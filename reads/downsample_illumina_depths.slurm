#!/usr/bin/env bash
# Create deterministic paired-end Illumina subsets at 20x, 15x, 10x and 5x.
# Coverage is calculated against the NIP-T2T reference used by the truth set.
# The original FASTQ files are never modified.
#
# Submit all samples:
#   sbatch /data/malong/benchmark/code/downsample_illumina_depths.slurm
# Submit one sample (0=9311, 1=IR64, 2=zhonghua11, 3=W0169, 4=W1735, 5=W3033):
#   sbatch --array=4 /data/malong/benchmark/code/downsample_illumina_depths.slurm

#SBATCH --job-name=rice_read_depth
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=7-00:00:00
#SBATCH --array=0-5
#SBATCH --output=/data/malong/benchmark/rice_read_depth_%A_%a.out
#SBATCH --error=/data/malong/benchmark/rice_read_depth_%A_%a.err

set -Eeuo pipefail

BASE=/data/malong/benchmark
READS="$BASE/Illumina"
OUTROOT="$BASE/Illumina_depth"
REF="$BASE/graph/pggb/rice36_by_sample/NIP-T2T.fa"
REF_FAI="$REF.fai"
SEQKIT=/data/malong/miniconda3/envs/benchmark/bin/seqkit
SEQTK=/data/malong/miniconda3/envs/benchmark/bin/seqtk
PIGZ=$(command -v pigz || true)
THREADS=${SLURM_CPUS_PER_TASK:-16}
TARGETS=(20 15 10 5)
SAMPLES=(9311 IR64 zhonghua11 W0169 W1735 W3033)

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]:-}
[[ -n "$sample" ]] || { echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-unset}" >&2; exit 2; }

R1="$READS/$sample/${sample}_1.fastq"
R2="$READS/$sample/${sample}_2.fastq"
[[ -s "$R1" ]] || { echo "Missing R1: $R1" >&2; exit 2; }
[[ -s "$R2" ]] || { echo "Missing R2: $R2" >&2; exit 2; }
[[ -x "$SEQKIT" ]] || { echo "Missing seqkit: $SEQKIT" >&2; exit 127; }
[[ -x "$SEQTK" ]] || { echo "Missing seqtk: $SEQTK" >&2; exit 127; }
[[ -n "$PIGZ" ]] || { echo "pigz is required for compressed output" >&2; exit 127; }

mkdir -p "$OUTROOT/$sample"
LOG="$OUTROOT/$sample/$sample.downsample.log"
exec > >(tee "$LOG") 2>&1

echo "[$(date -Is)] sample=$sample"
echo "[$(date -Is)] R1=$R1"
echo "[$(date -Is)] R2=$R2"

[[ -s "$REF_FAI" ]] || "$BASE/../miniconda3/envs/benchmark/bin/samtools" faidx "$REF"
REF_LEN=$(awk '{s+=$2} END {printf "%.0f\n", s}' "$REF_FAI")
[[ "$REF_LEN" -gt 0 ]] || { echo "Could not determine reference length" >&2; exit 1; }

# seqkit reports exact read count and total sequenced bases. Paired files are
# required to have the same number of records and the same read length sum.
STATS=$($SEQKIT stats -T -j "$THREADS" "$R1" "$R2" 2>/dev/null)
R1_N=$(awk 'NR==2 {print $4}' <<< "$STATS")
R1_BASES=$(awk 'NR==2 {print $5}' <<< "$STATS")
R2_N=$(awk 'NR==3 {print $4}' <<< "$STATS")
R2_BASES=$(awk 'NR==3 {print $5}' <<< "$STATS")
[[ "$R1_N" =~ ^[0-9]+$ && "$R2_N" =~ ^[0-9]+$ ]] || { echo "Unable to parse seqkit stats" >&2; echo "$STATS"; exit 1; }
[[ "$R1_N" == "$R2_N" ]] || { echo "R1/R2 record counts differ: $R1_N vs $R2_N" >&2; exit 1; }
TOTAL_BASES=$((R1_BASES + R2_BASES))
AVG_PAIR_BASES=$(awk -v b="$TOTAL_BASES" -v n="$R1_N" 'BEGIN {printf "%.6f", b/n}')
RAW_DEPTH=$(awk -v b="$TOTAL_BASES" -v r="$REF_LEN" 'BEGIN {printf "%.4f", b/r}')

echo "reference_length_bp=$REF_LEN"
echo "input_pairs=$R1_N"
echo "R1_bases=$R1_BASES"
echo "R2_bases=$R2_BASES"
echo "average_pair_bases=$AVG_PAIR_BASES"
echo "input_depth_x=$RAW_DEPTH"

for target in "${TARGETS[@]}"; do
  # Number of pairs is chosen from the actual average pair length, so the
  # achieved depth is close to the requested value even if read lengths vary.
  NPAIRS=$(awk -v t="$target" -v r="$REF_LEN" -v n="$R1_N" -v b="$TOTAL_BASES" \
    'BEGIN { x=t*r*n/b; printf "%.0f\n", int(x) }')
  (( NPAIRS > 0 )) || { echo "Computed zero pairs for ${target}x" >&2; exit 1; }
  (( NPAIRS <= R1_N )) || { echo "Requested ${target}x exceeds input depth (${RAW_DEPTH}x)" >&2; exit 1; }

  SEED=$((17000 + SLURM_ARRAY_TASK_ID * 100 + target))
  OD="$OUTROOT/$sample/${target}x"
  mkdir -p "$OD"
  O1="$OD/${sample}_${target}x_R1.fastq.gz"
  O2="$OD/${sample}_${target}x_R2.fastq.gz"
  META="$OD/${sample}_${target}x.metadata.tsv"

  if [[ ! -s "$O1" || ! -s "$O2" || ! -s "$META" ]]; then
    tmp1="$OD/.${sample}_${target}x_R1.fastq.gz.tmp"
    tmp2="$OD/.${sample}_${target}x_R2.fastq.gz.tmp"
    rm -f "$tmp1" "$tmp2"

    # seqtk uses the same deterministic RNG stream for two files with the
    # same record count and seed; sampling R1 and R2 separately therefore
    # preserves pairing without storing millions of read IDs in memory.
    "$SEQTK" sample -s"$SEED" "$R1" "$NPAIRS" \
      | "$PIGZ" -p "$THREADS" -c > "$tmp1"
    "$SEQTK" sample -s"$SEED" "$R2" "$NPAIRS" \
      | "$PIGZ" -p "$THREADS" -c > "$tmp2"
    mv "$tmp1" "$O1"
    mv "$tmp2" "$O2"

    ACHIEVED=$(awk -v n="$NPAIRS" -v b="$TOTAL_BASES" -v r="$R1_N" -v ref="$REF_LEN" \
      'BEGIN {printf "%.6f", (n*(b/r))/ref}')
    {
      printf 'sample\t%s\n' "$sample"
      printf 'target_depth_x\t%s\n' "$target"
      printf 'achieved_depth_x\t%s\n' "$ACHIEVED"
      printf 'input_R1\t%s\ninput_R2\t%s\n' "$R1" "$R2"
      printf 'output_R1\t%s\noutput_R2\t%s\n' "$O1" "$O2"
      printf 'reference\t%s\nreference_length_bp\t%s\n' "$REF" "$REF_LEN"
      printf 'input_pairs\t%s\nselected_pairs\t%s\n' "$R1_N" "$NPAIRS"
      printf 'seed\t%s\nmethod\tseqtk_sample_same_seed_on_R1_R2\n' "$SEED"
      printf 'sampling_note\tdeterministic independent subset for this target depth\n'
    } > "$META"
  else
    echo "[$(date -Is)] exists, skipping ${target}x"
  fi

  # Basic postcondition checks: both files must contain NPAIRS FASTQ records.
  C1=$(zcat "$O1" | awk 'END {print int(NR/4)}')
  C2=$(zcat "$O2" | awk 'END {print int(NR/4)}')
  [[ "$C1" == "$C2" && "$C1" == "$NPAIRS" ]] || {
    echo "Pair-count validation failed for ${target}x: R1=$C1 R2=$C2 expected=$NPAIRS" >&2
    exit 1
  }
  echo "[$(date -Is)] ${target}x complete: pairs=$C1 achieved_depth_x=$(awk -v n="$C1" -v b="$TOTAL_BASES" -v r="$R1_N" -v ref="$REF_LEN" 'BEGIN{printf "%.4f",n*(b/r)/ref}')"
done

echo "[$(date -Is)] completed sample=$sample"
