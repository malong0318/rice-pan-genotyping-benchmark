#!/usr/bin/env bash
set -Eeuo pipefail
ulimit -c 0
ROOT=/data/malong/benchmark
BIN=/data/malong/miniconda3/envs/benchmark/bin
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/performance_common.sh"
MODE=$(printenv PERF_MODE 2>/dev/null || true)
SAMPLE=$(printenv PERF_SAMPLE 2>/dev/null || true)
DEPTH=$(printenv PERF_DEPTH 2>/dev/null || true)
THREADS=$(printenv PERF_THREADS 2>/dev/null || true)
[ -n "$SAMPLE" ] || SAMPLE=9311
[ -n "$DEPTH" ] || DEPTH=5x
[ -n "$THREADS" ] || THREADS=128
case "$MODE" in
  vg_mc) tool=vg; graph=MC;;
  vg_minigraph) tool=vg; graph=minigraph;;
  pangenie_mc) tool=pangenie; graph=MC;;
  pangenie_pggb) tool=pangenie; graph=PGGB;;
  pangenie_minigraph) tool=pangenie; graph=minigraph;;
  graphtyper_mc) tool=graphtyper; graph=MC;;
  graphtyper_pggb) tool=graphtyper; graph=PGGB;;
  graphtyper_minigraph) tool=graphtyper; graph=minigraph;;
  bayestyper_mc) tool=bayestyper; graph=MC;;
  bayestyper_pggb) tool=bayestyper; graph=PGGB;;
  bayestyper_minigraph) tool=bayestyper; graph=minigraph;;
  *) perf_die 'invalid PERF_MODE';;
esac
OUT=$ROOT/results_performance/$tool/$graph/$SAMPLE/$DEPTH
perf_init "$OUT" "$0" "$@"
mkdir -p "$OUT"
exec > >(tee -a "$OUT/run.log") 2>&1
trap 'rc=$?; perf_log "FAILED line=$LINENO status=$rc"; exit "$rc"' ERR
PERF_START=$(date +%s)
REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
BCFTOOLS=$BIN/bcftools
perf_need_bin "$BCFTOOLS"
perf_need_file "$REF"
perf_need_file "$REF.fai"

if [ "$tool" = vg ]; then
  VG=$BIN/vg
  perf_need_bin "$VG"
  R1=$ROOT/Illumina_depth/$SAMPLE/$DEPTH/"$SAMPLE"_"$DEPTH"_R1.fastq.gz
  R2=$ROOT/Illumina_depth/$SAMPLE/$DEPTH/"$SAMPLE"_"$DEPTH"_R2.fastq.gz
  perf_need_file "$R1"; perf_need_file "$R2"
  if [ "$graph" = MC ]; then
    GRAPH=$ROOT/graph/cactus_new/rice36_mc.d2.gbz
    DIST=$ROOT/graph/cactus_new/rice36_mc.d2.dist
    MIN=$ROOT/graph/cactus_new/vg173_giraffe_index/rice36_mc.d2.shortread.withzip.min
    ZIP=$ROOT/graph/cactus_new/vg173_giraffe_index/rice36_mc.d2.shortread.zipcodes
    REF_ARGS=''; for c in $(seq 1 12); do REF_ARGS="$REF_ARGS -p NIP-T2T#0#Chr$c"; done
  else
    GRAPH=$ROOT/graph/minigraph/vg173_giraffe_index/rice36.minigraph.giraffe.gbz
    DIST=$ROOT/graph/minigraph/vg173_giraffe_index/rice36.minigraph.dist
    MIN=$ROOT/graph/minigraph/vg173_giraffe_index/rice36.minigraph.shortread.withzip.min
    ZIP=$ROOT/graph/minigraph/vg173_giraffe_index/rice36.minigraph.shortread.zipcodes
    REF_ARGS=''; for c in $(seq 1 12); do REF_ARGS="$REF_ARGS -p NIP-T2T#1#Chr$c"; done
  fi
  for f in "$GRAPH" "$DIST" "$MIN" "$ZIP"; do perf_need_file "$f"; done
  GAM=$OUT/$SAMPLE.$graph.vg.giraffe.gam
  PACK=$OUT/$SAMPLE.$graph.vg.pack
  VCF=$OUT/$SAMPLE.$graph.vg.vcf
  RESULT=$OUT/$SAMPLE.$graph.vg.vcf.gz
  perf_log "START vg $graph: giraffe, pack, call"
  "$VG" giraffe -Z "$GRAPH" -d "$DIST" -m "$MIN" -z "$ZIP" -t "$THREADS" -f "$R1" -f "$R2" -N "$SAMPLE" > "$GAM"
  perf_need_file "$GAM"
  "$VG" pack -x "$GRAPH" -g "$GAM" -Q 10 -t "$THREADS" -o "$PACK"
  perf_need_file "$PACK"
  "$VG" call "$GRAPH" -k "$PACK" -s "$SAMPLE" -t "$THREADS" $REF_ARGS > "$VCF"
  grep -q '^#CHROM' "$VCF" || perf_die 'vg call produced no VCF header'
  "$BIN/bgzip" -@ "$THREADS" -c "$VCF" > "$RESULT"
  "$BCFTOOLS" index -f -t "$RESULT"
elif [ "$tool" = pangenie ]; then
  READS=$ROOT/Illumina_depth/pangenie_fastq/"$SAMPLE"_"$DEPTH".fastq
  graphdir=$(printf '%s' "$graph" | tr '[:upper:]' '[:lower:]')
  [ "$graph" = MC ] && graphdir=cactus_new
  WORK=$ROOT/graph/$graphdir/pangenie_work
  INDEX_DIR=$WORK/index_all; PANEL_DIR=$WORK/vcf/by_chromosome
  perf_need_file "$READS"
  mkdir -p "$OUT/chromosomes" "$OUT/logs"
  LIST=$OUT/chromosomes.list; : > "$LIST"
  for c in $(seq 1 12); do
    chr=Chr$c; prefix=$INDEX_DIR/$chr.pangenie
    perf_need_file "$PANEL_DIR/$chr.vcf"
    perf_need_file "$prefix"_UniqueKmersMap.cereal
    perf_need_file "$prefix"_path_segments.fasta
    perf_need_file "$prefix"_"$chr"_Graph.cereal
    outprefix=$OUT/chromosomes/$SAMPLE.$DEPTH.$graph.$chr
    raw=$outprefix"_genotyping.vcf"
    if [ ! -s "$raw" ]; then
      perf_log "PanGenie $graph $chr: kmer_threads=$THREADS core_threads=1"
      "$BIN/PanGenie" -f "$prefix" -i "$READS" -o "$outprefix" -s "$SAMPLE" -j "$THREADS" -t 1 -u > "$OUT/logs/$chr.log" 2>&1
    fi
    perf_need_file "$raw"
    compressed=$raw.gz
    "$BIN/bgzip" -@ "$THREADS" -c "$raw" > "$compressed"
    "$BIN/tabix" -f -p vcf "$compressed"
    printf '%s\n' "$compressed" >> "$LIST"
  done
  ALL=$OUT/$SAMPLE.$DEPTH.$graph.pangenie.vcf.gz
  "$BCFTOOLS" concat -a -Oz -o "$ALL" -f "$LIST"
  "$BIN/tabix" -f -p vcf "$ALL"
  if [ "$graph" = minigraph ]; then
    RESULT=$OUT/$SAMPLE.$DEPTH.$graph.pangenie.sv.vcf.gz
    cp -f "$ALL" "$RESULT"; cp -f "$ALL.tbi" "$RESULT.tbi"
  else
    for kind in snp indel sv; do
      file=$OUT/$SAMPLE.$DEPTH.$graph.pangenie.$kind.vcf.gz
      if [ "$kind" = snp ]; then
        "$BCFTOOLS" view -v snps -Oz -o "$file" "$ALL"
      elif [ "$kind" = indel ]; then
        "$BCFTOOLS" view -i 'strlen(REF)!=strlen(ALT) && abs(strlen(REF)-strlen(ALT))<50' -Oz -o "$file" "$ALL"
      else
        "$BCFTOOLS" view -i 'strlen(REF)!=strlen(ALT) && abs(strlen(REF)-strlen(ALT))>=50' -Oz -o "$file" "$ALL"
      fi
      "$BIN/tabix" -f -p vcf "$file"
    done
    RESULT=$ALL
  fi
fi

if [ "$tool" = vg ] || [ "$tool" = pangenie ]; then
  perf_need_file "$RESULT"
  perf_finish "$OUT" "$tool" "$graph" "$RESULT"
  exit 0
fi

if [ "$tool" = graphtyper ]; then
  GT=$BIN/graphtyper
  perf_need_bin "$GT"
  SHARED=$ROOT/share/graphtyper
  REGIONS=$SHARED/regions.10000000bp.tsv
  BAMDIR=$SHARED/bam/$SAMPLE/$DEPTH
  BAM=$BAMDIR/$SAMPLE.$DEPTH.MC.sorted.bam
  BAMLIST=$BAMDIR/bams.list
  for f in "$REGIONS" "$BAM" "$BAM.bai" "$BAMLIST"; do perf_need_file "$f"; done
  if [ "$graph" = MC ]; then
    VCF=$ROOT/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype.snp_smallindel.vcf.gz
    SVVCF=$ROOT/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype_sv.vcf.gz
  elif [ "$graph" = PGGB ]; then
    VCF=$ROOT/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype.snp_smallindel.vcf.gz
    SVVCF=$ROOT/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype_sv.vcf.gz
  else
    SVVCF=$ROOT/graph/minigraph/graphtyper_work/rice36.minigraph.graphtyper.genotype_sv.vcf.gz
  fi
  perf_need_file "$SVVCF"
  [ "$graph" = minigraph ] || perf_need_file "$VCF"
  [ -s "$SVVCF.tbi" ] || "$BIN/tabix" -f -p vcf "$SVVCF"
  mkdir -p "$OUT/genotype" "$OUT/genotype_sv"
  if [ "$graph" != minigraph ]; then
    "$GT" genotype "$REF" --vcf "$VCF" --sams "$BAMLIST" --output "$OUT/genotype" --region_file "$REGIONS" --threads "$THREADS" --no_decompose --force_no_copy_reference
    find "$OUT/genotype" -type f -name '*.vcf.gz' -size +0c | sort -V > "$OUT/genotype.list"
    [ -s "$OUT/genotype.list" ] || perf_die 'no GraphTyper genotype fragments'
    GT_FINAL=$OUT/$SAMPLE.$DEPTH.$graph.graphtyper.snp_smallindel.vcf.gz
    "$BCFTOOLS" concat -a --threads 8 -Oz -o "$GT_FINAL" -f "$OUT/genotype.list"
    "$BIN/tabix" -f -p vcf "$GT_FINAL"
  fi
  "$GT" genotype_sv "$REF" "$SVVCF" --sams "$BAMLIST" --output "$OUT/genotype_sv" --region_file "$REGIONS" --threads "$THREADS" --force_no_copy_reference
  find "$OUT/genotype_sv" -type f -name '*.vcf.gz' -size +0c | sort -V > "$OUT/genotype_sv.list"
  [ -s "$OUT/genotype_sv.list" ] || perf_die 'no GraphTyper SV fragments'
  RAW=$OUT/$SAMPLE.$DEPTH.$graph.graphtyper.sv.raw.vcf.gz
  "$BCFTOOLS" concat -a --threads 8 -Oz -o "$RAW" -f "$OUT/genotype_sv.list"
  "$BIN/tabix" -f -p vcf "$RAW"
  RESULT=$OUT/$SAMPLE.$DEPTH.$graph.graphtyper.sv.vcf.gz
  "$BCFTOOLS" view -f PASS -i 'INFO/SVMODEL="AGGREGATED"' -Oz -o "$RESULT" "$RAW"
  "$BIN/tabix" -f -p vcf "$RESULT"
  perf_need_file "$RESULT"
  perf_finish "$OUT" "$tool" "$graph" "$RESULT"
  exit 0
fi

if [ "$tool" = bayestyper ]; then
  BT_SOFTWARE=$ROOT/code/04bayestyper/software
  BT=$BT_SOFTWARE/bayesTyper-debug
  BTT=$BT_SOFTWARE/bayesTyperTools-debug
  KMC=$BT_SOFTWARE/kmc-2.3/bin/kmc
  perf_need_bin "$BT"
  perf_need_bin "$BTT"
  perf_need_bin "$KMC"
  # The bundled BayesTyper binaries require Boost 1.85 at runtime.  The
  # matching libraries are installed in the benchmark Conda environment;
  # the BayesTyper executables themselves come only from BT_SOFTWARE.
  export LD_LIBRARY_PATH="$BIN/../lib:$BT_SOFTWARE:$BT_SOFTWARE/kmc-2.3/bin:${LD_LIBRARY_PATH:-}"
  # BayesTyper candidate panels use the graph-style NIP-T2T#1#ChrN names.
  # Use the matching FASTA rather than the ChrN-renamed FASTA used by the
  # other genotypers.
  BT_REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T.fa
  perf_need_file "$BT_REF"
  perf_need_file "$BT_REF.fai"
  WORK=$OUT/work; LOG=$OUT/log; mkdir -p "$WORK" "$LOG"
  CAND_RAW=$ROOT/results/$SAMPLE/$DEPTH/work/vcf/nonvg.cactus.bayestyper.vcf
  if [ "$graph" = PGGB ]; then CAND_RAW=$ROOT/results/$SAMPLE/$DEPTH/work/vcf/nonvg.pggb.bayestyper.vcf; fi
  if [ "$graph" = minigraph ]; then CAND_RAW=$ROOT/graph/minigraph/rice36.deconstruct/repaired/rice36.minigraph.sequence.vcf.gz; fi
  graphdir=$(printf '%s' "$graph" | tr '[:upper:]' '[:lower:]')
  bayesdir=$graph
  [ "$graph" = MC ] && bayesdir=MC
  [ "$graph" = PGGB ] && bayesdir=pggb
  KMC_PREFIX=$ROOT/results/$SAMPLE/$DEPTH/$bayesdir/bayestyper/work/"$SAMPLE"_"$DEPTH".kmc
  for f in "$CAND_RAW" "$KMC_PREFIX.kmc_pre" "$KMC_PREFIX.kmc_suf" "$KMC_PREFIX.bloomMeta" "$KMC_PREFIX.bloomData"; do perf_need_file "$f"; done
  CAND=$WORK/candidate.vcf.gz
  "$BCFTOOLS" norm --threads "$THREADS" -f "$BT_REF" -m +any -Oz -o "$CAND" "$CAND_RAW"
  "$BIN/tabix" -f -p vcf "$CAND"
  SAMPLES=$WORK/samples.tsv
  printf '%s\tF\t%s\n' "$SAMPLE" "$KMC_PREFIX" > "$SAMPLES"
  prefix=$(printf '%s' "$graph" | tr '[:upper:]' '[:lower:]')
  CLUSTER=$WORK/"$prefix"_cluster_data
  "$BT" cluster -v "$CAND" -s "$SAMPLES" -g "$BT_REF" -o "$WORK/$prefix" -p"$THREADS" > "$LOG/cluster.log" 2>&1
  find "$WORK/$prefix"_unit_*/variant_clusters.bin -type f 2>/dev/null | sort -V > "$WORK/units.list" || true
  [ -s "$WORK/units.list" ] || perf_die 'no BayesTyper units'
  : > "$WORK/vcfs.list"
  while IFS= read -r unit; do
    udir=$(dirname "$unit")
    uprefix=$udir/bayestyper
    "$BT" genotype -v "$unit" -c "$CLUSTER" -s "$SAMPLES" -g "$BT_REF" -o "$uprefix" -p"$THREADS" -z > "$LOG/$(basename "$udir").log" 2>&1
    bgz=$udir/bayestyper.bgzf.vcf.gz
    "$BCFTOOLS" view -Oz -o "$bgz" "$uprefix.vcf.gz"
    "$BIN/tabix" -f -p vcf "$bgz"
    printf '%s\n' "$bgz" >> "$WORK/vcfs.list"
  done < "$WORK/units.list"
  RESULT=$OUT/$SAMPLE.$DEPTH.$graph.bayestyper.vcf.gz
  "$BCFTOOLS" concat -a -f "$WORK/vcfs.list" -Oz -o "$RESULT"
  "$BIN/tabix" -f -p vcf "$RESULT"
  perf_need_file "$RESULT"
  perf_finish "$OUT" "$tool" "$graph" "$RESULT"
  exit 0
fi

perf_die 'unreachable mode'
