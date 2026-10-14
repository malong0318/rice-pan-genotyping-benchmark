#!/usr/bin/env bash
set -Eeuo pipefail

BASE=/data/malong/benchmark
if printenv BENCHMARK_BASE >/dev/null 2>&1; then BASE=$(printenv BENCHMARK_BASE); fi
PAPER=$BASE/paper
BIN=$BASE/../miniconda3/envs/benchmark/bin
REF=$BASE/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
BCF=$BIN/bcftools
PY=$BIN/python
mkdir -p "$PAPER"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
size_bytes() { stat -c '%s' "$1" 2>/dev/null || printf 'NA'; }
vcf_records() {
  local f=$1 n
  [[ -s "$f" ]] || { printf 'NA'; return; }
  if [[ -s "$f.tbi" || -s "$f.csi" ]]; then
    n=$("$BCF" index -n "$f" 2>/dev/null || true)
    [[ "$n" =~ ^[0-9]+$ ]] && { printf '%s' "$n"; return; }
  fi
  n=$("$BCF" view -H "$f" 2>/dev/null | wc -l | tr -d ' ' || true)
  [[ "$n" =~ ^[0-9]+$ ]] && printf '%s' "$n" || printf 'NA'
}
vcf_contigs() {
  local f=$1 n
  [[ -s "$f" ]] || { printf 'NA'; return; }
  n=$("$BCF" view -h "$f" 2>/dev/null | awk '/^##contig=/{n++} END{print n+0}' || true)
  [[ "$n" =~ ^[0-9]+$ ]] && printf '%s' "$n" || printf 'NA'
}
version_of() {
  local exe=$1 out
  [[ -x "$exe" ]] || { printf 'not_found'; return; }
  out=$("$exe" --version 2>&1 || true)
  [[ -n "$out" ]] || out=$("$exe" version 2>&1 || true)
  printf '%s\n' "$out" | sed -n '1p'
}
emit_candidate() {
  local graph=$1 tool=$2 kind=$3 description=$4 f=$5
  local status=missing bytes=NA records=NA contigs=NA
  if [[ -s "$f" ]]; then
    status=present; bytes=$(size_bytes "$f"); records=$(vcf_records "$f"); contigs=$(vcf_contigs "$f")
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$graph" "$tool" "$kind" "$status" "$f" "$bytes" "$records" "$contigs" "$description" >> "$CAND_TMP"
}
emit_truth() {
  local sample=$1 dir=$BASE/truth/$sample/Chr_vcf
  local snp=$dir/$sample.snp.confident.final.vcf.gz indel=$dir/$sample.indel.small.confident.final.vcf.gz
  local sv=$dir/$sample.sv.truth.PASS.noBND.vcf.gz bed=$dir/$sample.confident.Chr.bed assembly
  case "$sample" in
    W0169|W1735|W3033) assembly=$BASE/graph/pggb/rice36_by_sample/$sample.fa ;;
    9311) assembly=$BASE/assembly/9311/9311_pansn.fa ;;
    IR64) assembly=$BASE/assembly/IR64/IR64_pansn.fa ;;
    zhonghua11) assembly=$BASE/assembly/zhonghua11/zhonghua11_pansn.fa ;;
  esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" assembly-derived 'dipcall+confident FILTER=.' "$snp" "$(vcf_records "$snp")" 'dipcall+confident FILTER=.' "$indel" "$(vcf_records "$indel")" SVIM-asm_haploid_PASS_noBND "$sv" "$(vcf_records "$sv")" "$bed" "$assembly" >> "$TRUTH_TMP"
}
emit_eval() {
  local tool=$1 graph=$2 sample=$3 depth=$4 klass=$5
    local truth=$BASE/truth/$sample/Chr_vcf/$sample.$klass.confident.final.vcf.gz bed=$BASE/truth/$sample/Chr_vcf/$sample.confident.Chr.bed
    [[ "$klass" == indel ]] && truth=$BASE/truth/$sample/Chr_vcf/$sample.indel.small.confident.final.vcf.gz
  local call=NA status=present reason=none call_n=NA
  [[ "$klass" == sv ]] && truth=$BASE/truth/$sample/Chr_vcf/$sample.sv.truth.PASS.noBND.vcf.gz
  if [[ "$tool" == vg && "$graph" == PGGB ]]; then
    status=NA; reason='VG-PGGB VG index construction failed'
  elif [[ "$tool" == pangenie && "$graph" == minigraph && "$klass" != sv ]]; then
    status=NA; reason='minigraph PanGenie benchmark scope is SV only'
  elif [[ "$tool" == graphtyper && "$graph" == minigraph && "$klass" != sv ]]; then
    status=NA; reason='minigraph GraphTyper benchmark scope is SV only'
  elif [[ "$tool" == bayestyper && "$graph" == minigraph && "$klass" != sv ]]; then
    status=NA; reason='minigraph BayesTyper benchmark scope is SV only'
  else
    case "$tool" in
      vg) call=$BASE/results_vg2/$sample/$depth/$graph/vg/$sample.$graph.vg.PASS.$klass.vcf.gz ;;
      pangenie) call=$BASE/results_pangenie2/$sample/$depth/$graph/pangenie/$sample.$graph.pangenie.PASS.$klass.vcf.gz ;;
      graphtyper) call=$BASE/results_graphtyper3/$sample/$depth/$graph/graphtyper/$sample.$graph.graphtyper.PASS.$klass.vcf.gz ;;
      bayestyper) call=$BASE/results_bayestyper2/$sample/$depth/$graph/bayestyper/$sample.$graph.bayestyper.PASS.$klass.vcf.gz ;;
    esac
    if [[ -s "$call" ]]; then call_n=$(vcf_records "$call"); else status=missing; reason='expected final evaluation VCF is absent'; fi
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$graph" "$sample" "$depth" "$klass" "$truth" "$(vcf_records "$truth")" "$bed" "$call" "$call_n" "$status" "$reason" "$BASE/evaluation/$klass/$tool/$graph/$sample/$depth" >> "$EVAL_TMP"
}

[[ -x "$BCF" ]] || die "bcftools not found: $BCF"
[[ -x "$PY" ]] || die "python not found: $PY"
[[ -s "$REF" ]] || die "reference missing: $REF"

VERSIONS=$PAPER/software_versions.txt
{
  printf '# Software versions recorded on %s\n' "$(date -Is)"
  printf '# Benchmark conda environment: %s\n\n' "$BASE/../miniconda3/envs/benchmark"
  printf 'software\tpath\tversion\n'
  printf 'bcftools\t%s\t%s\n' "$BCF" "$(version_of "$BCF")"
  printf 'samtools\t%s\t%s\n' "$BIN/samtools" "$(version_of "$BIN/samtools")"
  printf 'minimap2\t%s\t%s\n' "$BIN/minimap2" "$(version_of "$BIN/minimap2")"
  printf 'bedtools\t%s\t%s\n' "$BIN/bedtools" "$(version_of "$BIN/bedtools")"
  printf 'vg\t%s\tv1.73.0 (use vg version)\n' "$BIN/vg"
  printf 'graphtyper\t%s\t2.7.7\n' "$BIN/graphtyper"
  printf 'PanGenie\t%s\t4.2.1\n' "$BIN/PanGenie"
  printf 'pggb\t%s\t%s\n' "$BIN/pggb" "$(version_of "$BIN/pggb")"
  printf 'minigraph\t%s\t%s\n' "$BIN/minigraph" "$(version_of "$BIN/minigraph")"
  printf 'gfatools\t%s\t%s\n' "$BIN/gfatools" "$(version_of "$BIN/gfatools")"
  printf 'rtg\t%s\t3.13\n' "$BIN/rtg"
  printf 'truvari\t%s\t5.4.0\n' "$BIN/python -m truvari"
  printf 'svim-asm\t%s\t%s\n' "$BIN/svim-asm" "$(version_of "$BIN/svim-asm")"
  printf 'kmc\t%s\t3.2.4 (conda binary)\n' "$BIN/kmc"
  printf 'BayesTyper\t%s\tv1.5 local binary\n' "$BASE/code/04bayestyper/software/bayesTyper-debug"
  printf 'BayesTyperTools\t%s\tv1.5 local binary\n' "$BASE/code/04bayestyper/software/bayesTyperTools-debug"
  printf 'KMC-local\t%s\t2.3.0 local binary\n' "$BASE/code/04bayestyper/software/kmc-2.3/bin/kmc"
  printf 'dipcall\t%s\t0.3 as recorded by truth workflow\n' run-dipcall
  printf '\n# Existing graph files are historical outputs; current environment versions and recommended commands are in graph_statistics.tsv and evaluation_parameters.txt.\n'
} > "$VERSIONS"

DESIGN=$PAPER/benchmark_design.tsv
printf 'study_id\tspecies\treference\tgraph\tsample\tdepth\tread_type\tploidy_model\tgraph_membership\tsample_source\tvg\tpangenie\tgraphtyper\tbayestyper\n' > "$DESIGN"
for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do
  membership=graph_out; [[ "$sample" == W0169 || "$sample" == W1735 || "$sample" == W3033 ]] && membership=graph_in
  case "$sample" in
    9311) source='NGDC_GSA_CRA023020; NGDC_GWH_Assembly_88980' ;;
    IR64) source='NCBI_Assembly_GCA_055690995.1_NARO_Osi_IR64_r1.0; BioProject_PRJNA1311239' ;;
    zhonghua11) source='Complete_genome_assembly_of_japonica_rice_variety_Zhonghua_11' ;;
    W0169|W1735|W3033) source='A_pangenome_reference_of_wild_and_cultivated_rice' ;;
  esac
  for depth in 5x 10x 15x 20x; do
    for graph in MC PGGB minigraph; do
      vg_status=run; pg_status=run; gt_status=run; bt_status=run
      [[ "$graph" == PGGB ]] && vg_status='NA: VG PGGB index failed'
      [[ "$graph" == minigraph ]] && pg_status='SV only'; [[ "$graph" == minigraph ]] && gt_status='SV only'; [[ "$graph" == minigraph ]] && bt_status='SV only'
      printf 'rice36_graph_genotyping_benchmark\trice\t%s\t%s\t%s\t%s\tIllumina_paired_end\tpseudo_diploid_single_consensus\t%s\t%s\t%s\t%s\t%s\t%s\n' "$REF" "$graph" "$sample" "$depth" "$membership" "$source" "$vg_status" "$pg_status" "$gt_status" "$bt_status" >> "$DESIGN"
    done
  done
done

READS=$PAPER/sample_read_manifest.tsv
printf 'sample\tdepth\tR1\tR2\tmetadata\ttarget_depth_x\tachieved_depth_x\tselected_pairs\tread_count\tread_length_bp\testimated_bases_bp\tR1_bytes\tR2_bytes\tselection_method\n' > "$READS"
for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do
  for depth in 5x 10x 15x 20x; do
    dir=$BASE/Illumina_depth/$sample/$depth; r1=$dir/$sample'_'$depth'_R1.fastq.gz'; r2=$dir/$sample'_'$depth'_R2.fastq.gz'; meta=$dir/$sample'_'$depth'.metadata.tsv'
    target=$(awk -F '\t' '$1=="target_depth_x"{print $2}' "$meta" 2>/dev/null || true); achieved=$(awk -F '\t' '$1=="achieved_depth_x"{print $2}' "$meta" 2>/dev/null || true); pairs=$(awk -F '\t' '$1=="selected_pairs"{print $2}' "$meta" 2>/dev/null || true)
    [[ -n "$target" ]] || target=NA; [[ -n "$achieved" ]] || achieved=NA; [[ -n "$pairs" ]] || pairs=NA
    readlen=NA
    if [[ -s "$r1" ]]; then
      readlen=$($PY - "$r1" <<'PY'
import gzip
import sys
with gzip.open(sys.argv[1], 'rt') as handle:
    handle.readline()
    print(len(handle.readline().strip()))
PY
)
    fi
    readcount=NA; bases=NA
    if [[ "$pairs" =~ ^[0-9]+$ && "$readlen" =~ ^[0-9]+$ ]]; then readcount=$((pairs * 2)); bases=$((readcount * readlen)); fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$depth" "$r1" "$r2" "$meta" "$target" "$achieved" "$pairs" "$readcount" "$readlen" "$bases" "$(size_bytes "$r1")" "$(size_bytes "$r2")" seqtk_deterministic_paired_subsampling >> "$READS"
  done
done

GRAPH_STATS=$PAPER/graph_statistics.tsv
printf 'graph\tsource\tn_input_genomes\tgraph_file\tgraph_bytes\tgfa_nodes\tgfa_edges\tgraph_vcf_source\tgraph_vcf_records\tdeconstruct_or_projection\trecommended_current_command\tstatus\tnote\n' > "$GRAPH_STATS"
gfa_counts() {
  local f=$1
  if [[ "$f" == *.gfa.gz && -s "$f" ]]; then
    gzip -cd -- "$f" | awk -F '\t' '$1=="S"{n++} $1=="L"{e++} END{printf "%s\t%s", n+0, e+0}'
  elif [[ "$f" == *.gfa && -s "$f" ]]; then
    awk -F '\t' '$1=="S"{n++} $1=="L"{e++} END{printf "%s\t%s", n+0, e+0}' "$f"
  else
    printf 'NA\tNA'
  fi
}
mc_graph=$BASE/graph/cactus_new/rice36_mc.gfa.gz; pggb_graph=$BASE/graph/pggb/rice36_graph/rice36.complete.smooth.final.gfa; mini_graph=$BASE/graph/minigraph/rice36.gfa
mc_vcf=$BASE/graph/cactus_new/rice36_mc.vcf.gz; pggb_vcf=$BASE/graph/pggb/rice36_graph/rice36.complete.smooth.final.NIP-T2T.vcf.gz; mini_vcf=$BASE/graph/minigraph/rice36.deconstruct/rice36.sv.vcf
mc_counts=$(gfa_counts "$mc_graph"); pggb_counts=$(gfa_counts "$pggb_graph"); mini_counts=$(gfa_counts "$mini_graph")
mc_nodes=$(printf '%s\n' "$mc_counts" | cut -f1); mc_edges=$(printf '%s\n' "$mc_counts" | cut -f2); pggb_nodes=$(printf '%s\n' "$pggb_counts" | cut -f1); pggb_edges=$(printf '%s\n' "$pggb_counts" | cut -f2); mini_nodes=$(printf '%s\n' "$mini_counts" | cut -f1); mini_edges=$(printf '%s\n' "$mini_counts" | cut -f2)
printf 'MC\tminigraph-cactus_cactus-pangenome\t36\t%s\t%s\t%s\t%s\t%s\t%s\tMC historical VCF retained; old exact projection command was not preserved\tcactus-pangenome JOBSTORE GENOMES.tsv --outDir mc_out --outName rice36_mc --reference NIP-T2T --giraffe --gfa --vcf\tpresent\tHistorical MC compressed GFA retained; graph_bytes is the on-disk .gfa.gz size.\n' "$mc_graph" "$(size_bytes "$mc_graph")" "$mc_nodes" "$mc_edges" "$mc_vcf" "$(vcf_records "$mc_vcf")" >> "$GRAPH_STATS"
printf 'PGGB\tPGGB\t36\t%s\t%s\t%s\t%s\t%s\t%s\tvg convert -g -p -t 128 merged.gfa > merged.pg; vg deconstruct -p NIP-T2T#1#Chr1 ... -p NIP-T2T#1#Chr12 -t 128 merged.pg\tpggb -i rice36.complete.fa -o rice36_pggb -n 36 -t 128 -p 90 -s 10k -l 5k -k 23 -w 10k\tpresent\tHistorical PGGB output retained; command template records the current workflow.\n' "$pggb_graph" "$(size_bytes "$pggb_graph")" "$pggb_nodes" "$pggb_edges" "$pggb_vcf" "$(vcf_records "$pggb_vcf")" >> "$GRAPH_STATS"
printf 'minigraph\tminigraph_rGFA\t36\t%s\t%s\t%s\t%s\t%s\t%s\tvg convert -g -r 35 -p -t 128 rice36.gfa > rice36.minigraph.allpaths.pg; vg deconstruct -p NIP-T2T#1#Chr1 ... -p NIP-T2T#1#Chr12 -t 128 rice36.minigraph.allpaths.pg\tminigraph -cxggs -t 64 NIP-T2T.fa 35_other_assemblies.fa > rice36.gfa\tpresent\tThe supplied VCF is gfatools bubble output; symbolic alleles were sequence-repaired downstream.\n' "$mini_graph" "$(size_bytes "$mini_graph")" "$mini_nodes" "$mini_edges" "$mini_vcf" "$(vcf_records "$mini_vcf")" >> "$GRAPH_STATS"

TRUTH=$PAPER/truth_manifest.tsv; TRUTH_TMP=$TRUTH.tmp.$$
printf 'sample\ttruth_kind\tsnp_method\tsnp_vcf\tsnp_records\tindel_method\tindel_vcf\tindel_records\tsv_method\tsv_vcf\tsv_records\tconfident_bed\tassembly\n' > "$TRUTH_TMP"
for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do emit_truth "$sample"; done
mv "$TRUTH_TMP" "$TRUTH"

CAND=$PAPER/candidate_vcf_statistics.tsv; CAND_TMP=$CAND.tmp.$$
printf 'graph\ttool\tcandidate_kind\tstatus\tpath\tbytes\trecords\tcontigs\tdescription\n' > "$CAND_TMP"
emit_candidate MC graph_panel all 'historical MC source VCF' "$mc_vcf"
emit_candidate MC pangenie panel 'PanGenie normalized candidate panel' "$BASE/graph/cactus_new/pangenie_work/vcf/rice36_mc.pangenie.panel.vcf.gz"
emit_candidate MC graphtyper all_levels 'GraphTyper capacity-compatible all-level panel' "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.all_levels.vcf.gz"
emit_candidate MC graphtyper snp_smallindel 'GraphTyper SNP and small-indel panel' "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype.snp_smallindel.vcf.gz"
emit_candidate MC graphtyper sv 'GraphTyper SV panel' "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype_sv.vcf.gz"
emit_candidate PGGB graph_panel all 'historical PGGB smooth-final source VCF' "$pggb_vcf"
emit_candidate PGGB pangenie panel 'PanGenie normalized candidate panel' "$BASE/graph/pggb/pangenie_work/vcf/rice36_pggb.pangenie.panel.vcf.gz"
emit_candidate PGGB graphtyper all_levels 'GraphTyper capacity-compatible all-level panel' "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.all_levels.vcf.gz"
emit_candidate PGGB graphtyper snp_smallindel 'GraphTyper SNP and small-indel panel' "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype.snp_smallindel.vcf.gz"
emit_candidate PGGB graphtyper sv 'GraphTyper SV panel' "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype_sv.vcf.gz"
emit_candidate minigraph graph_panel sv 'gfatools bubble VCF from minigraph rGFA' "$mini_vcf"
emit_candidate minigraph graph_panel sequence_resolved 'sequence-resolved repaired minigraph panel' "$BASE/graph/minigraph/rice36.deconstruct/repaired/rice36.minigraph.sequence.vcf.gz"
emit_candidate minigraph pangenie panel 'PanGenie normalized candidate panel' "$BASE/graph/minigraph/pangenie_work/vcf/rice36_minigraph.pangenie.panel.vcf.gz"
emit_candidate minigraph graphtyper all_levels 'GraphTyper capacity-compatible all-level panel' "$BASE/graph/minigraph/graphtyper_work/rice36.minigraph.graphtyper.all_levels.vcf.gz"
emit_candidate minigraph graphtyper sv 'GraphTyper SV panel' "$BASE/graph/minigraph/graphtyper_work/rice36.minigraph.graphtyper.genotype_sv.vcf.gz"
printf 'PGGB\tvg\tall\tNA\tNA\tNA\tNA\tNA\tVG-PGGB omitted because the PGGB VG index construction failed; this is not zero calls.\n' >> "$CAND_TMP"
mv "$CAND_TMP" "$CAND"

EVAL=$PAPER/evaluation_input_manifest.tsv; EVAL_TMP=$EVAL.tmp.$$
printf 'tool\tgraph\tsample\tdepth\tvariant_class\ttruth_vcf\ttruth_records\tconfident_bed\tcall_vcf\tcall_records\tstatus\treason\texpected_evaluation_output\n' > "$EVAL_TMP"
for tool in vg pangenie graphtyper bayestyper; do
  for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do
    for depth in 5x 10x 15x 20x; do
      for graph in MC PGGB minigraph; do
        for klass in snp indel sv; do emit_eval "$tool" "$graph" "$sample" "$depth" "$klass"; done
      done
    done
  done
done
mv "$EVAL_TMP" "$EVAL"

PARAM=$PAPER/evaluation_parameters.txt
cat > "$PARAM" <<'EOF'
Benchmark evaluation parameters
================================
Reference: /data/malong/benchmark/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
Reference contigs: Chr1..Chr12; NIP-T2T#1#ChrN was renamed to ChrN.
Truth: assembly-derived. SNP/small indel from dipcall confident FILTER=. records; SV from SVIM-asm haploid, PASS and non-BND.
Genotype model: one consensus assembly was supplied twice to dipcall, therefore the truth is pseudo-diploid/homozygous rather than a phased diploid truth.

Unified classification after PASS, non-reference genotype preparation, normalization and decomposition:
SNP = explicit DNA REF/ALT, strlen(REF)=1, strlen(ALT)=1, REF!=ALT.
Small indel = explicit DNA REF/ALT, strlen(REF)!=strlen(ALT), both allele lengths <50.
SV = explicit DNA REF/ALT, abs(strlen(REF)-strlen(ALT))>=50.
Symbolic alleles and BND are excluded from sequence-based SNP/indel classification; SV truth excludes BND.

SNP/indel evaluator:
rtg vcfeval -t NIP-T2T_Chr.sdf -b truth.vcf.gz -c call.vcf.gz -e sample.confident.Chr.bed --sample syndip,SAMPLE --sample-ploidy 2 --decompose -o OUTPUT -T 64

SV evaluator:
python -m truvari bench -b truth.vcf.gz -c call.vcf.gz -f NIP-T2T_Chr.fa -o OUTPUT --passonly --bSample SAMPLE --cSample SAMPLE --includebed sample.confident.Chr.bed --sizemin 50 --sizefilt 50 --sizemax -1 --refdist 500 --pctseq 0.7 --pctsize 0.7 --bnddist 100 --no-ref c

Scope: MC and PGGB are evaluated for SNP, small indel and SV; minigraph is evaluated for SV only.
VG-PGGB is NA because construction of the PGGB VG index failed. It is not zero calls.

PGGB deconstruct: vg convert -g -p -t 128 rice36.pggb.12chr.merged.gfa > rice36.pggb.12chr.merged.pg; vg deconstruct -p NIP-T2T#1#Chr1 ... -p NIP-T2T#1#Chr12 -t 128 rice36.pggb.12chr.merged.pg | bgzip -@ 128 -c > rice36.pggb.vgdeconstruct.vcf.gz.
minigraph deconstruct: vg convert -g -r 35 -p -t 128 rice36.gfa > rice36.minigraph.allpaths.pg; vg deconstruct -p NIP-T2T#1#Chr1 ... -p NIP-T2T#1#Chr12 -t 128 rice36.minigraph.allpaths.pg | bgzip -@ 128 -c > rice36.minigraph.vgdeconstruct.vcf.gz.
MC deconstruct/projection: the historical minigraph-cactus/cactus-pangenome VCF is retained; the old exact build release was not preserved. A current cactus-pangenome template is in graph_statistics.tsv.

See software_versions.txt for executable versions and paths. Existing graph files are historical outputs and are not silently relabeled as newly rebuilt graphs.

Sample and graph membership provenance
---------------------------------------
The 36 graph genomes are from the study "A pangenome reference of wild and cultivated rice".
Graph-in benchmark samples: W0169, W1735, W3033.
Graph-out benchmark samples: 9311, IR64, zhonghua11.
9311 sources: NGDC GSA CRA023020 (https://ngdc.cncb.ac.cn/gsa/browse/CRA023020) and NGDC GWH Assembly 88980 (https://ngdc.cncb.ac.cn/gwh/Assembly/88980/show).
IR64 source: NCBI Assembly GCA_055690995.1, NARO_Osi_IR64_r1.0, BioProject PRJNA1311239 (https://www.ncbi.nlm.nih.gov/assembly/GCA_055690995.1/).
Zhonghua 11 source: "Complete genome assembly of japonica rice variety Zhonghua 11".
EOF

PERF=$PAPER/performance_environment.txt
{
  printf 'Performance measurement environment\n====================================\n'
  printf 'host\t%s\ndate\t%s\nkernel\t%s\nlogical_cpus\t%s\ncpu_model\t%s\nmemory_kb\t%s\nperformance_root\t%s\nmeasurement\t/usr/bin/time -v around complete genotyping workflow; performance.raw.txt is authoritative\n\n' "$(hostname)" "$(date -Is)" "$(uname -srmo)" "$(nproc)" "$(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ //')" "$(awk '/MemTotal:/{print $2}' /proc/meminfo)" "$BASE/results_performance"
  printf 'tool\tgraph\tsample\tdepth\tconfigured_threads\telapsed\telapsed_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\taverage_cpu_cores\tmax_rss_kb\tmax_rss_gb\texit_status\traw_log\n'
} > "$PERF"
"$PY" - "$BASE/results_performance" <<'PY' >> "$PERF"
import pathlib, re, sys
root = pathlib.Path(sys.argv[1])
def sec(value):
    p = value.split(':')
    try:
        if len(p) == 3: return float(p[0])*3600 + float(p[1])*60 + float(p[2])
        if len(p) == 2: return float(p[0])*60 + float(p[1])
        return float(value)
    except ValueError: return float('nan')
for raw in sorted(root.glob('*/*/9311/5x/performance.raw.txt')):
    tool, graph, sample, depth = raw.relative_to(root).parts[:4]
    text = raw.read_text(errors='replace')
    def field(pattern):
        m = re.search(pattern, text, re.MULTILINE)
        return m.group(1).strip() if m else 'NA'
    user = field(r'^\s*User time \(seconds\):\s*(\S+)')
    system = field(r'^\s*System time \(seconds\):\s*(\S+)')
    cpu = field(r'^\s*Percent of CPU this job got:\s*(\S+)%')
    elapsed = field(r'^\s*Elapsed.*\):\s*(\S+)$')
    rss = field(r'^\s*Maximum resident set size \(kbytes\):\s*(\S+)')
    status = field(r'^\s*Exit status:\s*(\S+)')
    try: avg = float(cpu)/100.0
    except ValueError: avg = float('nan')
    try: gb = float(rss)/(1024.0*1024.0)
    except ValueError: gb = float('nan')
    print('\t'.join([tool, graph, sample, depth, '16', elapsed, f'{sec(elapsed):.3f}', user, system, cpu, f'{avg:.3f}', rss, f'{gb:.3f}', status, str(raw)]))
PY
printf '\nNotes\n-----\nVG-PGGB has no performance row because the PGGB VG index construction failed.\nCPU percent is GNU time output; average_cpu_cores is cpu_percent/100.\nThe performance test covered 9311 at 5x and 11 completed tool/graph combinations.\n' >> "$PERF"

CHECK=$PAPER/checksums.sha256
{
  printf '# SHA256 scope: metadata files, reference FASTA/index, graph files, final truth VCFs and candidate VCFs.\n'
  printf '# Evaluation call VCFs are listed in evaluation_input_manifest.tsv and are not hashed here because they are large and numerous.\n'
} > "$CHECK"
checksum_if_present() { local f=$1; [[ -f "$f" ]] || return 0; sha256sum "$f" >> "$CHECK"; }
for f in "$VERSIONS" "$DESIGN" "$READS" "$GRAPH_STATS" "$TRUTH" "$CAND" "$EVAL" "$PARAM" "$PERF"; do checksum_if_present "$f"; done
for f in "$REF" "$REF.fai" "$mc_graph" "$pggb_graph" "$mini_graph"; do checksum_if_present "$f"; done
for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do
  for f in "$BASE/truth/$sample/Chr_vcf/$sample.snp.confident.final.vcf.gz" "$BASE/truth/$sample/Chr_vcf/$sample.snp.confident.final.vcf.gz.tbi" "$BASE/truth/$sample/Chr_vcf/$sample.indel.small.confident.final.vcf.gz" "$BASE/truth/$sample/Chr_vcf/$sample.indel.small.confident.final.vcf.gz.tbi" "$BASE/truth/$sample/Chr_vcf/$sample.sv.truth.PASS.noBND.vcf.gz" "$BASE/truth/$sample/Chr_vcf/$sample.sv.truth.PASS.noBND.vcf.gz.tbi"; do checksum_if_present "$f"; done
done
for f in "$mc_vcf" "$BASE/graph/cactus_new/pangenie_work/vcf/rice36_mc.pangenie.panel.vcf.gz" "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.all_levels.vcf.gz" "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype.snp_smallindel.vcf.gz" "$BASE/graph/cactus_new/graphtyper_work/rice36_mc.graphtyper.genotype_sv.vcf.gz" "$pggb_vcf" "$BASE/graph/pggb/pangenie_work/vcf/rice36_pggb.pangenie.panel.vcf.gz" "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.all_levels.vcf.gz" "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype.snp_smallindel.vcf.gz" "$BASE/graph/pggb/graphtyper_work/rice36.pggb.graphtyper.genotype_sv.vcf.gz" "$mini_vcf" "$BASE/graph/minigraph/rice36.deconstruct/repaired/rice36.minigraph.sequence.vcf.gz" "$BASE/graph/minigraph/pangenie_work/vcf/rice36_minigraph.pangenie.panel.vcf.gz" "$BASE/graph/minigraph/graphtyper_work/rice36.minigraph.graphtyper.all_levels.vcf.gz" "$BASE/graph/minigraph/graphtyper_work/rice36.minigraph.graphtyper.genotype_sv.vcf.gz"; do
  checksum_if_present "$f"; checksum_if_present "$f.tbi"; checksum_if_present "$f.csi"
done

log 'metadata generation completed'
printf '%s\n' "$VERSIONS" "$DESIGN" "$READS" "$GRAPH_STATS" "$TRUTH" "$CAND" "$EVAL" "$PARAM" "$PERF" "$CHECK"
