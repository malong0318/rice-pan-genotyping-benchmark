#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/data/malong/benchmark/truth
REF=/data/malong/benchmark/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
REF_FAI=${REF}.fai
BCFTOOLS=${BCFTOOLS:-/data/malong/miniconda3/envs/benchmark/bin/bcftools}
TABIX=${TABIX:-/data/malong/miniconda3/bin/tabix}
BGZIP=${BGZIP:-/data/malong/miniconda3/bin/bgzip}

SAMPLES=(9311 IR64 W0169 W1735 W3033 zhonghua11)
TYPES=(snp indel sv)
MAP=${ROOT}/truth_to_chr.tsv
LOG=${ROOT}/prepare_truth_chr_vcfs.log
SUMMARY=${ROOT}/prepare_truth_chr_vcfs.summary.tsv

timestamp() { date "+%Y-%m-%dT%H:%M:%S%z"; }
log() { printf "[%s] %s\n" "$(timestamp)" "$*" | tee -a "$LOG"; }
die() { log "ERROR: $*"; exit 1; }
on_error() {
    local code=$?
    log "FAILED: line=${BASH_LINENO[0]:-unknown}, command=${BASH_COMMAND}, exit=${code}"
    exit "$code"
}
trap on_error ERR

[[ -x "$BCFTOOLS" ]] || die "bcftools not executable: $BCFTOOLS"
[[ -x "$TABIX" ]] || die "tabix not executable: $TABIX"
[[ -x "$BGZIP" ]] || die "bgzip not executable: $BGZIP"
[[ -s "$REF" ]] || die "reference missing or empty: $REF"
[[ -s "$REF_FAI" ]] || die "reference index missing: $REF_FAI"

mkdir -p "$ROOT"
if [[ ! -s "$MAP" ]]; then
    awk -F "\t" {'print "NIP-T2T#1#"$1"\t"$1'} "$REF_FAI" > "$MAP.tmp.$$"
    mv "$MAP.tmp.$$" "$MAP"
fi
for c in Chr{1..12}; do
    grep -F "NIP-T2T#1#$c" "$MAP" >/dev/null || die "chromosome mapping missing: $c"
done

printf "sample\ttype\tinput_records\tfinal_records\tinput\tfinal\n" > "$SUMMARY.tmp.$$"
log "START truth VCF preparation"
log "reference=$REF"
log "bcftools=$($BCFTOOLS --version | head -n 1)"
log "samples=${#SAMPLES[@]}, types=${#TYPES[@]}, total=18"

for sample in "${SAMPLES[@]}"; do
    outdir="$ROOT/$sample/Chr_vcf"
    mkdir -p "$outdir"
    mkdir -p "$outdir/.sort_tmp"
    for type in "${TYPES[@]}"; do
        prefix="$sample.$type"
        input="$ROOT/$sample/final/$prefix.vcf.gz"
        if [[ "$type" == "sv" ]]; then
            stem="$outdir/$prefix.bndfixed"
        else
            stem="$outdir/$prefix"
        fi
        renamed="$stem.renamed.vcf.gz"
        reheadered="$stem.reheadered.vcf.gz"
        sorted="$stem.sorted.vcf.gz"
        normed="$outdir/$prefix.norm.vcf.gz"
        final="$outdir/$prefix.final.vcf.gz"
        [[ -s "$input" ]] || die "input missing or empty: $input"
        log "[$sample/$type] input=$input"
        in_count=$($BCFTOOLS view -H "$input" | wc -l)
        log "[$sample/$type] input_records=$in_count"

        if [[ ! -s "$renamed" ]]; then
            log "[$sample/$type] 1/7 rename chromosome names"
            if [[ "$type" == "sv" ]]; then
                raw="$renamed.raw.tmp.$$"
                "$BCFTOOLS" annotate --rename-chrs "$MAP" -Oz -o "$raw" "$input"
                "$BGZIP" -cd "$raw" | sed "s/NIP-T2T#1#Chr/Chr/g" | "$BGZIP" -c > "$renamed.tmp.$$"
                rm -f "$raw"
                mv "$renamed.tmp.$$" "$renamed"
            else
                "$BCFTOOLS" annotate --rename-chrs "$MAP" -Oz -o "$renamed.tmp.$$" "$input"
                mv "$renamed.tmp.$$" "$renamed"
            fi
        else
            log "[$sample/$type] reuse renamed"
        fi

        if [[ ! -s "$reheadered" ]]; then
            log "[$sample/$type] 2/7 rewrite reference contig header"
            "$BCFTOOLS" reheader -f "$REF_FAI" -o "$reheadered.tmp.$$" "$renamed"
            mv "$reheadered.tmp.$$" "$reheadered"
        else
            log "[$sample/$type] reuse reheadered"
        fi

        if [[ ! -s "$sorted" ]]; then
            log "[$sample/$type] 3/7 sort records"
            "$BCFTOOLS" sort -T "$outdir/.sort_tmp" -Oz -o "$sorted.tmp.$$" "$reheadered"
            mv "$sorted.tmp.$$" "$sorted"
        else
            log "[$sample/$type] reuse sorted"
        fi
        if [[ ! -s "$sorted.tbi" ]]; then
            log "[$sample/$type] 4/7 index sorted VCF"
            "$TABIX" -f -p vcf "$sorted"
        fi

        if [[ ! -s "$normed" ]]; then
            if [[ "$type" == "sv" ]]; then
                log "[$sample/$type] 5/7 check SV anchors and normalize; no extra filtering"
                "$BCFTOOLS" norm -f "$REF" -c s -m -any -Oz -o "$normed.tmp.$$" "$sorted"
            else
                log "[$sample/$type] 5/7 strict REF check and normalize; no extra filtering"
                "$BCFTOOLS" norm -f "$REF" -c e -m -any -Oz -o "$normed.tmp.$$" "$sorted"
            fi
            mv "$normed.tmp.$$" "$normed"
        else
            log "[$sample/$type] reuse normalized"
        fi

        if [[ ! -s "$final" ]]; then
            log "[$sample/$type] 6/7 final sort"
            "$BCFTOOLS" sort -T "$outdir/.sort_tmp" -Oz -o "$final.tmp.$$" "$normed"
            mv "$final.tmp.$$" "$final"
        else
            log "[$sample/$type] reuse final"
        fi
        if [[ ! -s "$final.tbi" ]]; then
            log "[$sample/$type] 7/7 index final VCF"
            "$TABIX" -f -p vcf "$final"
        fi

        final_count=$($BCFTOOLS view -H "$final" | wc -l)
        [[ "$final_count" -gt 0 ]] || die "final VCF has no records: $final"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$type" "$in_count" "$final_count" "$input" "$final" >> "$SUMMARY.tmp.$$"
        log "[$sample/$type] DONE final_records=$final_count final=$final"
    done
done

mv "$SUMMARY.tmp.$$" "$SUMMARY"
log "SUCCESS: all 18 truth VCFs completed"
log "summary=$SUMMARY"
