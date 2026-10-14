#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/data/malong/benchmark
OUT=$ROOT/01classify_vcf
REF=$ROOT/graph/pggb/rice36_by_sample/NIP-T2T_Chr.fa
BCF=/data/malong/miniconda3/envs/benchmark/bin/bcftools
TABIX=/data/malong/miniconda3/envs/benchmark/bin/tabix
MAP=$OUT/contig_to_chr.tsv
SUMMARY=$OUT/summary.tsv
MISSING=$OUT/missing_inputs.tsv

THREADS=$(printenv BAYESMARK_CLASSIFY_THREADS 2>/dev/null || true)
SORT_MEM=$(printenv BAYESMARK_CLASSIFY_SORT_MEM 2>/dev/null || true)
ONLY_TOOL=$(printenv BAYESMARK_ONLY_TOOL 2>/dev/null || true)
ONLY_SAMPLE=$(printenv BAYESMARK_ONLY_SAMPLE 2>/dev/null || true)
ONLY_DEPTH=$(printenv BAYESMARK_ONLY_DEPTH 2>/dev/null || true)
ONLY_GRAPH=$(printenv BAYESMARK_ONLY_GRAPH 2>/dev/null || true)
FORCE=$(printenv BAYESMARK_FORCE 2>/dev/null || true)
[ -n "$THREADS" ] || THREADS=8
[ -n "$SORT_MEM" ] || SORT_MEM=4G

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*"
}

die() {
    log "ERROR: $*" >&2
    exit 1
}

on_error() {
    rc=$?
    log "FAILED line=$1 command=$2 exit=$rc" >&2
    exit "$rc"
}
trap 'on_error "$LINENO" "$BASH_COMMAND"' ERR

command -v "$BCF" >/dev/null 2>&1 || die "bcftools not found: $BCF"
command -v "$TABIX" >/dev/null 2>&1 || die "tabix not found: $TABIX"
[ -s "$REF" ] || die "reference not found: $REF"
mkdir -p "$OUT" "$OUT/vg" "$OUT/pangenie" "$OUT/bayestyper"

if [ ! -s "$MAP" ]; then
    printf '%s\n' \
        'NIP-T2T#1#Chr1 Chr1' \
        'NIP-T2T#1#Chr2 Chr2' \
        'NIP-T2T#1#Chr3 Chr3' \
        'NIP-T2T#1#Chr4 Chr4' \
        'NIP-T2T#1#Chr5 Chr5' \
        'NIP-T2T#1#Chr6 Chr6' \
        'NIP-T2T#1#Chr7 Chr7' \
        'NIP-T2T#1#Chr8 Chr8' \
        'NIP-T2T#1#Chr9 Chr9' \
        'NIP-T2T#1#Chr10 Chr10' \
        'NIP-T2T#1#Chr11 Chr11' \
        'NIP-T2T#1#Chr12 Chr12' > "$MAP"
fi

if [ ! -s "$OUT/classification_policy.txt" ]; then
    printf '%s\n' \
        '统一 VCF 分类规则（2026-09-28）' \
        '1. NIP-T2T#1#ChrN 重命名为 ChrN。' \
        '2. 只保留 FILTER=PASS。' \
        '3. 只保留完整非参考基因型：GT=alt 且 GT 不是 missing；因此排除 0/0、./.、0/.、./1 等。' \
        '4. 使用 NIP-T2T_Chr.fa 做 REF 一致性检查和标准化，并用 bcftools norm -m -any --multi-overlaps 0 拆分多等位记录。' \
        '5. SNP：REF、ALT 均为 A/C/G/T/N，长度均为 1 且 REF!=ALT。' \
        '6. small indel：REF、ALT 均为 A/C/G/T/N，长度不同，长度差小于 50 bp。' \
        '7. SV：REF、ALT 均为 A/C/G/T/N，长度差大于等于 50 bp。' \
        '8. MNP、symbolic ALT、BND 及其他不符合上述条件的记录写入 other，不参加 SNP/indel/SV 主评估。' \
        '9. minigraph 当前benchmark只评估SV；脚本仍输出完整分类，便于审计。' > "$OUT/classification_policy.txt"
fi

count_records() {
    if "$BCF" index -n "$1" >/dev/null 2>&1; then
        "$BCF" index -n "$1"
    else
        "$BCF" view -H "$1" | wc -l | tr -d ' '
    fi
}

index_vcf() {
    "$TABIX" -f -p vcf "$1"
}

is_selected() {
    value=$1
    selected=$2
    [ -z "$selected" ] || [ "$value" = "$selected" ]
}

write_class() {
    input=$1
    output=$2
    expr=$3
    if [ "$FORCE" = 1 ] || [ ! -s "$output" ] || [ ! -s "$output.tbi" ]; then
        tmp=$output.tmp.$$
        rm -f "$tmp" "$tmp.tbi"
        "$BCF" view --threads "$THREADS" -i "$expr" -Oz -o "$tmp" "$input"
        index_vcf "$tmp"
        mv -f "$tmp" "$output"
        mv -f "$tmp.tbi" "$output.tbi"
    fi
}

prepare_one() {
    tool=$1
    sample=$2
    depth=$3
    graph=$4
    input=$5
    dir=$OUT/$tool/$sample/$depth/$graph
    mkdir -p "$dir/tmp"

    all=$dir/$sample.$graph.$tool.PASS.nonref.split.sorted.vcf.gz
    snp=$dir/$sample.$graph.$tool.PASS.snp.vcf.gz
    indel=$dir/$sample.$graph.$tool.PASS.small_indel.vcf.gz
    sv=$dir/$sample.$graph.$tool.PASS.sv.vcf.gz
    other=$dir/$sample.$graph.$tool.PASS.other.vcf.gz

    if [ ! -s "$all" ] || [ ! -s "$all.tbi" ]; then
        log "$tool $sample $depth $graph: normalize and split"
        renamed=$dir/renamed.$$
        normalized=$dir/normalized.$$
        sorted=$dir/sorted.$$
        rm -f "$renamed" "$normalized" "$sorted" "$sorted.tbi"

        "$BCF" annotate --rename-chrs "$MAP" --threads "$THREADS" -Oz -o "$renamed" "$input"
        "$BCF" view --threads "$THREADS" -f PASS -i 'GT="alt" && GT!="mis"' -Ou "$renamed" |
            "$BCF" norm --threads "$THREADS" -f "$REF" -c e -m -any --multi-overlaps 0 -Ou |
            "$BCF" view --threads "$THREADS" -i 'GT="alt" && GT!="mis"' -Ou |
            "$BCF" sort -m "$SORT_MEM" -T "$dir/tmp" -Oz -o "$sorted"
        index_vcf "$sorted"
        mv -f "$sorted" "$all"
        mv -f "$sorted.tbi" "$all.tbi"
        rm -f "$renamed" "$normalized"
    else
        log "$tool $sample $depth $graph: reuse normalized VCF"
    fi

    # bcftools uses a POSIX-compatible regex here; spell one-or-more as X X*.
    dna='REF~"^[ACGTN][ACGTN]*$" && ALT~"^[ACGTN][ACGTN]*$"'
    snp_expr="$dna && strlen(REF)=1 && strlen(ALT)=1 && REF!=ALT"
    indel_expr="$dna && strlen(REF)!=strlen(ALT) && abs(strlen(REF)-strlen(ALT))<50"
    sv_expr="$dna && abs(strlen(REF)-strlen(ALT))>=50"
    other_expr="$snp_expr || $indel_expr || $sv_expr"

    write_class "$all" "$snp" "$snp_expr"
    write_class "$all" "$indel" "$indel_expr"
    write_class "$all" "$sv" "$sv_expr"
    if [ "$FORCE" = 1 ] || [ ! -s "$other" ] || [ ! -s "$other.tbi" ]; then
        tmp=$other.tmp.$$
        rm -f "$tmp" "$tmp.tbi"
        "$BCF" view --threads "$THREADS" -e "$other_expr" -Oz -o "$tmp" "$all"
        index_vcf "$tmp"
        mv -f "$tmp" "$other"
        mv -f "$tmp.tbi" "$other.tbi"
    fi

    input_records=$(count_records "$input")
    all_records=$(count_records "$all")
    snp_records=$(count_records "$snp")
    indel_records=$(count_records "$indel")
    sv_records=$(count_records "$sv")
    other_records=$(count_records "$other")
    total=$((snp_records + indel_records + sv_records + other_records))
    [ "$total" -eq "$all_records" ] || die "$tool $sample $depth $graph: class counts do not partition all records ($snp_records+$indel_records+$sv_records+$other_records != $all_records)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$tool" "$sample" "$depth" "$graph" "$input" "$input_records" "$all_records" \
        "$snp_records" "$indel_records" "$sv_records" "$other_records" "$dir" >> "$SUMMARY"
    rm -rf "$dir/tmp"
    log "$tool $sample $depth $graph: all=$all_records SNP=$snp_records small_indel=$indel_records SV=$sv_records other=$other_records"
}

source_for() {
    tool=$1
    sample=$2
    depth=$3
    graph=$4
    prefix=$sample"_"$depth
    case "$tool" in
        vg)
            printf '%s\n' "$ROOT/results/$sample/$depth/$graph/vg/$sample.$graph.vg.vcf.gz"
            ;;
        pangenie)
            if [ "$graph" = minigraph ]; then
                printf '%s\n' "$ROOT/results2/minigraph/02pangenie/$sample/$depth/$sample.$depth.minigraph.pangenie.sv.vcf.gz"
            else
                printf '%s\n' "$ROOT/results2/$graph/02pangenie/$sample/$depth/$sample.$depth.$graph.pangenie.vcf.gz"
            fi
            ;;
        bayestyper)
            bayes_graph=$graph
            [ "$graph" = minigraph ] && bayes_graph=MINIGRAPH
            printf '%s\n' "$ROOT/bayestyper_jiqun/$bayes_graph/04bayestyper/$sample/$depth/$prefix.$bayes_graph.bayestyper.vcf.gz"
            ;;
    esac
}

printf '%s\n' 'tool sample depth graph input input_records prepared_records SNP small_indel SV other prepared_dir' > "$SUMMARY"
printf '%s\n' 'tool sample depth graph input' > "$MISSING"
found=0
for tool in vg pangenie bayestyper; do
    is_selected "$tool" "$ONLY_TOOL" || continue
    for sample in 9311 IR64 W0169 W1735 W3033 zhonghua11; do
        is_selected "$sample" "$ONLY_SAMPLE" || continue
        for depth in 5x 10x 15x 20x; do
            is_selected "$depth" "$ONLY_DEPTH" || continue
            for graph in MC PGGB minigraph; do
                is_selected "$graph" "$ONLY_GRAPH" || continue
                input=$(source_for "$tool" "$sample" "$depth" "$graph")
                if [ ! -s "$input" ]; then
                    printf '%s\t%s\t%s\t%s\t%s\n' "$tool" "$sample" "$depth" "$graph" "$input" >> "$MISSING"
                    log "WARN missing input: $input"
                    continue
                fi
                found=$((found + 1))
                prepare_one "$tool" "$sample" "$depth" "$graph" "$input"
            done
        done
    done
done

[ "$found" -gt 0 ] || die 'no input VCF was found'
log "DONE: processed $found input VCFs"
log "summary: $SUMMARY"
log "missing inputs: $MISSING"
