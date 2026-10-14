#!/usr/bin/env bash
set -Eeuo pipefail

# Rebuild the 36-genome MC graph with Minigraph-Cactus.
# The historical Cactus version/options were not retained, so this is a
# documented reconstruction workflow, not a claim of byte-identical recovery.

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SEQFILE=${SEQFILE:-}
REFERENCE=${REFERENCE:-NIP-T2T}
OUTDIR=${OUTDIR:-$PWD/rice36_mc_rebuild}
JOBSTORE=${JOBSTORE:-$OUTDIR/jobstore}
OUTNAME=${OUTNAME:-rice36_mc}
THREADS=${THREADS:-128}
RESTART=${RESTART:-0}
CACTUS_PANGENOME=${CACTUS_PANGENOME:-cactus-pangenome}

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }
trap 'rc=$?; log "FAILED line=$LINENO command=$BASH_COMMAND exit=$rc"; exit "$rc"' ERR

[[ -n "$SEQFILE" ]] || die "set SEQFILE to a two-column Cactus seqfile (sample<TAB>FASTA path)"
[[ -s "$SEQFILE" ]] || die "seqfile missing or empty: $SEQFILE"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
[[ "$RESTART" == 0 || "$RESTART" == 1 ]] || die "RESTART must be 0 or 1"
command -v "$CACTUS_PANGENOME" >/dev/null 2>&1 || die "cannot find cactus-pangenome: $CACTUS_PANGENOME"

# Validate the seqfile while permitting the historical optional _MINIGRAPH_
# pseudo-sequence row. URLs other than local file paths are left to Cactus.
awk -v ref="$REFERENCE" '
  /^[[:space:]]*($|#)/ { next }
  NF < 2 { print "Malformed seqfile row " NR ": expected sample and FASTA path" > "/dev/stderr"; bad=1; next }
  $1 == ref { have_ref=1 }
  $1 == "_MINIGRAPH_" { minigraph++; next }
  { samples++; names[$1]++ }
  END {
    for (name in names) if (names[name] > 1) {
      print "Duplicate sample name: " name > "/dev/stderr"; bad=1
    }
    if (samples != 36) { print "Expected 36 genome rows (excluding _MINIGRAPH_), found " samples > "/dev/stderr"; bad=1 }
    if (!have_ref) { print "Reference sample not found in seqfile: " ref > "/dev/stderr"; bad=1 }
    if (minigraph > 1) { print "At most one _MINIGRAPH_ row is allowed" > "/dev/stderr"; bad=1 }
    exit bad
  }
' "$SEQFILE" || die "seqfile validation failed"

while read -r sample fasta rest; do
  [[ -n "${sample:-}" && "$sample" != \#* ]] || continue
  [[ -n "${fasta:-}" ]] || die "missing FASTA path for $sample"
  case "$fasta" in
    file://*) local_fasta=${fasta#file://}; [[ -s "$local_fasta" ]] || die "missing FASTA: $local_fasta" ;;
    /*) [[ -s "$fasta" ]] || die "missing FASTA: $fasta" ;;
    http://*|https://*) ;;
    *) [[ -s "$fasta" ]] || die "relative FASTA path not found from current directory: $fasta" ;;
  esac
done < "$SEQFILE"

if [[ -e "$JOBSTORE" && "$RESTART" != 1 ]]; then
  die "job store already exists: $JOBSTORE (set RESTART=1 to resume; it will not be deleted)"
fi
if [[ "$RESTART" == 1 && ! -e "$JOBSTORE" ]]; then
  die "RESTART=1 requested but job store does not exist: $JOBSTORE"
fi
if [[ -d "$OUTDIR" && "$RESTART" != 1 ]] && find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  die "output directory is not empty: $OUTDIR (choose a new OUTDIR or resume the existing job)"
fi

mkdir -p "$OUTDIR" "$(dirname "$JOBSTORE")"
LOG="$OUTDIR/${OUTNAME}.build.log"
exec > >(tee -a "$LOG") 2>&1

cmd=("$CACTUS_PANGENOME" "$JOBSTORE" "$SEQFILE"
  --outDir "$OUTDIR" --outName "$OUTNAME" --reference "$REFERENCE"
  --gfa --vcf --giraffe --gbz --batchSystem single_machine --maxCores "$THREADS")
if [[ "$RESTART" == 1 ]]; then cmd+=(--restart); fi

log "START Minigraph-Cactus graph build"
log "seqfile=$SEQFILE reference=$REFERENCE outdir=$OUTDIR jobstore=$JOBSTORE threads=$THREADS"
log "cactus-pangenome=$($CACTUS_PANGENOME --version 2>&1 | head -n 1 || true)"
log "command: ${cmd[*]}"
"${cmd[@]}"

mapfile -t GFA_FILES < <(find "$OUTDIR" -maxdepth 1 -type f \( -name '*.gfa' -o -name '*.gfa.gz' \) -print | sort)
mapfile -t VCF_FILES < <(find "$OUTDIR" -maxdepth 1 -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) -print | sort)
(( ${#GFA_FILES[@]} > 0 )) || die "Cactus completed without a GFA in $OUTDIR"
(( ${#VCF_FILES[@]} > 0 )) || die "Cactus completed without a VCF in $OUTDIR"

{
  printf 'artifact\tpath\tbytes\n'
  for f in "${GFA_FILES[@]}" "${VCF_FILES[@]}"; do
    printf '%s\t%s\t%s\n' "$(basename "$f")" "$f" "$(stat -c '%s' "$f")"
  done
} > "$OUTDIR/${OUTNAME}.artifacts.tsv"
log "SUCCESS: found ${#GFA_FILES[@]} GFA and ${#VCF_FILES[@]} VCF outputs"
log "artifact_manifest=$OUTDIR/${OUTNAME}.artifacts.tsv"
