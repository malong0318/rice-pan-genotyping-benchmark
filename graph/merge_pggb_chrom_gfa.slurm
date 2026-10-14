#!/usr/bin/env bash
# Merge 12 chromosome-level PGGB GFA 1.0 graphs with disjoint numeric IDs.
# Submit: sbatch /data/malong/benchmark/code/merge_pggb_chrom_gfa.slurm

#SBATCH --job-name=merge_pggb_gfa
#SBATCH --partition=AMD_9654
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/malong/benchmark/merge_pggb_gfa_%A.out
#SBATCH --error=/data/malong/benchmark/merge_pggb_gfa_%A.err

set -Eeuo pipefail
BASE=/data/malong/benchmark
ROOT="$BASE/graph/pggb/pggb_results"
OUTDIR="$BASE/graph/pggb/merged_12chr"
OUT="$OUTDIR/rice36.pggb.12chr.merged.gfa"
MANIFEST="$OUTDIR/rice36.pggb.12chr.merged.manifest.tsv"
THREADS=${SLURM_CPUS_PER_TASK:-32}
ODGI=/data/malong/miniconda3/envs/benchmark/bin/odgi
PYTHON=/data/malong/miniconda3/envs/benchmark/bin/python

[[ -x "$ODGI" ]] || { echo "Missing odgi: $ODGI" >&2; exit 127; }
[[ -x "$PYTHON" ]] || { echo "Missing Python: $PYTHON" >&2; exit 127; }

mkdir -p "$OUTDIR"
exec > >(tee "$OUTDIR/merge_pggb_gfa.log") 2>&1

declare -a GFA_FILES
for chr in $(seq 1 12); do
  mapfile -t matches < <(find "$ROOT" -mindepth 2 -maxdepth 2 -type f \
    -name "Chr${chr}.fa.gz.*.smooth.final.gfa" -print | sort)
  if [[ "${#matches[@]}" -ne 1 ]]; then
    echo "Expected exactly one GFA for Chr${chr}, found ${#matches[@]}" >&2
    printf '%s\n' "${matches[@]:-<none>}" >&2
    exit 2
  fi
  [[ -s "${matches[0]}" ]] || { echo "Empty GFA: ${matches[0]}" >&2; exit 2; }
  GFA_FILES[$chr]="${matches[0]}"
  echo -e "Chr${chr}\t${matches[0]}"
done

TMP="$OUTDIR/.rice36.pggb.12chr.merged.gfa.tmp"
rm -f "$TMP"
"$PYTHON" - "$TMP" "$MANIFEST" "$ROOT" "${GFA_FILES[@]:1}" <<'PY'
import re
import sys
from pathlib import Path

tmp = Path(sys.argv[1])
manifest = Path(sys.argv[2])
root = sys.argv[3]
files = sys.argv[4:]
if len(files) != 12:
    raise SystemExit(f"expected 12 input GFAs, got {len(files)}")
num = re.compile(r"^[0-9]+$")
ptok = re.compile(r"^([0-9]+)([+-])$")
wtok = re.compile(r"([<>])([0-9]+)")
max_ids, counts = [], []
for i, fn in enumerate(files, 1):
    mx = 0
    c = {k: 0 for k in "HSLPJW"}
    with open(fn, encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] in c: c[f[0]] += 1
            if f[0] == "S":
                if len(f) < 3 or not num.fullmatch(f[1]):
                    raise SystemExit(f"Chr{i}: invalid S node ID in {fn}: {line[:200]}")
                mx = max(mx, int(f[1]))
    if mx == 0: raise SystemExit(f"Chr{i}: no S records in {fn}")
    max_ids.append(mx); counts.append(c)
offsets, next_offset = [], 0
for mx in max_ids:
    offsets.append(next_offset); next_offset += mx
def nid(raw, off, fn):
    if not num.fullmatch(raw): raise SystemExit(f"invalid node reference {raw!r} in {fn}")
    return str(int(raw) + off)
with open(tmp, "w", encoding="utf-8") as out:
    out.write("H\tVN:Z:1.0\tZZ:Z:merged_from_12_PGGB_chromosome_GFAs\n")
    for fn, off in zip(files, offsets):
        with open(fn, encoding="utf-8") as fh:
            for line in fh:
                if not line or line[0] == "#": continue
                f = line.rstrip("\n").split("\t")
                tag = f[0]
                if tag == "H": continue
                if tag == "S":
                    f[1] = nid(f[1], off, fn)
                elif tag == "L":
                    if len(f) < 5: raise SystemExit(f"malformed L line in {fn}")
                    f[1] = nid(f[1], off, fn); f[3] = nid(f[3], off, fn)
                elif tag == "P":
                    if len(f) < 3: raise SystemExit(f"malformed P line in {fn}")
                    vals = []
                    for token in f[2].split(","):
                        if token == "*": vals.append(token); continue
                        m = ptok.fullmatch(token)
                        if not m: raise SystemExit(f"malformed P token {token!r} in {fn}")
                        vals.append(nid(m.group(1), off, fn) + m.group(2))
                    f[2] = ",".join(vals)
                elif tag == "W" and len(f) >= 7:
                    f[6] = wtok.sub(lambda m: m.group(1) + nid(m.group(2), off, fn), f[6])
                out.write("\t".join(f) + "\n")
with open(manifest, "w", encoding="utf-8") as mf:
    mf.write("chromosome\tinput_gfa\tmax_input_node_id\tnode_id_offset\n")
    for i, (fn, mx, off) in enumerate(zip(files, max_ids, offsets), 1):
        mf.write(f"Chr{i}\t{fn}\t{mx}\t{off}\n")
    mf.write(f"total_max_node_id\t{root}\t{next_offset}\t0\n")
    mf.write("method\tdisjoint cumulative numeric offsets; S/L/P/W references rewritten\n")
PY
mv "$TMP" "$OUT"
echo "[$(date -Is)] merged GFA: $OUT"
echo "[$(date -Is)] size: $(stat -c%s "$OUT") bytes"

if [[ -x "$ODGI" ]]; then
  "$ODGI" validate -i "$OUT" -t "$THREADS" > "$OUTDIR/odgi.validate.txt" 2>&1
  "$ODGI" build -g "$OUT" -o "$OUTDIR/rice36.pggb.12chr.merged.og" -t "$THREADS"
else
  echo "Missing odgi validator: $ODGI" >&2
  exit 127
fi

HCOUNT=$(awk '$1=="H"{n++} END{print n+0}' "$OUT")
[[ "$HCOUNT" == 1 ]] || { echo "Expected one H line, found $HCOUNT" >&2; exit 1; }
if awk -F '\t' '$1=="S"{print $2}' "$OUT" | sort -n | uniq -d | head -1 | grep -q .; then
  echo "Duplicate merged node ID detected" >&2
  exit 1
fi
echo "[$(date -Is)] duplicate-node check completed"
echo "[$(date -Is)] completed"
