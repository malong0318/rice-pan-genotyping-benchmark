#!/usr/bin/env bash
set -Eeuo pipefail

BASE=/data/malong/benchmark
ROOT=$BASE/results_performance
OUT=$BASE/paper/performance_statistics.tsv
mkdir -p "$BASE/paper"
PY=$BASE/../miniconda3/bin/python
[[ -x "$PY" ]] || PY=python3

"$PY" - "$ROOT" "$OUT" <<'PY'
import csv, pathlib, re, sys

root = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])
rows = []

def parse_kv(path):
    data = {}
    try:
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for r in reader:
                if r.get('metric') and r.get('value') is not None:
                    data[r['metric']] = r['value']
    except FileNotFoundError:
        pass
    return data

def field(text, pattern, default='NA'):
    m = re.search(pattern, text, re.MULTILINE)
    return m.group(1).strip() if m else default

def elapsed_seconds(value):
    if value == 'NA': return 'NA'
    parts = value.split(':')
    try:
        if len(parts) == 3: return '{:.3f}'.format(float(parts[0])*3600 + float(parts[1])*60 + float(parts[2]))
        if len(parts) == 2: return '{:.3f}'.format(float(parts[0])*60 + float(parts[1]))
        return '{:.3f}'.format(float(value))
    except ValueError:
        return 'NA'

def gb(value):
    try: return '{:.3f}'.format(float(value) / 1024 / 1024)
    except (TypeError, ValueError): return 'NA'

for raw in sorted(root.glob('*/*/9311/5x/performance.raw.txt')):
    rel = raw.relative_to(root).parts
    tool, graph, sample, depth = rel[:4]
    text = raw.read_text(errors='replace')
    meta = parse_kv(raw.parent / 'performance.tsv')
    elapsed = field(text, r'^\s*Elapsed.*\):\s*(\S+)$')
    user = field(text, r'^\s*User time \(seconds\):\s*(\S+)')
    system = field(text, r'^\s*System time \(seconds\):\s*(\S+)')
    cpu = field(text, r'^\s*Percent of CPU this job got:\s*(\S+)%')
    rss = field(text, r'^\s*Maximum resident set size \(kbytes\):\s*(\S+)')
    inputs = field(text, r'^\s*File system inputs:\s*(\S+)')
    outputs = field(text, r'^\s*File system outputs:\s*(\S+)')
    status = field(text, r'^\s*Exit status:\s*(\S+)')
    try: avg_cpu = '{:.3f}'.format(float(cpu) / 100)
    except ValueError: avg_cpu = 'NA'
    rows.append({
        'tool': tool, 'graph': graph, 'sample': sample, 'depth': depth,
        'configured_threads': meta.get('configured_threads', 'NA'),
        'status': 'complete' if status == '0' else 'failed',
        'elapsed': elapsed, 'elapsed_seconds': elapsed_seconds(elapsed),
        'user_seconds': user, 'system_seconds': system, 'cpu_percent': cpu,
        'average_cpu_cores': avg_cpu, 'max_rss_kb': rss, 'max_rss_gb': gb(rss),
        'filesystem_inputs': inputs, 'filesystem_outputs': outputs,
        'result_records': meta.get('result_records', 'NA'),
        'completed': meta.get('completed', 'NA'), 'raw_log': str(raw),
        'note': 'GNU time -v is authoritative; performance.tsv resource fields were recorded before raw timing finished.'
    })

# VG-PGGB is an unavailable benchmark combination, not a zero-valued result.
rows.append({
    'tool': 'vg', 'graph': 'PGGB', 'sample': '9311', 'depth': '5x',
    'configured_threads': '16', 'status': 'NA', 'elapsed': 'NA', 'elapsed_seconds': 'NA',
    'user_seconds': 'NA', 'system_seconds': 'NA', 'cpu_percent': 'NA',
    'average_cpu_cores': 'NA', 'max_rss_kb': 'NA', 'max_rss_gb': 'NA',
    'filesystem_inputs': 'NA', 'filesystem_outputs': 'NA', 'result_records': 'NA',
    'completed': 'NA', 'raw_log': 'NA',
    'note': 'VG-PGGB VG index construction failed; NA is not zero.'
})

fields = ['tool','graph','sample','depth','configured_threads','status','elapsed','elapsed_seconds','user_seconds','system_seconds','cpu_percent','average_cpu_cores','max_rss_kb','max_rss_gb','filesystem_inputs','filesystem_outputs','result_records','completed','raw_log','note']
rows.sort(key=lambda r: (r['tool'], r['graph']))
with out.open('w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
print('rows={}'.format(len(rows)))
print('complete={}'.format(sum(r['status'] == 'complete' for r in rows)))
print('na={}'.format(sum(r['status'] == 'NA' for r in rows)))
PY
printf '%s\n' "Wrote $OUT"

