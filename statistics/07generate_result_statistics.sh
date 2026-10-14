#!/usr/bin/env bash
set -Eeuo pipefail

BASE=/data/malong/benchmark
PAPER=$BASE/paper
OUT=$PAPER/result_statistics.tsv
AGG=$PAPER/result_statistics_aggregate.tsv
PY=$BASE/../miniconda3/bin/python
[[ -x "$PY" ]] || PY=python3

"$PY" - "$BASE" "$OUT" "$AGG" <<'PY'
import csv, json, pathlib, sys

base = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])
agg = pathlib.Path(sys.argv[3])
samples = ['9311', 'IR64', 'W0169', 'W1735', 'W3033', 'zhonghua11']
depths = ['5x', '10x', '15x', '20x']
tools = ['vg', 'pangenie', 'bayestyper', 'graphtyper']
classes = ['snp', 'indel', 'sv']
rows = []

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def text_or_na(x): return 'NA' if x is None else str(x)

def read_count(summary_tsv, sample, depth, graph):
    try:
        with summary_tsv.open() as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r.get('sample') == sample and r.get('depth') == depth and r.get('graph') == graph:
                    truth = r.get('truth_records') or r.get('truth_records_excluding_BND') or 'NA'
                    return truth, r.get('call_records') or 'NA'
    except FileNotFoundError: pass
    return 'NA', 'NA'

def parse_rtg(path):
    try: lines = path.read_text(errors='replace').splitlines()
    except FileNotFoundError: return None
    for line in lines:
        if not line.strip() or line.startswith('-') or line.lstrip().startswith('Threshold'): continue
        fields = line.split()
        if len(fields) >= 8:
            return {'tp_baseline': fields[1], 'tp_call': fields[2], 'fp': fields[3], 'fn': fields[4], 'precision': fields[5], 'recall': fields[6], 'f1': fields[7], 'evaluator': 'RTG-vcfeval'}
    return None

def parse_truvari(path):
    try: x = json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError): return None
    return {'tp_baseline': x.get('TP-base', 'NA'), 'tp_call': x.get('TP-comp', 'NA'), 'fp': x.get('FP', 'NA'), 'fn': x.get('FN', 'NA'), 'precision': x.get('precision', 'NA'), 'recall': x.get('recall', 'NA'), 'f1': x.get('f1', 'NA'), 'gt_concordance': x.get('gt_concordance', 'NA'), 'evaluator': 'Truvari'}

def scope(tool, graph, klass):
    if tool == 'vg' and graph == 'PGGB': return 'NA', 'VG-PGGB VG index construction failed'
    if graph == 'minigraph' and klass != 'sv': return 'NA', 'minigraph benchmark scope is SV only'
    return 'expected', ''

for tool in tools:
    for klass in classes:
        root = base / 'evaluation' / klass / tool
        summary_tsv = root / 'summary.tsv'
        for sample in samples:
            for depth in depths:
                for graph in ['MC', 'PGGB', 'minigraph']:
                    status, reason = scope(tool, graph, klass)
                    result = None if status == 'NA' else root / graph / sample / depth
                    truth, call = read_count(summary_tsv, sample, depth, graph)
                    rec = {'tool': tool, 'graph': graph, 'sample': sample, 'depth': depth, 'variant_class': klass, 'truth_records': truth, 'call_records': call, 'tp_baseline': 'NA', 'tp_call': 'NA', 'fp': 'NA', 'fn': 'NA', 'precision': 'NA', 'recall': 'NA', 'f1': 'NA', 'gt_concordance': 'NA', 'evaluator': 'NA', 'status': status, 'reason': reason, 'result_directory': str(result) if result else 'NA'}
                    if status != 'NA':
                        metric = parse_truvari(result / 'summary.json') if klass == 'sv' else parse_rtg(result / 'summary.txt')
                        if metric:
                            rec.update({k: text_or_na(v) for k, v in metric.items()}); rec['status'] = 'complete'; rec['reason'] = ''
                        else:
                            rec['status'] = 'missing_or_incomplete'; rec['reason'] = 'summary file missing or could not be parsed'
                    rows.append(rec)

fields = ['tool','graph','sample','depth','variant_class','truth_records','call_records','tp_baseline','tp_call','fp','fn','precision','recall','f1','gt_concordance','evaluator','status','reason','result_directory']
with out.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)

agg_fields = ['tool','graph','variant_class','n_complete','n_expected','mean_precision','mean_recall','mean_f1','mean_tp_baseline','mean_tp_call','mean_fp','mean_fn','note']
with agg.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=agg_fields, delimiter='\t', lineterminator='\n'); w.writeheader()
    for tool, graph, klass in sorted({(r['tool'], r['graph'], r['variant_class']) for r in rows}):
        subset = [r for r in rows if (r['tool'], r['graph'], r['variant_class']) == (tool, graph, klass)]
        complete = [r for r in subset if r['status'] == 'complete']; expected = [r for r in subset if r['status'] != 'NA']
        def mean(field):
            vals = [num(r[field]) for r in complete if num(r[field]) is not None]
            return 'NA' if not vals else '{:.6f}'.format(sum(vals) / len(vals))
        note = 'NA means out of scope or VG-PGGB index failure; not zero.' if any(r['status'] == 'NA' for r in subset) else ('some expected outputs are missing or incomplete.' if len(complete) < len(expected) else '')
        w.writerow({'tool': tool, 'graph': graph, 'variant_class': klass, 'n_complete': len(complete), 'n_expected': len(expected), 'mean_precision': mean('precision'), 'mean_recall': mean('recall'), 'mean_f1': mean('f1'), 'mean_tp_baseline': mean('tp_baseline'), 'mean_tp_call': mean('tp_call'), 'mean_fp': mean('fp'), 'mean_fn': mean('fn'), 'note': note})

print('rows={}'.format(len(rows)))
print('complete={}'.format(sum(r['status'] == 'complete' for r in rows)))
print('na={}'.format(sum(r['status'] == 'NA' for r in rows)))
print('missing_or_incomplete={}'.format(sum(r['status'] == 'missing_or_incomplete' for r in rows)))
PY
printf '%s\n' "Wrote $OUT" "Wrote $AGG"

