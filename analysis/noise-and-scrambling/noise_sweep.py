"""Noise-factor sweep for the frozen seed-7001 case of plan-a/REPRODUCE_FULL_CASE.py.

The original script hard-codes a noise improvement factor of 64 in two places:
  * relaxation: T1 = 400 us * R, T2 = 200 us * R
  * gate depolarizing error: p = (3e-4 for sx, 3e-3 for cz) / R
This script exposes them as T (coherence factor) and G (gate-error factor), keeps
the frozen 229.46 us schedule and the frozen receiver policies, and re-evaluates
the expected success probability. Readout assignment error (1-2%) is NOT scaled.

Usage: python -I noise_sweep.py <path to plan-a/REPRODUCE_FULL_CASE.py> [--quick] [--out FILE]
"""
import argparse
import csv
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

LABELS = ['oracle', 'n2048_d0_nominal', 'n2048_d1_nominal', 'n2048_d2_nominal']
NS = None


def build(script_path):
    src = Path(script_path).read_text(encoding='utf-8')
    reps = [
        ("def relaxation(rho,us,n,resource=64):",
         "def relaxation(rho,us,n,resource=None):\n    if resource is None:resource=CFG['T']"),
        ("p=(.0003 if g['name']=='sx' else .003)/64", "p=(.0003 if g['name']=='sx' else .003)/CFG['G']"),
        ("print('ALTERNATE full7001/64 completed_us',checkpoint/1000,flush=True)", "pass"),
    ]
    for old, new in reps:
        assert old in src, old
        src = src.replace(old, new)
    assert '/64' not in src.split('def transport')[1]
    ns = {'__file__': str(Path(script_path).resolve()), '__name__': 'sweep_ns', 'CFG': {'T': 64, 'G': 64}}
    exec(compile(src, script_path, 'exec'), ns)
    return ns


def init(script_path):
    global NS
    import base64, gzip
    NS = build(script_path)
    NS['_initial'], _ = NS['source'](7001)
    NS['_timeline'] = json.loads(gzip.decompress(base64.b64decode((NS['DATA'] / 'schedule_seed_7001.json.gz.b64').read_text())))
    NS['_policies'] = {l: json.loads((NS['DATA'] / f'policy_{l}_seed_7001.json').read_text()) for l in LABELS}


def run(cfg):
    series, T, G = cfg
    np = NS['np']
    NS['CFG']['T'], NS['CFG']['G'] = T, G
    t0 = time.perf_counter()
    full, _ = NS['transport'](NS['_initial'], NS['_timeline'])
    right = np.trace(full.reshape(16, 16, 16, 16), axis1=0, axis2=2)
    scores = {l: NS['receiver'](right, NS['_policies'][l])[0] for l in LABELS}
    return dict(series=series, T=T, G=G, **scores, seconds=time.perf_counter() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('script')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument('--out', default='noise_sweep_results.csv')
    a = ap.parse_args()
    inf = float('inf')
    if a.quick:
        cfgs = [('joint', 64, 64)]
    else:
        R = [1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128, 192, 256, 512, inf]
        cfgs = [('joint', r, r) for r in R]
        cfgs += [('coherence only', r, 64) for r in R if r != 64]
        cfgs += [('gates only', 64, r) for r in R if r != 64]
    # Resume: results.jsonl next to --out holds one finished configuration per line.
    jsonl = Path(a.out).with_name('results.jsonl')
    done = {}
    if jsonl.exists():
        for line in jsonl.read_text().splitlines():
            r = json.loads(line)
            done[(r['series'], float(r['T']), float(r['G']))] = r
    todo = [c for c in cfgs if (c[0], float(c[1]), float(c[2])) not in done]
    print(f'{len(cfgs) - len(todo)} done, {len(todo)} to run', flush=True)
    if todo:
        with ProcessPoolExecutor(max_workers=min(a.workers, len(todo)), initializer=init, initargs=(a.script,)) as ex:
            for r in ex.map(run, todo):
                with jsonl.open('a') as f:
                    f.write(json.dumps(r) + '\n')
                done[(r['series'], float(r['T']), float(r['G']))] = r
                print(json.dumps(r), flush=True)
    rows = [done[(c[0], float(c[1]), float(c[2]))] for c in cfgs]
    ref = json.loads((Path(a.script).parent / 'data' / 'results64_seed_7001.json').read_text())
    for r in rows:
        if r['T'] == 64 and r['G'] == 64:
            for l in LABELS:
                want = next(x for x in ref if x['label'] == l)['success']
                assert abs(r[l] - want) < 1e-9, (l, r[l], want)
            print('R=64 reproduces archived endpoints to <1e-9', flush=True)
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
