"""Plot noise_sweep_results.csv -> noise_sweep.png (mean of the three finite-calibration receivers)."""
import csv
import math
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src)))
NOM = ['n2048_d0_nominal', 'n2048_d1_nominal', 'n2048_d2_nominal']
INF_X = 1024  # where the noiseless point is drawn

SERIES = [  # fixed slot order from the reference palette (first three validate all-pairs)
    ('joint', 'Both scaled together', '#2a78d6'),
    ('coherence only', 'Only T1/T2 scaled (gates fixed at 64x)', '#eb6834'),
    ('gates only', 'Only gate error scaled (T1/T2 fixed at 64x)', '#1baf7a'),
]
SURF, INK, INK2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8984', '#e6e5e1'


def factor(r, series):
    v = float(r['T'] if series != 'gates only' else r['G'])
    return INF_X if math.isinf(v) else v


fig, ax = plt.subplots(figsize=(9, 5.4), dpi=150)
fig.patch.set_facecolor(SURF)
ax.set_facecolor(SURF)
by = {}
for key, _, _ in SERIES:
    pts = sorted((factor(r, key), sum(float(r[k]) for k in NOM) / 3 * 100)
                 for r in rows if r['series'] == key or (r['T'] == '64' and r['G'] == '64' and r['series'] == 'joint'))
    by[key] = pts

for y, lab in ((70, '70% target'), (50, 'chance (50%)')):
    ax.axhline(y, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(1.05, y + 0.35, lab, color=INK2, fontsize=9, va='bottom')
ax.axvline(64, color=GRID, lw=1, zorder=0)
ax.text(64, 48.6, ' assumed 64x', color=INK2, fontsize=9, va='bottom')

for key, label, color in SERIES:
    xs, ys = zip(*by[key])
    ax.plot(xs, ys, color=color, lw=2, zorder=3, label=label, solid_capstyle='round')
    ax.scatter(xs, ys, s=22, color=color, edgecolor=SURF, linewidth=1.5, zorder=4)

ax.set_xscale('log', base=2)
ticks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, INF_X]
ax.set_xticks(ticks)
ax.set_xticklabels([str(t) for t in ticks[:-1]] + ['noiseless'])
ax.set_xlabel('Noise improvement factor R over the baseline (T1 = 400 us x R, gate error / R)', color=INK2, fontsize=10)
ax.set_ylabel('Expected success probability (%)', color=INK2, fontsize=10)
ax.set_ylim(48, max(max(y for _, y in v) for v in by.values()) + 2)
ax.grid(axis='y', color=GRID, lw=0.8)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
for s in ('left', 'bottom'):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ax.set_title('Seed-7001 frozen case: success vs. noise factor\n(mean of three finite-calibration receivers; schedule and receivers frozen)',
             color=INK, fontsize=11, loc='left')
leg = ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.12), frameon=False, fontsize=9)
for t in leg.get_texts():
    t.set_color(INK)
fig.tight_layout()
fig.savefig(out, facecolor=SURF)
print('wrote', out)
