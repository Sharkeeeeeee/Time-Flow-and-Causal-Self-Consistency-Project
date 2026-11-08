"""Scrambling advantage: best SYK success over (t1, t2) at fixed coupling mu, minus the H=0 control
(pure coupling = partial swap, P_ctrl = 1/2 + 1/2 sin^2 mu). Both signs of mu."""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import minimize

sys.argv = sys.argv[:2]
import fit_protocol as fp

eH, wH = np.linalg.eigh(fp.H)


def helstrom(d):
    return 0.5 + 0.25 * float(np.abs(np.linalg.eigvalsh((d + d.conj().T) / 2)).sum())


def F(t):
    u = (wH * np.exp(-1j * eH * t)) @ wH.conj().T
    return np.kron(u, u)


def P(t1, mu, t2, Hzero=False):
    E = (fp.wV * np.exp(1j * mu * fp.eV)) @ fp.wV.conj().T
    W = E if Hzero else F(t2) @ E @ F(t1)
    return helstrom(fp.right_of(W))


grid = np.linspace(-4, 4, 17)
rows = []
for mu in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.76, 0.9, 1.1, 1.3]:
    rec = {'mu': mu, 'control': P(0, mu, 0, True)}
    for sign in (+1, -1):
        m = sign * mu
        pts = sorted(((P(a, m, b), (a, b)) for a, b in itertools.product(grid, grid)), key=lambda x: -x[0])
        best = max((minimize(lambda p: -P(p[0], m, p[1]), p0, method='Nelder-Mead', options=dict(maxiter=800)) for _, p0 in pts[:4]),
                   key=lambda r: -r.fun)
        rec['syk_+mu' if sign > 0 else 'syk_-mu'] = -best.fun
        rec['times_+mu' if sign > 0 else 'times_-mu'] = [round(float(x), 3) for x in best.x]
    rows.append(rec)
    print(json.dumps(rec), flush=True)
