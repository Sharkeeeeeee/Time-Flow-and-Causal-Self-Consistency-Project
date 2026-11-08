"""Noiseless Helstrom ceiling of the right-side output over the protocol family
F(t2) exp(i mu V) F(t1), F(t)=e^{-iHt} (x) e^{-iHt}, applied to the repo's seed-7001 source.
Compares the SYK H against an H=0 control (coupling only, no scrambling)."""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import minimize

sys.argv = sys.argv[:2]
import fit_protocol as fp  # reuses H, V, source, frozen-schedule target

src = fp.src


def helstrom(d):
    return 0.5 + 0.25 * float(np.abs(np.linalg.eigvalsh((d + d.conj().T) / 2)).sum())


def make(Hm):
    e, w = np.linalg.eigh(Hm)
    def F(t):
        u = (w * np.exp(-1j * e * t)) @ w.conj().T
        return np.kron(u, u)
    def P(p):
        t1, mu, t2 = p
        E = (fp.wV * np.exp(1j * mu * fp.eV)) @ fp.wV.conj().T
        W = F(t2) @ E @ F(t1)
        return helstrom(fp.right_of(W))
    return P


out = {'frozen_schedule_right_helstrom': helstrom(fp.target)}
for name, Hm in (('SYK H', fp.H), ('H=0 control', np.zeros((16, 16), complex))):
    P = make(Hm)
    grid = np.linspace(-4, 4, 17)
    pts = sorted(((P(p), p) for p in itertools.product(grid, np.linspace(-np.pi, np.pi, 25), grid)), key=lambda x: -x[0])
    best = max((minimize(lambda p: -P(p), p0, method='Nelder-Mead', options=dict(maxiter=3000)) for _, p0 in pts[:8]),
               key=lambda r: -r.fun)
    t1, mu, t2 = best.x
    mu_scan = {f'{m:.2f}': P((t1, m, t2)) for m in np.linspace(0, np.pi, 9)}
    out[name] = dict(max_helstrom=-best.fun, t1=t1, mu=mu, t2=t2, mu_scan_at_best_times=mu_scan)
print(json.dumps(out, indent=1))
