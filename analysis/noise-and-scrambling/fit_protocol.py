"""Fit the frozen schedule unitary to F(t2) exp(i mu V) F(t1), F(t) = e^{-iHt} (x) e^{-iHt}-variant,
using the same SYK H and coupling V that REPRODUCE_FULL_CASE.source(7001) builds.
Objective: distance between right-side difference operators produced from the source."""
import itertools
import json
import math
import sys

import numpy as np
from scipy.optimize import minimize

I = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1., -1.]).astype(complex)


def kron(items):
    out = np.ones((1, 1), complex)
    for it in items:
        out = np.kron(out, it)
    return out


gamma = [kron([Z] * q + [p] + [I] * (3 - q)) for q in range(4) for p in (X, Y)]
parity = kron([Z] * 4)
draws = np.random.default_rng(7001).normal(0, math.sqrt(6 / 8 ** 3), 70)
H = np.zeros((16, 16), complex)
for idx, j in zip(itertools.combinations(range(8), 4), draws):
    g = np.eye(16, dtype=complex)
    for q in idx:
        g = g @ gamma[q]
    H -= j * g / 4
V = sum(.5j * np.kron(g @ parity, g) for g in gamma)
V = (V + V.conj().T) / 2
eH, wH = np.linalg.eigh(H)
eV, wV = np.linalg.eigh(V)

U = np.load('U_sched.npy')
import noise_sweep as ns_mod
ns_mod.init(sys.argv[1])
src = ns_mod.NS['_initial']


def right_of(W):
    full = W @ src @ W.conj().T
    return np.trace(full.reshape(16, 16, 16, 16), axis1=0, axis2=2)


target = right_of(U)


def ev(t, conj_right):
    u = (wH * np.exp(-1j * eH * t)) @ wH.conj().T
    return np.kron(u, u.conj() if conj_right else u)


def cand(p, variant):
    t1, mu, t2 = p
    E = (wV * np.exp(1j * mu * eV)) @ wV.conj().T
    return ev(t2, variant) @ E @ ev(t1, variant)


def cost(p, variant):
    return float(np.abs(right_of(cand(p, variant)) - target).sum())


best = []
for variant in (False, True):
    grid = np.linspace(-3, 3, 13)
    for t1, mu, t2 in itertools.product(grid, grid, grid):
        best.append((cost((t1, mu, t2), variant), (t1, mu, t2), variant))
best.sort(key=lambda x: x[0])
results = []
for c0, p0, variant in best[:12]:
    r = minimize(cost, p0, args=(variant,), method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-12, maxiter=4000))
    W = cand(r.x, variant)
    fid = abs(np.trace(W.conj().T @ U)) / 256
    results.append(dict(cost=r.fun, t1=r.x[0], mu=r.x[1], t2=r.x[2], conj_right=variant, full_unitary_fidelity=fid))
results.sort(key=lambda d: d['cost'])
print('target L1 norm', float(np.abs(target).sum()))
for d in results[:5]:
    print(json.dumps(d))
