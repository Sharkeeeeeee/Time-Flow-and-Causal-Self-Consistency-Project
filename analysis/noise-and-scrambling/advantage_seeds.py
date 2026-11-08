"""Scrambling advantage of the SYK protocol over non-scrambling controls, across disorder seeds.

Protocol (noiseless, mirrors REPRODUCE_FULL_CASE.source then a transport family):
  TFD of H at inverse temperature beta  ->  free(t_ins)  ->  measure first right qubit (both outcomes kept)
  -> free(t_ins) -> reset first left qubit to |b>  ->  F(t1) exp(i mu V) F(t2)-family transport
  -> Helstrom success for b on the 4 right qubits (best possible receiver).
Hamiltonians: 'syk4' (the repo's), 'syk2' (random quadratic, same Frobenius norm), 'zero' (H = 0).
The control is applied to the WHOLE protocol (source and transport), not only to the transport.

Usage:
  python advantage_seeds.py check  <REPRODUCE_FULL_CASE.py>   # source(7001, beta=4, t_ins=1) == repo source
  python advantage_seeds.py seeds  out.jsonl                   # 12 seeds x mu x sign x H
  python advantage_seeds.py tune   out.jsonl MU                # beta x t_ins scan at coupling MU, 12 seeds
"""
import itertools
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import minimize

SEEDS = [5003, 5009, 5011, 5021, 5023, 5039, 7001, 7013, 7019, 7027, 7039, 7043]
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1., -1.]).astype(complex)


def kron(items):
    out = np.ones((1, 1), complex)
    for it in items:
        out = np.kron(out, it)
    return out


GAMMA = [kron([Z] * q + [p] + [I2] * (3 - q)) for q in range(4) for p in (X, Y)]
PARITY = kron([Z] * 4)
V = sum(.5j * np.kron(g @ PARITY, g) for g in GAMMA)
V = (V + V.conj().T) / 2
EV, WV = np.linalg.eigh(V)
_, BASIS = np.linalg.eigh(sum(.5j * np.kron(g @ PARITY, g) for g in GAMMA))
IDX = np.arange(256); RBIT = (IDX >> 3) & 1; LBIT = (IDX >> 7) & 1


def hamiltonian(seed, kind):
    draws = np.random.default_rng(seed).normal(0, math.sqrt(6 / 8 ** 3), 70)
    H4 = np.zeros((16, 16), complex)
    for idx, j in zip(itertools.combinations(range(8), 4), draws):
        g = np.eye(16, dtype=complex)
        for q in idx:
            g = g @ GAMMA[q]
        H4 -= j * g / 4
    if kind == 'syk4':
        return H4
    if kind == 'zero':
        return np.zeros((16, 16), complex)
    J = np.random.default_rng(seed + 900000).normal(size=28)
    H2 = sum(1j * j * GAMMA[a] @ GAMMA[b] for (a, b), j in zip(itertools.combinations(range(8), 2), J))
    H2 = (H2 + H2.conj().T) / 2
    return H2 * np.linalg.norm(H4) / np.linalg.norm(H2)


def source(H, beta=4.0, t_ins=1.0):
    e, w = np.linalg.eigh(H)
    u = (w * np.exp(-1j * e * t_ins)) @ w.conj().T
    filt = (w * np.exp(-beta / 2 * (e - e.min()))) @ w.conj().T
    vec = (filt @ BASIS[:, 0].reshape(16, 16)).reshape(256)
    vec /= np.linalg.norm(vec)
    free = lambda v: (u @ v.reshape(16, 16) @ u.T).reshape(256)
    before = free(vec)
    pair = np.zeros((2, 256, 256), complex)
    for m in (0, 1):
        branch = free(before * (RBIT == m))
        for b in (0, 1):
            retained = branch * (LBIT == b)
            jump = np.zeros(256, complex)
            t = IDX[LBIT == b]
            jump[t] = branch[t ^ 128]
            pair[b] += np.outer(retained, retained.conj()) + np.outer(jump, jump.conj())
    return pair[1] - pair[0]


def helstrom(d):
    return 0.5 + 0.25 * float(np.abs(np.linalg.eigvalsh((d + d.conj().T) / 2)).sum())


class Protocol:
    def __init__(self, H, beta=4.0, t_ins=1.0):
        self.src = source(H, beta, t_ins)
        self.e, self.w = np.linalg.eigh(H)

    def F(self, t):
        u = (self.w * np.exp(-1j * self.e * t)) @ self.w.conj().T
        return np.kron(u, u)

    def P(self, t1, mu, t2):
        E = (WV * np.exp(1j * mu * EV)) @ WV.conj().T
        W = self.F(t2) @ E @ self.F(t1)
        full = W @ self.src @ W.conj().T
        return helstrom(np.trace(full.reshape(16, 16, 16, 16), axis1=0, axis2=2))

    def best(self, mu, span=6.0, n=121):
        # F(t2) is local (u x u), so it cannot change the right-side Helstrom value: optimize t1 only.
        grid = np.linspace(-span, span, n)
        vals = [self.P(a, mu, 0.0) for a in grid]
        top = np.argsort(vals)[-3:]
        res = [minimize(lambda p: -self.P(p[0], mu, 0.0), [grid[k]], method='Nelder-Mead',
                        options=dict(maxiter=200, xatol=1e-5, fatol=1e-10)) for k in top]
        r = max(res, key=lambda r: -r.fun)
        return -r.fun, [float(r.x[0])]


def job_seed(args):
    seed, kind, mu = args
    pr = Protocol(hamiltonian(seed, kind))
    if kind == 'zero':
        return dict(seed=seed, kind=kind, mu=mu, plus=pr.P(0, mu, 0), minus=pr.P(0, -mu, 0))
    p, tp = pr.best(mu)
    m, tm = pr.best(-mu)
    return dict(seed=seed, kind=kind, mu=mu, plus=p, minus=m, t_plus=tp, t_minus=tm)


def job_tune(args):
    seed, kind, mu, beta, t_ins = args
    pr = Protocol(hamiltonian(seed, kind), beta, t_ins)
    if kind == 'zero':
        return dict(seed=seed, kind=kind, mu=mu, beta=beta, t_ins=t_ins, plus=pr.P(0, mu, 0), minus=pr.P(0, -mu, 0))
    p, _ = pr.best(mu, n=61)
    m, _ = pr.best(-mu, n=61)
    return dict(seed=seed, kind=kind, mu=mu, beta=beta, t_ins=t_ins, plus=p, minus=m)


def run(fn, jobs, out):
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(json.dumps(json.loads(line)['key']))
    todo = [j for j in jobs if json.dumps(list(j)) not in done]
    print(f'{len(jobs) - len(todo)} done, {len(todo)} to run', flush=True)
    with ProcessPoolExecutor(14) as ex:
        for j, r in zip(todo, ex.map(fn, todo, chunksize=1)):
            r['key'] = list(j)
            with open(out, 'a') as f:
                f.write(json.dumps(r) + '\n')
    print('finished', flush=True)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    mode = sys.argv[1]
    if mode == 'check':
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import noise_sweep as ns_mod
        ns_mod.init(sys.argv[2])
        mine = source(hamiltonian(7001, 'syk4'))
        print('max |mine - repo source| =', float(np.abs(mine - ns_mod.NS['_initial']).max()))
    elif mode == 'seeds':
        mus = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7]
        run(job_seed, [(s, k, m) for s in SEEDS for k in ('syk4', 'syk2', 'zero') for m in mus], sys.argv[2])
    elif mode == 'tune':
        mu = float(sys.argv[3])
        jobs = [(s, k, mu, b, t) for s in SEEDS for k in ('syk4', 'syk2', 'zero')
                for b in (1.0, 2.0, 4.0, 8.0) for t in (0.5, 1.0, 2.0, 3.0)]
        run(job_tune, jobs, sys.argv[2])
