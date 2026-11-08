"""Extract the noiseless 8-qubit unitary of the frozen schedule and probe its structure."""
import json
import sys

import numpy as np

import noise_sweep as ns_mod

ns_mod.init(sys.argv[1])
NS = ns_mod.NS
tl = NS['_timeline']
n = 8


def apply(U, u, sites):
    k = len(sites)
    A = U.reshape([2] * n + [U.shape[1]])
    G = u.reshape([2] * (2 * k))
    A = np.tensordot(G, A, axes=(list(range(k, 2 * k)), sites))
    A = np.moveaxis(A, list(range(k)), sites)
    return A.reshape(2 ** n, -1)


# order gates by start time (ties: virtual rz at the same instant act on disjoint or commuting qubits)
order = sorted(tl, key=lambda g: (g['start_ns'], 0 if g['name'] == 'rz' else 1, g['index']))
U = np.eye(256, dtype=complex)
for g in order:
    U = apply(U, NS['gate_matrix'](g), g['qubits'])
print('unitarity error', float(np.abs(U.conj().T @ U - np.eye(256)).max()))

# check against noiseless transport of the source
NS['CFG']['T'] = NS['CFG']['G'] = float('inf')
src = NS['_initial']
full, _ = NS['transport'](src, tl)
print('U reproduces noiseless transport:', float(np.abs(U @ src @ U.conj().T - full).max()))
np.save('U_sched.npy', U)

# operator Schmidt spectrum across left(0-3)|right(4-7)
T = U.reshape(16, 16, 16, 16).transpose(0, 2, 1, 3).reshape(256, 256)
s = np.linalg.svd(T, compute_uv=False)
s = s / np.linalg.norm(s)
print('operator Schmidt weights (top 12):', np.round(s[:12] ** 2, 4).tolist())
print('operator Schmidt rank (>1e-6):', int((s ** 2 > 1e-6).sum()))
