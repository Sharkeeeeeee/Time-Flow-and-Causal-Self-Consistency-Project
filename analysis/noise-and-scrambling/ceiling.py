"""Decompose the success ceiling of the seed-7001 frozen case into loss sources.

For R in {64, inf}: transport the source, then compare
  * frozen receivers (nominal d0-d2 and oracle), with and without readout error
  * Helstrom bound on the right-side difference operator (best possible receiver)
Also reports the Helstrom bound of the source itself (before transport).
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor

import noise_sweep as ns_mod

SCRIPT = sys.argv[1]


def helstrom(delta, np):
    # equal priors: P = 1/2 + 1/4 * ||rho1 - rho0||_1
    return 0.5 + 0.25 * float(np.abs(np.linalg.eigvalsh((delta + delta.conj().T) / 2)).sum())


def job(R):
    ns_mod.init(SCRIPT)
    NS = ns_mod.NS
    np = NS['np']
    NS['CFG']['T'] = NS['CFG']['G'] = R
    full, _ = NS['transport'](NS['_initial'], NS['_timeline'])
    right = np.trace(full.reshape(16, 16, 16, 16), axis1=0, axis2=2)
    left = np.trace(full.reshape(16, 16, 16, 16), axis1=1, axis2=3)
    out = {'R': R, 'helstrom_right': helstrom(right, np), 'helstrom_left_after': helstrom(left, np),
           'helstrom_both_after': helstrom(full, np)}
    real_assign = NS['assignment']
    for tag, assign in (('readout', real_assign), ('no_readout', lambda: np.eye(16))):
        NS['assignment'] = assign
        for l in ns_mod.LABELS:
            out[f'{l}|{tag}'] = NS['receiver'](right, NS['_policies'][l])[0]
    NS['assignment'] = real_assign
    np.save(f'right_R{R}.npy', right)
    return out


if __name__ == '__main__':
    ns_mod.init(SCRIPT)
    NS = ns_mod.NS
    np = NS['np']
    src = NS['_initial']
    r0 = np.trace(src.reshape(16, 16, 16, 16), axis1=0, axis2=2)
    l0 = np.trace(src.reshape(16, 16, 16, 16), axis1=1, axis2=3)
    print(json.dumps({'source_helstrom_8q': helstrom(src, np), 'source_helstrom_right': helstrom(r0, np),
                      'source_helstrom_left': helstrom(l0, np)}), flush=True)
    with ProcessPoolExecutor(2) as ex:
        for r in ex.map(job, [64, float('inf')]):
            print(json.dumps(r), flush=True)
