import os, sys
import advantage_seeds as a
if __name__ == '__main__':
    jobs = [(s, k, mu, b, t) for s in a.SEEDS for k in ('syk4',) for mu in (0.3, 0.4)
            for b in (0.25, 0.5, 1.0) for t in (0.1, 0.25, 0.5)]
    jobs += [(s, 'zero', mu, 1.0, 0.5) for s in a.SEEDS for mu in (0.3, 0.4)]
    a.run(a.job_tune, jobs, 'edge.jsonl')
