# Quantum information-transfer evidence

Steven Li · evidence subset v0.34 · 14 October 2026

This repository provides inspectable code and simulation data supporting two independent research proposals: complete implementation-cost comparisons for a finite SYK information-transfer model, and uncertainty/provenance/cost checks for small quantum experiments.

All current results are simulations, planning arithmetic or synthetic tests. There is no hardware execution or independent external reviewer response. Code and writing used AI assistance. The task is classical one-bit classification; no backward communication, physical wormhole, arbitrary-state teleportation or quantum advantage is established.

- [Plan A evidence and reproduction](plan-a/README.md)
- [Plan B evidence and reproduction](plan-b/README.md)
- [Per-file provenance and transformations](PROVENANCE.json)
- [License scope and exclusions](LICENSES.md)
- [Recorded local checks](LOCAL_CHECKS.json)

The public subset is deliberately separate from private administrative and funding records. Historical private packets remain unchanged. Frozen data reproduce the registered historical assumptions; they do not replace the proposed integration of reduced preparation, actual topology and measured noise.

The root integrity check is `python -I -S -B VERIFY_PUBLIC.py`. It verifies hashes only. Generated results are written under local results directories, which are excluded from the release manifest and are not testimonials.
