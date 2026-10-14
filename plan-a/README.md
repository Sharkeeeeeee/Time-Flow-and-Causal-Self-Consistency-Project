# Plan A: finite-SYK evidence

The archived 408-row endpoint table covers twelve models in two sequential cohorts, multiple resource/receiver conditions and three finite-calibration draws where specified. These are expected simulator success probabilities, not device counts. Three draws per model are not three new disorder models.

## Fast arithmetic check

```text
python -I -S -B plan-a/REPRODUCE_SUMMARY.py
```

It recomputes the factor-64, finite-2048 nominal aggregate: mean 71.262738%, minimum model mean 69.422806%, and 8/12 model means at least 70%. Assumptions include ideal preparation, all-to-all connectivity and strong joint noise improvements. No fresh full twelve-model propagation is performed by this helper.

## One complete frozen case

Use Python 3.12 and NumPy 2.3.5, for example in a virtual environment. Install NumPy from its official package distribution; this repository does not distribute a library binary or environment.

```text
python -m venv .venv
.venv/Scripts/python -m pip install numpy==2.3.5
.venv/Scripts/python -I -B plan-a/REPRODUCE_FULL_CASE.py
```

These environment paths are for Windows; other platforms use their virtual-environment Python path. Only the recorded Windows run is tested here. The complete historical seed-7001/factor-64 case reconstructs the source, propagates the frozen 229.46-us schedule and evaluates frozen receivers and calibration draws. It takes more time than the fast helper; do not interrupt it merely because its initial output is silent.

Expected finite-calibration endpoints are approximately 69.28%, 69.68% and 69.44%, all below 70%. The oracle endpoint is about 70.36%. Preserving this failure case is part of the evidence. Reproduction does not improve its score.

The alternate equations originate within this project. They are not an independent external review. The schedule and receiver policies are fixed inputs, rather than independently synthesized or optimized by this script. The public wrapper adds relative paths and decodes base64 versions of original NPZ/gzip inputs; decoded original hashes appear in PROVENANCE.json. NumPy NPZ reads disallow pickles.

## Resource ledger

data/source_resources.json records the original and reduced source/local-block counts and their timing limits. In the disclosed source0 comparison, PRx rotations fell from 36,479 to 2,250 and CZ gates from 2,594 to 1,178. The full-case reproduction uses the historical source, not this newly reduced source. The ledger is not the entire serialized gate-list dataset and does not complete a new noisy source-integration study.

Future work still needs matched state/information/intervention checks, preparation and routing on one clock, unused-model validation and measured hardware noise. Public evidence is not a hardware-feasibility certificate or large-system gravitational interpretation.
