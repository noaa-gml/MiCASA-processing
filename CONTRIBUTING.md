# Contributing

Thanks for your interest in the MiCASA-processing pipeline. This file
covers the repository layout, how to run the checks, and a few
conventions specific to this codebase.

## Branches

This repository has **two branches with unrelated histories**:

- **`main`** — active development; the v2 pipeline (PCHIP fitter
  default, the `verify_v2` suite). Work here.
- **`legacy`** — the historical MiCASA_v1 pipeline, preserved for
  archival reproducibility. Do **not** merge it into `main` — the two
  were initialised as separate git repositories.

## Environment

The pipeline runs on an HPC system (Orion) and needs:

- **Python** ≥ 3.9 with the packages in [`requirements.txt`](requirements.txt)
  (`numpy`, `xarray`, `netCDF4`, `pandas`). `lib/provenance.py` is
  standard-library only.
- **R** ≥ 4.0 with `ncdf4` (and `quadprog`, only for the optional
  `write_mss.r` fitter), plus the site `ct` helper library.
- **NCO** (`ncks`, `ncrcat`, `ncdump`) for the concatenate / split steps.

Install the Python side into an isolated environment:

```sh
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

## Running the checks

CI (`.github/workflows/ci.yml`) runs four jobs on every push: Python
byte-compile + notebook-sync, `bash -n` on every shell script, R
`parse()` on every R script, and the behaviour tests. Run the
behaviour tests locally with:

```sh
for t in tests/test_*.py; do python3 "$t"; done
for t in tests/test_*.r;  do Rscript   "$t"; done
```

Each test is self-contained, needs no cluster data, and exits non-zero
on failure. When you fix a bug, add a test that would have caught it —
CI learned the hard way that `bash -n` / `parse()` syntax checks are
not enough.

## The verify_v2 notebook is generated

`tests/verify_v2.ipynb` is a **derived artifact**; its source of truth is
`tests/build_verify_v2.py`. To change a check:

1. Edit `tests/build_verify_v2.py`.
2. Regenerate the notebook: `python3 tests/build_verify_v2.py`.
3. Commit **both** files.

CI fails the build if `tests/verify_v2.ipynb` is out of sync with
`tests/build_verify_v2.py`, so never hand-edit the `.ipynb`.

## Provenance / citation constants

`lib/provenance.conf` is the single source of truth for citation
metadata (institution, pipeline URL, DOI); both `lib/provenance.r` and
`lib/provenance.py` read it. The archival DOI ships as `PENDING` — when
it is minted, update `lib/provenance.conf` **and** `CITATION.cff`
(`grep -rl PENDING .` finds every spot).

## Conventions

- **Commit messages** — short, descriptive, imperative subject lines
  (e.g. *"diurnalize: decide climatology per month, not per year"*).
- **Documentation** — engineering changes go in `CHANGELOG.md`; design
  rationale goes in `docs/PROPOSALS.md` as a numbered ADR.
- **License** — this is a U.S. Government work, in the public domain in the
  United States (17 U.S.C. § 105; see [LICENSE](LICENSE)). By contributing you
  agree that your contribution is released as part of that public-domain work.
