# MiCASA Processing — Land Carbon Flux Pipeline

Scripts to take raw [MiCASA Land Carbon Flux](https://earth.gov/ghgcenter/data-catalog/micasa-carbonflux-grid-v1)
data (global 0.1° NPP, Rh, FIRE, FUEL fields from the MiCASA model v1) and
process it into hourly 1° NEE / NBE products consumed by NOAA's
[CarbonTracker](https://gml.noaa.gov/ccgg/carbontracker/documentation.php)
inverse-modelling system.

## Quick start

Required environment variables:

```sh
export MAIL_USER=you@example.org
export BASE_DIR=/path/to/GFED-CASA/tree
```

Then drive a year through the pipeline:

```sh
# Full pipeline for one year
./run_year.sh 2026

# Near-real-time stream
./run_year.sh 2026 vNRT

# Skip stages whose inputs already exist
./run_year.sh 2026 v1 --skip-download

# Show what would run, without running it
./run_year.sh 2026 --dry-run
```

Stage skip flags: `--skip-download`, `--skip-ingest`, `--skip-aggregate`,
`--skip-piqs`, `--skip-diurnalize`, `--skip-daysplit`. SBATCH stages are
submitted with `--wait` so the driver blocks until completion.

See [`docs/PIPELINE.md`](docs/PIPELINE.md) for full configuration and
stage details.

## Repository layout

This repo has two branches with **unrelated histories**:

- **`main`** — active development. PCHIP fitter (default), the `verify_v2`
  test suite (60 structural / sign / continuity / sanity checks),
  bake-off scripts, plus PIQS and MSS as selectable alternative fitters.
- **`legacy`** — historical MiCASA_v1 pipeline (classic PIQS only).
  Preserved for archival reproducibility; do not try to merge into `main`
  (separate git inits, unrelated histories).

Work on `main`. Use `legacy` only to reproduce a v1-vintage product.

## What's where

| Document | Contents |
|---|---|
| [`docs/V1_TO_V2_JUSTIFICATION.md`](docs/V1_TO_V2_JUSTIFICATION.md) | **Change register: an evidence-backed justification for every change from V1 -> V2**, classified behavior-preserving vs intentional improvement, with the verification (and its scope/limits) for each |
| [`docs/PIPELINE.md`](docs/PIPELINE.md) | Versions (v1 vs vNRT), configuration env vars, flowchart, every program in the tree, data-layout reference, output provenance metadata, NetCDF input schema |
| [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md) | PIQS / PCHIP / MSS fitter theory; diurnalization with ERA5; why NEE = Rh − NPP and not Rh − NPP − ATMC |
| [`docs/FITTER_COMPARISON.md`](docs/FITTER_COMPARISON.md) | Full sub-monthly-fitter bake-off (PIQS/PCHIP/PPM/minmod/MSS/ATP), equations, scorecard, uncertainty; why PCHIP is the default |
| [`docs/DIURNALIZATION_ALTERNATIVES.md`](docs/DIURNALIZATION_ALTERNATIVES.md) | Diurnal-redistribution survey + soil-temp / Lloyd-Taylor respiration prototypes (opt-in), with shadow-diff results |
| [`docs/PROPOSALS.md`](docs/PROPOSALS.md) | Architecture decision records: 18 numbered notes covering landed / proposed / rejected design changes, with rationale |
| [`CHANGELOG.md`](CHANGELOG.md) | Dated engineering entries: latent-bug sweep, performance tuning, ATMC integration arc, PCHIP promotion |
| [`README.notes`](README.notes) | Historical author log (Pera, Jacobson, Weir) — kept for provenance |

## Verification

The `verify_v2.ipynb` notebook runs 60 checks across the pipeline output:
schema, mass conservation across re-aggregation, sign-flip rates,
polar-night clipping, biome-cell sanity, climate-signal consistency
(NEE trend, El Niño anomaly, COVID impact), PCHIP fit invariants,
diurnalize timing, and output-provenance attributes.

```sh
# Build the notebook from source-of-truth
python3 build_verify_v2.py

# Execute (requires WORK_DIR set + ERA5 meteo accessible via $CARBONTRACKER)
python3 run_verify_v2.py verify_v2.ipynb
```

The summary cell prints `PASS=N FAIL=N WARN=N INFO=N` at the end.

## License

[CC0 1.0 Universal](LICENSE) — public domain dedication.

## Citation

If this pipeline supports your work:

- **MiCASA Land Carbon Flux v1** — Brad Weir et al., NASA GSFC.
  DOI [10.5067/ZBXSA1LEN453](https://doi.org/10.5067/ZBXSA1LEN453);
  data catalog at https://earth.gov/ghgcenter/data-catalog/micasa-carbonflux-grid-v1
- **CarbonTracker** — https://gml.noaa.gov/ccgg/carbontracker/
- **This processing code** — cite this repository (CC0; attribution
  appreciated, not required)

Authors: Ash Pera, Andy Jacobson, Brad Weir. See
[`README.notes`](README.notes) for the historical author log.
