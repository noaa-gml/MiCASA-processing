# Pipeline Reference

Versions, configuration knobs, stage flowchart, every program in the tree,
data layout, and the NetCDF input schema.

## Versions: v1 vs vNRT (the hybrid stream)

MiCASA publishes two parallel streams from the same upstream pipeline:

- **`v1`** — Final / authoritative stream. Lags the source data by some
  weeks; what you want for production-quality NEE.
- **`vNRT`** — Near-real-time. Available within days of the source data,
  but may be revised once `v1` lands and supersedes it.

Both streams use version-tagged basenames so they coexist in one tree:

```
portal.nccs.nasa.gov/daily/<YYYY>/<MM>/
    MiCASA_v1_flux_x3600_y1800_daily_<YYYYMMDD>.nc4    ← preferred
    MiCASA_vNRT_flux_x3600_y1800_daily_<YYYYMMDD>.nc4  ← fallback
    …_sha256.txt                                       (also versioned)
```

### Operational policy

Use `v1` wherever it exists; fall back to `vNRT` for the trailing window
where `v1` has not yet been published. Concretely:

1. **Download both streams for the current year:**
   ```sh
   MICASA_VERSION=both ./download.sh
   ```
   (`v1` is downloaded first; `--no-clobber` means `vNRT` only fills gaps.)

2. **Ingest with whichever version you intend to use *for the gap days*.**
   `v1` days that already exist on disk will be ingested as `v1`; `vNRT`
   days fill the rest. Set `MICASA_VERSION=v1` to ingest both as `v1`
   (the file basename then determines which you read), or
   `MICASA_VERSION=vNRT` to keep the vNRT-tagged outputs.

3. **After ingest, expose vNRT-tagged 1° outputs as v1-tagged for any
   downstream consumer that doesn't know about vNRT:**
   ```sh
   ./link_vNRT_to_v1.sh
   ```
   Skips days where `MiCASA_v1_*.nc` already exists, so `v1` always wins.

4. **When the upstream `v1` record catches up,** re-run ingest for those
   days with `MICASA_VERSION=v1` — the `v1` outputs replace the symlinks,
   and CarbonTracker silently switches over.

## Configuration (env-driven)

`config.sh` and `config.r` read the same environment variables, so any
knob can be set on the command line or in `run_year.sh`. **`MAIL_USER`
and `BASE_DIR` are required** (no defaults — set them in your shell
profile).

| Variable | Default | Purpose |
|---|---|---|
| `MAIL_USER` | — (required) | SBATCH `--mail-user` |
| `BASE_DIR` | — (required) | Parent of YYYY trees |
| `MICASA_YEAR` | `2025` | Single-year focus, used by SBATCH workers and `run_year.sh` |
| `MICASA_VERSION` | `v1` | `v1`, `vNRT`, or `both` (`download.sh` only) |
| `MICASA_YEAR_START` | `2001` | First year for multi-year stages |
| `MICASA_YEAR_END` | `${MICASA_YEAR}` | Last year for multi-year stages |
| `MICASA_MONTH_START` | `1` | First month for diurnalize / daysplit |
| `MICASA_MONTH_END` | `12` | Last month for diurnalize / daysplit |
| `MICASA_CLIM_YEARS` | `2000 $(date +%Y)` | Years `link_daily_clim.sh` fills with day-of-year climatology symlinks (diurnalize auto-detects clim per month, see below) |
| `WORK_DIR` | auto-detect | Set explicitly to point at a different checkout |
| `PORTAL_URL_BASE` | NCCS portal | Source download URL |
| `DAILY_1X1_DIR` | `daily_1x1` | Output/work locations — default to subdirs of `WORK_DIR`; **set absolute to redirect results** (see "Output locations" below) |
| `MONTHLY_1X1_DIR` | `monthly_1x1` | " |
| `ERA5_DIR` | `ERA5` | " — hourly `fluxes_*` **and** daily `*.nee.*` |
| `RAW_SRC_DIR` | `portal.nccs.nasa.gov` | " — raw 0.1° download mirror |
| `JOBS_DIR` | `jobs` | " — SLURM logs |

`MICASA_CLIM_YEARS` is consumed only by `link_daily_clim.sh` (which days
to fill with climatology symlinks); its default covers (a) years before
the record starts (≤ 2000) and (b) the current calendar year (NRT
phase). `diurnalize-ERA5.r` no longer uses it — it auto-detects
real-vs-climatology per month from monthly-file presence (proposal #14).
Independent of `MICASA_YEAR` so backfilling an earlier year
doesn't accidentally clim a fully-published year.

### Output locations

By default all results land in subdirectories of `WORK_DIR` (your checkout):
`$ERA5_DIR/` (hourly `fluxes_YYYYMM.nc` **and** the daily `MiCASA_*.nee.*.nc` —
the products CarbonTracker ingests), `$DAILY_1X1_DIR/` and `$MONTHLY_1X1_DIR/`
(the 1° ingested fields), `$JOBS_DIR/` (SLURM logs), and `$RAW_SRC_DIR/` (raw
0.1° downloads). A human-readable `PROVENANCE.txt` (streams, config, git commit,
timestamp) is written into `$ERA5_DIR` by `write_provenance.sh` — see Stage 7.

To send results elsewhere — e.g. heavy output on `/work2` while the repo stays on
`home` — **export absolute paths** for these layout knobs before running.
`config.sh` honours the environment (the defaults apply only when unset), and
`run_year.sh`'s SBATCH `--export=ALL` carries them to the R workers, which read
the same names (`config.r`) and accept absolute paths:

```sh
OUT=/work2/noaa/co2/$USER/micasa_out
mkdir -p "$OUT"/{ERA5,daily_1x1,monthly_1x1,jobs}     # create targets first
export ERA5_DIR="$OUT/ERA5"  DAILY_1X1_DIR="$OUT/daily_1x1" \
       MONTHLY_1X1_DIR="$OUT/monthly_1x1"  JOBS_DIR="$OUT/jobs"
./run_year.sh 2026
```

- **Pre-create the targets** (`mkdir -p`). Ingest creates its
  `daily_1x1`/`monthly_1x1` if absent, but `$JOBS_DIR` must exist before SBATCH
  writes its `%x.o%j` logs.
- These control **output**. The raw-download mirror is `RAW_SRC_DIR`; the ERA5
  **input** meteo is separate (`$CARBONTRACKER/METEO/…`, override via
  `MICASA_ERA5_DIR`).
- `WORK_DIR` is *not* the knob for this — the R scripts source `config.r` and
  `lib/` from it, so it must point at the checkout, not a data dir.

### Runtime-only knobs (driver-set; do not set manually)

| Variable | Purpose |
|---|---|
| `INGEST_YEAR` | Set by `ingest_byyear.r` driver |
| `diurn_year` | Set by `diurnalize-ERA5.r` driver |
| `RECOMPUTE_EXISTING` | `1` to force `ingest_monthly.r` to re-write existing outputs (default: skip them) |
| `MICASA_NO_BLESS_REFERENCE` | `1` to skip auto-blessing this year's downloaded file as next year's reference in `check_unchanged.sh` |
| `MICASA_PIQS_PAD_LEFT` / `_RIGHT` | PIQS edge padding (proposal #1) |
| `MICASA_STRICT_PIQS` | `1` to escalate "year past PIQS fit edge" warning to an error |
| `MICASA_DIURN_OUT_DIR` | A/B-test override for diurnalize output dir |
| `MICASA_DIURN_ONLY_MONTH` | Restrict diurnalize to a single month for shadow-output testing |
| `MICASA_ERA5_DIR` | Override the primary ERA5 meteo tree |
| `MICASA_ERA5_DIR_FALLBACK` | Override the FastTrack (`ea_0005`) fallback ERA5 meteo tree |

## Flowchart

```
run_year.sh
   |
   |--- link_old_micasa_raw.sh (auto-detect 2024 vs 2025+ layout)
   |
   |--- download.sh ---> portal.nccs.nasa.gov/{daily,monthly}/YYYY/...
   |       check_daily_downloads.r
   |       check_hashes.py            (year range from MICASA_YEAR_*)
   |       check_unchanged.sh         (auto-blesses next year's reference)
   |
   |--- ingest_monthly.r       ingest_byyear.r
   |       (both source lib/ingest_common.r — area-weighted aggregator
   |        bug-fixed 2026-04-26; see lib/test_aggregate.r)
   |
   |--- cat_monthly.sh                compute_daily_clim.sh
   |       check_bounds.sh                    |
   |                                  link_daily_clim.sh
   |       compute_clim.sh
   |       write_pchip.r          (or write_piqs.r / write_mss.r)
   |
   |--- diurnalize-ERA5.r
   |
   `--- daysplitter.sh
                |
        link_vNRT_to_v1.sh  (only if MICASA_VERSION=vNRT was used)
```

## Programs

### Drivers

- **`run_year.sh`** — Top-level driver. Sets `MICASA_YEAR`, sources
  `config.sh`, calls each pipeline stage in order. SBATCH stages
  submitted with `--wait`.
- **`produce_2025_2026.sh`** — Phase-2 NRT batch driver: ingest_monthly,
  v1/vNRT symlinking, climatology fill for the trailing partial month,
  `cat_monthly.sh`, `write_pchip.r`, then submits diurnalize drivers for
  2025 (full year) and 2026 Q1. Day-splitting is run separately
  (`daysplit_array.sbatch`).
- **`daysplit_array.sbatch`** — SLURM array wrapper around
  `daysplitter.sh` (one task per year). Useful for parallel backfills.
- **`config.sh` / `config.r`** — Single source of truth for env-driven
  knobs (sourced by every other script).

### Library

- **`lib/ingest_common.r`** — Shared helpers between `ingest_byyear.r`
  and `ingest_monthly.r`. Defines `archimedes()`, `compute.gca()`,
  `aggregate.to.1x1()`, dim/var helpers, `write.netcdf()`, plus the
  constants `micasa.tracers` and `EARTH_RADIUS_M`.
- **`lib/test_aggregate.r`** — Self-contained Rscript test for
  `aggregate.to.1x1`, with regression test against the pre-2026-04-26
  buggy implementation.
- **`lib/bench_compression_diurnal.r`** — Compression-level benchmark
  for the diurnalize output (see CHANGELOG entry).
- **`lib/era5_meteo.r`** — ERA5 hourly-meteo path resolution
  (`era5.relpath`, `resolve.era5.source`, `encode.day.runs`): the
  primary / FastTrack search order. Unit-tested by
  `tests/test_era5_meteo.r`.
- **`lib/pchip_fit.r`** — the PCHIP-on-cumulative fitter core
  (`pchip.fit.cell`, Fritsch-Carlson monotone cubic Hermite) that
  `write_pchip.r` runs per grid cell. Unit-tested by
  `tests/test_pchip_fit.r`.
- **`lib/diurnal.r`** — the diurnalize flux transform (`diurnal.flux`,
  `polar.night.clip`) that `diurnalize-ERA5.r` applies per time slot.
  Unit-tested by `tests/test_diurnal.r`.
- **`lib/mss_fit.r`** — the monotone-smoothing-spline QP fitter core
  (`mss.fit.setup`, `mss.fit.cell`) that `write_mss.r` runs per grid
  cell. Unit-tested by `tests/test_mss_fit.r`.
- **`lib/provenance.r`** / **`lib/provenance.py`** — build the CF/ACDD
  provenance global-attribute set (git commit, timestamp, input
  checksums, citation) for the R and Python netCDF writers
  respectively (see "Output provenance metadata" below). Unit-tested by
  `tests/test_provenance.{r,py}`.
- **`lib/provenance.conf`** — `KEY="VALUE"` citation constants
  (institution, pipeline URL, DOI) read by both helpers; the single
  place the archival DOI is set.
- **`lib/manifest.sh`** / **`lib/manifest.r`** — append `start` / `ok` /
  `fail` records to `jobs/run_manifest.tsv` (the structured run record
  verify_v2 reads instead of globbing logs). Failure-tolerant: a
  logging call never aborts the pipeline. See "Run manifest" below.
  Unit-tested by `tests/test_manifest.r`.

### Stage 1 — Download

- **`download.sh`** — `wget` MiCASA daily + monthly files for
  `$MICASA_YEAR` / `$MICASA_VERSION`. Supports `MICASA_VERSION=both`
  for the hybrid stream (v1 first, vNRT fills gaps via `--no-clobber`).
- **`check_daily_downloads.r`** — Verify NPP, Rh, FIRE, FUEL exist for
  every day in `[MICASA_YEAR_START, MICASA_YEAR_END]`.
- **`check_hashes.py`** — Verify SHA-256 of each downloaded file. Year
  range from env. Handles both `v1` and `vNRT` files in the same
  directory.
- **`check_unchanged.sh`** — Diff `ncdump -h` headers of new vs reference
  (previous year's tree). Catches silent provider-side metadata changes
  (e.g. the 2018 kg→g units flip). Auto-blesses this year's first daily
  / monthly as the *next* year's reference on a clean diff.

### Stage 2 — Ingest (0.1° → 1°)

- **`ingest_byyear.r`** — For a given `INGEST_YEAR`, aggregate every
  day's raw 0.1° NPP/Rh/FIRE/FUEL to 1° via
  `lib/ingest_common.r:aggregate.to.1x1`. Driver mode (no `INGEST_YEAR`)
  fans out one SBATCH per year in `[MICASA_YEAR_START, MICASA_YEAR_END]`.
- **`ingest_monthly.r`** — Plain year-loop monthly aggregator. Skips
  outputs that already exist unless `RECOMPUTE_EXISTING=1`.
- **`ingest.r`** — *Deprecated.* Superseded by the byyear/monthly split.
  Kept for reference.

### Stage 3 — Aggregate / climatology

- **`cat_monthly.sh`** — Concatenate per-month 1° monthly files into a
  single time-stacked `monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly.nc`,
  then stamp it with provenance via `stamp_provenance.py`. Runs
  `check_bounds.sh` afterward as a non-fatal sanity print.
- **`check_bounds.sh`** — Simple unweighted-area average sanity check.
  Not used in production aggregation (that's `aggregate.to.1x1`).
- **`compute_clim.sh`** / **`compute_clim.py`** — Modulo-month
  climatology (mean of each calendar month across all years) of the
  concatenated monthly file. Writes `monthly_1x1/{NPP,Rh}clim.nc`
  (variables `NPPCLIM`/`RHCLIM`). `compute_clim.sh` is a thin wrapper;
  the logic lives in `compute_clim.py` (xarray). It was a PyFerret
  script until 2026-05 — PyFerret is broken on Orion (NumPy ABI
  mismatch). See [proposal #13](PROPOSALS.md).
- **`compute_daily_clim.sh`** — `ncea` across-year average per day-of-year,
  writing `daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_0000<MMDD>.nc`.
- **`link_daily_clim.sh`** — For each year in `$MICASA_CLIM_YEARS`,
  symlink missing daily files to the `0000<MMDD>` clim.

### Stage 4 — PIQS / PCHIP / MSS fit

- **`write_pchip.r`** — *Production default since 2026-05-04.* Per grid
  cell, fit GPP and rtot with PCHIP-on-cumulative (Fritsch-Carlson
  monotone-cubic Hermite, R's `splinefun(method="monoH.FC")`), save to
  `fit.piqs.rda`. Sign-definite at the knots by Fritsch-Carlson
  construction (not everywhere — the derivative quadratic can dip
  mid-segment, ≤0.94% of cell-hours). Confirmed reduction in sub-monthly
  sign flips: GPP cell-hour mean 6.55% → 0.12%, max 14.70% → 0.94%.
  The per-cell fitter core `pchip.fit.cell` lives in `lib/pchip_fit.r`,
  unit-tested by `tests/test_pchip_fit.r`.
  See [`docs/METHODOLOGY.md`](METHODOLOGY.md) and
  [proposal #10](PROPOSALS.md).
- **`write_piqs.r`** — *Legacy alternative, retained.* Classic Piecewise
  Integral Quadratic Spline fit. PIQS overshoots zero in cells with
  sharp seasonality, which motivated the PCHIP switch. Still used for
  historical reproducibility on the `legacy` branch.
- **`write_mss.r`** — *Drop-in alternative.* Cubic spline on cumulative
  F minimizing ∫(F″)² subject to F(t_k) = F_k and f = F′ ≥ 0 at 8 test
  points per segment. Solved per-cell as a QP via the `quadprog`
  package. Recovers PIQS-level smoothness in cells where positivity
  isn't binding, drops sign-flip rate ~5–25× in cells where PIQS
  overshoots, but residual ~1% from constraint discretization.
  Slower (~30–450 ms/cell vs <1 ms for PCHIP/PIQS). Requires the R
  `quadprog` package.
- **`bakeoff_pchip.py`** / **`bakeoff_mss.py`** — Cell-level fitter
  comparisons used to evaluate write_pchip / write_mss before promotion.

### Stage 5 — Diurnalize (apply ERA5 hourly meteo)

- **`diurnalize-ERA5.r`** — Apply ERA5 hourly meteo (ssrd, t2m, stl1,
  swvl1) to the smoothed (PCHIP/PIQS/MSS) monthly fluxes to get hourly
  GPP/RESP/NEE per (year, month). Writes `ERA5/fluxes_<YYYYMM>.nc`.
  Driver mode (no `diurn_year`) fans out per year in
  `[MICASA_YEAR_START, MICASA_YEAR_END]`. Real-vs-climatology is decided
  per month by monthly-file presence: a month with a real monthly
  NPP/Rh file uses it; a month without one falls to the day-of-year
  climatology (`NPPclim.nc`, `Rhclim.nc`). See proposal #14.

  **Meteo source resolution.** Each day is resolved to the first ERA5
  tree that holds all four variables for it: the **primary** tree
  (`ec/ea/h06h18tr1/sfc/glb100x100`) is tried first, then the
  **FastTrack** fallback (`ec/ea_0005/...`), which is populated sooner
  during the NRT window. A day is read wholly from one tree, so
  provenance stays clean. Each output file records where its meteo
  came from in global attributes:

  | Attribute | Meaning |
  |---|---|
  | `meteo_source_primary` | Path to the primary tree |
  | `meteo_source_fasttrack` | Path to the FastTrack fallback tree |
  | `meteo_source_by_day` | Run-length per-day attribution, e.g. `primary:1-30 fasttrack:31` |
  | `meteo_source_directory` | Back-compat: the tree that supplied the most days |
  | `meteo_fallback_used` | `yes` if any day came from a non-primary tree, else `no` |

  Days absent from every tree are dropped from the hourly time axis and
  the output is marked `status = provisional` (see `meteo_partial`).

### Stage 6 — Daysplit

- **`daysplitter.sh`** — Split `ERA5/fluxes_<YYYYMM>.nc` into per-day
  `ERA5/MiCASA_v1.nee.<YYYYMMDD>.nc` files (NEE only). Range from
  `MICASA_YEAR_*` × `MICASA_MONTH_*`. `ncks` copies the source file's
  global attributes, so each daily file inherits the hourly file's
  provenance; `daily_split_from` markers are added (see "Output
  provenance metadata").

### Stage 7 — Provenance stamp

- **`write_provenance.sh`** — Drop a human-readable `PROVENANCE.txt` into the
  output dir (`$ERA5_DIR` by default; pass a dir as `$1` to override). It
  records what the data streams mean (**v1 / vNRT / FastTrack**), when/where/by
  whom it was generated, the resolved run configuration (version, year/month
  range, respiration driver / response function / polar-clip, fitter pointer,
  meteo root, output dirs), and the **git commit / describe / branch**. Reuses
  the citation constants from `lib/provenance.conf`. This is the top-level
  orientation file; it complements (does not replace) the per-file CF/ACDD
  netCDF attributes and `$JOBS_DIR/run_manifest.tsv`. `run_year.sh` runs it last;
  it is best-effort (always exits 0) and safe to run standalone (does not need
  `MAIL_USER` / `BASE_DIR`).

### Symlink helpers

- **`link_old_micasa_raw.sh`** — Auto-detect previous year's raw layout
  (legacy `from_weir/...` or current `portal.nccs.nasa.gov/...`) and
  absolute-path symlink daily/monthly into this year's tree.
- **`link_old_micasa_finals.sh`** — Same idea for the 1° outputs.
- **`link_vNRT_to_v1.sh`** — Symlink ingested `vNRT` daily files as
  `v1`-named files for the same year. `v1` always wins (skips days
  where `MiCASA_v1_*.nc` already exists).

### Verification

- **`build_verify_v2.py`** — Source-of-truth for the verify_v2 notebook;
  `python3 build_verify_v2.py` regenerates `verify_v2.ipynb` from this
  file. 24 sections, 60 distinct checks (Section 23 confirms output provenance;
  Section 24 reads the run manifest).
- **`stamp_provenance.py`** — CLI to write CF/ACDD provenance global
  attributes onto a netCDF; backs `cat_monthly.sh` and, with
  `--retrofit`, stamps pre-existing outputs (see "Output provenance
  metadata").
- **`verify_v2.ipynb`** — Generated. Run via `run_verify_v2.py`.
- **`run_verify_v2.py`** — Execute `verify_v2.ipynb` as a script.
- **`verify_pchip_invariants.r`** — Helper for verify_v2 §18 (PCHIP fit
  invariants: per-segment analytic sign + C¹ continuity at knots).
- **`verify_piqs_invariants.r`** — Helper for verify_v2 §1.2 / §2.1
  (PIQS metadata + integral preservation).
- **`verify_piqs_tail_horizon.r`** — Helper for verify_v2 §16.4 (tail-coef
  propagation horizon).
- **`verify_piqs_tail_stability.r`** — Helper for verify_v2 §11.2
  (tail-coef stability across re-fits).
- **`diag_v1_vNRT_handoff.r`** — Splice-continuity diagnostic ([proposal
  #3](PROPOSALS.md)).

### Deprecated / one-time scripts (kept for reference)

- `download_and_check.sh` — superseded by `run_year.sh` Stage 1
- `create_era5_move.py` — one-time data-move script
- `download-NRT.sh` — merged into `download.sh`

## Data layout

```
portal.nccs.nasa.gov/{daily/<YYYY>/<MM>,monthly/<YYYY>}/
    MiCASA_<VER>_flux_x3600_y1800_<freq>_<...>.nc4
        Created by download.sh. Both v1 and vNRT files coexist here.

daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_<YYYYMMDD>.nc
        Created by ingest_byyear.r (1° area-weighted aggregate of
        NPP, Rh, FIRE, FUEL).

daily_1x1/MiCASA_<VER>_flux_x360_y180_daily_0000<MMDD>.nc
        Day-of-year climatology, created by compute_daily_clim.sh.

monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly_<YYYYMM>.nc
        Created by ingest_monthly.r.

monthly_1x1/MiCASA_<VER>_flux_x360_y180_monthly.nc
        Concatenated multi-year monthly file (cat_monthly.sh).

monthly_1x1/{NPP,Rh}clim.nc
        Modulo-month climatology, created by compute_clim.sh / .py.

ERA5/fluxes_<YYYYMM>.nc
        Hourly diurnalized monthly file, created by diurnalize-ERA5.r.

ERA5/MiCASA_v1.nee.<YYYYMMDD>.nc
        Daily NEE-only files, created by daysplitter.sh.
```

## Output provenance metadata

Every netCDF the pipeline writes carries a CF/ACDD-style set of global
attributes describing how it was produced, so a downstream user can
trace any file back to an exact pipeline revision and its inputs. The
attributes are built by `lib/provenance.r` (R writers) and
`lib/provenance.py` (Python writers) from the shared citation file
`lib/provenance.conf`.

| Attribute | Meaning |
|---|---|
| `Conventions` | `CF-1.10, ACDD-1.3` |
| `title` / `summary` | Human-readable product description |
| `institution`, `creator_name`, `creator_url` | NOAA GML attribution |
| `source` | Pipeline + producing step |
| `references` | Pipeline URL; dataset DOI once registered |
| `license` | `CC0-1.0` |
| `date_created` | ISO-8601 UTC processing timestamp |
| `processing_pipeline` / `_url` | `MiCASA-processing` + repository URL |
| `processing_pipeline_commit` | Full git SHA of the producing code |
| `processing_pipeline_version` | `git describe --tags --always --dirty` |
| `processing_step` | Producing script |
| `processing_host` | Host the step ran on |
| `input_<name>` / `input_<name>_sha256` | Each input file + its SHA-256 |
| `flux_fit_method` | `pchip` / `piqs` / `mss` (hourly files) |
| `history` | CF audit line; NCO tools append their own |

Citation constants — institution, pipeline URL, and the archival
**DOI** — live in `lib/provenance.conf`. The DOI ships as the literal
`PENDING`; set `MICASA_DOI` (and `MICASA_LANDING_PAGE`) there and in
`CITATION.cff` once the archival record is minted (`grep -rl PENDING`
finds every spot).

Which step writes what:

- **`diurnalize-ERA5.r`** stamps `ERA5/fluxes_*.nc` directly, including
  the `input_*` checksums of the monthly NPP/Rh source (or climatology)
  and `fit.piqs.rda`, plus the FastTrack meteo attributes (Stage 5).
- **`compute_clim.py`** stamps `monthly_1x1/{NPP,Rh}clim.nc`.
- **`daysplitter.sh`** uses `ncks`, which copies the source file's
  global attributes — so each daily `MiCASA_*.nee.*.nc` inherits the
  hourly file's provenance — and adds `daily_split_from` /
  `daily_split_tool` markers.
- **`cat_monthly.sh`** stamps the concatenated monthly file via
  `stamp_provenance.py`.

**`stamp_provenance.py`** is a standalone CLI that writes the
provenance attributes onto an existing netCDF (it backs `cat_monthly.sh`).
Run with `--retrofit` to add the static citation subset to outputs
generated before provenance stamping existed; retrofit mode never
asserts a generating commit or input checksums it cannot recover, and
marks the file with a `provenance_note`. Outputs from a fresh pipeline
run carry the complete set. See [proposal #15](PROPOSALS.md).

## Run manifest

Each pipeline step appends a structured record to a single
tab-separated `jobs/run_manifest.tsv` as it runs, via the
`lib/manifest.sh` (shell) and `lib/manifest.r` (R) helpers. This is the
pipeline's machine-readable run log — `verify_v2` reads it instead of
regex-scraping `jobs/*.o*`. Columns:

```
timestamp  step  status  host  commit  elapsed_s  detail
```

`status` is `start` / `ok` / `fail`. The SBATCH-fanned steps
(`diurnalize-ERA5.r`, `daysplitter.sh`) record themselves — a `start`
when the worker begins and an `ok` (or `fail`, via an error handler or
`EXIT` trap) when it finishes, with the elapsed seconds and the year.
The orchestrators `run_year.sh` and `produce_2025_2026.sh` record every
stage they run. The helpers are failure-tolerant: a logging call never
aborts the pipeline, even under `set -e`.

`verify_v2` Check 22.1 reads diurnalize wall-times from the manifest,
and Section 24 verifies the manifest itself. The file is append-only
and lives under the git-ignored `jobs/`; it does not exist until the
instrumented pipeline has run once. See [proposal #16](PROPOSALS.md).

## NetCDF input schema (raw daily file)

For reference; this is what `download.sh` pulls from NCCS, and what
`ingest_byyear.r` reads:

```
netcdf MiCASA_v1_flux_x3600_y1800_daily_20130507 {
dimensions:
        lat = 1800 ;
        lon = 3600 ;
        time = UNLIMITED ; // (1 currently)
        nv = 2 ;
variables:
        double lat(lat) ;
                lat:units = "degrees_north" ;
                lat:long_name = "latitude" ;
        double lon(lon) ;
                lon:units = "degrees_east" ;
                lon:long_name = "longitude" ;
        double time(time) ;
                time:units = "days since 1980-01-01" ;
                time:long_name = "time" ;
                time:bounds = "time_bnds" ;
        double time_bnds(time, nv) ;
                time_bnds:units = "days since 1980-01-01" ;
                time_bnds:long_name = "time bounds" ;
        float NPP(time, lat, lon) ;
                NPP:units = "kg m-2 s-1" ;
                NPP:expressed_as = "carbon" ;
                NPP:long_name = "Net primary productivity" ;
        float Rh(time, lat, lon) ;
                Rh:units = "kg m-2 s-1" ;
                Rh:expressed_as = "carbon" ;
                Rh:long_name = "Heterotrophic respiration" ;
        float FIRE(time, lat, lon) ;
                FIRE:units = "kg m-2 s-1" ;
                FIRE:expressed_as = "carbon" ;
                FIRE:long_name = "Fire emission" ;
        float FUEL(time, lat, lon) ;
                FUEL:units = "kg m-2 s-1" ;
                FUEL:expressed_as = "carbon" ;
                FUEL:long_name = "Fuel wood emission" ;
        float ATMC(time, lat, lon) ;
                ATMC:units = "kg m-2 s-1" ;
                ATMC:expressed_as = "carbon" ;
                ATMC:long_name = "Atmospheric correction" ;
        float NEE(time, lat, lon) ;
                NEE:units = "kg m-2 s-1" ;
                NEE:expressed_as = "carbon" ;
                NEE:long_name = "Net ecosystem exchange" ;

// global attributes:
        :Conventions = "CF-1.9" ;
        :institution = "NASA Goddard Space Flight Center" ;
        :title = "MiCASA Daily NPP Rh ATMC NEE FIRE FUEL Fluxes 0.1 degree x 0.1 degree v1" ;
        :ShortName = "MICASA_FLUX_D" ;
        :VersionID = "1" ;
        :ProcessingLevel = "4" ;
        :IdentifierProductDOI = "10.5067/ZBXSA1LEN453" ;
        :ReadMeURL = "https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/MiCASA_README.pdf" ;
        :comment = "Positive NPP indicates uptake by vegetation. Positive Rh indicates emission to the atmosphere. NEE = Rh - NPP - ATMC, and NBE = NEE + FIRE + FUEL. ATMC adjusts net exchange to account for missing processes and better match long-term atmospheric budgets." ;
}
```

This pipeline does **not** subtract `ATMC` when computing NEE (the
file-level `:comment` notwithstanding) — see
[`docs/METHODOLOGY.md`](METHODOLOGY.md) and
[proposal #7](PROPOSALS.md#7-rejected-atmc-budget-closure-in-nee).

Note there is **no uncertainty field** in the schema — MiCASA is a single
deterministic realization, so we ship no native per-pixel σ. Any prior
uncertainty is *constructed*: a model-free ~3%-of-envelope band (0.1° sub-grid
heterogeneity + across-fitter structural spread), available as an opt-in `NEE_sd`
output, which is a lower bound (redistribution + representativeness) — not the
dominant monthly-NPP/Rh model error. See
[`docs/FITTER_COMPARISON.md`](FITTER_COMPARISON.md) §4.3.
