# Changelog

Dated engineering entries for the active (`main`) branch. Conceptual /
methodological reasoning lives in [`docs/PROPOSALS.md`](docs/PROPOSALS.md);
this file is for "what landed when, and what numbers it moved."

## 2026-05-16 — per-step run manifest

- **Pipeline steps now write a structured run manifest.** New helpers
  `lib/manifest.sh` and `lib/manifest.r` append `start` / `ok` / `fail`
  records — timestamp, host, git commit, elapsed seconds, detail — to
  `jobs/run_manifest.tsv`. `diurnalize-ERA5.r` and `daysplitter.sh`
  self-record (the SBATCH-fanned steps); `run_year.sh` and
  `produce_2025_2026.sh` record every stage they run. The helpers are
  failure-tolerant — a logging call never aborts the pipeline, even
  under `set -e`.
- **`verify_v2` reads the manifest instead of scraping logs.** Check
  22.1 (diurnalize wall-time) now reads `elapsed_s` from the manifest
  rather than regex-parsing `[R] session elapsed time` out of
  `d-*.o*` logs. New Section 24 verifies the manifest itself
  (24.1 integrity, 24.2 no failed steps). See
  [`docs/PROPOSALS.md` #16](docs/PROPOSALS.md).
- **`tests/test_manifest.r`** — 15 CI-run checks of the manifest record
  format and the never-error-the-caller guarantee.

## 2026-05-16 — PCHIP fitter unit test

- **`pchip.fit.cell` extracted to `lib/pchip_fit.r`.** The
  Fritsch-Carlson PCHIP-on-cumulative fitter core that `write_pchip.r`
  runs per grid cell is now a standalone pure-base-R helper, sourced by
  `write_pchip.r`. No behaviour change — the function body is unchanged.
- **`tests/test_pchip_fit.r`** — 12 CI-run checks pinning the fitter's
  contract on synthetic monthly series: integral / monthly-mean
  preservation (uniform and non-uniform knots), C¹ continuity at knots,
  non-negativity for single-signed input, the GPP sign-flip branch, and
  the constant / all-zero edge cases. Previously the fitter was only
  exercised by a post-production `verify_v2` run.

## 2026-05-16 — output provenance metadata + release prep

- **CF/ACDD provenance stamped into every netCDF the pipeline writes.**
  New helpers `lib/provenance.r` and `lib/provenance.py` build one
  standard global-attribute set — producing software with its git
  commit and `git describe` version, an ISO-8601 processing timestamp,
  the host, input files with SHA-256 checksums, the fitter, and
  citation metadata (`institution`, `references`, `license`,
  `creator_*`, CF/ACDD `Conventions`). `diurnalize-ERA5.r` (hourly
  `fluxes_*.nc`) and `compute_clim.py` (`{NPP,Rh}clim.nc`) call them
  directly; `daysplitter.sh` carries the attributes into the daily NEE
  files (copied by `ncks`) and tags each with `daily_split_from`;
  `cat_monthly.sh` stamps the concatenated monthly file. See
  [`docs/PROPOSALS.md` #15](docs/PROPOSALS.md).

- **`lib/provenance.conf`** — single source of truth for the citation
  constants (institution, pipeline URL, archival DOI), read by both
  helpers and shell-sourceable. The DOI ships as `PENDING`;
  `grep -rl PENDING` finds every spot to update once it is minted.

- **`stamp_provenance.py`** — CLI to write the provenance attributes
  onto an existing netCDF. Used in-pipeline by `cat_monthly.sh`, and
  with `--retrofit` to add the static citation subset to outputs that
  predate this change (it never asserts a generating commit it cannot
  know).

- **Tests + verify section.** `tests/test_provenance.r` and
  `tests/test_provenance.py` (26 checks each, CI-run) cover the conf
  parser, SHA-256, and the attribute builder. `verify_v2` gained
  Section 23 (3 checks) confirming production outputs carry the
  attributes.

- **Release prep.** Added `CITATION.cff` (GitHub "Cite this
  repository"), `CONTRIBUTING.md`, and pinned dependencies
  (`requirements.txt`); `.gitignore` widened to cover regenerated
  `fit.piqs.rda.*` / `verify_*.json` artifacts; tagged `v2.0.0`.

## 2026-05-16 — pipeline-robustness pass

- **`compute_clim` ported off PyFerret.** PyFerret is broken on Orion
  (NumPy ABI mismatch); `compute_clim.sh` is now a thin wrapper around
  `compute_clim.py` (xarray modulo-month mean). Algorithm validated
  exact against a hand-computed mean; see
  [`docs/PROPOSALS.md` #13](docs/PROPOSALS.md). Removes the last
  PyFerret dependency.

- **`compute_daily_clim.sh` hardened.** It built the day-of-year
  climatology with `fls=$(ls <glob>)`, which under `set -e` aborted
  the whole script with a cryptic `ls: cannot access` the moment a
  day's glob was empty (hit mid-aggregate in the 2026-05 v1 run).
  Now globs into a `nullglob` array and, on a genuinely empty day,
  fails with a clear message naming the missing pattern instead of a
  raw `ls` error.

- **`check_bounds` ported off NCO.** `check_bounds.sh` (the global-mean
  flux sanity print run by `cat_monthly.sh`) used `ncwa`, which hits an
  NCO chunking bug on the concatenated record — so `cat_monthly.sh`
  wrapped it in `|| true` and the check effectively never ran.
  `check_bounds.py` recomputes the same crude metric with xarray;
  `check_bounds.sh` is now a thin wrapper. With this and `compute_clim`,
  no pipeline step depends on NCO `ncwa`/PyFerret quirks anymore.

- **GitHub Actions CI added** (`.github/workflows/ci.yml`): three jobs
  — Python byte-compile + `verify_v2.ipynb`-in-sync check, `bash -n`
  on every shell script, R `parse()` on every R script. The data
  pipeline can't run in CI, but syntax/build regressions are now
  caught on every push. Present on both `main` and `legacy`.

- **Behavior tests + a CI `tests` job.** The CI above only checked
  *syntax* — which this session proved insufficient (the
  `compute_daily_clim.sh` quoted-glob bug passed `bash -n` cleanly;
  only running it caught the failure). Added a self-contained `tests/`
  suite that CI runs: `tests/test_compute_clim.py` (the modulo-month
  mean, numpy-only) and `tests/test_era5_meteo.r` (the FastTrack
  resolver + run-length encoder, base R). To make the latter testable,
  the ERA5 path helpers were extracted from `diurnalize-ERA5.r` into
  `lib/era5_meteo.r`; `compute_clim.py`'s core is now the pure-NumPy
  `modulo_month_mean`. diurnalize-ERA5.r was re-smoke-tested after the
  extraction (2026-02, FastTrack fallback, metadata intact).

- **Per-month climatology auto-detect** (`docs/PROPOSALS.md` #14).
  `diurnalize-ERA5.r` chose real-vs-climatology per *year* from
  `MICASA_CLIM_YEARS`, and the real branch had no file-existence check
  -- so a partially-published year forced a choice between
  climatologising its real months or crashing on its unpublished ones
  (the 2026-Q1 run was hand-scoped around this). It now decides per
  *month* by monthly-file presence: real file present -> use it, else
  fall to day-of-year climatology. `MICASA_CLIM_YEARS` is now solely
  `link_daily_clim.sh`'s knob; `produce_2025_2026.sh`'s 2026 step
  dropped its `MICASA_MONTH_END`/`MICASA_CLIM_YEARS` workaround.

- **CI dry-run caught two real bugs:**
  - `lib/bench_compression_diurnal.r` had an `if/else` split across
    lines (a syntax error introduced in the 2026-05-05 public-release
    cleanup) — fixed.
  - `lib/test_gca.r` was a broken, incomplete pre-refactor stub
    (`print()` missing a paren); `archimedes()` is maintained in
    `lib/ingest_common.r`. Removed.

- **verify_v2 Check 6.2** downgraded to INFO. It tested the PIQS
  `PAD_RIGHT=2` edge effect; PCHIP (the production fitter) uses local
  slopes and no edge padding, so the premise is obsolete. The
  v2-vs-v1 sanity invariant lives in Check 6.1.

- **verify_v2 Check 11.1** now scans only logs modified within
  `MICASA_VERIFY_LOG_AGE_DAYS` (default 14). Old experiment /
  superseded-run logs in `jobs/` were flagging the check forever.

`compute_clim.py`/`.sh` and the CI workflow were also ported to the
`legacy` (v1) branch.

## 2026-05-16 — 2026-Q1 production run; verify partial-year + fallback fixes

First production run using the FastTrack meteo fallback. Diurnalized
2026-01/02/03 for v2 and day-split them; the v2 archive now spans
2001-01..2026-03 (303 monthly hourly files; ~9220 daily NEE files).

Meteo provenance of the new months (`meteo_source_by_day`):

| Month | Source |
|---|---|
| 2026-01 | `primary:1-30 fasttrack:31` (mixed — primary ssrd ends Jan 30) |
| 2026-02 | `fasttrack:1-28` |
| 2026-03 | `fasttrack:1-31` |

All three read real monthly NPP/Rh via the PCHIP fit (not
climatology); GPP sub-monthly sign-flip 0.017-0.019%.

verify_v2 adjustments the extended archive surfaced:

- **Check 1.4** (ERA5 meteo coverage) probed only the primary tree,
  so the FastTrack-sourced 2026-02/03 files flagged as missing meteo.
  Now probes primary then fallback; `MET_BASE_FALLBACK` added.
- **Checks 5.1 / 5.2** (annual totals, YoY) summed monthly values per
  year and treated 2026's 3-month partial total as a full year
  (GPP -22.9 PgC/yr, -81% YoY). Now drop years with <12 months before
  the sanity comparison — the idiom 15.1 / 16.2 already use.

`produce_2025_2026.sh`: step 9 split into 9 (diurnalize 2025 full
year) and 10 (diurnalize 2026 Q1, `MICASA_MONTH_END=3`,
`MICASA_CLIM_YEARS=2000`). The stale "2026 intentionally not
diurnalized" header note is gone.

Re-verify after the run + fixes: 0 FAIL, 2 WARN (6.2 edge/interior
and 11.1 stale-log noise, both pre-existing). Check 3.1 spans 303
months, GPP cell-hour flip mean 0.11%.

## 2026-05-15 — FastTrack ERA5 meteo fallback

`diurnalize-ERA5.r` now consults two meteo trees instead of one:

- **primary** — `ec/ea/h06h18tr1/sfc/glb100x100`
- **FastTrack fallback** — `ec/ea_0005/h06h18tr1/sfc/glb100x100`,
  populated sooner during the NRT window (covers ~2 months further
  than the primary as of 2026-05).

Each day is resolved to the first tree holding all four variables
(t2m, ssrd, stl1, swvl1); a day is read wholly from one tree.
Provenance is written to the output file as global attributes
`meteo_source_primary`, `meteo_source_fasttrack`,
`meteo_source_by_day` (run-length, e.g. `primary:1-30 fasttrack:31`),
`meteo_fallback_used` (`yes`/`no`), and `meteo_source_directory`
(kept for back-compat, set to the dominant tree).

Both roots are overridable via `MICASA_ERA5_DIR` /
`MICASA_ERA5_DIR_FALLBACK`. Smoke-tested on 2026-02 (a month the
primary tree lacks entirely): resolved `fasttrack:1-28`, wrote a
clean `fluxes_202602.nc` with `meteo_fallback_used = "yes"`.

## 2026-05-05 — Public-release prep

- Added `LICENSE` (CC0 1.0 Universal) and `README.md` as the GitHub
  front page. Restructured the formerly 901-line README into:
  `README.md` (orientation), `docs/PIPELINE.md`, `docs/METHODOLOGY.md`,
  `docs/PROPOSALS.md`, `CHANGELOG.md`.
- Scrubbed personal email defaults from `config.r` / `config.sh` /
  `ingest.r`; `MAIL_USER` and `BASE_DIR` are now required from env.
- Parameterized site-specific paths in `produce_2025_2026.sh` and
  `lib/bench_compression_diurnal.r`.
- Untracked the regenerable `bakeoff_pchip.log`; added `*.log` to
  `.gitignore`.
- Pushed to `git@github.com:pera-noaa/MiCASA-processing.git`:
  `main` = v2 active dev, `legacy` = v1 historical pipeline (unrelated
  histories).

## 2026-05-04 — PCHIP promoted to production fitter

After full-record confirmation (300 months, 25 years), switched
`produce_2025_2026.sh` and `run_year.sh` from `write_piqs.r` to
`write_pchip.r`. PCHIP-on-cumulative is provably non-negative by
Fritsch-Carlson construction, eliminating sub-monthly sign flips
without requiring the polar-night clip.

verify_v2 Check 3.1 numbers:

| Metric | PIQS | PCHIP |
|---|---|---|
| GPP cell-hour mean | 6.55% | 0.12% (57× reduction) |
| GPP cell-hour max | 14.70% | 0.94% (16× reduction) |
| Rh cell-hour mean | 0.122% | 0.0000% |
| Rh cell-hour max | 0.444% | 0.002% (222× reduction) |

All Section 15 climate-signal checks unchanged: trend +0.0447 PgC/yr/yr,
El Niño anomaly +0.643, COVID effect -0.346 (all consistent with PIQS);
Section 5 globals GPP ∈ [-126.2, -119.8], resp ∈ [117.0, 123.9] PgC/yr.

Also fixed `build_verify_v2.py` Check 3.1 glob — was `d-*-MiCASA*.o*`,
which silently read stale PIQS logs after tagged reruns
(`d-*-pchip.o*`). New version picks the most-recent log per year.

`write_piqs.r` and `write_mss.r` remain in the tree as selectable
alternatives via direct invocation. README note (10) flipped
[PROPOSED] → [LANDED].

## 2026-04-30 — Add write_pchip.r and write_mss.r

Two drop-in alternative fitters next to `write_piqs.r`:

- `write_pchip.r` — Fritsch-Carlson monotone-cubic Hermite on cumulative
  F. R `splinefun(method="monoH.FC")`. Smoke test: 47 sec full grid,
  195 MB.
- `write_mss.r` — Monotone smoothing spline (cubic on cumulative F)
  solved as a per-cell QP via `quadprog`. Smoke test: 53 min full grid,
  239 MB.

Both produce `fit.piqs.rda` with `piqsfit.meta$fitter` recording which
one wrote it. `diurnalize-ERA5.r` consumes any of the three transparently.

Bake-off scripts: `bakeoff_pchip.py`, `bakeoff_mss.py` test on 6
representative cells (Manaus, Hyytiälä, Sahel, Arctic Tundra,
semi-arid Texas, AK Tundra). PCHIP gives 0% flip rate by construction
vs PIQS up to 30.91%, with absolute flux differences <2e-11.

## 2026-04-29 — ATMC integration tried and reverted; polar-night clip lands

**Tried (and reverted same day):** Subtracting MiCASA's `ATMC` field
from NEE per the file-level `:comment` formula (`NEE = Rh - NPP - ATMC`).
Reverted because ATMC is tuned to the global atmospheric CO₂ growth
rate and these fluxes feed an inversion that ALSO assimilates
atmospheric CO₂ — pre-correcting the prior with ATMC double-dips.
See [`docs/PROPOSALS.md` #7](docs/PROPOSALS.md) for the full reasoning.

Trend impact during the brief integration:

| | Slope | Mean NEE |
|---|---|---|
| Without ATMC | +0.0413 PgC/yr/yr | -2.45 PgC/yr |
| With ATMC | -0.0067 PgC/yr/yr | -5.99 PgC/yr |

Code state after revert: `lib/ingest_common.r` tracers list back to
`NPP/Rh/FIRE/FUEL`; `diurnalize-ERA5.r` computes `NEE = gpp + resp`;
`compute_clim.sh` produces only NPPclim/Rhclim. Existing monthly
files still carry ATMC (harmless leftover), and `ATMCclim.nc` sits
unused — both can stay; just aren't read.

**Landed:** Polar-night `gpp = 0` clip in `diurnalize-ERA5.r`. Without
this, the spline's quadratic `qmod.gpp - gpp.mn` term leaked a small
residual into hours where ssrd is identically 0 (~2.6% of cells in
`fluxes_202512.nc` with max |GPP| = 9.4e-9 mol m⁻² s⁻¹). The clip
zeros gpp at any cell-hour with `ssrd == 0` before NEE is summed.
verify_v2 Check 12.2 covers this; Check 2.2 threshold relaxed
1% → 5% to acknowledge the ~1.5% mass-conservation gap from the
clip at partial-polar-night latitudes.

## 2026-04-27 — PIQS edge padding, fit-window guard, sign-flip diag

Three additions to make the NRT cadence safer:

- **Edge padding** in `write_piqs.r`: `MICASA_PIQS_PAD_LEFT/RIGHT` env
  vars extend `x.time` with synthetic months (filled from same-month
  climatology), fit, then strip pad coefficients before saving. Output
  shape unchanged. Production setting: `RIGHT=2 LEFT=0`. See
  [proposal #1](docs/PROPOSALS.md).

- **Fit-window guard** in `diurnalize-ERA5.r`: prints fit window +
  padding metadata + active diurnalization year on startup; warns if
  the active year extends past the fit edge. `MICASA_STRICT_PIQS=1`
  escalates the warning to a hard error. See
  [proposal #2](docs/PROPOSALS.md).

- **Sub-monthly sign-flip log line** in `diurnalize-ERA5.r`: per-month
  count and percentage of cells / cell-hours where GPP > 0 or resp < 0.
  Drives verify_v2 Check 3.1. See
  [proposal #4](docs/PROPOSALS.md).

## 2026-04-26 — Latent-bug sweep (Tier-1 refactor)

Six bugs found and fixed during the post-refactor audit:

1. **`lib/ingest_common.r:aggregate.to.1x1`** — Latitude-area weights
   were being recycled column-major across a 10×10 sub-block, applying
   them along the **longitude** axis instead of latitude. Inner
   `for (inlon in inlons)` loop was also dead (×10 then ÷10).
   Magnitude depends on field gradient within a 1° block — typically
   <0.01% for smooth fields, growing toward the poles. Fix: build a
   flat length-100 weight vector that correctly assigns
   `gca[inlats[k]]` to every cell at lat-position k. See
   `lib/test_aggregate.r` for verification + regression test.

2. **`run_year.sh:sbatch_wait`** — `--export="ALL,${exports}"` produced
   a trailing comma when called with empty exports. Now passes `"ALL"`
   alone in that case.

3. **`write_piqs.r`** — `load.ncdf()` path was hardcoded to
   `MiCASA_v1_*.nc`. Now sources `config.r` and uses
   `micasa.out.monthly.cat(cfg)`, so it works under
   `MICASA_VERSION=vNRT` too.

4. **`check_hashes.py`** — Directory glob `202[4-5]` silently skipped
   any year outside 2024–2025. Now reads `MICASA_YEAR_START/END` from
   env and globs the requested years. Also added a missing-checksum-file
   warning.

5. **`link_old_micasa_raw.sh`** — Hardcoded `from_weir/...` legacy path,
   which only existed in the 2024 layout. Now auto-detects between
   legacy and 2025+ layouts via a `layout_candidates` array; uses
   absolute paths so the link survives `WORK_DIR` moves.

6. **`check_unchanged.sh`** — Used to silently warn-and-continue when
   the previous-year reference was missing; new years would slip
   through unchecked. Now: clearer warning with bootstrap instructions,
   and on a clean diff it auto-blesses the new year's file as next
   year's reference (the chain bootstraps itself once the initial 2024
   reference is in place).

Also: `link_daily_clim.sh` and `diurnalize-ERA5.r` default
`MICASA_CLIM_YEARS` to `"2000 $(date +%Y)"` instead of
`"2000 $MICASA_YEAR"` — climatology fallback should track *what's
missing on disk right now*, not which year you happen to be processing.

## 2026-04-26 — Performance: compression-level tuning

`diurnalize-ERA5.r` writes 12 ~660 MB files per year (9 hourly vars at
1°). Default deflate level was 9; bench results on a real
`fluxes_202401.nc` (`lib/bench_compression_diurnal.r`):

| Level | Time/file | Size MB | Per-year writes |
|---|---|---|---|
| 9 | 108 s | 632 | 1298 s (reference) |
| 6 | 72 s | 633 | 870 s (-33%) |
| 4 | 65 s | 634 | 786 s (-39%, +0.3% size) |
| 3 | 60 s | 646 | 715 s (-45%, +2.2% size) |
| 1 | 55 s | 646 | 654 s (-50%, +2.2% size) |

Chose **level 4**: nearly identical file size to level 9 (+0.3%) for
~9 min saved per year on the diurnalize stage. Levels 1-3 buy another
~2 min but cost +14 MB/file = +170 MB/year, not worth it for archived
output.

Ingest paths (`lib/ingest_common.r`, `ingest.r`) left at level 9 since
per-file output is only ~164 KB and the prior bench
(`lib/bench_compression.r`) showed ~9 s/year savings — not worth the
file-size cost for users who pull the daily 1° aggregates.

## 2026-04-26 — Performance: ingest_byyear skip-existing + read-only-needed

Two changes to `ingest_byyear.r` (and a smaller one to
`ingest_monthly.r`):

1. **Skip-existing** — `RECOMPUTE_EXISTING=1` to override (default off).
   mtime-aware: a day is re-ingested if the source `.nc4` is newer
   than the existing 1° output. NASA can republish source files
   (especially vNRT); a pure `file.exists` check would silently keep
   stale aggregates. `wget` in `download.sh` sets local mtime to
   download time, so a re-download of a republished file makes
   `mtime(src) > mtime(out)` and triggers re-ingest on the next
   pipeline pass.

   A daily NRT cycle that adds 1 new day previously deleted and rebuilt
   all 365 daily 1° outputs. Now: re-run skips finished days, processes
   only what's missing or stale.

2. **Read only the 4 needed tracers** (NPP, Rh, FIRE, FUEL) instead of
   the full 6-var raw file (which also has ATMC and NEE). Done by
   passing `vars = micasa.tracers` to `load.ncdf()`.

Measured impact on `ingest_byyear` 2024 (full year, 366 days):

| Run | Wall-time |
|---|---|
| Baseline (vectorized aggregator) | 610 s |
| + read-only-needed (`RECOMPUTE=1`) | 504 s (-17%) |
| + skip-existing (cached re-run) | 4 s (-99%) |

Output is bit-identical (`ncdiff` on 4 sample days × 4 tracers: max
|Δ| = 0). Only the `:history` attribute timestamp differs on rewrite,
as expected.

The vectorized aggregator (commit `ce1bccc`) was the big win that
collapsed `ingest_byyear` from 3.6 hr to ~10 min/year. These two
changes shave another ~17% of the throughput case and ~99% of the
NRT-rerun case.

Verified by:
- `lib/test_ingest_bitident.r` — read-path bit-identity
- `lib/profile_ingest_day.r` — per-step cost breakdown
- `lib/test_aggregate.r` — aggregator regression test (earlier)
