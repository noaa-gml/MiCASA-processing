# Proposals — Architecture Decision Records

Numbered design notes covering landed / staged / proposed / rejected
changes to the pipeline, with the rationale that drove each decision.

**Status legend:**

- `[LANDED]` — code is in-tree; behaviour-preserving by default; opt-in
  via env var.
- `[STAGED]` — diagnostic script is in-tree but must be run on the cluster
  to produce output.
- `[PROPOSED]` — no code change yet; documented for later work.
- `[CONSIDERED, NOT PURSUED]` / `[CONSIDERED, DOMINATED BY ...]` —
  evaluated and rejected with rationale.
- `[REJECTED]` — landed and reverted, or actively decided against.

## (1) [LANDED] Right-edge (and optional left-edge) padding for the NRT fit

The last quadratic piece has no future neighbour, so its curvature is
constrained only by the monthly mean and the slope inherited from the
previous piece. Every time `write_piqs.r` is re-run with one more month
of NRT data, the previous tail coefficients shift, and the published
fluxes for the last month or two of the record are revised.

`write_piqs.r` now accepts `MICASA_PIQS_PAD_RIGHT` and
`MICASA_PIQS_PAD_LEFT` (both default 0). When set to a positive integer
the script extends `x.time` by that many months at the corresponding
end, fills the synthetic months from the per-cell same-calendar-month
climatology of the unpadded data, fits, and strips the pad coefficients
before saving. Output dimensions and `piqsfit.time` are unchanged.
Padding settings are recorded in `piqsfit.meta` inside `fit.piqs.rda`
so downstream readers can detect them.

Recommended starting point for production:

```sh
MICASA_PIQS_PAD_RIGHT=2 MICASA_PIQS_PAD_LEFT=0 Rscript write_piqs.r
```

Note: PCHIP-on-cumulative ([proposal #10](#10-landed-2026-05-04-pchip-on-cumulative-as-the-production-fitter))
uses local Fritsch-Carlson slopes, so this padding is moot for the
PCHIP fitter. Kept here for the legacy PIQS fitter and for reference.

## (2) [LANDED] Active-year refit guard in diurnalize-ERA5.r

The climatology-fallback branch introduces a hard discontinuity at the
boundary between "in fit" and "outside fit" months.
`diurnalize-ERA5.r` prints the fit window and padding metadata on
startup, and — per month — warns if a month with **real** monthly data
lies past the fit edge (its sub-monthly shape would silently come from
coefficient-climatology). Set `MICASA_STRICT_PIQS=1` to escalate that
warning to a hard error — recommended for the NRT cadence.

Climatology months (no real monthly file) do not trigger the guard:
their coefficient-climatology is intentional. The check became
per-month with proposal #14.

## (3) [STAGED] v1 → vNRT handoff sanity check

`diag_v1_vNRT_handoff.r` reads the spliced monthly file, computes
area-weighted global monthly totals of NPP, Rh, and (if present) ATMC,
prints the values straddling the boundary so any step is immediately
visible, and saves a CSV plus a multi-panel PDF. Run from the working
directory after `monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc` has
been (re)built.

Optional env vars `MICASA_DIAG_FROM`, `MICASA_DIAG_TO`,
`MICASA_DIAG_BOUNDARY` (all `YYYYMM`) override the default plot window
202301–202612 with the boundary at 202501.

If a global step is visible, the 1-2 months of PIQS coefficients on
either side are biased and downstream fluxes inherit that.

## (4) [LANDED] Sub-monthly sign-flip logging in diurnalize-ERA5.r

The script does not have an explicit negative-GPP clip; the real risk
is that the spline overshoots above zero (GPP convention is
negative-for-uptake, `gpp = -2*NPP`) or below zero (respiration is
positive-typical) at sub-monthly resolution.

The script now prints a one-line per-month summary giving the count and
percentage of land cells that flipped sign at any hour, and the count
and percentage of cell-hours that flipped:

```
PIQS sign-flip [GPP > 0]:  N / M cells (X%), N / M cell-hours (Y%)
PIQS sign-flip [resp < 0]: N / M cells (X%), N / M cell-hours (Y%)
```

Use these counters to decide whether a fitter change is worth a
bake-off. After the PCHIP switch (proposal #10), the per-month rates
collapsed from ~6.55% to ~0.12% globally.

The verify_v2 suite parses these lines (Check 3.1) to aggregate across
the full record.

## (5) [PROPOSED] Document the original PIQS-vs-pils.2-vs-pics bake-off

`write_piqs.r` retains commented-out calls to `pils.2` and `pics`, but
the implementations of those alternatives were never imported into this
working tree (only `piqs()` landed, in
`ash-code/ccg_idl/john/general/piqs.r.txt`). If we want to redo the
bake-off now that the record is ~25 years instead of ~17, we need to
obtain `pils.2` and `pics` first. Until then the commented-out calls in
the per-cell loop are aspirational, not switchable.

## (6) [PROPOSED] CCGCRV as a diagnostic baseline

The NOAA-GML CCGCRV fit (Thoning/Tans: long-term polynomial + harmonics
+ smoothed residual) is available in-tree at `ash-code/ccgcrv` and
`ash-code/ccg_dataProcessing`. It does not preserve monthly mass like
PIQS / PCHIP, so it's not a drop-in replacement, but running it on a
handful of representative gridcells gives a useful sanity baseline —
particularly for right-edge behaviour, where its harmonic component
extrapolates cleanly and the smoothed residual fades to zero. Not
implemented; depends on whether (1) closes the gap on its own.

## (7) [REJECTED 2026-04-29] ATMC budget closure in NEE

NCCS publishes an "atmospheric correction" (ATMC) field alongside
NPP/Rh/FIRE/FUEL with the file-level comment `NEE = Rh - NPP - ATMC`.
Per Weir et al. 2021a (ACP, [doi:10.5194/acp-21-9609-2021](https://doi.org/10.5194/acp-21-9609-2021)),
ATMC is the Low-order Flux Inversion (LoFI) empirical sink: an additive
correction tuned **annually** so the global biospheric NBE matches the
observed atmospheric CO₂ growth rate. The Weir 2021a parameterization is

```
S_m = α_yr · max(T_m − T_{m-1}, 0)/10 · HR_m
```

with α scaled each year so the area-weighted global total of `S_m`
(added to the baseline NEE) closes against the NOAA-MBL CO₂ growth
rate. Spatially the correction is concentrated in the NH extratropics
during JJA via the dT⁺ weighting; magnitude typically ~3 PgC/yr global.
ATMC accounts for processes CASA does not represent (riverine/coastal
carbon export, CO₂/N fertilization, forest regrowth, Q10 effects on
warming-season respiration).

On 2026-04-29 we tried integrating ATMC — `diurnalize-ERA5.r`
subtracted `atmc.mn` from NEE, `lib/ingest_common.r` picked up ATMC,
`compute_clim.sh` built `ATMCclim.nc`. The verify_v2 Check 15.1 trend
impact was substantial:

| | Slope | Mean NEE |
|---|---|---|
| Without ATMC | +0.0413 PgC/yr/yr | −2.45 PgC/yr |
| With ATMC | −0.0067 PgC/yr/yr | −5.99 PgC/yr |

But the change was **reverted** the same day. These fluxes are consumed
as priors in a global atmospheric inversion that itself assimilates
atmospheric CO₂ measurements. ATMC was tuned to the same observation
class — the global atmospheric CO₂ growth rate. Pre-correcting the
prior with ATMC therefore smuggles observational information from the
data side into the prior, a classic data-leakage / double-dipping
problem: the inversion cannot then independently constrain the
long-term sink because ATMC has already used that constraint upstream.

The "right" picture in our usage: the inversion's atmospheric
assimilation **is** the place where the global growth-rate constraint
enters; the prior should reflect what the offline biospheric model
says **on its own**, and the inversion learns the bias correction from
data. This means we accept the +0.04 PgC/yr/yr long-term trend in
CASA-only NEE as a real feature of the prior — it's the inversion's
job to correct it.

Code state after revert: `lib/ingest_common.r` tracers list is
`NPP/Rh/FIRE/FUEL` (no ATMC); `diurnalize-ERA5.r` computes
`NEE = gpp + resp`; `compute_clim.sh` produces only NPPclim/Rhclim.
Existing `monthly_1x1/*.nc` files still carry the ATMC field (harmless
leftover from the brief integration), and `ATMCclim.nc` sits unused on
disk — both can stay; they just aren't read.

If MiCASA fluxes are ever used in a context other than an atmospheric
inversion (e.g., forward-model comparison vs. obs at site level, or as
an ensemble member without further optimization), the ATMC subtraction
may again be appropriate. For our current pipeline it isn't.

## (8) [LANDED 2026-04-29] Polar-night GPP=0 clip in diurnalize-ERA5.r

Without this clip, the PIQS quadratic component (`qmod.gpp - gpp.mn`)
leaked a small residual into hours where ssrd is identically 0 (~2.6%
of cells in `fluxes_202512.nc` with max |GPP| = 9.4e-9 mol m⁻² s⁻¹ in
the verify suite's Check 12.2).

The clip zeros `gpp` at any cell-hour with `ssrd == 0` before NEE is
summed; resp/qgpp/qresp are unaffected. Mass-conservation gap from
the clip is small (p99 GPP rel diff ~1.5%) and limited to
partial-polar-night latitudes; verify Check 2.2 threshold relaxed
1% → 5% to acknowledge this.

After the PCHIP switch ([#10](#10-landed-2026-05-04-pchip-on-cumulative-as-the-production-fitter))
this clip is technically redundant (PCHIP gives 0 by construction at
ssrd=0 cells), but kept as a defensive belt-and-suspenders.

## (9) [CONSIDERED, NOT PURSUED] Linear PIQS as a sign-flip remedy

Suggested in conversation 2026-04-30: drop the per-piece quadratic to
linear (each segment is a straight line, two coefficients,
integral-preserving by the trapezoidal relation
`y_{i+1} = 2*m_i - y_i`). Within-piece overshoot is impossible (a
straight line cannot bulge above either endpoint), so the U-shaped
sub-monthly sign flips Check 3.1 reports couldn't happen inside a
piece.

Tempting, but does not fix the actual problem cleanly:

a. The polar-night residual the (8) clip addresses is **not** a
   within-piece overshoot. It comes from the diurnalize formula
   `gpp = ssrd*gpp.mn/ssr.mn - gpp.mn + qmod`. At ssrd=0 hours this
   collapses to `qmod - gpp.mn`. Even with a linear `qmod`,
   `qmod(t)` within the segment is not exactly `gpp.mn` at every hour,
   so the residual remains. The clip is the right surgical fix.

b. The integral-preservation recursion `y_{i+1} = 2*m_i - y_i`
   **amplifies** knot-level oscillation when monthly means alternate
   small and large values (e.g., near-zero polar DJF flanked by months
   that aren't). PIQS-quadratic absorbs some of that alternation into
   the curvature degree of freedom; linear has to push it all onto
   the knots. So in exactly the cells we care about most (low-NPP,
   transition-heavy), linear may produce knot values of the wrong
   sign, undoing the within-piece guarantee.

c. C¹ continuity is lost: a kink at every month boundary, which
   downstream consumers (CT, etc.) take as hourly NEE — visible
   midnight-on-the-1st discontinuities in derived diel cycles.

d. Most non-polar sign flips are SSRD-redistribution-driven (the
   `ssrd*gpp.mn/ssr.mn` term, not `qmod`), so a smoother `qmod`
   doesn't help with the bulk of Check 3.1's count.

**Conclusion:** not a free win. The right alternative for "preserve
integral but no overshoot anywhere" is monotone-cubic Hermite on the
cumulative integral (PCHIP-on-cumulative) — see (10).

## (10) [LANDED 2026-05-04] PCHIP-on-cumulative as the production fitter

Build the cumulative monthly integral F at the knot times, apply
Fritsch-Carlson monotone-cubic Hermite interpolation to F (R's
`splinefun(method="monoH.FC")`; scipy's `PchipInterpolator` in the
Python bake-off), then differentiate analytically.

Properties:

- F is monotone non-decreasing (Rh) or non-increasing (negated GPP) by
  Fritsch-Carlson construction.
- The flux f = F′ is therefore non-negative (or non-positive)
  **everywhere** — knots and within pieces alike. No sign flips by
  construction, not by clipping.
- f is a piecewise quadratic (derivative of a piecewise cubic Hermite),
  so the storage layout is identical to PIQS — three coefficients per
  piece. `diurnalize-ERA5.r` needs no change.
- The Fritsch-Carlson slope rule is local (uses neighbouring monthly
  means only), no global solve, ~constant time per cell.
- C¹-smooth at knots (Hermite by construction).
- Mass-preserving by construction.

Slight cost vs PIQS-quadratic: in cells with smoothly-varying monthly
means (no near-zero pieces), PIQS's global smoothness solve can produce
a subtly smoother sub-monthly shape than the locally-determined PCHIP
slopes. This is a third-order aesthetic concern in non-pathological
cells; in the cells we actually care about (polar, semi-arid,
transition months), PCHIP is more sensible because it produces flat
segments at zero rather than oscillating through it.

Bake-off (`bakeoff_pchip.py`) on 6 representative cells confirmed PCHIP
gives 0% flip rate by construction vs PIQS up to 30.91% (AK Tundra),
with absolute flux differences <2e-11 invisible at hourly sampling.
Full-record diurnalize confirmation (25 years, 300 months) shows GPP
cell-hour mean flip rate 6.55% → 0.12% and Rh effectively zero.

`write_pchip.r` is now invoked by both `produce_2025_2026.sh` and
`run_year.sh`; the polar-night clip in `diurnalize-ERA5.r` (note 8) is
now redundant for new diurnalizes but kept as a defensive
belt-and-suspenders. `write_piqs.r` and `write_mss.r` remain in the
tree as selectable alternatives via direct invocation.

## (11) [CONSIDERED, DOMINATED BY PCHIP] Constrained-quadratic PIQS

Add a per-piece non-negativity constraint to PIQS-quadratic: enforce
that the quadratic's vertex value lies on the right side of zero (or
that the vertex is outside the piece interval). Solved as a small QP
per cell.

Rejected in favour of PCHIP-on-cumulative because:

- PCHIP guarantees no sign flip **globally** (within-piece + at knots).
  Constrained-quadratic only enforces within-piece; knot values are
  still set by PIQS's global solve and can still flip sign in cells
  with alternating monthly means.
- PCHIP is closed-form (Fritsch-Carlson is a local algebraic rule);
  constrained-quadratic needs a per-cell QP solve, much more expensive
  at 64,800 cells × 25 years.
- Constrained-quadratic has feasibility risk (non-negativity + integral
  + smoothness may be jointly infeasible in pathological cells),
  requiring a fallback path. PCHIP has no analogous risk.
- Storage layout is identical for both alternatives (piecewise
  quadratic for the flux), so PCHIP wins on equal terms with less work.

`write_mss.r` partially implements the constrained-quadratic idea on
the cumulative F (Monotone Smoothing Spline; see
[`docs/METHODOLOGY.md`](METHODOLOGY.md)). It's retained as a selectable
fitter for cells where one prefers PIQS-style smoothness over PCHIP's
local-slope determinism.

## (12) [LANDED 2026-05-15] FastTrack ERA5 meteo fallback

`diurnalize-ERA5.r` reads hourly ERA5 surface meteo (t2m, ssrd, stl1,
swvl1) to redistribute the smoothed monthly fluxes within each month.
It previously read from a single hardcoded tree:

```
$CARBONTRACKER/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100
```

That tree lags the NRT window — as of 2026-05 it stopped at
2026-01-30, blocking diurnalize of Feb/Mar 2026 even though the
MiCASA monthly product was current through 2026-03.

A second tree, the **FastTrack** product, carries the same data but is
populated sooner during the NRT window:

```
$CARBONTRACKER/METEO/tm5-nc/ec/ea_0005/h06h18tr1/sfc/glb100x100
```

It reaches 2026-03 against the primary's 2026-01; where the two
overlap the files are byte-identical.

`diurnalize-ERA5.r` now consults both. `resolve.era5.source()`
resolves each day to the first tree holding all four variables for
that day — the primary tree is always preferred, FastTrack fills the
trailing gap. A day is read wholly from one tree, so provenance stays
clean. Both roots are overridable via `MICASA_ERA5_DIR` /
`MICASA_ERA5_DIR_FALLBACK`.

Provenance is written to each `fluxes_<YYYYMM>.nc` as global
attributes:

| Attribute | Meaning |
|---|---|
| `meteo_source_primary` | path to the primary tree |
| `meteo_source_fasttrack` | path to the FastTrack tree |
| `meteo_source_by_day` | run-length per-day attribution, e.g. `primary:1-30 fasttrack:31` |
| `meteo_fallback_used` | `yes` if any day used a non-primary tree, else `no` |
| `meteo_source_directory` | kept for back-compat; the tree that supplied the most days |

**Why per-day, not per-month resolution.** Near a tree's coverage edge
a single month genuinely straddles both — 2026-01 resolved to
`primary:1-30 fasttrack:31` because the primary's ssrd ends Jan 30.
Per-day resolution keeps that exact, and the run-length encoding
records it without bloating the attribute.

Landed on `main` and backported to `legacy`. verify_v2 Check 1.4 was
updated to probe both trees. First production use: the 2026-Q1 run on
2026-05-16 (see [`CHANGELOG.md`](../CHANGELOG.md)).

## (13) [LANDED 2026-05-16] compute_clim ported off PyFerret

`compute_clim.sh` built the NPP/Rh modulo-month climatologies
(`monthly_1x1/{NPP,Rh}clim.nc`, consumed by diurnalize-ERA5.r's
climatology-fallback branch) by running PyFerret:

```
let <flux>clim = <flux>[d=1,GT=MONTH_IRREG@MOD]
```

PyFerret is broken on Orion — a NumPy 1.x/2.x ABI mismatch makes
`import pyferret` abort (`_ARRAY_API not found`). The clim files
survived only as a stale build from 2025-06; if they ever needed
regenerating, the pipeline could not do it.

The climatology is just the mean of each calendar month across every
year in the record. `compute_clim.py` computes it with xarray
(`groupby("time.month").mean("time")`) and writes the same schema —
variables `NPPCLIM`/`RHCLIM`, dimensions `(MONTH_IRREG, lat, lon)`,
units gC m-2 s-1 — so diurnalize-ERA5.r picks them up unchanged.
`compute_clim.sh` is now a thin wrapper.

Validation: the Python climatology restricted to 2001-2024 was
compared against the old PyFerret output. The algorithm is exact — a
hand-computed January mean matched `groupby` to 1e-12. The two builds
differ by ≤4.7% at a few cells (median 0.05%), and the difference
grows poleward — the signature of the area-weighted aggregator bug
fixed in the 2026-04 latent-bug sweep. The old PyFerret clim was
built from pre-bugfix monthly data; the new clim, built from the
corrected inputs, is the more accurate of the two.

Removes the last pipeline dependency on PyFerret.

## (14) [LANDED 2026-05-16] Per-month climatology auto-detect

`diurnalize-ERA5.r` decided real-vs-climatology per **year**, from the
`MICASA_CLIM_YEARS` env list (default `"2000 <current year>"`):

- a year in the list -> every month read `NPPclim`/`Rhclim`;
- a year not in the list -> every month read its real monthly file,
  with **no existence check** (`load.ncdf` crashed on a missing file).

Year granularity is wrong for the current NRT year, which is partly
real (early months published) and partly not. Neither blanket choice
works: listing the year climatologises the real months too; not
listing it crashes on the unpublished months. The 2026-Q1 production
run had to be hand-scoped (`MICASA_MONTH_END=3` + `MICASA_CLIM_YEARS=2000`)
to work around this.

`diurnalize-ERA5.r` now decides **per month, by file presence**:

```
real.monthly <- <monthly_1x1>/MiCASA_<ver>_..._<YYYYMM>.nc
use.clim     <- !file.exists(real.monthly)
```

A month with a real monthly file uses it; a month without one falls to
the day-of-year climatology. No year list, no crash, no manual scoping
-- a partially-published year just works (real early months,
climatology for the rest, skipped where ERA5 is also absent).

`MICASA_CLIM_YEARS` is no longer read by `diurnalize-ERA5.r`; it
remains `link_daily_clim.sh`'s knob (which days to fill with
climatology symlinks). The fit-edge guard (proposal #2) became
per-month to match -- it warns only when a month with real data lies
past the PIQS fit window.

Validated: diurnalizing all of 2026 with no `MICASA_MONTH_END` /
`MICASA_CLIM_YEARS` -- Jan-Mar read their real monthly files, Apr-Dec
fell to climatology and were skipped (no ERA5), with no crash.

## (15) [LANDED 2026-05-16] Output provenance metadata

Before this change a `fluxes_*.nc` recorded only a one-line `history`
("Created on &lt;date&gt; by script 'diurnalize-ERA5.r'") and the
FastTrack meteo-source attributes (#12). A downstream user — or a
future maintainer — could not tell which revision of this pipeline
produced a file, from which inputs, or when. For a public data product
that is a real gap.

`lib/provenance.r` and `lib/provenance.py` now build one standard
CF/ACDD global-attribute set:

- **Software** — `processing_pipeline`, `processing_pipeline_url`,
  `processing_pipeline_commit` (full git SHA), `processing_pipeline_version`
  (`git describe --tags --always --dirty`), `processing_step`.
- **Run** — `date_created` (ISO-8601 UTC), `processing_host`.
- **Inputs** — `input_<name>` and `input_<name>_sha256` for each input
  file (the monthly NPP/Rh source or climatology, and the coefficient
  fit `fit.piqs.rda`).
- **Method** — `flux_fit_method` (pchip / piqs / mss, from
  `piqsfit.meta`), `micasa_version`, `flux_from_climatology`.
- **Citation** — `institution`, `references`, `license`, `creator_*`,
  `Conventions = "CF-1.10, ACDD-1.3"`.
- **`history`** — one CF audit line; downstream NCO tools append theirs.

The citation constants live in `lib/provenance.conf` — one
`KEY="VALUE"` file, shell-sourceable and parsed by both helpers — so
the DOI is defined in exactly one place. It ships as `PENDING` until
the archival record is minted.

Writers: `diurnalize-ERA5.r` and `compute_clim.py` call the helper
directly. `daysplitter.sh` relies on `ncks` copying the source file's
global attributes into each daily NEE file (so they inherit the hourly
file's provenance) and adds `daily_split_from` markers. `cat_monthly.sh`
stamps the concatenated monthly file via `stamp_provenance.py`.

`stamp_provenance.py` doubles as a retrofit tool: `--retrofit` adds the
**static** citation subset to a file that predates this change, with an
explicit `provenance_note` — it never fabricates a generating commit or
input checksums it cannot recover. Outputs from the next full pipeline
run carry the complete set.

Verified: R/ncdf4 and Python/xarray netCDF round-trips, 26 + 26 unit
checks (`tests/test_provenance.{r,py}`), and `verify_v2` Section 23.

## (16) [LANDED 2026-05-16] Per-step run manifest

`verify_v2` answered "did each pipeline step run, succeed, and how long
did it take?" by globbing `jobs/*.o*` and regex-scraping the logs.
That is fragile: Check 3.1 once read a stale PIQS log, Check 11.1 hit a
self-referential recursion (a verify log quotes the error strings it
found), and Check 22.1's timing depended on the exact wording of an
`[R]` banner line. Logs are unstructured and accumulate forever.

Pipeline steps now append a structured record to a single
tab-separated `jobs/run_manifest.tsv` via `lib/manifest.sh` (shell) and
`lib/manifest.r` (R). Each record is

    timestamp  step  status  host  commit  elapsed_s  detail

with `status` one of `start` / `ok` / `fail`. `diurnalize-ERA5.r` and
`daysplitter.sh` -- the steps fanned out as separate SBATCH jobs, which
the orchestrator cannot time inline -- record themselves: a `start`
when the worker begins and an `ok` (or `fail`, via an R error handler /
shell `EXIT` trap) when it finishes, carrying the elapsed seconds and
the year. The orchestrators `run_year.sh` and `produce_2025_2026.sh`
record every stage they run.

The helpers are deliberately failure-tolerant: a logging call must
never abort a pipeline run, so `manifest_record` runs its body in a
guarded subshell and `manifest.record` wraps everything in `tryCatch`.
Both are safe to call under `set -e`.

`verify_v2` now reads the manifest. Check 22.1 (diurnalize wall-time)
takes `elapsed_s` straight from the manifest instead of parsing logs;
new Section 24 verifies the manifest itself (integrity; no failed
steps within the age window). Check 11.1's log error-scan is kept --
it catches a step that died so hard it could not record its own
`fail` (segfault, OOM-kill), which a manifest cannot capture; that is
defense-in-depth, not fragility.

The manifest does not exist until the instrumented pipeline runs once,
so 22.1 / 24.x report INFO until the next run -- there is nothing to
retrofit (per-run timing of past runs is not recoverable in structured
form).

Verified: shell and R round-trips, 15 unit checks
(`tests/test_manifest.r`), and the Section 24 / 22.1 parse logic on a
synthetic manifest.


## (17) [INVESTIGATED 2026-06-18] PPM/linear fitters evaluated; PCHIP retained as default

Prompted by a request to reconsider the V1→V2 (PIQS→PCHIP) switch and adopt an
"integral-preserving linear" smoother to avoid overshoot. Full method survey,
equations, citations, and an empirical scorecard on the 2001–2026 record live in
[`docs/FITTER_COMPARISON.md`](FITTER_COMPARISON.md). Summary of what was added
and concluded:

- **New selectable fitters** (drop-in via `MICASA_FIT_RDA`, `(a,b,c)` storage,
  no default change): `write_ppm.r` (PPM, Colella & Woodward 1984, limited
  piecewise-parabolic) and `write_linmm.r` (minmod/MUSCL integral-preserving
  linear). Cores in `lib/ppm_fit.r`, `lib/linmm_fit.r`, vectorized grid path,
  unit-tested (`tests/test_ppm_fit.r`, `tests/test_linmm_fit.r`).
- **`diurnalize-ERA5.r`** now reads the fit from `MICASA_FIT_RDA` (default
  `fit.piqs.rda`), so fitters can be A/B'd without clobbering the production fit.
- **Findings.** Mean-preservation is exact for all candidates (annual 2020 budget
  PCHIP −2.617 vs PPM −2.612 PgC, 0.2%). Overshoot: PCHIP bounded 1.5×, PPM and
  linear 1.0× (none). 2020 product GPP sign-flips: PCHIP 0.1–0.7%, PPM 0%.
  Daily-fidelity vs MiCASA's own daily product: PPM < PCHIP < minmod < flat.
  NRT-revision footprint: PCHIP 0 / minmod ≤1 / PPM ≤2 months, **PIQS all 302**.
- **Continuous integral-preserving linear is unstable** (trapezoidal recursion
  `y_{i+1}=2m_i−y_i`, pole at Nyquist) — confirms PROPOSAL #9.
- **Decision (2026-06-18):** PCHIP **remains** the V2 default. The real fix —
  moving off PIQS to a *local, sign-definite* method — already happened at V2
  (PCHIP). A paired same-cell test shows PCHIP vs PPM daily fidelity is a
  **statistical tie** (PPM better in 54% of cell-months, 49% in boreal; the
  paired-difference IQR straddles zero), and PCHIP **wins on continuity** (zero
  flux jumps vs PPM's small jumps at most month edges) -- and continuity is the
  smoother's whole purpose. So PCHIP stays default; `write_ppm.r` and
  `write_linmm.r` are selectable via `MICASA_FIT_RDA`. PIQS is unsuitable for an
  NRT product (28% wrong-sign GPP cell-months; global-solve non-locality); MSS
  likewise (overshoots + ~300x slower).

## (18) [IMPLEMENTED 2026-06-18] Area-to-point kriging fitter (with uncertainty)

Implements the one out-of-domain method (survey, FITTER_COMPARISON.md (4.3))
that adds what every spline lacks — a principled prior-uncertainty — while
keeping exact mass conservation. Globality is acceptable for most use cases, so
the global-solve concern that gated the probabilistic methods does not bind.

- **`lib/atpk_fit.r`** — 1-D temporal area-to-point kriging (Kyriakidis 2004; Yoo
  & Kyriakidis 2006). Block data = monthly means; each month's sub-points are
  kriged from a +-W-month window with an exponential covariance, then represented
  as a mass-preserving quadratic so the on-disk `(a,b,c)` layout is unchanged,
  PLUS a per-piece kriging **variance**. Properties (unit-tested,
  `tests/test_atpk_fit.r`, 14 checks): coherence/mass exact to ~1e-16; variance
  >= 0 and scales with flux magnitude; **sign-safe** via a selective per-piece
  flat fallback (0 wrong-sign); dormant (~0-variance) cells handled. Two scaling
  fixes vs a naive build: solve on the UNIT covariance (real flux variance ~1e-13
  otherwise wrecks system scaling) and fit the per-piece quadratic in the
  normalized fraction u in [0,1]. Kriging weights are data-independent, so
  `atpk.window.weights` precomputes them once per window geometry and
  `atpk.apply.series` applies them as a fast filter (windowed result matches the
  full-series solve to ~6e-5; +-6-month window is effectively exact and makes the
  method NRT-local, footprint <= W).
- **`write_atpk.r`** — driver; outputs `fit.piqs.rda` format with `$var` arrays.
  Knobs MICASA_ATPK_{W,NS,RANGE}. NOTE: fitting the covariance range from monthly
  autocorrelation is INAPPROPRIATE (that is the seasonal cycle, ~9 mo => singular
  + over-smooth); a fixed short range (1.5 mo default) is used.
- **`diurnalize-ERA5.r`** — guarded hook: when the fit carries `$var` AND
  MICASA_WRITE_FLUX_SD is set, emits an `NEE_sd` field (1-sigma prior uncertainty
  on the smoothed sub-monthly NEE flux, sqrt(var_GPP+var_RESP), constant within a
  month). Default off; no effect on the PCHIP/PPM/PIQS path.

**Status:** core + driver + hook implemented and unit-tested; the point estimate
matches PCHIP (prototype RMS/env 0.003-0.035), so the value over the default is
the uncertainty, not the central flux. PCHIP remains the deterministic default;
ATP kriging is selected via `MICASA_FIT_RDA=fit.atpk.rda` (+ MICASA_WRITE_FLUX_SD
for the band). Pending Orion e2e run (needs ct + monthly cat). Open refinements:
per-biome variogram/range, and a true selective-QP positivity (quadprog) instead
of the flat fallback in pieces where it triggers.
