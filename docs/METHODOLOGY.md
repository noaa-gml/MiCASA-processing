# Methodology

How monthly NPP and Rh become hourly NEE: smoothing the monthly time
series with a positivity-preserving spline, then redistributing within
each month using ERA5 hourly meteo.

## Smoothing the monthly time series

For each grid cell we fit a smooth function through the multi-year
monthly time series of GPP and total respiration (Rh), preserving the
monthly means by construction. Several fitters live in the tree, all
producing the same on-disk format (`fit.piqs.rda` with three
coefficients per piece per cell); see [`FITTER_COMPARISON.md`](FITTER_COMPARISON.md) for the full comparison:

- **`write_pchip.r`** — **production default** (Fritsch-Carlson monotone cubic;
  local; sign-definite at knots — 16–60× fewer sub-monthly sign flips than PIQS,
  not zero; a bounded ~1.5x within-piece bump)
- **`write_ppm.r`** — selectable alternative (PPM limited parabolic; zero
  overshoot but reintroduces small month-edge discontinuities; daily fidelity
  statistically tied with PCHIP — see (17))
- **`write_linmm.r`** — selectable (minmod/MUSCL integral-preserving linear)
- **`write_atpk.r`** — selectable; area-to-point kriging — exact coherence + a
  kriging-variance **prior-uncertainty** (`$var`); point estimate ≈ PCHIP (see (18))
- **`write_piqs.r`** — legacy; CT2022-documented but overshoots and its global
  solve rewrites the whole record on any NRT revision (unsuitable for NRT; (17))
- **`write_mss.r`** — alternative; slow and overshoots (see (17))

`diurnalize-ERA5.r` consumes whichever wrote `fit.piqs.rda` last (the
fitter records itself in `piqsfit.meta$fitter`).

### PCHIP-on-cumulative (default)

Fritsch-Carlson monotone-cubic Hermite interpolation applied to the
cumulative integral F(t) at the knot times, then differentiated
analytically to get the flux f = F′ as a piecewise quadratic.

Properties:

- F is monotone non-decreasing (Rh) or non-increasing (negated GPP) by
  Fritsch-Carlson construction — at the **knots**.
- The flux f = F′ is therefore sign-definite **at the knots** and
  overwhelmingly so in the interiors, but **not everywhere by
  construction**: the derivative quadratic can dip the wrong way
  mid-segment even on strictly single-signed input (reproduced — worst
  interior flux −0.042 on positive data; see
  `fitter_diagnostics/pchip_sign_definiteness.r`). In practice this is a
  **16–60× reduction** in sub-monthly sign flips vs PIQS, not an
  elimination — the small residual (≤0.94% of GPP cell-hours) is cleared
  by the polar-night clip / is negligible elsewhere.
- f is a piecewise quadratic (derivative of a piecewise cubic Hermite),
  so the storage layout matches PIQS — three coefficients per piece
  per cell.
- The Fritsch-Carlson slope rule is local (uses neighbouring monthly
  means only), no global solve, ~constant time per cell.
- C¹-smooth at knots (Hermite by construction).
- Mass-preserving by construction.

Confirmed reduction in sub-monthly sign-flip rates over the full
2001-01..2026-03 record (verify_v2 Check 3.1):

| Metric | PIQS | PCHIP |
|---|---|---|
| GPP cell-hour mean | 6.55% | 0.11% |
| GPP cell-hour max | 14.70% | 0.94% |
| Rh cell-hour mean | 0.122% | 0.0000% |
| Rh cell-hour max | 0.444% | 0.002% |

References:
- Fritsch & Carlson 1980, *Monotone Piecewise Cubic Interpolation*,
  SIAM J. Numer. Anal. 17(2) pp 238-246
- R: `stats::splinefun(method="monoH.FC")`
- Python: `scipy.interpolate.PchipInterpolator`
- See [`bakeoff_pchip.py`](../bakeoff_pchip.py) for the cell-level
  diagnostic and [proposal #10](PROPOSALS.md) for the rationale.

### PIQS — Piecewise Integral Quadratic Splines (legacy)

Rasmussen 1991, *Piecewise Integral Splines of Low Degree*, Computers &
Geosciences 17(9) pp 1255-1263. Implementation in
`ash-code/ccg_idl/john/general/piqs.r.txt`. The fit is computed in
`write_piqs.r`, once per gridcell, across the entire multi-year monthly
record (`x.time` has nmon+1 knots spanning every month from the start
of the record through the most recently ingested month). It is **not**
re-fit per calendar year — transitions between calendar years are
already smooth because the December piece and the following January
piece see each other.

For each month i a quadratic piece

```
f_i(t) = a_i (t − t_i)² + b_i (t − t_i) + c_i
```

is stored. The "Integral" in PIQS means the per-piece integral from t_i
to t_{i+1} equals the monthly mean flux exactly. Adjacent pieces share
their endpoint value at the month boundary (C⁰ continuity at the knots);
derivative continuity is **not** separately enforced.

PIQS overshoots zero in cells with sharp seasonality (verify_v2 Check
3.1: up to ~30% sub-monthly sign flips in boreal/tundra). That's the
reason for the PCHIP switch.

Three coefficient arrays `piqsfit.gpp$a/b/c` (and the same for
`piqsfit.resp`) are saved to `fit.piqs.rda`, together with `piqsfit.time`
(left-edge knot times) and `piqsfit.meta` (padding settings — see
[proposal #1](PROPOSALS.md)).

### MSS — Monotone Smoothing Spline (alternative)

Cubic spline on cumulative F minimizing ∫(F″)² subject to F(t_k) = F_k
and f = F′ ≥ 0 at 8 test points per segment. Solved per-cell as a QP
via the R `quadprog` package.

- Recovers PIQS-level smoothness in cells where positivity isn't
  binding.
- **Caveat (measured 2026-06-18):** the non-negativity constraint binds only
  at the interior test points, NOT at the knots, so MSS still overshoots
  (peak/envelope median ~1.35, max ~1.57) and ~24% of land cells carry a
  wrong-sign GPP knot — it is not the overshoot remedy its name suggests.
  See [`docs/FITTER_COMPARISON.md`](FITTER_COMPARISON.md) §4.1.
- Slower (~180–370 ms/cell vs <1 ms for PIQS / PCHIP), ~hours for the full
  grid. The QP's banded Hessian keeps it NRT-local (footprint ≤1 month).

See [`bakeoff_mss.py`](../bakeoff_mss.py) for the cell-level diagnostic
and [proposal #11](PROPOSALS.md) for why we did not adopt this as the
default.

## Diurnalization

The fit gives a smooth monthly NPP and Rh per grid cell. Within each
month, `diurnalize-ERA5.r` redistributes those totals across the 24
hours of each day using ERA5 surface meteo.

The hourly meteo (t2m, ssrd, stl1, swvl1) is resolved per day from the
primary ERA5 tree, falling back to the FastTrack (`ea_0005`) tree for
the NRT trailing window where the primary is not yet populated. Each
output file records the per-day source in `meteo_source_*` global
attributes — see [proposal #12](PROPOSALS.md).

For each month being diurnalized:

- If `t_i` lies within `[min(piqsfit.time), max(piqsfit.time)]` the
  gridcell's `(a, b, c)` for that month are used directly.
- Otherwise (any month past the right end of the current fit) a
  climatology of the coefficients is used: the per-cell mean of
  `(a, b, c)` across every same-calendar-month entry in the fit.

A month's monthly *mean* (the NPP/Rh the fit is evaluated against) is
likewise taken from the real monthly file when it exists, or the
day-of-year climatology (`NPPclim.nc` / `Rhclim.nc`) when it does not
— decided per month by file presence (proposal #14).

Within the month:

- **GPP** is redistributed in time using ERA5 surface solar downward
  radiation (`ssrd`):

  ```
  gpp[t] = ssrd[t] · gpp.mn / ssr.mn − gpp.mn + qmod.gpp(t)
  ```

  where `gpp.mn` is the monthly mean GPP, `ssr.mn` is the monthly mean
  ssrd, and `qmod.gpp(t)` is the spline's within-month deviation from
  its mean.

- **Total respiration** is redistributed using a Q10 function of 2-m
  temperature (`t2m`), with the monthly mean rescaled to match
  `rtot.mn`.

- **Polar-night clip** (proposal #8): at any cell-hour where
  `ssrd == 0`, `gpp` is forced to 0 before NEE is summed. This
  handles the residual `qmod.gpp − gpp.mn` term that would otherwise
  leak into pure-darkness hours.

- **NEE** is then computed as `nee = gpp + resp` — equivalently
  `Rh − NPP`. We do **not** subtract MiCASA's ATMC field; see below.

Fire and fuel-wood emissions bypass the smoothing entirely and are
taken straight from the MiCASA daily product.

### Why not diurnalize the daily product directly?

MiCASA ships a daily 0.1° NPP/Rh product. A natural question is why we
aggregate it to monthly, fit a sub-monthly spline, and re-impose
structure — rather than diurnalizing the daily fields directly. There
are two layers.

**The daily product can't feed the inversion as-is.** It is a daily
total with no diurnal cycle, and CarbonTracker assimilates afternoon-
biased observations against hourly transport. The day/night structure of
NEE — GPP only in daylight, respiration around the clock — is first-
order: a flat daily flux aliases badly against the sampling. So a
diurnalization step to hourly is mandatory regardless of fitter.

**Given that, why the monthly-mean intermediate rather than diurnalizing
the daily fields?** For the historical record (where the daily data
exist) one *could* diurnalize them directly, preserving daily totals and
imposing only the diurnal shape. That is a defensible alternative; we do
not take it, for three reasons:

1. **Meteorological consistency.** MiCASA's day-to-day variability is
   driven by *its* meteorology (MERRA-2). CarbonTracker's transport runs
   on a *different* meteorology (ERA5). Baking MiCASA's daily wiggles
   into the prior introduces a weather inconsistency the inversion cannot
   reconcile — e.g. an NPP dip on a MERRA-2-cloudy day that ERA5 thinks
   was clear. Taking the robust monthly mean and re-imposing sub-monthly
   and diurnal structure from the assimilation system's *own* meteo keeps
   the prior internally consistent with transport.

2. **NRT homogeneity.** The current month's daily data do not exist at
   production time; the fit's purpose is to extend the product to "now"
   from monthly means (PCHIP's ~1-month footprint). Building the
   historical record by a *different* method than the NRT present would
   plant a discontinuity exactly where trend detection is most sensitive.
   One method across the whole record keeps it homogeneous.

3. **The monthly mean is the trusted quantity.** It integrates out
   MiCASA's internal daily model noise; the inversion needs the seasonal-
   plus-diurnal *shape* layered on reliable monthly totals, not MiCASA's
   day-to-day model weather.

This is also why the daily-truth comparison (evaluating the fit at daily
resolution against MiCASA's native daily product;
[`V1_TO_V2_JUSTIFICATION.md`](V1_TO_V2_JUSTIFICATION.md) §1) is a fair
*check on the fitter's fidelity* rather than an argument to abandon it:
it confirms the reconstruction recovers MiCASA's sub-monthly shape
without importing MiCASA's meteorology into the prior.

### Why calendar-month means, not a rolling mean?

Calendar months are an arbitrary 28–31-day discretization, so a sliding
(rolling) monthly mean as the quantity the fit targets is a natural
thing to consider. We don't, for three reasons:

1. **Integral preservation / mass conservation.** Calendar-month bins
   *partition* time, so "the fit's per-piece integral equals the bin
   mean" is well-posed, and `verify_v2` checks the product by
   re-aggregating to those same months. A rolling mean is a low-pass
   filter over *overlapping* windows that do not partition time, so
   "preserve the rolling means" is over-determined and no longer
   corresponds to a clean conserved monthly total — losing the exact
   mass-conservation guarantee the inversion and the verify harness both
   depend on.

2. **NRT.** A *centered* rolling mean needs data on both sides of each
   point, which does not exist at the trailing edge (the current month) —
   exactly the case the fit is built to handle; a *trailing* rolling mean
   lags by half a window. Calendar means + PCHIP have a clean ~1-month
   NRT footprint with no lag.

3. **No new information, and it partly defeats the purpose.** A rolling
   mean is derived from the same daily data; retaining more sub-monthly
   structure would smuggle *some* of MiCASA's MERRA-2 sub-monthly timing
   back into a prior whose transport runs on ERA5 — the same meteo
   inconsistency the monthly-mean step exists to remove (see above).

The concern a rolling mean usually addresses — jagged month-to-month
steps or a Dec→Jan boundary kink — is already handled: PCHIP-on-
cumulative is C¹ at the knots by construction and uses true month-edge
times, so it accounts for unequal month lengths and gives a smooth curve
between month means without a sliding window.

## Why NEE = Rh − NPP, not Rh − NPP − ATMC

NCCS publishes an "atmospheric correction" (`ATMC`) field alongside
NPP/Rh/FIRE/FUEL with the file-level comment `:comment = "...NEE = Rh
- NPP - ATMC..."`. Per Weir et al. 2021a (ACP, doi:[10.5194/acp-21-9609-2021](https://doi.org/10.5194/acp-21-9609-2021)),
ATMC is the Low-order Flux Inversion (LoFI) empirical sink: an additive
correction tuned **annually** so the global biospheric NBE matches the
observed atmospheric CO₂ growth rate. The Weir 2021a parameterization is

```
S_m = α_yr · max(T_m − T_{m-1}, 0)/10 · HR_m
```

with α scaled each year so the area-weighted global total of S_m closes
against the NOAA-MBL CO₂ growth rate. Spatially, the correction is
concentrated in the NH extratropics during JJA via the dT⁺ weighting;
magnitude typically ~3 PgC/yr global. ATMC accounts for processes CASA
does not represent (riverine/coastal carbon export, CO₂/N fertilization,
forest regrowth, Q10 effects on warming-season respiration).

We tried integrating ATMC on 2026-04-29 — and reverted it the same day.
**These fluxes are consumed as priors in a global atmospheric inversion
that itself assimilates atmospheric CO₂ measurements.** ATMC was tuned
to the same observation class (the global atmospheric CO₂ growth rate).
Pre-correcting the prior with ATMC therefore smuggles observational
information from the data side into the prior — a classic
data-leakage / double-dipping problem: the inversion cannot then
independently constrain the long-term sink because ATMC has already
used that constraint upstream.

The "right" picture in our usage: the inversion's atmospheric
assimilation is the place where the global growth-rate constraint
enters; the prior should reflect what the offline biospheric model
says **on its own**, and the inversion learns the bias correction from
data. This means we ship the +0.04 PgC/yr/yr long-term trend in CASA-only NEE
(verify_v2 Check 15.1) as a **property of the CASA prior — not asserted to be a
real climate feature**: whether it is real (CO₂-fertilization/greening) or a CASA
structural bias, it's the inversion's job to correct it from independent
atmospheric data, and pre-closing it with ATMC would be circular.

If MiCASA fluxes are ever used in a context **other** than an
atmospheric inversion (e.g., forward-model comparison vs obs at site
level, or as an ensemble member without further optimization), the ATMC
subtraction may again be appropriate. For our current pipeline it
isn't.

See [proposal #7](PROPOSALS.md) for the full integrate-and-revert
chronology.

## Uncertainty

MiCASA provides **no native per-pixel uncertainty** — it is a single
deterministic realization (the raw file carries only `NPP/Rh/FIRE/FUEL/ATMC/NEE`).
Any prior uncertainty is therefore *constructed*. We quantify a model-free band of
~3% of the local flux envelope — 0.1° sub-grid heterogeneity (~3.5%, biome-
dependent) plus the across-fitter structural spread (~3%) — emitted as an opt-in
`NEE_sd` field (from the ATP-kriging fitter's variance; see
[`FITTER_COMPARISON.md`](FITTER_COMPARISON.md) §4.3–4.4). This is a **lower bound**
covering sub-monthly redistribution and 1° representativeness; it does **not**
include the dominant term — the model error in the monthly NPP/Rh itself (tens of
%) — which MiCASA does not provide and which the downstream inversion's prior
error covariance carries separately.
