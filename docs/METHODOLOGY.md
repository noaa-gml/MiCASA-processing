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
  local + sign-definite by construction; a bounded ~1.5x within-piece bump)
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
  Fritsch-Carlson construction.
- The flux f = F′ is therefore non-negative (or non-positive)
  **everywhere** — knots and within pieces alike. No sign flips by
  construction, not by clipping.
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
| GPP cell-hour mean | 6.55% | 0.12% |
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
data. This means we accept the +0.04 PgC/yr/yr long-term trend in
CASA-only NEE (verify_v2 Check 15.1) as a real feature of the prior
— it's the inversion's job to correct it.

If MiCASA fluxes are ever used in a context **other** than an
atmospheric inversion (e.g., forward-model comparison vs obs at site
level, or as an ensemble member without further optimization), the ATMC
subtraction may again be appropriate. For our current pipeline it
isn't.

See [proposal #7](PROPOSALS.md) for the full integrate-and-revert
chronology.
