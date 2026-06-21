# Sub-monthly flux reconstruction: method comparison

**Status:** decision document · **Date:** 2026-06-18 · **Default fitter: PCHIP**
(unchanged). This doc compares the alternatives and records why PCHIP is kept
and PIQS is unsuitable. · **Scope:** the monthly→sub-monthly flux smoother in
`diurnalize-ERA5.r` (the `fit.piqs.rda` coefficients).

This explains every reconstruction method (equations + citations), their pros
and cons, and the empirical scorecard measured on the real 2001–2026 record and
a full-year 2020 diurnalize. **Headline:** the consequential decision was moving
off PIQS to a *local, sign-definite* fitter — which already happened at V2
(PCHIP). Among the local methods the differences are second-order; PCHIP and PPM
are near-identical on fidelity (PCHIP significantly but **negligibly** better — see §4.6); the decision rests on PCHIP being local, sign-definite, and closed-form.

## Executive summary

This document records a full investigation of the monthly→sub-monthly flux
smoother, prompted by a concern about the V1→V2 (PIQS→PCHIP) switch and a
proposal to use an integral-preserving *linear* fit to avoid overshoot. Every
candidate was benchmarked on the real 2001–2026 record and a full-year 2020
diurnalize; a cross-domain literature survey and four uncertainty analyses were
added. The findings:

1. **The consequential fix already happened at V2.** What matters is moving off
   PIQS to a *local, sign-definite* fitter — done at V2 (PCHIP). PIQS is
   unsuitable for an NRT product on two measured grounds: overshoot → ~11% of
   GPP cell-hours wrong-sign (2020), and its **global solve rewrites all 303
   historical months** on any revision (§5).
2. **Among local methods the differences are second-order.** PCHIP vs PPM is a
   near-identical on product fidelity — a by-cell bootstrap finds PCHIP
   significantly but **negligibly** better (Δ≈0.3%, 95% CI excludes 0; §4.6).
   PCHIP is kept as default for the *sound* reasons (local, sign-definite,
   closed-form, native `(a,b,c)`), not the fidelity margin. PPM and minmod are
   selectable, not improvements.
   Steffen matches PCHIP on max overshoot (1.5×); MSS overshoots (24% wrong-sign) + ~300×
   slower; the bare/continuous integral-preserving *linear* either explodes
   (continuous, §2.4a) or is dominated by PPM (minmod).
3. **"PIQS + positivity" is exactly MSS** and is over-constrained by a
   degrees-of-freedom argument (§2.6); PCHIP/PPM are positivity done right (a
   representation where f≥0 is automatic, not a constraint).
4. **Fitting at native 0.1° then averaging gives no benefit** (1° or 4°×6°, §4.2)
   — for smooth seasonal data the fitter is nearly linear, so fit-then-average
   and average-then-fit nearly commute.
5. **Uncertainty.** The splines are point estimates. A cross-domain survey (§4.3:
   econometric temporal disaggregation, area-to-point kriging, GP-with-integral-
   constraints, penalized composite link, SCOP-splines) found the principled-
   variance methods are all *global*; **area-to-point kriging** is the standout
   (exact mass + QP positivity + native variance) and is now **implemented and
   verified** (§4.4, PROPOSALS #18). The defensible *model-free* uncertainties
   are 0.1° sub-grid spatial heterogeneity (~3.5%) and across-fitter structural
   spread (~3%); the ATP kriging "band" (9–52%) is an assumption-driven prior
   whose magnitude is set by a hand-chosen covariance range (a 3× swing; §4.6),
   not a measured indeterminacy.

**Decision.** PCHIP stays the deterministic default. `write_ppm.r`,
`write_linmm.r`, `write_piqs.r`, and `write_atpk.r` are selectable via
`MICASA_FIT_RDA` (no change to the default). For a principled prior-uncertainty,
use `write_atpk.r` (+ `MICASA_WRITE_FLUX_SD` for the `NEE_sd` field). PIQS and MSS
are rejected.

---

## 1. The problem and where the fit sits

MiCASA delivers **monthly** mean NPP and heterotrophic respiration per 1° grid
cell. Hourly NEE is built with the **Olsen & Randerson (2004)** scheme (the
CASA-GFED / CarbonTracker standard): split NEE into gross fluxes, redistribute
GPP within the month with ERA5 shortwave and respiration with a Q₁₀ of
temperature, conserving the monthly mean:

```
GPP(t) = GPP_mean · SSRD(t)/SSRD_mean
RE(t)  = RE_mean  · Q10(t)/Q10_mean ,   Q10(t) = 1.5^((T2m−273.15)/10)
```

Used as-is this leaves **abrupt steps at month boundaries** (each month is a flat
mean). To remove those steps, CarbonTracker fits a **mean-preserving smooth
curve** through the monthly series and uses its sub-monthly deviation (`qmod`).
That smoother is the subject here. It writes per-piece quadratic coefficients
`(a,b,c)` to `fit.piqs.rda`; `diurnalize-ERA5.r` evaluates, in month `i` of
width `hᵢ` (28–31 days; **all methods use the true non-uniform `hᵢ`** — see §2):

```
qmod(t) = a·(t−tᵢ)² + b·(t−tᵢ) + c
```

so any method producing `(a,b,c)` is a drop-in swap via `MICASA_FIT_RDA`.

### The constraint trilemma — and which property each method relaxes

For interval-mean reconstruction it is **well-established in practice** (Bartlein
`mp-interp` notes; JULES temporal-interpolation docs; and the smoothing-spline
literature, e.g. histosplines, Boneva et al. 1971) — though not a single named
theorem — that:

> **exact mass-conservation + boundedness/no-overshoot + global C⁰/C¹ smoothness
> cannot all hold at once.** A mean-preserving *smooth* fit must overshoot near
> sharp turning points; forcing no-overshoot breaks smoothness/continuity there.

Every method keeps mass and relaxes *something else* — and **which** thing it
relaxes is the whole story:

| Method | mass | relaxes | keeps |
|---|---|---|---|
| PIQS | ✓ | **boundedness** (overshoots, incl. wrong-sign, unbounded) | global smoothness, C⁰ |
| PCHIP | ✓ | strict boundedness → a **bounded** ≤1.5× bump | continuity (C⁰ flux), sign-definite |
| PPM | ✓ | **global continuity** (small jumps at most edges) | no overshoot, smooth (parabolic) |
| minmod-linear | ✓ | **continuity + curvature** | no overshoot |
| piecewise-const | ✓ | **all sub-monthly structure + continuity** | no overshoot |

So PCHIP and PPM sit on opposite sides of the same coin: PCHIP keeps continuity
and accepts a small bounded bump; PPM kills the bump and accepts small
discontinuities. This is the crux of the PCHIP-vs-PPM choice (§4, §6).

Three requirements are non-negotiable for our use: **(1) mass conservation**,
**(2) no wrong-sign flux** (GPP must not be a source), **(3) NRT stability** (a
revised/appended recent month must not rewrite the historical record).

---

## 2. The methods

Per cell, the monthly means `m₁…mₙ` at knot times `x₀…xₙ`, widths `hᵢ = xᵢ−xᵢ₋₁`,
cumulative `Fₖ = Σ_{j≤k} mⱼhⱼ`. **All in-tree fitters use the true `hᵢ`**; the
PPM edge formula below is the uniform-grid form (the ~10% month-length variation
makes the edge interpolation slightly approximate, but mass is exact regardless
— the parabola passes through the cell mean for any edge pair).

### 2.1 Piecewise-constant (mass-conserving baseline)
`f(t)=mᵢ` (linear interpolation of `F`). Mass-conserving, never overshoots, but
`qmod≡0` (no structure) and a hard step at every boundary. The null hypothesis.

### 2.2 PIQS — Piecewise Integral Quadratic Splines (legacy)
**Rasmussen (1991)**, *Computers & Geosciences* 17(9):1255–1263,
doi:[10.1016/0098-3004(91)90027-B](https://doi.org/10.1016/0098-3004(91)90027-B);
the method documented for CarbonTracker CT2022. Per month a quadratic with (i)
exact integral, (ii) C⁰ continuity at knots, (iii) the free DOF fixed by a
**global** smoothness solve (a tridiagonal system coupling all knots).
**Pro:** smooth, mass-conserving, citable standard. **Con:** overshoots through
zero (unphysical sign flips), and its global solve is unstable for NRT (§5).

### 2.3 PCHIP-on-cumulative (current default)
**Fritsch & Carlson (1980)**, *SIAM JNA* 17(2):238–246,
doi:[10.1137/0717021](https://doi.org/10.1137/0717021). Monotone cubic Hermite on
`F`, differentiated → piecewise-quadratic flux. Knot slopes `dₖ` limited so `F`
is monotone (`dₖ=0` at sign changes; else `|dₖ|≤3·min(|mₖ₋₁|,|mₖ|)`). On segment
`k`, `s=(t−xₖ)/hₖ`, `uₖ=mₖ`:
```
f(s) = (6s−6s²)·uₖ + (3s²−4s+1)·dₖ + (3s²−2s)·dₖ₊₁
```
`f` is sign-definite **at the knots** by construction (and overwhelmingly so in
the interiors) — a 16–60× reduction in sign flips vs PIQS, **not zero by
construction**: the derivative quadratic can still dip mid-segment (reproduced,
worst −0.042 on single-signed input; `fitter_diagnostics/pchip_sign_definiteness.r`).
The slope rule is **local**. **Pro:** mass-conserving, strongly (not perfectly)
sign-preserving, **C⁰-continuous (zero flux jumps)**, local, smooth. **Con:** a
bounded ≤1.5× within-piece bump (the `6u·s(1−s)` Hermite term peaks at `1.5u` when
both knot slopes vanish — a month flanked by near-zero neighbours) and a small
residual sign-flip rate (≤0.94% GPP cell-hours) in near-zero cells (§4).

### 2.4 Integral-preserving linear (MUSCL / slope-limited)
**van Leer (1979)**, *JCP* 32(1):101–136,
doi:[10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1);
TVD theory **Harten (1983)**, *JCP* 49(3):357–393,
doi:[10.1016/0021-9991(83)90136-5](https://doi.org/10.1016/0021-9991(83)90136-5).
Two variants, and the distinction is decisive:

- **(a) Continuous linear** forces endpoints to meet → the mass recursion
  `yᵢ₊₁ = 2mᵢ − yᵢ`, a one-parameter family with a pole at the Nyquist frequency:
  month-to-month alternation **resonates and diverges** (measured ~10⁹× the
  envelope, 65–72% sign-flips). **Unusable** (PROPOSALS #9).
- **(b) Discontinuous, slope-limited (minmod):** `f(t)=mᵢ+σᵢ(t−tᵢᶜ)`,
  `σᵢ=minmod(Sᵢ₋₁,Sᵢ)`, `minmod(a,b)=½(sgn a+sgn b)·min(|a|,|b|)`. Mass-conserving
  for any slope; minmod ⇒ **no overshoot (TVD)**; **discontinuous at every edge**.

### 2.5 PPM — Piecewise Parabolic Method (selectable alternative)
**Colella & Woodward (1984)**, *JCP* 54(1):174–201,
doi:[10.1016/0021-9991(84)90143-8](https://doi.org/10.1016/0021-9991(84)90143-8).
Per month a parabola through its mean with **continuous shared edge values** from
monotonized neighbour differences, then a monotonicity limiter (flatten at
extrema; else steepen one edge):
```
δmᵢ = minmod( (mᵢ₊₁−mᵢ₋₁)/2 , 2(mᵢ−mᵢ₋₁) , 2(mᵢ₊₁−mᵢ) )
aᵢ₊½ = mᵢ + ½(mᵢ₊₁−mᵢ) − (δmᵢ₊₁−δmᵢ)/6            (continuous before limiting)
f(ξ) = a_L + ξ(a_R − a_L + a₆(1−ξ)) ,  a₆ = 6(mᵢ − ½(a_L+a_R)) ,  ξ∈[0,1]
```
**Pro:** mass-conserving, **zero overshoot**, piecewise-quadratic (fits `(a,b,c)`),
local. **Con / honest correction:** it is **not** "continuous except at rare
turning points." Measured on this record the limiter resets an edge at **~70% of
month boundaries**, so most edges carry a small flux jump (median 1.8% of the
local envelope; §4). PPM is *approximately* continuous (small jumps), not C⁰.

### 2.6 Other methods evaluated (none adopted)
- **MSS — monotone smoothing spline** (in-tree, `write_mss.r`): cubic on `F`
  minimising `∫(F″)²` s.t. `F′≥0` at interior test points, per-cell QP. **Measured:
  overshoots** (peak/env median 1.35, max 1.57) and ~24% of land cells carry a
  wrong-sign GPP knot (the constraint binds at test points, not knots); ~180–370
  ms/cell (~100–300× PPM/PCHIP). Banded QP ⇒ NRT-local (≤1 mo). Rejected. **Note:** MSS *is* "PIQS with
positivity enforcement" (PROPOSAL #11). It cannot do better, by a
degrees-of-freedom argument: a quadratic piece has 3 DOF, and fixing its
integral + both C0 endpoint values leaves **zero** freedom to enforce
non-negativity — so positivity must be bought by moving knots (the QP), which
at a near-zero-flanked month is over-constrained/infeasible (a non-negative
quadratic cannot integrate to ~0 between positive endpoints). PCHIP/PPM are
positivity-enforcement done correctly — via a representation (monotone cubic
on the cumulative / limited parabola) where f>=0 is automatic, not constrained.
- **Steffen (1990)** monotone cubic, ADS:[1990A&A...239..443S](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S):
  overshoot **identical to PCHIP (1.50)** — the bump is intrinsic to
  monotone-cubic-on-cumulative, not the FC rule. No gain.
- **Unlimited parabolic** (PPM edges, limiter off): overshoots ≤1.25× (vs PPM
  1.00) — shows the limiter does real work. Dominated by PPM.
- **Unconstrained cubic histospline** (natural cubic on `F`; the Boneva/PIQS
  *parent class*, [Boneva, Kendall & Stefanov 1971](https://academic.oup.com/jrsssb/article/33/1/1/7027042)):
  **measured peak/env median 11.4, max 5.7×10⁵, 67.7% wrong-sign** (Runge
  oscillation). The usable histospline members are exactly PIQS/PCHIP/MSS/PPM,
  which we did test; the unconstrained foundation is dominated by everything.
- **Bounded iterative mean-preserving** (Rymes & Myers 2001,
  doi:[10.1016/S0038-092X(01)00052-4](https://doi.org/10.1016/S0038-092X(01)00052-4);
  Wang & Bartlein 2022,
  doi:[10.1175/JTECH-D-21-0154.1](https://doi.org/10.1175/JTECH-D-21-0154.1)):
  iterated 3-point moving average + per-interval mean-restoration + bound clip.
  **Measured 2026-06-18** (`fitter_diagnostics/bounded_iterative.r`, Rymes-Myers,
  daily sub-resolution): mass-preserving, **0 sign-flips** (bound clip), smooth,
  continuous, and — contrary to the a-priori worry — **NRT-local** (footprint
  ≤2 months even at 300 iterations; the per-interval mean-restoration confines a
  revision). It does carry a bounded overshoot that grows with iterations
  (peak/env ~1.24 → 1.45 for 30 → 300 iters), i.e. **PCHIP-like overshoot plus
  guaranteed sign-definiteness**. Genuinely competitive — *not* dominated — but
  not clearly better than PCHIP either, and it is iterative (niter tuning, no
  closed form) and produces point values rather than the native `(a,b,c)`
  quadratic, so adopting it would need a fit/convert step in `diurnalize`.

---

## 3. Context: why sub-monthly structure matters (not itself a fitter argument)

The rectifier literature shows resolving sub-monthly/diurnal structure is
first-order for an inversion — **but that structure comes from ERA5**, not the
smoother. Denning, Fung & Randall (1995), *Nature* 376:240–243,
doi:[10.1038/376240a0](https://doi.org/10.1038/376240a0); Munassar et al. (2025),
*ACP* 25:639, doi:[10.5194/acp-25-639-2025](https://doi.org/10.5194/acp-25-639-2025)
(removing diurnal structure biases fluxes up to ~48% regionally); Gourdji et al.
(2010), *ACP* 10:6151–6167, doi:[10.5194/acp-10-6151-2010](https://doi.org/10.5194/acp-10-6151-2010).

This justifies *having* sub-monthly structure (so piecewise-constant is a real
downgrade), and it justifies the Olsen–Randerson ERA5 redistribution. The
**fitter's** narrower job is only to (a) remove month-boundary steps and (b) add
a physically-plausible smooth seasonal `qmod` **without injecting artifacts**
(wrong-sign flux, spurious extrema). No rectifier result discriminates *between*
the mass-preserving fitters; their differences are second-order corrections on
top of the ERA5 signal.

---

## 4. Empirical scorecard

Production fit + full-year **2020** diurnalize, 1° grid, ~4.4 M land cell-months.
Reproducible from `fitter_diagnostics/` (§7). Uncertainty shown where it matters.

| Metric | PIQS (V1) | **PCHIP (default)** | PPM | minmod-lin | flat |
|---|---|---|---|---|---|
| Mass-conserving | ✓ | ✓ | ✓ | ✓ | ✓ |
| Overshoot peak/env (med / max) | 0.93 / **~10¹⁸** | 0.83 / **1.50** | 0.78 / **1.00** | / 1.00 | / 1.00 |
| GPP sign-flips, 2020 product (cell-hours, max mo) | **~11%** | 0.1–0.7% | **0%** | 0% | 0% |
| Daily-fidelity RMSE/env, GPP — med [IQR] (mean)¹ | 0.086 (18.6)² | 0.081 [.041,.148] (0.151) | 0.079 [.039,.147] (0.149) | 0.094 (0.159) | (0.181) |
| Daily-fidelity RMSE/env, RESP (mean) | 0.128 | 0.128 | 0.125 | 0.130 | 0.145 |
| Within-month structure (grad/env, med) | high | 0.111 | 0.060 | 0.031 | 0 |
| Flux value-continuity (jump/env, med; % edges jumping) | C⁰ | **0 ; 0%** | 0.018 ; **~70%** | 0.10 ; ~93% | 0.24 ; ~100% |
| Interior extrema (% of months) | high | 54–71% | 5–7% | 0% | 0% |
| **NRT-revision footprint** (months rewritten)³ | **all 303** | **0** | ≤2 | ≤1 | 0 |
| Annual global NEE budget, 2020⁴ | ~same | −2.617 PgC | −2.612 PgC | ~same | ~same |
| Month-boundary aliasing vs PCHIP | — | 1.0× | 0.44–1.0× | worse | worst |
| Fit cost | fast | fast | fast | fast | trivial |
| Lineage | Rasmussen 1991 / CT2022 | Fritsch-Carlson 1980 | Colella-Woodward 1984 | van Leer 1979 | — |

¹ Reconstruction at day-midpoints vs MiCASA's **own** daily 1° product, RMS / local
monthly-mean envelope; `fitter_diagnostics/{fidelity_daily,uncertainty_fidelity}.r`.
² PIQS GPP fidelity **median (0.086) is comparable**, but the **mean (18.6)** is
wrecked by a tail of cells where the global solve diverges to ~10¹⁸× the envelope
— fine in most cells, catastrophic in pathological ones (28% of GPP cell-months
carry a wrong-sign knot). All PIQS numbers measured 2026-06-18 on a regenerated
PIQS fit (`write_piqs.r`→`fit.piqs_v1.rda`), on the **same** record/diurnalize as
the others (`fitter_diagnostics/piqs_score.r`, `ERA5_2020_piqs`).
³ Perturb the latest monthly mean +10%, refit, count prior months moving >1%.
PIQS's global solve couples the whole record; MSS (banded QP) is ≤1.
⁴ The 0.2% PCHIP/PPM budget difference is **not** a fitter mass-leak (the fits
conserve to ~1e-16). It is the **diurnalize** discretization: `qmod` sampled at
discrete hours + the polar-night `GPP=0` clip + the SSRD/Q10 weighting perturb
the product's monthly mean ~0.1–0.7%/month. Same effect for all fitters.

**Two framing clarifications the table needs:**
- **"Overshoot" is not one thing.** The enemy is **wrong-sign** flux and
  **unbounded** blow-up (PIQS tail ~10¹⁸; MSS; histospline 5.7×10⁵). A **bounded
  magnitude** bump (PCHIP's ≤1.5×) is *not* an artifact — MiCASA's own daily data
  exceeds the monthly-mean envelope routinely (daily-max/env median ≈ 1.0, 90th ≈
  1.3–1.4), so a sub-monthly peak above the envelope is physically real.
- **PCHIP vs PPM fidelity is near-identical (superseded — see §4.6).** The pooled
  cell-month numbers here (median Δ 0.0006, IQR straddling 0) conflated population
  spread with inferential uncertainty; the proper by-cell bootstrap on the real
  *product* (§4.6) gives a CI that **excludes 0** with PCHIP significantly but
  negligibly better (~0.3%). Either way PPM does not meaningfully beat PCHIP.

### 4.2 Order of operations: fit-at-0.1deg then aggregate (evaluated, no benefit)

Tested whether fitting at native 0.1deg and area-averaging the `(a,b,c)`
coefficients to the target grid beats aggregating to the target grid first then
fitting (the fit is nonlinear, so they differ). Both derived from the same
0.1deg blocks (`fitter_diagnostics/refine_then_average{,_4x6}.r`):

| target | fit-coarse overshoot (med/90/max) | fit-0.1-then-avg | fine lower in |
|---|---|---|---|
| 1deg (10x10 subcells) | 1.07 / 1.24 / 1.50 | 1.07 / 1.29 / 1.50 | 26% of cells |
| 4x6 (2400 subcells) | 1.07 / 1.18 / 1.38 | 1.07 / 1.18 / 1.36 | 7% of cells |

**No overshoot reduction at either scale.** The hoped-for bump-cancellation
(heterogeneous sub-pixel phenology averaging the 1.5x bump down) does not
materialise: for smooth monthly seasonal cycles the fitter is *nearly linear*,
so fit-then-average and average-then-fit nearly commute. Fit-0.1-then-average
does preserve sign-definiteness, but PCHIP-at-1deg already achieves that to
≤0.94% residual (§4). Not
worth the 100x-2400x fit cost + a 0.1deg monthly cat. (Could still matter for a
flux field that is genuinely sub-cell heterogeneous in *shape*, not just level.)

---

### 4.3 Out-of-domain methods & uncertainty quantification

MiCASA provides **no native per-pixel uncertainty** (a single deterministic
realization — vars `NPP/Rh/FIRE/FUEL/ATMC/NEE` only), so any prior σ must be
constructed. A broad cross-domain survey (verified deep-research, 2026-06-18)
asked whether any field outside climate/numerics beats PCHIP/PPM on the five
requirements, especially by adding the **uncertainty estimate** the deterministic
splines (and MiCASA itself) lack.

| Family (field) | mass-exact | positivity | NRT-local | smooth | uncertainty |
|---|---|---|---|---|---|
| Temporal disaggregation — Denton/Denton-Cholette, Chow-Lin, Fernández, Litterman, Proietti state-space (econometrics; R `tempdisagg`) | ✅ sum/avg/first/last | ✗ | ✗ global (GLS/Kalman) | ✅ | ◑ regression SEs |
| **Area-to-point kriging** — Kyriakidis 2004; Yoo & Kyriakidis 2006 (geostatistics) | ✅ coherent/pycnophylactic | ✅ selective QP | ❓ untested (neighborhood) | ✅ | ✅ kriging variance |
| GP w/ integral + inequality constraints — Da Veiga & Marrel 2012; Maatouk & Bay 2017 | ✅ (linear constraint) | ◑→✅ (M&B exact) | ✗ global | ✅ | ✅ posterior variance |
| Penalized composite link / ungrouping — Rizzi, Gampe & Eilers 2015 (R `ungroup`) | ◑ ~0.1% only | ✅ exp-param | ✗ global | ✅ | ◑ SEs (lower bound) |
| SCOP-splines — Pya & Wood 2015 (R `scam`); mboost — Hofner et al. 2016 | ✗ regression | ✅ / ◑ soft | ✗ global | ✅ | ✅ approx-Bayes / bootstrap CIs |

**Findings.** (1) The genuine gain over the splines is **uncertainty
quantification** — the probabilistic methods all add a variance/CI. (2)
**Locality is the binding constraint**: every probabilistic method is a *global*
solve, so a revised recent month perturbs the historical reconstruction — the
NRT property PCHIP has trivially (footprint 0). **No surveyed method is
simultaneously exact-mass + hard-positive + local + smooth + uncertainty-bearing.**
(3) The standout is **area-to-point kriging** (exact mass + QP positivity + native
variance), but its locality is untested (global system unless forced local). (4)
Refuted (3-0): the PCLM "exactly re-aggregates" claim — it conserves mass only to
~0.1%. (5) Gap: hydrology (method-of-fragments, random cascades) and ML-downscaling
angles yielded no surviving verified claims (they produce stochastic realizations,
not a smooth mean-preserving curve).

**A cheap LOCAL uncertainty band (prototype, `fitter_diagnostics/uncertainty_bands.r`).**
Rather than adopt a global method, the uncertainty PCHIP lacks can be added while
*keeping* locality, two ways (measured on the real record / 2020 daily):

- **Structural** — spread across the mass-preserving fitters {PCHIP, PPM, minmod}
  at each sub-monthly point: **median 3% of the local flux envelope, 90th 11%,
  99th 28%** (roughly uniform across biomes). The "which-smoother" ambiguity
  moves the prior only a few %.
- **Bootstrap-PCHIP** — resample days within each month → bootstrap monthly means
  → refit PCHIP → 5–95% band: **median ~1% (tropics) to ~6% (boreal), tails
  12–25%** at sharp-transition months.

Both bands are per-cell + windowed, so they **preserve NRT locality** (unlike the
global probabilistic methods), and they are modest — a few % in the bulk, 10–28%
at sharp seasonal transitions (largest in boreal spring/fall) — which itself shows
the sub-monthly prior is fairly well-constrained.

**Sub-grid spatial heterogeneity (0.1° within 1°; `fitter_diagnostics/subgrid_uncertainty.r`).**
A model-free, data-driven uncertainty: a 1° cell averages ~100 0.1° pixels, and
the spread across them measures how representative the 1° value is of its actual
landscape. Measured: **median ~3.5% of the flux envelope**, strongly biome-
dependent — **temperate ~10%** (mixed crop/forest/urban), tropics ~6%, **boreal
~1%** (uniform tundra/taiga). Same magnitude for the sub-monthly spread (~3%).

**Uncertainty hierarchy (the useful summary).** The four quantified prior-
uncertainty sources are roughly independent and span an order of magnitude:

| source | magnitude (/envelope) | assumption |
|---|---|---|
| ATP kriging variance (sub-monthly *temporal* indeterminacy) | ~9–52% **but range-set** (0.62→0.20 across range 0.5→6 mo; §4.6) | covariance model — a *chosen* range |
| sub-grid spatial heterogeneity (0.1°) | ~3.5% (1% boreal → 10% temperate) | **none — data-driven** |
| structural (across mass-preserving fitters) | ~3% | none |
| bootstrap-PCHIP (monthly-mean sampling) | ~1–6% | resampling model |

So the *defensible* (model-free) prior-uncertainty is the 0.1° sub-grid spatial
term (~3.5%) and the across-fitter structural term (~3%). The ATP "band" can be
made larger or smaller than these at will by changing the covariance range
(§4.6), so it should be read as a *chosen prior assumption*, not a measured
"dominant" uncertainty. Largest spatial term is temperate (heterogeneous landscapes), largest
temporal term is boreal sharp transitions — they peak in different biomes.

**Recommendation.** Keep PCHIP; if a prior-uncertainty is wanted, use the **local
bootstrap/structural band** above (cheap, keeps locality). Pursue a global
probabilistic method (area-to-point kriging is the best candidate) only if a
*principled posterior variance* is required AND a local-neighborhood variant can
be shown to retain the NRT property. Globality is acceptable for most of our
use cases, so this candidate was implemented & verified — see §4.4.

---

### 4.4 Area-to-point kriging — implemented & verified (globality accepted)

Because globality is acceptable for most use cases, the standout candidate —
area-to-point (ATP) kriging — was prototyped in 1-D time
(`fitter_diagnostics/atp_kriging.r`): block data = monthly means; predict
sub-monthly points by ordinary kriging with an exponential covariance (range
1.5 mo), 6 points/month over a 36-month window. Results across biomes:

| property | result |
|---|---|
| **Mass (coherence)** | **exact to ~1e-16** — block-average of point predictions = monthly mean (pycnophylactic, verified) |
| **Point estimate** | **≈ PCHIP**: RMS(kriging − PCHIP)/env = 0.003–0.035 — the central curve is essentially identical |
| **Uncertainty** | **native kriging variance**: ±1.96σ ≈ **9–52% of the flux envelope** (median ~0.4), wider than the bootstrap band (1–6%) because it quantifies true sub-monthly *indeterminacy*, not just monthly-mean sampling. Width is set by the covariance range (a modeling choice). |
| **Positivity** | **not automatic**: 0% (tropics) → 5% (temperate) → **30–37% (boreal dormant)** wrong-sign → the production code uses a selective per-piece flat fallback, which fires on 12%% of land cell-months (boreal 22%%, §4.6) |
| **Robustness** | the ordinary-kriging system is **ill-conditioned for near-dormant (≈0-variance) cells** → needs a ridge/nugget or a skip |

**Verdict.** ATP kriging delivers what no spline does — exact mass **plus a
principled posterior variance** — with a point estimate that *matches PCHIP*, so
adopting it does not change the central flux; the value is purely the
uncertainty. Costs: a covariance/variogram model (range → band width), a QP
positivity step (material in boreal), regularization for dormant cells, and a
per-cell windowed linear solve. **Recommendation:** keep PCHIP as the
deterministic point estimate; if a principled prior-uncertainty is required
(globality OK), ATP kriging is the route — wrap it around the same monthly means,
fit the variogram per biome, add the selective-QP positivity step, and regularize
dormant cells. Open refinements: per-cell/per-biome variogram fitting (the 1.5-mo
range was assumed here) and the QP positivity enforcement.

**Implemented (2026-06-18, PROPOSALS #18).** A production module now exists:
`lib/atpk_fit.r` (core, 14 unit tests), `write_atpk.r` (driver → `fit.piqs.rda`
format + `$var` arrays; knobs `MICASA_ATPK_{W,NS,RANGE}`), and a guarded
`diurnalize-ERA5.r` hook (`MICASA_WRITE_FLUX_SD` → an `NEE_sd` prior-uncertainty
field; default off). Windowed kriging with precomputed data-independent weights
makes it tractable and NRT-local (footprint ≤ W; the windowed result matches the
full-series solve to ~6e-5). Coherence is exact to ~1e-16 and the reconstruction
is sign-safe via a selective per-piece flat fallback. The flat fallback stands in
for a true selective-QP positivity (a follow-up), and the covariance range is a
fixed 1.5 mo (per-biome variogram fitting is a follow-up — fitting it from the
monthly autocorrelation is *wrong*, as that is the seasonal cycle). Select via
`MICASA_FIT_RDA=fit.atpk.rda`; PCHIP remains the deterministic default. **Verified e2e on Orion (2026-06-20):** `write_atpk.r` produced `fit.atpk.rda` (mass 1.7e-15 vs PCHIP, `$var` on 100% of land months, 0 wrong-sign, point estimate RMS/env 0.04–0.06), and a diurnalize with `MICASA_WRITE_FLUX_SD=1` emitted `NEE_sd` (0 GPP sign-flips, sd 0–2.8e-6 mol m⁻² s⁻¹).

---

### 4.5 PIQS-with-linear-fallback-on-overshoot (evaluated, not adopted)

Tested whether PIQS's superior *global* smoothness could be salvaged by keeping
its quadratic where sign-safe and patching the overshooting pieces with a
sign-safe minmod-linear (`fitter_diagnostics/piqs_hybrid.r`):

| property | result |
|---|---|
| PIQS overshoot rate | **29.3%** of land cell-months need patching |
| smoothness (knot deriv-jump·D/env) | **PIQS 0.000 vs PCHIP 0.290** — PIQS's global solve is near-C¹ in the flux (the genuine motivation) |
| sign-safety after fallback | **0.000% wrong-sign** |
| fidelity (2020 daily RMSE/env) | hybrid 0.079 med / 0.139 mean ≈ PCHIP 0.081 / 0.141 (fixes PIQS's 15.7 mean blow-up; tied) |
| **discontinuity at patched edges** | **median ~5× the envelope** — PIQS overshoots by *large* amounts, so the patch sits far from its neighbour values → severe jumps at the ~30% patched edges |
| NRT locality | **inherits PIQS's global solve** — a revised month rewrites the whole record |

**Verdict: not adopted.** It is sign-safe and recovers PIQS's smoothness in the
easy ~70% of pieces, but injects *large* discontinuities (median ~5× envelope)
at exactly the hard transition pieces — far worse than PPM's small edge jumps —
and still inherits PIQS's global non-locality, while only tying PCHIP on fidelity.
Dominated by PCHIP. The PIQS smoothness advantage is real but cannot be captured
without either the overshoot (PIQS) or the patch discontinuities (this hybrid).

---

### 4.6 Second adversarial review — corrections and follow-ups (2026-06-21)

A second hostile review challenged the fidelity metric, the uncertainty
hierarchy, the statistics, and specifics. Re-analyses
(`fitter_diagnostics/{fidelity_product, atpk_range_sweep, critique_followups}.r`)
and the corrections they forced:

- **Fidelity re-scored on the real product (C1).** §4's original daily fidelity
  scored the *bare* fitter quadratic, omitting the ERA5 redistribution the
  product applies. Re-scoring the **diurnalized NEE product** vs MiCASA daily NEE
  (2020): PCHIP 1.186, PPM 1.190, PIQS 1.184 (RMSE / per-cell RMS scale) —
  near-identical. The fitter's `qmod` controls **~95%** of daily-mean NEE
  variance (the seasonal cycle), but all fitters reproduce that 95% identically
  (shared monthly means); the fitter *choice* affects only the sub-monthly shape
  — a sliver — so no daily metric separates them.
- **"Statistical tie" withdrawn; replaced by a CI (M1).** A by-**cell** block
  bootstrap (N=15,724 cells, B=2000) of PCHIP−PPM product fidelity gives
  Δ = **−0.004, 95% CI [−0.0043, −0.0037] — excludes zero**: PCHIP is
  *significantly but negligibly* better (~0.3% of RMSE). The earlier pooled-
  cell-month "tie" (which conflated population spread with median uncertainty)
  is corrected.
- **ATP variance demoted — a chosen prior, not a measurement (C2).** The band is
  set by the covariance range, a hand-set knob: median band/env = **0.62 / 0.53 /
  0.38 / 0.27 / 0.20** for range = 0.5 / 0.75 / 1.5 / 3.0 / 6.0 mo (a 3× swing),
  and the variance is **geometry×sill** (two cells with 60× different sill have
  identical √var/√sill = 0.374). So the "9–52%" is √sill/env re-expressed, not an
  independently measured indeterminacy. **Corrected ranking:** the *model-free*
  terms — 0.1° sub-grid (~3.5%) and structural (~3%) — are the defensible ones;
  the ATP band is an assumption-driven prior whose size is a free parameter.
- **ATP flat-fallback rate (M2):** the selective fallback flattens 12% of land
  cell-months (tropics 1.4%, temperate 10%, **boreal 22%**); there the point
  estimate is flat and `$var` is the kriging spread (legitimate uncertainty,
  punted point estimate) — documented in §4.4.
- **Envelope-drop transparency (M3):** ~11% of land cell-months have a near-zero
  envelope and are excluded from the strict overshoot medians; with a *floored*
  envelope the PIQS overshoot median is 0.41 vs 0.89 strict — the dropped cells
  have *lower* relative overshoot, so the strict medians do not understate.
- **Locality robust but synthetic (M4):** re-checking the NRT footprint on the
  **full flux (a,b,c)**: PCHIP ≤1, PPM ≤2, minmod ≤1 prior months move >1% under
  a +10% last-month perturbation — locality holds, but this is a *synthetic*
  perturbation; no diff of two real consecutive vNRT releases was done.
- **Minor:** "Steffen ≡ PCHIP" → "matches on max overshoot"; the production ATP
  point-estimate RMS/env is 0.04–0.06 (the flat fallback regressed it from the
  0.003–0.035 prototype); fidelity scripts pin explicit `fit.<method>.rda`.

**Net effect.** The decision is unchanged and better-grounded: PCHIP is local,
sign-definite, closed-form, `(a,b,c)`-native, and *significantly* (if negligibly)
the best on product fidelity; ATP is a *legitimate but assumption-parameterized*
uncertainty option, not a measured one. The overstatements the review flagged —
"statistical tie," "dominant uncertainty," "principled variance" — are corrected
above.

---

## 5. Why PIQS is unsuitable for this product (and why that fix already happened)

PIQS is the citable CT2022 standard and is genuinely smooth, so the argument is
specific — and note the fix is **historical, not pending**: V2 already moved off
PIQS to PCHIP. Two PIQS properties are disqualifying for an NRT product:

1. **Overshoot through zero.** Measured 2026-06-18 on the same 2020 diurnalize:
   up to **~11% of GPP cell-hours wrong-sign** (plant emitting CO₂), and a
   fidelity mean wrecked by ~10¹⁸ overshoot-tail cells. The local sign-definite
   methods remove this: PCHIP 0.1–0.7%, PPM/minmod 0%.
2. **Global solve ⇒ NRT non-reproducibility.** PIQS couples every knot to every
   month, so appending/revising the latest vNRT month (which MiCASA's stream does
   in place) re-solves and changes **all 303 historical months** — the published
   hourly record would silently change each cycle. The local methods confine a
   revision to **0–2 months** (PCHIP 0).

These are the grounds; both are already addressed by the V2 (PCHIP) default.

---

## 6. Recommendation

1. **Keep PCHIP as the default.** It is local, sign-definite, mass-conserving,
   and the **only** method with zero flux jumps (continuity is the smoother's
   raison d'être). Its blemishes — a bounded ≤1.5× bump (shown physical) and a
   ~0.1–0.7% residual sign-flip rate — are minor.
2. **PPM and minmod are selectable alternatives, not improvements.** PPM trades
   PCHIP's small residual for zero overshoot, but **reintroduces month-edge
   discontinuities at ~70% of boundaries** (the steps the smoother exists to
   remove) and is a **fidelity tie** (§4). Use PPM only if strict zero-overshoot
   is required and edge discontinuities are acceptable.
3. **Reject PIQS and MSS.** PIQS: overshoot + global non-locality. MSS:
   overshoot + ~300× slower.
4. **Honest framing.** The consequential change was PIQS→a local method, done at
   V2. PCHIP↔PPM is a refinement-scale tradeoff with no clear winner; the doc
   does **not** claim PPM is better.

All fitters are selectable via `MICASA_FIT_RDA` (`write_pchip.r`, `write_ppm.r`,
`write_linmm.r`, `write_piqs.r`) with no change to the default; cores are
unit-tested (`tests/test_{pchip,ppm,linmm,mss}_fit.r`).

**Field closed (2026-06-18):** the last untested candidate, bounded iterative
mean-preserving (Rymes-Myers), was measured (§2.6) — it is bounded, sign-definite,
smooth, continuous, and NRT-local (≤2 mo), i.e. competitive with PCHIP but not
clearly superior, and costs an iterative solve + a point-value→`(a,b,c)` convert.
No evaluated method dominates PCHIP for this product; PCHIP stays the default,
PPM/minmod/Rymes-Myers are documented selectable alternatives, PIQS/MSS rejected.

---

## 7. Reproducibility

Scripts are tracked in `fitter_diagnostics/`; run from the repo working dir.
The 2020 product comparisons require the shadow diurnalize runs first:

```
# regenerate the 2020 PCHIP/PPM/PIQS shadow products (one sbatch each):
diurn_year=2020 MICASA_MONTH_START=1 MICASA_MONTH_END=12 MICASA_VERSION=v1 \
  MICASA_FIT_RDA=fit.pchip.rda MICASA_DIURN_OUT_DIR=ERA5_2020_pchip Rscript diurnalize-ERA5.r
#   ... =fit.ppm.rda → ERA5_2020_ppm ;  =fit.piqs_v1.rda → ERA5_2020_piqs
```

| Result | Script |
|---|---|
| Overshoot / gradient / boundary-jump bake-off | `fitter_diagnostics/bakeoff_full.r` |
| Integral-preserving linear families (minmod / trapezoidal) | `fitter_diagnostics/linear_compare.r` |
| Smoothness, spurious-extrema, NRT-locality | `fitter_diagnostics/metrics_extra.r` |
| Full-year 2020 budget + aliasing | `fitter_diagnostics/analyze_2020.r` |
| Daily fidelity vs MiCASA daily | `fitter_diagnostics/fidelity_daily.r` |
| Uncertainty / paired / per-biome fidelity | `fitter_diagnostics/uncertainty_fidelity.r` |
| MSS / Steffen / unlimited-parabolic | `fitter_diagnostics/other_methods.r` |
| Unconstrained histospline | `fitter_diagnostics/histospline_check.r` |
| PIQS apples-to-apples score | `fitter_diagnostics/piqs_score.r` |
| Bounded iterative (Rymes-Myers) | `fitter_diagnostics/bounded_iterative.r` |
| Fit-then-aggregate order test (1deg / 4x6) | `fitter_diagnostics/refine_then_average{,_4x6}.r` |
| Out-of-domain survey + local uncertainty bands | `fitter_diagnostics/uncertainty_bands.r` |
| Area-to-point kriging prototype (1-D temporal) | `fitter_diagnostics/atp_kriging.r` |
| Sub-grid (0.1deg) heterogeneity uncertainty | `fitter_diagnostics/subgrid_uncertainty.r` |
| ATP production-fit verification | `fitter_diagnostics/verify_atpk.r` |
| PIQS-with-linear-fallback hybrid | `fitter_diagnostics/piqs_hybrid.r` |
| Product-level fidelity + by-cell bootstrap (C1/M1) | `fitter_diagnostics/fidelity_product.r` |
| ATP variance range-sensitivity (C2) | `fitter_diagnostics/atpk_range_sweep.r` |
| Review follow-ups: fallback rate / env-drop / footprint (M2–M4) | `fitter_diagnostics/critique_followups.r` |
| **Production fitter drivers** | `write_{pchip,ppm,linmm,piqs,mss,atpk}.r` (repo root) |
| Fitter cores + tests | `lib/{pchip,ppm,linmm,mss,atpk}_fit.r`, `tests/test_*_fit.r` |

---

## 8. References

- Bartlein, P. J., *mp-interp* — mean-preserving interpolation reference code
  (Epstein 1991 / Harzallah 1995 / Killworth 1996 / Rymes & Myers 2001 + a
  bounded `enforce_mean`): https://github.com/pjbartlein/mp-interp.
- Boneva, Kendall & Stefanov (1971), *Spline transformations*,
  [JRSS-B 33(1):1–70](https://academic.oup.com/jrsssb/article/33/1/1/7027042).
- Colella & Woodward (1984), *The Piecewise Parabolic Method*, JCP 54(1):174–201, doi:10.1016/0021-9991(84)90143-8.
- Denning, Fung & Randall (1995), *Nature* 376:240–243, doi:10.1038/376240a0.
- Fritsch & Carlson (1980), *Monotone Piecewise Cubic Interpolation*, SIAM JNA 17(2):238–246, doi:10.1137/0717021.
- Gourdji et al. (2010), *ACP* 10:6151–6167, doi:10.5194/acp-10-6151-2010.
- Harten (1983), *High resolution schemes...*, JCP 49(3):357–393, doi:10.1016/0021-9991(83)90136-5.
- JULES land-surface model, *Temporal interpolation* (conservative Sheng & Zwiers
  option; documents curve-fitting overshoot at turning points):
  https://jules-lsm.github.io/latest/input/temporal-interpolation.html.
- Munassar et al. (2025), *ACP* 25:639, doi:10.5194/acp-25-639-2025.
- Olsen & Randerson (2004), *Differences between surface and column atmospheric CO2...*, JGR 109:D02301, doi:10.1029/2003JD003968.
- Rasmussen (1991), *Piecewise integral splines of low degree*, Computers & Geosciences 17(9):1255–1263, doi:10.1016/0098-3004(91)90027-B.
- Rymes & Myers (2001), *Mean preserving algorithm...*, Solar Energy 71(4):225–231, doi:10.1016/S0038-092X(01)00052-4.
- Sheng & Zwiers (1998), *An improved scheme for time-dependent boundary conditions...*, Climate Dynamics 14:609–613, doi:10.1007/s003820050244.
- Steffen (1990), *A simple method for monotonic interpolation...*, A&A 239:443, ADS:[1990A&A...239..443S](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S).
- van Leer (1979), *Towards the ultimate conservative difference scheme V*, JCP 32(1):101–136, doi:10.1016/0021-9991(79)90145-1.
- Wang & Bartlein (2022), *A Fast Mean-Preserving Spline...*, JTECH 39(4):503–512, doi:10.1175/JTECH-D-21-0154.1.
- Weir et al. (2021), *Bias-correcting carbon fluxes...* (LoFI; the inversion
  context these fluxes feed, motivating requirement (3)), *ACP* 21:9609–9628,
  doi:10.5194/acp-21-9609-2021.
- Sax & Steiner (2013), *Temporal Disaggregation of Time Series*, R Journal 5(2):80–88, doi:10.32614/RJ-2013-028.
- Proietti (2006), *Temporal disaggregation by state space methods*, Econometrics J. 9(3):357–372, doi:10.1111/j.1368-423X.2006.00189.x.
- Kyriakidis (2004), *A geostatistical framework for area-to-point spatial interpolation*, Geographical Analysis 36(3):259–289, doi:10.1111/j.1538-4632.2004.tb01135.x.
- Yoo & Kyriakidis (2006), *Area-to-point kriging*, J. Geographical Systems, doi:10.1007/s10109-006-0036-7.
- Da Veiga & Marrel (2012), *Gaussian process modeling with inequality constraints*, Ann. Fac. Sci. Toulouse 21(3):529–555, doi:10.5802/afst.1344.
- Maatouk & Bay (2017), *Gaussian process emulators ... with inequality constraints*, Math. Geosci. 49:557–582, doi:10.1007/s11004-017-9673-2.
- Rizzi, Gampe & Eilers (2015), *Efficient estimation of smooth distributions from coarsely grouped data*, Am. J. Epidemiol. 182(2):138–147, doi:10.1093/aje/kwv020.
- Pya & Wood (2015), *Shape constrained additive models*, Statistics and Computing 25(3):543–559, doi:10.1007/s11222-013-9448-7.
- Hofner, Kneib & Hothorn (2016), *A unified framework of constrained regression*, Statistics and Computing 26, doi:10.1007/s11222-014-9520-y.
