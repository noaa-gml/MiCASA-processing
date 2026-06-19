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
are a **statistical tie** on fidelity, and PCHIP wins on continuity.

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
`f` keeps one sign by construction; the slope rule is **local**. **Pro:**
mass-conserving, sign-definite, **C⁰-continuous (zero flux jumps)**, local,
smooth. **Con:** a bounded ≤1.5× within-piece bump (the `6u·s(1−s)` Hermite term
peaks at `1.5u` when both knot slopes vanish — a month flanked by near-zero
neighbours) and a small residual sign-flip rate in near-zero cells (§4).

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
  ms/cell (~100–300× PPM/PCHIP). Banded QP ⇒ NRT-local (≤1 mo). Rejected.
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
| **NRT-revision footprint** (months rewritten)³ | **all 302** | **0** | ≤2 | ≤1 | 0 |
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
- **PCHIP vs PPM fidelity is a statistical tie.** Paired, same-cell test
  (`uncertainty_fidelity.r`): median Δ(PCHIP−PPM) = **0.0006**, IQR **[−0.006,
  +0.009]** (straddles zero), PPM better in **54%** of cell-months — and only
  **49%** in boreal/polar (PCHIP better there). The 0.149-vs-0.151 mean gap is
  within noise. PPM does **not** measurably beat PCHIP on fidelity.

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
   in place) re-solves and changes **all 302 historical months** — the published
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
| Fitter cores + tests | `lib/{pchip,ppm,linmm,mss}_fit.r`, `tests/test_*_fit.r` |

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
