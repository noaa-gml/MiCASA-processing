# Sub-monthly flux reconstruction: method comparison and the case for retiring PIQS

**Status:** decision document · **Date:** 2026-06-18 · **Scope:** the monthly→sub-monthly
flux smoother in `diurnalize-ERA5.r` (the `fit.piqs.rda` coefficients).

This document explains every reconstruction method in (and adjacent to) the
tree, with equations and citations; lays out their pros and cons; presents the
empirical scorecard measured on the real 2001–2026 record; and argues that
**PIQS — the original V1 method — should not be used for this product**, with a
recommendation among the overshoot-free alternatives.

---

## 1. The problem and where the fit sits

MiCASA delivers **monthly** mean NPP and heterotrophic respiration per 1° grid
cell. The hourly NEE product is built with the **Olsen & Randerson (2004)**
scheme, the community standard across the CASA-GFED / CarbonTracker lineage:
decompose NEE into gross fluxes, redistribute GPP within the month with ERA5
shortwave and respiration with a Q₁₀ function of temperature, conserving the
monthly mean:

```
GPP(t)  = GPP_mean · SSRD(t)/SSRD_mean
RE(t)   = RE_mean  · Q10(t)/Q10_mean ,   Q10(t) = 1.5^((T2m−273.15)/10)
NEE(t)  = GPP(t) + RE(t)            (sign convention: positive = source to atm)
```

Used as written, this imposes realistic *intra*-month structure but leaves
**abrupt steps at month boundaries** (each month carries a single flat mean).
To remove those steps, CarbonTracker fits a **mean-preserving smooth curve**
through the monthly series and uses its sub-monthly deviation (`qmod`) in place
of the flat mean. That smoother is the subject of this document. In our code
the chosen method writes per-piece quadratic coefficients `(a,b,c)` to
`fit.piqs.rda`; `diurnalize-ERA5.r` evaluates

```
qmod(t) = a·(t−t_i)² + b·(t−t_i) + c        within month i
```

so **any method that produces `(a,b,c)` per month is a drop-in swap** via
`MICASA_FIT_RDA`.

### The constraint trilemma

For a reconstruction of interval-mean data there is a well-known impossibility
(see [Bartlein's `mp-interp` notes](https://github.com/pjbartlein/mp-interp); the [JULES temporal-interpolation docs](https://jules-lsm.github.io/latest/input/temporal-interpolation.html)):

> **Exact mass-conservation + boundedness/positivity + global smoothness cannot
> all hold simultaneously.** Any mean-preserving *smooth* fit must overshoot
> near sharp turning points; that is precisely the over/undershoot that cancels
> to reproduce the mean.

Every method below is a different choice of *which* property to relax. Three
requirements are non-negotiable for our use:

1. **Mass conservation** — the monthly mean must be preserved exactly (it is the
   only quantity MiCASA actually reports, and the inversion budgets depend on it).
2. **No wrong-sign flux** — GPP must not become a source, RE must not become a
   sink, at any sub-monthly time (unphysical; pollutes the prior).
3. **NRT stability** — revising or appending the most recent month must not
   silently rewrite the historical record (the vNRT stream *is* revised in
   place; see the 2026-06 download-staleness incident).

---

## 2. The methods

All operate per grid cell on the monthly-mean series `m_1…m_n` at knot times
`x_0…x_n` (piece widths `h_i = x_i − x_{i−1}`). The cumulative integral is
`F_k = Σ_{j≤k} m_j h_j`.

### 2.1 Piecewise-constant (mass-conserving baseline)

The flux is the monthly mean itself, `f(t) = m_i` on month `i` (equivalently,
linear interpolation of the *cumulative* `F`). Mass-conserving, never
overshoots, but `qmod ≡ 0` — it contributes **no** sub-monthly structure and
restores a hard step at every boundary. It is the null hypothesis: any fitter
worth using must beat it.

### 2.2 PIQS — Piecewise Integral Quadratic Splines (V1)

**Citation:** Rasmussen, L. A. (1991), *Piecewise integral splines of low
degree*, Computers & Geosciences 17(9):1255–1263,
doi:[10.1016/0098-3004(91)90027-B](https://doi.org/10.1016/0098-3004(91)90027-B).
This is the method documented for **CarbonTracker CT2022**.

Each month gets a quadratic `f_i(t) = a_i(t−t_i)² + b_i(t−t_i) + c_i` subject to:

```
(i)  exact integral:   ∫_{x_{i-1}}^{x_i} f_i dt = m_i · h_i           (mass)
(ii) C0 continuity:    f_i(x_i) = f_{i+1}(x_i)                        (shared knot value)
(iii) smoothness:      remaining DOF fixed by a GLOBAL solve that
                       minimises the summed squared ordinate change,
                       giving a tridiagonal system coupling all knots.
```

Derivative continuity is **not** separately enforced (C⁰, not C¹). The defining
feature — and the source of its problems here — is that step (iii) is a
**global** solve: every knot value depends on every monthly mean.

**Pros:** smooth (minimises a global roughness objective); mass-conserving;
C⁰; the documented CarbonTracker standard, so maximally citable. **Cons:** see
§4 — it overshoots through zero (unphysical sign flips) and its global solve
makes it unstable for an NRT product.

### 2.3 PCHIP-on-cumulative (V2, current default)

**Citation:** Fritsch, F. N. & Carlson, R. E. (1980), *Monotone Piecewise Cubic
Interpolation*, SIAM J. Numer. Anal. 17(2):238–246,
doi:[10.1137/0717021](https://doi.org/10.1137/0717021). (R `splinefun(method="monoH.FC")`,
SciPy `PchipInterpolator`.)

Interpolate the **cumulative** `F` with a Fritsch-Carlson monotone cubic
Hermite, then differentiate analytically — the flux `f = F′` is piecewise
quadratic. The knot slopes `d_k` (= flux at knots) are limited so `F` stays
monotone:

```
secants  S_k = (F_{k+1}−F_k)/h_{k+1} = m_{k+1}
d_k = 0                       if S_{k-1}, S_k have opposite sign (a turning point)
d_k = limited harmonic mean,  |d_k| ≤ 3·min(|S_{k-1}|, |S_k|)      otherwise
```

On segment `k` with `s = (t−x_k)/h_k` and `u_k = m_k`:

```
f(s) = (6s−6s²)·u_k + (3s²−4s+1)·d_k + (3s²−2s)·d_{k+1}
```

Because `F` is monotone, `f = F′` keeps a single sign — GPP stays ≤ 0, RE ≥ 0
**by construction**, with no clipping. The slope rule is **local** (uses only
neighbouring means).

**Pros:** mass-conserving; sign-definite by construction; local; C⁰ flux,
smooth. **Cons:** the flux quadratic can bulge *within* a piece up to **1.5×**
the monthly mean (the Hermite bump `6u·s(1−s)` peaks at `1.5u` when both knot
slopes are zero — a month flanked by near-zero neighbours); a small residual
sign-flip rate survives in near-zero cells (§4).

### 2.4 Integral-preserving linear (MUSCL / slope-limited)

**Citations:** van Leer, B. (1979), *Towards the ultimate conservative
difference scheme V*, J. Comput. Phys. 32(1):101–136,
doi:[10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1);
Harten, A. (1983), *High resolution schemes for hyperbolic conservation laws*,
J. Comput. Phys. 49(3):357–393,
doi:[10.1016/0021-9991(83)90136-5](https://doi.org/10.1016/0021-9991(83)90136-5)
(TVD theory).

There are **two** integral-preserving piecewise-linear fluxes, and the
distinction is decisive:

**(a) Continuous linear** — force the line endpoints to meet at knots. Mass
conservation then *forces* the recursion

```
y_{i+1} = 2·m_i − y_i
```

a one-parameter family (pick the seed; the rest follow). It has a pole at the
Nyquist frequency, so any month-to-month *alternation* in the means — exactly
what low-flux polar cells show — **resonates and diverges**. Measured on the
real record it blows up to ~10⁹× the local envelope and flips sign on 65–72% of
pieces. **Not usable.** (This is the "linear PIQS" considered and rejected in
PROPOSALS #9.)

**(b) Discontinuous, slope-limited (MUSCL/minmod)** — give each month its own
line through its mean with a minmod-limited slope:

```
f_i(t) = m_i + σ_i·(t − t_i^c) ,   t_i^c = month centre
σ_i = minmod(S_{i-1}, S_i) / Δ ,   minmod(a,b) = ½(sgn a + sgn b)·min(|a|,|b|)
```

The line passes through the mean at the centre → **mass-conserving for any
slope**; minmod caps the slope at the smaller neighbouring secant and zeroes it
at turning points → **provably no overshoot (TVD)** and no new extrema. The
price is a **discontinuity at each month edge** (finite-volume reconstructions
are cell-local). Less-diffusive limiters (van Leer, MC, superbee) retain more
slope at the cost of larger edge jumps.

**Pros:** mass-conserving; zero overshoot; local; dead simple. **Cons:**
discontinuous flux at month boundaries; carries the least within-month
structure of the smooth methods.

### 2.5 PPM — Piecewise Parabolic Method (proposed)

**Citation:** Colella, P. & Woodward, P. R. (1984), *The Piecewise Parabolic
Method (PPM) for gas-dynamical simulations*, J. Comput. Phys. 54(1):174–201,
doi:[10.1016/0021-9991(84)90143-8](https://doi.org/10.1016/0021-9991(84)90143-8).

The conservative-reconstruction literature's recommended answer to "mean-
preserving, smoother than linear, but no overshoot." Each month gets a
**parabola** through its mean, with **continuous shared edge values** computed
from monotonized neighbour differences, then a monotonicity limiter:

```
monotonized diff:  δm_i = minmod( (m_{i+1}−m_{i-1})/2 , 2(m_i−m_{i-1}) , 2(m_{i+1}−m_i) )
edge value:        a_{i+½} = m_i + ½(m_{i+1}−m_i) − (δm_{i+1} − δm_i)/6     (continuous)
parabola (ξ∈[0,1]): f(ξ) = a_L + ξ·(a_R − a_L + a_6·(1−ξ)) ,  a_6 = 6(m_i − ½(a_L+a_R))
limiter:           flatten (a_L=a_R=m_i) at a local extremum; otherwise steepen
                   one edge so the parabola has no interior extremum.
```

Mass is exact (the parabola's mean is `m_i` for any edge pair). It is **C⁰
everywhere except where the limiter resets an edge** (i.e. only at genuine
turning points), and it is **local** (uses ±2 neighbours).

**Pros:** mass-conserving; zero overshoot; piecewise-quadratic (fits the
existing `(a,b,c)` storage exactly); mostly continuous; local; in the same
conservative finite-volume lineage the literature endorses. **Cons:** slightly
less within-month structure than PCHIP; tiny edge discontinuities at turning
points (far smaller than minmod's).

### 2.6 Other mean-preserving options (not pursued)

- **MSS — monotone smoothing spline** (in-tree, `write_mss.r`): cubic on `F`
  minimising `∫(F″)²` subject to `F′≥0` at test points, solved per-cell as a QP.
  Measured on the real record (2026-06-18) it is **not** overshoot-free (peak/env median 1.35, max 1.57) and ~24% of land cells carry a wrong-sign GPP knot — the non-negativity constraint binds only at interior test points, not at knots. It is also ~180–370 ms/cell (~100–300× slower than PPM/PCHIP). The QP's banded Hessian does keep it NRT-local (footprint ≤1 month). Not adopted: overshoots like PIQS and far costlier (see §4.1).
- **Bounded iterative mean-preserving**: Rymes & Myers (2001), Solar Energy
  71(4):225–231, doi:[10.1016/S0038-092X(01)00052-4](https://doi.org/10.1016/S0038-092X(01)00052-4);
  Wang & Bartlein (2022), J. Atmos. Oceanic Technol. 39(4):503–512,
  doi:[10.1175/JTECH-D-21-0154.1](https://doi.org/10.1175/JTECH-D-21-0154.1).
  Mean-preserving with explicit lower bounds (clip-then-redistribute). Worth
  knowing as the "force positivity into the fit" family; heavier than PPM and
  not needed once overshoot is gone.
- **Histospline foundations**: [Boneva, Kendall & Stefanov (1971)](https://academic.oup.com/jrsssb/article/33/1/1/7027042), JRSS-B 33(1):1–70. **Conservative remapping analogue**: Jones (1999), MWR
  127(9):2204–2210, doi:[10.1175/1520-0493(1999)127<2204:FASOCR>2.0.CO;2](https://doi.org/10.1175/1520-0493(1999)127<2204:FASOCR>2.0.CO;2).

---

## 3. Why sub-monthly structure matters at all

Resolving sub-monthly/diurnal flux structure is not cosmetic — the
**rectifier effect** (covariance of biosphere flux with boundary-layer mixing)
produces first-order CO₂ signals:

- Denning, Fung & Randall (1995), *Nature* 376:240–243,
  doi:[10.1038/376240a0](https://doi.org/10.1038/376240a0) — the diurnal
  rectifier rivals the seasonal one.
- Munassar et al. (2025), *ACP* 25:639, doi:[10.5194/acp-25-639-2025](https://doi.org/10.5194/acp-25-639-2025)
  — removing diurnal flux structure biases inferred fluxes ~2% globally but up
  to ~48–51% regionally.
- Gourdji et al. (2010), *ACP* 10:6151–6167,
  doi:[10.5194/acp-10-6151-2010](https://doi.org/10.5194/acp-10-6151-2010) —
  inversions must be able to adjust diurnally or incur temporal-aggregation error.

So piecewise-constant (no structure) is genuinely worse; the question is *which*
smooth, mass-preserving, non-overshooting fitter.

---

## 4. Empirical scorecard

Measured on the production fit and a full-year 2020 diurnalize (PCHIP vs PPM vs
minmod), 1° global grid, ~4.4 M land cell-months. Reproducible from the scripts
in `jobs/` (see §6).

| Metric | PIQS (V1) | PCHIP (V2, current) | **PPM (proposed)** | minmod-linear | flat |
|---|---|---|---|---|---|
| Mass-conserving | ✓ | ✓ | ✓ | ✓ | ✓ |
| Overshoot, max peak/envelope | large | 1.50 | **1.00** | 1.00 | 1.00 |
| % pieces overshooting envelope | high | ~19% | **0%** | 0% | 0% |
| GPP sign-flips, 2020 product (cell-hours) | ~6.5% mean, ≤30% cells¹ | 0.1–0.7% | **0%** | 0% | 0% |
| Daily-fidelity RMSE/env, GPP (mean)² | — | 0.151 | **0.149** | 0.159 | 0.181 |
| Daily-fidelity RMSE/env, RESP (mean)² | — | 0.128 | **0.125** | 0.130 | 0.145 |
| Within-month structure (gradient/env, med) | high | 0.111 | 0.060 | 0.031 | 0 |
| Flux value-continuity (jump/env, med) | 0 | **0** | 0.018 | 0.10 | 0.24 |
| Interior extrema (% of months)³ | high | 54–71% | 5–7% | 0% | 0% |
| Annual global NEE budget, 2020 | — | −2.617 PgC | **−2.612 PgC (0.2%)** | ~same | ~same |
| Month-boundary aliasing vs PCHIP⁴ | — | 1.0× | **0.44–1.0×** | worse | worst |
| **NRT-revision footprint** (months rewritten)⁵ | **all 302** | **0** | ≤2 | ≤1 | 0 |
| Fit cost | fast (global solve) | fast | fast | fast | trivial |
| Citable lineage | CT2022 / Rasmussen 1991 | Fritsch-Carlson 1980 | Colella-Woodward 1984 | van Leer 1979 | — |

¹ PIQS sign-flip rates from `verify_v2` Check 3.1 / METHODOLOGY.md and the
original `bakeoff_pchip.py` (AK Tundra 30.9%).
² Reconstruction sampled at day-midpoints vs MiCASA's **own** daily 1° product,
2020, RMS normalised by the local monthly-mean envelope; lower = more faithful.
³ Fraction of cell-months in which the reconstructed flux has a within-month
extremum not implied by the monthly means.
⁴ Power in the 1-month spectral band of the hourly NEE, ratio to PCHIP.
⁵ Perturb the latest monthly mean +10%, refit, count prior months whose flux
moves >1%. PIQS's global solve couples the whole record.

**Two findings that reframe "overshoot":** (1) MiCASA's *own* daily data
exceeds the monthly-mean envelope routinely — daily-max/env median ≈ 1.0, 90th
≈ 1.3–1.4 — so a sub-monthly peak above the envelope is *physically real*;
bounded overshoot is not the enemy, **wrong-sign and unbounded oscillatory
overshoot is.** (2) On fidelity to that daily truth the ranking is consistent:
**PPM < PCHIP < minmod < flat** — i.e. PCHIP's extra wiggle (interior extrema in
54–71% of months) does *not* buy better fidelity; PPM's disciplined shape
matches the daily data marginally better. (These fit-level differences are
small because ERA5 supplies the within-month weather in the real pipeline.)

### 4.1 Other methods evaluated (2026-06-18)

Four further candidates were scored on the same record so the field is complete; **none beats PPM**:

| Method | overshoot peak/env (med / max) | wrong-sign cells | cost | verdict |
|---|---|---|---|---|
| **MSS** (in-tree QP, `write_mss.r`) | 1.35 / 1.57 | ~24% | ~180–370 ms/cell | overshoots *and* slow — reject |
| **Steffen-on-cumulative** (1990) | 0.83 / 1.50 | — | fast | identical to PCHIP (same monotone-cubic bump cap) — no gain |
| **Unlimited parabolic** (PPM, limiter off) | 0.84 / 1.25 | bounded | fast | fully continuous but overshoots 1.25× — dominated by PPM |
| **van Leer / MC / superbee** linear | 1.00 / 1.00 | 0 | fast | gradient & jumps sit between minmod and PPM — dominated by PPM |

Takeaways: (1) MSS's `F′≥0` test-point constraint does not bind at knots, so it overshoots and flips sign like PIQS — the in-tree “alternative” is not a fix. (2) A different monotone-cubic slope rule ([Steffen 1990, A&A 239:443](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S)) yields the *same* 1.5× cap as Fritsch-Carlson — PCHIP's bump is intrinsic to monotone-cubic-on-cumulative, not a quirk of the FC rule. (3) PPM's limiter does real work: it removes the 1.25× overshoot the unlimited parabolic still shows.

---

## 5. The case for retiring PIQS

PIQS is the documented CarbonTracker standard and is genuinely smooth, so the
argument must be specific. Two properties make it the **wrong choice for this
NRT product**:

### 5.1 It overshoots through zero — unphysical sign flips

PIQS prioritises global smoothness over boundedness, so in high-seasonality
cells its quadratic dips across zero: **up to ~30% of cell-hours show positive
GPP (a plant emitting CO₂) or negative respiration** in boreal/tundra regions
(Check 3.1; AK Tundra 30.9%). This is not a tuning artefact — it is the
trilemma (§1) biting. The pipeline's polar-night `GPP=0` clip (PROPOSALS #8) was
a band-aid for exactly this. The local sign-definite methods remove the disease,
not the symptom: PCHIP drops it to 0.1–0.7%, and **PPM and minmod to exactly
0%** (verified across all 8.4 M cell-hours of 2020).

### 5.2 Its global solve is unstable for an NRT product

This is the decisive, less-obvious point. PIQS fixes its free degrees of freedom
with a **global** tridiagonal solve, so **every knot value depends on every
monthly mean.** When the latest vNRT month is appended *or revised in place*
(which MiCASA's vNRT stream routinely does — see the 2026-06 staleness incident
where NASA had silently re-issued April), re-running PIQS changes the
coefficients for **all 302 historical months**. The published hourly product
would therefore be **non-reproducible**: every NRT cycle rewrites the entire
past record by small amounts, defeating provenance and perturbing an already-
assimilated prior.

The local methods confine a revision's footprint to **0–2 months** (PCHIP: 0;
minmod: ≤1; PPM: ≤2). For an operationally-republished NRT product consumed by
an atmospheric inversion, locality is not a nicety — it is correctness.

### 5.3 What PIQS is *not* penalised for

In fairness: PIQS conserves mass, is C⁰, is fast, and its smoothness objective
is principled. If MiCASA were a one-shot, final, whole-record product with no
near-zero cells, PIQS would be defensible. It is neither — it is NRT and global,
and §5.1–5.2 are disqualifying in that setting.

---

## 6. Recommendation

1. **Do not use PIQS** for this product (overshoot + global non-locality).
2. **PPM is the recommended fitter.** It is the only option that is
   simultaneously mass-conserving, zero-overshoot, zero-sign-flip, smooth,
   local (≤2-month NRT footprint), best-on-daily-fidelity, and in the
   conservative finite-volume lineage the literature endorses. It fits the
   existing `(a,b,c)` storage, so `diurnalize-ERA5.r` consumes it unchanged.
3. **PCHIP (current V2) is an acceptable status quo.** Its only blemishes are a
   bounded 1.5× bump (shown to be physically plausible against the daily data)
   and a ~0.1–0.7% residual sign-flip rate. Switching PCHIP→PPM is a refinement,
   not a fix; switching PIQS→anything-local is a fix.
4. **The integral-preserving linear that was requested (minmod) works** — it is
   honest and overshoot-free — **but it is dominated by PPM** on fidelity,
   structure, and continuity for the same guarantees, and the *continuous*
   linear variant is numerically unstable (§2.4a).

All four candidate fitters are selectable today via `MICASA_FIT_RDA`
(`write_pchip.r`, `write_ppm.r`, `write_linmm.r`, `write_piqs.r`), with no change
to the default; the cores are unit-tested (`tests/test_{pchip,ppm,linmm,mss}_fit.r`).

---

## 7. Reproducibility

| Result | Script |
|---|---|
| Overshoot / gradient / boundary-jump bake-off (all methods) | `jobs/bakeoff_full.r` |
| Integral-preserving linear families (minmod / trapezoidal) | `jobs/linear_compare.r` |
| Smoothness, spurious-extrema, NRT-locality | `jobs/metrics_extra.r` |
| Full-year 2020 product comparison (budget, aliasing) | `jobs/analyze_2020.r` |
| Fidelity vs MiCASA daily | `jobs/fidelity_daily.r` |
| Other methods (MSS, Steffen, unlimited parabolic) | `jobs/other_methods.r` |
| Fitter cores + tests | `lib/{pchip,ppm,linmm,mss}_fit.r`, `tests/test_*_fit.r` |

---

## 8. References

- Bartlein, P. J., *mp-interp* — mean-preserving interpolation reference code (Epstein 1991 / Harzallah 1995 / Killworth 1996 / Rymes & Myers 2001 methods + a bounded `enforce_mean`): https://github.com/pjbartlein/mp-interp.
- Boneva, Kendall & Stefanov (1971), *Spline transformations*, [JRSS-B 33(1):1–70](https://academic.oup.com/jrsssb/article/33/1/1/7027042).
- Colella & Woodward (1984), *The Piecewise Parabolic Method*, JCP 54(1):174–201, doi:10.1016/0021-9991(84)90143-8.
- Denning, Fung & Randall (1995), *Nature* 376:240–243, doi:10.1038/376240a0.
- Fritsch & Carlson (1980), *Monotone Piecewise Cubic Interpolation*, SIAM JNA 17(2):238–246, doi:10.1137/0717021.
- Gourdji et al. (2010), *ACP* 10:6151–6167, doi:10.5194/acp-10-6151-2010.
- Harten (1983), *High resolution schemes...*, JCP 49(3):357–393, doi:10.1016/0021-9991(83)90136-5.
- Jones (1999), *Conservative remapping...*, MWR 127(9):2204–2210, doi:10.1175/1520-0493(1999)127<2204:FASOCR>2.0.CO;2.
- JULES land-surface model, *Temporal interpolation* (conservative Sheng & Zwiers 1998 option; documents curve-fitting overshoot at turning points): https://jules-lsm.github.io/latest/input/temporal-interpolation.html.
- Munassar et al. (2025), *ACP* 25:639, doi:10.5194/acp-25-639-2025.
- Olsen & Randerson (2004), *Differences between surface and column atmospheric CO2...*, JGR 109:D02301, doi:10.1029/2003JD003968.
- Rasmussen (1991), *Piecewise integral splines of low degree*, Computers & Geosciences 17(9):1255–1263, doi:10.1016/0098-3004(91)90027-B.
- Rymes & Myers (2001), *Mean preserving algorithm...*, Solar Energy 71(4):225–231, doi:10.1016/S0038-092X(01)00052-4.
- Steffen (1990), *A simple method for monotonic interpolation...*, A&A 239:443, ADS:[1990A&A...239..443S](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S).
- van Leer (1979), *Towards the ultimate conservative difference scheme V*, JCP 32(1):101–136, doi:10.1016/0021-9991(79)90145-1.
- Wang & Bartlein (2022), *A Fast Mean-Preserving Spline...*, JTECH 39(4):503–512, doi:10.1175/JTECH-D-21-0154.1.
- Weir et al. (2021), *Bias-correcting carbon fluxes...*, *ACP* 21:9609–9628, doi:10.5194/acp-21-9609-2021.
- CarbonTracker CT2022 documentation, NOAA GML: https://gml.noaa.gov/ccgg/carbontracker/CT2022/documentation.php
