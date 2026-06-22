# Switching MiCASA V1 → V2: the case for adopting V2

**Status:** decision / review document · **Date:** 2026-06-21 · **Purpose:** make the
case that it is worth switching the MiCASA pipeline feeding CarbonTracker from **V1**
— the long-running, verified production line — to **V2** (`main`, tagged `v2.1.0`).
The default position is "stay on V1": it works
and has years of track record, so the burden is on V2 to earn the change. Every
change from V1 is justified below with its rationale, quantified impact, and
verification, so the case can be audited change by change.

> **Decision requested:** adopt **V2** as the MiCASA pipeline feeding CarbonTracker,
> replacing V1. V2 is built, tagged (`v2.1.0`), and verified; everything below is the
> case that the switch is worth it, with the receipts to check it. *(V2 also flips the
> respiration driver default to **soil** temperature, §2 — a small (+2.3% NEE diurnal
> amplitude), validated change, reversible to the bit; independent of the fitter switch.)*

**This document is self-contained** — the
load-bearing scorecard, equations, figures, and references are inlined here; the
repo docs (`FITTER_COMPARISON.md`, `DIURNALIZATION_ALTERNATIVES.md`,
`METHODOLOGY.md`, `PROPOSALS.md`, `CHANGELOG.md`) hold the fuller bake-offs and
dated logs but are not needed to follow the argument below.

## Executive summary

**The case in brief — V1 is proven; here is why V2 is worth the switch:**

- **What you keep.** V1's verified climate-signal history transfers — and I
  **measured it**, not just argued it: a direct diff of the full PIQS and PCHIP
  products (§0) preserves the long-term **trend to Δ 2×10⁻⁵ PgC/yr/yr** and the
  ENSO/COVID anomalies to **< 0.001 PgC** (absolute annual budget within ≤0.5%), and
  the 0.1°→1° aggregation fix shifts **no** band-level mass (0.04%, §20.1). Every
  fitter is integral-preserving (§0); the non-fitter changes are proven no-ops (§4).
- **What you gain.** V2 fixes three genuine V1 defects — PIQS's unphysical wrong-sign
  sub-monthly fluxes (§1), PIQS's global solve that **rewrites V1's published record
  on every NRT update** (§1), and a latitude-weight **bug** in V1's 0.1°→1°
  aggregation (§3.1) — plus operational hardening (provenance, run manifest, ERA5
  FastTrack) and a 60-check + 153-test verification base V1 never had (§4, §6).
- **What it risks.** Essentially only revalidation effort — mitigated by the evidence
  this document is built on: bit-identical proofs for the no-ops (§4), the
  budget-invariance guarantee (§0), and the verification base (§6). V2 *departs* from
  V1's history in two deliberate places: **correcting** the 0.1°→1° aggregation bug
  (<0.01% typical, §3.1) and **defaulting the respiration driver to soil** temperature
  (+2.3% NEE diurnal amplitude, §2) — both validated, and the latter reversible to the
  bit (`MICASA_RESP_DRIVER=airtemp`). Neither adds open-ended risk.

**Bottom line: V2 preserves everything V1 got right and corrects the little it got
wrong — and the cost of switching is quantified here, not asserted.**

The rest of this summary follows the document's structure, one bullet per section:

- **Framing (§0).** Every V1→V2 change is either a **proven no-op (B)** —
  bit-identical / exact-equivalent — or a **quantified improvement (A)** carrying a
  measured impact + a physical or statistical justification + a guarding check. The
  change register above lists them all.
- **§1 — Fitter PIQS → PCHIP fixes two real PIQS defects, at zero risk to the
  long-term signal.** PIQS produced unphysical sub-monthly fluxes (wrong-sign
  overshoot — GPP as a *source* in **~11%** of GPP cell-hours) and, as a single global
  solve, **rewrote the entire ~303-month published record on every near-real-time
  (NRT) revision**. PCHIP fixes both (overshoot **→ ≤0.9%**; NRT footprint **→ 0**).
  This matters because the **sub-monthly shape is what the inversion ingests**; it is
  *safe* because every fitter is integral-preserving, so the monthly-and-longer budget
  — annual totals, trend, IAV, ENSO/COVID — is **identical by construction** (§0).
- **§2 — Diurnalization framework is V1's, unchanged; V2's one diurnal change is to
  default the respiration driver to soil temperature.** Driving respiration off soil
  rather than 2-m air is physically motivated and free (`stl1` already loaded), and damps
  the imposed respiration diurnal cycle (amplitude ratio **0.80**, boreal **0.61**;
  NEE-level effect small, **+2.3%**). The eddy-covariance gate, run properly, splits
  into two questions: soil is the better **seasonal** respiration driver (12/13
  AmeriFlux sites, p=0.003), and the **within-day** relationship the diurnalization
  actually sets is a **tie** (R²≈0.003 for both, below the EC noise floor). So the flip
  has no measured within-day downside and is seasonally + mechanistically correct.
  The legacy air path stays selectable (`MICASA_RESP_DRIVER=airtemp`), reversible to the
  bit (byte-identical to the full legacy product together with `MICASA_POLAR_CLIP=plain`,
  §3.2). Independent of the V1→V2 fitter switch.
- **§3 — The other number-moving changes are small and correct:** the 0.1°→1°
  aggregation latitude-weight bug fix (V1 mis-weighted; <0.01% typical, larger toward
  poles), the polar-night clip now **mass-conserving by default** (closes the ~0.16%
  median high-latitude GPP gap; `MICASA_POLAR_CLIP=plain` for the legacy zero-clip),
  the ERA5 FastTrack fallback (NRT
  trailing months only), and per-month climatology auto-detect.
- **§4 — Everything else is a proven no-op:** 15 library / refactor / compression /
  provenance / manifest / packaging changes, each verified bit-identical or
  exact-equivalent and CI-guarded.
- **§5 — Considered and rejected:** subtracting MiCASA's **ATMC** correction (it
  would double-dip the inversion's own atmospheric constraint — data leakage; and it
  shifts the prior's mean sink −2.45→−5.99 PgC/yr and flips its trend), and
  **"PIQS-then-revert-to-linear"** (the stakeholder-preferred alternative — dominated:
  doesn't fix non-locality, injects a **0.97 mol m⁻² s⁻¹** discontinuity vs PCHIP's
  exact **0**, only *ties* on daily fidelity).
- **§6 — Standing verification base:** `verify_v2` — 60 checks, committed clean at
  **54 PASS / 8 INFO / 0 FAIL / 0 WARN**, re-run on the shipped soil product
  (`verify_v2_summary_soil_20260622.txt`; identical tally to the pre-flip
  `verify_v2_summary_20260621.txt`): every
  product / science / provenance check passes, including the §20 cross-product checks
  and the §24 run-manifest audit. Plus `tests/` (153 on Orion; 143 reproduced green
  locally, 10 `quadprog`-gated). V1 has no comparable harness.

## Change register — the whole scope at a glance

Every V1→V2 change, classified **(A)** moves output numbers or **(B)** proven no-op.
The two changes that move the default product beyond V1's record — the fitter swap and
the respiration-driver flip — are bold.

| # | Change | Type | Impact / status | Where |
|---|---|---|---|---|
| 1 | Fitter **PIQS → PCHIP** | A | budget unchanged (integral-preserving); sub-monthly sign-flips ~11% → ≤0.9%; NRT-local | §1 |
| 2 | **Respiration driver: 2-m air → soil temp** | A (shipped) | **default flipped to soil** — soil better *seasonally* (12/13, p=0.003); within-day (what diurnalization sets) a tie (R²≈0.003 both) so no downside; +2.3% NEE diurnal amp; airtemp selectable & byte-identical | §2 |
| 3.1 | 0.1°→1° aggregation latitude-weight bug fix | A | V1 mis-weighted; <0.01% typical, larger toward poles | §3.1 |
| 3.2 | Polar-night clip — **mass-conserving by default** | A | closes the ~0.16% median high-latitude GPP gap (redistributes onto lit hours); legacy zero-clip = `MICASA_POLAR_CLIP=plain` | §3.2 |
| 3.3 | ERA5 FastTrack dual-tree fallback | A | NRT trailing months only; per-day provenance | §3.3 |
| 3.4 | Per-month climatology auto-detect | A | real-vs-climatology decided per month by file presence | §3.4 |
| — | 15 library / refactor / compression / provenance / manifest / packaging changes | B | proven no-ops (bit-identical or exact-equivalent) | §4 |
| — | ATMC subtraction; PIQS+linear-fallback; MSS; PPM-as-default; … | rejected | diligence, not changes | §5 |

**How this document is organized.** The executive summary and the register above are
the whole story at a glance; the rest is for auditing. **§§1–3** detail each change
that moves numbers (including the respiration-driver default flip, §2); **§4** the proven
no-ops; **§5** what was considered and rejected; **§6–7** the evidence and
limitations; **Appendix A** the fitter equations. Skip the equations and the full
scorecard (§1) if you trust the headline.

---

## 0. How to read this — two categories and one invariant

Every change is exactly one of:

- **(B) Behavior-preserving** — the output product is *provably unchanged*:
  bit-identical (`ncdiff` max |Δ| = 0), an exact-equivalent algorithm (matched to
  floating point), or purely additive metadata / observability / packaging.
  These need no *scientific* defense because they move no number the science
  sees; the defense is the verification that proves it (§4).
- **(A) Intentional improvement** — output numbers *do* change. Each carries:
  what V1 did and why it was worse, what V2 does, the quantified impact, why the
  new behavior is correct (physics / statistics / citation), the verification
  guarding it, and the residual risk (§1–§3).

**The master invariant (why the headline change is safe).** Every fitter in the
tree is **integral-preserving by construction**: each month's per-piece integral
equals that month's MiCASA mean. Therefore *the monthly-and-longer carbon budget
— annual totals, the long-term trend, interannual variability, the ENSO and
COVID signals — is identical across any two fitters* **at the fit level**. Two
scope caveats:
- This is a property of the *fitter*. The **shipped** hourly NEE additionally
  applies a polar-night clip (§3.2). Under the legacy *plain* zero-clip this opened a
  small **GPP monthly-mean gap** at high latitudes (Check 2.2: ~0.16% median cell-month,
  ~2% p99) — a property of the clip, not the fitter. **V2's default conserve clip closes
  it**: it redistributes the clipped dark-hour uptake onto each cell's lit hours,
  restoring the monthly mean exactly. So "mass-preserving" is exact for the fit and, on
  the conserve default, for the delivered GPP too (by construction). (`MICASA_POLAR_CLIP=plain`
  reverts to the legacy gap; the production record is re-diurnalized onto the conserve
  default for v2.2.0 — see CHANGELOG.)
- Mass-preservation is exact; **sign-preservation is not** (see §1, claim 2) —
  the two are independent.

The fitter can only change the **sub-monthly shape**. The budget invariance rests on
two legs:
- **Deductive:** verify_v2 Check 2.1 confirms each fitter preserves the per-piece
  integral (max-abs < 1e-9, max-rel < 1e-6). Equal monthly means ⇒ equal
  monthly-and-longer budget, so the PIQS↔PCHIP equality of the trend / ENSO / COVID
  signals follows by construction.
- **Measured — diffed both ways, directly.** I ran the full PIQS product
  (`MiCASA_v1_piqs`) and the shipped PCHIP product (`MiCASA_v2`) through the *same*
  global-annual-NEE computation, 2001–2024
  (`fitter_diagnostics/piqs_vs_pchip_section15.py`). The published climate signals
  match to high precision: **trend +0.0299 (PCHIP) vs +0.0299 (PIQS), Δ = 2×10⁻⁵
  PgC/yr/yr**; 2015-16 El Niño anomaly **+0.091 vs +0.090**; 2020 COVID **−0.042 vs
  −0.041** (Δ < 0.001 PgC). The absolute annual NEE agrees to **≤0.5% (worst year),
  ~0.4% in the mean** — that residual was the *plain* polar-night clip (§3.2) removing
  slightly different dark-hour GPP from the two fitters' different sub-monthly shapes,
  *not* a fit-level budget difference (and the V2 default conserve clip removes it
  entirely, since it restores each monthly mean exactly). (The +0.0299 trend over this window matches verify_v2
  Check 16.2's independently-computed v1-only slope, validating the calculation; the
  full-record headline trend is +0.0447, §15.1.)

So switching the fitter leaves the **published climate signal intact to <0.1%** and
the absolute annual budget to ≤0.5% — it changes the sub-monthly shape, removing most
of an unphysical artifact. This is now a *measurement*, not an argument.

**Evidence base.** Two independent harnesses back every claim below — `verify_v2`
(60 checks) and `tests/` (153) — with the full tallies, the §24 manifest-artifact
explanation, and the per-change → check map in **§6**.

---

## 1. Headline change — fitter PIQS → PCHIP (A)

**Why switch at all, if the budget is unchanged?** Because the long-term
budget-invariance (§0) does *not* protect against PIQS's two real defects — both
live in the **sub-monthly shape**, which is precisely what an atmospheric inversion
ingests (resolving sub-monthly / diurnal structure is first-order for inversions —
the rectifier effect; `FITTER_COMPARISON.md` §3). PIQS **(a)** produces unphysical
**wrong-sign sub-monthly fluxes** (GPP appearing as a *source*) and **(b)** — being a
single global solve — **rewrites the entire published record on every NRT update**.
PCHIP fixes both. The switch is *safe* exactly because every fitter is
integral-preserving: it cannot move the annual / trend / ENSO budget (§0). So the net
is **a real improvement where it matters (the sub-monthly product + NRT stability),
provably harmless where it doesn't (the long-term signal)** — which is why "it can't
move the science signal" is the reassurance, not the reason.

This change prompted the V1↔V2 concern, so it gets the fullest defense — the decisive
comparison is the scorecard and constraint trilemma below. (The fuller per-method
bake-off, incl. PPM / minmod / MSS / ATP-kriging, is in `FITTER_COMPARISON.md`, not
needed to follow this section.)

**V1 — PIQS** (Piecewise Integral Quadratic Splines, Rasmussen 1991;
CT2022-documented). Per-cell quadratic pieces, each preserving the monthly
integral, C⁰ at knots. Two disqualifying problems for an NRT product:
1. **Overshoot → unphysical sign flips.** In sharply seasonal cells the quadratic
   overshoots through zero, producing positive (source) GPP and negative
   respiration sub-monthly. Measured rate on a regenerated PIQS fit
   (`fitter_diagnostics/piqs_score.r`, full 2001–2026 record): **6.55%** of GPP
   cell-hours mean, **14.70%** max (the per-month max is ~11% on the 2020 product);
   Rh 0.122% / 0.444%. (The PCHIP product's own Check 3.1 is the 0.11%/0.94% below.)
2. **Non-locality.** PIQS is a single global solve over the entire record, so any
   NRT revision **rewrites the entire ~303-month record** — the published past
   changes every cycle.

**V2 — PCHIP-on-cumulative** (Fritsch & Carlson 1980; `splinefun(method="monoH.FC")`
/ `scipy PchipInterpolator`). Monotone-cubic Hermite interpolation of the
cumulative integral F(t), differentiated analytically to the flux f = F′ as a
piecewise quadratic (same `(a,b,c)` storage as PIQS).

Both store the same per-piece quadratic `(a,b,c)` with the same mass-preservation
identity; they differ only in how the free coefficient is set. **In one line:** PIQS
fixes it by a *global* C⁰ continuity solve (couples every month → non-local,
sign-unconstrained → overshoots in sharply seasonal cells); PCHIP sets it from
*local* Fritsch-Carlson knot slopes (~1-month revision footprint, sign-definite *at
the knots* by the monotonicity limiter). Both reproduce the identical monthly mean
`ȳᵢ` — hence the budget-invariance (§0). Full equations in **Appendix A**.

### The constraint trilemma — why every fitter relaxes *something*

For interval-mean reconstruction there is a **well-established practical tension**
— documented in Bartlein's *mp-interp* notes, the JULES temporal-interpolation
docs, and the smoothing-spline literature, though **not as a single named
theorem** — between **exact mass, strict no-overshoot, and global smoothness**: a
mean-preserving *smooth* fit tends to overshoot near sharp turning points, and
forcing strict no-overshoot tends to break smoothness or continuity there. It is a
strong empirical regularity, not a proof, and the axis that gives is *strict*
boundedness: tolerating a small **bounded** overshoot (no sign flip) buys back
mass, continuity, smoothness, *and* locality together. Every method keeps mass and
relaxes one axis — *which* one is the whole argument:

| Method | mass | relaxes | keeps |
|---|---|---|---|
| **PIQS (V1)** | ✓ | **boundedness** — overshoots, incl. wrong-sign, *unbounded* | global smoothness, C⁰ |
| **PCHIP (V2)** | ✓ | strict boundedness → a **bounded ≤1.5×** bump | C⁰ flux, sign-definite at knots, **local** |
| Rymes–Myers (bounded-iterative) | ✓ | strict boundedness → bounded ≤1.45× bump (**same axis as PCHIP**) | sign-**definite**, C⁰, smooth, **local** — competitive (`FITTER_COMPARISON.md` §2.6); not default (iterative, no closed form, emits point values not native `(a,b,c)`) |
| PPM | ✓ | **global continuity** (small jumps at ~70% of edges) | no overshoot, smooth |
| minmod-linear | ✓ | **continuity + curvature** | no overshoot |
| PIQS + linear-fallback | ✓ | **continuity** (large jumps where it patches) — **and keeps PIQS's non-locality** | no overshoot |

Three requirements are **non-negotiable** for a CO₂-inversion NRT prior:
**(1) mass conservation, (2) no wrong-sign flux** (GPP must not be a source),
**(3) NRT stability** (a revised recent month must not rewrite the published
record). PIQS fails (2) and (3); the linear-fallback fails (3) and breaks
continuity hard (§5.1); PCHIP meets all three, at the cost of only a bounded,
physically-real sub-monthly bump (MiCASA's own daily data exceeds the monthly-mean
envelope routinely, so a peak above it is real, not an artifact). That is the
whole fitter case in one table. **PCHIP is not the *only* fit that reaches this
corner:** the bounded-iterative Rymes–Myers scheme keeps mass, sign, continuity,
smoothness *and* locality too — relaxing the same strict-boundedness axis as PCHIP
(`FITTER_COMPARISON.md` §2.6, measured 2026-06-18). PCHIP is the default over it on **engineering** grounds, not a
trilemma one: PCHIP is closed-form and emits the native `(a,b,c)` quadratic
directly, whereas Rymes–Myers is iterative (no closed form, a `niter` to tune) and
produces point values that would need a fit/convert step in `diurnalize`. So the
table makes a *decisive* case against PIQS and the linear fallback, and an
implementation-broken tie for PCHIP over Rymes–Myers.

**Claims, stated to their exact scope:**
1. **Budget-invariant at the fit level** — monthly+ means identical by
   construction (each piece's integral = the monthly mean; mass-preserving). The
   climate signal is therefore fitter-invariant (master invariant above; Check
   2.1, Section 15). Globals unchanged: GPP ∈ [−126.2, −119.8], resp ∈
   [117.0, 123.9] PgC/yr (Check 5.1; reproduced 2026-06-21). The PCHIP↔PIQS
   equality of these globals is by the integral-preserving invariant (§0), argued
   not separately diffed.
2. **A large sub-monthly improvement — a ~16–60× reduction in sign flips, *not*
   elimination by construction.** PCHIP fits a Fritsch-Carlson *monotone* cubic to
   the cumulative integral, so the flux f = F′ is sign-definite **at the knots**
   and overwhelmingly so in the interiors. It is **not** sign-definite *everywhere*
   by construction: Fritsch-Carlson constrains the cubic's *knot* slopes, and the
   derivative quadratic can still dip mid-segment even on strictly single-signed
   input. I reproduced this — worst interior flux **−0.042 on strictly positive
   monthly means** (0.1% of 20,000 synthetic series carry any wrong-sign dip), see
   [`fitter_diagnostics/pchip_sign_definiteness.r`](../fitter_diagnostics/pchip_sign_definiteness.r)
   and its committed output `pchip_sign_definiteness_20260621.txt`.
   What PCHIP buys is a 1–2 order-of-magnitude *reduction* vs PIQS, leaving a small
   bounded residual: GPP **6.55% → 0.11%** mean (~60×), 14.70% → **0.94%** max
   (16×); Rh 0.122% → 0.0000% mean, 0.444% → 0.002% max (Check 3.1,
   `verify_v2_summary_20260621.txt`). Check 18.2 confirms **C⁰ flux continuity**
   (flux-value |jump| ≤ 1e-12 at knots — i.e. C¹ of the cumulative F; PCHIP is *not*
   C¹ in the flux, see §5.1). Check 18.1 (INFO)
   finds **0.646% of GPP *segments*** carry a wrong-sign interior
   point (max 1.24e-6) — a *different* denominator (segments, not cell-hours), so
   consistent in order of magnitude rather than a strict cross-check, but it
   likewise confirms the residual is real and that "sign-definite everywhere" would
   be false.
3. **Reduction is rule-based, not tuned** — the knot-level sign-definiteness and
   the interior reduction come from the Fritsch-Carlson monotonicity rule, not a
   fitted parameter; the small residual interior dips (and any dark-hour GPP) are
   then removed by the polar-night clip (§3.2). PIQS's overshoot, by contrast, was
   an order of magnitude larger and *not* removable without a clip that would
   distort the bulk flux.
4. **NRT-local** — Fritsch-Carlson slopes use only neighbouring monthly means, so
   a revision's footprint is ~1 month, vs PIQS rewriting the whole record. This is
   a correctness requirement for a published NRT product, independent of the
   physics. (Locality follows from the slope formula; it is argued, not separately
   diff-tested.)

**Verification:** Checks 2.1, 3.1, 6.1, 18.1, 18.2; `tests/test_pchip_fit.r`
(12 checks, green); `bakeoff_pchip.py` (6 biome cells: 0% flips *on those cells*
vs PIQS up to 30.91% — the full-grid residual is the ≤0.94% in claim 2, |Δ flux| < 2e-11).

### Empirical scorecard — production fit + full-year 2020 diurnalize (~4.4 M land cell-months)

The trilemma table above is conceptual — *which property* each method sacrifices.
This scorecard is the **measured** version of the same comparison on the production
fit. Throughout, **env** = the local monthly-mean flux magnitude (the natural scale
for normalizing a sub-monthly excursion).

| Metric | PIQS (V1) | **PCHIP (V2)** | PPM | minmod | PIQS+lin-fallback |
|---|---|---|---|---|---|
| Mass-conserving | ✓ | ✓ | ✓ | ✓ | ✓ |
| Overshoot peak/env (med / max) | 0.93 / **~10¹⁸** (diverges, see ²) | 0.83 / **1.50** | 0.78 / 1.00 | – / 1.00 | – / 1.00 |
| GPP wrong-sign (cell-hours, max month) | **~11%** (2020) | **0.1–0.9%** | 0% | 0% | 0% |
| Daily-fidelity RMSE/env, GPP (**median**) ² | 0.086 | **0.081** | 0.079 | 0.094 | 0.079 |
| Flux continuity (jump/env med ; % edges) | C⁰ | **0 ; 0%** | 0.018 ; ~70% | 0.10 ; ~93% | **0.25⁴ ; 29%** |
| **NRT footprint** (months rewritten, +10% revision)³ | **all 303** | **0** | ≤2 | ≤1 | **all 303** |
| Lineage | Rasmussen 1991 | Fritsch-Carlson 1980 | Colella-Woodward 1984 | van Leer 1979 | — |

² I report the **median** RMSE/env because the *mean* is tail-sensitive: PIQS's
GPP mean is **18.6**, wrecked by cells where the global solve diverges to ~10¹⁸× the
envelope (28% of GPP cell-months carry a wrong-sign knot), and even the local
methods' means (PCHIP 0.151 / PPM 0.149 / minmod 0.159; committed in
`fitter_diagnostics/uncertainty_fidelity_20260621.txt`) are noisier than their
medians. On the robust median all local methods sit within ~0.015, and **PIQS's own
median (0.086) is fine** — its disqualifiers are the overshoot tail and
non-locality, not median fidelity. PIQS numbers measured 2026-06-18 on a regenerated
fit, same record/diurnalize (`piqs_score.r`). ³ Perturb the latest monthly mean +10%,
refit, count prior months moving > 1% (PIQS's global solve couples the whole
record). ⁴ The PIQS+lin-fallback continuity entry is the **finite-envelope median**
jump/env (0.25) over the 29% of cell-months that trigger the patch; 38% of patched
edges fall in near-zero-envelope transition months where jump/env is undefined, so
the honest cross-method statement is the **absolute** discontinuity budget — 0.97
mol m⁻² s⁻¹ vs PCHIP's exact 0 (§5.1). **Among the columns shown, PCHIP is the only one that is sign-safe,
C⁰-continuous, *and* NRT-local** — the good corner (see the §5.1 tradeoff scatter).
(The bounded-iterative Rymes–Myers scheme, not tabulated here, also reaches that
corner — see the trilemma note above; PCHIP wins on closed form + native format.)

**Validated against MiCASA's own daily product.** The cleanest possible fitter test:
MiCASA *ships* daily NPP/Rh (`daily_1x1`) — the sub-monthly truth the monthly→fit step
discards — so I evaluated both production fits at daily resolution and compared to
that actual daily product, in the fitter's own space, over all 1° land (2020;
**5.97 M cell-days**; `fitter_diagnostics/piqs_vs_pchip_daily_truth.r`). Crucially
this is **not circular**: `daily_1x1` is MiCASA's *native* daily model output regridded
0.1°→1° by the ingest, **not** interpolated from monthly means — and my monthly means
(which the fit ingests) are the monthly *means* of that same daily stream, so it is a
genuine independent target. It is also smooth within-month (~4% day-to-day, not
weather noise), so the comparison is a clean interpolation test, free of scale,
quantity, weather, or partitioning confounds — the model's own output at the
production scale. Result: the two fitters reconstruct the daily NEE with **equal
fidelity** (median RMSE within 0.5%; PCHIP closer at 59% of cells) — PCHIP costs no
accuracy — but **PIQS reaches that fidelity only by reconstructing wrong-sign (source)
GPP at 12.7% of cell-days, vs PCHIP's 0.12%** (~108× fewer). The true daily flux
essentially never shows source GPP; PIQS's overshoot invents it, PCHIP doesn't. So the
PCHIP case, against ground truth, is **sign-physicality at equal accuracy**.

**Selectable alternatives** (not defaults): PPM, minmod/MUSCL, ATP-kriging, MSS,
PIQS all remain selectable via `MICASA_FIT_RDA`; PPM was briefly defaulted then
reverted (continuity — see §5). The on-disk format and all monthly+ budgets are
identical across them.

---

## 2. Diurnalization — framework unchanged; V2 defaults the respiration driver to soil

**The diurnal *framework* is V1's, unchanged** — GPP ∝ ERA5 shortwave, respiration ∝
Q10 of temperature (Olsen & Randerson 2004). **The one V2 change is the respiration
driver variable: V2 evaluates that Q10 on 0–7 cm soil temperature (`stl1`) rather than
2-m air temperature (`t2m`).** This is a deliberate, validated change to the default
product (justified below), not the byte-identical no-op the rest of §2–§4 are; its
effect is small at the NEE level (+2.3% diurnal amplitude) and concentrated in the
boreal cold season where the air-temperature proxy is least physical.

The full hourly product has been regenerated on the soil default: the entire record
(`ERA5/fluxes_YYYYMM.nc`, **303 months, 2001-01 … 2026-03**) was re-diurnalized in place
on 2026-06-21, and every file now carries `respiration_temperature_driver = soiltemp`.
The regeneration is mass-preserving (GPP monthly mean unchanged to the bit; resp/NEE
monthly means conserved to ~1e-6 vs the airtemp reference, resp July diurnal amplitude
ratio **0.831**) and leaves the global annual NEE budget unchanged (§0).

The legacy diurnal product remains fully reproducible by selecting the two legacy
defaults together — `MICASA_RESP_DRIVER=airtemp MICASA_POLAR_CLIP=plain` — which is
**byte-identical to the V1 canonical product**, verified by a committed `ncdiff` run
([`fitter_diagnostics/bytecheck_resp_driver_default.txt`](../fitter_diagnostics/bytecheck_resp_driver_default.txt):
max |Δ| = 0 for GPP/resp/NEE, the airtemp-selected code vs the canonical
`ERA5_2020_pchip/fluxes_202007.nc`), run-and-diffed, not argued from source. So the V2
diurnal changes (the soil driver here, the conserve clip in §3.2) are **reversible to
the bit** by env var — both are new defaults because the evidence supports them, but
nothing about V1's behaviour is lost.
(The **Lloyd-Taylor** response function `MICASA_RESP_TEMPFUN` stays **opt-in,
default-off** — its within-day effect is unvalidated; see that doc §5.3, §5.1–5.2 for
the measured shadow-diffs.)

**Decision (shipped in this release): soil temperature is the default respiration
driver.** The case rests on the *seasonal* eddy-covariance result plus mechanism, with
the within-day relationship a wash — so the flip is at worst neutral and at best correct
(below). The
*measurable* effect of the driver, on the matched full-year-2019 PCHIP air-vs-soil pair
(all 12 months; `fitter_diagnostics/resp_driver_blockboot.py`, committed output
`fitter_diagnostics/resp_driver_blockboot_2019.txt`), is a **damping and phase-lag of
the imposed respiration diurnal cycle**:

- The global area-weighted **respiration** amplitude ratio soil/air is **0.80, 95% CI
  [0.78, 0.83]** (10° spatial block bootstrap, annual), pulled down by the **boreal
  band 0.61 [0.58, 0.63]** — snow-insulated/frozen soil decoupled from swinging winter
  air. (SH-temperate 0.94 [0.86, 1.07] is the one band whose CI spans 1.) The **NEE**
  effect is small: amplitude ratio **1.023, 95% CI [1.021, 1.024]** (every month
  excludes 1; survives a conservative 20°-block CI). *(The 10° block respects spatial
  autocorrelation; a naive i.i.d.-cell resample gives a ~16×-tighter, invalid CI —
  block-bootstrap width 0.0032 vs i.i.d. 0.0002, both from the committed output.)*

**Crucially, this amplitude ratio is not evidence that soil is the *right* driver** —
respiration is a monotone function of its driver, so the ratio merely tracks the
`stl1`/`t2m` amplitude ratio (a self-consistency check on the implementation, not an
observation). The independent test is eddy covariance, and I have now run it
(`fitter_diagnostics/ec_resp_driver_validation.py`, **14** AmeriFlux sites after a
u\*>0.2 turbulence filter and raw — non-gap-filled — flux only; at night NEE ≈
ecosystem respiration). It separates two questions the first pass conflated:

- **Seasonal driver** (which temperature tracks the *seasonal magnitude* of
  respiration): soil wins decisively at **12 of 13 decisive sites** (1 air, 1 tie;
  binomial p = 0.003; median ΔR² +0.042), robust to the u\* filter and a by-night
  block bootstrap.
- **Within-day driver** (which temperature drives the *sub-daily shape* — the only
  thing the diurnalization controls, since it rescales respiration to the monthly
  mean): **neither.** Removing each night's mean, the within-night respiration anomaly
  is explained at **R² ≈ 0.003 by both** air and soil temperature (soil 2 / air 2 /
  tie 10 across sites; p = 1.0). The sub-daily temperature–respiration signal is
  **below the eddy-covariance noise floor.**

So the EC data confirm soil is the better *seasonal* driver, and at the within-day
timescale the driver actually controls — the diurnalization **preserves the monthly
mean** — the two are statistically indistinguishable. The earlier "soil wins 16/20"
was the seasonal cycle in disguise (a metric-vs-use mismatch an adversarial review
flagged); the honest within-day verdict is a tie, not a soil win.

**Why soil is the default.** The within-day tie means switching to soil carries
**no measured penalty** to the diurnal shape — the only thing the driver sets. On every
axis where the two drivers *can* be distinguished, soil is at least as good: it is the
seasonally better descriptor of the temperature–respiration relationship (12/13,
p=0.003), it is the mechanistically correct variable (decomposition responds to soil,
not air, temperature), and its damped, lagged imposed cycle (0.80 amplitude ratio,
+1 h) is the more conservative choice given respiration shows no strong within-night
temperature response of either kind. The switch is free (`stl1` already loaded), mass-
conserving, and reversible to the bit (`MICASA_RESP_DRIVER=airtemp`, plus
`MICASA_POLAR_CLIP=plain` for the full legacy product). So V2 defaults to
soil. I do **not** overclaim a within-day improvement — there isn't one in the tower
data; the argument is "principled default, with no measured downside." Caveats, all in the script output: r(TA,TS)=0.87 /
VIF 4.1, so the competitive-regression betas are variance-inflated and count as **one**
evidence line with the seasonal separate-fit test, not two; **3 of 14** sites show an
unphysical soil Q10 > 3.5 (the signature of seasonal-range compression in the
whole-record fit, not gap-fill — this run is raw-only); and **108** candidate sites
were dropped for lacking a raw soil-temperature sensor, biasing the sample toward
soil-instrumented towers. **Lloyd-Taylor** stays opt-in — unlike the driver variable,
it materially changes respiration amplitude (1.5–3.7×) on an uncertain low-T
sensitivity, and its within-day effect is likewise unvalidated, so it needs its own
gate before any flip.

![ERA5 forcing: 0–7 cm soil temp lags & (per-cell) damps vs 2-m air](figures/resp_forcing_t2m_vs_stl1.png)

![Respiration diurnal cycle: soil driver damps & lags the air driver](figures/resp_diurnal_air_vs_soil.png)

![Per-cell respiration diurnal amplitude ratio soil/air (< 1 = damped)](figures/resp_amplitude_ratio_hist.png)

---

## 3. Other product-number changes (A) — each justified

### 3.1 Aggregation latitude-weight bug fix — V1 was wrong, V2 is correct
V1's 0.1°→1° aggregator (`lib/ingest_common.r:aggregate.to.1x1`) recycled the
cos-latitude area weights **column-major**, applying them along the *longitude*
axis instead of latitude (with a dead ×10/÷10 inner loop). V2 builds a flat
length-100 weight vector that assigns each sub-cell its correct latitude weight.
This is a **genuine bug**: V1 area-weighted the wrong axis. Impact is small for
smooth fields (typically < 0.01%) and grows toward the poles where the cos-lat
gradient across a 1° block is largest. **Verification:** `tests/test_aggregate.r`
(regression test) pins the corrected weighting against the analytic spherical
area; `lib/test_ingest_bitident.r` confirms the read path. Justification is not
"I prefer V2" but "V1 mis-weighted; V2 matches the analytic cos-lat area."

### 3.2 Polar-night GPP = 0 clip — now mass-conserving by default
Physical: no incoming shortwave ⇒ no photosynthesis. The clip zeros GPP wherever
`ssrd == 0`, removing the small residual the sub-monthly quadratic otherwise
leaks into dark hours (a spot check of `fluxes_202512.nc`: ~2.6% of cells touched,
max |GPP| = 9.4e-9 mol m⁻² s⁻¹ — illustrative, not a verify_v2 check; Check 12.2
verifies >75 N GPP = 0). The *plain* zero-clip drops that dark-hour flux outright,
opening a small **GPP monthly-mean gap** — **not** a global 1.5%: Check 2.2 measured
the GPP monthly-mean rel-diff at **p50 ≈ 0.16%, p99 ≈ 2.0%** (the tail is the
partial-polar-night high-latitude cells; the "~1.5%" elsewhere in this doc is a rough
high-latitude figure, refined here).

**V2 closes that gap by default.** The default (`MICASA_POLAR_CLIP=conserve`) invokes
`polar.night.renorm` (`lib/diurnal.r`): after zeroing dark-hour GPP, it
**redistributes the clipped uptake onto each cell's remaining lit hours** (a uniform
per-cell rescale, preserving the ssrd-proportional shape) so the monthly mean is
**restored exactly** — precisely at the partial-polar-night cells where the gap lived.
A cell that is dark *all month* (full polar night) has no lit hours to redistribute
onto, but its monthly-mean GPP should be ≈0 there anyway, so it stays zeroed (the
residual removed is the fit's spurious leak). So on this default the delivered
(post-diurnalize) field is mass-preserving at high latitudes too — by construction —
closing the one place §0's integral-preservation did not previously hold in the
delivered field. Setting
`MICASA_POLAR_CLIP=plain` restores the legacy zero-clip (byte-identical to the
pre-conserve product). Unit-tested (`tests/test_diurnal.r`: mean restored to ~1e-12,
dark hours stay 0, full-dark cells don't blow up). **Verification:** Checks 2.2, 12.2, 17.1.
*(Staged for v2.2.0: the code default is `conserve`; the production record is
re-diurnalized onto it on the next pass — see CHANGELOG.)*

### 3.3 ERA5 dual-tree FastTrack fallback
Only affects NRT trailing months the primary ERA5 tree has not yet populated;
those days fall to the lower-latency `ea_0005` FastTrack stream (the *same* ERA5
product, earlier release). For any month the primary covers, the path is
unchanged. Per-day provenance is written (`meteo_source_by_day`, e.g.
`primary:1-30 fasttrack:31`). **Verification:** Checks 1.4, 10.1; first
production use 2026-Q1 (2026-02/03 wholly FastTrack), clean files.

### 3.4 Per-month climatology auto-detect
V1 chose real-vs-climatology per *year* from a hand-set `MICASA_CLIM_YEARS` list,
with no file-existence check — a partially-published year forced either
climatologising real months or crashing on unpublished ones. V2 decides per
*month* by file presence: real monthly file present ⇒ use it, else day-of-year
climatology. For fully-published months the path is identical. **Verification:**
Check 1.4; 2026-Q1 run (Jan–Mar real via PCHIP, later months climatology, no
crash).

---

## 4. Behavior-preserving changes (B) — proven no-ops on the product

Each item below changes *no flux value*; the proof is in the right column.

| Change | Proof it preserves the product |
|---|---|
| `diurnal.flux` / `polar.night.clip` extracted to `lib/diurnal.r` | Byte-for-byte identical on random arrays; `tests/test_diurnal.r` (21) |
| Fitter cores extracted to `lib/{pchip,mss,ppm,linmm}_fit.r` | Function bodies unchanged; unit tests `test_{pchip,mss,ppm,linmm}_fit.r` (12/10/13/11) |
| ERA5 path helpers → `lib/era5_meteo.r` | `tests/test_era5_meteo.r` (11) on resolver + run-length encoder |
| Grid-area fns `archimedes`/`compute.gca` made pure | `tests/test_ingest_geometry.r` (20) vs analytic 4πR² |
| `compute_clim` PyFerret → xarray | Algorithm exact to 1e-12 vs hand-computed mean; `tests/test_compute_clim.py` |
| `check_bounds` NCO `ncwa` → xarray | Pure `flux_to_tgc_per_year`; `tests/test_check_bounds.py` (7). NCO version never actually ran (guarded by `|| true`). |
| Ingest skip-existing + read-only-needed | `ncdiff` 4 days × 4 tracers max \|Δ\| = 0; `lib/test_ingest_bitident.r`; 610→504→4 s |
| Compression deflate 9 → 4 (**diurnalize output only**; ingest stays at 9, `lib/ingest_common.r:149`) | Lossless codec ⇒ data bit-identical; only size +0.3% / time −39% (`lib/bench_compression_diurnal.r`). Codec argument, not ncdiff-run. |
| Provenance CF/ACDD attributes | Additive global attributes only; `tests/test_provenance.{r,py}` (26 ea); Checks 23.1–23.3 |
| Per-step run manifest | Additive `jobs/run_manifest.tsv`; never aborts caller; `tests/test_manifest.r` (15); Checks 22.1, 24.1–24.2 |
| Sub-monthly sign-flip logging | Log lines only; drives Check 3.1; no flux touched |
| Download verify scoped to year | Verifies *which* files, not their content; same files for a given year |
| Op bug fixes: `sbatch_wait` comma, `check_hashes` glob, `compute_daily_clim` nullglob, hardcoded paths | Fix crashes/skips, not numbers; `tests/test_check_hashes.py` (12); 2026-Q1 multi-scenario run |
| verify_v2 harness edits (6.2→INFO, 11.1 log-age, 5.1/5.2 partial-year, 1.4 dual-tree) | Change what is *checked*, not what is *produced*; each justified in CHANGELOG 2026-05-16 |
| Public-release packaging (LICENSE CC0, README split, CITATION.cff, CI) | No pipeline effect |

The CI (`.github/workflows/ci.yml`) byte-compiles Python, `bash -n`s every shell
script, `parse()`s every R script, and runs the behavior tests on every push — so
these refactors cannot silently regress.

---

## 5. Considered and rejected (diligence, not changes)

Documenting what was *not* changed, and why, is part of the justification:

- **ATMC budget closure (NEE = Rh − NPP, *not* − ATMC)** — tried 2026-04-29,
  reverted same day. This is the most consequential "rejected" choice — it changes
  the sign of the prior's long-term trend — so it gets its own treatment in **§5.2**.
- **PPM as default** — briefly defaulted 2026-06-18, reverted: daily fidelity is a
  near-tie that the two bootstraps *split* — on the **pooled cell-month** metric PPM
  is marginally ahead (paired Δ ≈ 0.7% of the median level, PPM better in 54% of
  cell-months), while the more appropriate **by-cell product bootstrap**
  (`FITTER_COMPARISON.md` §4.6) puts PCHIP significantly but negligibly ahead (~0.3%
  of RMSE, CI excludes 0). The margin is immaterial either way; PPM is reverted
  because it reintroduces month-edge **discontinuities at ~70% of edges** — the steps
  the smoother exists to remove (CHANGELOG 2026-06-18).
- **MSS** (overshoots despite the name, ~24% wrong-sign GPP knots, ~hours/grid),
  **linear-recursion PIQS** (unstable), **constrained-quadratic PIQS** (dominated
  by PCHIP), **CCGCRV** (not pursued) — FITTER_COMPARISON.md §4.1/§5, PROPOSALS
  #9/#11/#6.

---

### 5.1 Why not "PIQS, then revert to linear on overshoot"

This hybrid — keep PIQS's smooth global-solve quadratic where it is sign-safe and
patch **only** the overshooting pieces with a sign-safe integral-preserving linear
— was the stakeholder-preferred alternative to PCHIP, and it has a **genuine
motivation I state plainly**: PIQS's global solve makes its flux **near-C¹**
(continuous *slope*), whereas PCHIP-on-cumulative is only **C⁰ in the flux** — the
flux carries a small slope kink at each month knot. Measured
(`FITTER_COMPARISON.md` §4.5): the knot derivative-jump (×width/env) is **PIQS
0.000 vs PCHIP 0.290**. A *selective* fallback would therefore keep PIQS's superior
derivative-smoothness on the ~71% of sign-safe pieces and patch only the rest. That
is the strongest case for it, and it is real.

I implemented and measured it (`fitter_diagnostics/piqs_hybrid.r`,
`linear_fallback_quantify.r`; committed output
`linear_fallback_quantify_20260621.txt`) and **did not adopt** it, because that C¹
advantage is inconsequential for *this* product and is outweighed on locality and
continuity:

0. **The C¹ edge does not reach the delivered prior.** The shipped hourly NEE is
   the fit's monthly mean redistributed by ERA5 hourly meteo; the smoother sets
   only the small within-month *deviation* term, and the hourly field is sampled,
   not differentiated. A C⁰-vs-C¹ distinction in the flux *slope* at month knots is
   below the ERA5 redistribution it rides on — which is exactly why the daily
   fidelity is a tie (point 3). The smoothness PIQS preserves is aesthetic here,
   not a measurable property of the prior.
1. **Non-locality is not fixed.** The linear fallback is applied *post-hoc* to
   PIQS's already-solved coefficients; it does not decouple the knots. A revised
   NRT month still re-solves PIQS and **rewrites the entire ~303-month record**
   (the **NRT footprint** row of the §1 scorecard) — the exact disqualifier PCHIP
   avoids (footprint 0). This applies to the *selective* fallback, i.e. Andy's
   actual proposal, not a strawman.
2. **It injects a genuine discontinuity where PCHIP injects none — reported
   denominator-free.** **29.3%** of land cell-months trigger the fallback, and
   patching breaks PIQS's C⁰ continuity there. The honest, denominator-free
   statement: the hybrid's **total absolute discontinuity budget is 0.97 mol m⁻²
   s⁻¹** summed over 1.44 M patched land edges, versus **exactly 0** for PCHIP (C⁰
   by construction). An envelope-normalized framing (“52% exceed 3× env”) is
   misleading here: the overshooting pieces are near-zero-transition months, so
   **38% of patched edges have envelope ≈ 0**, where “jump/env” blows up by dividing
   by ~0, not by being physically large (among the 62% of edges with a well-defined
   envelope the median jump is only **0.25× env**). The absolute budget is the
   defensible cross-method number.
3. **No fidelity gain.** On the robust **median**, hybrid daily RMSE/env (2020) is
   **0.079** vs PCHIP's **0.081** — a tie; the means (0.139 vs 0.151) are
   tail-sensitive and not a meaningful difference. The preserved smoothness buys no
   reconstruction accuracy, consistent with point 0.

(*A distinct, weaker proposal — “use continuous linear **everywhere**”, PROPOSALS
#9, which is **not** the selective fallback above — is worse than PIQS on its own
terms: the recursion `yᵢ₊₁ = 2·mᵢ − yᵢ` flips sign at 36.9% of interior knots and
rings with unbounded resonance (knot/env p99 ≈ 2.6×10⁵, a Nyquist pole). I note it
only to close the option; it is not Andy's proposal.*)

So against the *selective* fallback, PCHIP wins on locality and continuity while
conceding a real but immaterial C¹-flux advantage; it is sign-safe *and*
C⁰-continuous *and* local *and* closed-form — the "good corner" below:

![Sign-safety vs continuity: PCHIP in the good corner](figures/fitter_tradeoff_scatter.png)

![PIQS+linear-fallback patch-discontinuity distribution](figures/linear_fallback_discontinuity.png)

### 5.2 ATMC, and the sign of the prior's long-term trend

This is the choice most likely to be contested.

**What ATMC is.** NCCS publishes an "atmospheric correction" `ATMC` field with the
file comment `NEE = Rh − NPP − ATMC`. Per Weir et al. (2021) it is the Low-order
Flux Inversion (LoFI) empirical sink, `S_m = α_yr·max(T_m−T_{m-1},0)/10·HR_m`, with
α scaled **each year so the global biospheric total matches the observed
atmospheric CO₂ growth rate** (~3 PgC/yr, concentrated NH-extratropics JJA).

**The stakes are not cosmetic — and they hit both the level and the trend.**
Subtracting ATMC more than *doubles the mean biospheric sink*, from **−2.45 to
−5.99 PgC/yr** (a 3.5 PgC/yr shift; the ATMC field itself is ~3 PgC/yr), and *flips the sign of the long-term
trend*: CASA-only NEE trends **+0.0413 PgC/yr/yr**, with ATMC **−0.0067** (i.e.
essentially flat). The trend alone compounds to ≈ **+1.1 PgC/yr** of drift over the
25-yr record — itself of order the mean sink. So ATMC is a first-order change to
both the magnitude and the time-evolution of the prior, not a rounding term.

**Why I still do not subtract it — and this holds whether the trend is real or a
CASA bias.** These fluxes are a **prior to a CO₂ inversion that itself assimilates
atmospheric CO₂**, and ATMC was tuned to *that same observation class* (the global
growth rate). Subtracting it pre-loads the prior with the very constraint the
inversion exists to apply — **data leakage / double-dipping** — after which the
inversion can no longer *independently* constrain the long-term sink, because the
answer is already baked in. The growth-rate constraint belongs in the inversion's
assimilation step, applied once, against a prior that reports what the offline
biosphere model says **on its own**. Crucially this argument does **not** require
me to claim the +0.0447 trend is physically correct:
- *If* it is real (e.g. CO₂-fertilization / greening strengthening NPP, which CASA
  represents through satellite-APAR forcing — Zhu et al. 2016), the inversion keeps
  it and the data confirm it.
- *If* it is a CASA structural bias (e.g. warming-driven respiration outrunning
  NPP), the inversion corrects it from the atmospheric data.

Either way the prior's job is to carry CASA's *own* estimate; baking in ATMC
forecloses the correction in both cases. So I ship the +0.0447 trend as a property
of the CASA prior — **not** asserting it is a "real climate feature," only that
pre-correcting it would be circular.

**When ATMC *would* be appropriate.** If these fluxes are ever used *outside* an
inversion — forward site-level comparison against obs, or as a fixed ensemble
member with no further optimization — there is no double-dipping and the ATMC
subtraction is the right choice. For the current prior-to-inversion use it is not.
(Trend figures: the **+0.0413** in the 2026-04-29 ATMC table is the *PIQS-era*
CASA-only trend; **+0.0447** is the later *PCHIP* Section-15 value (2026-05-04 run);
**+0.04** elsewhere is the same number rounded — one ~+0.04 PgC/yr/yr trend, three
runs. The +0.0413→+0.0447 difference is between two **different-date pipeline runs**
with other V2 changes intervening — *not* a fitter effect; the budget-invariance of
§0 is a statement about a **pure fitter swap on a fixed pipeline**, which is why I
do not claim the cross-date runs match to floating point.)

## 6. Evidence matrix

**Validation harness — `verify_v2` (60 distinct checks / 24 sections).** Phase 1
structural (1.1–1.4); Phase 2 transformation + sanity (2.1–2.4 mass/integral,
5.1–5.3 global/YoY/seasonal); Phase 3 cross-boundary + spatial-vs-v1 + provenance
(4.1, 6.1–6.2, 7.1–7.4, 8.1–8.3, 9.1–9.2, 10.1, 11.1–11.2); Phase 4 edge cases +
biome cells + trends (12.1–12.2, 13.1–13.2, 14.1–14.3, 15.1–15.3 trend/ENSO/COVID,
16.x diagnostics); Sections 17 diurnal integrity, 18 PCHIP invariants, 19
additional biomes, 20 cross-product, 21 robustness, 22 performance, 23 provenance,
24 manifest (§24 = observability meta-checks on the run log, *not* product
assertions). Committed run `verify_v2_summary_20260621.txt`: **54 PASS / 8 INFO / 0 FAIL / 0
WARN** — every product / science / provenance check passes; the rest are INFO
context. The checks that earlier needed attention are now clean:
- **§3.1 / §20.1** — v2-vs-v1 per-band annual NEE agrees to **0.04%** over 2001–2024
  (`fitter_diagnostics/check_20_crossproduct.py`): the §3.1 aggregation fix shifts no
  band-level mass (boreal 0.04%, the largest, matches its poleward-growing
  prediction; diurnalize preserves monthly means, so this is the shipped per-band
  annual NEE to the polar-clip residual).
- **§20.2** — MiCASA global NBE **+0.99 PgC/yr** (2001–2024), ~3.6 PgC/yr off the GCB
  land sink ≈ the ATMC term: CASA-only does not self-close the growth-rate budget, by
  design (§5.2). Budget *context*, not a closure.
- **§11.1** (job-log error scan) clean. The ATP-kriging diagnostic crash-logs that
  once flagged it (a singular kriging system on dormant cells) are superseded — the
  production `lib/atpk_fit.r` guards that case (`solve` → ridge → flat **dormant**
  fallback, never halting; unit-tested in `test_atpk_fit.r`), so it cannot recur.
- **§24** (run-manifest meta-checks on the working-directory log, *not* product
  assertions) clean.

**Unit tests — all green on Orion (R 4.4.0 / Python, 2026-06-21); the 143
non-`quadprog` R checks reproduced green locally (R 4.6.0) for this revision:**

| Suite | Checks | Guards |
|---|---|---|
| test_pchip_fit.r | 12 | PCHIP mass / C⁰ flux continuity (= C¹ of F; not C¹ flux) / sign-flip-rate (not sign-definiteness — see §1) |
| test_diurnal.r | 21 | diurnalize transform + q10/lt factors |
| test_atpk_fit.r | 14 | ATP coherence/variance/sign |
| test_ppm_fit.r / test_linmm_fit.r | 13 / 11 | PPM & minmod mass/limiter |
| test_mss_fit.r | 10 | MSS QP fit — **requires `quadprog`; runs on Orion, SKIPs without it** |
| test_ingest_geometry.r | 20 | spherical area weights |
| test_era5_meteo.r | 11 | FastTrack resolver |
| test_manifest.r | 15 | manifest format / no-abort |
| test_provenance.r / .py | 26 / 26 | CF/ACDD attributes |
| test_check_hashes.py / test_check_bounds.py / test_compute_clim.py | 12 / 7 / 10 | hashing / unit conv / clim mean |

R total: 143 host-portable + 10 `quadprog`-gated (MSS) = 153.

**Per-change → guard map** (numbers-changing items): fitter → 2.1, 3.1, 18.1,
18.2 + test_pchip_fit; polar-night → 12.2, 17.1; aggregation fix → test_aggregate
+ test_ingest_bitident; FastTrack → 1.4, 10.1; per-month clim → 1.4. Every
behavior-preserving item maps to a proof in the §4 table.

---

## 7. Known residual limitations

- **+0.04 PgC/yr/yr** long-term trend in CASA-only NEE is shipped as a property of
  the CASA prior — *not* asserted to be a real climate feature. Whether real
  (CO₂-fertilization/greening) or a CASA structural bias, the inversion corrects it
  from independent atmospheric data, and pre-closing it with ATMC would be circular
  (§5.2). Sign-of-the-trend stakes: with ATMC the trend is −0.0067 (flat). The
  trend is also **non-stationary** (Check 16.2: first-half 2001–2012 +0.0274,
  second-half 2013–2025 +0.1031 — ~3.8× steeper), which only strengthens the case
  for letting the inversion, not a baked-in correction, resolve it.
- **Polar-night clip** — the legacy plain zero-clip left a small high-latitude GPP
  monthly-mean gap (Check 2.2: ~0.16% median, ~2% p99 cell-month). V2 now defaults to
  the **mass-conserving clip** (`MICASA_POLAR_CLIP=conserve`), which redistributes the
  clipped dark-hour uptake to restore the monthly mean exactly, so the *shipped* product
  is mass-preserving at high latitudes too (§0/§3.2). `MICASA_POLAR_CLIP=plain` reverts
  to the legacy zero-clip.
- **PCHIP is not sign-definite everywhere** — it cuts sub-monthly sign flips
  16–60× vs PIQS but leaves a small bounded residual (≤0.94% of GPP cell-hours;
  reproduced in `fitter_diagnostics/pchip_sign_definiteness.r`), mopped up by the
  clip. "Eliminated by construction" would be an overstatement (§1).
- **Respiration driver — soil is now the V2 default.** Its eddy-covariance gate was run
  properly: soil is the better *seasonal* respiration driver (12/13 AmeriFlux sites,
  p=0.003), and the *within-day* relationship the diurnalization actually sets is a tie
  (R²≈0.003 both, below the EC noise floor; §2). With no measured within-day downside
  and soil seasonally + mechanistically correct, V2 defaults to it — and the validation
  is explicitly that it is a **principled default with no measured downside**, not a
  demonstrated within-day improvement (the honest limit of what the tower data show).
  The legacy air path stays selectable and byte-identical (`MICASA_RESP_DRIVER=airtemp`).
  **Lloyd-Taylor** stays opt-in pending its own gate (it materially moves respiration
  amplitude on an uncertain low-T sensitivity).
- **Prior uncertainty is constructed, not native.** MiCASA ships **no per-pixel
  uncertainty** (a single deterministic realization — vars `NPP/Rh/FIRE/FUEL/ATMC/NEE`
  only), so any prior σ is one I build. I can bound **two distinct, small
  components — which should not be summed into a single "~3%"**: (i) 0.1° sub-grid
  heterogeneity within a 1° cell, which is strongly **biome-dependent, ~1% (boreal)
  to ~10% (temperate mosaic)**, median ~3.5%; and (ii) across-fitter structural
  spread, **~3%**. Both are emitted via the opt-in `NEE_sd` field (from the
  ATP-kriging variance; `FITTER_COMPARISON.md §4.3`). Together they are a **lower
  bound on the sub-monthly-redistribution + 1°-representativeness error *only*** —
  they explicitly **exclude the dominant term**, the model error in MiCASA's
  monthly NPP/Rh itself (tens of %), which the product does not carry and which the
  inversion's prior error covariance must supply. So this is a floor on two minor
  components, **not** an uncertainty on the prior as a whole — do not read "~3%" as
  "the prior is good to 3%."
- **Archival DOI** ships as `PENDING` (`grep -rl PENDING` finds every spot).

---

## Appendix A — Fitter equations (PIQS vs PCHIP)

Both fitters store, per cell and month *i*, a quadratic on `t ∈ [tᵢ, tᵢ₊₁]` of
width `hᵢ`, and both impose **mass preservation** — the piece integral equals the
MiCASA monthly mean `ȳᵢ` (this *is* the master invariant of §0):

```
fᵢ(t) = aᵢ (t−tᵢ)² + bᵢ (t−tᵢ) + cᵢ
(1/hᵢ) ∫[tᵢ→tᵢ₊₁] fᵢ dt = aᵢhᵢ²/3 + bᵢhᵢ/2 + cᵢ = ȳᵢ
```

**PIQS** (Rasmussen 1991) fixes the remaining freedom by a **single global solve**:
each piece preserves its integral *and* adjacent pieces share the knot value (C⁰),
`fᵢ(tᵢ₊₁) = fᵢ₊₁(tᵢ₊₁)`. That continuity system couples *every* month to every
other → non-local; nothing constrains the quadratic's sign, so it overshoots
through zero in sharply seasonal cells.

**PCHIP-on-cumulative** (Fritsch & Carlson 1980; `lib/pchip_fit.r`) instead works
on the cumulative integral and is **local**:

```
Fₖ = Σ_{i<k} ȳᵢ hᵢ           (F₀ = 0; monotone when the ȳᵢ share a sign)
secants     mₖ = ȳₖ
F-C knot slopes dₖ:  dₖ = 0 at a secant sign change,
                     else |dₖ| ≤ 3·min(|mₖ₋₁|, |mₖ|)   ← monotonicity limiter
```

The flux is the derivative of the monotone cubic Hermite on `F`; on segment *k*
with `s = (t−xₖ)/hₖ`,

```
f(s) = (6s−6s²)·mₖ + (3s²−4s+1)·dₖ + (3s²−2s)·dₖ₊₁
```

which in the stored `(a,b,c)` form is, with `Q = −6mₖ+3dₖ+3dₖ₊₁`,
`L = 6mₖ−4dₖ−2dₖ₊₁`, `K = dₖ`:

```
aᵢ = Q/hₖ²,   bᵢ = L/hₖ,   cᵢ = K          (signs negated for GPP ≤ 0)
```

Mass is automatic (`∫₀¹ f ds = mₖ = ȳₖ`). PIQS sets `(a,b,c)` by a *global* C⁰
system; PCHIP from *local* Fritsch-Carlson knot slopes `dₖ`. Both yield the
identical `ȳᵢ` — hence the budget-invariance (§0).

## 8. References

**Sub-monthly fitter**
- Rasmussen (1991), *Piecewise integral splines of low degree*, Computers & Geosciences 17(9):1255–1263, doi:10.1016/0098-3004(91)90027-B.
- Fritsch & Carlson (1980), *Monotone Piecewise Cubic Interpolation*, SIAM J. Numer. Anal. 17(2):238–246, doi:10.1137/0717021.
- Colella & Woodward (1984), *The Piecewise Parabolic Method (PPM)*, JCP 54(1):174–201, doi:10.1016/0021-9991(84)90143-8.
- van Leer (1979), *Towards the ultimate conservative difference scheme V*, JCP 32(1):101–136, doi:10.1016/0021-9991(79)90145-1.
- Boneva, Kendall & Stefanov (1971), *Spline transformations*, JRSS-B 33(1):1–70.
- Wang & Bartlein (2022), *A Fast Mean-Preserving Spline*, JTECH 39(4):503–512, doi:10.1175/JTECH-D-21-0154.1.
- Rymes & Myers (2001), *Mean-preserving algorithm for smoothly interpolating averaged data*, Solar Energy 71(4):225–231, doi:10.1016/S0038-092X(01)00052-4.
- Kyriakidis (2004), *A geostatistical framework for area-to-point interpolation*, Geographical Analysis 36(3):259–289, doi:10.1111/j.1538-4632.2004.tb01135.x.
- Bartlein, *mp-interp* (mean-preserving interpolation reference code): https://github.com/pjbartlein/mp-interp.
- JULES, *Temporal interpolation* docs: https://jules-lsm.github.io/latest/input/temporal-interpolation.html.

**Diurnalization & ecosystem respiration**
- Olsen & Randerson (2004), *Differences between surface and column atmospheric CO₂…*, JGR 109:D02301, doi:10.1029/2003JD003968.
- Potter et al. (1993), *Terrestrial ecosystem production (CASA)*, Global Biogeochem. Cycles 7(4):811–841, doi:10.1029/93GB02725.
- Lloyd & Taylor (1994), *On the temperature dependence of soil respiration*, Functional Ecology 8(3):315–323, doi:10.2307/2389824.
- Davidson, Janssens & Luo (2006), *…moving beyond Q10*, GCB 12:154–164, doi:10.1111/j.1365-2486.2005.01065.x.
- Reichstein et al. (2005), *On the separation of NEE into assimilation and respiration*, GCB 11:1424–1439, doi:10.1111/j.1365-2486.2005.001002.x.
- Lasslop et al. (2010), *…light response curve approach*, GCB 16:187–208, doi:10.1111/j.1365-2486.2009.02041.x.
- Haynes et al. (2019), *SiB4*, JAMES 11:4423–4439, doi:10.1029/2018MS001540.
- Hersbach et al. (2020), *The ERA5 global reanalysis*, QJRMS 146:1999–2049, doi:10.1002/qj.3803.

**Inversion context**
- Denning, Fung & Randall (1995), *Latitudinal gradient of atmospheric CO₂ … (the rectifier)*, Nature 376:240–243, doi:10.1038/376240a0.
- Weir et al. (2021), *Bias-correcting carbon fluxes* (LoFI / ATMC), ACP 21:9609–9628, doi:10.5194/acp-21-9609-2021.
- Zhu et al. (2016), *Greening of the Earth and its drivers*, Nature Climate Change 6:791–795, doi:10.1038/nclimate3004.
- Friedlingstein et al. (2023), *Global Carbon Budget 2023*, Earth Syst. Sci. Data 15:5301–5369, doi:10.5194/essd-15-5301-2023.
