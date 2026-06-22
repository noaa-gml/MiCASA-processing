# Diurnalization: method comparison and alternatives

**Status:** investigation document · **Date:** 2026-06-20 · **Default scheme:
Olsen & Randerson (2004) unchanged.** This doc records what the current
diurnalization is, the physical assumptions it embeds, and the ranked
alternatives. · **Scope:** the *diurnal* (sub-daily) redistribution in
`diurnalize-ERA5.r` — distinct from the monthly→sub-monthly *fitter*
(`fit.piqs.rda`), which is covered in [FITTER_COMPARISON.md](FITTER_COMPARISON.md).

This is the companion question to the fitter investigation. The fitter sets the
*sub-monthly* shape `qmod(t)`; the **diurnal** shape is imposed entirely by the
ERA5 meteo driver. For an atmospheric-inversion prior the diurnal shape is not
cosmetic — it is the **diurnal rectifier** (Denning et al. 1995): the covariance
between the surface flux and the depth of the boundary layer that ventilates it.
A diurnal cycle with the wrong amplitude or phase biases the rectifier, and the
inversion cannot fully undo a biased prior covariance structure. So the metric
for "is there a better way" is: *does the alternative get the diurnal amplitude
and phase of GPP and respiration closer to what eddy-covariance shows, while
still conserving each month's MiCASA total?*

## 1. What the current scheme is

Both components are **driver-proportional redistribution** (Olsen & Randerson
2004), mass-preserving against the monthly mean. Per cell, per hour `t` within a
month (`lib/diurnal.r` :: `diurnal.flux`):

```
GPP(t)  = GPP_mon  · ( SSRD(t)/SSRD_mon − 1 ) + qmod_gpp(t)
RESP(t) = RTOT_mon · ( Q10(t)/Q10_mon  − 1 ) + qmod_resp(t)
Q10(t)  = 1.5 ^ ( (t2m(t) − 273.15) / 10 )
```

- **GPP** (here negative-for-uptake, `GPP = −2·NPP`) is redistributed
  proportionally to ERA5 surface downward shortwave `ssrd`.
- **RESP** (`= Rh + Rauto`, positive-typical) is redistributed by a `Q10 = 1.5`
  factor on **2-m air temperature** `t2m`.
- `qmod` is the fitter's within-month deviation from the monthly mean; averaged
  over the month each component returns its MiCASA monthly total (exact mass).
- A **polar-night clip** zeros GPP wherever `ssrd == 0`.

The diurnal *shape* comes 100% from the driver ratio. `qmod` only modulates the
slowly-varying monthly envelope.

### 1.1 A structural fact: two of four ERA5 inputs are unused by the flux

`diurnalize-ERA5.r:236` reads **`t2m, ssrd, stl1, swvl1`** and writes all four to
every output file, but **`stl1` (0–7 cm soil temperature) and `swvl1` (0–7 cm
volumetric soil moisture) never enter the flux** — they are carried as
diagnostics only. Respiration, which in MiCASA is overwhelmingly *heterotrophic
soil decomposition*, is driven by **air** temperature. The physically-correct
driver for soil decomposition is already on disk and already loaded.

## 2. The assumptions, and which are questionable

A key algebraic point frames everything below: the redistribution is a **ratio**
normalized to the driver's own monthly mean. Therefore **any *linear* rescaling
of a driver cancels exactly** and has *zero* effect — e.g. converting `ssrd` to
PAR (`PAR ≈ 0.45 · 4.6 · ssrd`) leaves `SSRD(t)/SSRD_mon` unchanged. Only
**nonlinear** driver physics, or **changing which variable** drives a component,
changes the diurnal shape. This rules out a whole class of cosmetic "fixes" and
points at the ones that actually matter.

| # | Assumption | Status | Why |
|---|---|---|---|
| A1 | **Rh driven by air temp (`t2m`)**, not soil temp | **Questionable — highest leverage** | Air temp swings ~10–15 K diurnally and peaks ~14 h; 0–7 cm soil temp is damped (~½ amplitude) and lagged. Driving Rh off `t2m` *overstates* the respiration diurnal amplitude and runs its phase ~2–4 h early. Because Q10 is nonlinear, this changes shape, not just scale. |
| A2 | **Fixed `Q10 = 1.5`, single global value** | Improvable | Literature soil Q10 ~2.0 (Davidson et al. 2006); a Lloyd & Taylor (1994) Arrhenius-type curve has the correct increasing sensitivity at low T. The low 1.5 partly (and accidentally) compensates the over-amplitude of using air temp. |
| A3 | **No soil-moisture limitation on Rh** | Minor *diurnally* | `swvl1` is ~flat within a day → small diurnal effect; matters more sub-monthly. |
| A4 | **GPP linear in shortwave** (no light saturation) | Improvable, not a bug | Real canopy LUE drops at high light → linear-in-light makes the GPP peak too sharp at noon, too low at the shoulders. *But* CASA computes monthly NPP itself as `ε·APAR` (linear LUE; Potter et al. 1993), so linear redistribution is internally consistent with the parent model. A saturating (rectangular-hyperbola) response is more realistic at the hourly scale but needs a saturation parameter. |
| A5 | **GPP and autotrophic respiration share the light/temperature split** | Minor | `Rauto = NPP` is lumped into `RESP` and driven by Q10(Tair); leaf maintenance respiration partly tracks recent assimilation. Second-order. |
| A6 | **No PFT/biome dependence** of Q10 or light response | Realistic but heavy | Needs a PFT map; many parameters. Defer. |

## 3. Ranked alternatives

### Tier 1 — within-framework, defensible, low risk

**(1) Respiration on soil temperature (`stl1`) instead of air temperature.**
Replace `Q10(t) = 1.5^((t2m−273.15)/10)` with the same Q10 evaluated on `stl1`.
One-line change, **zero new inputs** (already loaded). Effect: damps and
phase-lags the respiration diurnal cycle toward what eddy-covariance and soil
chambers show; in particular it keeps overnight efflux flatter instead of
collapsing with the fast-cooling night air. Because nighttime NEE = Rh alone
(GPP clipped) and the nighttime PBL is shallow, this is precisely the term the
rectifier is most sensitive to. **This is prototype #1.**

**(2) Lloyd & Taylor (1994) instead of a fixed Q10.** Reco temperature response
`R = R_ref · exp[ E0 · ( 1/(T_ref−T0) − 1/(T−T0) ) ]` with `T0 = −46.02 °C`,
`E0 ≈ 309 K`. Standard in FLUXNET nighttime partitioning (Reichstein et al.
2005); correct curvature at low temperature where a fixed Q10 misbehaves.
Combines naturally with (1) on soil temperature.

### Tier 2 — within-framework, more effort

**(3) Soil-moisture scalar on Rh** using `swvl1` (e.g. a parabolic or
Davidson-type moisture function). Small diurnal effect; mostly a sub-monthly
refinement. Cheap to add once (1) is in.

**(4) Light-saturating GPP.** Rectangular hyperbola in PAR
`GPP = α·PAR·β / (α·PAR + β)` (Lasslop et al. 2010 daytime model) replacing the
linear `∝ ssrd`. Broadens the midday peak. Needs a saturation parameter `β`
(ideally per-PFT); without a PFT map a single global `β` is a guess. Defensible
improvement, not a correctness fix (see A4).

### Tier 3 — heavy / deferred

**(5) PFT-dependent Q10 and light parameters.** Most realistic, most parameters,
needs a land-cover map co-registered to the 1×1° grid. Defer until Tier-1
results justify it.

## 4. Out-of-framework options considered and rejected

These replace MiCASA's diurnal cycle with *another product's* rather than
redistributing MiCASA's own monthly totals:

- **Native sub-daily LSM output** (SiB4 — Haynes et al. 2019; ORCHIDEE; JULES —
  Best et al. 2011; or CASA-GFED's own hourly stream).
- **Machine-learning hourly upscaling** (FLUXCOM-X / X-BASE — Nelson et al.
  2024; FLUXCOM — Jung et al. 2020), trained directly on eddy-covariance.
- **FLUXNET PFT×month diurnal climatology** applied as a fixed shape (Falge et
  al. 2001-style mean diurnal cycles).

All three break the **"MiCASA-conserving prior"** identity: the product would no
longer be MiCASA redistributed in time but a blend of MiCASA monthly totals with
a different model's or dataset's sub-daily physics. Worse, the ML and FLUXNET
options pull in an observational basis (eddy covariance) that overlaps what the
downstream CO₂ inversion assimilates — the same **double-dipping** objection that
made us reject the ATMC correction (METHODOLOGY.md, "Why NEE = Rh − NPP"). Not
recommended.

## 5. Recommendation

Keep Olsen & Randerson (2004) as the framework. Pursue Tier-1 improvements in
order, each behind a default-off env flag so the canonical product is unchanged
until a shadow-diff quantifies the effect:

1. **`MICASA_RESP_DRIVER=soiltemp`** — Rh on `stl1` (prototype #1). Quantify the
   diurnal Rh/NEE amplitude and phase shift vs the canonical air-temp run on a
   representative month.
2. If (1) looks right, add **Lloyd-Taylor** and the **soil-moisture scalar** as
   further opt-in refinements.
3. Light-saturating GPP only if a PFT map / saturation parameter can be
   justified.

The shadow-diff uses the existing `MICASA_DIURN_OUT_DIR` / `MICASA_DIURN_ONLY_MONTH`
test harness — no change to the production path to measure the candidate.

## 5.1 Prototype #1 result — July 2020 shadow-diff

> **Note:** §5.1–5.2 are the original *single-month* (July/Jan 2020) prototype
> diagnostics. The headline numbers used for the recommendation are the **full-year
> 2019 spatial block-bootstrap** values in §5.4 (e.g. boreal resp ratio 0.61, not
> the 0.83 quoted here for July) — these single-month figures are kept for the
> shadow-diff record but are superseded by §5.4.

Implemented as `MICASA_RESP_DRIVER={airtemp|soiltemp}` (default `airtemp`,
byte-identical to legacy; `q10.factor()` extracted to `lib/diurnal.r` and
unit-tested). July 2020 was diurnalized both ways into shadow dirs with the
identical fit, and the global-land area-weighted diurnal cycles compared
(`fitter_diagnostics/resp_driver_shadowdiff.r`). Findings:

- **Pure redistribution confirmed.** Per-cell monthly means are unchanged to
  floating point: `resp` max |air−soil| = 3.8e-14 (rel 2.8e-9), `NEE`
  7.6e-14 (rel 1.4e-8), `GPP` exactly 0. The change touches only the
  sub-daily *shape*, not any monthly total.
- **Respiration diurnal cycle damps and lags, as predicted.** Per-cell
  amplitude ratio soil/air = **0.86** (median, area-weighted; boreal 0.83,
  tropics 0.85, SH-temperate 0.79); global-land mean diurnal range ratio
  **0.81** (a ~19% reduction). Phase shift **+1 h** (peak 15→16 UTC): soil
  respiration stays elevated into the evening instead of collapsing with the
  fast-cooling night air — the physically-expected behaviour.
- **Net NEE effect is small but in the right direction.** Because GPP
  dominates the NEE diurnal swing (~12× the respiration swing), the
  respiration change moves NEE diurnal amplitude by only **+2%** (median
  ratio 1.02), with no phase shift. So the switch is a *refinement* of the
  rectifier-relevant overnight respiration, not a disruption of the NEE the
  inversion sees.

This is the basis for a defensible recommendation to Andy: the soil-temp
driver is more physical, costs nothing (data already loaded), conserves every
monthly mean exactly, and changes the product modestly (~2% NEE diurnal
amplitude) in the correct direction. Default remains `airtemp` pending
sign-off on the §5.4 recommendation below.

## 5.2 Cold-season contrast — January 2020 shadow-diff

The air-vs-soil difference is strongly **seasonal** and concentrated in the
winter hemisphere's high latitudes, where snow-insulated / frozen soil is most
decoupled from the still-swinging air. Repeating the shadow-diff for January
2020 (`fluxes_202001.nc`):

| Region | Rh amplitude ratio soil/air, July | January |
|---|---|---|
| Boreal 50-70N | 0.83 | **0.35** |
| NH temperate 25-50N | 0.96 | 0.59 |
| Tropics 25S-25N | 0.85 | 0.85 |
| SH temperate 25-50S | 0.79 | 1.06 |
| Global median | 0.86 | 0.75 |

- **Boreal winter is where the legacy driver is most wrong.** Soil-temp
  respiration diurnal amplitude there is only **35%** of the air-temp version
  (vs 83% in July): frozen / snow-insulated 0-7 cm soil barely cycles while
  2-m air still swings ~10 K. Air temp grossly overstates the winter boreal
  respiration diurnal cycle; soil temp correctly flattens it.
- **Tropics are aseasonal** (0.85 both months); **SH temperate** in austral
  summer is ~neutral-to-slightly-amplified (1.06).
- **NEE consequence is localized but real.** Globally NEE diurnal amplitude
  still moves only ~+2% (the boreal winter flux is small in the global mean),
  but because GPP ~ 0 in polar winter, NEE there is respiration alone, so the
  **boreal-January NEE diurnal amplitude is halved** (ratio 0.51, +1 h phase) —
  in exactly the high-latitude cold-season regime the rectifier is sensitive to.

Takeaway: the soil-temp driver's largest, most defensible corrections land in
the cold winter hemisphere high latitudes — the seasons/regions where an
air-temperature proxy for soil decomposition is least physical.

## 5.3 Prototype #2 result — Lloyd-Taylor temperature response

Implemented as `MICASA_RESP_TEMPFUN={q10|lloydtaylor}`, orthogonal to the
driver-variable knob (so all four combinations are selectable). `q10` + `airtemp`
is **byte-identical** to legacy — verified by `ncdiff` (max |Δ| = 0 for GPP, resp,
NEE on July 2020), not just argued from source. Under the ratio-normalization the
L&T `R_ref` and `T_ref` constants cancel, leaving only the shape
`exp(−E0/(T−T0))`; `lt.factor()` is in `lib/diurnal.r` with 5 unit tests (21 total
diurnal-transform checks pass).

Holding the driver fixed at **soil temp**, Lloyd-Taylor vs fixed Q10=1.5
changes the *respiration* diurnal amplitude substantially (ratio LT/Q10):

| Region | July | January |
|---|---|---|
| Boreal 50-70N | 1.99 | **3.70** |
| NH temperate 25-50N | 1.48 | 3.33 |
| Tropics 25S-25N | 1.49 | 1.44 |
| SH temperate 25-50S | 2.22 | 1.42 |
| Global median range | 1.51 | 1.49 |

- **LT amplifies the respiration cycle**, because its apparent Q10 (~2.5–4 across
  the relevant range) exceeds the fixed 1.5, and the amplification is **strongest
  in the cold** — LT's defining steep low-T sensitivity — reaching 3.3–3.7× in
  NH-winter high latitudes. Phase is unchanged (same temperature variable).
- **Net NEE is again negligible globally** (LT/Q10 ≈ 0.99 both months): GPP
  dominates the NEE diurnal swing. The exception is, once more, **boreal winter**
  — GPP ≈ 0 there, so NEE ≈ Rh and the boreal-January NEE diurnal amplitude moves
  by **3.6×**.

### Synthesis across both prototypes

The respiration *treatment* — both the temperature variable (air↔soil) and the
response function (Q10↔Lloyd-Taylor) — strongly controls the respiration diurnal
amplitude, and the two effects **partly oppose**: soil-temp damps, Lloyd-Taylor
re-amplifies. The July global-land respiration diurnal range spans
1.68–2.53 ×10⁻⁷ mol m⁻² s⁻¹ (×1.5) across {Q10·air, Q10·soil, LT·soil}.

But the **NEE the inversion actually sees is robust to all of it** (~1–2%
globally), because GPP redistribution dominates the NEE diurnal cycle — *except*
in the polar / boreal cold season, where GPP → 0 makes NEE track respiration
directly and the choice matters (2–3.6×).

**Recommendation.** Adopt prototype #1 (soil-temp driver) as the
default-candidate: it is the physically-correct variable, costs nothing, and
moves NEE only ~2% (in the right direction). Keep Lloyd-Taylor implemented and
**opt-in but not default**: it changes respiration amplitude a lot yet NEE almost
nowhere except boreal winter, and its steep low-T sensitivity is the most
uncertain piece. The test that would actually discriminate Q10=1.5 from
Lloyd-Taylor is the **observed diurnal amplitude of ecosystem respiration**
(eddy covariance) — a validation, not another model run — has now been run and
supports soil (§5.4(4)).

## 5.4 Decision: recommend soil-temp as the default driver

With a full-year-2019 spatial-block-bootstrap analysis (all 12 months, matched
PCHIP air-vs-soil pair), the case for flipping the default from air to soil
temperature is solid, and the eddy-covariance validation gate has now been run and
**supports soil** (4). Recommendation: **make `MICASA_RESP_DRIVER=soiltemp` the
default; keep Lloyd-Taylor opt-in**, pending sign-off.
(`fitter_diagnostics/resp_driver_blockboot.py`; committed output
`fitter_diagnostics/resp_driver_blockboot_2019.txt`.)

**(1) The effect is real, sign-correct, and robust across seasons** — *spatial*
block bootstrap (resampling unit = 10° block, so the CI respects the field's
spatial autocorrelation; B=2000; full-year 2019):

| Quantity | soil/air ratio | 95% CI (spatial block, 10°) |
|---|---|---|
| Respiration diurnal amplitude | 0.80 | [0.78, 0.83] |
| NEE diurnal amplitude | 1.023 | [1.021, 1.024] |

The NEE CI excludes 1 — and **every one of the 12 months excludes 1 individually**
(1.016–1.027), surviving a conservative 20°-block CI [1.020, 1.025] — but the change
is ~2.3%: significant *and* small, the profile of a defensible refinement. By band
(resp amplitude ratio, annual, 10° block): boreal **0.61 [0.58, 0.63]**, NH-temp
0.81 [0.78, 0.85], tropics 0.85 [0.81, 0.90], SH-temp 0.94 [0.86, 1.07]. *(The
block-bootstrap NEE CI has width 0.0032; an earlier i.i.d.-cell resample gave
[1.0225, 1.0228] — width ~0.0002, ~16× tighter — because it treated ~15.6k
autocorrelated land cells as independent. The block-bootstrap CIs above are the
correct ones; the conclusion is unchanged.)*

**(2) Forcing consistency — confirms the implementation, not the physics.** The
per-cell diurnal amplitude ratio of the *driver* (`stl1`/`t2m`, annual 2019) is
**0.80**, matching the resulting respiration amplitude ratio (0.80) — but this match
is **near-tautological**: respiration is a monotone function of its temperature
driver, so its diurnal amplitude *must* track the driver's. It confirms the code
genuinely uses soil temperature (not a tuning); it is **not** independent evidence
that soil temperature is the *correct* driver. That rests on the physics
(decomposition responds to soil, not air, temperature) and on (4).

![ERA5 forcing: 0-7cm soil temp damps & lags 2-m air](figures/resp_forcing_t2m_vs_stl1.png)

![Respiration diurnal cycle, air vs soil driver](figures/resp_diurnal_air_vs_soil.png)

![Per-cell respiration amplitude ratio distribution](figures/resp_amplitude_ratio_hist.png)

**(3) Largest, most defensible where air-temp is least physical.** The damping is
strongest in the boreal band (resp amplitude 0.61 of air, annual; lower still in
deep winter, when snow-insulated or frozen soil is most decoupled from swinging
air) — exactly where the air-temp proxy is least appropriate. Since GPP ≈ 0 there,
the NEE amplitude is correspondingly reduced.

**(4) Cost and risk.** Zero new inputs (`stl1` already loaded), every monthly
total conserved exactly, default-off byte-identical to the current product
(`fitter_diagnostics/bytecheck_resp_driver_default.txt`, max |Δ| = 0). **The
eddy-covariance validation — the one independent test that soil is the *correct*
driver — has now been run** (`fitter_diagnostics/ec_resp_driver_validation.py`,
AmeriFlux half-hourly): at night (NEE ≈ respiration, no GPP, no partitioning model),
soil temperature explains nighttime respiration better than air at **11/14 sites
(79%)** with soil-temp + nighttime flux (median R² 0.315 vs 0.295; ΔR² +0.02–0.05),
and the partitioned-RECO diurnal amplitude is damped toward soil (0.86× the air-Q10).
It **supports soil** — modestly, consistent with the small NEE effect; a fuller
FLUXNET2015 test would sharpen it. The gate is addressed and points the right way;
default flip pending sign-off.

**Lloyd-Taylor stays opt-in:** it swings respiration amplitude 1.5–3.7× but NEE
only ~1% (§5.3), and its steep low-T sensitivity is the uncertain piece — flip it
only after an eddy-covariance amplitude check.

**Implementation:** a one-line default change (`MICASA_RESP_DRIVER` default
`airtemp`→`soiltemp`), pending sign-off; until then the production default is
unchanged.

## 6. References

- Best et al. (2011), *The Joint UK Land Environment Simulator (JULES)*, GMD 4:677–699, doi:10.5194/gmd-4-677-2011.
- Davidson, Janssens & Luo (2006), *On the variability of respiration in terrestrial ecosystems: moving beyond Q10*, GCB 12:154–164, doi:10.1111/j.1365-2486.2005.01065.x.
- Denning, Fung & Randall (1995), *Latitudinal gradient of atmospheric CO2 due to seasonal exchange with land biota*, Nature 376:240–243, doi:10.1038/376240a0.
- Falge et al. (2001), *Gap filling strategies for defensible annual sums of net ecosystem exchange*, Ag. For. Meteorol. 107:43–69, doi:10.1016/S0168-1923(00)00225-2.
- Haynes et al. (2019), *Representing grasslands using dynamic prognostic phenology in SiB4*, JAMES 11:4423–4439, doi:10.1029/2018MS001540.
- Hersbach et al. (2020), *The ERA5 global reanalysis*, QJRMS 146:1999–2049, doi:10.1002/qj.3803.
- Jung et al. (2020), *Scaling carbon fluxes from eddy covariance sites to globe (FLUXCOM)*, Biogeosciences 17:1343–1365, doi:10.5194/bg-17-1343-2020.
- Lasslop et al. (2010), *Separation of net ecosystem exchange into assimilation and respiration using a light response curve approach*, GCB 16:187–208, doi:10.1111/j.1365-2486.2009.02041.x.
- Lloyd & Taylor (1994), *On the temperature dependence of soil respiration*, Functional Ecology 8(3):315–323, doi:10.2307/2389824.
- Nelson et al. (2024), *X-BASE: a global high-resolution carbon and water flux product (FLUXCOM-X)*, Biogeosciences 21:5079–5115, doi:10.5194/bg-21-5079-2024.
- Olsen & Randerson (2004), *Differences between surface and column atmospheric CO2 and implications for carbon cycle research*, JGR 109:D02301, doi:10.1029/2003JD003968.
- Potter et al. (1993), *Terrestrial ecosystem production: a process model based on global satellite and surface data (CASA)*, Global Biogeochem. Cycles 7(4):811–841, doi:10.1029/93GB02725.
- Reichstein et al. (2005), *On the separation of net ecosystem exchange into assimilation and ecosystem respiration*, GCB 11:1424–1439, doi:10.1111/j.1365-2486.2005.001002.x.
