# The climatology prior

*Climatological fluxes, real meteorology.*

A MiCASA product whose monthly **magnitude** is a multi-year climatology, but
whose day-to-day and hour-to-hour **structure** comes from each target day's own
ERA5 meteorology. Interannual variability and trend are removed by
construction; weather is not.

Requested by Andy Jacobson, 2026-08-10 (CT2026 Issues):

> use the 2000-2020 climatology as the prior. This means computing the mean
> (climatology) from monthly data, then running it through the diurnalization
> process. I.e. climatological fluxes but diurnalized using real meteorology.

## Why you might want one

The MiCASA bio flux carries a strong trend over 2021-2025: the global land net
(NEE = Rh − NPP) weakens from −2.18 to −0.42 PgC/yr, a drift of about
+0.4 PgC/yr per year. Used as an inversion prior, that trend is a *growing*
prior error, and a flux-only filter with no state for accumulated burden cannot
remove it. A climatology prior holds the annual magnitude fixed at the
long-term mean, so the prior error stops growing, while keeping the
weather-driven sub-monthly structure an inversion needs to place flux in time.

It is a **diagnostic / counterfactual** product, not an observational one. Say
so whenever it is used.

## Running it

```bash
./run_climatology_prior.sh \
    --outdir    /work2/noaa/co2/ash/micasa-clim \
    --source    /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc \
    --baseline  2001-2020 \
    --years     2021-2025 \
    --reference /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/ERA5
```

The driver blocks on its SLURM stages, so run it from a batch job, not a login
shell you are about to close. `--dry-run` prints the plan and touches nothing.

Stages:

| # | Stage | What it does |
|---|-------|--------------|
| 1 | Series | `make_climatology_series.py` — per-month + concatenated climatological monthly files |
| 2 | Fit | the sub-monthly smoother (`write_pchip.r` by default), **rebuilt on that series** |
| 3 | Diurnalize | `climatology_diurnalize.sbatch`, one array task per year, real ERA5 meteo |
| 4 | Day-split | `climatology_daysplit.sbatch` — stamp provenance, split to per-day NEE |
| 5 | Verify | `tests/verify_climatology_prior.py` |

Output layout mirrors an ordinary MiCASA tree, so a consumer needs only a path
change:

```
<outdir>/monthly_1x1/MiCASA_v1_flux_x360_y180_monthly[_YYYYMM].nc
<outdir>/fit.piqs.rda
<outdir>/ERA5/fluxes_YYYYMM.nc
<outdir>/ERA5/MiCASA_v1.nee.YYYYMMDD.nc      <- what CarbonTracker reads
<outdir>/.micasa-climatology-prior           <- marker; the driver refuses to
                                                write into a tree without it
```

## ⚠ The trap: the fit sets the monthly mean, not the monthly-mean array

**This is the one thing to understand before modifying any of this.**

`diurnalize-ERA5.r` builds its hourly flux as (`lib/diurnal.r::diurnal.flux`):

```
f(t) = driver(t) * mean / mean_driver  -  mean  +  qmod(t)
```

The first two terms form a zero-mean hourly anomaly whose amplitude is set by
the monthly-mean array. The **monthly mean of `f` is `mean(qmod)`** — and
`qmod` is the piecewise quadratic evaluated from `fit.piqs.rda`. PCHIP-on-
cumulative is exactly mean-preserving, so the delivered monthly mean is the
monthly mean *the fit was built from*, and the monthly-mean array only sets the
diurnal amplitude.

Consequently:

> Feeding climatological monthly means to a fit built on the **real** monthly
> series reinstates the real interannual signal through `qmod`, while producing
> a complete and entirely healthy-looking set of output files.

Nothing errors. Every file is present, every artefact count is right, every
header is correct. The product is simply not a climatology. That is why stage 2
exists and why `tests/verify_climatology_fit.r` asserts, per delivered month,
that the fit's analytic month-mean reproduces the climatological mean *and*
that the same calendar month carries identical coefficients across years.

The negative control is worth keeping in mind as the scale of the effect: run
that fit check against the production `fit.piqs.rda` and the flatness statistic
is ~2.6 (order one, i.e. wholly interannual); against a climatological fit it
is ~7e-14 (floating-point round-off in `cumsum`). Fourteen orders of magnitude.

## Scope: only bio is climatologized

The **hourly** product is bio only: `NEE = Rh − NPP`.

The **monthly** product carries the full MiCASA variable set — but only NPP and
Rh are climatological. FIRE, FUEL and ATMC are copied **unchanged** from the
source archive and remain **real, year-specific** fields
(`--carry-imposed`, default `FIRE,FUEL,ATMC`; pass `''` for a strictly bio-only
product).

That mix is deliberate and matches reality: CarbonTracker takes bio from this
product and imposes fire/fuel *separately* from the unmodified daily tree, so a
file combining climatological bio with real fire is exactly the prior in use.
It is also the obvious way to mislead someone, so the distinction is
machine-readable:

```
NPP:climatology  = "yes"     Rh:climatology   = "yes"
FIRE:climatology = "no"      FUEL, ATMC       = "no"
:carried_imposed_variables = "ATMC,FIRE,FUEL"
```

`verify_climatology_prior.py` check 1.7 asserts **both halves** of that claim —
the climatological variables identical in every year, the carried ones not —
reading the file's own labels rather than a hardcoded list. A carried variable
that turned out identical across years (accidentally climatologized), or a
climatological one that varied, fails.

Two edges to know: months past the end of the source record carry no imposed
variables, and the **concatenation carries none at all** — it is the fitter's
input, which reads only NPP and Rh, and a partially populated FIRE there would
read as real data with silent holes.

### Why fire is not climatologized

CarbonTracker reads bio and fire through **two separate rc keys pointing at two
separate directories** (`ct.emissions.rc`):

```
ct.bio       == m  ->  bio.input.dir   : .../MiCASA.0/ERA5        prefix MiCASA_v1.nee
ct.wildfire  == m  ->  fires.input.dir : .../MiCASA.0/daily_1x1   vars FIRE, FUEL
```

So adopting the climatology prior is a **one-key change** — repoint
`bio.input.dir` — and fire/fuel keep reading the unmodified `daily_1x1` tree.
A climatological fire channel would be unused by that configuration and
misleading to anyone who found it. Fire also has its own separate history (the
2025 vNRT spike and its remediation) that a climatology would silently paper
over. `fire_flux_imp` is imposed, never optimized, so it is not part of the
prior-error problem this product addresses.

If you ever *do* want a climatological fire channel, it is a deliberate,
separate decision — not a default of this pipeline.

## Baseline availability

The MiCASA 1-degree monthly record **starts in 2001-01**. A request for a
"2000-2020" climatology therefore resolves to **2001-2020** (20 years).
`make_climatology_series.py` warns when a requested baseline year is absent and
records both the requested and the actual window in the file attributes
(`climatology_baseline`, `climatology_baseline_requested`), so the product can
never silently claim a window it did not use.

## Leap years and February

The climatology is a mean of monthly **rates** (gC m⁻² s⁻¹), so February
averages 28- and 29-day Februaries together as rates. A leap-year February then
carries the climatological February rate across 29 days — slightly more mass
than a common year, which is the physically sensible reading and needs no
special case. 29 February is a genuine day driven by its own meteorology, not a
copy of 28 February; `verify_climatology_prior.py` check 3.f asserts exactly
that.

The knock-on effect is visible and small: with a 2001-2020 baseline the
delivered annual net is −2.4756 PgC/yr in common years and −2.4633 in 2024,
because the extra February day is a weak-sink month.

## Padding the series

The default span (`--span 2020-01:2026-12`) is a year wider either side of the
delivered window. The PCHIP slope at a knot depends on its neighbours, so
without padding the first and last delivered months would inherit an edge slope
from real data. Check 1.3 warns if the span is not padded.

## Verification

```bash
python3 tests/verify_climatology_prior.py --selftest      # gates prove they can fail
python3 tests/verify_climatology_prior.py \
    --outdir /work2/noaa/co2/ash/micasa-clim --years 2021-2025 \
    --reference /work2/noaa/co2/GFED-CASA/2026/MiCASA.0/ERA5
```

Unlike `tests/verify_v2.py`, this battery **exits non-zero on failure**. The
product looks exactly like a real one from the outside, so a verdict that is
printed but not wired to an exit status would be a comment.

Sections: (1) the series is genuinely climatological over the recorded
baseline; (2) the fit spans every delivered month, is mean-preserving, and is
flat across years; (3) completeness — file counts, no zero-byte files, 24 slots
each, time axis at 00:30–23:30 UTC, real leap days; (4) format parity with the
reference product, with provenance attributes allowed to differ and everything
else required identical; (5) round-trip (produced dailies average back to the
climatology), trend removal against the raw product measured with the *same*
instrument, and liveness (the arms actually differ).

Unit tests, no cluster needed:

```bash
python3 tests/test_climatology.py
python3 tests/test_budget.py
```

## Configuration reference

Set in `config.sh`; the driver's flags override them.

| Variable | Default | Meaning |
|---|---|---|
| `MICASA_CLIM_BASELINE_START` | 2001 | first baseline year |
| `MICASA_CLIM_BASELINE_END` | 2020 | last baseline year |
| `MICASA_CLIM_SPAN_START` | 2020-01 | first month of the synthetic series |
| `MICASA_CLIM_SPAN_END` | 2026-12 | last month of the synthetic series |
| `MICASA_CLIM_YEAR_START` | 2021 | first delivered year |
| `MICASA_CLIM_YEAR_END` | 2025 | last delivered year |
| `MICASA_CLIM_SOURCE` | *(driver-derived)* | concatenated monthly file to average |

The diurnalization runs at production settings — `MICASA_RESP_DRIVER=airtemp`,
`MICASA_RESP_TEMPFUN=q10`, `MICASA_POLAR_CLIP=conserve`, `--fitter pchip` —
stated explicitly by the driver rather than inherited, so the climatology is
the only difference from the ordinary product. `MICASA_STRICT_PIQS=1` is forced
on: a delivered month falling outside the fit window would otherwise silently
take coefficient-climatology instead of its own coefficients.

## Related

- `docs/DIURNALIZATION_ALTERNATIVES.md` — the diurnal transform and driver choices
- `docs/V1_TO_V2_JUSTIFICATION.md` §5.2 — why ATMC is not subtracted
- `compute_clim.py` — the *other* climatology: the NPP/Rh day-of-year fallback
  `diurnalize-ERA5.r` uses for months with no published monthly file. Different
  purpose, shared core in `lib/climatology.py`.
