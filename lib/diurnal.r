## lib/diurnal.r -- the diurnalize flux transform.
##
## diurnalize-ERA5.r turns a monthly-mean flux into an hourly one by
## scaling it with an hourly meteo driver (ERA5 shortwave for GPP, a Q10
## temperature factor for respiration) and inserting the smooth
## sub-monthly shape from the PCHIP/PIQS coefficient fit. These are the
## pure-arithmetic cores of that transform, extracted so its invariants
## can be unit-tested (tests/test_diurnal.r) without the netCDF I/O.
##
## Pure base R; no ct / ncdf4 dependency.

## Hourly flux at one time slot.
##
##   driver        hourly meteo driver at this slot
##   monthly.mean  the cell's monthly-mean flux
##   mean.driver   the driver's monthly mean (so driver / mean.driver has
##                 monthly mean 1)
##   qmod          the fitted smooth sub-monthly flux value for this slot
##
## driver * monthly.mean / mean.driver has monthly mean `monthly.mean`
## (when mean.driver is the driver's true monthly mean); subtracting
## monthly.mean and adding qmod swaps that flat monthly mean for the
## fitted sub-monthly shape. Mean-preserving: averaged over a month the
## result equals mean(qmod). Operates element-wise -- arguments may be
## scalars or conformable arrays.
diurnal.flux <- function(driver, monthly.mean, mean.driver, qmod) {
  driver * monthly.mean / mean.driver - monthly.mean + qmod
}

## Q10 temperature factor for the respiration driver.
##
##   temp.K    temperature in Kelvin (2-m air t2m, or 0-7cm soil stl1)
##   q10.base  factor per 10 K (default 1.5, the legacy Olsen & Randerson value)
##   ref.K     reference temperature (default 273.15 K = 0 degC)
##
## Returns q10.base ^ ((temp.K - ref.K) / 10). Pure element-wise; temp.K may be
## a scalar or a conformable array. The *choice of temperature variable* (air
## vs soil) is made by the caller -- see diurnalize-ERA5.r MICASA_RESP_DRIVER.
## Because the diurnalize transform normalizes by this factor's monthly mean,
## the diurnal SHAPE depends on the nonlinearity (q10.base) and the driving
## variable, but not on ref.K (a constant ref shifts q10 by a constant factor
## that cancels in q10(t)/mean(q10)).
q10.factor <- function(temp.K, q10.base = 1.5, ref.K = 273.15) {
  q10.base ^ ((temp.K - ref.K) / 10.0)
}

## Lloyd & Taylor (1994) temperature response for the respiration driver.
##
##   temp.K  temperature in Kelvin (air t2m or soil stl1)
##   E0      activation-energy parameter (K); 308.56 is the L&T 1994 value
##   T0      lower temperature limit (K); 227.13 = -46.02 degC (L&T 1994)
##
## The full L&T form is R = R_ref * exp(E0*(1/(T_ref-T0) - 1/(T-T0))). Under the
## diurnalize ratio-normalization f(t)/mean(f) the R_ref and the constant
## exp(E0/(T_ref-T0)) cancel, so only the temperature-dependent SHAPE survives:
## exp(-E0/(T-T0)). Its apparent Q10 RISES as T falls (steeper low-temperature
## sensitivity than a fixed Q10), which is the point of using it. temp.K is
## clamped just above T0 so the curve stays finite and monotone for any input
## (below ~T0 respiration underflows to ~0, i.e. frozen). Element-wise.
lt.factor <- function(temp.K, E0 = 308.56, T0 = 227.13) {
  exp(-E0 / (pmax(temp.K, T0 + 0.1) - T0))
}

## Polar-night clip: GPP must be zero where the shortwave driver is zero
## (no incoming light => no photosynthesis). Without this the sub-monthly
## quadratic leaks a small residual into dark cell-hours. `gpp` and
## `driver` must be the same shape.
polar.night.clip <- function(gpp, driver) {
  gpp[driver == 0] <- 0
  gpp
}

## Mass-conserving alternative to the plain polar-night clip. After dark-hour GPP
## has been zeroed (per-slot, above), redistribute the clipped uptake onto each
## cell's remaining LIT hours so the monthly mean is restored to gpp.mn -- a uniform
## per-cell rescale, which preserves the ssrd-proportional diurnal shape. Cells with
## no lit hours all month (full polar night) have nothing to redistribute onto and
## stay zeroed; their monthly-mean GPP should be ~0 there anyway, so the residual
## removed is the fit's spurious dark-hour leak. This is the V2 default
## (MICASA_POLAR_CLIP=conserve); MICASA_POLAR_CLIP=plain falls back to the legacy
## byte-identical zero-clip (polar.night.clip alone).
##   gpp    [lat, lon, nslot] hourly GPP, already dark-clipped
##   gpp.mn [lat, lon]        target monthly-mean GPP
polar.night.renorm <- function(gpp, gpp.mn) {
  nslot         <- dim(gpp)[3]
  clipped.total <- rowSums(gpp, dims = 2)        # per-cell sum over slots (lit only)
  target.total  <- gpp.mn * nslot                # intended per-cell total
  scale         <- target.total / clipped.total
  ## rescale only where there is lit-hour mass of the matching sign; else no-op
  scale[!is.finite(scale) | clipped.total == 0 |
        sign(clipped.total) != sign(target.total)] <- 1
  sweep(gpp, c(1, 2), scale, "*")
}
