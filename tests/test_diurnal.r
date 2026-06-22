#!/usr/bin/env Rscript
## Unit tests for lib/diurnal.r (base R only, CI-runnable).
##
## diurnal.flux and polar.night.clip are the pure cores of the
## diurnalize-ERA5.r transform. These checks pin the transform's
## invariants -- monthly-mean preservation, driver proportionality,
## polar-night zeroing -- on synthetic data.
##
## Run:  Rscript tests/test_diurnal.r
## Exits non-zero on any failure.

.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "diurnal.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}
close.all <- function(a, b, tol = 1e-9) max(abs(a - b)) <= tol

## ---- diurnal.flux: monthly-mean preservation -----------------------------
## With mean.driver set to the driver's true monthly mean, the hourly flux
## averages to mean(qmod) -- the fitted sub-monthly shape carries the mean.
set.seed(1)
driver <- runif(720, 0, 900)            # a month of hourly shortwave-like values
qmod   <- rnorm(720, mean = 4, sd = 1)
mn     <- 4.0
flux   <- diurnal.flux(driver, mn, mean(driver), qmod)
check("mean-preserving: mean(flux) == mean(qmod)",
      abs(mean(flux) - mean(qmod)) < 1e-9)

## a flat qmod == the monthly mean -> the hourly flux averages to it
flux.flat <- diurnal.flux(driver, mn, mean(driver), rep(mn, length(driver)))
check("flat qmod == monthly mean -> mean(flux) == monthly mean",
      abs(mean(flux.flat) - mn) < 1e-9)

## ---- diurnal.flux: driver proportionality --------------------------------
## With qmod flat at the monthly mean, flux = driver * mn / mean.driver,
## i.e. the hourly flux tracks the driver exactly.
prop <- diurnal.flux(driver, mn, mean(driver), rep(mn, length(driver)))
check("flat qmod -> flux is proportional to the driver",
      close.all(prop, driver * mn / mean(driver)))

## a constant driver equal to mean.driver -> the flux is exactly qmod
check("constant driver == mean.driver -> flux == qmod",
      close.all(diurnal.flux(rep(7, 10), mn, 7, qmod[1:10]), qmod[1:10]))

## ---- diurnal.flux: edge cases --------------------------------------------
check("zero monthly mean -> flux == qmod",
      close.all(diurnal.flux(driver, 0, mean(driver), qmod), qmod))

## negative monthly mean (the GPP convention: gpp = -2*NPP <= 0)
fg <- diurnal.flux(driver, -4, mean(driver), -qmod)
check("negative monthly mean is still mean-preserving",
      abs(mean(fg) - mean(-qmod)) < 1e-9)

## operates element-wise on a matrix (a grid slice)
dm <- matrix(runif(12, 1, 5), 3, 4)
qm <- matrix(rnorm(12), 3, 4)
fm <- diurnal.flux(dm, 2.0, 3.0, qm)
check("operates element-wise on a matrix",
      is.matrix(fm) && all(dim(fm) == c(3, 4)) &&
      close.all(fm, dm * 2.0 / 3.0 - 2.0 + qm))

## ---- polar.night.clip ----------------------------------------------------
g <- c(1, 2, 3, 4, 5)
check("clip zeros gpp where the driver is 0, leaves the rest",
      identical(polar.night.clip(g, c(0, 9, 0, 9, 0)), c(0, 2, 0, 4, 0)))
check("clip leaves gpp untouched when the driver is non-zero everywhere",
      identical(polar.night.clip(g, rep(1, 5)), g))
check("clip zeros everything when the driver is 0 everywhere",
      all(polar.night.clip(g, rep(0, 5)) == 0))
check("clip operates on a matrix",
      identical(polar.night.clip(matrix(c(1, 2, 3, 4), 2, 2),
                                 matrix(c(0, 5, 5, 0), 2, 2)),
                matrix(c(0, 2, 3, 0), 2, 2)))

## ---- q10.factor (respiration temperature driver) -------------------------
## Legacy formula: q10.factor(T) must reproduce 1.5^((T-273.15)/10) exactly,
## so the q10 + airtemp combination is byte-identical to the pre-prototype code
## (V2 now defaults the driver to soiltemp; airtemp remains selectable).
Tk <- 273.15 + c(-20, -5, 0, 10, 15, 25, 35)
check("q10.factor reproduces legacy 1.5^((T-273.15)/10)",
      close.all(q10.factor(Tk), 1.5 ^ ((Tk - 273.15) / 10.0), tol = 1e-12))
check("q10.factor == 1 at the reference temperature",
      abs(q10.factor(273.15) - 1) < 1e-12)
check("q10.factor rises a factor q10.base per +10 K",
      abs(q10.factor(283.15) / q10.factor(273.15) - 1.5) < 1e-12)
## ref.K cancels under the diurnalize normalization q10(t)/mean(q10): two
## reference choices give factors differing only by a constant, so the
## normalized diurnal shape is identical.
qa <- q10.factor(Tk, ref.K = 273.15); qb <- q10.factor(Tk, ref.K = 280.0)
check("ref.K cancels in the normalized ratio q10(t)/mean(q10)",
      close.all(qa / mean(qa), qb / mean(qb), tol = 1e-12))
## operates element-wise on a grid slice
Tm <- matrix(273.15 + runif(12, -10, 30), 3, 4)
check("q10.factor operates element-wise on a matrix",
      is.matrix(q10.factor(Tm)) && all(dim(q10.factor(Tm)) == c(3, 4)))

## ---- lt.factor (Lloyd & Taylor 1994 respiration response) -----------------
## Monotonic increasing in temperature.
Tk2 <- 273.15 + c(-30, -10, 0, 10, 20, 30)
check("lt.factor strictly increases with temperature",
      all(diff(lt.factor(Tk2)) > 0))
## Apparent Q10 RISES as T falls: the +10 K ratio is larger in the cold.
cold <- lt.factor(283.15) / lt.factor(273.15)
warm <- lt.factor(303.15) / lt.factor(293.15)
check("lt.factor apparent Q10 is higher at low T than high T", cold > warm)
## Finite and ~0 below the T0 frozen limit (no NaN/Inf from the clamp).
fr <- lt.factor(227.13 + c(-50, -10, -0.05))
check("lt.factor stays finite and ~0 below T0 (frozen clamp)",
      all(is.finite(fr)) && all(fr >= 0) && max(fr) < 1e-30)
## ref/scale constants cancel under the diurnalize normalization f(t)/mean(f):
## an overall multiplicative constant leaves the normalized shape unchanged.
lt1 <- lt.factor(Tk2); lt2 <- 7.3 * lt1
check("lt.factor: constant prefactor cancels in f/mean(f)",
      close.all(lt1/mean(lt1), lt2/mean(lt2), tol = 1e-12))
check("lt.factor operates element-wise on a matrix",
      is.matrix(lt.factor(Tm)) && all(dim(lt.factor(Tm)) == c(3, 4)))

## ---- polar.night.renorm: mass-conserving polar-night clip (opt-in) ------------
.ns <- 24L; .nd <- 8L
.g0 <- array(-abs(seq_len(2L*3L*.ns)) * 1e-9, c(2, 3, .ns))   # GPP < 0 (uptake)
.mn <- apply(.g0, c(1, 2), mean)                              # true monthly mean (pre-clip)
.gc <- .g0; .gc[, , 1:.nd] <- 0                               # dark-clip first nd hours
.gr <- polar.night.renorm(.gc, .mn)
check("polar.night.renorm restores the monthly mean (mass-conserving)",
      close.all(apply(.gr, c(1, 2), mean), .mn, tol = 1e-12))
check("polar.night.renorm keeps dark hours zero", all(.gr[, , 1:.nd] == 0))
.gd <- array(0, c(1, 1, .ns)); .md <- array(-1e-9, c(1, 1))   # full polar night, all clipped
check("polar.night.renorm leaves full-dark cells zeroed (no blow-up)",
      all(polar.night.renorm(.gd, .md) == 0))

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall diurnal transform tests passed\n")
