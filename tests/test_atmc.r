#!/usr/bin/env Rscript
## test_atmc.r -- assertions about the ATMC removal (lib/diurnal.r ::
## atmc.removal, driven by MICASA_ATMC in diurnalize-ERA5.r).
##
## Usage:  Rscript tests/test_atmc.r
##
## ATMC is MiCASA's GMAO atmospheric-closure term. The product defines
## NEE = Rh - NPP - ATMC; the CT prep ingests Rh - NPP, so the term is dropped.
## MICASA_ATMC subtracts it back.
##
## The property that matters is MEAN PRESERVATION. diurnalize's hourly flux has
## monthly mean mean(qmod) -- taken from the piqs fit, not from the monthly
## array -- so a removal applied to the monthly array alone would change the
## diurnal amplitude and leave the budget untouched: a complete, plausible,
## wrong set of files. Subtracting a mean-preserving hourly field instead drops
## the delivered mean by exactly ATMC, with no refit and no dependence on the
## fit at all.
##
##   MEAN      both modes average to atmc.mn over the month, cell by cell
##   SHAPE     q10 mode reproduces the respiration channel's own weighting,
##             so the removal cannot drift from the flux it is subtracted from
##   FLAT      flat mode is uniform across the day
##   MODES     an unknown mode is refused rather than silently defaulted
##   FALSIFY   a plausible WRONG implementation fails the same MEAN assertion,
##             proving the check has teeth
suppressWarnings(source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "diurnal.r")))

fails <- 0L
ok <- function(label, cond, detail = "") {
  cat(sprintf("%-58s %s%s\n", label, if (isTRUE(cond)) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!isTRUE(cond)) fails <<- fails + 1L
}

## A synthetic cell-month: 24 slots of a realistic diurnal temperature swing,
## plus a cold cell where the Q10 factor is small and a warm one where it is not.
set.seed(1)
nslot <- 24L
hours <- seq_len(nslot) - 0.5
tempK <- cbind(cold = 263 + 6 * sin(2 * pi * (hours - 9) / 24),
               warm = 298 + 9 * sin(2 * pi * (hours - 9) / 24))
q10   <- q10.factor(tempK)                    # [slot, cell]
q10mn <- colMeans(q10)
atmc  <- c(cold = 1.7e-7, warm = 4.2e-7)      # mol m-2 s-1, monthly mean

apply.mode <- function(mode, driver = q10, mean.driver = q10mn) {
  t(vapply(seq_len(nslot),
           function(i) atmc.removal(mode, atmc, driver[i, ], mean.driver),
           numeric(length(atmc))))           # [slot, cell]
}

## ---- MEAN: both modes integrate to exactly the monthly ATMC ---------------
for (mode in c("q10", "flat")) {
  applied <- apply.mode(mode)
  err <- max(abs(colMeans(applied) - atmc))
  ok(sprintf("MEAN     %-4s averages to atmc.mn", mode),
     err < 1e-18, sprintf("max|d| = %.2e", err))
}

## ---- SHAPE: q10 mode carries the respiration channel's own weighting ------
## The respiration channel is weighted by q10/q10.mn; so is the removal, using
## the same arrays. Ratio removal/atmc must equal that weighting exactly.
applied <- apply.mode("q10")
wexp <- sweep(q10, 2, q10mn, "/")
ok("SHAPE    q10 removal matches the resp weighting",
   max(abs(sweep(applied, 2, atmc, "/") - wexp)) < 1e-12)
ok("SHAPE    q10 removal is not flat (it has a diurnal cycle)",
   min(apply(sweep(applied, 2, atmc, "/"), 2, function(x) max(x) - min(x))) > 0.1)

## ---- FLAT: uniform across the day ----------------------------------------
applied <- apply.mode("flat")
ok("FLAT     flat removal is constant over the day",
   max(apply(applied, 2, function(x) diff(range(x)))) == 0)

## ---- MODES: refuse the unknown -------------------------------------------
ok("MODES    unknown mode is an error",
   inherits(try(atmc.removal("lloyd", atmc, q10[1, ], q10mn), silent = TRUE), "try-error"))
ok("MODES    q10 without a driver is an error",
   inherits(try(atmc.removal("q10", atmc), silent = TRUE), "try-error"))

## ---- FALSIFY: the MEAN assertion must be able to fail ---------------------
## The plausible mistake is normalising by something other than the driver's
## own monthly mean -- e.g. by the mean over land cells, or by a stale value
## from another month. That still produces a complete diurnal field, and it is
## exactly what the MEAN check exists to catch.
wrong <- apply.mode("q10", mean.driver = rep(mean(q10mn), length(atmc)))
werr <- max(abs(colMeans(wrong) - atmc))
ok("FALSIFY  mis-normalised removal FAILS the MEAN check",
   werr > 1e-12, sprintf("max|d| = %.2e (must be large)", werr))

cat(sprintf("\n%s: %d check(s) failed\n", if (fails == 0L) "OK" else "FAILURES", fails))
quit(status = if (fails == 0L) 0L else 1L)
