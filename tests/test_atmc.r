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

## ---- NEGATIVITY: the diagnostic that can see what the other one cannot ----
## resp.negativity (lib/diurnal.r) exists because diurnalize's own
## `PIQS sign-flip [resp < 0]` counter tests qmod.resp -- the FITTED monthly
## value, evaluated before the ATMC subtraction -- and so printed identical
## lines for the off, q10 and flat arms of 2021-07 despite different physics.
## These checks assert the replacement is not blind in the same way.
resp0 <- array(c(3e-7, 1e-8, -2e-9, 5e-7), dim = c(2, 2, 4))
r0 <- resp.negativity(resp0)
ok("NEG      counts negatives in the delivered field",
   r0$n == 4L && r0$total == 16L && abs(r0$frac - 0.25) < 1e-12,
   sprintf("n=%d total=%d", r0$n, r0$total))
ok("NEG      reports the worst value and the field scale",
   isTRUE(all.equal(r0$worst, -2e-9)) && isTRUE(all.equal(r0$scale, 5e-7)))

clean <- resp.negativity(array(abs(rnorm(64)) + 1e-9, dim = c(4, 4, 4)))
ok("NEG      an all-positive field reports nothing",
   clean$n == 0L && clean$worst == 0 && clean$worst.frac == 0)

withna <- resp0; withna[1, 1, 1] <- NA; withna[2, 2, 4] <- NaN
rna <- resp.negativity(withna)
ok("NEG      NA/NaN excluded from numerator AND denominator",
   rna$total == 14L && rna$n == 4L, sprintf("total=%d n=%d", rna$total, rna$n))

## The load-bearing one: subtracting a mean-preserving removal must MOVE this
## number. A synthetic cell whose hourly respiration grazes zero is exactly the
## situation the real field is in -- 1,156 cell-hours were already negative in
## 2021-07 with no removal at all.
set.seed(2)
base <- array(rep(2e-7 * (1 + 0.98 * sin(2 * pi * (seq_len(nslot) - 0.5) / nslot)),
                  each = 1), dim = c(1, 1, nslot))
mn   <- mean(base)
w    <- q10[, "warm"]                       # a real diurnal weighting
removal <- array(vapply(seq_len(nslot),
                        function(i) atmc.removal("q10", 0.35 * mn, w[i], mean(w)),
                        numeric(1)), dim = c(1, 1, nslot))
before <- resp.negativity(base)
after  <- resp.negativity(base - removal)
ok("NEG      a mean-preserving removal CHANGES the count",
   after$n > before$n,
   sprintf("%d -> %d negative cell-hours", before$n, after$n))
## ...while the monthly mean -- what the old counter tests -- barely moves and
## never changes sign. That contrast IS the bug this replaces.
ok("NEG      ...where the monthly-mean test sees nothing",
   mean(base) > 0 && mean(base - removal) > 0,
   sprintf("monthly mean %.3e -> %.3e, still positive",
           mean(base), mean(base - removal)))

cat(sprintf("\n%s: %d check(s) failed\n", if (fails == 0L) "OK" else "FAILURES", fails))
quit(status = if (fails == 0L) 0L else 1L)
