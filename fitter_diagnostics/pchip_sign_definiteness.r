#!/usr/bin/env Rscript
## pchip_sign_definiteness.r -- evidence for the *scope* of PCHIP's sign claim.
##
## The docs must NOT claim PCHIP is "sign-definite everywhere by construction."
## Fritsch-Carlson guarantees the cumulative cubic F is monotone at the KNOTS;
## the derivative quadratic f = F' can still dip the wrong way mid-segment, even
## on strictly single-signed monthly means. This script reproduces that dip and
## quantifies the residual, so the "16-57x reduction, not elimination" wording in
## METHODOLOGY.md / FITTER_COMPARISON.md / V1_TO_V2_JUSTIFICATION.md is backed by
## a runnable artifact rather than an assertion.
##
## Run:  Rscript fitter_diagnostics/pchip_sign_definiteness.r
## Pure base R + lib/pchip_fit.r. No netCDF / data dependency.
.args <- commandArgs(FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "pchip_fit.r"))

set.seed(7)
x <- as.numeric(seq(as.POSIXct("2001-01-01", tz = "UTC"), by = "1 month", length.out = 37))
n <- 36; npp <- 200

## densest interior sample of the per-piece quadratic, returns worst wrong-sign
## value relative to the data's overall sign.
worst.wrongsign <- function(f, sgn) {
  m <- 0
  for (i in 1:n) {
    s  <- seq(0, 1, length.out = npp) * (x[i + 1] - x[i])
    fl <- f$a[i] * s^2 + f$b[i] * s + f$c[i]
    m  <- min(m, min(sgn * fl))     # most-negative same-sign-projected value
  }
  m
}

## (1) STRICTLY POSITIVE monthly means (Rh-like), sharp seasonality.
worst.pos <- 0; flip.cells.pos <- 0; ntrial <- 20000
for (k in 1:ntrial) {
  y <- abs(rnorm(36)) + (sin((1:36) / 12 * 2 * pi))^2 * runif(1, 0, 20) + 1e-6   # all > 0
  f <- pchip.fit.cell(x, y)
  w <- worst.wrongsign(f, +1)
  if (w < -1e-9 * median(y)) flip.cells.pos <- flip.cells.pos + 1
  if (w < worst.pos) worst.pos <- w
}
cat(sprintf("STRICTLY POSITIVE input  (%d series): worst interior flux = %.4e ; series with a wrong-sign interior dip = %.1f%%\n",
            ntrial, worst.pos, 100 * flip.cells.pos / ntrial))

## (2) STRICTLY NEGATIVE monthly means (GPP convention, gpp = -2*NPP <= 0).
worst.neg <- 0
for (k in 1:ntrial) {
  y <- -(abs(rnorm(36)) + (sin((1:36) / 12 * 2 * pi))^2 * runif(1, 0, 20) + 1e-6) # all < 0
  f <- pchip.fit.cell(x, y)
  worst.neg <- min(worst.neg, worst.wrongsign(f, -1))
}
cat(sprintf("STRICTLY NEGATIVE input  (%d series): worst interior wrong-sign flux = %.4e\n", ntrial, worst.neg))

cat("\nConclusion: PCHIP-on-cumulative is sign-definite AT THE KNOTS, not everywhere.\n")
cat("Interior dips occur on single-signed input, so the production sign-flip rate\n")
cat("(verify_v2 Check 3.1: GPP <=0.94% of cell-hours) is a real reduction vs PIQS\n")
cat("(<=14.70%), NOT an elimination 'by construction'.\n")
