#!/usr/bin/env Rscript
## Unit tests for lib/pchip_fit.r :: pchip.fit.cell (base R only, CI-runnable).
##
## pchip.fit.cell is the PCHIP-on-cumulative fitter core (Fritsch-Carlson
## monotone cubic Hermite) that write_pchip.r runs per grid cell. These
## checks pin its contract on synthetic monthly series, so a regression in
## the fitter is caught by CI rather than by a post-production verify_v2 run.
##
## Run:  Rscript tests/test_pchip_fit.r
## Exits non-zero on any failure.

.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "pchip_fit.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}

## piece k: f(s) = a*s^2 + b*s + c on s in [0, h_k]
peval <- function(fit, k, s) fit$a[k] * s^2 + fit$b[k] * s + fit$c[k]
## analytic integral of piece k over s in [0, h]
pint  <- function(fit, k, h) fit$a[k]*h^3/3 + fit$b[k]*h^2/2 + fit$c[k]*h

## ---- output shape ---------------------------------------------------------
fit <- pchip.fit.cell(0:5, c(1, 4, 2, 8, 3))
check("returns a list with a, b, c", all(c("a", "b", "c") %in% names(fit)))
check("coefficient vectors have length n (= length(ybar))",
      length(fit$a) == 5L && length(fit$b) == 5L && length(fit$c) == 5L)

## ---- all-zero input -> all-zero coefficients ------------------------------
z <- pchip.fit.cell(0:3, c(0, 0, 0))
check("all-zero input -> zero coefficients",
      all(z$a == 0) && all(z$b == 0) && all(z$c == 0))

## ---- constant input -> constant fit ---------------------------------------
cf <- pchip.fit.cell(0:4, rep(5, 4))
check("constant input -> quadratic term a == 0", max(abs(cf$a)) < 1e-9)
check("constant input -> linear term b == 0",    max(abs(cf$b)) < 1e-9)
check("constant input -> constant term c == value", max(abs(cf$c - 5)) < 1e-9)

## ---- integral / monthly-mean preservation (the core invariant) ------------
## The piecewise fit must integrate back to ybar[k] * h on every piece.
ybar <- c(1, 4, 2, 8, 3)
fit  <- pchip.fit.cell(0:5, ybar)
errs <- sapply(seq_along(ybar), function(k) abs(pint(fit, k, 1) - ybar[k]))
check("uniform knots: each piece integrates to ybar[k]", max(errs) < 1e-9)

## non-uniform knot spacing
x2    <- c(0, 1, 3, 4, 7)
yb2   <- c(2, 5, 1, 6)
fit2  <- pchip.fit.cell(x2, yb2)
h2    <- diff(x2)
errs2 <- sapply(seq_along(yb2), function(k) abs(pint(fit2, k, h2[k]) - yb2[k]*h2[k]))
check("non-uniform knots: each piece integrates to ybar[k]*h", max(errs2) < 1e-9)

## ---- C1 continuity at interior knots --------------------------------------
## right-limit of piece k-1 at s = h must equal c[k] (left-limit of piece k).
jumps <- sapply(2:length(ybar),
                function(k) abs(peval(fit, k - 1, 1) - fit$c[k]))
check("C1 continuous at interior knots", max(jumps) < 1e-9)

## ---- non-negativity for non-negative input (Fritsch-Carlson guarantee) ----
## Sharp seasonality -- the case where the old PIQS-quadratic overshot zero.
pos      <- pchip.fit.cell(0:5, c(0.1, 5, 0.2, 8, 0.1))
fmin_pos <- min(sapply(1:5, function(k)
  min(peval(pos, k, seq(0, 1, length.out = 50)))))
check("non-negative input -> fit >= 0 everywhere", fmin_pos > -1e-9)

## ---- sign-flip branch: all-negative input (the GPP case) ------------------
## GPP = -2*NPP is <= 0; the fitter negates, fits, then negates the result.
negbar   <- c(-0.1, -5, -0.2, -8, -0.1)
neg      <- pchip.fit.cell(0:5, negbar)
fmax_neg <- max(sapply(1:5, function(k)
  max(peval(neg, k, seq(0, 1, length.out = 50)))))
check("non-positive input -> fit <= 0 everywhere", fmax_neg < 1e-9)
errs_neg <- sapply(seq_along(negbar), function(k) abs(pint(neg, k, 1) - negbar[k]))
check("non-positive input: integral preserved through the sign-flip branch",
      max(errs_neg) < 1e-9)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall pchip_fit tests passed\n")
