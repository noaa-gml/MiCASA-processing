#!/usr/bin/env Rscript
## Unit tests for lib/mss_fit.r :: mss.fit.cell -- the QP (monotone
## smoothing-spline) fitter that write_mss.r runs per grid cell.
##
## SKIPPED (exit 0) where the quadprog package is unavailable: without it
## mss.fit.cell silently falls back to PCHIP, which would make these checks
## test the wrong fitter. quadprog is installed on Orion (write_mss.r needs
## it); a CI image without quadprog skips cleanly.
##
## Run:  Rscript tests/test_mss_fit.r
## Exits non-zero on any failure.

.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))

if (!requireNamespace("quadprog", quietly = TRUE)) {
  cat("SKIP: quadprog package not available -- MSS QP fitter not testable here\n")
  quit(status = 0L)
}
source(file.path(.repo, "lib", "mss_fit.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}

## piece k: f(s) = a*s^2 + b*s + c on s in [0, h_k]
peval <- function(fit, k, s) fit$a[k] * s^2 + fit$b[k] * s + fit$c[k]
pint  <- function(fit, k, h) fit$a[k]*h^3/3 + fit$b[k]*h^2/2 + fit$c[k]*h

## ---- setup ---------------------------------------------------------------
x  <- 0:6                                   # 6 unit segments
st <- mss.fit.setup(x)
check("setup: h.seg has one entry per segment", length(st$h.seg) == 6L)
check("setup: n.var = n.seg + 1", st$n.var == st$n.seg + 1L)
check("setup: smoothness Hessian G is n.var x n.var",
      all(dim(st$G) == c(st$n.var, st$n.var)))

## ---- output shape & all-zero --------------------------------------------
fit <- mss.fit.cell(x, c(1, 4, 2, 8, 3, 5), st)
check("mss.fit.cell returns a, b, c each of length n.seg",
      length(fit$a) == 6L && length(fit$b) == 6L && length(fit$c) == 6L)
z <- mss.fit.cell(x, rep(0, 6), st)
check("all-zero input -> zero coefficients",
      all(z$a == 0) && all(z$b == 0) && all(z$c == 0))

## ---- integral / monthly-mean preservation -------------------------------
## The Q/L/K coefficient formulas preserve int = ybar[k]*h on every piece,
## whatever knot slopes the QP returns.
ybar <- c(1, 4, 2, 8, 3, 5)
fit  <- mss.fit.cell(x, ybar, st)
errs <- sapply(seq_along(ybar), function(k) abs(pint(fit, k, 1) - ybar[k]))
check("each piece integrates back to ybar[k]", max(errs) < 1e-8)

x2    <- c(0, 1, 3, 4, 7, 9)                # non-uniform knots
st2   <- mss.fit.setup(x2)
yb2   <- c(2, 5, 1, 6, 3)
fit2  <- mss.fit.cell(x2, yb2, st2)
h2    <- diff(x2)
errs2 <- sapply(seq_along(yb2),
                function(k) abs(pint(fit2, k, h2[k]) - yb2[k]*h2[k]))
check("non-uniform knots: each piece integrates to ybar[k]*h",
      max(errs2) < 1e-8)

## ---- non-negativity at the QP test points -------------------------------
## The QP enforces f >= 0 at n.test.per.segment points per segment. Sharp
## seasonality -- the case the MSS fitter exists for.
ts   <- (seq_len(st$n.test.per.segment) - 0.5) / st$n.test.per.segment
pos  <- mss.fit.cell(x, c(0.1, 6, 0.2, 9, 0.1, 4), st)
fmin <- min(sapply(seq_len(st$n.seg), function(k) min(peval(pos, k, ts))))
check("non-negative input -> fit >= 0 at every QP test point", fmin > -1e-7)

## ---- sign-flip branch: all-negative input (the GPP case) ----------------
negbar <- c(-0.1, -6, -0.2, -9, -0.1, -4)
neg    <- mss.fit.cell(x, negbar, st)
fmax   <- max(sapply(seq_len(st$n.seg), function(k) max(peval(neg, k, ts))))
check("non-positive input -> fit <= 0 at every QP test point", fmax < 1e-7)
errsn  <- sapply(seq_along(negbar), function(k) abs(pint(neg, k, 1) - negbar[k]))
check("non-positive input: integral preserved through the sign-flip branch",
      max(errsn) < 1e-8)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall mss_fit tests passed\n")
