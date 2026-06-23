#!/usr/bin/env Rscript

## verify_piqs_invariants.r
##
## Phase 1 verification helper for verify_v2.py.
##
## Loads fit.piqs.rda and the multi-year monthly cat file, computes the
## PIQS integral-preservation residual per (i, j, month) for both gpp and
## rtot, and writes a small JSON summary to verify_piqs_invariants.json.
##
## The invariant: for a piecewise integral quadratic f_m(t) = a (t-t_m)^2
## + b (t-t_m) + c on segment [t_m, t_{m+1}] of width delta, the integral
## is delta * ybar where ybar is the segment monthly mean (this is what
## "Integral" in PIQS means -- mass conservation per piece). So
##   integral_per_segment / delta_per_segment  ==  input monthly mean
## up to floating-point error.
##
## Output JSON contains: n_cells_checked, n_segments_checked, residual
## stats per quantity (max abs, max rel, median abs, fraction > 1e-9),
## plus the piqsfit.meta object verbatim.
##
## Usage: Rscript verify_piqs_invariants.r [output.json]

suppressPackageStartupMessages({
  library(ncdf4)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
out.json <- if (length(args) >= 1) args[1] else "verify_piqs_invariants.json"

## --- Load PIQS fit ---------------------------------------------------------
load("fit.piqs.rda")
stopifnot(exists("piqsfit.gpp"), exists("piqsfit.resp"), exists("piqsfit.time"))

dim.gpp  <- dim(piqsfit.gpp$a)
dim.resp <- dim(piqsfit.resp$a)
nmon     <- length(piqsfit.time)

## --- Load multi-year monthly cat file -------------------------------------
nc <- nc_open("monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc")
NPP <- ncvar_get(nc, "NPP")    # gC m-2 s-1, dim (lon, lat, time)
Rh  <- ncvar_get(nc, "Rh")
nc_close(nc)
stopifnot(dim(NPP)[3] == nmon)
gpp.in  <- -2 * NPP            # same convention as write_piqs.r
rtot.in <- Rh + NPP

## --- Compute segment widths in seconds ------------------------------------
## piqsfit.time stores left-edge knot times; compute delta to next knot.
## For the last segment, use median delta as a placeholder (would only
## bias the LAST month residual; we'll exclude it from stats).
delta <- diff(piqsfit.time)
delta <- c(delta, median(delta))   # length nmon
last  <- nmon                      # the last segment uses approx delta

## --- Per-segment integral computed from coefficients -----------------------
##
## For PIQS, integral over segment of width delta equals
##   integral = a/3 * delta^3 + b/2 * delta^2 + c * delta
## (since f(t) = a (t-t_m)^2 + b (t-t_m) + c, integrated from 0 to delta).
##
## Then mean-over-segment = integral / delta should equal the input monthly
## mean for that cell-month, by the PIQS construction.

compute.residual <- function(piqsfit, ybar.in) {
  a <- piqsfit$a; b <- piqsfit$b; c <- piqsfit$c
  d <- array(rep(delta, each = prod(dim(a)[1:2])), dim = dim(a))
  integral <- (a/3) * d^3 + (b/2) * d^2 + c * d
  mean.from.fit <- integral / d
  resid <- mean.from.fit - ybar.in
  abs.resid <- abs(resid)
  rel.resid <- abs.resid / pmax(abs(ybar.in), .Machine$double.eps)

  ## Drop last segment from stats (delta is approximate).
  abs.resid.s <- abs.resid[, , -last]
  rel.resid.s <- rel.resid[, , -last]
  ybar.s      <- ybar.in[, , -last]

  ## Restrict to "active" cells (any |ybar| > 0 anywhere over the record).
  active <- apply(abs(ybar.in), c(1, 2), max) > 0
  active.3d <- array(rep(active, dim(abs.resid.s)[3]), dim = dim(abs.resid.s))

  list(
    n_segments_active     = sum(active.3d, na.rm = TRUE),
    max_abs_residual      = max(abs.resid.s[active.3d], na.rm = TRUE),
    median_abs_residual   = median(abs.resid.s[active.3d], na.rm = TRUE),
    max_rel_residual      = max(rel.resid.s[active.3d], na.rm = TRUE),
    median_rel_residual   = median(rel.resid.s[active.3d], na.rm = TRUE),
    frac_abs_resid_gt_1e9 = mean(abs.resid.s[active.3d] > 1e-9, na.rm = TRUE),
    frac_rel_resid_gt_1e6 = mean(rel.resid.s[active.3d] > 1e-6, na.rm = TRUE)
  )
}

stats.gpp  <- compute.residual(piqsfit.gpp,  gpp.in)
stats.rtot <- compute.residual(piqsfit.resp, rtot.in)

## --- Assemble summary ------------------------------------------------------
out <- list(
  generated_at      = format(Sys.time(), tz = "UTC", usetz = TRUE),
  fit_dims_gpp      = dim.gpp,
  fit_dims_resp     = dim.resp,
  nmon              = nmon,
  fit_window        = c(min(piqsfit.time), max(piqsfit.time)),
  piqsfit_meta      = if (exists("piqsfit.meta")) piqsfit.meta else NA,
  gpp               = stats.gpp,
  rtot              = stats.rtot
)

write_json(out, out.json, pretty = TRUE, auto_unbox = TRUE, na = "string")
cat(sprintf("wrote %s\n", out.json))
cat(sprintf("  GPP : max_abs=%.3e, max_rel=%.3e, frac_abs_resid_gt_1e9=%.4f\n",
            stats.gpp$max_abs_residual, stats.gpp$max_rel_residual,
            stats.gpp$frac_abs_resid_gt_1e9))
cat(sprintf("  Rtot: max_abs=%.3e, max_rel=%.3e, frac_abs_resid_gt_1e9=%.4f\n",
            stats.rtot$max_abs_residual, stats.rtot$max_rel_residual,
            stats.rtot$frac_abs_resid_gt_1e9))
