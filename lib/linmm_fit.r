## lib/linmm_fit.r -- minmod-limited integral-preserving piecewise-LINEAR
## fitter core (finite-volume / MUSCL reconstruction).
##
## linmm.fit.cell() is the per-grid-cell smoother that write_linmm.r runs
## over all 360x180 cells. Given a monthly-mean flux series it reconstructs,
## for each month, a straight line through that month's mean with a
## minmod-limited slope, and returns it in the same per-piece (a, b, c)
## layout as write_piqs.r / write_pchip.r (here a == 0, the flux is linear).
##
## Properties (the reason this exists -- see docs/PROPOSALS.md (17)):
##   - Integral-preserving by construction: the line passes through the
##     monthly mean at the month centre, so its integral over the month is
##     exactly mean * width, for ANY slope.
##   - No overshoot: the minmod limiter caps the slope at the smaller of the
##     two neighbouring secants (and zeros it at a sign change / turning
##     point), so a month's edge values never exceed the envelope of its
##     neighbouring monthly means. Provably non-oscillatory.
##   - Carries a within-month gradient (unlike piecewise-constant), but is
##     DISCONTINUOUS at month boundaries (the finite-volume cost).
##
## Contrast with a CONTINUOUS integral-preserving line (y_{i+1}=2*m_i-y_i):
## that recursion has a pole at the Nyquist frequency and is numerically
## unstable on real seasonal data -- see PROPOSALS.md (9). minmod trades the
## continuity for stability + the no-overshoot guarantee.
##
## Extracted for standalone unit testing (tests/test_linmm_fit.r); pure base
## R, no ct / ncdf4 dependency.

minmod <- function(a, b) 0.5 * (sign(a) + sign(b)) * pmin(abs(a), abs(b))

## minmod-limited linear reconstruction: per-cell, return list(a,b,c) of
## length nmon such that piece i is f(t) = a*(t-x[i])^2 + b*(t-x[i]) + c.
## a is identically 0 (the flux is linear within each month).
linmm.fit.cell <- function(x, ybar) {
  n <- length(x) - 1
  D  <- diff(x)                          # piece widths
  tc <- (x[1:n] + x[2:(n + 1)]) / 2      # month centres
  if (all(ybar == 0))
    return(list(a = rep(0, n), b = rep(0, n), c = rep(0, n)))

  slope <- numeric(n)                    # flux per second within each month
  if (n >= 3) {
    sec <- diff(ybar) / diff(tc)         # length n-1: secant between i and i+1
    ## interior months: minmod of the two adjoining secants
    slope[2:(n - 1)] <- minmod(sec[1:(n - 2)], sec[2:(n - 1)])
  }
  ## record-boundary months (i=1, i=n) keep slope 0: one-sided extrapolation
  ## there could overshoot on the outer edge, and a flat first/last month is
  ## cosmetically negligible.

  ## f(t) = ybar_i + slope_i*(t - tc_i) = slope_i*(t - x_i) + [ybar_i - slope_i*D_i/2]
  list(a = rep(0, n),
       b = slope,
       c = ybar - slope * D / 2)
}

## Vectorized grid version of linmm.fit.cell: U is an [Ncell, nmon] matrix of
## monthly means (one cell per row); returns a/b/c each [Ncell, nmon].
## Bit-for-bit equivalent to applying linmm.fit.cell row-by-row (verified by
## tests/test_linmm_fit.r). ~6x faster than the per-cell loop in write_linmm.r.
linmm.fit.grid <- function(x, U) {
  N <- nrow(U); M <- ncol(U)
  D  <- diff(x)
  tc <- (x[1:M] + x[2:(M + 1)]) / 2
  slope <- matrix(0, N, M)
  if (M >= 3) {
    dtc <- diff(tc)                                   # length M-1
    sec <- (U[, 2:M, drop = FALSE] - U[, 1:(M - 1), drop = FALSE]) /
           matrix(dtc, N, M - 1, byrow = TRUE)        # [N, M-1]
    slope[, 2:(M - 1)] <- minmod(sec[, 1:(M - 2), drop = FALSE],
                                 sec[, 2:(M - 1), drop = FALSE])
  }
  Dm <- matrix(D, N, M, byrow = TRUE)
  list(a = matrix(0, N, M), b = slope, c = U - slope * Dm / 2)
}
