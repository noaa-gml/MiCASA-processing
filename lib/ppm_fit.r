## lib/ppm_fit.r -- PPM (Colella & Woodward 1984) limited piecewise-PARABOLIC
## integral-preserving reconstruction. The literature's recommended
## "conservative + monotone (no overshoot) + smoother-than-linear" option
## (research synthesis 2026-06-18; Colella & Woodward, J. Comput. Phys. 54).
##
## Reconstructs, per month i, a parabola through that month's mean with
## edge values shared between adjacent months (so continuous where the
## monotonicity limiter does not bite), monotonicity-limited so it never
## overshoots the neighbouring monthly-mean envelope. Returns the same
## per-piece (a, b, c) quadratic layout as write_pchip.r / write_piqs.r, so
## diurnalize-ERA5.r consumes it transparently.
##
## Properties:
##   - Integral-preserving by construction (parabola mean == month mean).
##   - No overshoot: PPM limiting flattens at local extrema and steepens
##     edges so the parabola stays within the monthly-mean envelope.
##   - Piecewise QUADRATIC -> the (a,b,c) storage is exact, not a==0.
##   - C0-continuous at month edges EXCEPT where the limiter resets an edge
##     value (i.e. at turning points) -- strictly smoother than minmod-linear.

minmod3 <- function(a, b, cc) {
  s <- sign(a)
  ifelse(s == sign(b) & s == sign(cc), s * pmin(abs(a), abs(b), abs(cc)), 0)
}

## PPM per-cell reconstruction. x: knot times length n+1; ybar: means length n.
ppm.fit.cell <- function(x, ybar) {
  n <- length(x) - 1
  D <- diff(x)
  u <- ybar
  if (all(u == 0)) return(list(a = rep(0, n), b = rep(0, n), c = rep(0, n)))

  ## monotonized central differences dm_i (van Leer / PPM eq. 1.7-1.8)
  um1 <- c(u[1], u[1:(n - 1)])      # u_{i-1}
  up1 <- c(u[2:n], u[n])            # u_{i+1}
  dc  <- (up1 - um1) / 2
  dm  <- minmod3(dc, 2 * (u - um1), 2 * (up1 - u))

  ## continuous edge values a_{i+1/2} (PPM eq. 1.6); aedge[i] = edge i|i+1
  dmp1 <- c(dm[2:n], dm[n])
  aedge <- u + 0.5 * (up1 - u) - (1 / 6) * (dmp1 - dm)   # length n

  aL <- c(u[1], aedge[1:(n - 1)])   # left edge of cell i  = aedge_{i-1}
  aR <- aedge                       # right edge of cell i = aedge_i
  ## record-boundary cells: flatten to the mean (no outside neighbour)
  aL[1] <- u[1]; aR[1] <- u[1]; aL[n] <- u[n]; aR[n] <- u[n]

  ## PPM monotonicity limiter (eq. 1.10): flatten at extrema, else steepen
  ext <- (aR - u) * (u - aL) <= 0
  aL[ext] <- u[ext]; aR[ext] <- u[ext]

  d6 <- 6 * (u - 0.5 * (aL + aR))
  dd <- aR - aL
  ## overshoot left  -> reset aL ; overshoot right -> reset aR
  hi <- (!ext) & (dd * (u - 0.5 * (aL + aR)) >  dd^2 / 6)
  aL[hi] <- 3 * u[hi] - 2 * aR[hi]
  lo <- (!ext) & (-(dd^2) / 6 > dd * (u - 0.5 * (aL + aR)))
  aR[lo] <- 3 * u[lo] - 2 * aL[lo]

  ## parabola f(s)=aL + s*(aR-aL+a6) - a6*s^2, s=(t-x_i)/D, a6=6(u-(aL+aR)/2)
  a6 <- 6 * (u - 0.5 * (aL + aR))
  list(a = -a6 / D^2,
       b = (aR - aL + a6) / D,
       c = aL)
}

## Vectorized grid version of ppm.fit.cell: U is an [Ncell, nmon] matrix of
## monthly means (one cell per row); returns a/b/c each [Ncell, nmon].
## Equivalent to applying ppm.fit.cell row-by-row to FP tolerance (verified by
## tests/test_ppm_fit.r). ~6x faster than the per-cell loop in write_ppm.r.
ppm.fit.grid <- function(x, U) {
  N <- nrow(U); M <- ncol(U); D <- diff(x); Dm <- matrix(D, N, M, byrow = TRUE)
  uL <- cbind(U[, 1], U[, 1:(M - 1), drop = FALSE])
  uR <- cbind(U[, 2:M, drop = FALSE], U[, M])
  dL <- U - uL; dR <- uR - U
  dc <- (uR - uL) / 2
  dm <- minmod3(dc, 2 * dL, 2 * dR)
  dmp1 <- cbind(dm[, 2:M, drop = FALSE], dm[, M])
  aedge <- U + 0.5 * dR - (dmp1 - dm) / 6
  aL <- cbind(U[, 1], aedge[, 1:(M - 1), drop = FALSE]); aR <- aedge
  aL[, 1] <- U[, 1]; aR[, 1] <- U[, 1]; aL[, M] <- U[, M]; aR[, M] <- U[, M]
  ext <- (aR - U) * (U - aL) <= 0; aL[ext] <- U[ext]; aR[ext] <- U[ext]
  dd <- aR - aL; midd <- U - 0.5 * (aL + aR)
  hi <- (!ext) & (dd * midd >  dd^2 / 6); aL[hi] <- 3 * U[hi] - 2 * aR[hi]
  lo <- (!ext) & (-(dd^2) / 6 > dd * midd); aR[lo] <- 3 * U[lo] - 2 * aL[lo]
  a6 <- 6 * (U - 0.5 * (aL + aR))
  list(a = -a6 / Dm^2, b = (aR - aL + a6) / Dm, c = aL)
}
