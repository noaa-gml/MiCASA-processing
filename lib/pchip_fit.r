## lib/pchip_fit.r -- PCHIP-on-cumulative fitter core.
##
## pchip.fit.cell() is the per-grid-cell smoother that write_pchip.r runs
## over all 360x180 cells. Given a monthly-mean flux series it fits a
## Fritsch-Carlson monotone cubic Hermite spline to the cumulative
## integral and returns the piecewise-quadratic derivative as per-piece
## (a, b, c) coefficients -- the smoothed sub-monthly flux.
##
## Extracted from write_pchip.r so it can be unit-tested standalone
## (tests/test_pchip_fit.r) without the netCDF I/O. Pure base R
## (stats::splinefun); no ct / ncdf4 dependency.
##
## Fritsch & Carlson 1980, "Monotone Piecewise Cubic Interpolation",
## SIAM J. Numer. Anal. 17(2) pp 238-246.

## PCHIP-on-cumulative: per-cell, return list(a,b,c) of length nmon
## such that piece i is f(t) = a*(t - x[i])^2 + b*(t - x[i]) + c.
pchip.fit.cell <- function(x, ybar) {
  n <- length(x) - 1
  delta <- diff(x)
  ## Fritsch-Carlson is monotone iff input is monotone; cumulative is
  ## monotone iff all ybar same sign. We pass the SIGN of the data
  ## explicitly: PCHIP on cum_F where F is non-decreasing for ybar>=0.
  ## For GPP (ybar <= 0), the cumulative is non-increasing -- we negate
  ## before fitting so monoH.FC sees a non-decreasing sequence, then
  ## negate the resulting derivative coefficients.
  if (all(ybar == 0)) {
    return(list(a = rep(0, n), b = rep(0, n), c = rep(0, n)))
  }
  ## Determine direction. If data is mixed-sign (not the GPP/Rh case
  ## here, both are unsigned-magnitude), splinefun(monoH.FC) won't
  ## guarantee monotone -- but it will still return SOMETHING. We
  ## handle the unsigned case (flip sign if needed).
  if (mean(ybar) < 0) {
    ybar.s <- -ybar
    sign.flip <- -1
  } else {
    ybar.s <- ybar
    sign.flip <- 1
  }
  cum.F <- c(0, cumsum(ybar.s * delta))   # length n+1
  ## R's splinefun with monoH.FC returns a function of t. Internal
  ## representation can be queried via the function's environment for
  ## the per-piece slopes.
  fn <- splinefun(x, cum.F, method = "monoH.FC")
  ## Sample F' at the knots to get slopes m_k. monoH.FC uses
  ## Fritsch-Carlson per-knot slopes; calling fn(x, deriv = 1) returns
  ## the slope EXACTLY at each knot.
  m <- fn(x, deriv = 1)
  ## Build the per-piece quadratic. Cubic Hermite on segment k:
  ##   f(s) = (6s - 6s^2) u_k + (3s^2 - 4s + 1) m_k + (3s^2 - 2s) m_{k+1}
  ## where s = (t - x[k])/h_k, u_k = (F[k+1] - F[k])/h_k = ybar.s[k],
  ## h_k = delta[k]. Convert to t coordinates:
  ##   f(t) = (Q/h^2) (t - x[k])^2 + (L/h) (t - x[k]) + K
  ## with Q = -6 u + 3 m_k + 3 m_{k+1},
  ##      L =  6 u - 4 m_k - 2 m_{k+1},
  ##      K = m_k.
  u <- ybar.s
  Q <- -6 * u + 3 * m[1:n] + 3 * m[2:(n + 1)]
  L <-  6 * u - 4 * m[1:n] - 2 * m[2:(n + 1)]
  K <- m[1:n]
  list(a = sign.flip * Q / delta^2,
       b = sign.flip * L / delta,
       c = sign.flip * K)
}
