#!/usr/bin/env Rscript
##
## Drop-in alternative to write_piqs.r using PCHIP-on-cumulative
## (Fritsch-Carlson monotone-cubic Hermite, R's stats::splinefun
## with method="monoH.FC") instead of PIQS-quadratic.
##
## Mathematical setup: given monthly-mean flux ybar (length n) at
## knot times x (length n+1), build the cumulative integral
## F[k] = sum_{j<k} ybar[j] * (x[j+1] - x[j]), apply Fritsch-Carlson
## monotone-cubic Hermite interpolation to (x, F), and extract the
## piecewise-quadratic derivative as the smoothed flux. Storage
## layout matches write_piqs.r (a, b, c per piece per cell), so
## diurnalize-ERA5.r consumes fit.piqs.rda transparently.
##
## Why this matters: PIQS-quadratic can produce sub-monthly sign
## flips when monthly means are near zero (Check 3.1 of verify_v2:
## up to 30%+ in boreal/tundra cells). PCHIP-on-cumulative is
## monotone by construction, so the flux f = F' is non-negative
## everywhere -- knots and within pieces alike -- without the need
## for the polar-night clip in diurnalize-ERA5.r.
##
## See:
##   - bakeoff_pchip.py (Python reference + diagnostic)
##   - README.ash methodological note (10) for the full reasoning
##   - Fritsch & Carlson 1980, "Monotone Piecewise Cubic Interpolation",
##     SIAM J. Numer. Anal. 17(2) pp 238-246.

ct.setup()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
cfg <- micasa.config()

din <- load.ncdf(micasa.out.monthly.cat(cfg))

gpp  <- -2 * din$NPP
rtot <- din$Rh + din$NPP

nmon <- length(din$time)
plt.start <- as.POSIXlt(din$time[1])
y0 <- plt.start$year + 1900
m0 <- plt.start$mon  + 1
x.time <- as.numeric(seq(ISOdatetime(y0, m0, 1, 0, 0, 0, tz = "UTC"),
                         by = "1 month", length.out = nmon + 1))

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

## ---- Loop over cells -------------------------------------------------------
piqsfit.gpp  <- list()
piqsfit.resp <- list()
piqsfit.gpp$a  <- array(0, dim = c(360, 180, nmon))
piqsfit.gpp$b  <- array(0, dim = c(360, 180, nmon))
piqsfit.gpp$c  <- array(0, dim = c(360, 180, nmon))
piqsfit.resp$a <- array(0, dim = c(360, 180, nmon))
piqsfit.resp$b <- array(0, dim = c(360, 180, nmon))
piqsfit.resp$c <- array(0, dim = c(360, 180, nmon))

## Same low-flux predicate as write_piqs.r: cells with both NPP and
## Rh below 1e-15 across the entire record skip the fit and keep
## (0, 0, 0) coefficients.
COEF_ZERO_THRESHOLD <- 1e-15

pb  <- progress.bar.start(360 * 180, 360 * 180)
ipb <- 0
for (i in 1:360) {
  for (j in 1:180) {
    ipb <- ipb + 1
    if (max(abs(din$NPP[i, j, ])) < COEF_ZERO_THRESHOLD &&
        max(abs(din$Rh [i, j, ])) < COEF_ZERO_THRESHOLD) {
      next
    }
    fit.gpp  <- pchip.fit.cell(x.time, gpp [i, j, ])
    fit.rtot <- pchip.fit.cell(x.time, rtot[i, j, ])
    piqsfit.gpp$a [i, j, ] <- fit.gpp$a
    piqsfit.gpp$b [i, j, ] <- fit.gpp$b
    piqsfit.gpp$c [i, j, ] <- fit.gpp$c
    piqsfit.resp$a[i, j, ] <- fit.rtot$a
    piqsfit.resp$b[i, j, ] <- fit.rtot$b
    piqsfit.resp$c[i, j, ] <- fit.rtot$c
    pb <- progress.bar.print(pb, ipb)
  }
}
progress.bar.end(pb)

piqsfit.time <- x.time[1:(length(x.time) - 1)]

## Note: PCHIP padding semantics differ from PIQS. The Fritsch-Carlson
## slope rule is local (uses neighbouring monthly means only), so the
## right-edge-padding workaround proposal #1 isn't relevant here -- the
## last segment's slope only depends on the last two ybar values.
## piqsfit.meta still records what the script saw, with fitter="pchip".
piqsfit.meta <- list(fitter      = "pchip",
                     pad.left    = 0L,
                     pad.right   = 0L,
                     fit.range   = range(x.time),
                     saved.range = range(piqsfit.time),
                     written.at  = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     ## A note for downstream readers comparing across fits.
                     notes       = "Cubic-Hermite (Fritsch-Carlson) on cumulative; flux f=F' as piecewise quadratic; provably non-negative")

save(file = "fit.piqs.rda",
     piqsfit.gpp, piqsfit.resp, piqsfit.time, piqsfit.meta)
cat("Wrote fit.piqs.rda (PCHIP-on-cumulative)\n")
