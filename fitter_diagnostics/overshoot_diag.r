#!/usr/bin/env Rscript
## Quantify MAGNITUDE overshoot of the production PCHIP fit:
## how far the instantaneous sub-monthly flux peaks above the local
## monthly-mean envelope, per land cell-piece. Compares against the
## two mass-conserving "linear" alternatives Andy is asking about.
##
## All quantities derived from the shipped fit.piqs.rda coefficients.
## Monthly mean recovered exactly as the piece integral (mass-preserving):
##   u_i = a*D^2/3 + b*D/2 + c       (D = piece width in seconds)
suppressWarnings(load("fit.piqs.rda"))

## Reconstruct the 304 monthly knot edges the fit was built on.
t0  <- as.POSIXct("2001-01-01", tz = "UTC")
edges <- as.numeric(seq(t0, by = "1 month", length.out = length(piqsfit.time) + 1))
D <- diff(edges)                      # 303 piece widths (sec); 28-31 days

analyse <- function(fit, label) {
  N <- 360 * 180; M <- length(D)
  a <- matrix(fit$a, N, M); b <- matrix(fit$b, N, M); c0 <- matrix(fit$c, N, M)
  Dm <- matrix(D, N, M, byrow = TRUE)

  ## flux at the two knots and the monthly mean (mass-preserving integral)
  fL <- c0                                   # f at left knot
  fR <- a * Dm^2 + b * Dm + c0               # f at right knot
  u  <- a * Dm^2 / 3 + b * Dm / 2 + c0       # monthly mean

  ## within-piece vertex (quadratic extremum), only if interior to [0,D]
  sv  <- ifelse(a != 0, -b / (2 * a), -1)    # vertex offset from left knot (sec)
  intr <- a != 0 & sv > 0 & sv < Dm
  fV  <- ifelse(intr, c0 - b^2 / (4 * a), 0) # vertex flux

  peak <- pmax(abs(fL), abs(fR), ifelse(intr, abs(fV), 0))  # peak |flux| in piece
  knot <- pmax(abs(fL), abs(fR))                            # peak |flux| at knots

  ## local monthly-mean envelope: max |mean| over {i-1, i, i+1}
  um <- abs(u)
  uprev <- cbind(um[, 1], um[, 1:(M - 1)])
  unext <- cbind(um[, 2:M], um[, M])
  env <- pmax(uprev, um, unext)

  ## land cell-pieces only: cell active somewhere AND a real local signal
  active <- rowSums(um) > 1e-15
  keep   <- active & env > 1e-12
  r.peak <- (peak / env)[keep]               # within-piece overshoot ratio
  r.knot <- (knot / env)[keep]               # knot overshoot ratio

  qs <- c(.5, .9, .99, .999, 1)
  cat(sprintf("\n==== %s : %d land cell-pieces ====\n", label, length(r.peak)))
  cat("within-piece peak/envelope ratio:\n")
  cat(sprintf("  pctl 50/90/99/99.9/max = %.2f / %.2f / %.2f / %.2f / %.2f\n",
              quantile(r.peak, qs[1]), quantile(r.peak, qs[2]),
              quantile(r.peak, qs[3]), quantile(r.peak, qs[4]), max(r.peak)))
  for (thr in c(1.0, 1.25, 1.5, 2.0, 3.0))
    cat(sprintf("  %% pieces peak > %4.2f x envelope : %6.3f%%\n",
                thr, 100 * mean(r.peak > thr)))
  cat(sprintf("knot peak/envelope: pctl 99 = %.2f, max = %.2f (FC bound ~3)\n",
              quantile(r.knot, .99), max(r.knot)))

  ## ---- mass-conserving LINEAR alternatives, same cells ----
  ## (A) linear-on-cumulative = piecewise-CONSTANT flux: f==u, ratio==1 always.
  ## (B) linear-flux trapezoidal (proposal #9 "linear PIQS"):
  ##     knot recursion y_{i+1} = 2 u_i - y_i, seeded y_1 = u_1.
  yk <- matrix(0, N, M + 1); yk[, 1] <- u[, 1]
  for (i in 1:M) yk[, i + 1] <- 2 * u[, i] - yk[, i]
  ## within-piece peak of a straight line = max(|y_i|,|y_{i+1}|)
  pk.lin <- pmax(abs(yk[, 1:M]), abs(yk[, 2:(M + 1)]))
  ## sign flip: a GPP/Rh segment whose endpoints straddle 0 (knot wrong sign)
  signflip <- (yk[, 1:M] * u) < 0 | (yk[, 2:(M + 1)] * u) < 0
  r.lin <- (pk.lin / env)[keep]
  cat("linear-flux trapezoidal (proposal #9): \n")
  cat(sprintf("  peak/envelope pctl 99/max = %.2f / %.2f ; %% knot SIGN FLIPS = %.3f%%\n",
              quantile(r.lin, .99), max(r.lin), 100 * mean(signflip[keep])))
  cat("piecewise-constant (linear-on-cumulative): peak/envelope = 1.00 by construction; qmod=0 (no within-month gradient)\n")

  ## worst within-piece overshoot location
  idx <- which(keep)[which.max(r.peak)]
  cell <- ((idx - 1) %% N) + 1; mon <- ((idx - 1) %/% N) + 1
  ii <- ((cell - 1) %% 360) + 1; jj <- ((cell - 1) %/% 360) + 1
  lat <- -90 + (jj - 0.5); lon <- -180 + (ii - 0.5)
  cat(sprintf("worst piece: lon %.1f lat %.1f, %s, ratio=%.2f\n",
              lon, lat, format(as.POSIXct(edges[mon], origin = "1970-01-01", tz = "UTC"), "%Y-%m"),
              max(r.peak)))
}

analyse(piqsfit.gpp,  "GPP (-2*NPP)")
analyse(piqsfit.resp, "RESP (Rh+NPP)")
