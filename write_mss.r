#!/usr/bin/env Rscript
##
## Drop-in alternative to write_piqs.r using a monotone smoothing
## spline (MSS): a cubic spline F on cumulative-integral knot values
## that minimizes int(F'')^2 over the record subject to F(t_k) = F_k
## (integral preservation) AND F'(t) >= 0 (cumulative monotone, i.e.
## flux f = F' >= 0 everywhere). Solved per cell as a QP on knot
## slopes m_0, ..., m_n via the quadprog package.
##
## Storage layout matches write_piqs.r (a, b, c per piece per cell)
## so diurnalize-ERA5.r consumes fit.piqs.rda transparently.
##
## Properties (vs PIQS):
##   - Same smoothness optimum WHEN the monotonicity constraints are
##     not binding (smooth-seasonality cells: tropical evergreens,
##     mid-latitudes outside winter). Recovers PIQS coefficients
##     exactly in those cases.
##   - On rough cells where PIQS overshoots zero (boreal/tundra
##     winters, semi-arid transitions): drops sign-flip rate ~5-25x
##     while preserving PIQS-level smoothness. Bake-off shows
##     residual flips ~0.6-1.2% from constraint discretization
##     (constraints enforced at 8 test points per segment); 0%
##     would require the full QCQP variant or a denser grid.
##
## Properties (vs PCHIP):
##   - Smoother in absolute terms (matches PIQS); PCHIP is ~50%
##     rougher in max|df| but absolute differences are tiny (<2e-11)
##     and invisible at hourly sampling.
##   - More expensive: ~30-450 ms per cell (QP solve) vs ~1 ms
##     for PCHIP. ~30 min single-threaded for full grid.
##   - NOT provably zero sign flips (PCHIP is); test-point
##     discretization is leaky between sample points.
##
## See:
##   - bakeoff_mss.py (Python prototype with quadprog)
##   - README.ash methodological note (10) for the broader context
##   - Pya & Wood 2015, "Shape constrained additive models" / He & Shi
##     1998, "Monotone B-spline smoothing", for related literature.

suppressPackageStartupMessages(library(quadprog))

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

## Pre-build the constant smoothness Hessian and constraint template.
## The Hessian and the "shape" of constraints depend only on x.time
## (the knot positions), not on ybar values. Precompute once outside
## the cell loop.
n.test.per.segment <- 8
h.seg <- diff(x.time)            # length n
n.seg <- length(h.seg)
n.var <- n.seg + 1               # number of m_k (knot slopes)

build.G <- function(h) {
  G <- matrix(0, n.var, n.var)
  for (k in seq_len(n.seg)) {
    inv.h <- 1 / h[k]
    G[k,     k]     <- G[k,     k]     + 8 * inv.h
    G[k + 1, k + 1] <- G[k + 1, k + 1] + 8 * inv.h
    G[k,     k + 1] <- G[k,     k + 1] + 4 * inv.h
    G[k + 1, k]     <- G[k + 1, k]     + 4 * inv.h
  }
  ## tiny ridge for numerical stability of solve.QP (D must be PD)
  G + 1e-20 * diag(n.var)
}
G.const <- build.G(h.seg)

## Constraint matrix C is partly constant (the test-point coefficients
## on m_k, m_{k+1}) and partly cell-dependent (the b vector includes u_k).
## Pre-build the dense (n.var, n.constraints) C matrix.
n.constraints <- n.seg * n.test.per.segment + 2
C.const <- matrix(0, n.var, n.constraints)
test.s  <- (seq_len(n.test.per.segment) - 0.5) / n.test.per.segment
j <- 0
for (k in seq_len(n.seg)) {
  for (s in test.s) {
    j <- j + 1
    C.const[k,     j] <- 3 * s * s - 4 * s + 1
    C.const[k + 1, j] <- 3 * s * s - 2 * s
  }
}
j <- j + 1; C.const[1,     j] <- 1
j <- j + 1; C.const[n.var, j] <- 1
## Per-segment, the s-dependent constant in the constraint b is
##   b_j = -(6s - 6s^2) * u_k
## We'll build b per cell using a (n.constraints, n.seg) coefficient
## matrix that maps u to b (boundary rows have zero coefficient).
b.coef <- matrix(0, n.constraints, n.seg)
j <- 0
for (k in seq_len(n.seg)) {
  for (s in test.s) {
    j <- j + 1
    b.coef[j, k] <- -(6 * s - 6 * s * s)
  }
}

mss.fit.cell <- function(x, ybar) {
  if (all(ybar == 0)) {
    return(list(a = rep(0, n.seg), b = rep(0, n.seg), c = rep(0, n.seg)))
  }
  ## Like write_pchip.r, allow the data to be uniformly negative by
  ## flipping sign and unflipping the final coefficients.
  if (mean(ybar) < 0) {
    u <- -ybar
    sign.flip <- -1
  } else {
    u <- ybar
    sign.flip <- 1
  }
  inv.h <- 1 / h.seg
  ## Linear term in 0.5 m^T G m - a^T m form: a[k] = +12 u[k] / h[k]
  a.lin <- numeric(n.var)
  for (k in seq_len(n.seg)) {
    a.lin[k]     <- a.lin[k]     + 12 * u[k] * inv.h[k]
    a.lin[k + 1] <- a.lin[k + 1] + 12 * u[k] * inv.h[k]
  }
  ## b vector for inequality constraints: b[j] = sum_k b.coef[j,k] * u[k]
  bvec <- as.vector(b.coef %*% u)
  sol <- tryCatch(
    solve.QP(Dmat = G.const, dvec = a.lin, Amat = C.const, bvec = bvec, meq = 0),
    error = function(e) NULL)
  if (is.null(sol)) {
    ## Fall back to PCHIP on this cell if QP fails (rare; very pathological
    ## inputs only).
    fn <- splinefun(x, c(0, cumsum(u * h.seg)), method = "monoH.FC")
    m <- fn(x, deriv = 1)
  } else {
    m <- sol$solution
  }
  Q <- -6 * u + 3 * m[1:n.seg] + 3 * m[2:n.var]
  L <-  6 * u - 4 * m[1:n.seg] - 2 * m[2:n.var]
  K <- m[1:n.seg]
  list(a = sign.flip * Q / h.seg^2,
       b = sign.flip * L / h.seg,
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
    fit.gpp  <- mss.fit.cell(x.time, gpp [i, j, ])
    fit.rtot <- mss.fit.cell(x.time, rtot[i, j, ])
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
piqsfit.meta <- list(fitter      = "mss",
                     pad.left    = 0L,
                     pad.right   = 0L,
                     fit.range   = range(x.time),
                     saved.range = range(piqsfit.time),
                     written.at  = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     notes       = sprintf("Monotone smoothing spline via QP on knot slopes; %d test points/segment; falls back to PCHIP per-cell if QP fails",
                                           n.test.per.segment))

save(file = "fit.piqs.rda",
     piqsfit.gpp, piqsfit.resp, piqsfit.time, piqsfit.meta)
cat("Wrote fit.piqs.rda (MSS via quadprog)\n")
