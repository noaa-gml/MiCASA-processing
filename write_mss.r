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
##   - tests/bakeoff_mss.py (Python prototype with quadprog)
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

## QP fitter core (mss.fit.setup, mss.fit.cell) -- lib/mss_fit.r,
## unit-tested standalone by tests/test_mss_fit.r. The QP smoothness
## Hessian and constraint matrices depend only on the knot positions,
## so build them once here and reuse across every cell.
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "mss_fit.r"))
mss.setup <- mss.fit.setup(x.time)

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
    fit.gpp  <- mss.fit.cell(x.time, gpp [i, j, ], mss.setup)
    fit.rtot <- mss.fit.cell(x.time, rtot[i, j, ], mss.setup)
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
