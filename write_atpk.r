#!/usr/bin/env Rscript
##
## Area-to-point (ATP) kriging fitter driver. Produces the same per-piece
## (a,b,c) layout as write_pchip.r so diurnalize-ERA5.r consumes it via
## MICASA_FIT_RDA -- PLUS a per-piece kriging VARIANCE (piqsfit.gpp$var,
## piqsfit.resp$var): a principled prior-uncertainty on the sub-monthly flux,
## which the deterministic splines (PCHIP/PPM/PIQS) do not provide.
##
## Block data = monthly means; each month's Ns sub-points are kriged from a
## +-W-month window with an exponential covariance, then represented as a
## mass-preserving quadratic (exact coherence) with a selective sign-safe
## fallback for one-signed quantities. Kriging weights are data-independent, so
## they are precomputed once and applied as a fast filter to every cell.
##
## Output: $MICASA_FIT_OUT (default fit.piqs.rda). Knobs: MICASA_ATPK_W (window
## half-width, default 6), MICASA_ATPK_NS (sub-points/month, 6), MICASA_ATPK_RANGE
## (covariance range in months, 1.5). See lib/atpk_fit.r, docs/FITTER_COMPARISON.md
## (4.4) and PROPOSALS (18).
ct.setup()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
cfg <- micasa.config()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib", "atpk_fit.r"))

din  <- load.ncdf(micasa.out.monthly.cat(cfg))
gpp  <- -2 * din$NPP
rtot <- din$Rh + din$NPP
nmon <- length(din$time)
plt.start <- as.POSIXlt(din$time[1]); y0<-plt.start$year+1900; m0<-plt.start$mon+1
x.time <- as.numeric(seq(ISOdatetime(y0,m0,1,0,0,0,tz="UTC"), by="1 month", length.out=nmon+1))
h <- diff(x.time)

W     <- as.integer(Sys.getenv("MICASA_ATPK_W",     "6"))
Ns    <- as.integer(Sys.getenv("MICASA_ATPK_NS",    "6"))
range <- as.numeric(Sys.getenv("MICASA_ATPK_RANGE", as.character(ATPK_DEFAULT_RANGE)))
cat(sprintf("ATP kriging: W=%d Ns=%d range=%.2f mo, %d months\n", W, Ns, range, nmon))
ww <- atpk.window.weights(nmon, W = W, Ns = Ns, range = range)

mk <- function() list(a=array(0,c(360,180,nmon)), b=array(0,c(360,180,nmon)),
                      c=array(0,c(360,180,nmon)), var=array(0,c(360,180,nmon)))
piqsfit.gpp <- mk(); piqsfit.resp <- mk()
COEF_ZERO_THRESHOLD <- 1e-15

pb <- progress.bar.start(360*180, 360*180); ipb <- 0
for (i in 1:360) for (j in 1:180) {
  ipb <- ipb + 1
  if (max(abs(din$NPP[i,j,])) < COEF_ZERO_THRESHOLD &&
      max(abs(din$Rh [i,j,])) < COEF_ZERO_THRESHOLD) next
  fg <- atpk.apply.series(gpp [i,j,], h, ww)
  fr <- atpk.apply.series(rtot[i,j,], h, ww)
  piqsfit.gpp$a[i,j,]<-fg$a; piqsfit.gpp$b[i,j,]<-fg$b; piqsfit.gpp$c[i,j,]<-fg$c; piqsfit.gpp$var[i,j,]<-fg$var
  piqsfit.resp$a[i,j,]<-fr$a; piqsfit.resp$b[i,j,]<-fr$b; piqsfit.resp$c[i,j,]<-fr$c; piqsfit.resp$var[i,j,]<-fr$var
  pb <- progress.bar.print(pb, ipb)
}
progress.bar.end(pb)

piqsfit.time <- x.time[1:(length(x.time)-1)]
piqsfit.meta <- list(fitter="atpk", pad.left=0L, pad.right=0L, fit.range=range(x.time),
  saved.range=range(piqsfit.time), written.at=format(Sys.time(),tz="UTC",usetz=TRUE),
  atpk=list(W=W, Ns=Ns, range=range),
  notes="Area-to-point kriging (Kyriakidis 2004); windowed exponential covariance; exact coherence; sign-safe; per-piece kriging variance in $var")
out <- Sys.getenv("MICASA_FIT_OUT","fit.piqs.rda")
save(file=out, piqsfit.gpp, piqsfit.resp, piqsfit.time, piqsfit.meta)
cat(sprintf("Wrote %s (ATP kriging; with $var uncertainty arrays)\n", out))
