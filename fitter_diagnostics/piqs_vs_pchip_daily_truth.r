## PIQS vs PCHIP, validated against MiCASA's OWN daily product.
##
## The fitter's job: reconstruct the sub-monthly NPP/Rh shape from MONTHLY means.
## MiCASA already ships that sub-monthly shape daily (daily_1x1) -- the ground truth
## the monthly->fit step discards. So: evaluate each PRODUCTION fit (PIQS fit.piqs_v1,
## PCHIP fit.pchip) at daily resolution and compare to the actual daily product, in
## the fitter's own (gpp,rtot,NEE) space.  NEE=(Rh-NPP)/12, gpp=-2NPP/12, rtot=(Rh+NPP)/12.
## No scale/quantity/weather/partitioning confounds; every 1-deg land cell, full year.
suppressMessages(library(ncdf4))
P <- new.env(); load("fit.piqs_v1.rda", envir=P)
C <- new.env(); load("fit.pchip.rda",  envir=C)
ptime <- get("piqsfit.time", P)
Pg<-get("piqsfit.gpp",P); Pr<-get("piqsfit.resp",P)
Cg<-get("piqsfit.gpp",C); Cr<-get("piqsfit.resp",C)
nmon <- length(ptime)
cat("ptime: n=",nmon," median dt=",median(diff(ptime))," (≈2.6e6 -> seconds)\n",sep="")

evalc <- function(co,k,dt) co$a[,,k]*dt^2 + co$b[,,k]*dt + co$c[,,k]
dim_of <- function(y,m){ nm<-if(m==12) as.Date(sprintf("%d-01-01",y+1)) else as.Date(sprintf("%d-%02d-01",y,m+1))
                         as.integer(nm - as.Date(sprintf("%d-%02d-01",y,m))) }

YEAR <- 2020
k0 <- (YEAR-2001)*12
files <- Sys.glob(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_%d*.nc",YEAR))
sseP<-sseC<-nobs<-array(0,c(360,180))
gppwrongP<-gppwrongC<-ncell<-0
for(f in files){
  nc<-nc_open(f); npp<-ncvar_get(nc,"NPP"); rh<-ncvar_get(nc,"Rh"); nc_close(nc)
  ymd<-regmatches(f,regexpr("[0-9]{8}",f)); m<-as.integer(substr(ymd,5,6)); d<-as.integer(substr(ymd,7,8))
  k<-k0+m; if(k<1||k>=nmon) next
  h<-ptime[k+1]-ptime[k]; dt<-((d-0.5)/dim_of(YEAR,m))*h
  gtru<- -2*npp/12; rtru<-(rh+npp)/12; netru<-(rh-npp)/12
  gP<-evalc(Pg,k,dt); neP<-gP+evalc(Pr,k,dt)
  gC<-evalc(Cg,k,dt); neC<-gC+evalc(Cr,k,dt)
  land <- is.finite(netru) & is.finite(neP) & is.finite(neC) & (abs(npp)+abs(rh) > 0)
  e2P<-(neP-netru)^2; e2C<-(neC-netru)^2
  sseP[land]<-sseP[land]+e2P[land]; sseC[land]<-sseC[land]+e2C[land]; nobs[land]<-nobs[land]+1
  ## overshoot: reconstructed gpp should be <= 0 (uptake). count wrong-sign on land.
  gppwrongP<-gppwrongP+sum(gP[land]>0); gppwrongC<-gppwrongC+sum(gC[land]>0); ncell<-ncell+sum(land)
}
ok <- nobs>0
rmseP<-sqrt(sseP/pmax(nobs,1)); rmseC<-sqrt(sseC/pmax(nobs,1))
rP<-rmseP[ok]; rC<-rmseC[ok]
cat(sprintf("\n=== PIQS vs PCHIP vs MiCASA daily truth (%d, %d land cells, %d cell-days) ===\n",YEAR,sum(ok),sum(nobs)))
cat(sprintf("daily NEE reconstruction RMSE (fitter vs true daily), median over land cells:\n"))
cat(sprintf("  PIQS  : %.4e\n  PCHIP : %.4e\n  PCHIP/PIQS ratio : %.3f  (<1 = PCHIP closer to the true daily shape)\n",
            median(rP),median(rC),median(rC)/median(rP)))
cat(sprintf("  cells where PCHIP RMSE < PIQS RMSE : %d/%d (%.0f%%)\n",
            sum(rC<rP),length(rC),100*mean(rC<rP)))
cat(sprintf("  mean RMSE (area-unweighted): PIQS %.4e  PCHIP %.4e\n", mean(rP),mean(rC)))
cat(sprintf("\nwrong-sign GPP in the reconstruction (should be 0): PIQS %.3f%% of cell-days ; PCHIP %.3f%%\n",
            100*gppwrongP/ncell, 100*gppwrongC/ncell))
saveRDS(list(rmseP=rmseP,rmseC=rmseC,nobs=nobs), "fitter_diagnostics/piqs_vs_pchip_daily_truth.rds")
cat("DONE\n")
