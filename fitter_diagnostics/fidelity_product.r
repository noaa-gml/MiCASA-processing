## C1 + M1 + M3 fix. Scores the REAL diurnalized NEE product (meteo included) vs
## MiCASA daily NEE -- not the bare fitter quadratic -- and:
##  (C1) reports the fraction of daily-NEE variance the fitter's qmod controls
##       (var(daily qmod)/var(daily NEE)); if tiny, no daily metric can separate
##       fitters because the fitter governs little of the scored variance.
##  (M1) compares PCHIP vs PPM with a by-CELL block bootstrap + 95% CI (not pooled
##       autocorrelated cell-months).
##  (M3) normalizes by a per-cell RMS flux scale (no divide-by-near-zero envelope).
suppressWarnings(suppressMessages(library(ncdf4)))
N<-360*180; dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31); ndays<-sum(dpm)
daymean<-function(A){n<-dim(A)[3]/24; d<-matrix(0,N,n); for(k in 1:n) d[,k]<-rowMeans(matrix(A[,,((k-1)*24+1):(k*24)],N,24)); d}
rdvar<-function(dir,m,v){f<-nc_open(sprintf("%s/fluxes_2020%02d.nc",dir,m)); x<-ncvar_get(f,v); nc_close(f); x}
## accumulate daily-mean series [N, ndays]
NEEp<-NEEpp<-NEEq<-QMOD<-TRUTH<-matrix(0,N,ndays); col<-0
for(m in 1:12){ d<-dpm[m]; idx<-(col+1):(col+d); col<-col+d
  NEEp[,idx]<-daymean(rdvar("ERA5_2020_pchip",m,"NEE"))
  NEEpp[,idx]<-daymean(rdvar("ERA5_2020_ppm",m,"NEE"))
  NEEq[,idx]<-daymean(rdvar("ERA5_2020_piqs",m,"NEE"))
  QMOD[,idx]<-daymean(rdvar("ERA5_2020_pchip",m,"QGPP")) + daymean(rdvar("ERA5_2020_pchip",m,"qresp"))
  ## truth biosphere daily NEE = (Rh - NPP)/12 (same as product nee=gpp+resp)
  for(dd in 1:d){ f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%02d%02d.nc",m,dd))
    TRUTH[,idx[dd]]<-(as.vector(ncvar_get(f,"Rh"))-as.vector(ncvar_get(f,"NPP")))/12; nc_close(f) } }
scl<-sqrt(rowMeans(TRUTH^2)); land<-scl>1e-9                       # per-cell RMS scale (robust)
j<-((1:N-1)%/%360)+1; latc<- -90+j-0.5; band<-ifelse(abs(latc)<=23.5,"tropics",ifelse(abs(latc)<=50,"temperate","boreal/polar"))
rmse<-function(P) sqrt(rowMeans((P-TRUTH)^2))/scl
rp<-rmse(NEEp); rpp<-rmse(NEEpp); rq<-rmse(NEEq)
## (C1) fitter variance fraction
vfrac<-apply(QMOD,1,var)/apply(TRUTH,1,var)
cat("== C1: product-level daily-NEE fidelity (RMSE / per-cell RMS scale) ==\n")
cat(sprintf("  land cells %d ; PCHIP mean %.3f  PPM %.3f  PIQS %.3f\n", sum(land), mean(rp[land]), mean(rpp[land]), mean(rq[land])))
cat(sprintf("  FITTER variance fraction var(qmod_daily)/var(NEE_daily): median %.3f  90th %.3f  (=> fitter controls ~%.0f%% of daily NEE variance)\n",
    median(vfrac[land],na.rm=T), quantile(vfrac[land],.9,na.rm=T), 100*median(vfrac[land],na.rm=T)))
## (M1) PCHIP vs PPM, by-CELL block bootstrap (resample whole cells)
d<-(rp-rpp)[land]; lc<-which(land); set.seed(1); B<-2000; bs<-numeric(B)
for(b in 1:B){ s<-sample(length(d),replace=TRUE); bs[b]<-mean(d[s]) }
ci<-quantile(bs,c(.025,.975))
cat(sprintf("\n== M1: PCHIP - PPM product fidelity, by-cell bootstrap (N=%d cells, B=2000) ==\n", sum(land)))
cat(sprintf("  mean Delta(PCHIP-PPM) = %.5f  95%% CI [%.5f, %.5f]  (>0 => PCHIP worse)\n", mean(d), ci[1], ci[2]))
cat(sprintf("  => %s\n", if(ci[1]<0 & ci[2]>0) "CI straddles 0: indistinguishable (a real tie, with a CI)" else "CI excludes 0: distinguishable"))
cat("\nper-biome (mean RMSE PCHIP/PPM/PIQS ; var-fraction):\n")
for(bb in c("tropics","temperate","boreal/polar")){ m<-land & band==bb
  cat(sprintf("  %-13s %.3f / %.3f / %.3f ; vfrac %.3f\n", bb, mean(rp[m]), mean(rpp[m]), mean(rq[m]), median(vfrac[m],na.rm=T))) }
