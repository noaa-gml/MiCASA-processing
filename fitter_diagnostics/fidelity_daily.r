#!/usr/bin/env Rscript
## METRIC: fidelity of each fitter's sub-monthly reconstruction to MiCASA's
## OWN daily 1-degree data (2020). Decides whether PCHIP's overshoot bump is
## real sub-monthly structure (daily exceeds the monthly-mean envelope) or an
## artifact (daily stays bounded -> PPM/linear faithful).
suppressWarnings(suppressMessages(library(ncdf4)))
P<-new.env(); load("fit.piqs.rda",envir=P)
M<-new.env(); load("fit.ppm.rda",envir=M)
L<-new.env(); load("fit.linmm.rda",envir=L)
t0<-as.POSIXct("2001-01-01",tz="UTC")
edges<-as.numeric(seq(t0,by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(edges)
i2020<-(2020-2001)*12 + (1:12)                    # fit month indices for 2020
dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31)

## monthly means per cell (from PCHIP integral; identical across fits)
Ng<-360*180
um_g<-um_r<-matrix(0,Ng,12)
for(k in 1:12){ idx<-i2020[k]
  ag<-as.vector(P$piqsfit.gpp$a[,,idx]); bg<-as.vector(P$piqsfit.gpp$b[,,idx]); cg<-as.vector(P$piqsfit.gpp$c[,,idx])
  ar<-as.vector(P$piqsfit.resp$a[,,idx]);br<-as.vector(P$piqsfit.resp$b[,,idx]);cr<-as.vector(P$piqsfit.resp$c[,,idx])
  um_g[,k]<-ag*D[idx]^2/3+bg*D[idx]/2+cg; um_r[,k]<-ar*D[idx]^2/3+br*D[idx]/2+cr }
env_g<-pmax(cbind(abs(um_g[,1]),abs(um_g[,1:11])),abs(um_g),cbind(abs(um_g[,2:12]),abs(um_g[,12])))
env_r<-pmax(cbind(abs(um_r[,1]),abs(um_r[,1:11])),abs(um_r),cbind(abs(um_r[,2:12]),abs(um_r[,12])))

evalfit<-function(fit,comp,idx,dt) {            # flux at offset dt (sec) in month idx
  a<-as.vector(fit[[comp]]$a[,,idx]); b<-as.vector(fit[[comp]]$b[,,idx]); c0<-as.vector(fit[[comp]]$c[,,idx])
  a*dt^2+b*dt+c0 }

## accumulators: truth daily-max/env, and per-fitter daily RMSE & recon-peak/env
acc<-function() list(dmax=c(), rmseP=c(), rmseM=c(), rmseL=c(), rmseC=c(), pkP=c(), pkM=c(), pkL=c(), keep=c())
AG<-acc(); AR<-acc()
for(k in 1:12){
  nd<-dpm[k]; mm<-sprintf("%02d",k)
  g<-array(0,c(Ng,nd)); r<-array(0,c(Ng,nd))
  for(dd in 1:nd){ f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%s%02d.nc",mm,dd))
    npp<-as.vector(ncvar_get(f,"NPP")); rh<-as.vector(ncvar_get(f,"Rh")); nc_close(f)
    g[,dd]<- -2*npp; r[,dd]<- rh + npp }   # gC units, matching fit coeffs
  ## diurnalize uses gpp=-2*NPP/12, rtot=(Rh - 0.5*gpp)/12 with gpp negative => rtot=Rh/12 + NPP/12; align:
  ## recompute r consistently: rtot = Rh/12 - 0.5*gpp ; gpp=-2NPP/12 => -0.5*gpp = NPP/12 => rtot=(Rh+NPP)/12  OK
  idx<-i2020[k]
  dts<-((1:nd)-0.5)*86400                       # day midpoints from month start
  ## per fitter recon: build [Ng,nd]
  reconP_g<-reconM_g<-reconL_g<-matrix(0,Ng,nd); reconP_r<-reconM_r<-reconL_r<-matrix(0,Ng,nd)
  for(dd in 1:nd){ dt<-dts[dd]
    reconP_g[,dd]<-evalfit(P,"piqsfit.gpp",idx,dt); reconM_g[,dd]<-evalfit(M,"piqsfit.gpp",idx,dt); reconL_g[,dd]<-evalfit(L,"piqsfit.gpp",idx,dt)
    reconP_r[,dd]<-evalfit(P,"piqsfit.resp",idx,dt);reconM_r[,dd]<-evalfit(M,"piqsfit.resp",idx,dt);reconL_r[,dd]<-evalfit(L,"piqsfit.resp",idx,dt) }
  rms<-function(x,truth) sqrt(rowMeans((x-truth)^2))
  keepg<- env_g[,k]>1e-12 & rowSums(abs(g))>1e-15
  keepr<- env_r[,k]>1e-12 & rowSums(abs(r))>1e-15
  AG$dmax<-c(AG$dmax,(apply(abs(g),1,max)/env_g[,k])[keepg])
  AG$rmseP<-c(AG$rmseP,(rms(reconP_g,g)/env_g[,k])[keepg]); AG$rmseM<-c(AG$rmseM,(rms(reconM_g,g)/env_g[,k])[keepg])
  AG$rmseL<-c(AG$rmseL,(rms(reconL_g,g)/env_g[,k])[keepg]); AG$rmseC<-c(AG$rmseC,(rms(matrix(um_g[,k],Ng,nd),g)/env_g[,k])[keepg])
  AG$pkP<-c(AG$pkP,(apply(abs(reconP_g),1,max)/env_g[,k])[keepg]); AG$pkM<-c(AG$pkM,(apply(abs(reconM_g),1,max)/env_g[,k])[keepg]); AG$pkL<-c(AG$pkL,(apply(abs(reconL_g),1,max)/env_g[,k])[keepg])
  AR$dmax<-c(AR$dmax,(apply(abs(r),1,max)/env_r[,k])[keepr])
  AR$rmseP<-c(AR$rmseP,(rms(reconP_r,r)/env_r[,k])[keepr]); AR$rmseM<-c(AR$rmseM,(rms(reconM_r,r)/env_r[,k])[keepr])
  AR$rmseL<-c(AR$rmseL,(rms(reconL_r,r)/env_r[,k])[keepr]); AR$rmseC<-c(AR$rmseC,(rms(matrix(um_r[,k],Ng,nd),r)/env_r[,k])[keepr])
  cat(sprintf("month %02d done (land g=%d r=%d)\n",k,sum(keepg),sum(keepr)))
}
rep<-function(A,lab){ q<-function(x)sprintf("%.3f/%.3f/%.3f",median(x),quantile(x,.9),quantile(x,.99))
  cat(sprintf("\n==== %s ====\n",lab))
  cat(sprintf("TRUTH daily-max/env (med/90/99): %s   <- real sub-monthly excursion\n", q(A$dmax)))
  cat(sprintf("recon-peak/env  PCHIP %s  PPM %s  minmod %s\n", q(A$pkP),q(A$pkM),q(A$pkL)))
  cat(sprintf("daily RMSE/env  PCHIP %s\n                PPM   %s\n                minmod%s\n                flat  %s\n",
              q(A$rmseP),q(A$rmseM),q(A$rmseL),q(A$rmseC)))
  cat(sprintf("mean daily RMSE/env: PCHIP %.4f  PPM %.4f  minmod %.4f  flat %.4f\n",
              mean(A$rmseP),mean(A$rmseM),mean(A$rmseL),mean(A$rmseC))) }
rep(AG,"GPP fidelity vs MiCASA daily 2020")
rep(AR,"RESP fidelity vs MiCASA daily 2020")
