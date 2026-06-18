## Uncertainty + per-biome + PAIRED fidelity test (addresses the "no spread,
## 0.149 vs 0.151 needs an IQR" critique). Same cells across methods, so the
## honest test is the per-cell PAIRED difference PCHIP_rmse - PPM_rmse.
## NOTE: fit.piqs.rda is now PPM (default switched); read PCHIP from fit.pchip.rda.
suppressWarnings(suppressMessages(library(ncdf4)))
PC<-new.env(); load("fit.pchip.rda",envir=PC)
PP<-new.env(); load("fit.ppm.rda",envir=PP)
LM<-new.env(); load("fit.linmm.rda",envir=LM)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(PC$piqsfit.time)+1)); D<-diff(x)
N<-360*180; i2020<-(2020-2001)*12+(1:12); dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31)
lat<- -90 + ((rep(1:180,each=360)-1)%/%360 + ((1:N-1)%/%360))   # placeholder
j<- ((1:N-1)%/%360)+1; latc<- -90 + j - 0.5
band<- ifelse(abs(latc)<=23.5,"tropics",ifelse(abs(latc)<=50,"temperate","boreal/polar"))
ev<-function(E,comp,idx,dt){a<-as.vector(E[[comp]]$a[,,idx]);b<-as.vector(E[[comp]]$b[,,idx]);c0<-as.vector(E[[comp]]$c[,,idx]);a*dt^2+b*dt+c0}
um<-matrix(0,N,12)
for(k in 1:12){idx<-i2020[k]; um[,k]<-as.vector(PC$piqsfit.gpp$a[,,idx])*D[idx]^2/3+as.vector(PC$piqsfit.gpp$b[,,idx])*D[idx]/2+as.vector(PC$piqsfit.gpp$c[,,idx])}
env<-pmax(cbind(abs(um[,1]),abs(um[,1:11])),abs(um),cbind(abs(um[,2:12]),abs(um[,12])))
rmsePC<-rmsePP<-rmseLM<-bnd<-c()
for(k in 1:12){nd<-dpm[k];mm<-sprintf("%02d",k);idx<-i2020[k];dts<-((1:nd)-0.5)*86400
  g<-matrix(0,N,nd)
  for(dd in 1:nd){f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%s%02d.nc",mm,dd));g[,dd]<- -2*as.vector(ncvar_get(f,"NPP"));nc_close(f)}
  rc<-sapply(dts,function(dt)ev(PC,"piqsfit.gpp",idx,dt)); rp<-sapply(dts,function(dt)ev(PP,"piqsfit.gpp",idx,dt)); rl<-sapply(dts,function(dt)ev(LM,"piqsfit.gpp",idx,dt))
  rm<-function(X)sqrt(rowMeans((X-g)^2))
  keep<-env[,k]>1e-12 & rowSums(abs(g))>1e-15
  rmsePC<-c(rmsePC,(rm(rc)/env[,k])[keep]); rmsePP<-c(rmsePP,(rm(rp)/env[,k])[keep]); rmseLM<-c(rmseLM,(rm(rl)/env[,k])[keep]); bnd<-c(bnd,band[keep])}
q<-function(v)sprintf("med %.3f  IQR[%.3f,%.3f]  mean %.3f",median(v),quantile(v,.25),quantile(v,.75),mean(v))
cat("GPP daily RMSE/env (cell-months):\n")
cat(sprintf("  PCHIP  %s\n  PPM    %s\n  minmod %s\n", q(rmsePC),q(rmsePP),q(rmseLM)))
cat("\nPAIRED PCHIP-PPM per cell-month (positive => PPM better):\n")
d<-rmsePC-rmsePP
cat(sprintf("  median Δ %.4f  IQR[%.4f,%.4f]  %% cell-months PPM better %.1f%%  (median |rmse| ~%.2f so Δ is ~%.1f%% of level)\n",
    median(d),quantile(d,.25),quantile(d,.75),100*mean(d>0),median(rmsePC),100*median(d)/median(rmsePC)))
cat("\nPer-biome median RMSE/env (PCHIP / PPM / minmod) and PPM-better %:\n")
for(bb in c("tropics","temperate","boreal/polar")){m<-bnd==bb
  cat(sprintf("  %-13s n=%6d  %.3f / %.3f / %.3f   PPM-better %.1f%%\n",
      bb,sum(m),median(rmsePC[m]),median(rmsePP[m]),median(rmseLM[m]),100*mean((rmsePC-rmsePP)[m]>0)))}
