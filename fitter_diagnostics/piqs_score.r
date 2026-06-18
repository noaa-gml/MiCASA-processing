## Fill the PIQS gap in the scorecard: overshoot peak/env (from coeffs) and
## daily-fidelity RMSE/env vs MiCASA daily 2020 (same basis as fidelity_daily.r).
suppressWarnings(suppressMessages(library(ncdf4)))
V<-new.env(); load("fit.piqs_v1.rda",envir=V)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(V$piqsfit.time)+1)); D<-diff(x)
N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
peak_q<-function(a,b,c0,DD){fL<-c0;fR<-a*DD^2+b*DD+c0;sv<-ifelse(a!=0,-b/(2*a),-1);intr<-a!=0&sv>0&sv<DD;fV<-ifelse(intr,c0-b^2/(4*a),0);pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))}
## ---- overshoot + sign-flip from coeffs ----
for(comp in c("piqsfit.gpp","piqsfit.resp")){
  a<-matrix(V[[comp]]$a,N,M); b<-matrix(V[[comp]]$b,N,M); c0<-matrix(V[[comp]]$c,N,M)
  u<-a*Dm^2/3+b*Dm/2+c0; um<-abs(u)
  env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M])); keep<-rowSums(um)>1e-15
  pk<-peak_q(a,b,c0,Dm); r<-(pk/env)[keep,]; r<-r[is.finite(r)]
  ## sign-flip: GPP flux should be <=0, RESP >=0. count cell-months with wrong-sign anywhere in piece
  wrong<- if(comp=="piqsfit.gpp") (pmax(c0,a*Dm^2+b*Dm+c0)>1e-12*env) else (pmin(c0,a*Dm^2+b*Dm+c0) < -1e-12*env)
  cat(sprintf("PIQS %s: overshoot peak/env med %.2f 99th %.2f max %.2f ; %% cell-months wrong-sign %.2f%%\n",
      comp, median(r), quantile(r,.99), max(r), 100*mean(wrong[keep,],na.rm=TRUE)))
}
## ---- daily fidelity vs MiCASA daily 2020 ----
i2020<-(2020-2001)*12+(1:12); dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31)
ev<-function(comp,idx,dt){a<-as.vector(V[[comp]]$a[,,idx]);b<-as.vector(V[[comp]]$b[,,idx]);c0<-as.vector(V[[comp]]$c[,,idx]);a*dt^2+b*dt+c0}
um_g<-um_r<-matrix(0,N,12)
for(k in 1:12){idx<-i2020[k]
  um_g[,k]<-as.vector(V$piqsfit.gpp$a[,,idx])*D[idx]^2/3+as.vector(V$piqsfit.gpp$b[,,idx])*D[idx]/2+as.vector(V$piqsfit.gpp$c[,,idx])
  um_r[,k]<-as.vector(V$piqsfit.resp$a[,,idx])*D[idx]^2/3+as.vector(V$piqsfit.resp$b[,,idx])*D[idx]/2+as.vector(V$piqsfit.resp$c[,,idx])}
env_g<-pmax(cbind(abs(um_g[,1]),abs(um_g[,1:11])),abs(um_g),cbind(abs(um_g[,2:12]),abs(um_g[,12])))
env_r<-pmax(cbind(abs(um_r[,1]),abs(um_r[,1:11])),abs(um_r),cbind(abs(um_r[,2:12]),abs(um_r[,12])))
rg<-rr<-c()
for(k in 1:12){nd<-dpm[k];mm<-sprintf("%02d",k);idx<-i2020[k];dts<-((1:nd)-0.5)*86400
  g<-matrix(0,N,nd);r<-matrix(0,N,nd)
  for(dd in 1:nd){f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%s%02d.nc",mm,dd))
    npp<-as.vector(ncvar_get(f,"NPP"));rh<-as.vector(ncvar_get(f,"Rh"));nc_close(f);g[,dd]<- -2*npp;r[,dd]<-rh+npp}
  rePg<-sapply(dts,function(dt)ev("piqsfit.gpp",idx,dt)); rePr<-sapply(dts,function(dt)ev("piqsfit.resp",idx,dt))
  rmse<-function(X,T)sqrt(rowMeans((X-T)^2))
  kg<-env_g[,k]>1e-12&rowSums(abs(g))>1e-15; kr<-env_r[,k]>1e-12&rowSums(abs(r))>1e-15
  rg<-c(rg,(rmse(rePg,g)/env_g[,k])[kg]); rr<-c(rr,(rmse(rePr,r)/env_r[,k])[kr])}
cat(sprintf("PIQS daily RMSE/env: GPP mean %.4f median %.3f ; RESP mean %.4f median %.3f\n",
    mean(rg),median(rg),mean(rr),median(rr)))
cat("   (compare: GPP mean PPM 0.149 PCHIP 0.151 minmod 0.159 flat 0.181 ; RESP PPM 0.125 PCHIP 0.128)\n")
