## Prototype prior-uncertainty bands for the sub-monthly reconstruction (the
## thing the deterministic splines lack). Two cheap, LOCALITY-PRESERVING bands:
##  (A) STRUCTURAL: spread across the mass-preserving fitters {PCHIP,PPM,minmod}
##      -- "how much does the (arbitrary) smoother choice move the prior?"
##  (B) BOOTSTRAP-PCHIP: resample days within each month -> bootstrap monthly
##      means -> refit PCHIP (B times) -> 5-95% envelope (input/data uncertainty).
## Both stay per-cell + windowed, so they keep the NRT locality PCHIP has.
suppressWarnings(suppressMessages(library(ncdf4)))
source("lib/pchip_fit.r")
PC<-new.env(); load("fit.pchip.rda",envir=PC); PP<-new.env(); load("fit.ppm.rda",envir=PP); LM<-new.env(); load("fit.linmm.rda",envir=LM)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(PC$piqsfit.time)+1)); D<-diff(x); N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
a<-matrix(PC$piqsfit.gpp$a,N,M); b<-matrix(PC$piqsfit.gpp$b,N,M); c0<-matrix(PC$piqsfit.gpp$c,N,M); U<-a*Dm^2/3+b*Dm/2+c0
latof<-function(cl) -90 + ((cl-1)%/%360)+1 - 0.5
band<-function(l) ifelse(abs(l)<=23.5,"tropics",ifelse(abs(l)<=50,"temperate","boreal/polar"))
land<-which(rowSums(abs(U))>1e-15); ll<-latof(land)
set.seed(8); samp<-c(sample(land[abs(ll)<23.5],10),sample(land[abs(ll)>=23.5&abs(ll)<50],10),sample(land[abs(ll)>=50],10)); ns<-length(samp)
Ns<-6; ss<-(seq_len(Ns)-0.5)/Ns
gv<-function(E,i3,j3,idx,dt) E[[ "piqsfit.gpp" ]]$a[i3,j3,idx]*dt^2 + E[[ "piqsfit.gpp" ]]$b[i3,j3,idx]*dt + E[[ "piqsfit.gpp" ]]$c[i3,j3,idx]

## ---- (A) structural spread across 3 fitters ----
bw<-bnd<-c()
for(s in 1:ns){ cl<-samp[s]; bbc<-band(latof(cl)); i3<-((cl-1)%%360)+1; j3<-((cl-1)%/%360)+1
  for(idx in 6:(M-6)){ env<-max(abs(U[cl,(idx-1):(idx+1)])); if(env<1e-12)next
    for(fr in ss){ dt<-fr*D[idx]
      vals<-c(gv(PC,i3,j3,idx,dt),gv(PP,i3,j3,idx,dt),gv(LM,i3,j3,idx,dt))
      bw<-c(bw,(max(vals)-min(vals))/env); bnd<-c(bnd,bbc) } } }
cat("== (A) STRUCTURAL band (max-min across PCHIP/PPM/minmod) / local envelope ==\n")
cat(sprintf("  overall: median %.3f  90th %.3f  99th %.3f\n",median(bw),quantile(bw,.9),quantile(bw,.99)))
for(bb in c("tropics","temperate","boreal/polar")) cat(sprintf("  %-13s median %.3f  90th %.3f\n",bb,median(bw[bnd==bb]),quantile(bw[bnd==bb],.9)))

## ---- (B) bootstrap-PCHIP on 2020 daily ----
dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31); nb<-12
edges20<-as.numeric(seq(as.POSIXct("2020-01-01",tz="UTC"),by="1 month",length.out=nb+1)); D20<-diff(edges20)
cellsB<-samp[c(3,7,13,17,23,27)]                       # 2 tropics, 2 temperate, 2 boreal
DAY<-vector("list",nb)
for(m in 1:nb){ arr<-array(0,c(length(cellsB),dpm[m]))
  for(dd in 1:dpm[m]){ f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%02d%02d.nc",m,dd)); npp<-ncvar_get(f,"NPP"); nc_close(f)
    for(ci in seq_along(cellsB)){ i3<-((cellsB[ci]-1)%%360)+1; j3<-((cellsB[ci]-1)%/%360)+1; arr[ci,dd]<- -2*npp[i3,j3] } }
  DAY[[m]]<-arr }
B<-200; intm<-3:10
cat("\n== (B) BOOTSTRAP-PCHIP 5-95% band (resample days->monthly mean->refit PCHIP), 2020 ==\n")
for(ci in seq_along(cellsB)){ l<-latof(cellsB[ci])
  curves<-matrix(0,B,nb*Ns)
  for(bk in 1:B){ mm<-sapply(1:nb,function(m){d<-DAY[[m]][ci,]; mean(d[sample(length(d),replace=TRUE)])})
    f<-pchip.fit.cell(edges20,mm); pt<-1
    for(m in 1:nb) for(fr in ss){ curves[bk,pt]<-f$a[m]*(fr*D20[m])^2+f$b[m]*(fr*D20[m])+f$c[m]; pt<-pt+1 } }
  cmean<-sapply(1:nb,function(m)mean(DAY[[m]][ci,])); env<-max(abs(cmean)); if(env<1e-12){cat(sprintf("  lat %+5.1f: ~zero flux\n",l));next}
  idxint<-as.vector(sapply(intm,function(m)((m-1)*Ns+1):(m*Ns)))
  lo<-apply(curves[,idxint],2,quantile,.05); hi<-apply(curves[,idxint],2,quantile,.95)
  cat(sprintf("  lat %+5.1f (%-12s): band/env  median %.3f  max %.3f\n",l,band(l),median((hi-lo)/env),max((hi-lo)/env))) }
