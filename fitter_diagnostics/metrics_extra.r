#!/usr/bin/env Rscript
## Extra rating metrics (coefficient-based, fast):
##  (A) smoothness: knot value-jump (C0), knot first-derivative jump (C1),
##      within-piece curvature; (B) spurious wiggle: % pieces with an interior
##      extremum; (C) NRT-revision locality: perturb last month, footprint.
source("lib/pchip_fit.r"); source("lib/ppm_fit.r"); source("lib/linmm_fit.r")
P<-new.env(); load("fit.piqs.rda",envir=P); M<-new.env(); load("fit.ppm.rda",envir=M); L<-new.env(); load("fit.linmm.rda",envir=L)
t0<-as.POSIXct("2001-01-01",tz="UTC"); edges<-as.numeric(seq(t0,by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(edges)
N<-360*180; nm<-length(D)

smooth_extrema<-function(fit,comp,lab){
  a<-matrix(fit[[comp]]$a,N,nm); b<-matrix(fit[[comp]]$b,N,nm); c0<-matrix(fit[[comp]]$c,N,nm); Dm<-matrix(D,N,nm,byrow=TRUE)
  u<-a*Dm^2/3+b*Dm/2+c0; um<-abs(u)
  env<-pmax(cbind(um[,1],um[,1:(nm-1)]),um,cbind(um[,2:nm],um[,nm])); env[env<1e-12]<-NA; keep<-rowSums(um)>1e-15
  fRk<-a*Dm^2+b*Dm+c0                          # f at right knot of piece k
  vjump<-abs(fRk[,1:(nm-1)]-c0[,2:nm])/env[,1:(nm-1)]          # C0: value jump
  fpR<-2*a*Dm+b; fpL<-b                                         # f' at right / left
  djump<-abs((fpR[,1:(nm-1)]-fpL[,2:nm])*Dm[,1:(nm-1)])/env[,1:(nm-1)]  # C1: deriv jump * D /env
  curv<-abs(2*a*Dm^2)/env                                      # within-piece curvature (flux units)
  sv<- -b/(2*a); intr<- a!=0 & sv>0 & sv<Dm                    # interior extremum?
  g<-function(x){x<-x[keep[1:N],]; x<-x[is.finite(x)]; x}
  vj<-vjump[keep[1:N],]; vj<-vj[is.finite(vj)]; dj<-djump[keep[1:N],]; dj<-dj[is.finite(dj)]
  cv<-curv[keep,]; cv<-cv[is.finite(cv)]; iv<-intr[keep,]
  cat(sprintf("%-8s  Cjump(val)/env med %.3f  C1 deriv-jump*D/env med %.3f  curvature/env med %.3f  %%pieces interior-extremum %.1f\n",
      lab, median(vj), median(dj), median(cv), 100*mean(iv,na.rm=TRUE)))
}
cat("== (A) smoothness + (B) spurious wiggle, GPP ==\n")
smooth_extrema(P,"piqsfit.gpp","PCHIP"); smooth_extrema(M,"piqsfit.gpp","PPM"); smooth_extrema(L,"piqsfit.gpp","minmod")
cat("== smoothness + spurious wiggle, RESP ==\n")
smooth_extrema(P,"piqsfit.resp","PCHIP"); smooth_extrema(M,"piqsfit.resp","PPM"); smooth_extrema(L,"piqsfit.resp","minmod")

## ---- (C) NRT-revision locality: perturb LAST month +10%, footprint ----
set.seed(7)
um_all<-matrix(P$piqsfit.gpp$a,N,nm)*matrix(D,N,nm,byrow=TRUE)^2/3 + matrix(P$piqsfit.gpp$b,N,nm)*matrix(D,N,nm,byrow=TRUE)/2 + matrix(P$piqsfit.gpp$c,N,nm)
land<-which(rowSums(abs(um_all))>1e-12 & abs(um_all[,nm])>1e-12)
samp<-sample(land, min(1500,length(land)))
x<-edges
foot<-function(fitfun){
  fps<-integer(length(samp))
  for(s in seq_along(samp)){ yb<-um_all[samp[s],]
    f0<-fitfun(x,yb); yb2<-yb; yb2[nm]<-yb2[nm]*1.10; f1<-fitfun(x,yb2)
    Dm<-D; e<-max(abs(yb)); if(e<=0){fps[s]<-0;next}
    ## per-month max |Δflux| at knots, normalized; footprint = # months back changed >1%
    dc<-abs(f1$c-f0$c)/e; chg<-which(dc>0.01)
    fps[s]<- if(length(chg)) nm-min(chg) else 0 }   # months before the revised one that changed
  fps }
cat("\n== (C) NRT-revision locality: perturb last month +10%, # PRIOR months whose flux moves >1% ==\n")
for(nmlab in c("PCHIP","PPM","minmod")){
  ff<-switch(nmlab, PCHIP=pchip.fit.cell, PPM=ppm.fit.cell, minmod=linmm.fit.cell)
  fp<-foot(ff); cat(sprintf("%-8s  footprint months: median %d  90th %d  max %d  (PIQS global solve = all %d)\n",
      nmlab, median(fp), quantile(fp,.9), max(fp), nm-1)) }
