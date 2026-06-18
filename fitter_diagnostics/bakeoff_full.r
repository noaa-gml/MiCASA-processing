#!/usr/bin/env Rscript
## Comprehensive fitter bake-off on the REAL monthly means (recovered exactly
## from the shipped PCHIP fit). Scores every integral-preserving candidate on
## overshoot, within-month gradient retained, and month-edge discontinuity.
suppressWarnings(load("fit.piqs.rda"))
t0<-as.POSIXct("2001-01-01",tz="UTC")
edges<-as.numeric(seq(t0,by="1 month",length.out=length(piqsfit.time)+1)); D<-diff(edges)

mm  <- function(a,b) 0.5*(sign(a)+sign(b))*pmin(abs(a),abs(b))            # minmod
mm3 <- function(a,b,cc){s<-sign(a); ifelse(s==sign(b)&s==sign(cc), s*pmin(abs(a),abs(b),abs(cc)),0)}
maxmod <- function(a,b) 0.5*(sign(a)+sign(b))*pmax(abs(a),abs(b))

peakvals <- function(a,b,c0,Dm){                 # peak |f| within each piece
  fL<-c0; fR<-a*Dm^2+b*Dm+c0
  sv<-ifelse(a!=0,-b/(2*a),-1); intr<-a!=0&sv>0&sv<Dm; fV<-ifelse(intr,c0-b^2/(4*a),0)
  pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))
}

score <- function(fit, label){
  N<-360*180; M<-length(D)
  a<-matrix(fit$a,N,M); b<-matrix(fit$b,N,M); c0<-matrix(fit$c,N,M); Dm<-matrix(D,N,M,byrow=TRUE)
  u <- a*Dm^2/3 + b*Dm/2 + c0
  um<-abs(u); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M]))
  keep <- (rowSums(um)>1e-15) & (env>1e-12)
  uL<-cbind(u[,1],u[,1:(M-1)]); uR<-cbind(u[,2:M],u[,M])
  dL<-u-uL; dR<-uR-u

  methods <- list()
  ## edge-value pairs (Ledge, Redge) per method -> derive a,b,c (linear: a=0)
  lin <- function(slopehalf){ list(L=u-slopehalf, R=u+slopehalf, a=matrix(0,N,M)) }
  methods[["const"]]    <- lin(0*u)
  methods[["minmod"]]   <- lin(0.5*mm(dL,dR))
  methods[["vanLeer"]]  <- lin(0.5*ifelse(dL*dR>0, 2*dL*dR/(dL+dR+1e-300), 0))
  methods[["MC"]]       <- lin(0.5*mm3((dL+dR)/2, 2*dL, 2*dR))
  methods[["superbee"]] <- lin(0.5*maxmod(mm(dL,2*dR), mm(2*dL,dR)))

  ## PPM (vectorized)
  dc<-(uR-uL)/2; dm<-mm3(dc,2*dL,2*dR); dmp1<-cbind(dm[,2:M],dm[,M])
  aedge<-u+0.5*dR-(dmp1-dm)/6
  aLp<-cbind(u[,1],aedge[,1:(M-1)]); aRp<-aedge
  aLp[,1]<-u[,1];aRp[,1]<-u[,1];aLp[,M]<-u[,M];aRp[,M]<-u[,M]
  ext<-(aRp-u)*(u-aLp)<=0; aLp[ext]<-u[ext]; aRp[ext]<-u[ext]
  dd<-aRp-aLp; midd<-u-0.5*(aLp+aRp)
  hi<-(!ext)&(dd*midd> dd^2/6); aLp[hi]<-3*u[hi]-2*aRp[hi]
  lo<-(!ext)&(-(dd^2)/6> dd*midd); aRp[lo]<-3*u[lo]-2*aLp[lo]
  a6<-6*(u-0.5*(aLp+aRp))
  methods[["PPM"]] <- list(L=aLp, R=aRp, a=-a6/Dm^2, ppmb=(aRp-aLp+a6)/Dm, ppmc=aLp)

  cat(sprintf("\n==== %s (%d land cell-pieces) ====\n", label, sum(keep)))
  cat(sprintf("%-9s %8s %8s %9s %10s %10s %8s\n","method","max o/s","%os>1","grad/env","jump med","jump99","%disc"))
  ## PCHIP reference (shipped quadratic coeffs)
  pk<-peakvals(a,b,c0,Dm); r<-(pk/env)[keep]
  ## PCHIP knot discontinuity is 0 (C1). gradient = within-piece deviation amplitude
  grP<-(abs(pk-abs(u))/env)[keep]
  cat(sprintf("%-9s %8.3f %8.3f %9.3f %10s %10s %8s\n","PCHIP",max(r),100*mean(r>1.0001),median(grP,na.rm=T),"0","0","0"))

  for(nm in names(methods)){
    mth<-methods[[nm]]
    if(is.null(mth$ppmb)){ aa<-mth$a; bb<-(mth$R-mth$L)/Dm; cc<-mth$L
    } else { aa<-mth$a; bb<-mth$ppmb; cc<-mth$ppmc }
    pk<-peakvals(aa,bb,cc,Dm); r<-(pk/env)[keep]
    grad<-(abs(mth$R-mth$L)/2/env)[keep]                       # half edge-to-edge span / env
    jL<-mth$R[,1:(M-1)]; jR<-mth$L[,2:M]; jump<-abs(jL-jR)/env[,1:(M-1)]
    jk<-jump[keep]; jk<-jk[is.finite(jk)]
    cat(sprintf("%-9s %8.3f %8.3f %9.3f %10.3f %10.3f %8.1f\n",
        nm, max(r), 100*mean(r>1.0001), median(grad,na.rm=T),
        median(jk), quantile(jk,.99), 100*mean(jk>1e-6)))
  }
}
score(piqsfit.gpp,"GPP (-2*NPP)")
score(piqsfit.resp,"RESP (Rh+NPP)")
