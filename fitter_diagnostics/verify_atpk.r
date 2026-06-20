## Verify the production ATP fit (fit.atpk.rda): structure, coherence/mass,
## variance arrays, sign-safety, and point-estimate agreement with PCHIP.
A<-new.env(); load("fit.atpk.rda",envir=A); P<-new.env(); load("fit.pchip.rda",envir=P)
N<-360*180; M<-length(A$piqsfit.time)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=M+1)); D<-diff(x); Dm<-matrix(D,N,M,byrow=TRUE)
cat("fitter:",A$piqsfit.meta$fitter,"  has $var:",!is.null(A$piqsfit.gpp$var),"  atpk knobs:",unlist(A$piqsfit.meta$atpk),"\n")
cat("dims gpp$a:",paste(dim(A$piqsfit.gpp$a),collapse="x"),"  var:",paste(dim(A$piqsfit.gpp$var),collapse="x"),"\n")
for(comp in c("piqsfit.gpp","piqsfit.resp")){
  a<-matrix(A[[comp]]$a,N,M); b<-matrix(A[[comp]]$b,N,M); c0<-matrix(A[[comp]]$c,N,M); v<-matrix(A[[comp]]$var,N,M)
  ap<-matrix(P[[comp]]$a,N,M); bp<-matrix(P[[comp]]$b,N,M); cp<-matrix(P[[comp]]$c,N,M)
  u <- a*Dm^2/3+b*Dm/2+c0; up<- ap*Dm^2/3+bp*Dm/2+cp        # monthly means from each fit
  um<-abs(up); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M])); keep<-rowSums(um)>1e-15
  ## coherence/mass: ATP fit's piece integral == PCHIP's (both == cat monthly mean)
  massrel<-max((abs(u-up)/pmax(env,1e-30))[keep,])
  ## variance
  vk<-v[keep,]; sdband<-(sqrt(v)/env)[keep,]
  ## point estimate vs PCHIP at sub-monthly points
  ss<-((1:6)-0.5)/6; rms<-0; np<-0
  for(q in ss){ dt<-q*Dm; fa<-a*dt^2+b*dt+c0; fp<-ap*dt^2+bp*dt+cp; rms<-rms+sum(((fa-fp)/env)[keep,]^2,na.rm=TRUE); np<-np+sum(keep)*M }
  rms<-sqrt(rms/np)
  ## sign-safety: GPP should be <=0, RESP >=0
  sgn<- if(comp=="piqsfit.gpp") -1 else 1
  fmax<-pmax(c0, a*Dm^2+b*Dm+c0); wrong<- if(sgn<0) (fmax>1e-9*env) else (pmin(c0,a*Dm^2+b*Dm+c0) < -1e-9*env)
  cat(sprintf("\n%s:\n  coherence(mass) max|u_atpk-u_pchip|/env = %.2e\n  var>=0:%s  median sd/env=%.3f  90th=%.3f  %% land-months var>0=%.0f%%\n  point-est RMS vs PCHIP /env = %.4f\n  sign-safe: wrong-sign land-months = %.3f%%\n",
      comp, massrel, all(v>=0), median(sdband,na.rm=T), quantile(sdband,.9,na.rm=T), 100*mean(vk>0), rms, 100*mean(wrong[keep,])))
}
