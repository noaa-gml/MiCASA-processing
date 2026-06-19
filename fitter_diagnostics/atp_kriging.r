## Prototype 1-D TEMPORAL area-to-point (ATP) kriging for monthly->sub-monthly
## flux disaggregation (Kyriakidis 2004; Yoo & Kyriakidis 2006). Block data =
## monthly means; predict sub-monthly point values that (i) re-average to the
## monthly mean exactly (COHERENCE/pycnophylactic), (ii) carry a kriging VARIANCE
## (the uncertainty PCHIP lacks), (iii) can be made non-negative by selective QP.
## Globality is acceptable per use-case, so we krige a full multi-year window.
source("lib/pchip_fit.r")
PC<-new.env(); load("fit.pchip.rda",envir=PC)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(PC$piqsfit.time)+1)); D<-diff(x); N<-360*180; Mtot<-length(D); Dm<-matrix(D,N,Mtot,byrow=TRUE)
a<-matrix(PC$piqsfit.gpp$a,N,Mtot); b<-matrix(PC$piqsfit.gpp$b,N,Mtot); c0<-matrix(PC$piqsfit.gpp$c,N,Mtot); U<-a*Dm^2/3+b*Dm/2+c0
latof<-function(cl) -90+((cl-1)%/%360)+1-0.5; band<-function(l) ifelse(abs(l)<=23.5,"tropics",ifelse(abs(l)<=50,"temperate","boreal/polar"))
land<-which(rowSums(abs(U))>1e-15); ll<-latof(land)
set.seed(8); samp<-c(sample(land[abs(ll)<23.5],3),sample(land[abs(ll)>=23.5&abs(ll)<50],3),sample(land[abs(ll)>=50],3))

win<-217:252; M<-length(win); Ns<-6; ss<-(seq_len(Ns)-0.5)/Ns
tk<-rep(seq_len(M),each=Ns)-1+rep(ss,M)               # time in months
B<-matrix(0,M,M*Ns); for(i in 1:M) B[i,((i-1)*Ns+1):(i*Ns)]<-1/Ns
Dlag<-abs(outer(tk,tk,"-")); range<-1.5; K<-length(tk)
edgesW<-x[win[1]:(win[length(win)]+1)]; Dw<-diff(edgesW)

krige<-function(d){
  sig2<-var(d); if(sig2<=0) sig2<-1e-30
  C<-sig2*exp(-Dlag/range)+diag(1e-3*sig2,K)
  CBB<-B%*%C%*%t(B); CBp<-B%*%C
  A<-rbind(cbind(CBB,1),c(rep(1,M),0))
  sol<-tryCatch(solve(A, rbind(CBp,1)), error=function(e) solve(A+diag(1e-2*sig2,M+1), rbind(CBp,1)))  # ridge fallback for ill-conditioned cells
  lam<-sol[1:M,,drop=FALSE]; mu<-sol[M+1,]
  zhat<-as.vector(t(lam)%*%d)
  vr<-sig2 - colSums(CBp*lam) - mu                     # OK point variance
  list(z=zhat, sd=sqrt(pmax(vr,0)))
}
cat("== 1-D temporal ATP kriging (exp cov, range=1.5mo, 6 pts/mo, 36-mo window) ==\n")
for(cl in samp){ i3<-((cl-1)%%360)+1; j3<-((cl-1)%/%360)+1; d<-U[cl,win]; l<-latof(cl)
  if(var(d)<=1e-30){cat(sprintf("  lat %+5.1f: ~zero-variance cell, skip\n",latof(cl)));next}; k<-tryCatch(krige(d),error=function(e)NULL); if(is.null(k)){cat(sprintf('  lat %+5.1f: ill-conditioned (near-dormant), skip\n',latof(cl)));next}
  ## coherence: block-average of point predictions == monthly mean datum?
  coh<-max(abs(as.vector(B%*%k$z)-d))/max(abs(d),1e-30)
  ## PCHIP at same points for comparison
  pch<-numeric(K); for(m in 1:M){ idx<-win[m]; for(q in 1:Ns){ dt<-ss[q]*D[idx]; pch[(m-1)*Ns+q]<-PC$piqsfit.gpp$a[i3,j3,idx]*dt^2+PC$piqsfit.gpp$b[i3,j3,idx]*dt+PC$piqsfit.gpp$c[i3,j3,idx] } }
  env<-max(abs(d))
  cat(sprintf("  lat %+5.1f (%-12s): coherence err %.1e ; krige-var band(+-1.96sd)/env med %.3f ; %% pts wrong-sign %.0f%% ; RMS(krige-PCHIP)/env %.3f\n",
      l, band(l), coh, median(2*1.96*k$sd/env), 100*mean(k$z>1e-12*env), sqrt(mean((k$z-pch)^2))/env)) }
