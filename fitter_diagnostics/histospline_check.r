## Histospline = unconstrained cubic-on-cumulative (the parent class; PCHIP is
## its monotone-limited member, MSS its positivity-constrained member). Measure
## its overshoot / sign-flip on a sample to confirm it is dominated.
P<-new.env(); load("fit.ppm.rda",envir=P)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(x)
N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
a<-matrix(P$piqsfit.gpp$a,N,M); b<-matrix(P$piqsfit.gpp$b,N,M); c0<-matrix(P$piqsfit.gpp$c,N,M)
U<-a*Dm^2/3+b*Dm/2+c0
land<-which(rowSums(abs(U))>1e-15)
set.seed(5); samp<-sample(land,300)
os<-flip<-numeric(length(samp))
for(si in seq_along(samp)){
  yb<-U[samp[si],]; Fk<-c(0,cumsum(yb*D))
  fn<-splinefun(x,Fk,method="natural")
  tt<-unlist(lapply(1:M,function(k)seq(x[k],x[k+1],length.out=20)))
  fl<-fn(tt,deriv=1)
  e<-pmax(c(abs(yb)[1],abs(yb)[1:(M-1)]),abs(yb),c(abs(yb)[2:M],abs(yb)[M]))
  mth<-pmin(M,findInterval(tt,x,rightmost.closed=TRUE)); ok<-e[mth]>1e-12
  os[si]<-max(abs(fl[ok])/e[mth][ok])
  flip[si]<-as.integer(max(fl) > 1e-12*max(abs(yb)))
}
cat(sprintf("Cubic histospline (natural, unconstrained): overshoot peak/env median %.2f max %.2f ; %% cells wrong-sign %.1f%%\n",
            median(os),max(os),100*mean(flip)))
cat("   vs PCHIP 1.50 (monotone-limited member) and MSS 1.35 / 24%% (constrained member)\n")
