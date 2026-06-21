## C2 fix: the ATP kriging variance band magnitude is set by the covariance
## RANGE, a hand-chosen knob with no data-driven way to set it. Sweep range and
## show the "9-52%" band moves with it; also confirm the variance is essentially
## data-independent geometry x sill (so the band is a chosen prior, not a measured
## sub-monthly indeterminacy).
source("lib/atpk_fit.r")
P<-new.env(); load("fit.pchip.rda",envir=P)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(x); N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
a<-matrix(P$piqsfit.gpp$a,N,M); b<-matrix(P$piqsfit.gpp$b,N,M); c0<-matrix(P$piqsfit.gpp$c,N,M); U<-a*Dm^2/3+b*Dm/2+c0
land<-which(rowSums(abs(U))>1e-15); set.seed(8); samp<-sample(land,300)
env_cell<-function(cl) max(abs(U[cl,]))                    # per-cell scale (robust)

cat("== C2: ATP kriging variance band vs covariance range (the hand-set knob) ==\n")
cat("   range(mo) | median band(+-1.96sd)/env  | 90th\n")
for(r in c(0.5,0.75,1.5,3.0,6.0)){
  ww<-atpk.window.weights(M, W=6, Ns=6, range=r); bw<-c()
  for(cl in samp){ yb<-U[cl,]; if(var(yb)<=0) next
    f<-atpk.apply.series(yb,D,ww); e<-env_cell(cl); if(e<1e-12) next
    bw<-c(bw, median(2*1.96*sqrt(f$var))/e) }
  cat(sprintf("   %7.2f   |   %.3f                  | %.3f\n", r, median(bw), quantile(bw,.9)))
}
cat("\n=> the 'dominant' uncertainty magnitude is a free parameter: a shorter range\n")
cat("   inflates it, a longer range collapses it toward the spline. There is no\n")
cat("   data-driven way to set it (monthly autocorrelation = the seasonal cycle).\n")
## confirm near data-independence: unit-variance reduction is geometry-only;
## var = sill * geometry. Show two cells with very different sill have the SAME
## band/sqrt(sill) shape.
cl1<-samp[1]; cl2<-samp[which.max(sapply(samp,function(c)var(U[c,])))]
ww<-atpk.window.weights(M,W=6,Ns=6,range=1.5); f1<-atpk.apply.series(U[cl1,],D,ww); f2<-atpk.apply.series(U[cl2,],D,ww)
cat(sprintf("\nvariance is sill x geometry: cell A sill=%.2e cell B sill=%.2e ; sqrt(var)/sqrt(sill) median A=%.3f B=%.3f (≈equal => geometry-driven)\n",
    var(U[cl1,]),var(U[cl2,]), median(sqrt(f1$var)/sqrt(var(U[cl1,]))), median(sqrt(f2$var)/sqrt(var(U[cl2,])))))
