## Investigate the bounded iterative mean-preserving family (Rymes & Myers 2001):
## iterated 3-point moving average + per-interval mean-restoration + bound clip.
## Decisive question: NRT-locality (iteration couples neighbours -> how far does a
## last-month revision propagate?), plus overshoot / sign-flip / mass.
P<-new.env(); load("fit.ppm.rda",envir=P)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(x)
N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
a<-matrix(P$piqsfit.gpp$a,N,M); b<-matrix(P$piqsfit.gpp$b,N,M); c0<-matrix(P$piqsfit.gpp$c,N,M)
U<-a*Dm^2/3+b*Dm/2+c0
land<-which(rowSums(abs(U))>1e-15)
set.seed(9); samp<-sample(land,200)

## Rymes-Myers MPA on a 1-cell monthly-mean series; nsub points/month; GPP<=0 clip.
mpa<-function(means,nsub=30,niter=100,upper=0){
  Mn<-length(means); y<-rep(means,each=nsub); Nn<-length(y)
  for(it in 1:niter){
    ys<-y; ys[2:(Nn-1)]<-(y[1:(Nn-2)]+y[2:(Nn-1)]+y[3:Nn])/3
    ys[1]<-(y[1]+y[2])/2; ys[Nn]<-(y[Nn-1]+y[Nn])/2; y<-ys
    for(m in 1:Mn){idx<-((m-1)*nsub+1):(m*nsub); y[idx]<-y[idx]+(means[m]-mean(y[idx]))}
    if(!is.na(upper)) y[y>upper]<-upper
  }
  matrix(y,nrow=nsub)        # [nsub, Mn]
}
foot<-function(yb,niter){
  f0<-mpa(yb,niter=niter); yb2<-yb; yb2[M]<-yb2[M]*1.10; f1<-mpa(yb2,niter=niter)
  permon<-apply(abs(f1-f0),2,max)/max(abs(yb),1e-30)   # per-month max rel change
  chg<-which(permon>0.01); if(length(chg)) M-min(chg) else 0
}
for(niter in c(30,100,300)){
  ov<-flip<-mass<-fp<-numeric(length(samp))
  for(si in seq_along(samp)){ yb<-U[samp[si],]; fc<-mpa(yb,niter=niter)
    mn<-colMeans(fc); mass[si]<-max(abs(mn-yb)/pmax(abs(yb),1e-30))
    e<-pmax(c(abs(yb)[1],abs(yb)[1:(M-1)]),abs(yb),c(abs(yb)[2:M],abs(yb)[M]))
    pk<-apply(abs(fc),2,max); ov[si]<-max((pk/pmax(e,1e-30))[e>1e-12])
    flip[si]<-as.integer(max(fc)>1e-12*max(abs(yb)))
    fp[si]<-foot(yb,niter) }
  cat(sprintf("Rymes-Myers niter=%3d: mass-err med %.1e ; overshoot peak/env med %.2f max %.2f ; %%sign-flip %.0f%% ; NRT footprint months med %d max %d\n",
      niter, median(mass), median(ov), max(ov), 100*mean(flip), median(fp), max(fp)))
}
cat("  (compare NRT footprint: PCHIP 0, PPM <=2, minmod <=1, PIQS all 302)\n")
