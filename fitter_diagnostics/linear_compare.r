#!/usr/bin/env Rscript
## Integral-preserving LINEAR flux reconstructions vs PCHIP, on the real
## monthly means (recovered exactly from the shipped PCHIP fit).
##
##   (B) CONTINUOUS at knots -> forced trapezoidal recursion y_{i+1}=2u_i-y_i.
##       Pole at z=-1 resonates with month-to-month alternation -> unstable.
##   (C) DISCONTINUOUS slope-limited (MUSCL/minmod): per-cell line through the
##       mean, slope = minmod(dL,dR). Integral-preserving, provably NO
##       overshoot, keeps a within-month gradient; jump at each month edge.
suppressWarnings(load("fit.piqs.rda"))
t0 <- as.POSIXct("2001-01-01", tz="UTC")
edges <- as.numeric(seq(t0, by="1 month", length.out=length(piqsfit.time)+1))
D <- diff(edges)
minmod <- function(a,b) 0.5*(sign(a)+sign(b))*pmin(abs(a),abs(b))

cmp <- function(fit, label) {
  N<-360*180; M<-length(D)
  a<-matrix(fit$a,N,M); b<-matrix(fit$b,N,M); c0<-matrix(fit$c,N,M)
  Dm<-matrix(D,N,M,byrow=TRUE)
  u <- a*Dm^2/3 + b*Dm/2 + c0                   # monthly means (exact)
  um<-abs(u)
  env<-pmax(cbind(um[,1],um[,1:(M-1)]), um, cbind(um[,2:M],um[,M]))
  keep <- (rowSums(um) > 1e-15) & (env > 1e-12)

  ## PCHIP within-piece peak
  fL<-c0; fR<-a*Dm^2+b*Dm+c0
  sv<-ifelse(a!=0,-b/(2*a),-1); intr<-a!=0&sv>0&sv<Dm; fV<-ifelse(intr,c0-b^2/(4*a),0)
  pk.pchip<-pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))

  ## (C) minmod-limited integral-preserving linear
  dL <- u - cbind(u[,1], u[,1:(M-1)])
  dR <- cbind(u[,2:M], u[,M]) - u
  half <- 0.5*minmod(dL,dR)
  Redge<-u+half; Ledge<-u-half
  pk.mm<-pmax(abs(Ledge),abs(Redge))
  jump <- abs(Redge[,1:(M-1)] - Ledge[,2:M])    # edge discontinuity, [N,M-1]
  grad <- abs(half)                             # within-month gradient (qmod amp)

  ## (B) continuous trapezoidal
  yk<-matrix(0,N,M+1); yk[,1]<-u[,1]
  for(i in 1:M) yk[,i+1]<-2*u[,i]-yk[,i]
  pk.tr<-pmax(abs(yk[,1:M]),abs(yk[,2:(M+1)]))

  r.pchip<-(pk.pchip/env)[keep]; r.mm<-(pk.mm/env)[keep]; r.tr<-(pk.tr/env)[keep]
  jvals<-jump/env[,1:(M-1)]; jland<-jvals[keep]; jland<-jland[is.finite(jland)]
  gland<-(grad/env)[keep]

  cat(sprintf("\n==== %s : %d land cell-pieces ====\n", label, sum(keep)))
  cat(sprintf("PCHIP (current)        peak/env: med %.2f  99th %.2f  max %.2f\n",
              quantile(r.pchip,.5),quantile(r.pchip,.99),max(r.pchip)))
  cat(sprintf("minmod linear (C)      peak/env: med %.2f  99th %.2f  max %.2f   [no overshoot]\n",
              quantile(r.mm,.5),quantile(r.mm,.99),max(r.mm)))
  cat(sprintf("trapezoidal linear (B) peak/env: 99th %.3g  max %.3g   [unstable]\n",
              quantile(r.tr,.99),max(r.tr)))
  cat(sprintf("minmod month-edge JUMP / env: med %.2f  90th %.2f  99th %.2f  (PCHIP=0, continuous)\n",
              quantile(jland,.5),quantile(jland,.9),quantile(jland,.99)))
  cat(sprintf("minmod within-month gradient / env: med %.2f  (piecewise-constant=0)\n",
              quantile(gland,.5)))
}
cmp(piqsfit.gpp,"GPP (-2*NPP)")
cmp(piqsfit.resp,"RESP (Rh+NPP)")
