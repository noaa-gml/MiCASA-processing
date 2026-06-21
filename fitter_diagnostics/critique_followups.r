## M2/M3/M4 follow-ups to the adversarial review.
source("lib/pchip_fit.r"); source("lib/ppm_fit.r"); source("lib/linmm_fit.r")
N<-360*180
A<-new.env(); load("fit.atpk.rda",envir=A); P<-new.env(); load("fit.pchip.rda",envir=P); Q<-new.env(); load("fit.piqs_v1.rda",envir=Q)
M<-length(A$piqsfit.time); x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=M+1)); D<-diff(x); Dm<-matrix(D,N,M,byrow=TRUE)
latof<-function(cl) -90+((cl-1)%/%360)+1-0.5; bandf<-function(l) ifelse(abs(l)<=23.5,"tropics",ifelse(abs(l)<=50,"temperate","boreal/polar"))
pa<-matrix(P$piqsfit.gpp$a,N,M); pb<-matrix(P$piqsfit.gpp$b,N,M); pc<-matrix(P$piqsfit.gpp$c,N,M); U<-pa*Dm^2/3+pb*Dm/2+pc
land<-rowSums(abs(U))>1e-15; lat<-latof(1:N); bnd<-bandf(lat)

## --- M2: ATP flat-fallback rate per biome ---
## fallback sets a==0 & b==0 in a varying cell with a nonzero mean.
aa<-matrix(A$piqsfit.gpp$a,N,M); ab<-matrix(A$piqsfit.gpp$b,N,M)
flat <- (aa==0 & ab==0) & (abs(U)>1e-15)            # piece is flat despite nonzero mean
cat("== M2: ATP flat-fallback rate (flat piece in a non-zero-mean land cell) ==\n")
for(bb in c("tropics","temperate","boreal/polar")){ m<-land & bnd==bb
  cat(sprintf("  %-13s %.1f%% of land cell-months flattened\n", bb, 100*mean(flat[m,]))) }
cat(sprintf("  overall %.1f%% ; (in flattened pieces the point estimate is flat; the $var there is the kriging spread, not the spread of a structured curve)\n", 100*mean(flat[land,])))

## --- M3: dropped-cell (non-finite ratio) rate for the overshoot metric ---
qa<-matrix(Q$piqsfit.gpp$a,N,M); qb<-matrix(Q$piqsfit.gpp$b,N,M); qc<-matrix(Q$piqsfit.gpp$c,N,M)
um<-abs(U); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M]))
peak<-pmax(abs(qc),abs(qa*Dm^2+qb*Dm+qc)); ratio<-peak/env
dropped<- land & !is.finite(ratio[,1]) | (env<1e-12)   # would be dropped
cat(sprintf("\n== M3: envelope-normalization drop rate ==\n  %% land cell-months with env<1e-12 (dropped from overshoot medians): %.2f%%\n",
    100*mean((env<1e-12)[land,])))
## floored-envelope version: env_floor = max(env, global median flux)
gf<-median(env[land,][env[land,]>0]); envf<-pmax(env,gf)
cat(sprintf("  PIQS overshoot peak/env: median(strict, drops near-zero)=%.2f  median(floored env)=%.2f  (floored includes the near-zero cells)\n",
    median((peak/env)[land,][is.finite((peak/env)[land,])]), median((peak/envf)[land,])))

## --- M4: NRT footprint checking the FULL flux (a,b,c), not just c ---
set.seed(3); samp<-sample(which(land & abs(U[,M])>1e-12), 1500); xe<-x; ss<-((1:6)-0.5)/6
foot<-function(fitfun){ fp<-integer(length(samp))
  for(s in seq_along(samp)){ yb<-U[samp[s],]; f0<-fitfun(xe,yb); yb2<-yb; yb2[M]<-yb2[M]*1.10; f1<-fitfun(xe,yb2)
    e<-max(abs(yb)); chg<-rep(FALSE,M)
    for(q in ss){ dt<-q*D; fl0<-f0$a*dt^2+f0$b*dt+f0$c; fl1<-f1$a*dt^2+f1$b*dt+f1$c; chg<-chg | (abs(fl1-fl0)/e>0.01) }
    w<-which(chg); fp[s]<-if(length(w)) M-min(w) else 0 }
  fp }
cat("\n== M4: NRT footprint via full-flux (a,b,c) check, +10% last-month perturb ==\n")
for(nm in c("PCHIP","PPM","minmod")){ ff<-switch(nm,PCHIP=pchip.fit.cell,PPM=ppm.fit.cell,minmod=linmm.fit.cell)
  fp<-foot(ff); cat(sprintf("  %-8s prior months changed >1%%: median %d  90th %d  max %d\n", nm, median(fp), quantile(fp,.9), max(fp))) }
cat("  (PIQS global tridiagonal solve couples all 302 months by construction; not refit here)\n")
