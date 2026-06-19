#!/usr/bin/env Rscript
## Unit tests for lib/atpk_fit.r :: area-to-point kriging fitter.
## Pins the contract: exact coherence (mass), non-negative variance, sign-safe
## reconstruction for a one-signed quantity, dormant-cell handling, and the
## small-magnitude scaling that the unit-covariance solve fixes.
## Run:  Rscript tests/test_atpk_fit.r   (exits non-zero on failure)
.args<-commandArgs(FALSE); .fa<-grep("^--file=",.args,value=TRUE)
.dir<-if(length(.fa)) dirname(sub("^--file=","",.fa[1])) else "."
.repo<-normalizePath(file.path(.dir,"..")); source(file.path(.repo,"lib","atpk_fit.r"))
.fail<-0L
check<-function(name,ok){cat(sprintf("  %s  %s\n",if(isTRUE(ok))"PASS" else "FAIL",name)); if(!isTRUE(ok)).fail<<-.fail+1L}

x <- as.numeric(seq(as.POSIXct("2018-01-01",tz="UTC"), by="1 month", length.out=37)); h<-diff(x); n<-36
seas <- -(pmax(sin((1:n)/12*2*pi - 1.2),0))^1.5 * 5e-7        # near-zero-flanked GPP peak
massrel <- function(f,yb) max(abs((f$a*h^3/3 + f$b*h^2/2 + f$c*h)/h - yb))/max(abs(yb),1e-300)
wrongsign <- function(f,yb){ sgn<-sign(mean(yb)); w<-0
  for(i in 1:n){ ss<-((1:25)-0.5)/25*h[i]; fl<-f$a[i]*ss^2+f$b[i]*ss+f$c[i]; if(any(sgn*fl < -1e-9*max(abs(yb)))) w<-w+1 }; w }

f <- atpk.fit.cell(x, seas)
check("returns a,b,c,var,dormant", all(c("a","b","c","var","dormant")%in%names(f)))
check("lengths == n", length(f$a)==n && length(f$var)==n)
check("not dormant for a real seasonal series", !f$dormant)

## (1) COHERENCE: per-piece integral reproduces the monthly mean exactly
check("mass-preserving (GPP seasonal) to FP", massrel(f,seas) < 1e-10)
yb2 <- abs(seas)+3e-8; f2<-atpk.fit.cell(x,yb2)
check("mass-preserving (positive RESP) to FP", massrel(f2,yb2) < 1e-10)
x3 <- as.numeric(seq(as.POSIXct("2018-01-01",tz="UTC"),by="1 month",length.out=37))  # real (non-uniform) month widths
check("mass-preserving holds on non-uniform months", massrel(atpk.fit.cell(x3,seas),seas) < 1e-10)

## (2) VARIANCE: non-negative, present, scales with data magnitude
check("variance non-negative", all(f$var>=0))
fa<-atpk.fit.cell(x, seas); fb<-atpk.fit.cell(x, seas*10)
check("variance scales ~100x with 10x data", abs(median(fb$var[fb$var>0])/median(fa$var[fa$var>0]) - 100) < 1)

## (3) SIGN-SAFE reconstruction (one-signed quantity), incl. the hard cases
check("no wrong-sign: GPP seasonal", wrongsign(f,seas)==0L)
check("no wrong-sign: alternating", wrongsign(atpk.fit.cell(x,rep(c(-1e-7,0),18)),rep(c(-1e-7,0),18))==0L)
check("no wrong-sign: positive RESP", wrongsign(f2,yb2)==0L)

## (4) DORMANT-cell handling
fz<-atpk.fit.cell(x, rep(0,n))
check("all-zero -> dormant, flat, zero var", fz$dormant && all(fz$a==0)&&all(fz$b==0)&&all(fz$var==0))

## (5) SMALL-MAGNITUDE scaling (the unit-covariance solve): tiny but nonzero
ft<-atpk.fit.cell(x, seas*1e-6)
check("tiny-magnitude series still fits (not spuriously dormant)", !ft$dormant && massrel(ft,seas*1e-6)<1e-10)

## (6) points->abc is itself mass-preserving (unit check)
z<-rep(seas,each=6)+rnorm(n*6,0,1e-9); ss<-((1:6)-0.5)/6
abc<-atpk.points.to.abc(z,ss,6,n,h,seas,sgn=-1)
check("points.to.abc preserves each piece mean", max(abs((abc$a*h^3/3+abc$b*h^2/2+abc$c*h)/h - seas))/max(abs(seas)) < 1e-10)

if(.fail>0L){cat(sprintf("\n%d FAILED\n",.fail)); quit(status=1L)}
cat("\nall atpk_fit tests passed\n")
