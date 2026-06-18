#!/usr/bin/env Rscript
## Unit tests for lib/ppm_fit.r :: ppm.fit.cell + ppm.fit.grid.
## PPM (Colella & Woodward 1984) limited piecewise-parabolic integral-
## preserving fitter that write_ppm.r runs over the grid. Pins the contract +
## the vectorized grid==cell equivalence (write_ppm.r relies on the grid path).
## Run:  Rscript tests/test_ppm_fit.r   (exits non-zero on failure)
.args<-commandArgs(FALSE); .fa<-grep("^--file=",.args,value=TRUE)
.dir<-if(length(.fa)) dirname(sub("^--file=","",.fa[1])) else "."
.repo<-normalizePath(file.path(.dir,"..")); source(file.path(.repo,"lib","ppm_fit.r"))
.fail<-0L
check<-function(name,ok){cat(sprintf("  %s  %s\n",if(isTRUE(ok))"PASS" else "FAIL",name)); if(!isTRUE(ok)).fail<<-.fail+1L}
peval<-function(fit,k,s) fit$a[k]*s^2+fit$b[k]*s+fit$c[k]
pint <-function(fit,k,h) fit$a[k]*h^3/3+fit$b[k]*h^2/2+fit$c[k]*h

fit<-ppm.fit.cell(0:5,c(1,4,2,8,3))
check("returns list a,b,c",all(c("a","b","c")%in%names(fit)))
check("length n",length(fit$a)==5L&&length(fit$b)==5L&&length(fit$c)==5L)

z<-ppm.fit.cell(0:3,c(0,0,0)); check("all-zero -> zero coeffs",all(z$a==0)&&all(z$b==0)&&all(z$c==0))
cf<-ppm.fit.cell(0:4,rep(5,4)); check("constant -> a==0",max(abs(cf$a))<1e-12)
check("constant -> b==0",max(abs(cf$b))<1e-12); check("constant -> c==value",max(abs(cf$c-5))<1e-12)

## integral preservation (uniform + non-uniform) -- the core invariant
yb<-c(1,4,2,8,3); fit<-ppm.fit.cell(0:5,yb)
check("uniform knots: piece integrates to ybar[k]",
      max(sapply(seq_along(yb),function(k)abs(pint(fit,k,1)-yb[k])))<1e-12)
x2<-c(0,1,3,4,7); yb2<-c(2,5,1,6); fit2<-ppm.fit.cell(x2,yb2); h2<-diff(x2)
check("non-uniform knots: piece integrates to ybar[k]*h",
      max(sapply(seq_along(yb2),function(k)abs(pint(fit2,k,h2[k])-yb2[k]*h2[k])))<1e-12)

## NO OVERSHOOT including the alternating + spike cases that break naive fits
for(nm in c("sharp","alt","spike")){
  yb<-switch(nm, sharp=c(0.1,5,0.2,8,0.1), alt=c(1,0,2,0,3,0,2,0), spike=c(0,0,9,0,0))
  f<-ppm.fit.cell(seq(0,length(yb)),yb); n<-length(yb)
  env<-pmax(c(yb[1],yb[1:(n-1)]),yb,c(yb[2:n],yb[n]))
  pk<-sapply(1:n,function(k)max(abs(peval(f,k,seq(0,1,length.out=40)))))
  check(sprintf("no overshoot: %s peak/env <= 1",nm), max((pk/pmax(env,1e-30))[env>0])<=1+1e-12)
}
## negative (GPP) input
neg<-ppm.fit.cell(0:5,-c(0.1,5,0.2,8,0.1))
check("no overshoot (negative input)",
      max(sapply(1:5,function(k)max(abs(peval(neg,k,seq(0,1,length.out=40))))))<=max(c(0.1,5,0.2,8,0.1))+1e-12)

## GRID == CELL to FP (the vectorization write_ppm.r depends on)
set.seed(13); x<-0:12; U<-matrix(runif(7*12,-3,5),7,12); U[3,]<-0
g<-ppm.fit.grid(x,U); mx<-0; sc<-max(abs(U))
for(r in 1:7){ cc<-ppm.fit.cell(x,U[r,])
  mx<-max(mx,abs(g$a[r,]-cc$a)*1,abs(g$b[r,]-cc$b),abs(g$c[r,]-cc$c)) }
check("grid == cell (FP tolerance)", mx < 1e-12*max(sc,1))

if(.fail>0L){cat(sprintf("\n%d FAILED\n",.fail)); quit(status=1L)}
cat("\nall ppm_fit tests passed\n")
