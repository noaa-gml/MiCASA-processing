#!/usr/bin/env Rscript
## Unit tests for lib/linmm_fit.r :: linmm.fit.cell + linmm.fit.grid.
## minmod-limited integral-preserving piecewise-linear (MUSCL) fitter that
## write_linmm.r runs over the grid. Pins the contract + the vectorized
## grid==cell equivalence (write_linmm.r relies on the grid path).
## Run:  Rscript tests/test_linmm_fit.r   (exits non-zero on failure)
.args<-commandArgs(FALSE); .fa<-grep("^--file=",.args,value=TRUE)
.dir<-if(length(.fa)) dirname(sub("^--file=","",.fa[1])) else "."
.repo<-normalizePath(file.path(.dir,"..")); source(file.path(.repo,"lib","linmm_fit.r"))
.fail<-0L
check<-function(name,ok){cat(sprintf("  %s  %s\n",if(isTRUE(ok))"PASS" else "FAIL",name)); if(!isTRUE(ok)).fail<<-.fail+1L}
peval<-function(fit,k,s) fit$a[k]*s^2+fit$b[k]*s+fit$c[k]
pint <-function(fit,k,h) fit$a[k]*h^3/3+fit$b[k]*h^2/2+fit$c[k]*h

fit<-linmm.fit.cell(0:5,c(1,4,2,8,3))
check("returns list a,b,c",all(c("a","b","c")%in%names(fit)))
check("length n",length(fit$a)==5L&&length(fit$b)==5L&&length(fit$c)==5L)
check("linear: a==0 everywhere",all(fit$a==0))

z<-linmm.fit.cell(0:3,c(0,0,0)); check("all-zero -> zero coeffs",all(z$a==0)&&all(z$b==0)&&all(z$c==0))
cf<-linmm.fit.cell(0:4,rep(5,4)); check("constant -> slope 0",max(abs(cf$b))<1e-12)
check("constant -> c==value",max(abs(cf$c-5))<1e-12)

## integral preservation (uniform + non-uniform)
yb<-c(1,4,2,8,3); fit<-linmm.fit.cell(0:5,yb)
check("uniform knots: piece integrates to ybar[k]",
      max(sapply(seq_along(yb),function(k)abs(pint(fit,k,1)-yb[k])))<1e-12)
x2<-c(0,1,3,4,7); yb2<-c(2,5,1,6); fit2<-linmm.fit.cell(x2,yb2); h2<-diff(x2)
check("non-uniform knots: piece integrates to ybar[k]*h",
      max(sapply(seq_along(yb2),function(k)abs(pint(fit2,k,h2[k])-yb2[k]*h2[k])))<1e-12)

## NO OVERSHOOT: peak within local monthly-mean envelope (the whole point)
sharp<-c(0.1,5,0.2,8,0.1); sf<-linmm.fit.cell(0:5,sharp)
env<-pmax(c(sharp[1],sharp[1:4]),sharp,c(sharp[2:5],sharp[5]))
pk<-sapply(1:5,function(k)max(abs(peval(sf,k,seq(0,1,length.out=40)))))
check("no overshoot: peak/envelope <= 1",max(pk/env)<=1+1e-12)
## works for negative (GPP) input too
neg<-linmm.fit.cell(0:5,-sharp)
check("no overshoot (negative input)", max(sapply(1:5,function(k)max(abs(peval(neg,k,seq(0,1,length.out=40)))))/env)<=1+1e-12)

## GRID == CELL (the vectorization write_linmm.r depends on)
set.seed(11); x<-0:12; U<-matrix(runif(7*12,-3,5),7,12); U[3,]<-0
g<-linmm.fit.grid(x,U); mx<-0
for(r in 1:7){ cc<-linmm.fit.cell(x,U[r,])
  mx<-max(mx,abs(g$a[r,]-cc$a),abs(g$b[r,]-cc$b),abs(g$c[r,]-cc$c)) }
check("grid == cell (bit-identical)", mx==0)

if(.fail>0L){cat(sprintf("\n%d FAILED\n",.fail)); quit(status=1L)}
cat("\nall linmm_fit tests passed\n")
