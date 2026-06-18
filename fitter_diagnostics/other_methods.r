#!/usr/bin/env Rscript
## Investigate the remaining candidate fitters against the same yardsticks:
##  - MSS  (in-repo QP smoothing spline): overshoot, sign-flip, NRT-locality, cost
##  - unlimited parabolic (PPM edges, NO limiter): shows what PPM's limiter prevents
##  - Steffen-on-cumulative (alt monotone-cubic to Fritsch-Carlson): overshoot
suppressWarnings(suppressMessages(library(quadprog)))
source("lib/mss_fit.r"); source("lib/ppm_fit.r")
P<-new.env(); load("fit.ppm.rda",envir=P)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=length(P$piqsfit.time)+1)); D<-diff(x)
N<-360*180; M<-length(D); Dm<-matrix(D,N,M,byrow=TRUE)
a<-matrix(P$piqsfit.gpp$a,N,M); b<-matrix(P$piqsfit.gpp$b,N,M); c0<-matrix(P$piqsfit.gpp$c,N,M)
U<-a*Dm^2/3+b*Dm/2+c0
um<-abs(U); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M]))
land<-which(rowSums(um)>1e-15 & env[,1]>-1)
peak_q<-function(aa,bb,cc,DD){fL<-cc;fR<-aa*DD^2+bb*DD+cc;sv<-ifelse(aa!=0,-bb/(2*aa),-1);intr<-aa!=0&sv>0&sv<DD;fV<-ifelse(intr,cc-bb^2/(4*aa),0);pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))}

## ---------- (1) MSS on a sample ----------
set.seed(11); samp<-sample(land,400)
setup<-mss.fit.setup(x)
t0<-Sys.time(); osMSS<-c(); flipMSS<-c(); foot<-integer(length(samp))
for(si in seq_along(samp)){ s<-samp[si]; yb<-U[s,]
  f<-mss.fit.cell(x,yb,setup)
  pk<-peak_q(f$a,f$b,f$c,D); e<-pmax(c(abs(yb)[1],abs(yb)[1:(M-1)]),abs(yb),c(abs(yb)[2:M],abs(yb)[M]))
  osMSS<-c(osMSS,max((pk/pmax(e,1e-30))[e>1e-12]))
  ## sign-flip: GPP should be <=0; flux>0 anywhere?
  fmin<-pmin(f$c, f$a*D^2+f$b*D+f$c); flipMSS<-c(flipMSS, as.integer(max(f$c, f$a*D^2+f$b*D+f$c) > 1e-12*max(abs(yb))))
  yb2<-yb; yb2[M]<-yb2[M]*1.10; f2<-mss.fit.cell(x,yb2,setup)
  dc<-abs(f2$c-f$c)/max(abs(yb),1e-30); chg<-which(dc>0.01); foot[si]<-if(length(chg))M-min(chg) else 0
}
tMSS<-as.numeric(Sys.time()-t0,units="secs")/length(samp)*1000
cat(sprintf("MSS (n=400): overshoot peak/env med %.3f max %.3f ; %% cells with flux sign-flip %.1f%% ; NRT footprint med %d max %d ; %.1f ms/cell\n",
            median(osMSS),max(osMSS),100*mean(flipMSS),median(foot),max(foot),tMSS))

## ---------- (2) unlimited parabolic (PPM edges, NO limiter) ----------
mm3<-function(a,b,cc){s<-sign(a);ifelse(s==sign(b)&s==sign(cc),s*pmin(abs(a),abs(b),abs(cc)),0)}
uL<-cbind(U[,1],U[,1:(M-1)]); uR<-cbind(U[,2:M],U[,M]); dc<-(uR-uL)/2; dm<-mm3(dc,2*(U-uL),2*(uR-U))
dmp1<-cbind(dm[,2:M],dm[,M]); aedge<-U+0.5*(uR-U)-(dmp1-dm)/6
aLp<-cbind(U[,1],aedge[,1:(M-1)]); aRp<-aedge; aLp[,1]<-U[,1];aRp[,1]<-U[,1];aLp[,M]<-U[,M];aRp[,M]<-U[,M]
a6<-6*(U-0.5*(aLp+aRp)); au<- -a6/Dm^2; bu<-(aRp-aLp+a6)/Dm; cu<-aLp   # NO ext/hi/lo limiting
pk<-peak_q(au,bu,cu,Dm); r<-(pk/env)[land,]; r<-r[is.finite(r)]
cat(sprintf("Unlimited parabolic: overshoot peak/env med %.3f 99th %.3f max %.3f  (vs PPM 1.00) -> what the limiter prevents\n",
            median(r),quantile(r,.99),max(r)))

## ---------- (3) Steffen-on-cumulative (monotone cubic, alt slope rule) ----------
## cumulative secant over segment k = U_k. per-row sign flip so secants positive.
sgn<-ifelse(rowMeans(U)<0,-1,1); Us<-U*sgn          # positive means
h<-Dm
## knot slopes d_k (k=1..M+1). interior: Steffen 1990
sL<-Us; sR<-cbind(Us[,2:M],Us[,M])                  # s_{k-1}=Us[,k], s_k=Us[,k+1] per knot k(2..M)
# build knot arrays length M+1
dk<-matrix(0,N,M+1)
for(k in 2:M){ sm1<-Us[,k-1]; sk<-Us[,k]; hk1<-h[,k-1]; hk<-h[,k]
  p<-(sm1*hk+sk*hk1)/(hk1+hk)
  dk[,k]<-(sign(sm1)+sign(sk))*pmin(abs(sm1),abs(sk),0.5*abs(p)) }
dk[,1]<-Us[,1]; dk[,M+1]<-Us[,M]                    # one-sided ends
Q<- -6*Us+3*dk[,1:M]+3*dk[,2:(M+1)]; L<-6*Us-4*dk[,1:M]-2*dk[,2:(M+1)]; K<-dk[,1:M]
aS<-sgn*Q/h^2; bS<-sgn*L/h; cS<-sgn*K
pk<-peak_q(aS,bS,cS,Dm); r<-(pk/env)[land,]; r<-r[is.finite(r)]
cat(sprintf("Steffen-on-cumulative: overshoot peak/env med %.3f 99th %.3f max %.3f  (vs PCHIP 1.50)\n",
            median(r),quantile(r,.99),max(r)))
