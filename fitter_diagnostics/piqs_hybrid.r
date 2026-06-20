## Investigate PIQS-with-linear-fallback-on-overshoot: keep PIQS's global-solve
## quadratic where it is sign-safe, patch overshooting pieces with a sign-safe
## minmod-linear. Questions: (1) sign-safety achieved? (2) is PIQS measurably
## SMOOTHER than PCHIP (the only motivation)? (3) patch rate + discontinuity
## introduced? (4) daily fidelity vs PCHIP? (Still inherits PIQS global solve =>
## NRT non-local, footprint = whole record.)
suppressWarnings(suppressMessages(library(ncdf4)))
Q<-new.env(); load("fit.piqs_v1.rda",envir=Q)      # PIQS (global-solve quadratics)
P<-new.env(); load("fit.pchip.rda",envir=P)        # PCHIP
N<-360*180; M<-length(Q$piqsfit.time)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=M+1)); D<-diff(x); Dm<-matrix(D,N,M,byrow=TRUE)
minmod<-function(a,b)0.5*(sign(a)+sign(b))*pmin(abs(a),abs(b))

mat<-function(E,comp,k) matrix(E[[comp]][[k]],N,M)
comp<-"piqsfit.gpp"
qa<-mat(Q,comp,"a"); qb<-mat(Q,comp,"b"); qc<-mat(Q,comp,"c")
pa<-mat(P,comp,"a"); pb<-mat(P,comp,"b"); pc<-mat(P,comp,"c")
u <- pa*Dm^2/3+pb*Dm/2+pc                            # monthly means (PCHIP integral)
um<-abs(u); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M])); keep<-rowSums(um)>1e-15
sgn<- -1                                             # GPP <= 0

## (1) PIQS overshoot per piece: quadratic violates sign on [0,D]?
qend<-qa*Dm^2+qb*Dm+qc; sv<-ifelse(qa!=0,-qb/(2*qa),-1); intr<-qa!=0&sv>0&sv<Dm
qvtx<-ifelse(intr,qc-qb^2/(4*qa),0)
overshoot<- (qc>1e-9*env) | (qend>1e-9*env) | (intr & qvtx>1e-9*env)   # any wrong-sign point
cat(sprintf("PIQS overshoot pieces (land cell-months): %.1f%%\n", 100*mean(overshoot[keep,])))

## (2) smoothness: knot first-derivative jump |f'(end_i)-f'(start_{i+1})|*D/env
##     (smaller = smoother sub-monthly transitions). PIQS vs PCHIP.
djump<-function(a,b){ fpR<-2*a*Dm+b; fpL<-b; abs((fpR[,1:(M-1)]-fpL[,2:M])*Dm[,1:(M-1)])/env[,1:(M-1)] }
dq<-djump(qa,qb); dp<-djump(pa,pb)
kk<-keep[1:N]
cat(sprintf("knot deriv-jump*D/env (smoothness; lower=smoother):  PIQS med %.3f   PCHIP med %.3f\n",
    median(dq[kk,],na.rm=T), median(dp[kk,],na.rm=T)))
## curvature |2a|*D^2/env
cq<-(abs(2*qa*Dm^2)/env); cp<-(abs(2*pa*Dm^2)/env)
cat(sprintf("within-piece curvature/env (median):                  PIQS %.3f   PCHIP %.3f\n",
    median(cq[keep,]),median(cp[keep,])))

## (3) build hybrid: PIQS where safe, minmod-linear fallback where overshoot
dL<-u-cbind(u[,1],u[,1:(M-1)]); dR<-cbind(u[,2:M],u[,M])-u
slope<-minmod(dL,dR)/Dm; hc<-u-slope*Dm/2
Ha<-qa; Hb<-qb; Hc<-qc
Ha[overshoot]<-0; Hb[overshoot]<-slope[overshoot]; Hc[overshoot]<-hc[overshoot]
## hybrid sign-flips
hend<-Ha*Dm^2+Hb*Dm+Hc; hsv<-ifelse(Ha!=0,-Hb/(2*Ha),-1); hintr<-Ha!=0&hsv>0&hsv<Dm; hvtx<-ifelse(hintr,Hc-Hb^2/(4*Ha),0)
hwrong<-(Hc>1e-9*env)|(hend>1e-9*env)|(hintr&hvtx>1e-9*env)
cat(sprintf("hybrid wrong-sign pieces after fallback: %.3f%% (target ~0)\n",100*mean(hwrong[keep,])))
## (3b) discontinuity the patches introduce: PIQS is C0 (shares knot values); a
## patched piece breaks that. value-jump at patched-piece left edge / env.
qval_endprev<- qend[,1:(M-1)]            # PIQS right-edge of piece i
hval_start<- Hc[,2:M]                    # hybrid left value of piece i+1
patched_next<- overshoot[,2:M]
vj<-abs(qval_endprev - hval_start)/env[,1:(M-1)]
cat(sprintf("discontinuity at patched edges / env: median %.3f (only at the %.1f%% patched pieces)\n",
    median(vj[patched_next & keep[1:N]],na.rm=T), 100*mean(overshoot[keep,])))

## (4) daily fidelity vs MiCASA daily 2020: PCHIP vs PIQS vs hybrid
i2020<-(2020-2001)*12+(1:12); dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31)
ev<-function(a,b,c0,idx,dt) a[,idx]*dt^2+b[,idx]*dt+c0[,idx]
rms<-list(pchip=c(),piqs=c(),hybrid=c())
for(k in 1:12){ idx<-i2020[k]; nd<-dpm[k]; dts<-((1:nd)-0.5)*86400; mm<-sprintf("%02d",k)
  g<-matrix(0,N,nd)
  for(dd in 1:nd){ f<-nc_open(sprintf("daily_1x1/MiCASA_v1_flux_x360_y180_daily_2020%s%02d.nc",mm,dd)); g[,dd]<- -2*as.vector(ncvar_get(f,"NPP")); nc_close(f) }
  e<-env[,idx]; kp<-e>1e-12 & rowSums(abs(g))>1e-15
  rp<-sapply(dts,function(dt)ev(pa,pb,pc,idx,dt)); rq<-sapply(dts,function(dt)ev(qa,qb,qc,idx,dt)); rh<-sapply(dts,function(dt)ev(Ha,Hb,Hc,idx,dt))
  rmse<-function(X)sqrt(rowMeans((X-g)^2))
  rms$pchip<-c(rms$pchip,(rmse(rp)/e)[kp]); rms$piqs<-c(rms$piqs,(rmse(rq)/e)[kp]); rms$hybrid<-c(rms$hybrid,(rmse(rh)/e)[kp]) }
cat(sprintf("\ndaily RMSE/env 2020 (median / mean): PCHIP %.3f/%.3f  PIQS %.3f/%.3f  hybrid %.3f/%.3f\n",
    median(rms$pchip),mean(rms$pchip),median(rms$piqs),mean(rms$piqs),median(rms$hybrid),mean(rms$hybrid)))
