## Close the quantitative gaps in the case against "PIQS + linear-fallback-on-
## overshoot" (the stakeholder-preferred hybrid). Extends piqs_hybrid.r with:
##   (A) FULL distribution of the patch discontinuity (not just the median),
##   (B) continuous integral-preserving LINEAR knot-sign-flip count (the
##       y_{i+1}=2*m_i - y_i recursion that PROPOSALS #9 warns amplifies
##       alternation onto the knots),
##   (C) a discontinuity-vs-overshoot "error budget" (L1 norm).
## NRT footprint is measured separately by a perturbed PIQS refit (see the
## launching job). Reads fit.piqs_v1.rda + fit.pchip.rda. Pure base R + ncdf4.
suppressWarnings(suppressMessages(library(ncdf4)))
Q<-new.env(); load("fit.piqs_v1.rda",envir=Q); P<-new.env(); load("fit.pchip.rda",envir=P)
N<-360*180; M<-length(Q$piqsfit.time)
x<-as.numeric(seq(as.POSIXct("2001-01-01",tz="UTC"),by="1 month",length.out=M+1)); D<-diff(x); Dm<-matrix(D,N,M,byrow=TRUE)
minmod<-function(a,b)0.5*(sign(a)+sign(b))*pmin(abs(a),abs(b))
mat<-function(E,c,k) matrix(E[[c]][[k]],N,M); comp<-"piqsfit.gpp"
qa<-mat(Q,comp,"a");qb<-mat(Q,comp,"b");qc<-mat(Q,comp,"c"); pa<-mat(P,comp,"a");pb<-mat(P,comp,"b");pc<-mat(P,comp,"c")
u<-pa*Dm^2/3+pb*Dm/2+pc; um<-abs(u); env<-pmax(cbind(um[,1],um[,1:(M-1)]),um,cbind(um[,2:M],um[,M])); keep<-rowSums(um)>1e-15
qend<-qa*Dm^2+qb*Dm+qc; sv<-ifelse(qa!=0,-qb/(2*qa),-1); intr<-qa!=0&sv>0&sv<Dm; qvtx<-ifelse(intr,qc-qb^2/(4*qa),0)
overshoot<-(qc>1e-9*env)|(qend>1e-9*env)|(intr&qvtx>1e-9*env)
dL<-u-cbind(u[,1],u[,1:(M-1)]); dR<-cbind(u[,2:M],u[,M])-u; slope<-minmod(dL,dR)/Dm; hc<-u-slope*Dm/2
Ha<-qa;Hb<-qb;Hc<-qc; Ha[overshoot]<-0;Hb[overshoot]<-slope[overshoot];Hc[overshoot]<-hc[overshoot]

cat(sprintf("PIQS overshoot / fallback-trigger rate (land cell-months): %.1f%%\n",100*mean(overshoot[keep,])))

## (A) FULL discontinuity distribution at patched edges, value-jump/env
qval<-qend[,1:(M-1)]; hval<-Hc[,2:M]; patched<-overshoot[,2:M] & keep[1:N]
ajump<-abs(qval-hval)[patched]                # absolute jump (un-normalised)
vj<-abs(qval-hval)/env[,1:(M-1)]; vjp<-vj[patched]
fin<-is.finite(vjp); ninf<-sum(!fin)
qs<-quantile(vjp[fin],c(.1,.25,.5,.75,.9,.99),na.rm=T)
cat(sprintf("(A) patch discontinuity / env (finite, n=%d)  p10=%.2f p25=%.2f MED=%.2f p75=%.2f p90=%.2f p99=%.2f\n",
    sum(fin),qs[1],qs[2],qs[3],qs[4],qs[5],qs[6]))
cat(sprintf("    %.1f%% of patched edges exceed 1x env ; %.1f%% exceed 3x env ; %.1f%% have env~0 so the jump EXCEEDS the entire local monthly flux\n",
    100*mean(vjp>1,na.rm=T),100*mean(vjp>3,na.rm=T),100*ninf/length(vjp)))
write.csv(data.frame(jump_over_env=round(pmin(vjp[fin],50),4)),
          "fitter_diagnostics/linear_fallback_discontinuity.csv",row.names=FALSE)

## (B) continuous integral-preserving LINEAR: y_{i+1}=2*m_i - y_i, seed y_1=m_1.
## Count cells/knots where the knot series flips sign vs the (single-signed) data.
sgn<-sign(rowMeans(u)); y<-matrix(0,N,M+1); y[,1]<-u[,1]
for(i in 1:M) y[,i+1]<-2*u[,i]-y[,i]
land<-keep
flipknot<- (sgn*y[,2:M] < -1e-9*env[,1:(M-1)])      # interior knots wrong sign
cells.any<-rowSums(flipknot[land,,drop=FALSE]>0)>0
cat(sprintf("(B) continuous-linear knot recursion: %.1f%% of land cells get >=1 wrong-sign knot; %.2f%% of interior knots flip (vs PIQS overshoot %.1f%%)\n",
    100*mean(cells.any), 100*mean(flipknot[land,]), 100*mean(overshoot[keep,])))
ry<-(abs(y[,2:M])/env[,1:(M-1)])[land,]; ry<-ry[is.finite(ry)]
cat(sprintf("    knot amplification |knot|/env (continuous linear): p99=%.1f  p99.9=%.1f  (unbounded resonance, PROPOSALS #9 Nyquist pole; PIQS quadratic absorbs this into curvature)\n",
    quantile(ry,.99), quantile(ry,.999)))

## (C) error budget: L1 norm of hybrid discontinuity field vs PCHIP (0, C0).
##     Sum of |value-jump| over all patched edges, area-agnostic, /env-normalised.
cat(sprintf("(C) total ABSOLUTE discontinuity budget (sum |jump| over %d patched land edges, mol m-2 s-1): hybrid=%.3e  PCHIP=0 (C0 by construction)\n",
    length(ajump), sum(ajump,na.rm=T)))
