## Sub-grid spatial heterogeneity as a prior-uncertainty source: a 1deg cell is
## an area-average of 100 0.1deg pixels; the spread ACROSS those pixels measures
## how representative the 1deg value is of its actual landscape -- a data-driven
## uncertainty with NO model assumption (unlike ATP covariance / bootstrap model
## / across-fitter structural spread). Two flavours: (1) monthly-mean spread,
## (2) sub-monthly spread (PCHIP each sub-cell, spread of the 100 curves).
suppressWarnings(suppressMessages(library(ncdf4)))
source("lib/pchip_fit.r")
PC<-new.env(); load("fit.pchip.rda",envir=PC); a<-matrix(PC$piqsfit.gpp$a,360*180,length(PC$piqsfit.time))
land<-which(rowSums(abs(a))>0); latof<-function(cl) -90+((cl-1)%/%360)+1-0.5
band<-function(l) ifelse(abs(l)<=23.5,"tropics",ifelse(abs(l)<=50,"temperate","boreal/polar"))
ll<-latof(land); set.seed(8)
samp<-c(sample(land[abs(ll)<23.5],10),sample(land[abs(ll)>=23.5&abs(ll)<50],10),sample(land[abs(ll)>=50],10)); ns<-length(samp)
i1<-((samp-1)%%360)+1; j1<-((samp-1)%/%360)+1

yrs<-2018:2021; mons<-as.vector(sapply(yrs,function(y)sprintf("%04d%02d",y,1:12))); nM<-length(mons)
edges<-as.numeric(seq(as.POSIXct("2018-01-01",tz="UTC"),by="1 month",length.out=nM+1)); h<-diff(edges)
fpath<-function(mm)sprintf("portal.nccs.nasa.gov/monthly/%s/MiCASA_v1_flux_x3600_y1800_monthly_%s.nc4",substr(mm,1,4),mm)
GPP<-array(0,c(ns,100,nM))
for(m in 1:nM){ f<-nc_open(fpath(mons[m])); npp<-ncvar_get(f,"NPP"); nc_close(f)
  for(s in 1:ns){ li<-((i1[s]-1)*10+1):(i1[s]*10); lj<-((j1[s]-1)*10+1):(j1[s]*10); GPP[s,,m]<- -2*as.vector(npp[li,lj]) } }

intM<-6:(nM-6); Ns<-6; ss<-(seq_len(Ns)-0.5)/Ns
hetM<-subM<-bnd<-nact<-c()
for(s in 1:ns){ sub<-GPP[s,,]; bb<-band(latof(samp[s]))
  active<-which(rowSums(abs(sub))>1e-15); na<-length(active); if(na<3) next
  m1<-colMeans(sub[active,,drop=FALSE]); env<-max(abs(m1))
  ## (1) monthly-mean heterogeneity: sd across active sub-cells / env
  hm<-sapply(intM,function(m) sd(sub[active,m])/env)
  ## (2) sub-monthly heterogeneity: PCHIP each active sub-cell, spread of curves
  curves<-matrix(0,na,nM*Ns)
  for(k in seq_along(active)){ f<-pchip.fit.cell(edges,sub[active[k],]); pt<-1
    for(m in 1:nM) for(q in 1:Ns){ dt<-ss[q]*h[m]; curves[k,pt]<-f$a[m]*dt^2+f$b[m]*dt+f$c[m]; pt<-pt+1 } }
  idxint<-as.vector(sapply(intM,function(m)((m-1)*Ns+1):(m*Ns)))
  sm<-apply(curves[,idxint,drop=FALSE],2,sd)/env
  hetM<-c(hetM,median(hm)); subM<-c(subM,median(sm)); bnd<-c(bnd,bb); nact<-c(nact,na) }

cat("== Sub-grid (0.1deg-within-1deg) heterogeneity as prior uncertainty ==\n")
cat(sprintf("overall: monthly-mean spread/env median %.3f [%.3f,%.3f] ; sub-monthly spread/env median %.3f ; active subcells med %d\n",
    median(hetM),quantile(hetM,.25),quantile(hetM,.75),median(subM),median(nact)))
for(bb in c("tropics","temperate","boreal/polar")){ m<-bnd==bb
  cat(sprintf("  %-13s monthly-mean spread/env %.3f ; sub-monthly %.3f ; active subcells %d (n=%d)\n",
      bb,median(hetM[m]),median(subM[m]),round(median(nact[m])),sum(m))) }
cat("\ncompare other bands: structural(across-fitter) ~0.03 ; bootstrap(monthly-mean) ~0.01-0.06 ; ATP kriging var ~0.10-0.52\n")
