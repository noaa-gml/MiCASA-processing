## Order-of-operations test: fit at 0.1deg then area-average to 1deg  vs
## aggregate to 1deg then fit. Both derived from the SAME 0.1deg blocks, so the
## only difference is fit-order (the fit is nonlinear -> they differ). Measures
## whether fit-fine-then-average lowers overshoot + guarantees sign-definiteness.
suppressWarnings(suppressMessages(library(ncdf4)))
source("lib/pchip_fit.r")
yrs<-2018:2021; mons<-as.vector(sapply(yrs,function(y)sprintf("%04d%02d",y,1:12))); nM<-length(mons)
edges<-as.numeric(seq(as.POSIXct(sprintf("%d-01-01",yrs[1]),tz="UTC"),by="1 month",length.out=nM+1))
fpath<-function(mm)sprintf("portal.nccs.nasa.gov/monthly/%s/MiCASA_v1_flux_x3600_y1800_monthly_%s.nc4",substr(mm,1,4),mm)

## sample 1deg land cells across latitude bands (from the existing 1deg fit)
P<-new.env(); load("fit.pchip.rda",envir=P); a<-matrix(P$piqsfit.gpp$a,360*180,length(P$piqsfit.time))
land1<-which(rowSums(abs(a))>0); set.seed(4)
j1<-((land1-1)%/%360)+1; latc<- -90+j1-0.5
samp<-c(sample(land1[abs(latc)<23.5],14),sample(land1[abs(latc)>=23.5&abs(latc)<50],14),sample(land1[abs(latc)>=50],14))
ns<-length(samp); i1<-((samp-1)%%360)+1; j1<-((samp-1)%/%360)+1

## accumulate 0.1deg blocks [ns,100,nM] for gpp
GPP<-array(0,c(ns,100,nM))
for(m in 1:nM){ f<-nc_open(fpath(mons[m])); npp<-ncvar_get(f,"NPP"); nc_close(f)  # [lon3600, lat1800]
  for(s in 1:ns){ li<-((i1[s]-1)*10+1):(i1[s]*10); lj<-((j1[s]-1)*10+1):(j1[s]*10)
    GPP[s,,m]<- -2*as.vector(npp[li,lj]) } }
## area weights per 0.1deg row (cos lat), per 1deg cell
peak_q<-function(a,b,c0,D){fL<-c0;fR<-a*D^2+b*D+c0;sv<-ifelse(a!=0,-b/(2*a),-1);intr<-a!=0&sv>0&sv<D;fV<-ifelse(intr,c0-b^2/(4*a),0);pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))}
D<-diff(edges); intM<-6:(nM-6)                 # interior months w/ full context
os_c<-os_f<-signf<-numeric(ns)
for(s in 1:ns){ lat10<- -90 + ((j1[s]-1)*10 + (1:10) - 0.5)*0.1
  w<-rep(cos(lat10*pi/180),each=1); w<-rep(cos(lat10*pi/180),times=10)  # 100 cells: lon x lat
  w<-w/sum(w)
  sub<-GPP[s,,]                                 # [100, nM]
  ## coarse: aggregate then fit
  gm<-colSums(w*sub); fc<-pchip.fit.cell(edges,gm)
  ## fine: fit each subcell then area-average coeffs
  af<-bf<-cf<-numeric(nM)
  for(k in 1:100){ fk<-pchip.fit.cell(edges,sub[k,]); af<-af+w[k]*fk$a; bf<-bf+w[k]*fk$b; cf<-cf+w[k]*fk$c }
  um<-abs(gm); env<-pmax(c(um[1],um[1:(nM-1)]),um,c(um[2:nM],um[nM]))
  pc<-peak_q(fc$a,fc$b,fc$c,D); pf<-peak_q(af,bf,cf,D)
  ok<-env>1e-12 & (1:nM)%in%intM
  os_c[s]<-max((pc/env)[ok]); os_f[s]<-max((pf/env)[ok])
  ## sign-definiteness of the fine (averaged) GPP flux: should be <=0 everywhere
  fl<-af[intM]*0; for(mm in intM) fl<-c(fl, af[mm]*(D[mm]/2)^2*c(0,1)+bf[mm]*(D[mm]/2)*c(0,1)+cf[mm]) # sample
  signf[s]<-as.integer(max(c(cf[intM], af[intM]*D[intM]^2+bf[intM]*D[intM]+cf[intM]))>1e-12*max(um))
}
cat(sprintf("fit-at-1deg (coarse):  overshoot peak/env  median %.2f  90th %.2f  max %.2f\n",median(os_c),quantile(os_c,.9),max(os_c)))
cat(sprintf("fit-0.1deg-then-avg:   overshoot peak/env  median %.2f  90th %.2f  max %.2f\n",median(os_f),quantile(os_f,.9),max(os_f)))
cat(sprintf("per-cell reduction (coarse-fine): median %.3f  ; fine lower in %.0f%% of cells\n",median(os_c-os_f),100*mean(os_f<os_c)))
cat(sprintf("fine-avg GPP wrong-sign cells: %.0f%% (PCHIP avg of <=0 should be 0)\n",100*mean(signf)))
