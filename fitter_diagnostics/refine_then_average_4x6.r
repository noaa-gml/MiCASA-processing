## Order-of-operations test at a COARSE 4degx6deg target grid (TM5-like), where
## within-cell phenological heterogeneity is large. fit-0.1deg-then-average vs
## aggregate-to-4x6-then-fit, both from the same 0.1deg blocks (60 lon x 40 lat
## = 2400 subcells per coarse cell). Proper cos(lat) area weights.
suppressWarnings(suppressMessages(library(ncdf4)))
source("lib/pchip_fit.r")
yrs<-2018:2021; mons<-as.vector(sapply(yrs,function(y)sprintf("%04d%02d",y,1:12))); nM<-length(mons)
edges<-as.numeric(seq(as.POSIXct(sprintf("%d-01-01",yrs[1]),tz="UTC"),by="1 month",length.out=nM+1))
fpath<-function(mm)sprintf("portal.nccs.nasa.gov/monthly/%s/MiCASA_v1_flux_x3600_y1800_monthly_%s.nc4",substr(mm,1,4),mm)
nlon<-60; nlat<-40; nsub<-nlon*nlat                     # 4deg lat x 6deg lon
## candidate coarse cells over land lat/lon; filter by flux after a probe read
set.seed(2); cells<-list()
f<-nc_open(fpath(mons[7])); probe<- -2*ncvar_get(f,"NPP"); nc_close(f)
for(I in 1:60) for(J in 1:45){ li<-((I-1)*nlon+1):(I*nlon); lj<-((J-1)*nlat+1):(J*nlat)
  if(mean(abs(probe[li,lj]))>1e-9) cells[[length(cells)+1]]<-c(I,J) }
cells<-cells[sample(length(cells),min(15,length(cells)))]; ns<-length(cells)
GPP<-array(0,c(ns,nsub,nM))
for(m in 1:nM){ f<-nc_open(fpath(mons[m])); npp<-ncvar_get(f,"NPP"); nc_close(f)
  for(s in 1:ns){ I<-cells[[s]][1]; J<-cells[[s]][2]; li<-((I-1)*nlon+1):(I*nlon); lj<-((J-1)*nlat+1):(J*nlat)
    GPP[s,,m]<- -2*as.vector(npp[li,lj]) } }
peak_q<-function(a,b,c0,D){fL<-c0;fR<-a*D^2+b*D+c0;sv<-ifelse(a!=0,-b/(2*a),-1);intr<-a!=0&sv>0&sv<D;fV<-ifelse(intr,c0-b^2/(4*a),0);pmax(abs(fL),abs(fR),ifelse(intr,abs(fV),0))}
D<-diff(edges); intM<-6:(nM-6); os_c<-os_f<-numeric(ns)
for(s in 1:ns){ J<-cells[[s]][2]; lat40<- -90 + ((J-1)*nlat + (1:nlat) - 0.5)*0.1
  w<-rep(cos(lat40*pi/180),each=nlon); w<-w/sum(w)     # cell order: lon fastest -> rep each=nlon
  sub<-GPP[s,,]
  gm<-colSums(w*sub); fc<-pchip.fit.cell(edges,gm)
  af<-bf<-cf<-numeric(nM)
  for(k in 1:nsub){ fk<-pchip.fit.cell(edges,sub[k,]); af<-af+w[k]*fk$a; bf<-bf+w[k]*fk$b; cf<-cf+w[k]*fk$c }
  um<-abs(gm); env<-pmax(c(um[1],um[1:(nM-1)]),um,c(um[2:nM],um[nM]))
  ok<-env>1e-12 & (1:nM)%in%intM
  os_c[s]<-max((peak_q(fc$a,fc$b,fc$c,D)/env)[ok]); os_f[s]<-max((peak_q(af,bf,cf,D)/env)[ok])
}
cat(sprintf("4x6 fit-coarse:        overshoot peak/env median %.2f 90th %.2f max %.2f\n",median(os_c),quantile(os_c,.9),max(os_c)))
cat(sprintf("4x6 fit-0.1-then-avg:  overshoot peak/env median %.2f 90th %.2f max %.2f\n",median(os_f),quantile(os_f,.9),max(os_f)))
cat(sprintf("per-cell reduction (coarse-fine): median %.3f ; fine lower in %.0f%% of %d cells\n",median(os_c-os_f),100*mean(os_f<os_c),ns))
