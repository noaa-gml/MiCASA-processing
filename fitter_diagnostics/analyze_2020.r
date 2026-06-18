#!/usr/bin/env Rscript
## Product-level head-to-head: PCHIP vs PPM diurnalized NEE, full-year 2020.
## (1) annual global NEE budget (mass check), (2) per-month global mean,
## (3) spatial annual-mean NEE difference, (4) month-boundary aliasing via FFT.
suppressWarnings(suppressMessages(library(ncdf4)))
DP<-"ERA5_2020_pchip"; DM<-"ERA5_2020_ppm"
dpm<-c(31,29,31,30,31,30,31,31,30,31,30,31)   # 2020 leap

## grid-cell area (m^2), 1deg
R<-6.371e6; lat<-(-89.5):89.5; dlat<-pi/180
area_lat <- R^2*(pi/180)*(sin((lat+0.5)*pi/180)-sin((lat-0.5)*pi/180))  # per 1deg lon ring cell
A <- matrix(rep(area_lat,each=360),360,180)    # [lon,lat]

gC<-12.011; secyr<-sum(dpm)*86400
budP<-budM<-0; mP<-mM<-numeric(12)
annP<-array(0,c(360,180)); annM<-array(0,c(360,180))
for(m in 1:12){
  fP<-nc_open(sprintf("%s/fluxes_2020%02d.nc",DP,m)); fM<-nc_open(sprintf("%s/fluxes_2020%02d.nc",DM,m))
  neeP<-ncvar_get(fP,"NEE"); neeM<-ncvar_get(fM,"NEE"); nc_close(fP); nc_close(fM)
  mmP<-apply(neeP,c(1,2),mean); mmM<-apply(neeM,c(1,2),mean)   # month-mean NEE per cell
  w<-dpm[m]/sum(dpm)
  annP<-annP+w*mmP; annM<-annM+w*mmM
  ## global month-mean flux (area-weighted), mol m-2 s-1
  mP[m]<-sum(mmP*A)/sum(A); mM[m]<-sum(mmM*A)/sum(A)
}
## annual budget PgC/yr = sum(annual-mean NEE [mol/m2/s] * area * secyr * gC) *1e-15
budP<-sum(annP*A)*secyr*gC*1e-15
budM<-sum(annM*A)*secyr*gC*1e-15
cat(sprintf("ANNUAL GLOBAL NEE BUDGET 2020:  PCHIP %.4f PgC/yr   PPM %.4f PgC/yr   diff %.2e PgC/yr (%.3f%%)\n",
            budP,budM,budM-budP,100*(budM-budP)/abs(budP)))
cat("Per-month global-mean NEE (mol m-2 s-1), PCHIP vs PPM (rel diff):\n")
for(m in 1:12) cat(sprintf("  2020-%02d  %.4e  %.4e  (%.2e)\n",m,mP[m],mM[m],abs(mM[m]-mP[m])/max(abs(mP[m]),1e-30)))

## spatial annual-mean difference
land<-abs(annP)>1e-9 | abs(annM)>1e-9
d<-(annM-annP)[land]
cat(sprintf("\nAnnual-mean NEE difference PPM-PCHIP (land, mol m-2 s-1):\n  rms %.3e  max|.| %.3e  median|.| %.3e\n",
            sqrt(mean(d^2)),max(abs(d)),median(abs(d))))
reldiff<-abs(annM-annP)/pmax(abs(annP),1e-12)
cat(sprintf("  median rel|diff| over land annual-mean: %.3f%%\n",100*median(reldiff[land])))

## ---- month-boundary aliasing: FFT of full-year hourly NEE at 3 cells ----
read_cell_year<-function(dir,i,j){ v<-c(); for(m in 1:12){ f<-nc_open(sprintf("%s/fluxes_2020%02d.nc",dir,m))
  x<-ncvar_get(f,"NEE",start=c(i,j,1),count=c(1,1,-1)); nc_close(f); v<-c(v,as.numeric(x)) }; v }
## pick cells where annual diff largest + a typical mid-lat
ord<-order(abs(annM-annP)*land,decreasing=TRUE)
cells<-c(ord[1], which(land)[which.min(abs(reldiff[which(land)]-median(reldiff[land])))])
cat("\nMonth-boundary aliasing (FFT of hourly NEE): power at period=1 month / total power\n")
for(idx in cells){ i<-((idx-1)%%360)+1; j<-((idx-1)%/%360)+1
  yp<-read_cell_year(DP,i,j); ym<-read_cell_year(DM,i,j); n<-length(yp)
  ## remove diurnal (24h) + annual mean, focus on sub-seasonal; periodogram
  fp<-Mod(fft(yp-mean(yp)))^2; fm<-Mod(fft(ym-mean(ym)))^2
  freq<-(0:(n-1))/n                       # cycles per hour
  ## monthly band ~ 1/(30.4*24) cph; sum power in +-20% band
  fmonth<-1/(30.4*24); band<-which(freq>0 & abs(freq-fmonth)<0.2*fmonth)
  cat(sprintf("  lon %.1f lat %.1f:  PCHIP monthly-band frac %.3e   PPM %.3e   (ratio PPM/PCHIP %.2f)\n",
              -180+i-0.5,-90+j-0.5, sum(fp[band])/sum(fp[-1]), sum(fm[band])/sum(fm[-1]),
              (sum(fm[band])/sum(fm[-1]))/(sum(fp[band])/sum(fp[-1]))))
}
