#!/usr/bin/env Rscript

ct.setup()

# Time-stamp: <~/Documents/co2-orion/GFED-CASA/2025/MiCASA_v1/check_daily_downloads.r: 13 May 2025 16:03:49 -06:00>"

srcdir <- 'portal.nccs.nasa.gov/daily'

for (year in 2001:2024) {
  nms <- c("NPP","Rh","FIRE","FUEL")

  dpm <- days.in.month(year)
  pb <- progress.bar.start(message=sprintf("%d: %d days",year,sum(dpm)),nx=sum(dpm))

  iday <- 0
  for (month in 1:12) {
    
    for (day in 1:dpm[month]) {
      
      iday <- iday + 1
      
      this.date <- ISOdatetime(year,month,day,0,0,0,tz="UTC")+86400/2

      src_name <- sprintf('%s/%d/%02d/MiCASA_v1_flux_x3600_y1800_daily_%d%02d%02d.nc',srcdir,year,month,year,month,day)
      src_name4 <- sprintf('%s/%d/%02d/MiCASA_v1_flux_x3600_y1800_daily_%d%02d%02d.nc4',srcdir,year,month,year,month,day)
      if(!file.exists(src_name) && !file.exists(src_name4)) {
        cat(sprintf("no file %s\n",src_name))
      }
      pb <- progress.bar.print(pb,iday)
      
    } # day
    
  } # month
  progress.bar.end(pb)  
  #} # year

} # if-else

