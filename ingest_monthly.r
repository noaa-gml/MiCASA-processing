#!/usr/bin/env Rscript

#SBATCH --account co2
#SBATCH --time 8:00:00
#SBATCH --ntasks 1
#SBATCH --mem 20g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.pera@noaa.gov

ct.setup()
ingest.weir.timestamp <- "Time-stamp: <orion-login-3.hpc.msstate.edu:/work2/noaa/co2/GFED-CASA/2025/MiCASA_v1/MiCASA_v1/ingest_monthly.r: 2025-06-06 12:24:22 MT>"

recompute.existing <- FALSE

srcdir <- 'portal.nccs.nasa.gov/monthly'

outdir <- "monthly_1x1"

if(!dir.exists(outdir)) {
  cat(sprintf("Creating output dir \"%s\"\n",outdir))
  dir.create(outdir,recursive=TRUE,showWarnings=TRUE)
}

start_year <- 2025
end_year <- 2025

year <- 2001
month <- 1
regs <- load.ncdf(sprintf('%s/%d/MiCASA_v1_flux_x3600_y1800_monthly_%d%02d.nc4',srcdir,year,year,month))

nms <- c("NPP","Rh","FIRE","FUEL")

# compute the area of a single grid cell with corners at lons and lats
archimedes <- function(lons,lats) {
  if(length(lons)!=2) {
    stop("Lons vector length not 2")
  }
  if(length(lats)!=2) {
    stop("Lats vector length not 2")
  }
  if(any(abs(range(lons))>pi)) {
    stop("abs(lons) vector exceeds pi")
  }
  if(any(abs(range(lats))>(pi/2))) {
    stop("abs(lats) vector exceeds pi/2")
  }
  R <- 6371007.2 # mean radius of the earth, in meters
  return((sin(lats[2]) - sin(lats[1])) * (lons[2]-lons[1]) * R^2)
}

aggregate.to.1x1 <- function(fld,gca) {
  retval <- matrix(0,360,180)
  for (jlat in 1:180) {
    inlats <- 1:10 + 10*(jlat-1)
    for(ilon in 1:360) {
      inlons <- 1:10 + 10*(ilon-1)
      retval[ilon,jlat] <- 0
      for (inlon in inlons) {
        retval[ilon,jlat] <- retval[ilon,jlat]+weighted.mean(fld[inlons,inlats],weights=gca[inlats],na.rm=T)
      }
      retval[ilon,jlat] <- retval[ilon,jlat]/10
    }
  }
  return(retval)
}

#
write.netcdf <- function(ncout,vars,vals,srcnm) {
  if(file.exists(ncout)) {
    file.remove(ncout)
  }
  ncf <- nc_create(ncout,vars=vars)
  
  # add attributes
  ncatt_put(ncf,0,"history",
            attval=sprintf("Created on %s\nby script '%s'",
                           format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),ingest.weir.timestamp),
            prec="text")

  ncatt_put(ncf,0,"Source",attval=srcnm,prec="text")

  for (nm in names(vars)) {
    ncvar_put(ncf,vars[[nm]],vals[[nm]])
  }

  nc_close(ncf)
}

gca <- rep(NA,1800)


lon.rad <- c(-0.05,0.05)*(pi/180)

for (ilat in 1:1800) {
  lat.rad <- (pi/180)*(regs$lat[ilat]+c(-0.05,0.05))
  gca[ilat] <- archimedes(lon.rad,lat.rad)
}

dim.lon <- ncdim_def("longitude","degrees_east",vals=seq(-179.5,179.5,1))
dim.lat <- ncdim_def("latitude","degrees_north",vals=seq(-89.5,89.5,1))

epoch <- ISOdatetime(1970,1,1,0,0,0,tz="UTC")
timeunits <- "seconds"
timeunits.difftime <- "secs"

vals <- list()

for (year in start_year:end_year) {

  for (month in 1:12) {
    
    this.date <- seq.midmon(year,month,year,month)

    # subtract the epoch to make timeunits-since
    time.vals <- as.numeric(difftime(this.date,epoch,units=timeunits.difftime)) 
    
    dim.time <- ncdim_def("time",
                          sprintf("%s since %s",timeunits,format(epoch,format="%Y-%m-%d %H:%M:%S UTC")),
                          vals=time.vals,unlim=TRUE)
    
    vars <- list()

    srcnm <- sprintf('%s/%d/MiCASA_vNRT_flux_x3600_y1800_monthly_%d%02d.nc4',srcdir,year,year,month)
    cat(sprintf("Processing %s...",basename(srcnm)))
    
    ncout <- sprintf('%s/MiCASA_vNRT_flux_x360_y180_monthly_%d%02d.nc',outdir,year,month)

    if(!recompute.existing & file.exists(ncout)) {
      cat(sprintf("Skipping existing \"%s\"\n",ncout))
      next
    }

    ncin <- load.ncdf(srcnm)

    for (nm in nms) {
      vars[[nm]] <- ncvar_def(name=nm,units="gC m^-2 s^-1",
                              dim=list(dim.lon,dim.lat,dim.time),
                              missval=-1e34,compression=9,
                              longname=attributes(ncin[[nm]])$long_name,
                              prec="float")
      
      vals[[nm]] <- aggregate.to.1x1(ncin[[nm]])*1000 # from kgC m-2 s-1 to gC m-2 s-1
    }

    write.netcdf(ncout,vars,vals,srcnm)
    cat('\n')
    
  } # month

} # year

