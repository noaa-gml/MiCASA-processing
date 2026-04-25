#!/usr/bin/env Rscript

#SBATCH --account co2
#SBATCH --time 8:00:00
#SBATCH --ntasks 1
#SBATCH --mem 40g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ashley.pera@noaa.gov

# ALERT: This is set to FastTrack data, and ONLY after April. Beware.

script.name <- "Time-stamp: <orion-login-3.hpc.msstate.edu:/work2/noaa/co2/GFED-CASA/2025/MiCASA_v1/diurnalize-ERA5.r: 09 Jul 2025 (Wed) 19:15:22 UTC>"

gfed.version <- "MiCASA_v1"
product.name <- "MiCASA_v1_flux_x360_y180_monthly"

# Unlike previous versions, this one gets the relevant meteo (ssrd,
# t2m, stl1, swvl1) from TM5 netCDF files in the METEO directory.
# These have hourly 1x1 resolution.

era5dir <- sprintf("%s/METEO/tm5-nc/ec/ea_0005/h06h18tr1/sfc/glb100x100", CARBONTRACKER)
era5template <- "YYYY/MM/VVV_YYYYMMDD_00p01.nc"
metstr <- "ERA5"

yr <- Sys.getenv("diurn_year")

if(nchar(yr)==0) {
  for (yr in 2025) {
    system(sprintf("sbatch -J d-%d-%s --export=ALL,gfed_version=%s,diurn_year=%d diurnalize-ERA5.r",
                   yr,gfed.version,gfed.version,yr))
  }
} else {

  yr <- as.integer(yr)
  
  clim.yrs <- c(2000,2025)
  # TODO: make this one launch-directory dependant
  in.dir <- sprintf("%s/GFED-CASA/2025/%s/monthly_1x1",CARBONTRACKER,gfed.version)
  out.dir <- sprintf("%s/GFED-CASA/2025/%s/%s",CARBONTRACKER,gfed.version,metstr)
  
  ct.setup()

  # load PIQS coefficients to smooth month-month variability
  load(sprintf("%s/fit.piqs.rda",'.'))
  piqsfit.time <- epoch.seconds.to.POSIX(piqsfit.time)
  piqsfit.lts <- as.POSIXlt(piqsfit.time)
  
  lon.dim <- ncdim_def("longitude","degrees_east",vals=seq(-179.5,179.5,1))
  lat.dim <- ncdim_def("latitude","degrees_north",vals=seq(-89.5,89.5,1))

  epoch <- ISOdatetime(1970,1,1,0,0,0,tz="UTC")
  timeunits <- "days"
  
  dir.create(out.dir,showWarnings = FALSE, recursive = TRUE)

  for (mon in 4:12) {

    cat(sprintf("%d/%02d\n",yr,mon))
    monstr <- sprintf("%d%02d",yr,mon)

    current.time <- ISOdatetime(yr,mon,1,0,0,0,tz="UTC")
    
    ncname.out <- sprintf("%s/fluxes_%s.nc",out.dir,monstr)
    
    if( yr %in% clim.yrs) {
      fname <- sprintf("%s/NPPclim.nc",in.dir)
      if(!file.exists(fname)) {
        stop(sprintf("%d-%02d:  %s does not exist.",yr,mon,fname))
      }
      foo <- load.ncdf(fname)
      # change sign so that negative is a sink,
      # change units from gC m-2 s-1 to mol m-2 s-1
      gpp.clim <- -2*foo$NPPCLIM/12
      
      fname <- sprintf("%s/Rhclim.nc",in.dir)
      if(!file.exists(fname)) {
        stop(sprintf("%d-%02d:  %s does not exist.",yr,mon,fname))
      }
      foo <- load.ncdf(fname)
      
      # change units from gC m-2 s-1 to mol m-2 s-1
      rh.clim <- foo$RHCLIM/12
      rtot.clim <- rh.clim - 0.5*gpp.clim # total respiration in mol m-2 s-1 is heterotrophic plus autotrophic (half of GPP, equal to NPP)
      
      # The climatological fields will have a third dimension with 12 months
      rtot.mn <- rtot.clim[,,mon]
      gpp.mn <- gpp.clim[,,mon]
      
      rm(foo)
    } else {
      fname <- sprintf("%s/%s_%s.nc",in.dir,product.name,monstr)
      foo <- load.ncdf(fname)
      gpp.mn <- -2*foo$NPP/12 # from gC m-2 s-1 to mol m-2 s-1
      rh.mn <- foo$Rh/12 # from gC m-2 s-1 to mol m-2 s-1
      rtot.mn <- rh.mn - 0.5*gpp.mn # total respiration in mol m-2 s-1 is heterotrophic plus autotrophic (half of GPP, equal to NPP)
      rm(foo)
    }          
    cat(sprintf("Finished reading %s...\n",fname))
    
    # get ssr and q10
    varnms <- c("t2m","ssrd","stl1","swvl1")
    mets <- list()
    dpm <- days.in.month(yr)[mon]
    times <- rep(NA,dpm*24)
    for (day in 1:dpm) {
      for(varnm in varnms) {
        e5nm <- gsub("YYYY",sprintf("%d",yr),era5template)
        e5nm <- gsub("MM",sprintf("%02d",mon),e5nm)
        e5nm <- gsub("DD",sprintf("%02d",day),e5nm)
        e5nm <- gsub("VVV",varnm,e5nm)
        ncname.in <- sprintf("%s/%s",era5dir,e5nm)
        foo <- load.ncdf(ncname.in)
        if(is.null(mets[[varnm]])) {
          mets[[varnm]] <- array(NA,dim=c(360,180,24*dpm))
        }
        k0 <- 1+(day-1)*24
        k1 <- 23+k0
        mets[[varnm]][,,k0:k1] <- foo[[varnm]]
        if(varnm==varnms[1]) {
          # This fails between 2022-04-30 and 2024-06-29 due
          # to faulty time axes in the ERA5 files.
          #times[k0:k1] <- foo$time
          times[k0:k1] <- seq(ISOdatetime(yr,mon,day,0,0,0),
                              ISOdatetime(yr,mon,day,23,0,0),
                              by="1 hour")
                              
        }
      }
    }
    times <- epoch.seconds.to.POSIX(times)
    times <- times+1800 # add 30 minutes to get to center of hour
    nslots <- length(times)
    
    q10 <- 1.5^((mets$t2m-273.15)/10.0)

    # compute monthly means of ssr and q10
    q10.mn <- apply(q10,c(1,2),mean)
    ssr.mn <- apply(mets$ssrd,c(1,2),mean)

    lx <- which(ssr.mn==0)
    if(length(lx)>0) {
      # This used to be ssr.mn[lx] <- NA
      # but ERA5 north of about 70N is
      # identically zero, and ERA-i was 
      # 1e-16 W/m2. This causes trouble
      # for the GPP calculation.
      ssr.mn[lx] <- 1e-16
    }
    cat("GPP mean summary:\n")
    print(summary(as.vector(gpp.mn)))
    cat("RTOT mean summary:\n")
    print(summary(as.vector(rtot.mn)))
    factor.gpp <- gpp.mn/ssr.mn
    factor.resp <- rtot.mn/q10.mn
    
    lx <- which(is.na(factor.gpp))
    if(length(lx)>0) {
      factor.gpp[lx] <- 0
    }

    gpp <- array(NA,dim=dim(mets$ssrd))
    resp <- array(NA,dim=dim(mets$ssrd))
    nee <- array(NA,dim=dim(mets$ssrd))
    qgpp <- array(NA,dim=dim(mets$ssrd))
    qresp <- array(NA,dim=dim(mets$ssrd))


    # now substract mean monthly fluxes and insert the smoothed PIQS fit
    if((current.time >= min(piqsfit.time)) & (current.time <= max(piqsfit.time))) {
      imon <- which(piqsfit.time == current.time)
      gpp.a <- piqsfit.gpp$a[,,imon]
      gpp.b <- piqsfit.gpp$b[,,imon]
      gpp.c <- piqsfit.gpp$c[,,imon]
      resp.a <- piqsfit.resp$a[,,imon]
      resp.b <- piqsfit.resp$b[,,imon]
      resp.c <- piqsfit.resp$c[,,imon]
    } else {
      monseq <- which((piqsfit.lts$mon + 1) == mon)
      gpp.a <- apply(piqsfit.gpp$a[,,monseq],c(1,2),mean)
      gpp.b <- apply(piqsfit.gpp$b[,,monseq],c(1,2),mean)
      gpp.c <- apply(piqsfit.gpp$c[,,monseq],c(1,2),mean)
      resp.a <- apply(piqsfit.resp$a[,,monseq],c(1,2),mean)
      resp.b <- apply(piqsfit.resp$b[,,monseq],c(1,2),mean)
      resp.c <- apply(piqsfit.resp$c[,,monseq],c(1,2),mean)
    }
    for(islot in 1:nslots) {

      dt <- as.numeric(times[islot])-as.numeric(times[1])

      qmod.gpp <- (gpp.a*(dt)^2+gpp.b*(dt)+gpp.c)/12
      qmod.resp <- (resp.a*(dt)^2+resp.b*(dt)+resp.c)/12

      #        gpp[,,islot] <- mets$ssrd[,,islot]*qmod.gpp/mets$ssrd.mn
      #        resp[,,islot] <- q10[,,islot]*qmod.resp/q10.mn
      gpp[,,islot] <- mets$ssrd[,,islot]*gpp.mn/ssr.mn
      resp[,,islot] <- q10[,,islot]*rtot.mn/q10.mn
      
      gpp[,,islot] <- gpp[,,islot]-gpp.mn+(gpp.a*(dt)^2+gpp.b*(dt)+gpp.c)/12
      resp[,,islot] <- resp[,,islot]-rtot.mn+(resp.a*(dt)^2+resp.b*(dt)+resp.c)/12
      nee[,,islot] <- gpp[,,islot]+resp[,,islot]
      qgpp[,,islot] <- qmod.gpp
      qresp[,,islot] <- qmod.resp
    }

    #            if(any(is.na(gpp[80,160,]))) {
    # When switching from ERA-i to ERA5, had land GPPs
    # which were NA.  This was a diagnostic for that
    # problem. Issue was ERA5 SSRD identically zero, when
    # ERA-i set it to 1e-16. NOTE: only works for
    # GFED4.1s; suspect coastal gridpoint that is ocean in
    # GFEDCMS.
    #                stop("GPP NAs")
    #            }
    lx <- which(is.na(nee))
    if(length(lx)>0) {
      nee[lx] <- 0
    }
    cat(sprintf("  %s\n",ncname.out))
    
    date.vals <- as.numeric(difftime(times,epoch,units=timeunits)) # subtract the epoch to make days-since
    decimal.date <- POSIX.to.decimal(times)

    date.dim <- ncdim_def("time",sprintf("%s since %s",timeunits,format(epoch,format="%Y-%m-%d %H:%M:%S UTC")),vals=date.vals,unlim=TRUE)

    vars <- list()
    vars$dd <- ncvar_def(name="decimal_date",units="years",dim=list(date.dim),
                         missval=-1e34,compression=9,longname="decimal_date",prec="double")

    vars$gpp <- ncvar_def(name="GPP",units="mol m-2 s-1",dim=list(lon.dim,lat.dim,date.dim),
                          missval=-1e34,compression=9,longname="gross_primary_production, twice the modeled NPP, positive is source to atm (contrary to conventional definition)")

    vars$resp <- ncvar_def(name="resp",units="mol m-2 s-1",dim=list(lon.dim,lat.dim,date.dim),
                           missval=-1e34,compression=9,longname="ecosystem_respiration, as sum of Rhetero and Rauto, positive is source to atm")

    vars$nee <- ncvar_def(name="NEE",units="mol m-2 s-1",dim=list(lon.dim,lat.dim,date.dim),
                          missval=-1e34,compression=9,longname="NEE=GPP+RESP, positive is source to atm, as is each component")

    vars$qgpp <- ncvar_def(name="QGPP",units="mol m-2 s-1",dim=list(lon.dim,lat.dim,date.dim),
                           missval=-1e34,compression=9,longname="gross_primary_production model")

    vars$qresp <- ncvar_def(name="qresp",units="mol m-2 s-1",dim=list(lon.dim,lat.dim,date.dim),
                            missval=-1e34,compression=9,longname="ecosystem_respiration model")

    vars$ssr <- ncvar_def(name="ssr",units="W/m2",dim=list(lon.dim,lat.dim,date.dim),
                          missval=-1e34,compression=9,longname="ERA5 surface shortwave radiation downward")

    vars$t2m <- ncvar_def(name="t2m",units="K",dim=list(lon.dim,lat.dim,date.dim),
                          missval=-1e34,compression=9,longname="ERA5 2-meter air temperature")

    vars$stl1 <- ncvar_def(name="stl1",units="K",dim=list(lon.dim,lat.dim,date.dim),
                           missval=-1e34,compression=9,longname="ERA5 soil level 1 temperature (0-7 cm)")

    vars$swvl1 <- ncvar_def(name="swvl1",units="m3/m3",dim=list(lon.dim,lat.dim,date.dim),
                            missval=-1e34,compression=9,longname="ERA5 soil level 1 volumetric moisture content (0-7 cm)")

    if(file.exists(ncname.out)) {
      cat(sprintf("Removing existing output file \"%s\"\n",ncname.out))
      file.remove(ncname.out)
    }
    ncf <- nc_create(ncname.out,vars=vars)
    
    ncatt_put(ncf,0,"history",
              attval=sprintf("Created on %s\nby script '%s'",
                             format(Sys.time(), "%a %b %d %Y %H:%M:%S %Z"),script.name),
              prec="text")

    ncatt_put(ncf,0,"meteo_source_directory",
              attval=era5dir,prec='text')

    ncvar_put(ncf,vars$dd,vals=decimal.date)
    ncvar_put(ncf,vars$gpp,vals=gpp)
    ncvar_put(ncf,vars$resp,vals=resp)
    ncvar_put(ncf,vars$nee,vals=nee)
    ncvar_put(ncf,vars$qgpp,vals=qgpp)
    ncvar_put(ncf,vars$qresp,vals=qresp)

    ncvar_put(ncf,vars$ssr,vals=mets$ssrd)
    ncvar_put(ncf,vars$t2m,vals=mets$t2m)
    ncvar_put(ncf,vars$stl1,vals=mets$stl1)
    ncvar_put(ncf,vars$swvl1,vals=mets$swvl1)
    
    nc_close(ncf)


    #      cat(sprintf("  %s\n",hdfname.facs.out))
    
    #     cmd <- sprintf("%s %s|%s -b -o %s",ncdump,ncname.facs.out,hdfgen,hdfname.facs.out)
    #     system(cmd)
    
    cat("\n")

  } # month loop
  #    } # year loop
} # else
