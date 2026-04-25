README written by Ash Pera 2025-05-13 17:03:37

##########################
# Overview
##########################

This directory holds scripts to take raw MiCASA data and process it for use in CarbonTracker
https://gml.noaa.gov/ccgg/carbontracker/documentation.php

MiCASA Land Carbon Flux
Global, daily and monthly mean 0.1 degree resolution carbon fluxes from net primary production (NPP),
heterotrophic respiration (Rh), wildfire emissions (FIRE), fuel wood burning emissions (FUEL),
net ecosystem exchange (NEE), and net biosphere exchange (NBE) derived from the MiCASA model, version 1
https://earth.gov/ghgcenter/data-catalog/micasa-carbonflux-grid-v1

##########################
# Flowchart
##########################

symlink_old_micasa.sh
        |
download_and_check.sh--------\
(download.sh)                |
(check_daily_downloads.r)    |
(check_hashes.py)            |
(check_unchanged.sh)         |
        |                    |
 ingest_monthly.r      ingest_byyear.r
        |                    |
  cat_monthly.sh    compute_daily_clim.sh 
 (check_bounds.sh)           |
        |            link_daily_clim.sh   
 compute_clim.sh
        |      
  write_piqs.r 
        |
diurnalize-ERA5.r
        |                 test_gca.r  
  daysplitter.sh

##########################
# Programs
##########################

download_and_check.sh:
    Calls download.sh, check_daily_downloads.sh, check_hashes.py, and check_unchanged.sh

test_gca.r:
    Compute the lat-lon grid cell areas, Load CARBONTRACKER/tools/shared/aux/regions.nc, and print the relative error.
    Should be on the order of 2.7e-6, consistent with single-precision numerics.

download.sh:
    Make the from_weir dir, and wget https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/daily to
    get daily MiCASA data (in nc4 format). Can also get monthlies.
    Should (but doesn't) check the hash.

download.sh-orig:
    Same, but from https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/MICASA_D_FLUX/ instead. Depricated.

check_daily_downloads.r:
    Verify that the data downloaded to from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/daily/
    have NPP, Rh, FIRE, and FUEL for every day from 2001-2023.

check_hashes.py:
   Runs through all the monthly and daily data and checks their sha256 hash against the ones provided.

check_unchanged.sh:
   Checks what has changed in the headers from a referance file (January 2012) 

check_bounds.sh:
   Computes a simple average of data, and print out. Does NOT correctly average by cell areas. Called by cat_monthly.sh

ingest_byyear.r:
    Take raw from_weir/~/daily/*.nc4 files for a given year and--for every month and day--aggregate NPP, Rh, FIRE, and FUEL
    from the original 0.1 deg grid to a 1 deg grid. Write the results to daily_1x1/MiCASA_v1_flux_x360_y180_daily_YYYYMMDD.nc

    Parallel job; if INGEST_YEAR is not in env, launch 23 recursive jobs where it's set to [2001:2023].

ingest_monthly.r:
    Take raw from_weir/~/monthly/*.nc4 files from [2001:2023] for every month and aggregate NPP, Rh, FIRE, and FUEL
    from the original 0.1 deg grid to a 1 deg grid. Write the results to monthly_1x1/MiCASA_v1_flux_x360_y180_daily_MMDD.nc

ingest_monthly_special_201801.r:
    Ingest only January 2018, load from from_weir/~/monthly/*.nc instead of .nc4.
    Depricated, 1-time bug fix.

ingest.r:
    For every day of every month in the years [2001:2023], take raw from_weir/~/monthly/*.nc4 files and aggregate NPP, Rh,
    FIRE, and FUEL from the original 0.1 deg grid to a 1 deg grid; skipping files that exist already.
    Write the results to monthly_1x1/MiCASA_v1_flux_x360_y180_daily_MMDD.nc

    Hooked job--submit a dependency recursive job to continue if it didn't finish.
    Deprecated, replaced with ingest_monthly and ingest_byyear.

cat_monthly.sh :
    Takes monthly data and combines to monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc. Calls check_bounds.sh

compute_clim.sh:
    Using Ferret, load monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc and take a modulo average of each
    month of NPP and Rh, writing the outputs to monthly_1x1/Rhclim.nc and monthly_1x1/NPPclim.nc

compute_daily_clim.sh:
    Use ncea to average every day of the year in daily_1x1/MiCASA_v1_flux_x360_y180_dailyYYYYMMDD.nc across all years,
    and save output as daily_1x1/MiCASA_v1_flux_x360_y180_daily_0000MMDD.nc

link_daily_clim.sh:
    For every day [2000:2024], if there is no MiCASA_v1_flux_x360_y180_daily_YYYYMMDD.nc, link it to the 
    MiCASA_v1_flux_x360_y180_daily_0000MMDD.nc climate average for that day of the year.

write_piqs.r:
    Load monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc and for every grid cell do a piecewise integral quadratic splines (PIQS)
    fit of GPP and rtot, and save the results to fit.piqs.rda. Can also make a pdf of plots of TSER quadradic fits from 2000-2016.
    Skips cells where nee2 = 0. 
    Time Series Extrinsic Regression (TSER), Gross Primary Production (GPP = -2*NPP), total respiration (rtot = Rh + NPP)
    nee1 (Rh - NPP), and nee2 (gpp + rtot). 
    https://gml.noaa.gov/ccgg/carbontracker/documentation.php#tth_sEc2.2

diurnalize-ERA5.r:
    For a year, load PIQS coefficients (from fit.piqs.rda) to smooth month-month variability. Load climatic averages of
    NPP and respiration (from NPP and Rhclim.nc). Load ERA5 surface solar radiation (ssrd), volumetric soil water layer (swvl1),
    and average 2-meter temperature (t2m) from CARBONTRACKER/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100/YYYY/MM/VVV_YYYYMMDD_00p01.nc,
    to get ssr and q10, and find the monthly means. Subtract those, and insert the smoothed PIQS fit. Use to diernalize gpp, resp, nee,
    qgpp and qresp. Write GPP, resp, NEE, QGPP, qresp, ssr, t2m, stl1 and swvl1  to ERA5/fluxes_YYYYMM.nc.

    Parallel job; if diurn_year is not in env, launch 25 recursive jobs where it's set to [2000:2024].

daysplitter.sh:
    For every year [2000:2004] and month in ERA5/fluxes_YYYYMM.nc, split into ERA5/MiCASA_v1.nee.YYYYMMDD.nc dailies,
    keeping only NEE. NEE=GPP+RESP (positive is source to atm). GPP=gross_primary_production, twice the modeled NPP,
    RESP=ecosystem_respiration, sum of Rhetero and Rauto.

    Use ncdump, grep and sed to find the number of timeslices, divide by 24 for the days, and split/subset with ncks.




##########################
# Data
##########################

from_weir/portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/netcdf/daily/YYYY/MM/MiCASA_v1_flux_x3600_y1800_daily_YYYYMMDD.nc4
    Created by download.sh, intial data. See end of this file for header dump.

daily_1x1/MiCASA_v1_flux_x360_y180_dailyYYYYMMDD.nc
    Created by ingest_byyear.rm from inital daily data

daily_1x1/MiCASA_v1_flux_x360_y180_daily_0000MMDD.nc
    Created by compute_daily_clim.sh from daily_1x1/MiCASA_v1_flux_x360_y180_dailyYYYYMMDD.nc

monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_YYYYMM.nc
    Created by ingest_monthly.r from inital monthly data (TODO: where's that from?)

monthly_1x1/NPPclim.nc monthly_1x1/NPPclim.nc 
    Created by compute_clim.sh (Ferret) from monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc

monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc
    [MISSING, but uses NCO from monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_YYYYMM.nc]

ERA5/MiCASA_v1.nee.YYYYMMDD.nc
    Created by daysplitter.sh from ERA5/fluxes_YYYYMM.nc

ERA5/fluxes_YYYYMM.nc
    Created by diurnalize-ERA5.r from and 
    CARBONTRACKER/METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100/YYYY/MM/VVV_YYYYMMDD_00p01.nc and NPPclim.nc and fit.piqs.rda


##########################
# Extra Notes
##########################


  - Climatology is used when "real" years are not available. This includes 2000-2002  and could also be used in 2024 onwards,
  whenever real data run out.

  - Diurnalize is described in the CT documentation. It needs monthly Rh and NPP, from which it generates temporally-downscaled
    GPP and total respiration.

  - Fire and (bio)fuel emissions are taken from the daily files provided by MiCASA. This is the only source for fairly 
    high-resolution-in-time emissions for those processes.

  - MiCASA provides temporally-downscaled (non-fire and -fuel) fluxes, using a method similar to ours, but with different meteorology
    (NASA's MERRA2 reanalysis). Our method is a little bit better due to the PIQS part of the scheme, which smooths out abrupt changes
    at monthly boundaries, and our meteo comes from ERA5. That makes the downscaling consistent with the atmospheric transport provided
    by TM5. So, we do not use the MiCASA temporally-downscaled fluxes; we start with the monthlies and apply the downscaling ourselves. 

https://nco.sourceforge.net/nco.pdf

ncea (netCDF Ensemble Average)
    performs gridpoint averages of variables across an arbitrary number (an ensemble) 
    of input files, with each file receiving an equal weight in the average

    -O, overwrite output if it exists

    Note: ncea is deprecated for nces (netCDF Ensemble Statistics)

ncks (netCDF Kitchen Sink)
    extracts (a subset of the) data from input-file, regrids it according to map-file if specified,
    then writes in netCDF format to output-file, and optionally writes it in flat binary format to fl_bnr,
    and optionally prints it to screen.


VSEM-ET (possibly unrelated, but nice picture)
    https://insightmaker.com/insight/6DkHwGgVTkedbnCviUX8bD/Clone-of-Very-Simple-Ecosystem-Model-with-Evapotranspiration-VSEM-ET


Original daily header dump:
netcdf MiCASA_v1_flux_x3600_y1800_daily_20130507 {
dimensions:
	lat = 1800 ;
	lon = 3600 ;
	time = UNLIMITED ; // (1 currently)
	nv = 2 ;
variables:
	double lat(lat) ;
		lat:units = "degrees_north" ;
		lat:long_name = "latitude" ;
	double lon(lon) ;
		lon:units = "degrees_east" ;
		lon:long_name = "longitude" ;
	double time(time) ;
		time:units = "days since 1980-01-01" ;
		time:long_name = "time" ;
		time:bounds = "time_bnds" ;
	double time_bnds(time, nv) ;
		time_bnds:units = "days since 1980-01-01" ;
		time_bnds:long_name = "time bounds" ;
	float NPP(time, lat, lon) ;
		NPP:units = "kg m-2 s-1" ;
		NPP:expressed_as = "carbon" ;
		NPP:long_name = "Net primary productivity" ;
	float Rh(time, lat, lon) ;
		Rh:units = "kg m-2 s-1" ;
		Rh:expressed_as = "carbon" ;
		Rh:long_name = "Heterotrophic respiration" ;
	float FIRE(time, lat, lon) ;
		FIRE:units = "kg m-2 s-1" ;
		FIRE:expressed_as = "carbon" ;
		FIRE:long_name = "Fire emission" ;
	float FUEL(time, lat, lon) ;
		FUEL:units = "kg m-2 s-1" ;
		FUEL:expressed_as = "carbon" ;
		FUEL:long_name = "Fuel wood emission" ;
	float ATMC(time, lat, lon) ;
		ATMC:units = "kg m-2 s-1" ;
		ATMC:expressed_as = "carbon" ;
		ATMC:long_name = "Atmospheric correction" ;
	float NEE(time, lat, lon) ;
		NEE:units = "kg m-2 s-1" ;
		NEE:expressed_as = "carbon" ;
		NEE:long_name = "Net ecosystem exchange" ;

// global attributes:
		:Conventions = "CF-1.9" ;
		:contact = "Brad Weir <brad.weir@nasa.gov>" ;
		:institution = "NASA Goddard Space Flight Center" ;
		:title = "MiCASA Daily NPP Rh ATMC NEE FIRE FUEL Fluxes 0.1 degree x 0.1 degree v1" ;
		:LongName = "MiCASA Daily NPP Rh ATMC NEE FIRE FUEL Fluxes 0.1 degree x 0.1 degree" ;
		:ShortName = "MICASA_FLUX_D" ;
		:VersionID = "1" ;
		:GranuleID = "MiCASA_v1_flux_x3600_y1800_daily_20130507.nc4" ;
		:Format = "netCDF" ;
		:ProcessingLevel = "4" ;
		:IdentifierProductDOIAuthority = "https://doi.org/" ;
		:IdentifierProductDOI = "10.5067/ZBXSA1LEN453" ;
		:ReadMeURL = "https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA/v1/MiCASA_README.pdf" ;
		:RangeBeginningDate = "2013-05-07" ;
		:RangeBeginningTime = "00:00:00.000000" ;
		:RangeEndingDate = "2013-05-07" ;
		:RangeEndingTime = "23:59:59.999999" ;
		:NorthernmostLatiude = "90.0" ;
		:WesternmostLongitude = "-180.0" ;
		:SouthernmostLatitude = "-90.0" ;
		:EasternmostLongitude = "180.0" ;
		:comment = "Positive NPP indicates uptake by vegetation. Positive Rh indicates emission to the atmosphere. NEE = Rh - NPP - ATMC, and NBE = NEE + FIRE + FUEL. ATMC adjusts net exchange to account for missing processes and better match long-term atmospheric budgets." ;
		:ProductionDateTime = "2024-09-23T01:58:21Z" ;
}

