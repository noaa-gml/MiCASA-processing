ct.setup()
regs <- load.ncdf(sprintf('%s/tools/shared/aux/regions.nc',CARBONTRACKER),vars="grid_cell_area")

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

gca <- matrix(NA,nrow=360,ncol=180)


for (irow in 1:360) {
  lon.rad <- (pi/180)*(regs$longitude[irow]+c(-0.5,0.5))
  for (icol in 1:180) {
    lat.rad <- (pi/180)*(regs$latitude[icol]+c(-0.5,0.5))
    gca[irow,icol] <- archimedes(lon.rad,lat.rad)
  }
}


print(summary(as.vector((gca-regs$grid_cell_area)/regs$grid_cell_area))

#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#1.749e-06 2.206e-06 2.274e-06 2.271e-06 2.352e-06 2.710e-06
#
# Relative error of 2.7e-6, consistent with single-precision numerics

      
      
