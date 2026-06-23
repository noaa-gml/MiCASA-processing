#!/usr/bin/env Rscript

# diag_v1_vNRT_handoff.r
#
# ARCHIVED (lives in archive/). A one-off splice-continuity diagnostic for the
# v1→vNRT (2024-12 / 2025-01) handoff. Paths are cwd-relative, so it still runs
# from the working directory regardless of the script's location.
#
# Proposal #3 in README.ash: sanity-check that splicing MiCASA v1 (final,
# through 2024-12) and vNRT (2025-01 onward) into a single monthly record
# does not introduce a step in absolute scale that PIQS would smooth across.
#
# Reads the existing spliced monthly file
#     monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc
# and reports per-month area-weighted global totals of NPP and Rh (and ATMC
# if present) around the Dec 2024 / Jan 2025 boundary. Saves a PDF plot and a
# CSV table to the working directory.
#
# Run on Orion from the MiCASA working directory:
#     Rscript archive/diag_v1_vNRT_handoff.r
#
# Optional env vars:
#   MICASA_DIAG_FROM=YYYYMM   start month for the plot window (default 202301)
#   MICASA_DIAG_TO=YYYYMM     end   month for the plot window (default 202612)
#   MICASA_DIAG_BOUNDARY=YYYYMM
#                             month of the v1 -> vNRT handoff (default 202501,
#                             i.e. the first vNRT month). The plotted vertical
#                             line falls between BOUNDARY-1 and BOUNDARY.

ct.setup()

monthly.nc <- "monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc"
if(!file.exists(monthly.nc)) {
  stop(sprintf("Cannot find %s -- run from the MiCASA working directory.", monthly.nc))
}

cat(sprintf("Reading %s ...\n", monthly.nc))
din <- load.ncdf(monthly.nc)

# Grid-cell area (m^2). regions.nc carries it precomputed; if unavailable we
# fall back to test_gca.r's analytic formula.
gca <- NULL
regions.path <- sprintf("%s/tools/shared/aux/regions.nc", CARBONTRACKER)
if(file.exists(regions.path)) {
  regs <- load.ncdf(regions.path, vars="grid_cell_area")
  gca <- regs$grid_cell_area
}
if(is.null(gca) || any(dim(gca) != c(360, 180))) {
  cat("Falling back to analytic grid-cell area...\n")
  R <- 6371007.2
  gca <- matrix(NA, 360, 180)
  for(irow in 1:360) {
    lon.rad <- (pi/180) * (din$longitude[irow] + c(-0.5, 0.5))
    for(icol in 1:180) {
      lat.rad <- (pi/180) * (din$latitude[icol] + c(-0.5, 0.5))
      gca[irow, icol] <- (sin(lat.rad[2]) - sin(lat.rad[1])) * (lon.rad[2] - lon.rad[1]) * R^2
    }
  }
}

# Convert kg C m-2 s-1 monthly mean -> Pg C / month total.
# Days in each calendar month (no leap-year handling -- 28 used for Feb is
# fine for the size of step we are looking for).
times.lt <- as.POSIXlt(din$time)
yrs  <- times.lt$year + 1900
mons <- times.lt$mon  + 1
dpm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
sec.per.month <- dpm[mons] * 86400

flux.total <- function(field) {
  # field has dim (lon, lat, time); gca has dim (lon, lat).
  # Result: per-time global total in Pg C / month.
  apply(field, 3, function(slab) sum(slab * gca, na.rm=TRUE)) * sec.per.month / 1e12
}

vars.present <- intersect(c("NPP", "Rh", "ATMC", "FIRE", "FUEL"), names(din))
cat(sprintf("Variables present in monthly file: %s\n", paste(vars.present, collapse=", ")))
totals <- list()
for(v in vars.present) totals[[v]] <- flux.total(din[[v]])

# Window selection.
parse.yyyymm <- function(s, default.year, default.month) {
  if(nchar(s) == 0) return(c(default.year, default.month))
  s <- as.integer(s)
  c(s %/% 100, s %% 100)
}
window.from <- parse.yyyymm(Sys.getenv("MICASA_DIAG_FROM"), 2023, 1)
window.to   <- parse.yyyymm(Sys.getenv("MICASA_DIAG_TO"),   2026, 12)
boundary    <- parse.yyyymm(Sys.getenv("MICASA_DIAG_BOUNDARY"), 2025, 1)

mask <- (yrs > window.from[1] | (yrs == window.from[1] & mons >= window.from[2])) &
        (yrs < window.to[1]   | (yrs == window.to[1]   & mons <= window.to[2]))
ix <- which(mask)
if(length(ix) == 0) stop("No months fall inside the requested plot window.")

# CSV.
out.csv <- "diag_v1_vNRT_handoff.csv"
hdr <- c("year", "month", names(totals))
csv.df <- data.frame(year=yrs[ix], month=mons[ix])
for(v in names(totals)) csv.df[[v]] <- totals[[v]][ix]
write.csv(csv.df, out.csv, row.names=FALSE)
cat(sprintf("Wrote %s (%d months)\n", out.csv, nrow(csv.df)))

# Print the values surrounding the boundary, so the size of any step is
# immediately obvious without opening the PDF.
boundary.t <- ISOdatetime(boundary[1], boundary[2], 1, 0, 0, 0, tz="UTC")
near <- which(abs(as.numeric(din$time[ix]) - as.numeric(boundary.t)) < 92*86400)
cat(sprintf("\nMonths near the v1 -> vNRT handoff (boundary at %d-%02d):\n",
            boundary[1], boundary[2]))
print(csv.df[near, ], row.names=FALSE)

# Plot.
out.pdf <- "diag_v1_vNRT_handoff.pdf"
pdf(out.pdf, width=10, height=3 * length(totals))
par(mfrow=c(length(totals), 1), mar=c(3, 4, 2, 1), las=1)
plot.times <- din$time[ix]
boundary.x <- boundary.t - 15.5 * 86400  # midpoint between boundary-1 and boundary

for(v in names(totals)) {
  yy <- totals[[v]][ix]
  plot(plot.times, yy, type="o", pch=20,
       xlab="", ylab=sprintf("%s (Pg C / month)", v),
       main=sprintf("MiCASA monthly global %s -- v1 -> vNRT handoff", v))
  abline(v=boundary.x, col="red", lty=2)
  text(boundary.x, par("usr")[4], "v1 | vNRT", adj=c(-0.05, 1.2), col="red", cex=0.8)
  grid()
}
dev.off()
cat(sprintf("Wrote %s\n", out.pdf))
