#!/usr/bin/env Rscript
## Self-contained tests for aggregate.to.1x1 in lib/ingest_common.r.
##
## Doesn't require ct.setup() or netCDF tooling — sources ingest_common.r
## (whose ncdim/ncvar wrappers are lazy and never invoked here) and
## exercises the geometry helpers on synthetic fields.
##
## Usage: Rscript lib/test_aggregate.r

work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "lib", "ingest_common.r"))

cat("--- aggregate.to.1x1 tests ---\n")

## Build the 1800-element latitude-area vector on the real MiCASA 0.1° grid.
lats <- seq(-89.95, 89.95, 0.1)
gca  <- compute.gca(lats)

pass <- function(label, ok, info = "") {
  cat(sprintf("[%s] %-32s %s%s\n",
              if (ok) "PASS" else "FAIL",
              label, info, if (!ok) " <<<" else ""))
  invisible(ok)
}

## --- Test 1: constant field aggregates to itself ----------------------------
fld <- matrix(7.5, 3600, 1800)
out <- aggregate.to.1x1(fld, gca)
pass("constant field",
     max(abs(out - 7.5)) < 1e-10,
     sprintf("max|err|=%.2e", max(abs(out - 7.5))))

## --- Test 2: lon-only field aggregates per-lon-bin mean (1..3600) -----------
fld <- matrix(rep(1:3600, 1800), nrow = 3600, ncol = 1800)
out <- aggregate.to.1x1(fld, gca)
expected <- 10 * (1:360) - 4.5
err <- max(abs(out - expected))
pass("lon-only field", err < 1e-9, sprintf("max|err|=%.2e", err))

## --- Test 3: lat-only sin(lat) field gets correct area-weighted mean --------
lat.rad <- lats * pi / 180
fld <- matrix(rep(sin(lat.rad), each = 3600), nrow = 3600, ncol = 1800)
out <- aggregate.to.1x1(fld, gca)

expected <- rep(NA_real_, 180)
for (jlat in 1:180) {
  inlats <- 1:10 + 10 * (jlat - 1)
  expected[jlat] <- sum(sin(lat.rad[inlats]) * gca[inlats]) / sum(gca[inlats])
}
err <- max(abs(out[1, ] - expected))
pass("sin(lat) field — area-weighted",
     err < 1e-12, sprintf("max|err|=%.2e", err))

## --- Test 4: regression check — old buggy formula must NOT pass test 3 -----
old.aggregate.to.1x1 <- function(fld, gca) {
  retval <- matrix(0, 360, 180)
  for (jlat in 1:180) {
    inlats <- 1:10 + 10 * (jlat - 1)
    for (ilon in 1:360) {
      inlons <- 1:10 + 10 * (ilon - 1)
      retval[ilon, jlat] <- 0
      for (inlon in inlons) {
        retval[ilon, jlat] <- retval[ilon, jlat] +
          weighted.mean(fld[inlons, inlats], weights = gca[inlats], na.rm = TRUE)
      }
      retval[ilon, jlat] <- retval[ilon, jlat] / 10
    }
  }
  retval
}
old.out <- old.aggregate.to.1x1(fld, gca)
old.err.abs <- max(abs(old.out[1, ] - expected))
old.err.rel <- max(abs((old.out[1, ] - expected) /
                       pmax(abs(expected), 1e-9)))
pass("regression: old formula is wrong",
     old.err.abs > 1e-6,
     sprintf("max|err|=%.2e (rel %.2e)", old.err.abs, old.err.rel))

cat("--- done ---\n")
