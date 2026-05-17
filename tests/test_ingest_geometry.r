#!/usr/bin/env Rscript
## Unit tests for lib/ingest_common.r geometry (base R only, CI-runnable).
##
## archimedes() and compute.gca() compute spherical grid-cell areas -- the
## weights the 0.1-degree -> 1-degree aggregation (aggregate.to.1x1) depends
## on. They are pure base-R functions; these checks pin them against
## analytic ground truth (a sphere has area 4*pi*R^2). out.is.fresh() is the
## mtime-based skip-existing gate sourced from the same file.
##
## Run:  Rscript tests/test_ingest_geometry.r
## Exits non-zero on any failure.

.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "ingest_common.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}
close.rel <- function(a, b, rtol = 1e-9) abs(a - b) <= rtol * abs(b)
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

SPH <- 4 * pi * EARTH_RADIUS_M^2          # area of the whole sphere

## ---- archimedes: analytic areas ------------------------------------------
check("whole sphere = 4*pi*R^2",
      close.rel(archimedes(c(-pi, pi), c(-pi/2, pi/2)), SPH))
check("northern hemisphere = half the sphere",
      close.rel(archimedes(c(-pi, pi), c(0, pi/2)), SPH / 2))
check("a 180-degree lon wedge, pole to pole = half the sphere",
      close.rel(archimedes(c(0, pi), c(-pi/2, pi/2)), SPH / 2))
check("a cell has positive area", archimedes(c(0, 1), c(0, 1)) > 0)

d <- 0.1
check("equator-symmetric bands have equal area",
      close.rel(archimedes(c(0, d), c(0.5, 0.5 + d)),
                archimedes(c(0, d), c(-0.5 - d, -0.5))))
check("an equatorial band is larger than an equal-width polar band",
      archimedes(c(0, d), c(0, d)) > archimedes(c(0, d), c(pi/2 - d, pi/2)))
check("splitting a cell in latitude conserves area",
      close.rel(archimedes(c(0, d), c(0.2, 0.6)),
                archimedes(c(0, d), c(0.2, 0.4)) +
                archimedes(c(0, d), c(0.4, 0.6))))

## ---- archimedes: input validation ----------------------------------------
check("rejects a lons vector not of length 2", errs(archimedes(c(0), c(0, 1))))
check("rejects a lats vector not of length 2", errs(archimedes(c(0, 1), c(0))))
check("rejects |lon| > pi",   errs(archimedes(c(-4, 4), c(0, 1))))
check("rejects |lat| > pi/2", errs(archimedes(c(0, 1), c(-2, 2))))

## ---- compute.gca: the 0.1-degree latitude cell-area vector ---------------
lats <- seq(-89.95, 89.95, length.out = 1800)   # 0.1-deg cell centres
gca  <- compute.gca(lats)
check("compute.gca returns one area per latitude", length(gca) == 1800L)
check("every cell area is positive", all(gca > 0))
check("the 0.1-deg grid tiles the sphere (sum x 3600 lon cells = 4*pi*R^2)",
      close.rel(sum(gca) * 3600, SPH, rtol = 1e-6))
check("cell areas are symmetric about the equator",
      max(abs(gca - rev(gca))) <= 1e-6 * max(gca))
check("the equatorial cell is the largest", which.max(gca) %in% c(900L, 901L))
check("a polar cell is far smaller than the equatorial cell",
      gca[1] < 0.01 * max(gca))

## ---- out.is.fresh: mtime-based skip-existing gate ------------------------
tmp <- file.path(tempdir(), paste0("fresh_", Sys.getpid()))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
src <- file.path(tmp, "src.nc")
out <- file.path(tmp, "out.nc")
writeLines("src", src)
writeLines("out", out)
Sys.setFileTime(src, Sys.time() - 100)          # source older than output
Sys.setFileTime(out, Sys.time())
check("fresh when the output is newer than the source",
      isTRUE(out.is.fresh(out, src)))
Sys.setFileTime(src, Sys.time())                # source now newer
Sys.setFileTime(out, Sys.time() - 100)
check("stale when the source is newer than the output",
      isFALSE(out.is.fresh(out, src)))
check("not fresh when the output is missing",
      isFALSE(out.is.fresh(file.path(tmp, "nonexistent.nc"), src)))
unlink(tmp, recursive = TRUE)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall ingest geometry tests passed\n")
