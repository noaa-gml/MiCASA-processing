#!/usr/bin/env Rscript
## Self-contained tests + benchmark for aggregate.to.1x1.
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
  cat(sprintf("[%s] %-44s %s%s\n",
              if (ok) "PASS" else "FAIL",
              label, info, if (!ok) " <<<" else ""))
  invisible(ok)
}

## ---------------------------------------------------------------------------
## Reference scalar implementation, kept here so the vectorized version in
## lib/ingest_common.r has something to compare against. This is the
## "correct, slow" formulation — bug-fixed in 2026-04-26 to weight lat-area
## along the latitude axis (not longitude) in the inner weighted mean.
## ---------------------------------------------------------------------------
aggregate.to.1x1.scalar <- function(fld, gca) {
  retval <- matrix(0, 360, 180)
  for (jlat in 1:180) {
    inlats <- 1:10 + 10 * (jlat - 1)
    w <- rep(gca[inlats], each = 10)
    for (ilon in 1:360) {
      inlons <- 1:10 + 10 * (ilon - 1)
      retval[ilon, jlat] <- weighted.mean(
        as.vector(fld[inlons, inlats]), w = w, na.rm = TRUE)
    }
  }
  retval
}

## ---------------------------------------------------------------------------
## Synthetic-field correctness tests
## ---------------------------------------------------------------------------

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

## --- Test 4: regression — pre-fix buggy formula must NOT match expected ----
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

## --- Test 5: vectorized matches scalar reference on random no-NA field ------
set.seed(42)
fld <- matrix(rnorm(3600 * 1800), 3600, 1800)
out_v <- aggregate.to.1x1(fld, gca)
out_s <- aggregate.to.1x1.scalar(fld, gca)
err <- max(abs(out_v - out_s))
pass("vectorized vs scalar — random field",
     err < 1e-10, sprintf("max|err|=%.2e", err))

## --- Test 6: vectorized matches scalar with sparse NAs ----------------------
set.seed(43)
fld_na <- fld
fld_na[sample(length(fld_na), length(fld_na) %/% 100)] <- NA  # ~1% NA
out_v <- aggregate.to.1x1(fld_na, gca)
out_s <- aggregate.to.1x1.scalar(fld_na, gca)
err <- max(abs(out_v - out_s), na.rm = TRUE)
pass("vectorized vs scalar — 1% NA",
     err < 1e-10, sprintf("max|err|=%.2e", err))

## --- Test 7: all-NA block becomes NaN, not 0 or NA --------------------------
fld_blk <- matrix(1, 3600, 1800)
fld_blk[1:10, 1:10] <- NA   # entire (lon-block 1, lat-block 1) is NA
out_v <- aggregate.to.1x1(fld_blk, gca)
out_s <- aggregate.to.1x1.scalar(fld_blk, gca)
pass("all-NA block → NaN, scalar agrees",
     is.nan(out_v[1, 1]) && is.nan(out_s[1, 1]) &&
       max(abs(out_v[-1, -1] - 1)) < 1e-12,
     sprintf("v[1,1]=%s s[1,1]=%s", out_v[1, 1], out_s[1, 1]))

## ---------------------------------------------------------------------------
## Benchmark: vectorized vs scalar on a realistic-size random field.
## ingest_byyear.r calls aggregate.to.1x1 ~1460 times per year (4 tracers
## × 365 days), so the per-call speedup compounds.
## ---------------------------------------------------------------------------

cat("\n--- benchmark (random 3600x1800 field, 1% NA) ---\n")
set.seed(44)
bench_fld <- matrix(rnorm(3600 * 1800), 3600, 1800)
bench_fld[sample(length(bench_fld), length(bench_fld) %/% 100)] <- NA

n_v <- 5
t_v <- system.time({ for (i in 1:n_v) out_v <- aggregate.to.1x1(bench_fld, gca) })
cat(sprintf("vectorized: %.3f s/call (mean of %d)\n", t_v["elapsed"] / n_v, n_v))

n_s <- 1
t_s <- system.time({ out_s <- aggregate.to.1x1.scalar(bench_fld, gca) })
cat(sprintf("scalar:     %.3f s/call (n=%d)\n", t_s["elapsed"] / n_s, n_s))

speedup <- (t_s["elapsed"] / n_s) / (t_v["elapsed"] / n_v)
cat(sprintf("speedup:    %.1fx\n", speedup))

## Project to ingest_byyear.r wall time.
calls_per_year <- 4 * 365
cat(sprintf("\nProjected per-year aggregate-only cost (ignoring I/O):\n"))
cat(sprintf("  scalar:     %.0f s = %.1f h\n",
            calls_per_year * (t_s["elapsed"] / n_s),
            calls_per_year * (t_s["elapsed"] / n_s) / 3600))
cat(sprintf("  vectorized: %.0f s = %.1f min\n",
            calls_per_year * (t_v["elapsed"] / n_v),
            calls_per_year * (t_v["elapsed"] / n_v) / 60))

cat("\n--- done ---\n")
