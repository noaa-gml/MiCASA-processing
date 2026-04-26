#!/usr/bin/env Rscript
## test_ingest_bitident.r
##
## Confirm that switching load.ncdf(srcnm) -> load.ncdf(srcnm, vars=micasa.tracers)
## produces bit-identical aggregated output for one daily file.

ct.setup()
work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
source(file.path(work.dir, "lib", "ingest_common.r"))
cfg <- micasa.config()

src <- micasa.raw.daily(cfg, 2024, 1, 1)
cat(sprintf("source: %s\n", src))

probe <- load.ncdf(src, quiet = TRUE)
gca   <- compute.gca(probe$lat)

ncin_full <- load.ncdf(src, quiet = TRUE)                           # old path
ncin_4    <- load.ncdf(src, vars = micasa.tracers, quiet = TRUE)    # new path

ok <- TRUE
for (nm in micasa.tracers) {
  a <- ncin_full[[nm]]
  b <- ncin_4[[nm]]
  if (!identical(dim(a), dim(b))) {
    cat(sprintf("FAIL %s: dim mismatch\n", nm)); ok <- FALSE; next
  }
  diff_max <- max(abs(a - b), na.rm = TRUE)
  na_match <- identical(is.na(a), is.na(b))
  cat(sprintf("  %-5s  max|a-b|=%.3e  NA-pattern-match=%s  identical=%s\n",
              nm, diff_max, na_match, identical(a, b)))
  if (!identical(a, b)) ok <- FALSE
}

# also compare aggregated outputs end-to-end
cat("\nAggregated outputs:\n")
for (nm in micasa.tracers) {
  va <- aggregate.to.1x1(ncin_full[[nm]], gca) * 1e3
  vb <- aggregate.to.1x1(ncin_4[[nm]],    gca) * 1e3
  diff_max <- max(abs(va - vb), na.rm = TRUE)
  cat(sprintf("  %-5s  max|agg(a)-agg(b)|=%.3e  identical=%s\n",
              nm, diff_max, identical(va, vb)))
  if (!identical(va, vb)) ok <- FALSE
}

# attribute pass-through (used by make.tracer.vars for longname)
cat("\nlong_name attributes:\n")
for (nm in micasa.tracers) {
  la <- attributes(ncin_full[[nm]])$long_name
  lb <- attributes(ncin_4[[nm]])$long_name
  cat(sprintf("  %-5s  full='%s'  4var='%s'  match=%s\n",
              nm, la, lb, identical(la, lb)))
  if (!identical(la, lb)) ok <- FALSE
}

cat(sprintf("\n%s\n", if (ok) "PASS: bit-identical" else "FAIL"))
if (!ok) quit(status = 1)
