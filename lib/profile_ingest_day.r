#!/usr/bin/env Rscript
## profile_ingest_day.r — break down per-day cost of ingest_byyear into
## (load.ncdf full file)  vs  (load only 4 needed vars)  vs  aggregate  vs  write.

ct.setup()
work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "config.r"))
source(file.path(work.dir, "lib", "ingest_common.r"))
cfg <- micasa.config()

library(ncdf4)

src <- micasa.raw.daily(cfg, 2024, 1, 1)
cat(sprintf("source: %s  (%.1f MB on disk)\n", src, file.info(src)$size / 1024^2))

probe <- load.ncdf(src)
gca   <- compute.gca(probe$lat)
rm(probe); gc(verbose = FALSE)

dim.lon <- micasa.dim.lon()
dim.lat <- micasa.dim.lat()
dim.time <- micasa.time.dim(ISOdatetime(2024, 1, 1, 12, 0, 0, tz = "UTC"))

timeit <- function(label, expr) {
  t0 <- proc.time()[3]
  r  <- eval(expr, envir = parent.frame())
  t1 <- proc.time()[3]
  cat(sprintf("  %-35s %6.2fs\n", label, t1 - t0))
  attr(r, ".elapsed") <- t1 - t0
  r
}

n <- 5
cat(sprintf("\nProfiling %d iterations on the same daily file:\n", n))

for (iter in seq_len(n)) {
  cat(sprintf("\n-- iter %d --\n", iter))

  # (A) full-file load via ct.setup helper (current behavior)
  t0 <- proc.time()[3]; ncin_full <- load.ncdf(src); t_full <- proc.time()[3] - t0
  cat(sprintf("  %-35s %6.2fs\n", "load.ncdf full (current)", t_full))

  # (B) read only 4 needed vars via raw ncdf4
  t0 <- proc.time()[3]
  ncf <- nc_open(src)
  vars4 <- list()
  for (nm in micasa.tracers) vars4[[nm]] <- ncvar_get(ncf, nm)
  nc_close(ncf)
  t_4vars <- proc.time()[3] - t0
  cat(sprintf("  %-35s %6.2fs   (saves %.1fs)\n",
              "ncvar_get 4 needed only", t_4vars, t_full - t_4vars))

  # (C) aggregate.to.1x1 ×4 (vectorized form)
  t0 <- proc.time()[3]
  vals <- list()
  for (nm in micasa.tracers) vals[[nm]] <- aggregate.to.1x1(vars4[[nm]], gca) * 1e3
  t_agg <- proc.time()[3] - t0
  cat(sprintf("  %-35s %6.2fs\n", "aggregate.to.1x1 x4", t_agg))

  # (D) write 1° output (level 9, current)
  ncout <- tempfile(fileext = ".nc")
  ncin_for_attrs <- ncin_full  # for long_name passthrough
  vars_def <- make.tracer.vars(ncin_for_attrs, dim.lon, dim.lat, dim.time)
  t0 <- proc.time()[3]
  write.netcdf(ncout, vars_def, vals, src, "profile_ingest_day.r")
  t_write <- proc.time()[3] - t0
  cat(sprintf("  %-35s %6.2fs   (file %.1f KB)\n",
              "write.netcdf level 9", t_write,
              file.info(ncout)$size / 1024))
  file.remove(ncout)

  cat(sprintf("  TOTAL (current path): %.2fs   |   if 4-var read: %.2fs (save %.0f%%)\n",
              t_full + t_agg + t_write,
              t_4vars + t_agg + t_write,
              100 * (t_full - t_4vars) / (t_full + t_agg + t_write)))

  rm(ncin_full, vars4, vals); gc(verbose = FALSE)
}
