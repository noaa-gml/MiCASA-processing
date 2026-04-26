#!/usr/bin/env Rscript
## test_mtime_skip.r
##
## Unit test for out.is.fresh() — the make-style skip helper used by
## ingest_byyear.r and ingest_monthly.r.

ct.setup()
work.dir <- Sys.getenv("WORK_DIR", getwd())
source(file.path(work.dir, "lib", "ingest_common.r"))

td  <- tempfile("mtime_test_"); dir.create(td)
src <- file.path(td, "src.nc")
out <- file.path(td, "out.nc")

ok <- TRUE
fail <- function(msg) { cat(sprintf("FAIL: %s\n", msg)); ok <<- FALSE }

# -- case 1: output absent -> not fresh
file.create(src); Sys.sleep(0.05)
if (out.is.fresh(out, src)) fail("missing output considered fresh")

# -- case 2: output written after source -> fresh
file.create(out); Sys.sleep(0.05)
if (!out.is.fresh(out, src)) fail("output newer than source not considered fresh")

# -- case 3: source touched after output -> stale, must re-ingest
Sys.sleep(0.05); Sys.setFileTime(src, Sys.time())
if (out.is.fresh(out, src)) fail("stale output incorrectly considered fresh")

# -- case 4: explicit equal mtime -> not fresh (we use strict >)
Sys.setFileTime(out, file.mtime(src))
if (out.is.fresh(out, src)) fail("equal mtimes should not count as fresh")

# -- case 5: output written after source again -> fresh
Sys.sleep(0.05); Sys.setFileTime(out, Sys.time())
if (!out.is.fresh(out, src)) fail("re-touched output not considered fresh")

unlink(td, recursive = TRUE)
cat(sprintf("%s\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
