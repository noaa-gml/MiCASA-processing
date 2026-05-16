#!/usr/bin/env Rscript
## Unit tests for lib/era5_meteo.r (base R only, CI-runnable).
##
## Run:  Rscript tests/test_era5_meteo.r
## Exits non-zero on any failure.

## Locate lib/era5_meteo.r relative to this script, so the test runs
## from any working directory.
.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
source(file.path(.dir, "..", "lib", "era5_meteo.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}

## ---- era5.relpath ---------------------------------------------------------
tmpl <- "YYYY/MM/VVV_YYYYMMDD_00p01.nc"
check("era5.relpath substitutes YYYY/MM/DD/VVV",
      era5.relpath(tmpl, 2026, 2, 15, "ssrd") == "2026/02/ssrd_20260215_00p01.nc")
check("era5.relpath zero-pads month and day",
      era5.relpath(tmpl, 2001, 1, 5, "t2m") == "2001/01/t2m_20010105_00p01.nc")

## ---- encode.day.runs ------------------------------------------------------
check("encode.day.runs: empty -> ''",
      encode.day.runs(integer(0), character(0)) == "")
## srcvec is indexed positionally by day number, so day 5 needs a
## length->=5 vector with position 5 set.
src.one <- rep(NA_character_, 5L); src.one[5L] <- "primary"
check("encode.day.runs: single day",
      encode.day.runs(5L, src.one) == "primary:5")
src31 <- setNames(rep("primary", 31), as.character(1:31))
check("encode.day.runs: full contiguous run",
      encode.day.runs(1:31, src31) == "primary:1-31")
## genuine mixed month: days 1-30 primary, day 31 fasttrack
srcmix <- setNames(c(rep("primary", 30), "fasttrack"), as.character(1:31))
check("encode.day.runs: two sources",
      encode.day.runs(1:31, srcmix) == "primary:1-30 fasttrack:31")
## non-contiguous run within one source
srcgap <- setNames(rep("primary", 7), as.character(1:7))
check("encode.day.runs: gap -> comma-separated runs",
      encode.day.runs(c(1L, 2L, 3L, 5L, 6L, 7L), srcgap) == "primary:1-3,5-7")

## ---- resolve.era5.source (synthetic meteo trees) --------------------------
root  <- file.path(tempdir(), paste0("era5test_", Sys.getpid()))
prim  <- file.path(root, "primary")
fast  <- file.path(root, "fasttrack")
vn    <- c("t2m", "ssrd", "stl1", "swvl1")
mkday <- function(treedir, yr, mon, day, vars = vn) {
  d <- file.path(treedir, sprintf("%04d/%02d", yr, mon))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (v in vars)
    file.create(file.path(d, sprintf("%s_%04d%02d%02d_00p01.nc", v, yr, mon, day)))
}
mkday(prim, 2025, 6, 15)                  # primary has all 4 vars
mkday(fast, 2026, 2, 15)                  # only fasttrack has this day
mkday(prim, 2026, 3, 1, vars = vn[1:3])   # primary missing swvl1
mkday(fast, 2026, 3, 1)                   # fasttrack complete
era5dirs <- c(primary = prim, fasttrack = fast)

check("resolve: primary preferred when it has the day",
      identical(resolve.era5.source(era5dirs, tmpl, 2025, 6, 15, vn), "primary"))
check("resolve: falls back to fasttrack",
      identical(resolve.era5.source(era5dirs, tmpl, 2026, 2, 15, vn), "fasttrack"))
check("resolve: incomplete primary -> fasttrack",
      identical(resolve.era5.source(era5dirs, tmpl, 2026, 3, 1, vn), "fasttrack"))
check("resolve: NA when no tree has the day",
      is.na(resolve.era5.source(era5dirs, tmpl, 2099, 1, 1, vn)))

unlink(root, recursive = TRUE)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall era5_meteo tests passed\n")
