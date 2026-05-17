#!/usr/bin/env Rscript
## Unit tests for lib/manifest.r :: manifest.record (base R only, CI-runnable).
##
## manifest.record appends structured rows to jobs/run_manifest.tsv; verify_v2
## reads that manifest instead of globbing job logs. These checks pin its
## output format and its "never error the caller" guarantee, using a temp
## WORK_DIR so nothing real is touched.
##
## Run:  Rscript tests/test_manifest.r
## Exits non-zero on any failure.

.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "manifest.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}

tmp  <- file.path(tempdir(), paste0("mtest_", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
path <- file.path(tmp, "jobs", "run_manifest.tsv")

## ---- first record creates the file with a header ------------------------
manifest.record("step-a", "start", detail = "hello", work.dir = tmp)
check("manifest file is created under jobs/", file.exists(path))
lines <- readLines(path)
check("first line is the 7-column header",
      startsWith(lines[1],
                 "# timestamp\tstep\tstatus\thost\tcommit\telapsed_s\tdetail"))
check("one data row after one record", length(lines) == 2L)

f <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
check("data row has 7 tab-separated columns", length(f) == 7L)
check("timestamp is ISO-8601 UTC",
      grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", f[1]))
check("step column is the step", identical(f[2], "step-a"))
check("status column is the status", identical(f[3], "start"))
check("elapsed column is '-' when not given", identical(f[6], "-"))
check("detail column is the detail", identical(f[7], "hello"))

## ---- second record appends (header not repeated) ------------------------
manifest.record("step-b", "ok", elapsed = 42, detail = "done", work.dir = tmp)
lines <- readLines(path)
check("second record appends a row", length(lines) == 3L)
f2 <- strsplit(lines[3], "\t", fixed = TRUE)[[1]]
check("elapsed column carries the integer", identical(f2[6], "42"))
check("status column carries 'ok'", identical(f2[3], "ok"))

## ---- tab / newline in detail are squashed to spaces ---------------------
manifest.record("step-c", "info", detail = "a\tb\nc", work.dir = tmp)
lines <- readLines(path)
check("one row per record -- embedded newline did not add a line",
      length(lines) == 4L)
f3 <- strsplit(lines[4], "\t", fixed = TRUE)[[1]]
check("tabs/newlines in detail squashed to spaces, row still 7 columns",
      length(f3) == 7L && f3[7] == "a b c")

## ---- never errors the caller, even when the manifest cannot be written --
ok <- tryCatch({
  suppressWarnings(manifest.record("x", "ok", work.dir = "/nonexistent/zz/qq"))
  TRUE
}, error = function(e) FALSE)
check("manifest.record never errors its caller", isTRUE(ok))

unlink(tmp, recursive = TRUE)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall manifest tests passed\n")
