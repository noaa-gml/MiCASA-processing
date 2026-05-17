## lib/manifest.r -- append structured run records to jobs/run_manifest.tsv.
##
##   manifest.record(step, status = "info", elapsed = NA, detail = "")
##
##     step       producing script / stage, e.g. "diurnalize-ERA5.r"
##     status     start | ok | fail | info
##     elapsed    integer seconds, or NA when not applicable
##     detail     free text (tabs / newlines are squashed to spaces)
##
## The manifest is the pipeline's structured run record -- verify_v2 reads it
## instead of globbing job logs. The file is tab-separated; the columns are
##     timestamp  step  status  host  commit  elapsed_s  detail
##
## manifest.record never errors out its caller: a logging call must not abort
## a pipeline run, so the whole body is wrapped in tryCatch.

manifest.record <- function(step, status = "info", elapsed = NA, detail = "",
                            work.dir = Sys.getenv("WORK_DIR", unset = getwd())) {
  tryCatch({
    jobs.dir <- file.path(work.dir, "jobs")
    dir.create(jobs.dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(jobs.dir, "run_manifest.tsv")
    if (!file.exists(path))
      cat("# timestamp\tstep\tstatus\thost\tcommit\telapsed_s\tdetail\n",
          file = path)
    ts   <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    host <- tryCatch(as.character(Sys.info()[["nodename"]]),
                     error = function(e) "unknown")
    commit <- tryCatch(
      suppressWarnings(system2("git", c("-C", work.dir, "rev-parse",
                                        "--short", "HEAD"),
                               stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0))
    commit <- if (length(commit) == 0 || !nzchar(commit[1])) "unknown"
              else commit[1]
    clean <- function(s) gsub("[\t\n]", " ", as.character(s))
    el <- if (length(elapsed) != 1 || is.na(elapsed) ||
              !nzchar(as.character(elapsed))) "-" else as.character(elapsed)
    cat(paste(ts, clean(step), clean(status), clean(host), commit,
              el, clean(detail), sep = "\t"), "\n",
        file = path, append = TRUE, sep = "")
  }, error = function(e) invisible(NULL))
  invisible(NULL)
}
