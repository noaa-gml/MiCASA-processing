## lib/provenance.r -- pipeline provenance for netCDF global attributes.
##
## Builds the CF/ACDD-style global-attribute set that every netCDF the
## MiCASA pipeline writes should carry: which software produced it (git
## commit + version), when, on what host, from which inputs (path +
## SHA-256), plus citation metadata. R writers (diurnalize-ERA5.r ...)
## call provenance.attrs() and nc.write.attrs(), or the nc.write.provenance()
## convenience wrapper.
##
## Citation constants (institution, DOI, ...) come from lib/provenance.conf
## -- a KEY="VALUE" file shared with lib/provenance.py and the shell
## stampers, so the DOI lives in exactly one place.
##
## Building the attribute list needs only base R (no ncdf4), so
## provenance.attrs() is unit-tested standalone -- tests/test_provenance.r.

## ---- citation config ------------------------------------------------------

## Parse lib/provenance.conf (KEY="VALUE" lines, '#' comments) into a named
## list. Missing file -> empty list (callers fall back to defaults).
prov.load.conf <- function(work.dir) {
  path <- file.path(work.dir, "lib", "provenance.conf")
  out <- list()
  if (!file.exists(path)) return(out)
  for (ln in readLines(path, warn = FALSE)) {
    ln <- trimws(ln)
    if (nchar(ln) == 0 || startsWith(ln, "#")) next
    eq <- regexpr("=", ln, fixed = TRUE)
    if (eq < 1) next
    key <- trimws(substr(ln, 1, eq - 1))
    val <- trimws(substr(ln, eq + 1, nchar(ln)))
    if (nchar(val) >= 2 &&
        substr(val, 1, 1) == substr(val, nchar(val), nchar(val)) &&
        substr(val, 1, 1) %in% c("\"", "'")) {
      val <- substr(val, 2, nchar(val) - 1)
    }
    out[[key]] <- val
  }
  out
}

## ---- git / host / checksum helpers ---------------------------------------

## Run `git` in repo.dir; return the trimmed first stdout line, or `default`
## if git is unavailable / the command fails. system2() execs directly (no
## shell), so `args` is a plain character vector -- no quoting needed.
prov.git <- function(repo.dir, args, default = "unknown") {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo.dir, args),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(out) == 0 || !nzchar(out[1])) default else trimws(out[1])
}

prov.git.commit  <- function(repo.dir) prov.git(repo.dir, c("rev-parse", "HEAD"))
prov.git.version <- function(repo.dir)
  prov.git(repo.dir, c("describe", "--tags", "--always", "--dirty"))

## SHA-256 of a file via the system `sha256sum` (GNU coreutils, present on
## Orion and ubuntu CI). Returns NA_character_ if the file is absent or the
## tool is unavailable.
prov.file.sha256 <- function(path) {
  if (length(path) != 1 || is.na(path) || !nzchar(path) || !file.exists(path))
    return(NA_character_)
  out <- tryCatch(
    suppressWarnings(system2("sha256sum", path, stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(out) == 0) return(NA_character_)
  sub("\\s.*$", "", out[1])
}

## ISO-8601 UTC timestamp, e.g. "2026-05-16T05:21:09Z".
prov.timestamp <- function()
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

prov.host <- function()
  tryCatch(as.character(Sys.info()[["nodename"]]), error = function(e) "unknown")

## ---- attribute builder ----------------------------------------------------

## Build the ordered named list of global attributes for a netCDF written
## by `step`. Every value is a character scalar (prec = "text" in ncdf4).
##
##   step      producing script, e.g. "diurnalize-ERA5.r"
##   work.dir  pipeline checkout dir (git repo root; holds lib/)
##   title     ACDD title       (optional; omitted when NA)
##   summary   ACDD summary     (optional; omitted when NA)
##   inputs    named list/character of input file paths; each emits
##             input_<name> and input_<name>_sha256
##   extra     named list of extra attributes, merged last (wins on key
##             collision) -- e.g. flux_fit_method, micasa_version
provenance.attrs <- function(step, work.dir,
                             title = NA, summary = NA,
                             inputs = list(), extra = list()) {
  conf <- prov.load.conf(work.dir)
  getc <- function(k, d) if (!is.null(conf[[k]])) conf[[k]] else d
  doi     <- getc("MICASA_DOI", "PENDING")
  landing <- getc("MICASA_LANDING_PAGE", "PENDING")
  pipe    <- getc("MICASA_PROV_PIPELINE", "MiCASA-processing")
  pipeurl <- getc("MICASA_PROV_PIPELINE_URL", "")
  commit  <- prov.git.commit(work.dir)
  version <- prov.git.version(work.dir)
  has.doi     <- nzchar(doi)     && !identical(doi, "PENDING")
  has.landing <- nzchar(landing) && !identical(landing, "PENDING")

  ## references: pipeline URL always; DOI / landing only once registered.
  refs <- sprintf("%s pipeline: %s", pipe, pipeurl)
  if (has.doi)     refs <- c(refs, sprintf("dataset DOI: https://doi.org/%s", doi))
  if (has.landing) refs <- c(refs, sprintf("dataset landing page: %s", landing))

  ts <- prov.timestamp()
  a <- list()
  a[["Conventions"]] <- getc("MICASA_PROV_CONVENTIONS", "CF-1.10, ACDD-1.3")
  if (!is.na(title))   a[["title"]]   <- title
  if (!is.na(summary)) a[["summary"]] <- summary
  a[["institution"]]                 <- getc("MICASA_PROV_INSTITUTION", "")
  a[["source"]]                      <- sprintf("%s pipeline, step %s", pipe, step)
  a[["references"]]                  <- paste(refs, collapse = " ; ")
  a[["license"]]                     <- getc("MICASA_PROV_LICENSE", "")
  a[["creator_name"]]                <- getc("MICASA_PROV_CREATOR_NAME", "")
  a[["creator_url"]]                 <- getc("MICASA_PROV_CREATOR_URL", "")
  a[["date_created"]]                <- ts
  a[["processing_pipeline"]]         <- pipe
  a[["processing_pipeline_url"]]     <- pipeurl
  a[["processing_pipeline_commit"]]  <- commit
  a[["processing_pipeline_version"]] <- version
  a[["processing_step"]]             <- step
  a[["processing_host"]]             <- prov.host()
  if (has.doi) a[["doi"]] <- doi

  ## input files: input_<name> and input_<name>_sha256
  nms <- names(inputs)
  for (i in seq_along(inputs)) {
    nm <- nms[i]
    p  <- inputs[[i]]
    if (is.null(nm) || !nzchar(nm) || is.null(p) || length(p) != 1 || is.na(p))
      next
    a[[sprintf("input_%s", nm)]] <- as.character(p)
    s <- prov.file.sha256(as.character(p))
    a[[sprintf("input_%s_sha256", nm)]] <- if (is.na(s)) "unavailable" else s
  }

  ## history: one CF audit line. On a freshly nc_create()d file the writer
  ## SETs this (no prior history == append); downstream NCO tools then
  ## append their own lines.
  a[["history"]] <- sprintf("%s: created by %s [%s %s, commit %s]",
                            ts, step, pipe, version, substr(commit, 1, 12))

  ## extra wins on key collision
  for (k in names(extra)) a[[k]] <- as.character(extra[[k]])
  a
}

## ---- netCDF writers -------------------------------------------------------

## Write a named attribute list onto an open ncdf4 file as global attributes
## (varid 0). NULL / NA / empty values are skipped.
nc.write.attrs <- function(ncf, attrs) {
  for (k in names(attrs)) {
    v <- attrs[[k]]
    if (is.null(v) || length(v) != 1 || is.na(v) || !nzchar(as.character(v)))
      next
    ncatt_put(ncf, 0, k, attval = as.character(v), prec = "text")
  }
  invisible(ncf)
}

## Convenience: build provenance.attrs(...) and write them in one call.
nc.write.provenance <- function(ncf, step, work.dir,
                                title = NA, summary = NA,
                                inputs = list(), extra = list()) {
  nc.write.attrs(ncf, provenance.attrs(step, work.dir, title = title,
                                       summary = summary, inputs = inputs,
                                       extra = extra))
}
