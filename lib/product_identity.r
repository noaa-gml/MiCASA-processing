## lib/product_identity.r -- refuse to mix incompatible physics in one tree.
##
## The pipeline's product-defining settings (fitter, respiration driver and
## response function, polar clip, ATMC removal, climatology mode) are supplied
## through the environment. Nothing used to stop a tree being assembled from
## runs with different settings: re-diurnalizing one month with a different
## MICASA_RESP_TEMPFUN, or with MICASA_ATMC on, wrote into the same paths and
## left no trace at the tree level. The result is a product whose months do not
## share physics -- an ~4 PgC/yr discontinuity in the ATMC case -- and which
## looks exactly like a homogeneous one from the outside.
##
## `run_climatology_prior.sh` already solved this for its own output with a
## `.micasa-climatology-prior` marker it refuses to write without. This
## generalizes that idea to every product-defining knob, and puts it in the
## WRITER rather than in a driver, so it holds however diurnalize is invoked
## (drivers, sbatch, or a bare Rscript).
##
##   first write into a tree  -> record the identity
##   later writes             -> compare; refuse on any disagreement
##   deliberate change        -> new output directory, or MICASA_ALLOW_MIXED=1
##
## The marker is plain `key: value` text so it can be read, diffed and grepped
## without tooling. It is advisory for humans and binding for the pipeline.

PRODUCT.MARKER <- ".micasa-product"

## Format the identity as sorted "key: value" lines (stable for diffing).
product.identity.format <- function(id) {
  keys <- sort(names(id))
  paste0(sprintf("%-22s %s", paste0(keys, ":"), unlist(id[keys])), collapse = "\n")
}

## Compare two identities; returns a character vector of disagreeing keys.
##
## Only keys the CURRENT run declares are compared. A marker written by an older
## version may carry keys this one no longer tracks -- that is a schema change,
## not a physics disagreement, and refusing on it would strand every existing
## tree at each upgrade. The corollary is deliberate: a key dropped from the
## identity stops being enforced, so drop one only when it stops defining the
## product.
product.identity.diff <- function(have, want) {
  keys <- names(want)
  keys[vapply(keys, function(k) {
    a <- if (is.null(have[[k]])) "<absent>" else as.character(have[[k]])
    b <- if (is.null(want[[k]])) "<absent>" else as.character(want[[k]])
    !identical(a, b)
  }, logical(1))]
}

product.identity.read <- function(dir) {
  path <- file.path(dir, PRODUCT.MARKER)
  if (!file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln)) & !startsWith(trimws(ln), "#")]
  kv <- regmatches(ln, regexpr("^[^:]+:", ln))
  if (!length(kv)) return(list())
  out <- as.list(trimws(sub("^[^:]+:", "", ln)))
  names(out) <- trimws(sub(":$", "", kv))
  out
}

## Check the tree's identity against this run's, and record it on first write.
##
##   dir          output directory (created if absent)
##   id           named list of this run's product-defining settings
##   allow.mixed  TRUE to warn instead of stopping (MICASA_ALLOW_MIXED=1)
##
## Returns invisibly: "created" | "match" | "mixed".
product.identity.enforce <- function(dir, id, allow.mixed = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, PRODUCT.MARKER)
  have <- product.identity.read(dir)

  if (is.null(have)) {
    ## An existing tree with output but no marker predates this check; adopt it
    ## rather than refuse, but say so -- we cannot vouch for what is already there.
    existing <- length(list.files(dir, pattern = "^fluxes_[0-9]{6}\\.nc$"))
    writeLines(c("# MiCASA product identity -- what this tree's files were made with.",
                 "# Written on first output; later runs must agree or be refused.",
                 "# See lib/product_identity.r.",
                 product.identity.format(id)), path)
    if (existing > 0)
      warning(sprintf("%s: adopting %d pre-existing file(s) into a new %s; their settings were not recorded and are assumed to match.",
                      dir, existing, PRODUCT.MARKER), immediate. = TRUE)
    return(invisible("created"))
  }

  bad <- product.identity.diff(have, id)
  if (!length(bad)) return(invisible("match"))

  detail <- paste(vapply(bad, function(k)
    sprintf("    %-20s tree=%-14s this run=%s", k,
            if (is.null(have[[k]])) "<absent>" else have[[k]],
            if (is.null(id[[k]]))   "<absent>" else id[[k]]),
    character(1)), collapse = "\n")
  msg <- sprintf("product identity mismatch in %s\n%s\n  The files already there were made with different physics. Writing here would produce a tree whose months do not share a configuration.\n  Use a different output directory, or set MICASA_ALLOW_MIXED=1 to override deliberately.",
                 dir, detail)
  if (allow.mixed) {
    warning(msg, immediate. = TRUE)
    return(invisible("mixed"))
  }
  stop(msg)
}
