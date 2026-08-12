#!/usr/bin/env Rscript
## test_product_identity.r -- the guard that stops a tree being assembled from
## runs with different physics (lib/product_identity.r).
##
## Usage:  Rscript tests/test_product_identity.r
##
## Before this guard, re-diurnalizing one month with a different
## MICASA_RESP_TEMPFUN -- or with MICASA_ATMC on -- wrote into the same paths
## and left no tree-level trace, producing a product whose months do not share
## physics and which is indistinguishable from a homogeneous one.
##
##   CREATE   first write records the identity
##   MATCH    an identical configuration is allowed through
##   REFUSE   any disagreeing key stops the run, and the message names the key
##   OVERRIDE MICASA_ALLOW_MIXED downgrades the refusal to a warning
##   ADOPT    a pre-existing tree with no marker is adopted, with a warning
##   ROUNDTRIP the marker reads back exactly what was written
suppressWarnings(source(file.path(Sys.getenv("WORK_DIR", getwd()), "lib",
                                  "product_identity.r")))

fails <- 0L
ok <- function(label, cond, detail = "") {
  cat(sprintf("%-58s %s%s\n", label, if (isTRUE(cond)) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!isTRUE(cond)) fails <<- fails + 1L
}
tmp <- function() { d <- tempfile("prodid"); dir.create(d); d }

ID <- list(micasa_version = "v1", flux_fit_method = "pchip",
           resp_driver = "airtemp", resp_tempfun = "q10",
           polar_clip = "conserve", atmc_removal = "off")

## ---- CREATE ---------------------------------------------------------------
d <- tmp()
r <- product.identity.enforce(d, ID)
ok("CREATE   first write records the identity",
   identical(r, "created") && file.exists(file.path(d, PRODUCT.MARKER)))

## ---- MATCH ----------------------------------------------------------------
ok("MATCH    an identical configuration is allowed",
   identical(product.identity.enforce(d, ID), "match"))

## ---- ROUNDTRIP ------------------------------------------------------------
back <- product.identity.read(d)
ok("ROUNDTRIP marker reads back what was written",
   length(product.identity.diff(back, ID)) == 0,
   paste(product.identity.diff(back, ID), collapse = ","))

## ---- REFUSE ---------------------------------------------------------------
## Every product-defining key must be load-bearing, not just the new one.
for (k in names(ID)) {
  changed <- ID
  changed[[k]] <- paste0(ID[[k]], "-X")
  res <- try(product.identity.enforce(d, changed), silent = TRUE)
  refused <- inherits(res, "try-error")
  named   <- refused && grepl(k, conditionMessage(attr(res, "condition")), fixed = TRUE)
  ok(sprintf("REFUSE   changing %-16s stops the run", k), refused && named)
}

## ---- SCHEMA: a marker with an extra legacy key must not refuse -------------
d3 <- tmp()
product.identity.enforce(d3, c(ID, list(retired_key = "old")))
ok("SCHEMA   marker key we no longer track is ignored",
   identical(product.identity.enforce(d3, ID), "match"))

## ---- OVERRIDE -------------------------------------------------------------
changed <- ID; changed$atmc_removal <- "q10"
res <- withCallingHandlers(
  product.identity.enforce(d, changed, allow.mixed = TRUE),
  warning = function(w) invokeRestart("muffleWarning"))
ok("OVERRIDE allow.mixed downgrades refusal to a warning",
   identical(res, "mixed"))

## ---- ADOPT ----------------------------------------------------------------
d2 <- tmp()
invisible(file.create(file.path(d2, "fluxes_202107.nc")))
warned <- FALSE
res <- withCallingHandlers(
  product.identity.enforce(d2, ID),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
ok("ADOPT    pre-existing unmarked tree adopted, with a warning",
   identical(res, "created") && warned)

## ---- the guard must not fire on a fresh directory --------------------------
ok("CREATE   a fresh tree is never refused",
   identical(product.identity.enforce(tmp(), ID), "created"))

cat(sprintf("\n%s: %d check(s) failed\n", if (fails == 0L) "OK" else "FAILURES", fails))
quit(status = if (fails == 0L) 0L else 1L)
