#!/usr/bin/env Rscript

## verify_piqs_tail_horizon.r
##
## Like verify_piqs_tail_stability.r but sweeps over k = 1..12 instead
## of fixing k = 3. For each k, computes median |delta/coef| over the
## last k common segments of fit.piqs.rda vs fit.piqs.rda.prev.
##
## Reports the largest k for which median |delta/coef| > 1% — the
## "propagation horizon" of input changes through the global PIQS solve.

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
out.json <- if (length(args) >= 1) args[1] else "verify_piqs_tail_horizon.json"

cur.path  <- "fit.piqs.rda"
prev.path <- "fit.piqs.rda.prev"

if (!file.exists(cur.path) || !file.exists(prev.path)) {
  write_json(list(status = "no_pair_to_compare"), out.json,
             pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no")
}

e.cur  <- new.env()
e.prev <- new.env()
load(cur.path,  envir = e.cur)
load(prev.path, envir = e.prev)

common <- sort(intersect(e.cur$piqsfit.time, e.prev$piqsfit.time))
n.common <- length(common)
if (n.common < 4) {
  write_json(list(status = "insufficient_overlap", n_common = n.common),
             out.json, pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no")
}

per.k <- list()
horizon <- 0L
for (k in 1:min(12, n.common - 1)) {
  tail.times <- tail(common, k)
  cur.idx  <- match(tail.times, e.cur$piqsfit.time)
  prev.idx <- match(tail.times, e.prev$piqsfit.time)

  abs.deltas <- c()
  rel.deltas <- c()
  for (which_set in list(list("gpp",  e.cur$piqsfit.gpp,  e.prev$piqsfit.gpp),
                         list("rtot", e.cur$piqsfit.resp, e.prev$piqsfit.resp))) {
    cur.set  <- which_set[[2]]
    prev.set <- which_set[[3]]
    for (coef in c("a", "b", "c")) {
      cur.arr  <- cur.set [[coef]][, , cur.idx]
      prev.arr <- prev.set[[coef]][, , prev.idx]
      active   <- !is.na(cur.arr) & !is.na(prev.arr)
      d <- abs(cur.arr - prev.arr)[active]
      m <- pmax(abs(prev.arr[active]), .Machine$double.eps)
      abs.deltas <- c(abs.deltas, d)
      rel.deltas <- c(rel.deltas, d / m)
    }
  }
  med.rel <- median(rel.deltas, na.rm = TRUE)
  per.k[[length(per.k) + 1]] <- list(k = k, median_rel = med.rel)
  if (med.rel > 0.01) horizon <- k
}

write_json(list(
  status   = "compared",
  n_common = n.common,
  horizon_1pct_median = horizon,
  per_k    = per.k
), out.json, pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("propagation horizon (median |delta/coef| > 1%%): last %d months\n",
            horizon))
