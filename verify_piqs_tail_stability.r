#!/usr/bin/env Rscript

## verify_piqs_tail_stability.r
##
## Phase 3 verification helper: compare the right-edge (last 3 months
## pre-pad) PIQS coefficients in the current fit.piqs.rda to those in a
## snapshot fit.piqs.rda.prev. Used to validate that PAD_RIGHT=2 actually
## stabilises the tail under NRT re-fit.
##
## On first run (no .prev file): snapshot the current .rda as .prev and
## emit status=snapshot_established.
##
## On subsequent runs: load both, restrict to the LAST 3 segments common
## between them, compute |delta|/|coef| per (i,j,seg) for each of a, b, c
## (gpp and resp), report max-rel and median-rel.
##
## Output JSON written to argv[1] (default verify_piqs_tail_stability.json).
##
## "Common" tail logic: if the previous fit ended at month T-3 and the
## current ends at T (one more month of NRT data ingested), the
## "last 3 months" for comparison are (T-3, T-2, T-1) -- the months that
## existed in BOTH fits but moved to non-edge position in the new fit.
## We compare those interior coefs in the new fit to their old edge-coef
## values in the prior fit.

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
out.json <- if (length(args) >= 1) args[1] else "verify_piqs_tail_stability.json"

cur.path  <- "fit.piqs.rda"
prev.path <- "fit.piqs.rda.prev"

if (!file.exists(cur.path)) {
  write_json(list(status = "no_current_fit"), out.json, pretty = TRUE, auto_unbox = TRUE)
  cat("no fit.piqs.rda; nothing to do\n")
  quit(save = "no")
}

if (!file.exists(prev.path)) {
  ## First run: just snapshot.
  file.copy(cur.path, prev.path, overwrite = FALSE)
  write_json(list(status        = "snapshot_established",
                  snapshot_path = prev.path,
                  taken_at      = format(Sys.time(), tz = "UTC", usetz = TRUE)),
             out.json, pretty = TRUE, auto_unbox = TRUE)
  cat(sprintf("baseline snapshot taken: %s\n", prev.path))
  quit(save = "no")
}

## Both files exist: load them as separate environments to avoid name clash.
e.cur  <- new.env()
e.prev <- new.env()
load(cur.path,  envir = e.cur)
load(prev.path, envir = e.prev)

cur.time  <- e.cur$piqsfit.time
prev.time <- e.prev$piqsfit.time

## Find the last 3 segments common to both (intersection of left-edge times,
## taking the rightmost 3).
common <- intersect(cur.time, prev.time)
if (length(common) < 3) {
  write_json(list(status = "insufficient_overlap",
                  n_common = length(common)),
             out.json, pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no")
}
common.tail <- tail(sort(common), 3)
cur.idx  <- match(common.tail, cur.time)
prev.idx <- match(common.tail, prev.time)

extract <- function(env, idxs) {
  list(gpp_a  = env$piqsfit.gpp$a [, , idxs],
       gpp_b  = env$piqsfit.gpp$b [, , idxs],
       gpp_c  = env$piqsfit.gpp$c [, , idxs],
       rtot_a = env$piqsfit.resp$a[, , idxs],
       rtot_b = env$piqsfit.resp$b[, , idxs],
       rtot_c = env$piqsfit.resp$c[, , idxs])
}
cur.coefs  <- extract(e.cur,  cur.idx)
prev.coefs <- extract(e.prev, prev.idx)

## Active mask: cells where both fits have non-NA coefs.
active <- !is.na(cur.coefs$gpp_a) & !is.na(prev.coefs$gpp_a)

rel.diff <- function(cur.arr, prev.arr) {
  d <- abs(cur.arr - prev.arr)
  m <- pmax(abs(prev.arr), .Machine$double.eps)
  rel <- d / m
  rel[!active] <- NA
  list(max_rel    = max(rel, na.rm = TRUE),
       median_rel = median(rel, na.rm = TRUE))
}

stats.gpp_a  <- rel.diff(cur.coefs$gpp_a,  prev.coefs$gpp_a)
stats.gpp_b  <- rel.diff(cur.coefs$gpp_b,  prev.coefs$gpp_b)
stats.gpp_c  <- rel.diff(cur.coefs$gpp_c,  prev.coefs$gpp_c)
stats.rtot_a <- rel.diff(cur.coefs$rtot_a, prev.coefs$rtot_a)
stats.rtot_b <- rel.diff(cur.coefs$rtot_b, prev.coefs$rtot_b)
stats.rtot_c <- rel.diff(cur.coefs$rtot_c, prev.coefs$rtot_c)

snapshot.age.hours <- as.numeric(
  difftime(Sys.time(), file.info(prev.path)$mtime, units = "hours"))

out <- list(
  status                = "compared",
  common_tail_times     = common.tail,
  n_segments_compared   = sum(active) %/% 3,
  snapshot_path         = prev.path,
  snapshot_age_hours    = snapshot.age.hours,
  max_rel_diff_gpp      = max(stats.gpp_a$max_rel, stats.gpp_b$max_rel, stats.gpp_c$max_rel),
  max_rel_diff_rtot     = max(stats.rtot_a$max_rel, stats.rtot_b$max_rel, stats.rtot_c$max_rel),
  median_rel_diff_gpp   = max(stats.gpp_a$median_rel, stats.gpp_b$median_rel, stats.gpp_c$median_rel),
  median_rel_diff_rtot  = max(stats.rtot_a$median_rel, stats.rtot_b$median_rel, stats.rtot_c$median_rel),
  per_coef = list(
    gpp_a  = stats.gpp_a,  gpp_b  = stats.gpp_b,  gpp_c  = stats.gpp_c,
    rtot_a = stats.rtot_a, rtot_b = stats.rtot_b, rtot_c = stats.rtot_c
  )
)

write_json(out, out.json, pretty = TRUE, auto_unbox = TRUE, na = "string")
cat(sprintf("compared %d cells across last 3 common segments; max_rel GPP=%.3e, Rtot=%.3e\n",
            out$n_segments_compared, out$max_rel_diff_gpp, out$max_rel_diff_rtot))
