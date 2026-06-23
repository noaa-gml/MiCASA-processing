#!/usr/bin/env Rscript
##
## verify_pchip_invariants.r — PCHIP fit invariants
##
## Loads fit.piqs.rda and validates two invariants the PCHIP fitter is
## supposed to guarantee by Fritsch-Carlson construction:
##
##   (a) per-segment sign: f(τ) = a + bτ + cτ² is ≤ 0 (GPP) or ≥ 0 (Rh)
##       on τ ∈ [0, L]. Endpoints + interior vertex are checked
##       analytically.
##   (b) C¹ continuity at interior knots: f from the right of segment k-1
##       must equal f from the left of segment k (i.e., a_k).
##
## Outputs JSON; consumed by verify_v2.py Checks 18.1 / 18.2.

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) > 0) args[1] else "verify_pchip_invariants.json"

load("fit.piqs.rda")

t_sec <- as.numeric(piqsfit.time)
nmon  <- length(t_sec)
# Segment lengths in seconds. Segment i is [t_i, t_{i+1}). The last
## stored segment (index nmon) extends to the next month boundary which
## isn't in piqsfit.time -- approximate it by the median segment length
## so we have a coverage value to use.
dt <- diff(t_sec)
dt <- c(dt, median(dt))
stopifnot(length(dt) == nmon)

eps_zero <- 1e-30   # treat coef == 0 as machine zero
tol_fp   <- 1e-12   # floating-point tolerance below which we don't flag

## ---- (a) per-segment sign check ----------------------------------------
##
## sign_required: -1 means f ≤ 0 expected (GPP); +1 means f ≥ 0 (Rh).
## Returns total land segments scanned and total violations.
check_sign <- function(coef, sign_required) {
  max_mag <- 0.0
  ## storage convention from write_pchip.r:
  ##   f(tau) = A * tau^2 + B * tau + C, tau in [0, L]
  ## (i.e., coef$a is the tau^2 coeff, coef$c is the constant.)
  a <- coef$a; b <- coef$b; c <- coef$c
  total <- 0L; viol <- 0L
  for (im in seq_len(nmon)) {
    A <- a[, , im]; B <- b[, , im]; C <- c[, , im]
    L <- dt[im]
    f0 <- C
    fL <- A * L^2 + B * L + C
    ## interior vertex when |A| not negligible; tau_v in (0, L)
    safe_a <- ifelse(abs(A) < eps_zero, NA_real_, A)
    tau_v  <- -B / (2 * safe_a)
    f_v    <- C - B^2 / (4 * safe_a)
    inside <- !is.na(tau_v) & tau_v > 0 & tau_v < L
    ## fmax / fmin over [0, L]
    fmax <- pmax(f0, fL, na.rm = TRUE)
    fmin <- pmin(f0, fL, na.rm = TRUE)
    fmax[inside] <- pmax(fmax[inside], f_v[inside])
    fmin[inside] <- pmin(fmin[inside], f_v[inside])
    ## land mask: any non-zero coefficient
    land <- (abs(A) + abs(B) + abs(C)) > eps_zero
    if (sign_required > 0) {
      v <- (fmin < -tol_fp) & land
      mag <- ifelse(v, -fmin, 0)  # how far below zero
    } else {
      v <- (fmax >  tol_fp) & land
      mag <- ifelse(v, fmax, 0)   # how far above zero
    }
    total   <- total + sum(land, na.rm = TRUE)
    viol    <- viol  + sum(v,    na.rm = TRUE)
    max_mag <- max(max_mag, max(mag, na.rm = TRUE))
  }
  list(total = total, viol = viol, max_mag = max_mag)
}

## ---- (b) C¹ continuity at interior knots -------------------------------
c1_jump <- function(coef) {
  ## storage convention: f(tau) = A tau^2 + B tau + C
  ## right-limit of segment k-1 at tau = L_{k-1}: A_{k-1} L^2 + B_{k-1} L + C_{k-1}
  ## left-limit of segment k at tau = 0:           C_k
  a <- coef$a; b <- coef$b; c <- coef$c
  max_jump <- 0.0
  for (im in 2:nmon) {
    L <- dt[im - 1]
    f_right_of_prev <- a[, , im - 1] * L^2 + b[, , im - 1] * L + c[, , im - 1]
    f_left_of_curr  <- c[, , im]
    jump <- abs(f_right_of_prev - f_left_of_curr)
    max_jump <- max(max_jump, max(jump, na.rm = TRUE))
  }
  max_jump
}

cat("Checking GPP per-segment sign (≤ 0)...\n")
gpp_sign <- check_sign(piqsfit.gpp,  sign_required = -1)
cat(sprintf("  GPP segments: %d, violations: %d, max |violation|: %.3e\n",
            gpp_sign$total, gpp_sign$viol, gpp_sign$max_mag))

cat("Checking Rh per-segment sign (≥ 0)...\n")
rh_sign  <- check_sign(piqsfit.resp, sign_required = +1)
cat(sprintf("  Rh segments:  %d, violations: %d, max |violation|: %.3e\n",
            rh_sign$total, rh_sign$viol, rh_sign$max_mag))

cat("Checking GPP C¹ continuity at interior knots...\n")
gpp_c1 <- c1_jump(piqsfit.gpp)
cat(sprintf("  GPP max |jump|: %.3e\n", gpp_c1))

cat("Checking Rh C¹ continuity at interior knots...\n")
rh_c1  <- c1_jump(piqsfit.resp)
cat(sprintf("  Rh max |jump|: %.3e\n", rh_c1))

result <- list(
  gpp_seg_total      = gpp_sign$total,
  gpp_seg_violations = gpp_sign$viol,
  gpp_seg_max_mag    = gpp_sign$max_mag,
  rh_seg_total       = rh_sign$total,
  rh_seg_violations  = rh_sign$viol,
  rh_seg_max_mag     = rh_sign$max_mag,
  gpp_c1_max_jump    = gpp_c1,
  rh_c1_max_jump     = rh_c1,
  fitter             = if (exists("piqsfit.meta")) piqsfit.meta$fitter else "unknown",
  written_at         = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), out_path)
} else {
  ## minimal hand-written JSON if jsonlite isn't available
  fmt <- function(x) {
    if (is.character(x)) sprintf('"%s"', x) else format(x, scientific = FALSE)
  }
  pairs <- vapply(names(result),
                  function(k) sprintf('  "%s": %s', k, fmt(result[[k]])),
                  character(1))
  writeLines(c("{", paste(pairs, collapse = ",\n"), "}"), out_path)
}
cat(sprintf("Wrote %s\n", out_path))
