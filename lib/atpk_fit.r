## lib/atpk_fit.r -- 1-D temporal area-to-point (ATP) kriging fitter core.
##
## Disaggregates a monthly-mean flux series into a sub-monthly reconstruction
## that (i) re-aggregates to each month's mean EXACTLY (coherence / pycnophylactic
## property), (ii) carries a kriging VARIANCE -- the principled prior-uncertainty
## the deterministic splines (PCHIP/PPM/PIQS) lack, (iii) is made non-negative for
## a one-signed quantity (GPP) by a selective coherent projection / QP.
##
## Block data = the monthly means; we krige Ns sub-monthly points per month with
## an exponential covariance, then represent each month as a MASS-PRESERVING
## quadratic fit to its points so the on-disk layout matches write_pchip.r
## ((a,b,c) per piece) and diurnalize-ERA5.r consumes it unchanged -- PLUS a
## per-piece variance array (the new output).
##
## Refs: Kyriakidis (2004) Geographical Analysis 36(3):259-289; Yoo & Kyriakidis
## (2006) J. Geographical Systems, DOI 10.1007/s10109-006-0036-7. See
## docs/FITTER_COMPARISON.md (4.4) and PROPOSALS (18).
##
## Pure base R; `quadprog` is used for the positivity QP when available, else a
## base-R alternating-projection (POCS) fallback (tests/test_atpk_fit.r).

## Default exponential-covariance range (months). NOTE: fitting the range from
## the monthly-mean autocorrelation is INAPPROPRIATE here -- that autocorrelation
## is dominated by the seasonal cycle (lag-1 ~0.9 => range ~9 mo), which both
## over-smooths and makes the kriging system singular. The sub-monthly process
## correlation cannot be estimated from monthly data, so we use a fixed short
## range (interpolate smoothly between adjacent months, no more); override per
## biome if desired. ATPK_DEFAULT_RANGE is in months.
ATPK_DEFAULT_RANGE <- 1.5

## Selective coherent non-negativity: nudge z to the nearest point that keeps the
## per-block means (B z == ybar) AND has sign(mean)*z >= 0. quadprog if present,
## else alternating projection onto the (affine) coherence set and the sign cone.
atpk.positive <- function(z, blk, ybar, sgn, qp = TRUE) {
  K <- length(z); n <- length(ybar)
  if (qp && requireNamespace("quadprog", quietly = TRUE)) {
    ## min ||z'-z||^2 s.t. B z'=ybar (meq=n), sgn*z' >= 0
    Amat <- t(rbind(blk, diag(sgn, K)))
    bvec <- c(ybar, rep(0, K))
    sol <- tryCatch(quadprog::solve.QP(diag(K), z, Amat, bvec, meq = n),
                    error = function(e) NULL)
    if (!is.null(sol)) return(sol$solution)
  }
  ## fallback: projections onto coherence (per-block mean restore) then sign cone
  for (it in 1:500) {
    z0 <- z
    z <- ifelse(sgn * z < 0, 0, z)                       # project onto sign cone
    for (i in seq_len(n)) {                              # restore each block mean
      idx <- which(blk[i, ] > 0); z[idx] <- z[idx] + (ybar[i] - mean(z[idx]))
    }
    if (max(abs(z - z0)) < 1e-14 * (abs(mean(ybar)) + 1e-30)) break
  }
  z
}

## Per piece: mass-preserving quadratic LS fit f(s)=a s^2+b s+c to the Ns kriged
## points (s = time offset from the piece start, seconds), constrained so the
## integral over the piece equals ybar_i * h_i exactly.
atpk.points.to.abc <- function(z, ss, Ns, n, h, ybar, sgn = 0) {
  ## Fit in the NORMALIZED fraction u = (t - x_i)/D in [0,1] (well-conditioned),
  ## then convert g(u)=A u^2 + B u + C to the seconds-based f(t)=a(t-x)^2+b(t-x)+c
  ## via a = A/D^2, b = B/D, c = C. Mass constraint in u: A/3 + B/2 + C = ybar.
  ## SELECTIVE positivity: where the mass-preserving quadratic would flip sign on
  ## [0,1], fall back to flat (A=B=0, C=ybar) -- sign-safe + mass-exact, triggers
  ## only in near-zero pieces where the sub-monthly shape is negligible.
  X <- cbind(ss^2, ss, 1); XtX <- crossprod(X); Xi <- solve(XtX)
  g <- c(1/3, 1/2, 1); gXig <- as.numeric(t(g) %*% Xi %*% g)
  a <- b <- c0 <- numeric(n)
  for (i in seq_len(n)) {
    zi <- z[((i - 1) * Ns + 1):(i * Ns)]; D <- h[i]
    bols <- Xi %*% crossprod(X, zi)
    beta <- bols + (Xi %*% g) * ((ybar[i] - sum(g * bols)) / gXig)
    A <- beta[1]; B <- beta[2]; C <- beta[3]
    if (sgn != 0) {                                   # check g(u)>=0 (signed) on [0,1]
      vals <- c(C, A + B + C)                         # endpoints
      uv <- if (A != 0) -B / (2 * A) else -1          # vertex
      if (uv > 0 && uv < 1) vals <- c(vals, C - B^2 / (4 * A))
      if (any(sgn * vals < -1e-12 * abs(ybar[i]))) { A <- 0; B <- 0; C <- ybar[i] }
    }
    a[i] <- A / D^2; b[i] <- B / D; c0[i] <- C
  }
  list(a = a, b = b, c = c0)
}

## ---- Windowed production path -------------------------------------------------
## Full-series kriging is O(n^3) per cell (a ~1800x1800 solve for 303 months) and
## intractable over 64800 cells. Because kriging WEIGHTS depend only on the
## covariance geometry (window size m, target-block position p), not the data, we
## precompute them once per (m,p) and apply them as a fixed linear filter to every
## cell. A +-W month window is effectively exact (exp(-W/range) is negligible for
## W >> range). This also makes the method NRT-local (footprint <= W).

## Precompute, for an n-month record, the kriging weights mapping the surrounding
## window's monthly means to each target month's Ns sub-points (+ unit variance).
atpk.window.weights <- function(n, W = 6L, Ns = 6L, range = ATPK_DEFAULT_RANGE,
                                nugget = 1e-3) {
  ss <- (seq_len(Ns) - 0.5) / Ns
  cache <- list()                                   # key "m_p" -> list(lam[Ns,m], uvar[Ns])
  geom <- function(m, p) {
    key <- paste(m, p, sep = "_"); if (!is.null(cache[[key]])) return(cache[[key]])
    tk <- rep(seq_len(m), each = Ns) - 1 + rep(ss, m); K <- length(tk)
    blk <- matrix(0, m, K); for (i in 1:m) blk[i, ((i-1)*Ns+1):(i*Ns)] <- 1/Ns
    C <- exp(-abs(outer(tk, tk, "-")) / range); diag(C) <- diag(C) + nugget
    CBB <- blk %*% C %*% t(blk); CBp <- blk %*% C
    A <- rbind(cbind(CBB, 1), c(rep(1, m), 0))
    cols <- ((p-1)*Ns+1):(p*Ns)
    sol <- tryCatch(solve(A, rbind(CBp[, cols, drop=FALSE], 1)),
                    error = function(e) solve(A + diag(c(rep(10*nugget,m),0)), rbind(CBp[,cols,drop=FALSE],1)))
    lam <- sol[1:m, , drop=FALSE]; mu <- sol[m+1, ]
    uvar <- pmax((1+nugget) - colSums(CBp[, cols, drop=FALSE] * lam) - mu, 0)
    cache[[key]] <<- list(lam = lam, uvar = uvar); cache[[key]]
  }
  lo <- pmax(1L, seq_len(n) - W); hi <- pmin(n, seq_len(n) + W)
  out <- vector("list", n)
  for (i in seq_len(n)) out[[i]] <- c(list(lo = lo[i], hi = hi[i]),
                                      geom(hi[i] - lo[i] + 1L, i - lo[i] + 1L))
  list(per = out, Ns = Ns, ss = ss)
}

## Apply precomputed window weights to one cell's monthly means -> (a,b,c,var).
atpk.apply.series <- function(ybar, h, ww, positive = TRUE) {
  n <- length(ybar); Ns <- ww$Ns; ss <- ww$ss
  sill <- var(ybar)
  if (!is.finite(sill) || sill <= 0 || all(ybar == 0))
    return(list(a = rep(0,n), b = rep(0,n), c = ybar, var = rep(0,n), dormant = TRUE))
  z <- numeric(n * Ns); vpc <- numeric(n)
  for (i in seq_len(n)) {
    w <- ww$per[[i]]; wm <- ybar[w$lo:w$hi]
    zi <- as.vector(t(w$lam) %*% wm)
    zi <- zi + (ybar[i] - mean(zi))                 # exact coherence (belt+braces)
    z[((i-1)*Ns+1):(i*Ns)] <- zi
    vpc[i] <- sill * mean(w$uvar)
  }
  sgn.out <- if (positive) { s0 <- sign(mean(ybar)); if (s0 == 0) 1 else s0 } else 0
  abc <- atpk.points.to.abc(z, ss, Ns, n, h, ybar, sgn = sgn.out)
  list(a = abc$a, b = abc$b, c = abc$c, var = vpc, dormant = FALSE)
}

## Fit one cell. x: knot times (length n+1, seconds). ybar: monthly means (n).
## Returns list(a,b,c,var,dormant): per-piece quadratic coeffs (mass-preserving),
## per-piece kriging variance, and a dormant flag for ~zero-variance cells.
atpk.fit.cell <- function(x, ybar, Ns = 6L, range = NULL, nugget = 1e-3,
                          positive = TRUE, qp = TRUE) {
  n <- length(ybar); h <- diff(x)
  sill <- var(ybar)
  if (!is.finite(sill) || sill <= 0 || all(ybar == 0)) {
    return(list(a = rep(0, n), b = rep(0, n), c = ybar, var = rep(0, n), dormant = TRUE))
  }
  if (is.null(range)) range <- ATPK_DEFAULT_RANGE
  ss <- (seq_len(Ns) - 0.5) / Ns
  tk <- rep(seq_len(n), each = Ns) - 1 + rep(ss, n)     # time in months
  K <- length(tk)
  blk <- matrix(0, n, K); for (i in 1:n) blk[i, ((i - 1) * Ns + 1):(i * Ns)] <- 1 / Ns
  ## Solve on the UNIT covariance (correlation structure only); kriging weights
  ## are sill-independent, and using sill here would wreck the system scaling
  ## (CBB ~ sill ~ 1e-13 vs the Lagrange row's 1). Variance is scaled by sill at
  ## the end.
  C <- exp(-abs(outer(tk, tk, "-")) / range); diag(C) <- diag(C) + nugget
  CBB <- blk %*% C %*% t(blk); CBp <- blk %*% C
  A <- rbind(cbind(CBB, 1), c(rep(1, n), 0))
  sol <- tryCatch(solve(A, rbind(CBp, 1)), error = function(e)
          tryCatch(solve(A + diag(c(rep(10 * nugget, n), 0)), rbind(CBp, 1)),
                   error = function(e2) NULL))
  if (is.null(sol)) return(list(a = rep(0, n), b = rep(0, n), c = ybar,
                                var = rep(0, n), dormant = TRUE))
  lam <- sol[1:n, , drop = FALSE]; mu <- sol[n + 1, ]
  z <- as.vector(t(lam) %*% ybar)
  v <- sill * pmax((1 + nugget) - colSums(CBp * lam) - mu, 0)
  for (i in 1:n) { idx <- ((i - 1) * Ns + 1):(i * Ns); z[idx] <- z[idx] + (ybar[i] - mean(z[idx])) }
  if (positive) {
    sgn <- sign(mean(ybar)); if (sgn == 0) sgn <- 1
    if (any(sgn * z < -1e-9 * sqrt(sill))) z <- atpk.positive(z, blk, ybar, sgn, qp = qp)
  }
  sgn.out <- if (positive) { s0 <- sign(mean(ybar)); if (s0 == 0) 1 else s0 } else 0
  abc <- atpk.points.to.abc(z, ss, Ns, n, h, ybar, sgn = sgn.out)
  vpc <- vapply(1:n, function(i) mean(v[((i - 1) * Ns + 1):(i * Ns)]), numeric(1))
  list(a = abc$a, b = abc$b, c = abc$c, var = vpc, dormant = FALSE)
}
