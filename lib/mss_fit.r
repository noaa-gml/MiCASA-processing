## lib/mss_fit.r -- monotone smoothing-spline (MSS) fitter core.
##
## The QP-based alternative fitter used by write_mss.r: a cubic spline on
## the cumulative integral that minimizes int(F'')^2 subject to integral
## preservation F(t_k) = F_k and non-negativity f = F' >= 0 (enforced at
## test points), solved per cell as a quadratic program on the knot slopes.
##
## mss.fit.setup(x) precomputes the smoothness Hessian and constraint
## matrices -- they depend only on the knot positions x, not on the data --
## so write_mss.r builds them once and reuses them across all 360x180
## cells. mss.fit.cell(x, ybar, setup) then fits one cell.
##
## Extracted from write_mss.r so the fitter can be unit-tested standalone
## (tests/test_mss_fit.r). Requires the `quadprog` package (solve.QP).

## Precompute the data-independent QP structures for knot positions `x`.
mss.fit.setup <- function(x, n.test.per.segment = 8) {
  h.seg <- diff(x)
  n.seg <- length(h.seg)
  n.var <- n.seg + 1

  ## Smoothness Hessian G (the integral of (F'')^2 as a quadratic form in
  ## the knot slopes), plus a tiny ridge so solve.QP sees a PD matrix.
  G <- matrix(0, n.var, n.var)
  for (k in seq_len(n.seg)) {
    inv.h <- 1 / h.seg[k]
    G[k,     k]     <- G[k,     k]     + 8 * inv.h
    G[k + 1, k + 1] <- G[k + 1, k + 1] + 8 * inv.h
    G[k,     k + 1] <- G[k,     k + 1] + 4 * inv.h
    G[k + 1, k]     <- G[k + 1, k]     + 4 * inv.h
  }
  G <- G + 1e-20 * diag(n.var)

  ## Constraint matrix C: per-segment test-point coefficients on
  ## (m_k, m_{k+1}) plus two boundary rows. b.coef maps the per-segment
  ## u_k onto the constraint right-hand side.
  n.constraints <- n.seg * n.test.per.segment + 2
  C      <- matrix(0, n.var, n.constraints)
  b.coef <- matrix(0, n.constraints, n.seg)
  test.s <- (seq_len(n.test.per.segment) - 0.5) / n.test.per.segment
  j <- 0
  for (k in seq_len(n.seg)) {
    for (s in test.s) {
      j <- j + 1
      C[k,     j]  <- 3 * s * s - 4 * s + 1
      C[k + 1, j]  <- 3 * s * s - 2 * s
      b.coef[j, k] <- -(6 * s - 6 * s * s)
    }
  }
  j <- j + 1; C[1,     j] <- 1
  j <- j + 1; C[n.var, j] <- 1

  list(h.seg = h.seg, n.seg = n.seg, n.var = n.var,
       G = G, C = C, b.coef = b.coef,
       n.test.per.segment = n.test.per.segment)
}

## Fit one cell: knot times `x`, monthly means `ybar`, and the precomputed
## `setup` from mss.fit.setup(x). Returns list(a, b, c) of per-piece
## quadratic coefficients, the same layout as pchip.fit.cell / write_piqs.r.
mss.fit.cell <- function(x, ybar, setup) {
  n.seg <- setup$n.seg
  n.var <- setup$n.var
  h.seg <- setup$h.seg
  if (all(ybar == 0)) {
    return(list(a = rep(0, n.seg), b = rep(0, n.seg), c = rep(0, n.seg)))
  }
  ## Uniformly-negative data (the GPP case): flip sign, fit, unflip the
  ## resulting coefficients.
  if (mean(ybar) < 0) {
    u <- -ybar
    sign.flip <- -1
  } else {
    u <- ybar
    sign.flip <- 1
  }
  inv.h <- 1 / h.seg
  ## Linear term of 0.5 m^T G m - a^T m:  a[k] = 12 u[k] / h[k].
  a.lin <- numeric(n.var)
  for (k in seq_len(n.seg)) {
    a.lin[k]     <- a.lin[k]     + 12 * u[k] * inv.h[k]
    a.lin[k + 1] <- a.lin[k + 1] + 12 * u[k] * inv.h[k]
  }
  bvec <- as.vector(setup$b.coef %*% u)
  sol <- tryCatch(
    quadprog::solve.QP(Dmat = setup$G, dvec = a.lin,
                       Amat = setup$C, bvec = bvec, meq = 0),
    error = function(e) NULL)
  if (is.null(sol)) {
    ## QP failed (very pathological input only) -- fall back to PCHIP.
    fn <- splinefun(x, c(0, cumsum(u * h.seg)), method = "monoH.FC")
    m  <- fn(x, deriv = 1)
  } else {
    m <- sol$solution
  }
  Q <- -6 * u + 3 * m[1:n.seg] + 3 * m[2:n.var]
  L <-  6 * u - 4 * m[1:n.seg] - 2 * m[2:n.var]
  K <- m[1:n.seg]
  list(a = sign.flip * Q / h.seg^2,
       b = sign.flip * L / h.seg,
       c = sign.flip * K)
}
