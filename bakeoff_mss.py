#!/usr/bin/env python3
"""bakeoff_mss.py — three-way bake-off: PIQS vs PCHIP vs MSS.

MSS (monotone smoothing spline): cubic smoothing spline on the cumulative
integral, minimizing integral(F'')^2 subject to F(t_k) = F_k (integral
preservation) and f = F' >= 0 everywhere (monotone cumulative -> non-
negative flux). Solved as a QP per cell on the knot slopes m_k via
scipy.optimize.minimize(method="trust-constr").

The flux is a piecewise quadratic — same storage layout as PIQS/PCHIP,
drop-in for diurnalize-ERA5.r.

Run from MiCASA_v2/ working directory after cat_monthly.sh.
"""
import sys
import time
import numpy as np
import netCDF4 as nc
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize, LinearConstraint
from scipy.sparse import csr_matrix, lil_matrix

MONTHLY = "monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc"

# ---- PIQS port (Rasmussen 1991 S2 smoothness) ------------------------------
def piqs(x, ybar):
    x = np.asarray(x, dtype=np.float64)
    ybar = np.asarray(ybar, dtype=np.float64)
    n = len(x) - 1
    delta = np.diff(x)
    A = np.zeros((n - 1, n - 1))
    for i in range(n - 1):
        if i + 1 < n - 1: A[i, i + 1] = delta[i]
        if i > 0:         A[i, i - 1] = delta[i + 1]
        A[i, i] = 2 * (delta[i] + delta[i + 1])
    g1 = np.zeros(n - 1); g1[0] = -delta[1]
    g2 = np.zeros(n - 1); g2[-1] = -delta[n - 2]
    g3 = 3 * (delta[:n-1] * ybar[1:n] + delta[1:n] * ybar[:n-1])
    for i in range(1, n - 1):
        z = A[i, i-1] / A[i-1, i-1]
        A[i, i] -= z * A[i-1, i]
        g1[i] -= z * g1[i-1]
        g3[i] -= z * g3[i-1]
    f1 = np.empty(n+1); f2 = np.empty(n+1); f3 = np.empty(n+1)
    f1[0]=1; f2[0]=0; f3[0]=0; f1[n]=0; f2[n]=1; f3[n]=0
    f1[n-1] = g1[n-2] / A[n-2, n-2]
    f2[n-1] = g2[n-2] / A[n-2, n-2]
    f3[n-1] = g3[n-2] / A[n-2, n-2]
    for i in range(n-2, 0, -1):
        f1[i] = (g1[i-1] - A[i-1, i] * f1[i+1]) / A[i-1, i-1]
        f2[i] = (g2[i-1] - A[i-1, i] * f2[i+1]) / A[i-1, i-1]
        f3[i] = (g3[i-1] - A[i-1, i] * f3[i+1]) / A[i-1, i-1]
    inv_d3 = 1.0 / delta**3
    sf1 = f1[:n] + f1[1:n+1]; sf2 = f2[:n] + f2[1:n+1]; sf3 = f3[:n] + f3[1:n+1]
    t1 = np.sum(inv_d3 * sf1 * (2*ybar - sf3))
    t2 = np.sum(inv_d3 * sf2 * (2*ybar - sf3))
    R11 = np.sum(inv_d3 * sf1 * sf1); R22 = np.sum(inv_d3 * sf2 * sf2)
    R12 = np.sum(inv_d3 * sf1 * sf2)
    z = R11 * R22 - R12 * R12
    y = np.empty(n + 1)
    y[0] = (R22 * t1 - R12 * t2) / z
    y[n] = (R11 * t2 - R12 * t1) / z
    for i in range(1, n):
        y[i] = f1[i] * y[0] + f2[i] * y[n] + f3[i]
    a = (3 * y[:n] + 3 * y[1:n+1] - 6 * ybar) / delta**2
    b = (-4 * y[:n] - 2 * y[1:n+1] + 6 * ybar) / delta
    c = y[:n].copy()
    return dict(a=a, b=b, c=c, y=y)

# ---- PCHIP-on-cumulative ---------------------------------------------------
def pchip_oncum(x, ybar):
    x = np.asarray(x, dtype=np.float64)
    ybar = np.asarray(ybar, dtype=np.float64)
    n = len(x) - 1
    delta = np.diff(x)
    F = np.concatenate(([0.0], np.cumsum(ybar * delta)))
    pp = PchipInterpolator(x, F, extrapolate=False)
    pp_d = pp.derivative()
    coefs = pp_d.c
    a = coefs[0].copy(); b = coefs[1].copy(); c = coefs[2].copy()
    y = pp_d(x)
    return dict(a=a, b=b, c=c, y=y)


# ---- Monotone Smoothing Spline (MSS) ---------------------------------------
def mss(x, ybar, n_test_per_segment=8, init_from_pchip=True):
    """Cubic smoothing spline on cumulative F minimizing int(F'')^2,
    subject to F(t_k) = F_k and f = F' >= 0. Solved as QP on knot slopes.

    Each segment's flux f(s) on s in [0,1]:
        f(s) = (6s - 6s^2) u_k + (3s^2 - 4s + 1) m_k + (3s^2 - 2s) m_{k+1}
    where u_k = (F_{k+1} - F_k) / h_k = ybar_k.

    Smoothness functional per segment (analytic):
        S_k = (1/h_k) [4 m_k^2 + 4 m_{k+1}^2 + 4 m_k m_{k+1}
                       + 12 u_k m_k + 12 u_k m_{k+1} + 12 u_k^2]
    """
    x = np.asarray(x, dtype=np.float64)
    ybar = np.asarray(ybar, dtype=np.float64)
    n = len(x) - 1
    h = np.diff(x)               # length n
    u = ybar.copy()              # length n; the secant slopes (= monthly means)

    # Build symmetric Hessian H (size (n+1) x (n+1)) and gradient g for
    # objective 0.5 m^T H m + g^T m.
    # Per-segment contributions:
    #   m_k^2   coef = 4/h_k     -> diag[k] += 8/h_k  (factor 2 to convert to 0.5 m^T H m form)
    #   m_{k+1}^2 coef = 4/h_k   -> diag[k+1] += 8/h_k
    #   m_k m_{k+1} coef = 4/h_k -> H[k, k+1] += 4/h_k symmetric
    #   m_k linear coef = 12 u_k / h_k
    #   m_{k+1} linear coef = 12 u_k / h_k
    H = np.zeros((n + 1, n + 1))
    g = np.zeros(n + 1)
    for k in range(n):
        inv_h = 1.0 / h[k]
        H[k, k]       += 8 * inv_h
        H[k+1, k+1]   += 8 * inv_h
        H[k, k+1]     += 4 * inv_h
        H[k+1, k]     += 4 * inv_h
        # Linear cross term: 2 * <A, B'> * u * m = 2 * (-6) * u * m = -12 u m
        # where A(s) = (6-12s) and B'(s) = (6s-4), C(s) = (6s-2). The cross
        # integrals are <A,B'> = <A,C> = -6, <B',C> = +2. Was +12 (sign bug).
        g[k]   += -12 * u[k] * inv_h
        g[k+1] += -12 * u[k] * inv_h
    # objective(m) = 0.5 m^T H m + g^T m

    # Linear constraints: f >= 0 at test points within each segment, plus
    # endpoints m_0, m_n >= 0 (which equal f at t_0 and t_n).
    # Test points: s_j = (j + 0.5) / n_test for j=0..n_test-1
    test_s = (np.arange(n_test_per_segment) + 0.5) / n_test_per_segment
    # Each constraint row: f(s, k) >= 0
    #   = (6s - 6s^2) u_k + (3s^2 - 4s + 1) m_k + (3s^2 - 2s) m_{k+1}
    # Solver expects A m >= lb, with lb_j = -(constant_term_j)
    # i.e. (3s^2 - 4s + 1) m_k + (3s^2 - 2s) m_{k+1} >= -(6s - 6s^2) u_k
    n_tests = n * n_test_per_segment + 2  # +2 for boundary (m_0 >= 0, m_n >= 0)
    A = lil_matrix((n_tests, n + 1))
    lb = np.zeros(n_tests)
    row = 0
    for k in range(n):
        for s in test_s:
            coef_mk   = 3*s*s - 4*s + 1
            coef_mk1  = 3*s*s - 2*s
            const_k   = (6*s - 6*s*s) * u[k]
            A[row, k]   = coef_mk
            A[row, k+1] = coef_mk1
            lb[row] = -const_k   # so A m >= lb  <=>  flux >= 0
            row += 1
    # Boundary endpoints
    A[row, 0] = 1.0; lb[row] = 0.0; row += 1
    A[row, n] = 1.0; lb[row] = 0.0; row += 1

    A_csr = csr_matrix(A)
    constraints = LinearConstraint(A_csr, lb=lb, ub=np.inf)

    # Initial guess: PCHIP slopes (always feasible since PCHIP gives f >= 0)
    if init_from_pchip:
        m0 = pchip_oncum(x, ybar)["y"]
    else:
        m0 = np.maximum(0.0, np.zeros(n + 1))  # all-zero fallback

    # Objective and gradient (for trust-constr)
    def fun(m): return 0.5 * m @ H @ m + g @ m
    def jac(m): return H @ m + g
    def hess(m): return H

    res = minimize(
        fun, m0, jac=jac, hess=hess,
        constraints=constraints,
        method="trust-constr",
        options=dict(maxiter=200, verbose=0, gtol=1e-10, xtol=1e-12),
    )
    if not res.success:
        # Fall back to PCHIP if QP didn't converge cleanly
        return pchip_oncum(x, ybar)

    m = res.x
    # Convert back to per-segment piecewise-quadratic (a, b, c) with
    # f(t) on [t_k, t_{k+1}] = a*(t-t_k)^2 + b*(t-t_k) + c
    # f(s) = (6s - 6s^2) u_k + (3s^2 - 4s + 1) m_k + (3s^2 - 2s) m_{k+1}
    # In s, f = -6 (u_k - 0) s^2 + ... let me expand:
    #   = u_k (6s - 6s^2) + m_k (3s^2 - 4s + 1) + m_{k+1} (3s^2 - 2s)
    #   = s^2 (-6 u_k + 3 m_k + 3 m_{k+1}) + s (6 u_k - 4 m_k - 2 m_{k+1}) + m_k
    # Convert to t-coordinate: t = t_k + s h_k => s = (t-t_k)/h_k
    # f(t) = a (t-t_k)^2 + b (t-t_k) + c
    #      = (s^2 coef) / h_k^2 * (t-t_k)^2 + (s coef) / h_k * (t-t_k) + (const)
    a = np.empty(n); b = np.empty(n); c = np.empty(n)
    for k in range(n):
        Q = -6 * u[k] + 3 * m[k] + 3 * m[k+1]
        L = 6 * u[k] - 4 * m[k] - 2 * m[k+1]
        K = m[k]
        a[k] = Q / (h[k] ** 2)
        b[k] = L / h[k]
        c[k] = K
    return dict(a=a, b=b, c=c, y=m)


# ---- Sampling and metrics --------------------------------------------------
def sample_hourly(fit, x_knot, t_hourly):
    a, b, c = fit["a"], fit["b"], fit["c"]
    seg = np.searchsorted(x_knot, t_hourly, side="right") - 1
    seg = np.clip(seg, 0, len(a) - 1)
    dt = t_hourly - x_knot[seg]
    return a[seg]*dt*dt + b[seg]*dt + c[seg]

def metrics(name, fit, x_knot, ybar, t_hourly, sign_target):
    f_h = sample_hourly(fit, x_knot, t_hourly)
    flips = (f_h * sign_target) < -1e-15
    flip_pct = 100.0 * flips.sum() / len(f_h)
    a, b, c = fit["a"], fit["b"], fit["c"]
    delta = np.diff(x_knot)
    integral = a/3 * delta**3 + b/2 * delta**2 + c * delta
    expected = ybar * delta
    int_max_rel = np.max(np.abs((integral - expected) / np.where(np.abs(expected) > 1e-30, expected, 1)))
    df_left  = np.abs(b)
    df_right = np.abs(2 * a * delta + b)
    smooth   = np.max(np.maximum(df_left, df_right))
    overshoot = 0.0
    for i in range(len(a)):
        if a[i] != 0:
            tv = -b[i] / (2 * a[i])
            if 0 < tv < delta[i]:
                fv = a[i]*tv*tv + b[i]*tv + c[i]
                if sign_target > 0 and fv < 0:
                    overshoot = max(overshoot, -fv)
                elif sign_target < 0 and fv > 0:
                    overshoot = max(overshoot, fv)
    # Sub-monthly smoothness as integrated |F''|^2 across record
    # F''(t) on segment k = 2 a_k * (t-t_k) + b_k for the FLUX f = F'.
    # Wait — f = F' so f' = F''. f' = 2 a (t-t_k) + b. Integrate (f')^2:
    # int_0^delta (2a*tau + b)^2 dtau = (4a^2/3)*delta^3 + 2a*b*delta^2 + b^2*delta
    smooth_int = float(np.sum(4*a*a/3 * delta**3 + 2*a*b * delta**2 + b*b * delta))
    return dict(name=name, flip_pct=flip_pct, int_max_rel=int_max_rel,
                smooth=smooth, smooth_int=smooth_int, overshoot=overshoot,
                fmin=f_h.min(), fmax=f_h.max())

def main():
    print(f"Loading {MONTHLY}...")
    ds = nc.Dataset(MONTHLY)
    npp = ds.variables["NPP"][:]
    rh  = ds.variables["Rh"][:]
    times = ds.variables["time"][:]
    lats = ds.variables["latitude"][:]
    lons = ds.variables["longitude"][:]
    import datetime
    nmon = len(times)
    epoch = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
    t0 = epoch + datetime.timedelta(seconds=float(times[0]))
    y0, m0 = t0.year, t0.month
    knot_times = []
    for k in range(nmon + 1):
        m = m0 + k
        y = y0 + (m - 1) // 12
        m = (m - 1) % 12 + 1
        knot_times.append((datetime.datetime(y, m, 1, tzinfo=datetime.timezone.utc) - epoch).total_seconds())
    x_knot = np.array(knot_times, dtype=np.float64)
    dt_h = 3600.0
    t_h = np.arange(x_knot[0], x_knot[-1], dt_h)

    targets = [
        ("Manaus_evergreen",   -3.0,   -60.0),
        ("Hyytiala_boreal",    61.5,    24.0),
        ("Sahel_savanna",      15.0,     0.0),
        ("Polar_night_Arctic", 80.5,   100.0),
        ("Semi-arid_Texas",    33.0,  -100.0),
        ("Boreal_Tundra_AK",   68.0,  -150.0),
    ]
    print(f"Knots: n={nmon} segments; hourly samples: {len(t_h)}")
    print()
    rows = []
    for name, lat, lon in targets:
        ilat = int(np.argmin(np.abs(lats - lat)))
        ilon = int(np.argmin(np.abs(lons - lon)))
        npp_cell = npp[:, ilat, ilon].astype(np.float64)
        rh_cell  = rh [:, ilat, ilon].astype(np.float64)
        if not np.all(np.isfinite(npp_cell)) or not np.all(np.isfinite(rh_cell)):
            print(f"-- {name} ({lat:+.1f}, {lon:+.1f}): skip (NaN)"); continue
        if np.max(np.abs(npp_cell)) < 1e-15 and np.max(np.abs(rh_cell)) < 1e-15:
            print(f"-- {name} ({lat:+.1f}, {lon:+.1f}): skip (zero)"); continue
        print(f"=== {name} ({lat:+.1f}, {lon:+.1f}) [grid {ilat},{ilon}] ===")
        for var_name, ybar, sign_target in [("NPP", npp_cell, +1), ("Rh", rh_cell, +1)]:
            try:    pf = piqs(x_knot, ybar)
            except Exception as e:
                print(f"   {var_name} PIQS FAILED: {e}"); continue
            ph = pchip_oncum(x_knot, ybar)
            t0 = time.perf_counter()
            pm = mss(x_knot, ybar, n_test_per_segment=8)
            t_mss = time.perf_counter() - t0
            mp = metrics("PIQS",  pf, x_knot, ybar, t_h, sign_target)
            mc = metrics("PCHIP", ph, x_knot, ybar, t_h, sign_target)
            mm = metrics("MSS",   pm, x_knot, ybar, t_h, sign_target)
            rows.append((name, var_name, mp, mc, mm))
            for m in (mp, mc, mm):
                print(f"   {var_name:3s} {m['name']:5s}: flip%={m['flip_pct']:6.3f}  int_rel={m['int_max_rel']:.2e}  "
                      f"overshoot={m['overshoot']:.2e}  max|df|={m['smooth']:.2e}  "
                      f"int(F'')^2={m['smooth_int']:.3e}  range=[{m['fmin']:.2e},{m['fmax']:.2e}]")
            print(f"   ({var_name:3s} MSS solve time: {t_mss*1000:.0f} ms)")
        print()
    print("=" * 100)
    print("Aggregate (median across cells × NPP/Rh):")
    print("=" * 100)
    for label, key in [("flip%", "flip_pct"), ("int_max_rel", "int_max_rel"),
                       ("max_overshoot", "overshoot"), ("max|df|", "smooth"),
                       ("int(F'')^2", "smooth_int")]:
        piqs_vals  = [r[2][key] for r in rows]
        pchip_vals = [r[3][key] for r in rows]
        mss_vals   = [r[4][key] for r in rows]
        print(f"  {label:14s}: PIQS={np.median(piqs_vals):.3e}  PCHIP={np.median(pchip_vals):.3e}  MSS={np.median(mss_vals):.3e}")

if __name__ == "__main__":
    main()
