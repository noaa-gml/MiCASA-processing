#!/usr/bin/env python3
"""bakeoff_pchip.py — PIQS vs PCHIP-on-cumulative spline fitter comparison.

Targeted bake-off on a handful of representative 1-degree cells.
Reports sign-flip rate, integral preservation, mass conservation,
sub-monthly smoothness, and polar-night residual for both methods,
sampled at hourly resolution over the full multi-year record.

Run from MiCASA_v2/ working directory after cat_monthly.sh.
"""
import sys
import numpy as np
import netCDF4 as nc
from scipy.interpolate import PchipInterpolator

MONTHLY = "monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc"

# ---- PIQS port (Rasmussen 1991 S2 smoothness, from piqs.r.txt) -------------
def piqs(x, ybar):
    """Piecewise integral quadratic spline.

    x    : knot times, length n+1
    ybar : segment means,  length n
    returns dict(a, b, c, y) where each segment is
        f_i(t) = a[i]*(t - x[i])**2 + b[i]*(t - x[i]) + c[i]
    """
    x    = np.asarray(x, dtype=np.float64)
    ybar = np.asarray(ybar, dtype=np.float64)
    n    = len(x) - 1
    if n < 3:
        raise ValueError("piqs requires at least 3 segments")
    delta = np.diff(x)

    # Tridiagonal A (n-1 by n-1)
    A = np.zeros((n - 1, n - 1))
    for i in range(n - 1):
        if i + 1 < n - 1:
            A[i, i + 1] = delta[i]
        if i > 0:
            A[i, i - 1] = delta[i + 1]
        A[i, i] = 2 * (delta[i] + delta[i + 1])

    g1 = np.zeros(n - 1); g1[0] = -delta[1]
    g2 = np.zeros(n - 1); g2[-1] = -delta[n - 2]
    g3 = 3 * (delta[:n - 1] * ybar[1:n] + delta[1:n] * ybar[:n - 1])

    # Forward elimination (Thomas)
    for i in range(1, n - 1):
        z = A[i, i - 1] / A[i - 1, i - 1]
        A[i, i] -= z * A[i - 1, i]
        g1[i]   -= z * g1[i - 1]
        g3[i]   -= z * g3[i - 1]

    f1 = np.empty(n + 1); f2 = np.empty(n + 1); f3 = np.empty(n + 1)
    f1[0] = 1.0; f2[0] = 0.0; f3[0] = 0.0
    f1[n] = 0.0; f2[n] = 1.0; f3[n] = 0.0
    f1[n - 1] = g1[n - 2] / A[n - 2, n - 2]
    f2[n - 1] = g2[n - 2] / A[n - 2, n - 2]
    f3[n - 1] = g3[n - 2] / A[n - 2, n - 2]
    for i in range(n - 2, 0, -1):
        f1[i] = (g1[i - 1] - A[i - 1, i] * f1[i + 1]) / A[i - 1, i - 1]
        f2[i] = (g2[i - 1] - A[i - 1, i] * f2[i + 1]) / A[i - 1, i - 1]
        f3[i] = (g3[i - 1] - A[i - 1, i] * f3[i + 1]) / A[i - 1, i - 1]

    # Smoothness condition (Eq 15)
    inv_d3 = 1.0 / delta**3
    sum_f1 = f1[:n] + f1[1:n + 1]
    sum_f2 = f2[:n] + f2[1:n + 1]
    sum_f3 = f3[:n] + f3[1:n + 1]
    t1 = np.sum(inv_d3 * sum_f1 * (2 * ybar - sum_f3))
    t2 = np.sum(inv_d3 * sum_f2 * (2 * ybar - sum_f3))
    R11 = np.sum(inv_d3 * sum_f1 * sum_f1)
    R22 = np.sum(inv_d3 * sum_f2 * sum_f2)
    R12 = np.sum(inv_d3 * sum_f1 * sum_f2)
    z = R11 * R22 - R12 * R12
    y = np.empty(n + 1)
    y[0] = (R22 * t1 - R12 * t2) / z
    y[n] = (R11 * t2 - R12 * t1) / z
    for i in range(1, n):
        y[i] = f1[i] * y[0] + f2[i] * y[n] + f3[i]

    a = (3 * y[:n] + 3 * y[1:n + 1] - 6 * ybar) / delta**2
    b = (-4 * y[:n] - 2 * y[1:n + 1] + 6 * ybar) / delta
    c = y[:n].copy()
    return dict(a=a, b=b, c=c, y=y)


# ---- PCHIP-on-cumulative ---------------------------------------------------
def pchip_oncum(x, ybar):
    """Monotone-cubic Hermite on the cumulative integral, differentiated.

    Returns the same a/b/c piecewise-quadratic structure as piqs(),
    where f_i(t) = a[i]*(t - x[i])**2 + b[i]*(t - x[i]) + c[i].
    """
    x    = np.asarray(x, dtype=np.float64)
    ybar = np.asarray(ybar, dtype=np.float64)
    n    = len(x) - 1
    delta = np.diff(x)

    # Cumulative integral at knots: F(x[0])=0, F(x[k+1]) = F(x[k]) + ybar[k]*delta[k]
    F = np.concatenate(([0.0], np.cumsum(ybar * delta)))

    # Fritsch-Carlson via scipy. PchipInterpolator gives a piecewise cubic.
    pp = PchipInterpolator(x, F, extrapolate=False)
    pp_d = pp.derivative()  # piecewise quadratic, c shape (3, n)
    # scipy PPoly convention: f(t) on segment k = sum c[i,k]*(t - x[k])**(deg - i)
    # For a degree-2 derivative, c has shape (3, n) -> a=c[0], b=c[1], c=c[2]
    coefs = pp_d.c  # shape (3, n)
    a = coefs[0].copy()
    b = coefs[1].copy()
    c = coefs[2].copy()
    y = pp_d(x)     # endpoint flux values; len n+1
    return dict(a=a, b=b, c=c, y=y)


def sample_hourly(fit, x_knot, t_hourly):
    """Sample a piecewise-quadratic fit at hourly times."""
    a, b, c = fit["a"], fit["b"], fit["c"]
    out = np.full_like(t_hourly, np.nan, dtype=np.float64)
    seg = np.searchsorted(x_knot, t_hourly, side="right") - 1
    seg = np.clip(seg, 0, len(a) - 1)
    dt = t_hourly - x_knot[seg]
    out = a[seg] * dt * dt + b[seg] * dt + c[seg]
    return out


def metrics(name, fit, x_knot, ybar, t_hourly, sign_target):
    """sign_target = +1 (expect non-negative) or -1 (expect non-positive)."""
    f_h = sample_hourly(fit, x_knot, t_hourly)

    # 1. Sign-flip rate: cell-hours with f * sign_target < 0
    flips = (f_h * sign_target) < -1e-15
    flip_pct = 100.0 * flips.sum() / len(f_h)

    # 2. Integral preservation per piece (analytic)
    a, b, c = fit["a"], fit["b"], fit["c"]
    delta = np.diff(x_knot)
    integral = a/3 * delta**3 + b/2 * delta**2 + c * delta
    expected = ybar * delta
    int_max_rel = np.max(np.abs((integral - expected) / np.where(np.abs(expected) > 1e-30, expected, 1)))

    # 3. Sub-monthly smoothness: max |df/dt| per piece
    # f' = 2 a (t-x_i) + b. Max within piece is at t-x_i = 0 or delta -> |b| or |2 a delta + b|
    df_left  = np.abs(b)
    df_right = np.abs(2 * a * delta + b)
    smooth   = np.max(np.maximum(df_left, df_right))

    # 4. Within-piece extremum (vertex): vertex at t-x_i = -b/(2a). Value at vertex.
    #    Compute amount by which fit overshoots the "wrong" sign.
    overshoot = 0.0
    for i in range(len(a)):
        if a[i] != 0:
            tv = -b[i] / (2 * a[i])
            if 0 < tv < delta[i]:
                fv = a[i] * tv * tv + b[i] * tv + c[i]
                # how far above zero (for negative-target) or below (for positive-target)
                if sign_target > 0 and fv < 0:
                    overshoot = max(overshoot, -fv)
                elif sign_target < 0 and fv > 0:
                    overshoot = max(overshoot, fv)

    return dict(name=name, flip_pct=flip_pct, int_max_rel=int_max_rel,
                smooth=smooth, overshoot=overshoot,
                fmin=f_h.min(), fmax=f_h.max())


def main():
    print(f"Loading {MONTHLY}...")
    ds = nc.Dataset(MONTHLY)
    npp  = ds.variables["NPP"][:]   # (time, lat, lon), gC m-2 s-1
    rh   = ds.variables["Rh"][:]
    times = ds.variables["time"][:] # seconds since epoch
    lats = ds.variables["latitude"][:]
    lons = ds.variables["longitude"][:]

    # Knot times: n+1 endpoints. The cat'd file has time at month centers.
    # Build knot times as month-start; same convention as write_piqs.r:
    #   x.time = seq(start_of_record, by="1 month", length.out=nmon+1)
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

    # Hourly times across the whole record
    dt_h = 3600.0
    t_h = np.arange(x_knot[0], x_knot[-1], dt_h)

    # Representative cells: (name, lat, lon, biome notes)
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
        npp_cell = npp[:, ilat, ilon].astype(np.float64)  # gC m-2 s-1
        rh_cell  = rh [:, ilat, ilon].astype(np.float64)

        # GPP convention: gpp = -2 * NPP / 12 (mol m-2 s-1) — but for the
        # bake-off we compare on the underlying NPP & Rh time series since
        # that's what PIQS/PCHIP fits operate on. Use gC m-2 s-1.
        if not np.all(np.isfinite(npp_cell)) or not np.all(np.isfinite(rh_cell)):
            print(f"-- {name} ({lat:+.1f}, {lon:+.1f}): skip (NaN in input)")
            continue

        # Skip cells where both fluxes are essentially zero (water/ice)
        if np.max(np.abs(npp_cell)) < 1e-15 and np.max(np.abs(rh_cell)) < 1e-15:
            print(f"-- {name} ({lat:+.1f}, {lon:+.1f}): skip (cell at machine zero, ocean/ice)")
            continue

        print(f"=== {name} ({lat:+.1f}, {lon:+.1f}) [grid {ilat},{ilon}] ===")
        print(f"   NPP range [{npp_cell.min():.3e}, {npp_cell.max():.3e}] gC m-2 s-1")
        print(f"   Rh  range [{rh_cell.min():.3e}, {rh_cell.max():.3e}] gC m-2 s-1")

        for var_name, ybar, sign_target in [
            ("NPP", npp_cell, +1),  # NPP is non-negative
            ("Rh",  rh_cell,  +1),  # Rh is non-negative
        ]:
            try:
                pf = piqs(x_knot, ybar)
            except Exception as e:
                print(f"   {var_name} PIQS FAILED: {e}")
                continue
            ph = pchip_oncum(x_knot, ybar)
            mp = metrics("PIQS",  pf, x_knot, ybar, t_h, sign_target)
            mc = metrics("PCHIP", ph, x_knot, ybar, t_h, sign_target)
            rows.append((name, var_name, mp, mc))
            for m in (mp, mc):
                print(f"   {var_name:3s} {m['name']:5s}: flip%={m['flip_pct']:6.3f}  int_rel={m['int_max_rel']:.2e}  "
                      f"max_overshoot={m['overshoot']:.2e}  smooth(max|df|)={m['smooth']:.2e}  "
                      f"fmin={m['fmin']:.2e}  fmax={m['fmax']:.2e}")
        print()

    # Aggregate summary
    print("=" * 88)
    print("Aggregate (median across cells × NPP/Rh):")
    print("=" * 88)
    for label, key in [("flip%", "flip_pct"),
                       ("int_max_rel", "int_max_rel"),
                       ("max_overshoot", "overshoot"),
                       ("smooth", "smooth")]:
        piqs_vals  = [r[2][key] for r in rows]
        pchip_vals = [r[3][key] for r in rows]
        print(f"  {label:14s}: PIQS median={np.median(piqs_vals):.3e}  PCHIP median={np.median(pchip_vals):.3e}")

if __name__ == "__main__":
    main()
