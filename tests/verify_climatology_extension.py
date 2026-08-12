#!/usr/bin/env python3
"""verify_climatology_extension.py -- check months appended to an existing
climatology-prior tree.

    python3 tests/verify_climatology_extension.py --outdir DIR \
            --year 2026 --months 1-3 --against 2025

tests/verify_climatology_prior.py iterates `for m in range(1, 13)` and so can
only verify whole years. An extension is by nature a PARTIAL year, and the
question it raises is not the one the full battery asks. The full battery asks
"is this tree a climatology prior?"; an extension asks "are these new months
the SAME product as the ones already delivered?"

That question has a sharp answer, because the series is climatological and the
fit is interior at both ends of the compared months:

  MEAN      the monthly mean of the new month must equal the monthly mean of
            the SAME CALENDAR MONTH in the reference year, cell by cell. The
            delivered monthly mean comes from the PCHIP fit, and the fit was
            built once over a periodic 2020-01..2026-12 series -- so a refit, a
            padding change, a trend leak or a wrong series file all break this
            and nothing else catches them. (Real meteorology redistributes flux
            WITHIN the month; it cannot move the mean.)

  LIVE      the HOURLY fields must nevertheless differ everywhere on land.
            Identical hourly fields would mean the year's real meteorology was
            never read -- the failure mode where a month is silently copied.

  FALSIFY   MEAN is compared against the WRONG calendar month as well. A check
            that cannot fail is not a check; if Jan-2026-vs-Feb-2025 also
            "passes", the comparison is measuring nothing.

## Why MEAN is not an exact-equality test, and how the threshold was set

The two months carry the same monthly mean ANALYTICALLY, but they are computed
from different hourly weights (different real meteorology) and stored as
float32, so the two sums differ in their last bits. The tolerance therefore has
to come from float32 resolution, not from a number chosen to make the test
pass:

    max|d|  <=  ULP_TOL * eps32 * max|reference|

Measured 2026 Q1 vs 2025 Q1 with the accumulation done in float64: **0.17 ULP**
at the peak value, with the global monthly integral agreeing to 7e-9 PgC/yr out
of 6.77. The wrong-month comparison differs by 2.4e-6 mol m-2 s-1 -- about
5e7 times larger. ULP_TOL = 10 sits four orders of magnitude below anything a
real defect produces, so the gate discriminates by an enormous margin rather
than by a tuned threshold.

⚠ Accumulate the monthly mean in FLOAT64. A `nanmean` over the float32 stack
makes the CHECKER the dominant error term -- it reported 3.4e-12 where the
honest answer is 4.8e-14, a 70x inflation, which is enough to make a clean
product look like a failure.

Structural gates (file count, 24-hour axes, zero-byte files) are already
enforced by climatology_daysplit.sbatch's post-gate; they are re-run here
because a gate that only ever runs inside the job that produces the artefact
cannot be re-run against the artefact later.
"""

import argparse
import calendar
import os
import sys

import numpy as np
from netCDF4 import Dataset

SEC_PER_YR = 365.0 * 86400.0
GC_PER_MOL = 12.0107
EPS32 = float(np.finfo(np.float32).eps)


def parse_months(s):
    if "-" in s:
        a, b = s.split("-", 1)
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def daily_paths(outdir, ver, yr, mon):
    n = calendar.monthrange(yr, mon)[1]
    return [f"{outdir}/ERA5/MiCASA_{ver}.nee.{yr}{mon:02d}{d:02d}.nc"
            for d in range(1, n + 1)]


def exists(p):
    return os.path.exists(p) and os.path.getsize(p) > 0


def read_month(paths, want_stack=True):
    """Return (monthly mean [lat,lon] in float64, hourly stack or None, day-1 attrs)."""
    acc = None
    nslot = 0
    hours = []
    attrs = None
    for p in paths:
        with Dataset(p) as nc:
            v = nc.variables["NEE"][:]
            if attrs is None:
                attrs = {k: nc.getncattr(k) for k in nc.ncattrs()}
            if v.shape[0] != 24:
                raise SystemExit(f"FATAL: {p} has {v.shape[0]} time slots, not 24")
            f = np.ma.filled(v, np.nan)
            # float64 accumulation: the checker must not be the error term.
            s = f.astype(np.float64).sum(axis=0)
            acc = s if acc is None else acc + s
            nslot += f.shape[0]
            if want_stack:
                hours.append(f)
    return acc / nslot, (np.concatenate(hours, axis=0) if want_stack else None), attrs


def area_grid(nlat=180, nlon=360):
    """Cell area [m2] on a 1-degree grid, latitude centres -89.5..89.5."""
    R = 6.371e6
    edges = np.deg2rad(np.linspace(-90.0, 90.0, nlat + 1))
    band = 2.0 * np.pi * R * R * (np.sin(edges[1:]) - np.sin(edges[:-1])) / nlon
    return np.repeat(band[:, None], nlon, axis=1)


def pgc_per_yr(mean_field, area):
    """mol m-2 s-1 -> PgC/yr, over the whole grid."""
    good = np.isfinite(mean_field)
    return float(np.nansum(mean_field[good] * area[good]) * GC_PER_MOL * SEC_PER_YR * 1e-15)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--months", default="1-3")
    ap.add_argument("--against", type=int, required=True,
                    help="reference year already verified in this tree")
    ap.add_argument("--version", default="v1")
    ap.add_argument("--ulp-tol", type=float, default=10.0,
                    help="MEAN tolerance in float32 ULP at the peak value (default 10)")
    a = ap.parse_args()

    months = parse_months(a.months)
    area = area_grid()
    fails = []

    def ok(label, cond, detail=""):
        print(f"  {label:<52} {'PASS' if cond else 'FAIL'}" + (f"   {detail}" if detail else ""))
        if not cond:
            fails.append(label)

    print(f"\nclimatology extension: {a.year} months {a.months}  vs reference year {a.against}")
    print(f"tree: {a.outdir}")
    print(f"MEAN tolerance: {a.ulp_tol:g} float32 ULP at the peak value\n")

    def mean_gate(label, new_mean, ref_mean, tag):
        d = np.abs(new_mean - ref_mean)
        worst = float(np.nanmax(d))
        peak = float(np.nanmax(np.abs(ref_mean)))
        ulp = worst / (EPS32 * peak) if peak > 0 else np.inf
        ok(label, ulp <= a.ulp_tol,
           f"max|d| = {worst:.3e} = {ulp:.2f} ULP   {tag}")
        return ulp

    for mon in months:
        new_paths = daily_paths(a.outdir, a.version, a.year, mon)
        ref_paths = daily_paths(a.outdir, a.version, a.against, mon)
        ndays = calendar.monthrange(a.year, mon)[1]
        print(f"--- {a.year}-{mon:02d} ({ndays} days) ---")

        missing = [p for p in new_paths if not exists(p)]
        ok(f"COUNT    {ndays} daily files present", not missing,
           "" if not missing else f"missing {len(missing)}")
        if missing:
            continue

        new_mean, new_stack, at = read_month(new_paths)

        # -- STAMP: the product attributes a downstream reader will trust ----
        ok("STAMP    flux_from_climatology=yes",
           at.get("flux_from_climatology") == "yes", str(at.get("flux_from_climatology")))
        ok("STAMP    climatology_baseline=2001-2020",
           at.get("climatology_baseline") == "2001-2020", str(at.get("climatology_baseline")))
        ok("STAMP    resp driver/tempfun = airtemp/q10",
           (at.get("respiration_temperature_driver"),
            at.get("respiration_temperature_function")) == ("airtemp", "q10"),
           f"{at.get('respiration_temperature_driver')}/{at.get('respiration_temperature_function')}")

        # -- METEO: disclosure, not a pass/fail ------------------------------
        print(f"  {'METEO    source by day':<52} {at.get('meteo_source_by_day')}"
              f"   (fallback used: {at.get('meteo_fallback_used')})")

        if not all(exists(p) for p in ref_paths):
            print(f"  (no {a.against}-{mon:02d} in this tree; skipping MEAN/LIVE)")
            continue
        ref_mean, ref_stack, _ = read_month(ref_paths)

        # -- MEAN: the load-bearing identity ---------------------------------
        mean_gate(f"MEAN     == {a.against}-{mon:02d} monthly mean, cell by cell",
                  new_mean, ref_mean, "")

        # -- LIVE: real meteorology actually differs -------------------------
        land = np.isfinite(new_mean) & (np.abs(new_mean) > 0)
        if new_stack.shape == ref_stack.shape:
            hd = np.abs(new_stack - ref_stack)
            frac = float(np.nanmean(hd[:, land] > 0)) if land.any() else 0.0
            hmax = float(np.nanmax(hd))
            ok("LIVE     hourly fields differ from the reference year",
               frac > 0.99 and hmax > 0,
               f"{100*frac:.2f}% of land cell-hours differ, max|d| = {hmax:.3e}")
        else:
            ok("LIVE     hourly fields differ from the reference year", False,
               f"shape {new_stack.shape} vs {ref_stack.shape}")

        # -- BUDGET: the same identity, integrated ---------------------------
        gn, gr = pgc_per_yr(new_mean, area), pgc_per_yr(ref_mean, area)
        print(f"  {'BUDGET   global bio flux':<52} "
              f"{gn:+.9f} PgC/yr   ({a.against}: {gr:+.9f}, d = {gn-gr:+.2e})")

    # -- FALSIFY: the MEAN gate must be able to fail -------------------------
    print("\n--- FALSIFY (the MEAN gate must be able to fail) ---")
    m0 = months[0]
    alt = [m for m in range(1, 13)
           if m != m0 and all(exists(p) for p in daily_paths(a.outdir, a.version, a.against, m))]
    if alt:
        m1 = alt[0] if alt[0] != m0 else alt[1]
        new_mean, _, _ = read_month(daily_paths(a.outdir, a.version, a.year, m0),
                                    want_stack=False)
        wrong_mean, _, _ = read_month(daily_paths(a.outdir, a.version, a.against, m1),
                                      want_stack=False)
        d = np.abs(new_mean - wrong_mean)
        worst = float(np.nanmax(d))
        peak = float(np.nanmax(np.abs(wrong_mean)))
        ulp = worst / (EPS32 * peak)
        ok(f"FALSIFY  {a.year}-{m0:02d} vs {a.against}-{m1:02d} does NOT match",
           ulp > a.ulp_tol,
           f"max|d| = {worst:.3e} = {ulp:.3e} ULP (must exceed {a.ulp_tol:g})")
    else:
        ok("FALSIFY  a wrong-month comparison was available", False,
           "no other complete reference month; gate not exercised")

    print(f"\n{'OK' if not fails else 'FAILURES'}: {len(fails)} check(s) failed")
    for f in fails:
        print(f"  - {f}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
