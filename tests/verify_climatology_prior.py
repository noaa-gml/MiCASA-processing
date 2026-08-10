#!/usr/bin/env python3
"""Verification battery for the MiCASA CLIMATOLOGY PRIOR.

    python3 tests/verify_climatology_prior.py --outdir DIR --years 2021-2025 \
        --source /path/to/MiCASA_v1_flux_x360_y180_monthly.nc \
        --baseline 2001-2020 [--reference /path/to/reference/ERA5]

Unlike tests/verify_v2.py, this EXITS NON-ZERO when a check fails. The product
is a synthetic prior that looks exactly like a real one from the outside; a
verdict that is printed but not wired to an exit status is a comment, and this
is precisely the artefact where that would go unnoticed.

Everything is gated on artefact counts and numbers, never on a stage's exit
code -- a diurnalization that silently produced eleven months, or that quietly
reused the production fit, exits 0.

Sections
    1  the climatological monthly series      (is the input what we think?)
    2  the sub-monthly fit                    (delegated to verify_climatology_fit.r)
    3  completeness of the delivered product
    4  format parity with the reference product
    5  science: round-trip, trend removal, liveness

--selftest builds synthetic known-good and known-bad files and asserts the bad
ones FAIL, so the battery has demonstrated it can fail before any verdict from
it is believed.
"""
import argparse
import calendar
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
from budget import (grid_area, clean, pgc, nee_from_npp_rh, mol_to_gc,
                    seconds_in_month, days_in_year)

try:
    import netCDF4
except ImportError:
    netCDF4 = None

PASS, FAIL, WARN, INFO = "PASS", "FAIL", "WARN", "INFO"
_RESULTS = []

# Attributes whose whole purpose is to differ between the climatology prior and
# the ordinary product. Anything else differing is a finding.
EXPECTED_ATTR_DIFF = {
    "title", "summary", "history", "date_created", "processing_host",
    "input_monthly_flux", "input_monthly_flux_sha256",
    "input_pchip_fit", "input_pchip_fit_sha256",
    "input_npp_climatology", "input_rh_climatology",
    "input_monthly_flux_series", "input_monthly_flux_series_sha256",
    "flux_from_climatology", "prior_variant", "climatology_baseline",
    "climatology_baseline_requested", "climatology_method",
    "climatology_years_per_month", "climatology_span",
    "scope_note", "provenance_note",
    "meteo_source_by_day", "processing_pipeline_version",
    "processing_pipeline_commit", "processing_step", "source",
    # Tool-version stamp written by ncks during the day-split.
    "NCO",
    # Documentation attributes added to the pipeline after the reference
    # product was built. Their presence in the new files and absence in an
    # older reference is a pipeline-version difference, not a data difference;
    # both products follow the same convention (NEE = Rh - NPP, ATMC not
    # subtracted).
    "nee_atmc_convention", "atmc_convention",
}
STRUCTURAL_VAR_ATTRS = ("units", "_FillValue", "long_name")

# Attributes the delivered files MUST carry, so a downstream reader can tell
# what this product is from the file alone.
REQUIRED_CLIM_ATTRS = ("flux_from_climatology", "prior_variant",
                       "climatology_baseline", "climatology_method",
                       "scope_note")


def record(check_id, name, status, detail=""):
    color = {PASS: "\x1b[32m", FAIL: "\x1b[31m",
             WARN: "\x1b[33m", INFO: "\x1b[34m"}.get(status, "")
    reset = "\x1b[0m" if color else ""
    _RESULTS.append((check_id, name, status, detail))
    print(f"  {color}{status:4s}{reset}  {check_id:<7s} {name}"
          f"{('  -- ' + detail) if detail else ''}")


def section(title):
    print(f"\n=== {title} ===")


# ------------------------------------------------------------------ helpers

def series_path(outdir, version, y=None, m=None):
    p = os.path.join(outdir, "monthly_1x1",
                     f"MiCASA_{version}_flux_x360_y180_monthly")
    return f"{p}_{y}{m:02d}.nc" if y else f"{p}.nc"


def daily_path(era5, version, y, m, d):
    return os.path.join(era5, f"MiCASA_{version}.nee.{y}{m:02d}{d:02d}.nc")


def read_bio(path, idx=0):
    ds = netCDF4.Dataset(path)
    npp = ds.variables["NPP"]
    rh = ds.variables["Rh"]
    npp = npp[idx] if npp.ndim == 3 else npp[:]
    rh = rh[idx] if rh.ndim == 3 else rh[:]
    ds.close()
    return nee_from_npp_rh(npp, rh)


def describe(path):
    ds = netCDF4.Dataset(path)
    d = {"dims": {k: (len(v), v.isunlimited()) for k, v in ds.dimensions.items()},
         "vars": {}, "gattrs": {k: ds.getncattr(k) for k in ds.ncattrs()},
         "coords": {}}
    for name, v in ds.variables.items():
        filt = v.filters() or {}
        try:
            chunk = v.chunking()
        except Exception:
            chunk = None
        d["vars"][name] = {
            "dtype": str(v.dtype), "dims": tuple(v.dimensions),
            "shape": tuple(v.shape),
            "attrs": {a: v.getncattr(a) for a in v.ncattrs()},
            "deflate": filt.get("complevel", 0),
            "shuffle": bool(filt.get("shuffle", False)),
            "chunking": chunk,
        }
        if name in ("latitude", "longitude", "time"):
            d["coords"][name] = np.asarray(v[:], dtype="float64")
    ds.close()
    return d


def structural_diff(new, ref, compare_time=True):
    """Structural differences that would make the files non-interchangeable."""
    out = []
    if new["dims"] != ref["dims"]:
        out.append(f"dimensions {new['dims']} vs {ref['dims']}")
    nv, rv = set(new["vars"]), set(ref["vars"])
    if nv != rv:
        out.append(f"variables only-new={sorted(nv - rv)} only-ref={sorted(rv - nv)}")
    for name in sorted(nv & rv):
        a, b = new["vars"][name], ref["vars"][name]
        for key in ("dtype", "dims", "shape", "deflate", "shuffle", "chunking"):
            if a[key] != b[key]:
                out.append(f"{name}.{key}: {a[key]} vs {b[key]}")
        for at in STRUCTURAL_VAR_ATTRS:
            av, bv = a["attrs"].get(at), b["attrs"].get(at)
            if (av is None) != (bv is None):
                out.append(f"{name}.{at} presence differs")
            elif av is not None:
                same = (bool(np.isclose(av, bv))
                        if isinstance(av, (float, np.floating)) else av == bv)
                if not same:
                    out.append(f"{name}.{at}: {av!r} vs {bv!r}")
    # Shape-safe: a truncated axis must be REPORTED, not raise. A gate that
    # crashes on malformed input reads as a broken tool and gets rerun without
    # it, which is worse than no gate at all.
    def coord_diff(c, label):
        a, b = new["coords"].get(c), ref["coords"].get(c)
        if a is None or b is None:
            if (a is None) != (b is None):
                out.append(f"coordinate {c} present in only one file")
            return
        if a.shape != b.shape:
            out.append(f"coordinate {c} length {a.shape[0]} vs {b.shape[0]}")
        elif not np.allclose(a, b):
            out.append(label)

    coord_diff("latitude", "coordinate latitude values differ")
    coord_diff("longitude", "coordinate longitude values differ")
    if compare_time:
        coord_diff("time", "time values differ for the same calendar date")
    for k in sorted(set(new["gattrs"]) | set(ref["gattrs"])):
        if k in EXPECTED_ATTR_DIFF:
            continue
        if new["gattrs"].get(k) != ref["gattrs"].get(k):
            out.append(f"global attr {k!r}: "
                       f"{new['gattrs'].get(k)!r} vs {ref['gattrs'].get(k)!r}")
    return out


def day_mean_nee_gc(path):
    """Time-mean NEE (gC m-2 s-1) over a daily file's slots, and slot count."""
    ds = netCDF4.Dataset(path)
    v = ds.variables["NEE"]
    n = v.shape[0]
    a = mol_to_gc(v[:])
    ds.close()
    return a.mean(axis=0), n


# ------------------------------------------------------------------ sections

def sec1_series(args):
    section("1. Climatological monthly series")
    cat = series_path(args.outdir, args.version)
    if not os.path.exists(cat):
        record("1.1", "series concatenation present", FAIL, f"missing {cat}")
        return None
    ds = netCDF4.Dataset(cat)
    tv = ds.variables["time"]
    dates = netCDF4.num2date(tv[:], tv.units)
    yrs = np.array([d.year for d in dates])
    mons = np.array([d.month for d in dates])
    gattrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
    ds.close()
    record("1.1", "series concatenation present", PASS,
           f"{len(yrs)} months, {yrs[0]}-{mons[0]:02d}..{yrs[-1]}-{mons[-1]:02d}")

    # The series must cover every delivered year plus a pad month either side,
    # or the fit's edge slopes come from nothing.
    need = set(range(args.y0, args.y1 + 1))
    have = set(yrs.tolist())
    record("1.2", "series covers every delivered year",
           PASS if need <= have else FAIL,
           f"delivered {args.y0}-{args.y1}; series {min(have)}-{max(have)}")

    padded = (min(have) < args.y0) and (max(have) > args.y1)
    record("1.3", "series padded either side of the delivery window",
           PASS if padded else WARN,
           "edge months have climatological neighbours" if padded
           else "unpadded: first/last delivered month inherits an edge slope")

    # Genuinely climatological: same calendar month identical across years.
    worst, worst_m = 0.0, None
    for m in range(1, 13):
        idx = np.where(mons == m)[0]
        if len(idx) < 2:
            continue
        ref = read_bio(cat, idx[0])
        for i in idx[1:]:
            d = float(np.abs(read_bio(cat, i) - ref).max())
            if d > worst:
                worst, worst_m = d, m
    record("1.4", "series is genuinely climatological",
           PASS if worst == 0.0 else FAIL,
           f"max|same-calendar-month difference| = {worst:.3e} gC m-2 s-1"
           + (f" (month {worst_m})" if worst_m else ""))

    base_actual = gattrs.get("climatology_baseline", "?")
    base_req = gattrs.get("climatology_baseline_requested", base_actual)
    counts = gattrs.get("climatology_years_per_month", "")
    even = len(set(counts.split(","))) == 1 if counts else False
    record("1.5", "baseline window recorded and even",
           PASS if even else WARN,
           f"actual {base_actual}, requested {base_req}, years/month {counts}")
    if base_actual != base_req:
        record("1.6", "requested baseline was available", WARN,
               f"requested {base_req} but the record only supports {base_actual}")
    else:
        record("1.6", "requested baseline was available", PASS, base_actual)

    # 1.7 -- mixed provenance. If imposed channels are carried, the file claims
    # two different things about two sets of variables, and BOTH claims have to
    # hold: the climatological ones identical in every year, the carried ones
    # not. A carried variable that turned out to be identical across years
    # would mean it had been climatologized by accident; a climatological one
    # that varied would mean the opposite. Read the file's own labels rather
    # than a hardcoded list, so the check follows the product.
    pm = series_path(args.outdir, args.version, args.y0, 7)
    if os.path.exists(pm):
        ds = netCDF4.Dataset(pm)
        labels = {n: (v.getncattr("climatology")
                      if "climatology" in v.ncattrs() else None)
                  for n, v in ds.variables.items()
                  if v.ndim == 3}
        ds.close()
        clim_vars = [n for n, l in labels.items() if l == "yes"]
        imp_vars = [n for n, l in labels.items() if l == "no"]

        def july_spread(name):
            ref = None
            worst = 0.0
            for y in range(args.y0, args.y1 + 1):
                p = series_path(args.outdir, args.version, y, 7)
                if not os.path.exists(p):
                    continue
                d = netCDF4.Dataset(p)
                if name not in d.variables:
                    d.close()
                    continue
                a = clean(d.variables[name][0])
                d.close()
                if ref is None:
                    ref = a
                else:
                    worst = max(worst, float(np.abs(a - ref).max()))
            return worst

        if not labels or all(l is None for l in labels.values()):
            record("1.7", "variables labelled climatological / imposed", WARN,
                   "no per-variable `climatology` attribute to check")
        else:
            bad = []
            for n in clim_vars:
                if july_spread(n) != 0.0:
                    bad.append(f"{n} labelled climatological but VARIES by year")
            for n in imp_vars:
                if july_spread(n) == 0.0:
                    bad.append(f"{n} labelled imposed but is IDENTICAL every year")
            record("1.7", "per-variable climatology labels are truthful",
                   PASS if not bad else FAIL,
                   (f"climatological {sorted(clim_vars)} identical across years; "
                    f"imposed {sorted(imp_vars)} year-specific")
                   if not bad else "; ".join(bad))
    return gattrs


def sec2_fit(args):
    section("2. Sub-monthly coefficient fit")
    fit = os.path.join(args.outdir, "fit.piqs.rda")
    if not os.path.exists(fit):
        record("2.0", "fit present", FAIL, f"missing {fit}")
        return
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "verify_climatology_fit.r")
    if not os.path.exists(helper):
        record("2.0", "fit helper present", FAIL, f"missing {helper}")
        return
    try:
        r = subprocess.run(["Rscript", helper, fit,
                            series_path(args.outdir, args.version),
                            str(args.y0), str(args.y1)],
                           capture_output=True, text=True, timeout=3600)
    except (OSError, subprocess.SubprocessError) as e:
        record("2.0", "fit checks ran", FAIL, f"{e}")
        return
    ids = {"fit.window": ("2.1", "every delivered month has its own coefficients"),
           "fit.mean_preserving": ("2.2", "fit is mean-preserving"),
           "fit.flat": ("2.3", "fit carries no interannual structure")}
    seen = set()
    for line in r.stdout.splitlines():
        if line.startswith("CHECK|"):
            _, cid, status, detail = line.split("|", 3)
            num, name = ids.get(cid, (cid, cid))
            record(num, name, PASS if status == "PASS" else FAIL, detail)
            seen.add(cid)
        elif line.startswith("INFO|"):
            print(f"  {INFO}  {line[5:]}")
    missing = set(ids) - seen
    if missing:
        record("2.9", "all fit checks reported", FAIL,
               f"missing {sorted(missing)}; stderr: {r.stderr[-300:]}")


def sec3_completeness(args):
    section("3. Completeness of the delivered product")
    era5 = os.path.join(args.outdir, "ERA5")
    total_expected = total_found = 0
    zero_byte, bad_slots, missing_days = [], [], []
    for y in range(args.y0, args.y1 + 1):
        exp = days_in_year(y)
        found = 0
        for m in range(1, 13):
            for d in range(1, calendar.monthrange(y, m)[1] + 1):
                p = daily_path(era5, args.version, y, m, d)
                if not os.path.exists(p):
                    missing_days.append(f"{y}{m:02d}{d:02d}")
                    continue
                found += 1
                if os.path.getsize(p) == 0:
                    zero_byte.append(p)
        total_expected += exp
        total_found += found
        record(f"3.{y - args.y0 + 1}", f"daily NEE files for {y}",
               PASS if found == exp else FAIL, f"{found}/{exp}")
    record("3.a", "total daily file count",
           PASS if total_found == total_expected else FAIL,
           f"{total_found}/{total_expected} files "
           f"({args.y0}-01-01 .. {args.y1}-12-31)")
    record("3.b", "no zero-byte files", PASS if not zero_byte else FAIL,
           f"{len(zero_byte)} zero-byte" if zero_byte else "none")
    if missing_days:
        record("3.c", "no missing days", FAIL,
               f"{len(missing_days)} missing, first: {missing_days[:5]}")
    else:
        record("3.c", "no missing days", PASS, "every calendar day present")

    # Slot count + time convention, sampled across the record (opening 1826
    # files is wasteful; a stratified sample catches a systematic fault).
    import datetime as dt
    sample = []
    for y in range(args.y0, args.y1 + 1):
        for m in (1, 2, 6, 7, 12):
            for d in (1, 15, calendar.monthrange(y, m)[1]):
                p = daily_path(era5, args.version, y, m, d)
                if os.path.exists(p):
                    sample.append((y, m, d, p))
    worst_t = 0.0
    for y, m, d, p in sample:
        ds = netCDF4.Dataset(p)
        n = ds.variables["NEE"].shape[0]
        t = np.asarray(ds.variables["time"][:], dtype="float64")
        ds.close()
        if n != 24:
            bad_slots.append(f"{p}:{n}")
            continue
        base = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
        day0 = dt.datetime(y, m, d, tzinfo=dt.timezone.utc)
        want = np.array([((day0 - base).total_seconds() + h * 3600 + 1800) / 86400.0
                         for h in range(24)])
        worst_t = max(worst_t, float(np.abs(t - want).max()))
    record("3.d", "24 hourly slots per file",
           PASS if not bad_slots else FAIL,
           f"{len(sample)} files sampled" if not bad_slots
           else f"{len(bad_slots)} bad: {bad_slots[:3]}")
    record("3.e", "time axis is 00:30..23:30 UTC of the nominal date",
           PASS if worst_t < 1e-9 else FAIL,
           f"max deviation {worst_t:.3e} days over {len(sample)} files")
    # Leap day must exist and must not be a copy of Feb 28.
    for y in range(args.y0, args.y1 + 1):
        if calendar.isleap(y):
            p = daily_path(era5, args.version, y, 2, 29)
            ok = os.path.exists(p)
            detail = "present" if ok else "MISSING"
            if ok:
                a, _ = day_mean_nee_gc(p)
                b, _ = day_mean_nee_gc(daily_path(era5, args.version, y, 2, 28))
                differs = float(np.abs(a - b).max()) > 0
                ok = differs
                detail = (f"present and distinct from Feb 28 "
                          f"(max|d| = {float(np.abs(a - b).max()):.3e} gC m-2 s-1)"
                          if differs else "present but IDENTICAL to Feb 28")
            record("3.f", f"leap day {y}-02-29", PASS if ok else FAIL, detail)


def sec4_parity(args):
    section("4. Format parity with the reference product")
    if not args.reference:
        record("4.0", "reference product supplied", WARN,
               "no --reference; parity not checked")
        return
    era5 = os.path.join(args.outdir, "ERA5")
    # A date that exists in both, mid-record.
    probe = None
    for y in range(args.y0, args.y1 + 1):
        cand = daily_path(era5, args.version, y, 7, 15)
        ref = daily_path(args.reference, args.version, y, 7, 15)
        if os.path.exists(cand) and os.path.exists(ref):
            probe = (cand, ref, f"{y}0715")
            break
    if probe is None:
        record("4.0", "counterpart file in the reference tree", FAIL,
               "no shared date found")
        return
    new_p, ref_p, tag = probe
    record("4.0", "counterpart file found", INFO, tag)
    diffs = structural_diff(describe(new_p), describe(ref_p))
    record("4.1", "structurally interchangeable with the reference",
           PASS if not diffs else FAIL,
           "dims, variables, dtypes, units, fill, chunking, compression, "
           "coordinates and time all identical" if not diffs
           else f"{len(diffs)} differences: {diffs[:4]}")

    ng = describe(new_p)["gattrs"]
    missing = [a for a in REQUIRED_CLIM_ATTRS if a not in ng]
    record("4.2", "delivered files are self-describing as a climatology",
           PASS if not missing else FAIL,
           f"carries {', '.join(REQUIRED_CLIM_ATTRS)}" if not missing
           else f"missing {missing}")
    if "flux_from_climatology" in ng:
        record("4.3", "flux_from_climatology says yes",
               PASS if ng["flux_from_climatology"] == "yes" else FAIL,
               f"{ng['flux_from_climatology']!r}")


def walk_tree(dirpath, version, y0, y1, area):
    """One pass over a daily-NEE tree.

    Returns (monthly, annual) where monthly[(y, m)] is the month's mean NEE
    field in gC m-2 s-1 and annual is [(year, PgC/yr, nfiles)]. Both the
    round-trip and the trend table are derived from this, so the two arms are
    never measured by two different code paths -- and each tree is read once.
    """
    monthly, annual = {}, []
    for y in range(y0, y1 + 1):
        tot, nfiles = 0.0, 0
        for m in range(1, 13):
            acc, slots = None, 0
            for d in range(1, calendar.monthrange(y, m)[1] + 1):
                p = daily_path(dirpath, version, y, m, d)
                if not os.path.exists(p):
                    continue
                dm, n = day_mean_nee_gc(p)
                acc = dm * n if acc is None else acc + dm * n
                slots += n
                tot += pgc(dm, area, n * 3600.0)
                nfiles += 1
            if acc is not None:
                monthly[(y, m)] = acc / slots
        annual.append((y, tot, nfiles))
    return monthly, annual


def roundtrip_worst(monthly, source_for):
    """Worst relative deviation between produced month-means and their input.

    monthly    {(year, month): produced mean field}
    source_for callable (year, month) -> path of the monthly file it was built
               from, or None if there isn't one.
    """
    worst, tag, n = 0.0, "", 0
    for (y, m), produced in sorted(monthly.items()):
        src = source_for(y, m)
        if not src or not os.path.exists(src):
            continue
        ref = read_bio(src, 0)
        denom = float(np.abs(ref).max()) or 1.0
        rel = float(np.abs(produced - ref).max()) / denom
        n += 1
        if rel > worst:
            worst, tag = rel, f"{y}-{m:02d}"
    return worst, tag, n


def sec5_science(args):
    section("5. Science: round-trip, trend removal, liveness")
    era5 = os.path.join(args.outdir, "ERA5")
    area = grid_area()

    new_monthly, new_rows = walk_tree(era5, args.version, args.y0, args.y1, area)

    # --- 5.1 round-trip: the produced dailies must average back to the
    # climatological monthly mean they were built from.
    #
    # The residual is NOT zero and is not supposed to be: diurnalize samples
    # the sub-monthly quadratic at mid-hour points and the conserve polar clip
    # renormalizes GPP onto the monthly-mean array, so a few parts in 10^3 of
    # the peak cell is the pipeline's own discretization. The meaningful
    # question is therefore comparative -- is it WORSE than what the
    # unmodified product already does? -- so when a reference is supplied the
    # same measurement is made on it and the ratio is the verdict.
    worst, tag, nmon = roundtrip_worst(
        new_monthly, lambda y, m: series_path(args.outdir, args.version, y, m))
    record("5.1", "round-trip: monthly means reproduce the climatology",
           PASS if worst < args.roundtrip_tol else FAIL,
           f"worst relative deviation {worst:.3e} ({tag}) over {nmon} months, "
           f"ceiling {args.roundtrip_tol:.0e}")

    ref_monthly = ref_rows = None
    if args.reference:
        ref_monthly, ref_rows = walk_tree(args.reference, args.version,
                                          args.y0, args.y1, area)
        ref_mon_dir = args.reference_monthly or os.path.join(
            os.path.dirname(os.path.normpath(args.reference)), "monthly_1x1")

        def ref_src(y, m):
            return os.path.join(
                ref_mon_dir,
                f"MiCASA_{args.version}_flux_x360_y180_monthly_{y}{m:02d}.nc")

        rworst, rtag, rn = roundtrip_worst(ref_monthly, ref_src)
        if rn:
            ratio = worst / rworst if rworst else float("inf")
            record("5.1b", "round-trip residual is inherited, not introduced",
                   PASS if ratio <= 2.0 else FAIL,
                   f"climatology {worst:.3e} vs unmodified product {rworst:.3e} "
                   f"over {rn} months -> {ratio:.2f}x (the residual is the "
                   f"diurnalization's own discretization)")
        else:
            record("5.1b", "round-trip control on the reference product", WARN,
                   f"no reference monthly files under {ref_mon_dir}")

    # --- 5.2 trend removal, measured with the same instrument on both arms.
    print("\n  Global annual land NET BIO flux (NEE = Rh - NPP), PgC/yr, "
          "positive = source to atmosphere")
    hdr = f"  {'year':>6} {'climatology':>13}"
    if ref_rows:
        hdr += f" {'raw MiCASA':>13} {'difference':>12}"
    print(hdr)
    for i, (y, v, n) in enumerate(new_rows):
        line = f"  {y:>6} {v:>13.4f}"
        if ref_rows:
            rv = ref_rows[i][1]
            line += f" {rv:>13.4f} {v - rv:>12.4f}"
        print(line)

    ys = np.array([r[0] for r in new_rows], dtype=float)
    vs = np.array([r[1] for r in new_rows], dtype=float)
    slope = float(np.polyfit(ys, vs, 1)[0]) if len(ys) > 1 else 0.0
    spread = float(vs.max() - vs.min())
    print(f"  {'mean':>6} {vs.mean():>13.4f}")
    print(f"  {'trend':>6} {slope:>13.4f} PgC/yr per yr   (spread {spread:.4f})")
    record("5.2", "trend removed from the delivered prior",
           PASS if abs(slope) < args.trend_tol else FAIL,
           f"least-squares trend {slope:+.4f} PgC/yr per yr, "
           f"year-to-year spread {spread:.4f} PgC/yr, tol {args.trend_tol}")

    if ref_rows:
        rv = np.array([r[1] for r in ref_rows], dtype=float)
        rslope = float(np.polyfit(ys, rv, 1)[0])
        record("5.3", "reference product does show the trend (control)",
               PASS if abs(rslope) > abs(slope) * 10 else WARN,
               f"raw trend {rslope:+.4f} vs climatology {slope:+.4f} PgC/yr per yr "
               f"({abs(rslope) / max(abs(slope), 1e-12):.0f}x)")

        # --- liveness: the two arms must actually differ.
        probe = None
        for y in range(args.y0, args.y1 + 1):
            a = daily_path(era5, args.version, y, 7, 15)
            b = daily_path(args.reference, args.version, y, 7, 15)
            if os.path.exists(a) and os.path.exists(b):
                probe = (a, b, f"{y}-07-15")
                break
        if probe:
            a, b, tag = probe
            da = netCDF4.Dataset(a); db = netCDF4.Dataset(b)
            x = clean(da.variables["NEE"][:]); yv = clean(db.variables["NEE"][:])
            da.close(); db.close()
            d = np.abs(x - yv)
            ndiff = int((d > 0).sum())
            record("5.4", "treatment reached the product (liveness)",
                   PASS if ndiff > 0 else FAIL,
                   f"{ndiff:,}/{d.size:,} cell-hours differ on {tag}, "
                   f"max|d| = {float(d.max()):.3e} mol m-2 s-1")

    # --- 5.5 the delivered values must be finite and physically plausible.
    bad = 0
    mx = 0.0
    for y in range(args.y0, args.y1 + 1):
        p = daily_path(era5, args.version, y, 7, 15)
        if not os.path.exists(p):
            continue
        ds = netCDF4.Dataset(p)
        v = ds.variables["NEE"][:]
        ds.close()
        raw = (v.filled(np.nan) if isinstance(v, np.ma.MaskedArray)
               else np.asarray(v))
        bad += int((~np.isfinite(raw)).sum() + (np.abs(raw) > 1e30).sum())
        mx = max(mx, float(np.abs(clean(v)).max()))
    record("5.5", "no non-finite or sentinel values in delivered NEE",
           PASS if bad == 0 else FAIL, f"{bad} bad cells")
    record("5.6", "peak |NEE| is physically plausible",
           PASS if 1e-8 < mx < 1e-3 else WARN,
           f"max |NEE| = {mx:.3e} mol m-2 s-1")


# ------------------------------------------------------------------ selftest

def _mkdaily(path, *, units="mol m-2 s-1", fill=-1e34, dtype="f4",
             nt=24, deflate=4, values=None, extra_var=False, lat_off=0.0,
             nlat=18, nlon=36, gattrs=None):
    ds = netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC")
    ds.createDimension("time", None)
    ds.createDimension("latitude", nlat)
    ds.createDimension("longitude", nlon)
    v = ds.createVariable("NEE", dtype, ("time", "latitude", "longitude"),
                          fill_value=np.dtype(dtype).type(fill),
                          zlib=deflate > 0, complevel=max(deflate, 1),
                          shuffle=True, chunksizes=(1, nlat, nlon))
    v.units = units
    v.long_name = "NEE=GPP+RESP, positive is source to atm, as is each component"
    v[:] = np.zeros((nt, nlat, nlon)) if values is None else values
    la = ds.createVariable("latitude", "f8", ("latitude",))
    la.units, la.long_name = "degrees_north", "latitude"
    la[:] = np.arange(nlat, dtype="float64") + lat_off
    lo = ds.createVariable("longitude", "f8", ("longitude",))
    lo.units, lo.long_name = "degrees_east", "longitude"
    lo[:] = np.arange(nlon, dtype="float64")
    tv = ds.createVariable("time", "f8", ("time",))
    tv.units, tv.long_name = "days since 1970-01-01 00:00:00 UTC", "time"
    tv[:] = np.arange(nt) / 24.0
    if extra_var:
        e = ds.createVariable("BONUS", "f4", ("time", "latitude", "longitude"))
        e.units, e.long_name = "1", "spurious"
        e[:] = 0.0
    ds.Conventions = "CF-1.10, ACDD-1.3"
    ds.institution = "NOAA Global Monitoring Laboratory"
    for k, val in (gattrs or {}).items():
        setattr(ds, k, val)
    ds.close()


def selftest():
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}"
              f"{('  -- ' + detail) if detail else ''}")

    tmp = tempfile.mkdtemp()
    ref = os.path.join(tmp, "ref.nc")
    _mkdaily(ref)

    good = os.path.join(tmp, "good.nc")
    _mkdaily(good)
    check("identical files are structurally interchangeable",
          not structural_diff(describe(good), describe(ref)))

    prov = os.path.join(tmp, "prov.nc")
    _mkdaily(prov, gattrs={"title": "different", "flux_from_climatology": "yes"})
    check("provenance-only differences do not fail parity",
          not structural_diff(describe(prov), describe(ref)))

    for name, kw in [("wrong units", dict(units="gC m-2 s-1")),
                     ("wrong fill value", dict(fill=-9999.0)),
                     ("wrong dtype", dict(dtype="f8")),
                     ("wrong deflate level", dict(deflate=1)),
                     ("extra variable", dict(extra_var=True)),
                     ("shifted latitudes", dict(lat_off=0.5)),
                     ("truncated time axis", dict(nt=12))]:
        bad = os.path.join(tmp, f"bad{abs(hash(name))}.nc")
        _mkdaily(bad, **kw)
        d = structural_diff(describe(bad), describe(ref))
        check(f"negative control: {name} FAILS parity", bool(d),
              d[0] if d else "NOT DETECTED")

    # Budget arithmetic reaches the right answer through this module's helpers.
    area = grid_area(18, 36)
    check("uniform flux integrates to the analytic value",
          abs(pgc(np.ones((18, 36)), area, seconds_in_month(2021, 1))
              - area.sum() * 31 * 86400.0 / 1e15) < 1e-6)
    check("negative control: unmasked sentinel is absurd",
          abs(pgc(np.full((18, 36), -1e34), area,
                  seconds_in_month(2021, 1))) > 1e6)

    # Round-trip arithmetic: a constant NEE file must average back exactly.
    vals = np.full((24, 18, 36), 2.0 / 12.011)   # -> 2.0 gC m-2 s-1
    live = os.path.join(tmp, "live.nc")
    _mkdaily(live, values=vals)
    dm, n = day_mean_nee_gc(live)
    check("day_mean_nee_gc converts mol -> gC and averages",
          n == 24 and abs(float(dm.mean()) - 2.0) < 1e-6,
          f"mean {float(dm.mean()):.6f} gC m-2 s-1")

    print(f"\nSELFTEST {'PASSED' if ok else 'FAILED'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir")
    ap.add_argument("--reference", help="ERA5 dir of the unmodified product")
    ap.add_argument("--source", help="concatenated monthly file that was averaged")
    ap.add_argument("--baseline", default="")
    ap.add_argument("--years", default="2021-2025")
    ap.add_argument("--version", default=os.environ.get("MICASA_VERSION", "v1"))
    ap.add_argument("--reference-monthly",
                    help="reference tree's monthly_1x1 dir (default: sibling "
                         "of --reference)")
    ap.add_argument("--roundtrip-tol", type=float, default=5e-3,
                    help="absolute ceiling on the round-trip residual. The "
                         "meaningful test is check 5.1b, which compares it "
                         "against the unmodified product's own residual "
                         "(~1.2e-3, the diurnalization's discretization).")
    ap.add_argument("--trend-tol", type=float, default=0.02,
                    help="PgC/yr per yr; the raw product's is ~0.4")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(selftest())
    if not args.outdir:
        ap.error("--outdir required (or --selftest)")
    if netCDF4 is None:
        sys.exit("netCDF4 unavailable")
    y0, _, y1 = args.years.partition("-")
    args.y0, args.y1 = int(y0), int(y1 or y0)

    print("=" * 72)
    print("MiCASA climatology prior -- verification")
    print(f"  outdir     {args.outdir}")
    print(f"  years      {args.y0}-{args.y1}")
    print(f"  baseline   {args.baseline or '(from file attributes)'}")
    print(f"  reference  {args.reference or '(none)'}")
    print("=" * 72)

    sec1_series(args)
    sec2_fit(args)
    sec3_completeness(args)
    sec4_parity(args)
    sec5_science(args)

    n_pass = sum(1 for r in _RESULTS if r[2] == PASS)
    n_fail = sum(1 for r in _RESULTS if r[2] == FAIL)
    n_warn = sum(1 for r in _RESULTS if r[2] == WARN)
    print("\n" + "=" * 72)
    print(f"SUMMARY  {n_pass} PASS  {n_warn} WARN  {n_fail} FAIL")
    if n_fail:
        print("\nFAILURES:")
        for cid, name, status, detail in _RESULTS:
            if status == FAIL:
                print(f"  {cid:<7s} {name} -- {detail}")
    print("=" * 72)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
