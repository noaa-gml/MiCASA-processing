#!/usr/bin/env python3
"""Build a CLIMATOLOGICAL MiCASA monthly series -- stage 1 of the climatology prior.

Requested by Andy Jacobson 2026-08-10 (CT2026 Issues): "use the 2000-2020
climatology as the prior. This means computing the mean (climatology) from
monthly data, then running it through the diurnalization process. I.e.
climatological fluxes but diurnalized using real meteorology."

The product is a prior whose MAGNITUDE carries no interannual variability or
trend, but whose day-to-day and hour-to-hour structure is entirely real, because
the stock diurnalization is driven by each target day's own ERA5 meteorology.

WHAT THIS SCRIPT EMITS, AND WHY IT IS A SERIES RATHER THAN 12 FIELDS
--------------------------------------------------------------------
Two artefacts, both required, in $MONTHLY_1X1_DIR:

    MiCASA_<ver>_flux_x360_y180_monthly_YYYYMM.nc   one per month of the span
    MiCASA_<ver>_flux_x360_y180_monthly.nc          the concatenation

diurnalize-ERA5.r reads the per-month file for the diurnal AMPLITUDE; the
fitter (write_pchip.r et al.) reads the concatenation to build fit.piqs.rda.

Both matter, and the second is the subtle one. diurnalize's hourly flux is

    f(t) = driver(t) * mean / mean_driver  -  mean  +  qmod(t)

whose monthly mean is mean(qmod) -- taken ENTIRELY from the coefficient fit,
not from the monthly-mean array. So a climatology prior that reused the
production fit.piqs.rda would silently reinstate the real interannual signal
through qmod while producing a complete, healthy-looking set of files. The fit
must be rebuilt on the climatological series, which is why the series exists.
run_climatology_prior.sh does that; see docs/CLIMATOLOGY_PRIOR.md.

The span is padded a year either side of the delivery window by default so the
PCHIP slopes at the first and last delivered month have climatological
neighbours rather than inheriting one from real data.

SCOPE: BIO ONLY. NPP and Rh are written; FIRE, FUEL and ATMC are deliberately
not carried. CarbonTracker reads fire and fuel from a separate rc key
(ct.wildfire -> fires.input.dir) pointing at the daily_1x1 tree, which this
product does not replace -- so a climatological fire channel would be both
unused and misleading. Their absence makes that scoping decision visible in the
files rather than resting on a convention.

Configuration (config.sh; all overridable on the command line):
    MICASA_CLIM_SOURCE          concatenated monthly file to average
    MICASA_CLIM_BASELINE_START  first baseline year   (default 2001)
    MICASA_CLIM_BASELINE_END    last  baseline year   (default 2020)
    MICASA_CLIM_SPAN_START      first month of output (default 2020-01)
    MICASA_CLIM_SPAN_END        last  month of output (default 2026-12)
    MONTHLY_1X1_DIR             where the series is written
    MICASA_VERSION              v1 | vNRT (filename component)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from provenance import provenance_attrs
from climatology import (modulo_month_mean, modulo_month_counts,
                         month_span, month_mid_epoch, parse_yyyymm)

FILL = -1.0e34
BIO_VARS = ("NPP", "Rh")

# Kept identical to what ingest_monthly.r stamps, so a reader cannot tell the
# difference in anything except the deliberately-added climatology attributes.
VAR_ATTRS = {
    "NPP": dict(
        units="gC m^-2 s^-1",
        long_name="Net primary productivity",
        sign_and_usage=("Stored POSITIVE, but this is an UPTAKE term: NEE is "
                        "formed as (Rh - NPP), so the stored sign is NOT this "
                        "term's sign in the atmospheric budget."),
    ),
    "Rh": dict(
        units="gC m^-2 s^-1",
        long_name="Heterotrophic respiration",
        sign_and_usage="Stored POSITIVE = source to the atmosphere.",
    ),
}

SCOPE_NOTE_BIO_ONLY = (
    "BIO CHANNEL ONLY (NEE = Rh - NPP). FIRE, FUEL and ATMC are deliberately "
    "absent: CarbonTracker reads fire and fuel from a separate rc key "
    "(ct.wildfire -> fires.input.dir) pointing at the daily_1x1 tree, which "
    "this product does not replace. Do not treat these files as a drop-in "
    "substitute for the 1-degree daily MiCASA product.")

SCOPE_NOTE_MIXED = (
    "MIXED PROVENANCE -- READ THE PER-VARIABLE `climatology` ATTRIBUTE. "
    "ONLY NPP and Rh are climatological. {imposed} are copied UNCHANGED from "
    "the source archive and are REAL, YEAR-SPECIFIC fields. This mirrors what "
    "CarbonTracker actually does: it takes bio from this product and imposes "
    "fire/fuel separately (ct.wildfire -> fires.input.dir) from the unmodified "
    "daily tree, so a file combining climatological bio with real fire is a "
    "faithful picture of the prior in use. It is NOT a fully climatological "
    "product, and nothing in it should be summed as though it were.")

PROV_NOTE = (
    "SYNTHETIC SERIES. Every instance of a given calendar month is identical by "
    "construction IN THE CLIMATOLOGICAL VARIABLES; interannual variability and "
    "trend are removed by design in those, and retained in any carried "
    "imposed variables. Not an observational product.")

CLIM_VAR_NOTE = ("CLIMATOLOGICAL: the mean of this calendar month over the "
                 "baseline years. Identical in every year of the series.")
IMPOSED_VAR_NOTE = ("NOT CLIMATOLOGICAL: copied unchanged from the source "
                    "archive for this specific month. Real, year-specific.")


def env(name, default=None):
    v = os.environ.get(name, "")
    return v if v else default


def read_imposed(args, year, month):
    """Imposed (non-climatological) variables for one month, copied unchanged.

    Returns {name: (values, attrs)} for whichever of args.carry_imposed the
    source per-month file actually has. A missing month or variable yields
    nothing -- the caller reports the gap rather than inventing a field.

    These are REAL, year-specific data. They are carried so the monthly product
    is usable by consumers who expect the full MiCASA variable set (Andy reads
    FIRE from it), and because CarbonTracker genuinely does combine
    climatological bio with real imposed fire. Every carried variable is
    stamped `climatology = "no"`.
    """
    if not args.carry_imposed:
        return {}
    import netCDF4
    src = os.path.join(args.imposed_source_dir,
                       "MiCASA_%s_flux_x360_y180_monthly_%d%02d.nc"
                       % (args.version, year, month))
    if not os.path.exists(src):
        return {}
    out = {}
    ds = netCDF4.Dataset(src)
    for name in args.carry_imposed:
        if name not in ds.variables:
            continue
        var = ds.variables[name]
        vals = var[0] if var.ndim == 3 else var[:]
        if np.ma.isMaskedArray(vals):
            vals = vals.filled(FILL)
        attrs = {a: var.getncattr(a) for a in var.ncattrs()
                 if a != "_FillValue"}
        attrs["climatology"] = "no"
        attrs["provenance"] = IMPOSED_VAR_NOTE
        attrs["copied_from"] = src
        out[name] = (np.asarray(vals, dtype="float64"), attrs)
    ds.close()
    return out


def build(args):
    import netCDF4

    src = args.source
    if not src or not os.path.exists(src):
        sys.exit("make_climatology_series: source monthly file not found: %r\n"
                 "                          set MICASA_CLIM_SOURCE or --source "
                 "(run cat_monthly.sh to build one)." % (src,))

    ds = netCDF4.Dataset(src)
    tv = ds.variables["time"]
    dates = netCDF4.num2date(tv[:], tv.units)
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])

    print("make_climatology_series: source %s" % src)
    print("  record            %d-%02d .. %d-%02d  (%d months)"
          % (years[0], months[0], years[-1], months[-1], len(years)))
    print("  baseline          %d-%d" % (args.baseline_start, args.baseline_end))

    counts = modulo_month_counts(months, years,
                                 args.baseline_start, args.baseline_end)
    present = sorted(set(years[(years >= args.baseline_start) &
                               (years <= args.baseline_end)].tolist()))
    if not present:
        sys.exit("make_climatology_series: no data in the baseline window "
                 "%d-%d" % (args.baseline_start, args.baseline_end))
    requested = set(range(args.baseline_start, args.baseline_end + 1))
    missing = sorted(requested - set(present))
    print("  baseline years    %d present (%d..%d)"
          % (len(present), present[0], present[-1]))
    if missing:
        print("  ** WARNING: requested baseline years absent from the record: %s"
              % (missing,))
        print("  ** the climatology is over %d-%d, NOT %d-%d as requested"
              % (present[0], present[-1], args.baseline_start, args.baseline_end))
    print("  months/calendar-month %s" % (counts.tolist(),))
    if len(set(counts.tolist())) != 1:
        print("  ** NOTE: uneven year counts per calendar month -- the "
              "climatology weights calendar months unequally.")

    fields = {}
    for name in BIO_VARS:
        if name not in ds.variables:
            sys.exit("make_climatology_series: variable %s missing from %s"
                     % (name, src))
        arr = ds.variables[name][:]
        if np.ma.isMaskedArray(arr):
            arr = arr.filled(np.nan)
        arr = np.asarray(arr, dtype="float64")
        # Defensive: the house -1e34 sentinel must never reach a mean.
        arr[np.abs(arr) > 1e30] = np.nan
        clim = modulo_month_mean(arr, months, years,
                                 args.baseline_start, args.baseline_end)
        # Cells with no data anywhere in the baseline stay 0, matching the
        # source product's convention (ocean is stored as exact 0.0, not fill).
        clim = np.where(np.isfinite(clim), clim, 0.0)
        fields[name] = clim
        print("  %-4s clim %s  global mean %.6e gC m-2 s-1"
              % (name, clim.shape, clim.mean()))
    ds.close()

    span = month_span(args.span_start, args.span_end)
    baseline_actual = "%d-%d" % (present[0], present[-1])
    extra = {
        "micasa_version": args.version,
        "prior_variant": "MiCASA climatology %s (bio only)" % baseline_actual,
        "flux_from_climatology": "yes",
        "climatology_baseline": baseline_actual,
        "climatology_baseline_requested": "%d-%d" % (args.baseline_start,
                                                     args.baseline_end),
        "climatology_method": (
            "modulo-month mean of the monthly-mean RATE (gC m-2 s-1) over the "
            "baseline years. February averages 28- and 29-day Februaries as "
            "rates, so a leap-year February carries the climatological rate "
            "over 29 days."),
        "climatology_years_per_month": ",".join(str(c) for c in counts),
        "climatology_span": "%s .. %s" % (args.span_start, args.span_end),
        "climatological_variables": ",".join(BIO_VARS),
        "scope_note": (SCOPE_NOTE_MIXED.format(
                           imposed=", ".join(args.carry_imposed))
                       if args.carry_imposed else SCOPE_NOTE_BIO_ONLY),
        "provenance_note": PROV_NOTE,
    }
    title = ("MiCASA %s climatological monthly land carbon flux (NPP, Rh%s), "
             "%s baseline"
             % (args.version,
                "; " + ",".join(args.carry_imposed) + " carried unchanged"
                if args.carry_imposed else "",
                baseline_actual))
    summary = (
        "Modulo-month climatology of MiCASA NPP and heterotrophic respiration: "
        "for each calendar month, the mean of that month across %s, replicated "
        "onto the target span. Built as the monthly input to the MiCASA "
        "diurnalization so the resulting hourly fluxes are climatological in "
        "magnitude but driven by each target day's real ERA5 meteorology. "
        "Sign convention: NEE = Rh - NPP, positive = source to atmosphere."
        % baseline_actual)

    gattrs = provenance_attrs(step="make_climatology_series.py",
                              work_dir=args.work_dir, title=title,
                              summary=summary,
                              inputs={"monthly_flux_series": src},
                              extra=extra)

    os.makedirs(args.outdir, exist_ok=True)
    prefix = "MiCASA_%s_flux_x360_y180_monthly" % args.version

    def coords(dset):
        dset.createDimension("longitude", 360)
        dset.createDimension("latitude", 180)
        dset.createDimension("time", None)
        lon = dset.createVariable("longitude", "f8", ("longitude",))
        lon.units, lon.long_name = "degrees_east", "longitude"
        lon[:] = np.arange(-179.5, 180.0, 1.0)
        lat = dset.createVariable("latitude", "f8", ("latitude",))
        lat.units, lat.long_name = "degrees_north", "latitude"
        lat[:] = np.arange(-89.5, 90.0, 1.0)
        t = dset.createVariable("time", "f8", ("time",))
        t.units = "seconds since 1970-01-01 00:00:00 UTC"
        t.long_name = "time"
        return t

    def bio_var(dset, name, **kw):
        v = dset.createVariable(name, "f4", ("time", "latitude", "longitude"),
                                fill_value=np.float32(FILL), **kw)
        for k, val in VAR_ATTRS[name].items():
            setattr(v, k, val)
        v.climatology = "yes"
        v.provenance = CLIM_VAR_NOTE
        return v

    def imposed_var(dset, name, attrs, **kw):
        v = dset.createVariable(name, "f4", ("time", "latitude", "longitude"),
                                fill_value=np.float32(FILL), **kw)
        for k, val in attrs.items():
            setattr(v, k, val)
        return v

    # Imposed channels, read once per month and reused for both writers.
    imposed = {}
    if args.carry_imposed:
        for y, m in span:
            got_month = read_imposed(args, y, m)
            if got_month:
                imposed[(y, m)] = got_month
        have = sorted({n for d_ in imposed.values() for n in d_})
        missing = [f"{y}-{m:02d}" for (y, m) in span if (y, m) not in imposed]
        print("  carried imposed (real, year-specific): %s over %d/%d months"
              % (",".join(have) if have else "NONE", len(imposed), len(span)))
        if missing:
            print("  ** %d month(s) have no source file, so carry NOTHING: %s%s"
                  % (len(missing), ", ".join(missing[:6]),
                     " ..." if len(missing) > 6 else ""))

    written = 0
    for y, m in span:
        path = os.path.join(args.outdir, "%s_%d%02d.nc" % (prefix, y, m))
        d = netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC")
        t = coords(d)
        t[0] = month_mid_epoch(y, m)
        for name in BIO_VARS:
            bio_var(d, name)[0, :, :] = fields[name][m - 1]
        for name, (vals, attrs) in sorted(imposed.get((y, m), {}).items()):
            imposed_var(d, name, attrs)[0, :, :] = vals
        for k, val in gattrs.items():
            setattr(d, k, val)
        setattr(d, "carried_imposed_variables",
                ",".join(sorted(imposed.get((y, m), {}))) or "none")
        d.close()
        written += 1
    print("  wrote %d per-month files to %s" % (written, args.outdir))

    cat = os.path.join(args.outdir, "%s.nc" % prefix)
    d = netCDF4.Dataset(cat, "w", format="NETCDF4_CLASSIC")
    t = coords(d)
    var = {}
    for name in BIO_VARS:
        v = bio_var(d, name, zlib=True, complevel=4)
        v.cell_methods = "time: mean"
        var[name] = v
    # Only carry a variable into the concatenation if EVERY month has it --
    # a partially populated field would read as real data with silent holes.
    cat_names = sorted(set.intersection(
        *[set(imposed[k]) for k in imposed]) ) if len(imposed) == len(span) else []
    for name in cat_names:
        attrs = imposed[span[0]][name][1]
        v = imposed_var(d, name, attrs, zlib=True, complevel=4)
        v.cell_methods = "time: mean"
        var[name] = v
    for i, (y, m) in enumerate(span):
        t[i] = month_mid_epoch(y, m)
        for name in BIO_VARS:
            var[name][i, :, :] = fields[name][m - 1]
        for name in cat_names:
            var[name][i, :, :] = imposed[(y, m)][name][0]
    for k, val in gattrs.items():
        setattr(d, k, val)
    setattr(d, "carried_imposed_variables", ",".join(cat_names) or "none")
    d.close()
    print("  wrote concatenation %s (%d months, imposed: %s)"
          % (cat, len(span), ",".join(cat_names) or "none"))

    # GATE: artefact count, not the exit code. A partially written span would
    # otherwise surface much later as a missing month inside diurnalize.
    got = len([f for f in os.listdir(args.outdir)
               if f.startswith(prefix + "_") and f.endswith(".nc")])
    if got != len(span):
        sys.exit("FATAL: expected %d per-month files in %s, found %d"
                 % (len(span), args.outdir, got))
    if not os.path.getsize(cat):
        sys.exit("FATAL: concatenation %s is empty" % cat)
    print("  GATE series-count: %d/%d per-month files, concatenation present"
          % (got, len(span)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=env("MICASA_CLIM_SOURCE"),
                    help="concatenated monthly file to average")
    ap.add_argument("--outdir", default=env("MONTHLY_1X1_DIR", "monthly_1x1"))
    ap.add_argument("--version", default=env("MICASA_VERSION", "v1"))
    ap.add_argument("--baseline-start", type=int,
                    default=int(env("MICASA_CLIM_BASELINE_START", "2001")))
    ap.add_argument("--baseline-end", type=int,
                    default=int(env("MICASA_CLIM_BASELINE_END", "2020")))
    ap.add_argument("--span-start", default=env("MICASA_CLIM_SPAN_START", "2020-01"))
    ap.add_argument("--span-end", default=env("MICASA_CLIM_SPAN_END", "2026-12"))
    ap.add_argument("--work-dir", default=env("WORK_DIR", os.getcwd()))
    ap.add_argument("--carry-imposed",
                    default=env("MICASA_CLIM_CARRY_IMPOSED", "FIRE,FUEL,ATMC"),
                    help="comma-separated variables copied UNCHANGED (real, "
                         "year-specific) from the source monthly product, so "
                         "the output carries the full MiCASA variable set. "
                         "Each is stamped climatology='no'. Pass an empty "
                         "string for a strictly bio-only product.")
    ap.add_argument("--imposed-source-dir",
                    default=env("MICASA_CLIM_IMPOSED_DIR"),
                    help="directory of per-month source files to copy the "
                         "imposed variables from (default: alongside --source)")
    args = ap.parse_args()

    args.carry_imposed = [s.strip() for s in args.carry_imposed.split(",")
                          if s.strip()]
    if args.imposed_source_dir is None:
        args.imposed_source_dir = os.path.dirname(os.path.abspath(args.source)) \
            if args.source else ""

    # normalise / validate the span early, before any I/O
    parse_yyyymm(args.span_start)
    parse_yyyymm(args.span_end)
    if args.baseline_end < args.baseline_start:
        ap.error("--baseline-end precedes --baseline-start")

    build(args)
    print("make_climatology_series: done")


if __name__ == "__main__":
    main()
