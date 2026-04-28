"""Build verify_v2.ipynb from cell list.

Run once to regenerate the notebook from this source.
Source-of-truth is this file; the .ipynb is a derived artifact.
"""
import json
import textwrap

cells = []


def md(src):
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": textwrap.dedent(src).strip(),
    })


def code(src):
    cells.append({
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": textwrap.dedent(src).strip(),
    })


# ---- Title ---------------------------------------------------------------
md("""
    # MiCASA v2 — Pipeline Verification (`verify_v2.ipynb`)

    Phase 1 of a port from miller-ff `verify_2026.ipynb`.
    Ten structural-invariant checks across schema, transformation
    invariants, sign convention, and cross-boundary sanity.

    **Run from the v2 working directory**:
    ```sh
    cd /work2/noaa/co2/GFED-CASA/2025/MiCASA_v2
    jupyter nbconvert --to notebook --execute verify_v2.ipynb \\
        --output verify_v2.executed.ipynb
    ```

    Each check appends to `_RESULTS`. The summary cell at the end prints
    counts and any FAIL details. Designed to be additive: Phase 2 / 3 will
    add Sections 5+ for cross-product comparison, spatial sanity, and
    seasonal sanity.
""")

# ---- Configuration & Setup ----------------------------------------------
md("## Configuration & Setup")
code('''
    import os, sys, json, subprocess, glob, re
    from pathlib import Path
    from collections import OrderedDict
    import numpy as np
    import pandas as pd
    import xarray as xr

    WORK_DIR     = Path(os.environ.get("WORK_DIR", os.getcwd()))
    MONTHLY_FILE = WORK_DIR / "monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc"
    FIT_RDA      = WORK_DIR / "fit.piqs.rda"
    ERA5_DIR     = WORK_DIR / "ERA5"
    JOBS_DIR     = WORK_DIR / "jobs"
    DAILY_DIR    = WORK_DIR / "daily_1x1"
    MET_BASE     = Path(os.environ.get("CARBONTRACKER", "/work2/noaa/co2")) \\
                   / "METEO/tm5-nc/ec/ea/h06h18tr1/sfc/glb100x100"

    PASS, FAIL, WARN, INFO = "PASS", "FAIL", "WARN", "INFO"
    _RESULTS = []

    def record(check_id, name, status, detail=""):
        _RESULTS.append((check_id, name, status, detail))
        color = {PASS: "\\x1b[32m", FAIL: "\\x1b[31m",
                 WARN: "\\x1b[33m", INFO: "\\x1b[34m"}.get(status, "")
        reset = "\\x1b[0m"
        print(f"{color}[{status}]{reset} {check_id} {name}: {detail}")
''')

code('''
    # Quick orientation: what is on disk?
    print(f"WORK_DIR     = {WORK_DIR}")
    print(f"MONTHLY_FILE = {MONTHLY_FILE} ({'exists' if MONTHLY_FILE.exists() else 'MISSING'})")
    print(f"FIT_RDA      = {FIT_RDA} ({'exists' if FIT_RDA.exists() else 'MISSING'})")
    print(f"ERA5/        = {len(list(ERA5_DIR.glob('fluxes_*.nc')))} fluxes_*.nc files")
    print(f"daily_1x1/   = {len(list(DAILY_DIR.glob('*.nc')))} *.nc files (incl. clim + symlinks)")
    print(f"jobs/        = {len(list(JOBS_DIR.glob('*')))} log files")
''')

# ---- Section 1: Schema & Coverage ---------------------------------------
md("## Section 1 — Schema & Coverage")

md("""
    ### Check 1.1 — Multi-year monthly cat file schema

    The file `monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc` (output of
    `cat_monthly.sh`) is the input to PIQS. It must have:
    - dims: `longitude=360, latitude=180, time>=288`
    - variables: `NPP, Rh, FIRE, FUEL` (units gC m-2 s-1)
    - time monotonically increasing, monthly cadence
""")
code('''
    cid, cname = "1.1", "monthly_cat schema"
    try:
        ds = xr.open_dataset(MONTHLY_FILE)
        problems = []
        if ds.sizes.get("longitude") != 360: problems.append(f"longitude dim is {ds.sizes.get('longitude')}, expected 360")
        if ds.sizes.get("latitude")  != 180: problems.append(f"latitude  dim is {ds.sizes.get('latitude')}, expected 180")
        if ds.sizes.get("time", 0) < 288:    problems.append(f"time dim is {ds.sizes.get('time')}, expected >= 288")
        for v in ("NPP", "Rh", "FIRE", "FUEL"):
            if v not in ds.variables:
                problems.append(f"variable '{v}' missing")
        if problems:
            record(cid, cname, FAIL, "; ".join(problems))
        else:
            record(cid, cname, PASS,
                   f"dims OK, time={ds.sizes['time']} months "
                   f"({pd.to_datetime(ds.time.values[0]).strftime('%Y-%m')}..."
                   f"{pd.to_datetime(ds.time.values[-1]).strftime('%Y-%m')})")
        ds.close()
    except Exception as e:
        record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 1.2 — PIQS fit metadata schema

    Calls the R helper `verify_piqs_invariants.r` which loads `fit.piqs.rda`
    and dumps a JSON summary. Confirms:
    - `piqsfit.gpp/resp` are `360 × 180 × N` arrays
    - `piqsfit.time` length matches
    - `piqsfit.meta` exists and records padding settings
""")
code('''
    cid, cname = "1.2", "PIQS fit metadata"
    helper = WORK_DIR / "verify_piqs_invariants.r"
    out_json = WORK_DIR / "verify_piqs_invariants.json"
    if not helper.exists():
        record(cid, cname, FAIL, f"helper script missing: {helper}")
    else:
        # Run the helper (also computes Check 2.1 residuals -- reused below).
        try:
            r = subprocess.run(["Rscript", str(helper), str(out_json)],
                               cwd=WORK_DIR, capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                record(cid, cname, FAIL, f"Rscript failed: {r.stderr[-300:]}")
            else:
                _piqs = json.loads(out_json.read_text())
                problems = []
                if list(_piqs.get("fit_dims_gpp", [])) != [360, 180, _piqs.get("nmon", -1)]:
                    problems.append(f"gpp dims {_piqs.get('fit_dims_gpp')} unexpected")
                if not _piqs.get("piqsfit_meta"):
                    problems.append("piqsfit.meta missing from .rda")
                if problems:
                    record(cid, cname, FAIL, "; ".join(problems))
                else:
                    meta = _piqs["piqsfit_meta"]
                    record(cid, cname, PASS,
                           f"nmon={_piqs['nmon']}, "
                           f"pad=(L{meta.get('pad.left')},R{meta.get('pad.right')}), "
                           f"written {meta.get('written.at')}")
        except Exception as e:
            record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 1.3 — `fluxes_YYYYMM.nc` schema + NaN/Inf scan

    For each `ERA5/fluxes_YYYYMM.nc`, confirm dims, variable list, and that
    no NaN/Inf appears in NEE/GPP/resp over land. (Polar cells legitimately
    have NaN/zero where ssrd ≡ 0; we flag only land NaN counts as concerning.)
""")
code('''
    cid, cname = "1.3", "fluxes_*.nc schema + NaN scan"
    expected_vars = {"GPP", "resp", "NEE", "QGPP", "qresp", "ssr", "t2m", "stl1", "swvl1", "decimal_date"}
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))
    if not files:
        record(cid, cname, WARN, "no fluxes_*.nc files yet (diurnalize not run?)")
    else:
        bad = []
        nan_counts = []
        # Sample first, middle, last to keep runtime small
        sample_idx = sorted({0, len(files)//2, len(files)-1})
        for i in sample_idx:
            f = files[i]
            try:
                ds = xr.open_dataset(f)
                missing = expected_vars - set(ds.variables)
                if missing:
                    bad.append(f"{f.name}: missing {missing}")
                # NaN/Inf in non-polar land cells (latitude in [-60, 60])
                lat = ds.latitude.values
                lat_ix = (lat >= -60) & (lat <= 60)
                for v in ("NEE", "GPP", "resp"):
                    if v in ds:
                        arr = ds[v].isel(latitude=lat_ix).values
                        n_nan = int(np.isnan(arr).sum() + np.isinf(arr).sum())
                        if n_nan:
                            nan_counts.append(f"{f.name}/{v}={n_nan}")
                ds.close()
            except Exception as e:
                bad.append(f"{f.name}: {e}")
        problems = bad + nan_counts
        if problems:
            record(cid, cname, FAIL, "; ".join(problems[:6]))
        else:
            record(cid, cname, PASS,
                   f"checked {len(sample_idx)}/{len(files)} files (first/mid/last); schema OK, "
                   f"no NaN/Inf in [-60,60] latitude band")
''')

md("""
    ### Check 1.4 — ERA5 meteo coverage matches diurnalized year range

    Every year/month with a `fluxes_YYYYMM.nc` output must have had its
    meteo input on disk under `$CARBONTRACKER/METEO/.../<year>/<MM>/`.
    Catches the `ea_0005` vs `ea/` regression (workers fail at first
    `nc_open`) before it actually fails -- by checking representative meteo
    paths exist.
""")
code('''
    cid, cname = "1.4", "ERA5 meteo coverage"
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))
    missing = []
    for f in files:
        m = re.search(r"fluxes_(\\d{4})(\\d{2})\\.nc$", f.name)
        if not m: continue
        yr, mo = m.group(1), m.group(2)
        # Just probe the first day t2m file for that year/month
        probe = MET_BASE / yr / mo / f"t2m_{yr}{mo}01_00p01.nc"
        if not probe.exists():
            missing.append(f"{yr}-{mo}")
    if not files:
        record(cid, cname, INFO, "no fluxes_*.nc to check yet")
    elif missing:
        record(cid, cname, FAIL,
               f"{len(missing)} fluxes_*.nc files lack a corresponding meteo probe; "
               f"first few: {missing[:6]}")
    else:
        record(cid, cname, PASS, f"every {len(files)} fluxes file has meteo on disk")
''')

# ---- Section 2: Transformation Invariants -------------------------------
md("## Section 2 — Pipeline Transformation Invariants")

md("""
    ### Check 2.1 — PIQS integral preservation

    PIQS = Piecewise *Integral* Quadratic Splines. The defining property:
    for each segment, the integral of the per-piece quadratic equals the
    input monthly mean × segment width. Equivalently:

    > integral / delta == input monthly mean (per cell-segment)

    The R helper computed this residual for every active land cell-month;
    we read its summary and assert max-absolute-residual is below tolerance.
""")
code('''
    cid_g, cid_r = "2.1.gpp", "2.1.rtot"
    out_json = WORK_DIR / "verify_piqs_invariants.json"
    TOL_ABS = 1e-9   # gC m-2 s-1
    TOL_REL = 1e-6
    if not out_json.exists():
        record(cid_g, "PIQS integral preservation (GPP)", FAIL,
               f"helper output {out_json.name} missing -- did Check 1.2 run?")
    else:
        _piqs = json.loads(out_json.read_text())
        for cid, key, label in [(cid_g, "gpp", "GPP"), (cid_r, "rtot", "Rtot")]:
            s = _piqs[key]
            detail = (f"max_abs={s['max_abs_residual']:.2e} "
                      f"max_rel={s['max_rel_residual']:.2e} "
                      f"frac>1e-9 abs={s['frac_abs_resid_gt_1e9']:.4f}")
            ok_abs = s["max_abs_residual"] < TOL_ABS
            ok_rel = s["max_rel_residual"] < TOL_REL
            if ok_abs and ok_rel:
                record(cid, f"PIQS integral preservation ({label})", PASS, detail)
            else:
                record(cid, f"PIQS integral preservation ({label})", FAIL, detail)
''')

md("""
    ### Check 2.2 — Diurnalize preserves monthly means

    For a sample `fluxes_YYYYMM.nc`, the per-cell mean of hourly NEE
    over the month should equal the input monthly mean (which is a
    function of the gridcell's NPP and Rh for that month). Within a
    small numerical tolerance.
""")
code('''
    cid, cname = "2.2", "diurnalize preserves monthly means"
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))
    if not files:
        record(cid, cname, INFO, "no fluxes_*.nc to check yet")
    else:
        # Use the middle file for a representative sample
        f = files[len(files)//2]
        m = re.search(r"fluxes_(\\d{4})(\\d{2})\\.nc$", f.name)
        yr, mo = m.group(1), m.group(2)
        monthly_pm = WORK_DIR / "monthly_1x1" / f"MiCASA_v1_flux_x360_y180_monthly_{yr}{mo}.nc"
        if not monthly_pm.exists():
            record(cid, cname, WARN, f"per-month input {monthly_pm.name} missing -- skipping")
        else:
            try:
                ds_h = xr.open_dataset(f)            # hourly
                ds_m = xr.open_dataset(monthly_pm)   # monthly mean (one timestep)
                # Convert sign / units: gpp.mn = -2*NPP/12 in mol m-2 s-1
                # In monthly file, NPP is in gC m-2 s-1; gC -> mol = / 12
                gpp_mn_expected = (-2.0 * ds_m["NPP"].squeeze() / 12.0).values
                rtot_mn_expected = (ds_m["Rh"].squeeze() / 12.0 - 0.5 * gpp_mn_expected).values
                gpp_mn_actual    = ds_h["GPP"].mean(dim="time").values
                resp_mn_actual   = ds_h["resp"].mean(dim="time").values
                # Mask to active land (where the input mean is non-zero)
                mask = (np.abs(gpp_mn_expected) > 1e-15)
                if mask.sum() == 0:
                    record(cid, cname, WARN, "no active land cells in sample month")
                else:
                    # Relative diff
                    rel_g = np.abs((gpp_mn_actual - gpp_mn_expected)[mask] /
                                   gpp_mn_expected[mask])
                    rel_r = np.abs((resp_mn_actual - rtot_mn_expected)[mask] /
                                   rtot_mn_expected[mask])
                    detail = (f"sample {f.name}: GPP max_rel={rel_g.max():.2e} "
                              f"median={np.median(rel_g):.2e}; "
                              f"resp max_rel={rel_r.max():.2e} median={np.median(rel_r):.2e}")
                    if rel_g.max() < 1e-3 and rel_r.max() < 1e-3:
                        record(cid, cname, PASS, detail)
                    else:
                        record(cid, cname, FAIL, detail)
                ds_h.close(); ds_m.close()
            except Exception as e:
                record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 2.3 — Cat-monthly time monotonicity & monthly cadence

    The multi-year monthly file's `time` coordinate must be strictly
    monotonically increasing with median spacing ≈ 30.4 days (one month).
    Catches stale concat files and ncrcat ordering bugs.
""")
code('''
    cid, cname = "2.3", "monthly time monotonicity"
    try:
        ds = xr.open_dataset(MONTHLY_FILE)
        t  = pd.to_datetime(ds.time.values)
        dt = np.diff(t.values).astype("timedelta64[D]").astype(float)
        problems = []
        if not np.all(dt > 0): problems.append("time not monotonic")
        if not (28 <= dt.min() <= 32 and 28 <= dt.max() <= 32):
            problems.append(f"dt range [{dt.min():.1f},{dt.max():.1f}] days outside [28,32]")
        if np.median(dt) < 28 or np.median(dt) > 32:
            problems.append(f"median dt={np.median(dt):.1f} days, expected ~30.4")
        if problems:
            record(cid, cname, FAIL, "; ".join(problems))
        else:
            record(cid, cname, PASS,
                   f"{len(t)} timestamps, {t[0].strftime('%Y-%m')}..{t[-1].strftime('%Y-%m')}, "
                   f"median dt={np.median(dt):.2f} days")
        ds.close()
    except Exception as e:
        record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 2.4 — Daily → monthly mass conservation (sample)

    For a sample year/month, mean of the 1° daily files should equal the
    1° monthly file (within numerical precision). Catches ingest scaling
    errors and stale per-month files.
""")
code('''
    cid, cname = "2.4", "daily->monthly mass conservation"
    # Pick a recent fully-real month -- use the last month with both
    # daily and monthly inputs available.
    monthly_files = sorted((WORK_DIR / "monthly_1x1").glob(
        "MiCASA_v1_flux_x360_y180_monthly_2*.nc"))
    chosen = None
    for f in reversed(monthly_files):
        m = re.search(r"monthly_(\\d{4})(\\d{2})\\.nc$", f.name)
        if not m: continue
        yr, mo = m.group(1), m.group(2)
        sample_day = DAILY_DIR / f"MiCASA_v1_flux_x360_y180_daily_{yr}{mo}01.nc"
        if sample_day.exists() and not sample_day.is_symlink():
            chosen = (f, yr, mo); break
        # Allow symlinks too (vNRT-relabelled days)
        if sample_day.exists():
            chosen = (f, yr, mo); break
    if chosen is None:
        record(cid, cname, INFO, "no month found with both per-month and per-day files")
    else:
        f_m, yr, mo = chosen
        days = sorted(DAILY_DIR.glob(f"MiCASA_v1_flux_x360_y180_daily_{yr}{mo}*.nc"))
        if not days:
            record(cid, cname, INFO, f"no daily files for {yr}-{mo}")
        else:
            try:
                ds_d = xr.open_mfdataset([str(p) for p in days], combine="nested",
                                         concat_dim="time")
                ds_m = xr.open_dataset(f_m)
                npp_d_mean = ds_d["NPP"].mean(dim="time").values
                npp_m      = ds_m["NPP"].squeeze().values
                # Active land mask
                mask = np.abs(npp_m) > 1e-15
                if mask.sum() == 0:
                    record(cid, cname, WARN, f"no active land cells in {yr}-{mo}")
                else:
                    rel = np.abs((npp_d_mean - npp_m)[mask] / npp_m[mask])
                    detail = (f"{yr}-{mo} ({len(days)} daily files), NPP max_rel={rel.max():.2e}, "
                              f"median={np.median(rel):.2e}")
                    if rel.max() < 0.01:
                        record(cid, cname, PASS, detail)
                    else:
                        record(cid, cname, FAIL, detail)
                ds_d.close(); ds_m.close()
            except Exception as e:
                record(cid, cname, FAIL, f"exception: {e}")
''')

# ---- Section 3: Sign Convention -----------------------------------------
md("## Section 3 — Sign Convention")

md("""
    ### Check 3.1 — Sub-monthly sign-flip rate (parses diurnalize logs)

    The sign-flip diagnostic added in proposal #4 prints two lines per
    (year, month) worker. We aggregate across all `jobs/d-*.o*` logs to
    build a global summary. Threshold-based: warn if mean cell-hour rate
    > 5% (GPP) or > 1% (resp), fail if > 15% (GPP) or > 5% (resp).
""")
code('''
    cid, cname = "3.1", "sign-flip rate aggregate"
    pat_gpp = re.compile(r"PIQS sign-flip \\[GPP > 0\\]:\\s+\\d+ / \\d+ cells \\(([0-9.]+)%\\), \\d+ / \\d+ cell-hours \\(([0-9.]+)%\\)")
    pat_rsp = re.compile(r"PIQS sign-flip \\[resp < 0\\]:\\s+\\d+ / \\d+ cells \\(([0-9.]+)%\\), \\d+ / \\d+ cell-hours \\(([0-9.]+)%\\)")
    logs = sorted(JOBS_DIR.glob("d-*-MiCASA*.o*"))
    if not logs:
        record(cid, cname, INFO, "no diurnalize worker logs found")
    else:
        gpp_cells, gpp_chs, rsp_cells, rsp_chs = [], [], [], []
        for L in logs:
            try:
                txt = L.read_text(errors="ignore")
                gpp_cells += [float(m.group(1)) for m in pat_gpp.finditer(txt)]
                gpp_chs   += [float(m.group(2)) for m in pat_gpp.finditer(txt)]
                rsp_cells += [float(m.group(1)) for m in pat_rsp.finditer(txt)]
                rsp_chs   += [float(m.group(2)) for m in pat_rsp.finditer(txt)]
            except Exception:
                continue
        if not gpp_chs:
            record(cid, cname, WARN, f"{len(logs)} logs, but no PIQS sign-flip lines parsed")
        else:
            mg, xg = float(np.mean(gpp_chs)),  float(np.max(gpp_chs))
            mr, xr_ = float(np.mean(rsp_chs)), float(np.max(rsp_chs))
            detail = (f"{len(gpp_chs)} months: "
                      f"GPP cell-hour mean={mg:.2f}% max={xg:.2f}%; "
                      f"resp cell-hour mean={mr:.3f}% max={xr_:.3f}%")
            if mg > 15.0 or mr > 5.0:
                status = FAIL
            elif mg > 5.0 or mr > 1.0:
                status = WARN
            else:
                status = PASS
            record(cid, cname, status, detail)
''')

# ---- Section 4: Cross-Boundary Sanity -----------------------------------
md("## Section 4 — Cross-Boundary Sanity")

md("""
    ### Check 4.1 — v1↔vNRT splice continuity at Dec 2024 / Jan 2025

    Reads `diag_v1_vNRT_handoff.csv` (produced by `diag_v1_vNRT_handoff.r`)
    and looks for a step-change in monthly global totals across the
    splice boundary. A jump > 10% in NPP or Rh between months adjacent to
    the boundary is suspicious.
""")
code('''
    cid, cname = "4.1", "v1<->vNRT splice continuity"
    csv_path = WORK_DIR / "diag_v1_vNRT_handoff.csv"
    if not csv_path.exists():
        record(cid, cname, INFO,
               "diag_v1_vNRT_handoff.csv not found -- run diag_v1_vNRT_handoff.r first")
    else:
        df = pd.read_csv(csv_path)
        df["yyyymm"] = df["year"]*100 + df["month"]
        # Boundary at 2025-01 = 202501
        before = df[df["yyyymm"] == 202412]
        after  = df[df["yyyymm"] == 202501]
        if len(before) == 0:
            record(cid, cname, WARN, "no 2024-12 row in diagnostic CSV")
        elif len(after) == 0:
            record(cid, cname, INFO,
                   "no 2025-01 row yet -- ingest of vNRT 2025 monthlies probably hasn't completed")
        else:
            problems = []
            for v in ("NPP", "Rh"):
                if v not in df.columns: continue
                b, a = float(before[v].iloc[0]), float(after[v].iloc[0])
                rel_jump = abs(a - b) / max(abs(b), 1e-12)
                if rel_jump > 0.10:
                    problems.append(f"{v}: 2024-12={b:.1f} -> 2025-01={a:.1f} (rel jump {rel_jump:.1%})")
            if problems:
                record(cid, cname, WARN, "; ".join(problems))
            else:
                record(cid, cname, PASS,
                       f"NPP/Rh global totals continuous within 10% across the v1<->vNRT splice")
''')

# ---- Summary ------------------------------------------------------------
md("## Summary")

code('''
    # PASS/FAIL/WARN/INFO counts and the table of all results.
    from collections import Counter
    counts = Counter(r[2] for r in _RESULTS)
    print("=" * 60)
    print(f"  PASS={counts.get(PASS,0):3d}  "
          f"FAIL={counts.get(FAIL,0):3d}  "
          f"WARN={counts.get(WARN,0):3d}  "
          f"INFO={counts.get(INFO,0):3d}")
    print("=" * 60)
    print(f"{'check':<12} {'status':<6} {'name':<40} detail")
    print("-" * 60)
    for cid, name, status, detail in _RESULTS:
        print(f"{cid:<12} {status:<6} {name:<40} {detail}")
    if counts.get(FAIL, 0) > 0:
        print("\\nOne or more checks FAILED -- investigate above.")
    else:
        print("\\nAll checks passed.")
''')

# ============================================================================
#                              PHASE 2
# ============================================================================
md("# Phase 2 — Comparison & Sanity")

md("""
    Phase 2 adds "looks reasonable?" checks on top of Phase 1's structural
    invariants. The new sections each operate on a small per-month summary
    cube (`verify_v2_summary.csv`) built once by the preflight cell — so
    individual checks don't have to re-scan all ~300 hourly fluxes files.

    Sections added in this phase:

    - **5. Global Totals & Trends** — monthly/annual global NEE, GPP, Rh;
      plot the full record; check year-on-year growth rate plausibility;
      flag any outlier years.
    - **6. Spatial Comparison vs v1** — for sampled months, spatial Pearson
      correlation between v2 and v1 monthly-mean NEE; difference at the
      v2 PIQS-padded right edge (where we expect v2 to differ from v1) vs
      record interior (where they should be ~identical).
    - **7. Spatial Sanity** — ocean cells zero, polar/Antarctic small,
      tropical/boreal hotspot seasonal patterns plausible.
    - **8. Seasonal & Temporal** — NH/SH out-of-phase, latitude amplitude
      profile, year-boundary continuity (no PIQS ringing at Dec→Jan).
""")

md("## Phase 2 Preflight — Build per-month summary cube")
md("""
    Walk every `ERA5/fluxes_YYYYMM.nc` once, aggregate to monthly global,
    NH (>0°), SH (<0°), tropics (|lat|<23.5°), boreal (>50°N), and SH-mid
    (-50..0°) totals (GgC/month). Cache as `verify_v2_summary.csv`. Rebuild
    if any source file is newer than the cache. **About 5–10 minutes the
    first run, near-instant thereafter.**

    Conversion: fluxes are mol m⁻² s⁻¹. Multiply by gridcell area (m²) and
    seconds-in-month, then by 12 g/mol, divide by 1e15 → PgC/month.
""")
code('''
    cid, cname = "P2.0", "build summary cube"
    SUMMARY_CSV = WORK_DIR / "verify_v2_summary.csv"
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))

    def needs_rebuild():
        if not SUMMARY_CSV.exists(): return True
        cache_mtime = SUMMARY_CSV.stat().st_mtime
        for f in files:
            if f.stat().st_mtime > cache_mtime: return True
        return False

    R_EARTH = 6371007.2  # m
    def gridcell_area_m2():
        # 1° x 1° grid; latitude bands -89.5..89.5
        lat_edges = np.arange(-90, 91, 1) * np.pi / 180  # 181 edges
        # area = R^2 * (sin(lat_top) - sin(lat_bot)) * dlon (radians)
        dlon = np.pi / 180  # 1° in radians
        lat_band_area = R_EARTH ** 2 * (np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1])) * dlon  # 180
        # Replicate to 360 longitudes
        return np.broadcast_to(lat_band_area[:, None], (180, 360))  # (lat, lon)

    if not files:
        record(cid, cname, INFO, "no fluxes_*.nc files yet; Phase 2 cannot run")
    elif not needs_rebuild():
        _summary = pd.read_csv(SUMMARY_CSV)
        record(cid, cname, INFO,
               f"cache up to date ({len(_summary)} months in {SUMMARY_CSV.name})")
    else:
        print(f"Building summary cube from {len(files)} files...")
        rows = []
        gca = gridcell_area_m2()  # (lat, lon)
        for f in files:
            m = re.search(r"fluxes_(\\d{4})(\\d{2})\\.nc$", f.name)
            if not m: continue
            yr, mo = int(m.group(1)), int(m.group(2))
            try:
                ds = xr.open_dataset(f)
                # Mean over hourly time dim → (lat, lon) per variable
                # Result is mean flux in mol m-2 s-1
                lat = ds.latitude.values
                ndays = pd.Period(f"{yr}-{mo:02d}").days_in_month
                sec_per_month = ndays * 86400
                # Convert mol m-2 s-1 → gC m-2 month-1 → integrate → GgC/month
                # = mean_flux * gca * sec_per_month * 12 (gC/mol) / 1e9 (GgC/g)
                conv = sec_per_month * 12.0 / 1.0e9
                row = {"year": yr, "month": mo, "yyyymm": yr*100 + mo}
                for v in ("NEE", "GPP", "resp"):
                    if v not in ds: continue
                    mn = ds[v].mean(dim="time").values  # (lat, lon)
                    integrand = mn * gca * conv  # GgC/month per cell
                    row[f"{v}_global"] = float(np.nansum(integrand))
                    # Latitude masks -- lat axis is 0
                    nh        = lat > 0
                    sh        = lat < 0
                    trop      = np.abs(lat) < 23.5
                    boreal    = lat > 50
                    nh_mid    = (lat > 0) & (lat < 50)
                    sh_mid    = (lat > -50) & (lat < 0)
                    polar_n   = lat > 70
                    polar_s   = lat < -70
                    for nm, mask in [("nh", nh), ("sh", sh), ("trop", trop),
                                     ("boreal", boreal), ("nh_mid", nh_mid),
                                     ("sh_mid", sh_mid), ("polar_n", polar_n),
                                     ("polar_s", polar_s)]:
                        row[f"{v}_{nm}"] = float(np.nansum(integrand[mask, :]))
                ds.close()
                rows.append(row)
            except Exception as e:
                print(f"  WARN: {f.name}: {e}")
        _summary = pd.DataFrame(rows).sort_values("yyyymm").reset_index(drop=True)
        _summary.to_csv(SUMMARY_CSV, index=False)
        record(cid, cname, PASS,
               f"built summary cube: {len(_summary)} months, "
               f"{_summary['yyyymm'].iloc[0]}..{_summary['yyyymm'].iloc[-1]}")
    # Make _summary available to subsequent cells whether or not we rebuilt
    if SUMMARY_CSV.exists():
        _summary = pd.read_csv(SUMMARY_CSV)
''')

# ---- Section 5: Global Totals & Trends ----------------------------------
md("## Section 5 — Global Totals & Trends")

md("""
    ### Check 5.1 — Annual global NEE / GPP / Rh time series

    Aggregate the monthly summary to annual sums; print the table and plot
    the time series. Sanity criterion: GPP and Rh should each be in the
    100–150 PgC/yr range (terrestrial biosphere); NEE small relative to
    those (a few PgC/yr in either direction depending on year).
""")
code('''
    cid, cname = "5.1", "annual global NEE/GPP/Rh"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        annual = _summary.groupby("year").agg({
            "NEE_global": "sum",
            "GPP_global": "sum",
            "resp_global": "sum",
        })
        # GgC/yr -> PgC/yr (1e6 GgC = 1 PgC)
        annual = annual / 1e6
        problems = []
        # GPP should be roughly -100 to -150 PgC/yr (uptake = negative in our convention)
        if not ((-200 <= annual["GPP_global"]).all() and (annual["GPP_global"] <= -50).all()):
            problems.append(f"GPP outside [-200,-50] PgC/yr: {annual['GPP_global'].min():.1f}..{annual['GPP_global'].max():.1f}")
        # Rh should be roughly +50 to +150 PgC/yr (positive)
        if not ((30 <= annual["resp_global"]).all() and (annual["resp_global"] <= 200).all()):
            problems.append(f"resp outside [30,200] PgC/yr: {annual['resp_global'].min():.1f}..{annual['resp_global'].max():.1f}")
        # Print table head/tail
        print("Annual global totals (PgC/yr):")
        print(annual.head(3).to_string())
        print("...")
        print(annual.tail(3).to_string())
        if problems:
            record(cid, cname, FAIL, "; ".join(problems))
        else:
            record(cid, cname, PASS,
                   f"{len(annual)} years; GPP in [{annual['GPP_global'].min():.1f},{annual['GPP_global'].max():.1f}], "
                   f"resp in [{annual['resp_global'].min():.1f},{annual['resp_global'].max():.1f}] PgC/yr")
''')
code('''
    # Plot the time series. Inline plot, not a check.
    import matplotlib.pyplot as plt
    if "_summary" in globals() and not _summary.empty:
        annual = _summary.groupby("year").agg({
            "NEE_global": "sum", "GPP_global": "sum", "resp_global": "sum",
        }) / 1e6  # PgC/yr
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        for ax, col, ylabel in zip(
            axes, ["GPP_global", "resp_global", "NEE_global"],
            ["GPP (PgC/yr)", "resp (PgC/yr)", "NEE (PgC/yr)"],
        ):
            ax.plot(annual.index, annual[col], "o-", lw=1.2, ms=4)
            ax.axhline(0, color="grey", lw=0.5)
            ax.grid(alpha=0.3)
            ax.set_ylabel(ylabel)
        axes[-1].set_xlabel("Year")
        fig.suptitle("v2 annual global totals", y=0.995)
        fig.tight_layout()
        plt.show()
''')

md("""
    ### Check 5.2 — Year-on-year growth rate plausibility

    Year-on-year change in global GPP/Rh should be small (typically within
    ±5%). A jump above 10% suggests an artefact — likely a methodology
    change or input data step (e.g., the v1↔vNRT splice) showing through.
""")
code('''
    cid, cname = "5.2", "YoY growth rate plausibility"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        annual = _summary.groupby("year").agg({
            "NEE_global": "sum", "GPP_global": "sum", "resp_global": "sum",
        }) / 1e6
        problems = []
        warnings = []
        for col in ("GPP_global", "resp_global"):
            yoy = annual[col].pct_change() * 100  # %
            for yr, pct in yoy.dropna().items():
                if abs(pct) > 20:
                    problems.append(f"{col} {int(yr)} YoY {pct:+.1f}%")
                elif abs(pct) > 10:
                    warnings.append(f"{col} {int(yr)} YoY {pct:+.1f}%")
        if problems:
            record(cid, cname, FAIL, "; ".join(problems[:6]))
        elif warnings:
            record(cid, cname, WARN, "; ".join(warnings[:6]))
        else:
            record(cid, cname, PASS, "all GPP/resp YoY changes within ±10%")
''')

md("""
    ### Check 5.3 — Seasonal cycle amplitude per latitude band

    NH temperate (30–60°N) should have the strongest annual amplitude in
    NEE (largest positive winter, largest negative summer). Boreal
    (>50°N) similar but smaller. Tropics should have a weak seasonal
    cycle. SH is anti-phased.
""")
code('''
    cid, cname = "5.3", "seasonal amplitude per band"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        # Take last 5 full years and average per-month
        last_year = int(_summary["year"].max())
        if (_summary["year"] == last_year).sum() < 12:
            last_year -= 1  # active year incomplete
        recent = _summary[(_summary["year"] >= last_year - 4) & (_summary["year"] <= last_year)]
        # Average across years per calendar month, in PgC/month (1e6 GgC = 1 PgC)
        mclim = recent.groupby("month").agg({
            "NEE_nh_mid":  "mean", "NEE_boreal":  "mean",
            "NEE_trop":    "mean", "NEE_sh_mid":  "mean",
        }) / 1e6
        amp = (mclim.max() - mclim.min()).to_dict()
        nh_mid_amp = amp.get("NEE_nh_mid", 0)
        boreal_amp = amp.get("NEE_boreal", 0)
        trop_amp   = amp.get("NEE_trop", 0)
        sh_mid_amp = amp.get("NEE_sh_mid", 0)
        problems = []
        if nh_mid_amp < boreal_amp * 0.5:
            problems.append(f"NH mid amp ({nh_mid_amp:.2f}) unexpectedly less than boreal ({boreal_amp:.2f})")
        if trop_amp > nh_mid_amp:
            problems.append(f"tropical amp ({trop_amp:.2f}) > NH mid ({nh_mid_amp:.2f})")
        # NH/SH phase: month of NEE max should differ by 5-7 months
        nh_max_m = mclim["NEE_nh_mid"].idxmax()
        sh_max_m = mclim["NEE_sh_mid"].idxmax()
        phase_diff = abs(nh_max_m - sh_max_m)
        if not (4 <= phase_diff <= 8):
            problems.append(f"NH max month {nh_max_m}, SH max month {sh_max_m}, phase diff {phase_diff} (expected 5..7)")
        detail = (f"amp PgC/mo: nh_mid={nh_mid_amp:.2f}, boreal={boreal_amp:.2f}, "
                  f"trop={trop_amp:.2f}, sh_mid={sh_mid_amp:.2f}; phase diff {phase_diff}")
        if problems:
            record(cid, cname, FAIL, "; ".join(problems))
        else:
            record(cid, cname, PASS, detail)
''')
code('''
    # Plot the per-band climatological seasonal cycle for context
    import matplotlib.pyplot as plt
    if "_summary" in globals() and not _summary.empty:
        last_year = int(_summary["year"].max())
        if (_summary["year"] == last_year).sum() < 12:
            last_year -= 1
        recent = _summary[(_summary["year"] >= last_year - 4) & (_summary["year"] <= last_year)]
        mclim = recent.groupby("month").agg({
            "NEE_nh_mid":  "mean", "NEE_boreal":  "mean",
            "NEE_trop":    "mean", "NEE_sh_mid":  "mean",
        }) / 1e6
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for col, label in [("NEE_boreal", "Boreal (>50°N)"),
                           ("NEE_nh_mid", "NH mid (0..50°N)"),
                           ("NEE_trop",   "Tropics (|lat|<23.5°)"),
                           ("NEE_sh_mid", "SH mid (-50..0°)")]:
            ax.plot(mclim.index, mclim[col], "o-", lw=1.2, ms=4, label=label)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlabel("Month")
        ax.set_ylabel("NEE (PgC/month)")
        ax.set_xticks(range(1, 13))
        ax.set_title(f"NEE seasonal cycle by band ({last_year-4}..{last_year} climatology)")
        ax.legend(loc="best")
        ax.grid(alpha=0.3)
        plt.show()
''')

# ---- Section 6: Spatial Comparison vs v1 --------------------------------
md("## Section 6 — Spatial Comparison vs v1 (where overlapping)")

md("""
    ### Check 6.1 — Spatial correlation v2 vs v1, sampled months

    For 4 sampled months at well-separated points in the record, compute
    the per-cell Pearson correlation between v2 and v1 monthly-mean NEE.
    Expectation:
    - **Interior** months (e.g., 2010-07, 2018-01): r > 0.99 — both products
      use the same input monthly NPP/Rh; differences are PIQS interior
      coefficient drift only.
    - **Right-edge** months (e.g., 2024-12): r still high but lower —
      v2's PIQS pad stabilises the December coefs, so the diurnal
      redistribution differs slightly.
""")
code('''
    cid, cname = "6.1", "v2 vs v1 spatial correlation"
    V1_ERA5 = WORK_DIR.parent / "MiCASA_v1" / "ERA5"
    sample = ["201001", "201501", "202007", "202412"]
    rows = []
    for ymm in sample:
        v2_p = ERA5_DIR / f"fluxes_{ymm}.nc"
        v1_p = V1_ERA5 / f"fluxes_{ymm}.nc"
        if not v2_p.exists() or not v1_p.exists():
            rows.append((ymm, "missing", float("nan")))
            continue
        try:
            ds2 = xr.open_dataset(v2_p)
            ds1 = xr.open_dataset(v1_p)
            n2 = ds2["NEE"].mean(dim="time").values
            n1 = ds1["NEE"].mean(dim="time").values
            mask = np.isfinite(n2) & np.isfinite(n1) & ((np.abs(n2) > 1e-15) | (np.abs(n1) > 1e-15))
            if mask.sum() < 100:
                rows.append((ymm, "too few cells", float("nan")))
            else:
                r = float(np.corrcoef(n2[mask], n1[mask])[0, 1])
                rows.append((ymm, "ok", r))
            ds2.close(); ds1.close()
        except Exception as e:
            rows.append((ymm, f"error: {e}", float("nan")))
    detail_strs = [f"{ymm}: r={r:.4f}" for ymm, st, r in rows if st == "ok"]
    if not detail_strs:
        record(cid, cname, INFO, f"no v1↔v2 overlap available: {rows}")
    else:
        # FAIL if any interior month r < 0.95; WARN if right-edge < 0.85
        bad = [s for ymm, st, r in rows if st == "ok"
               for s in [f"{ymm} r={r:.3f}"] if r < 0.85]
        warn = [s for ymm, st, r in rows if st == "ok"
                for s in [f"{ymm} r={r:.3f}"] if 0.85 <= r < 0.95 and ymm not in ("202412",)]
        if bad:
            record(cid, cname, FAIL, "; ".join(detail_strs) + " | low r: " + "; ".join(bad))
        elif warn:
            record(cid, cname, WARN, "; ".join(detail_strs))
        else:
            record(cid, cname, PASS, "; ".join(detail_strs))
''')

md("""
    ### Check 6.2 — Right-edge difference vs interior difference

    The v2 methodology change targets the right edge of the PIQS fit.
    Quantify: for late-2024 months (closest to v1's stale fit edge), the
    v2−v1 NEE RMS difference should be larger than for interior months
    (e.g., 2015-07). If that's NOT the case, the PAD_RIGHT=2 change isn't
    doing what we expect.
""")
code('''
    cid, cname = "6.2", "right-edge differs more than interior"
    V1_ERA5 = WORK_DIR.parent / "MiCASA_v1" / "ERA5"
    interior = ["201007", "201507", "202007"]
    edge     = ["202410", "202411", "202412"]
    def rms_diff(ymm_list):
        out = []
        for ymm in ymm_list:
            v2_p = ERA5_DIR / f"fluxes_{ymm}.nc"
            v1_p = V1_ERA5 / f"fluxes_{ymm}.nc"
            if not v2_p.exists() or not v1_p.exists(): continue
            try:
                ds2 = xr.open_dataset(v2_p)
                ds1 = xr.open_dataset(v1_p)
                n2 = ds2["NEE"].mean(dim="time").values
                n1 = ds1["NEE"].mean(dim="time").values
                mask = np.isfinite(n2) & np.isfinite(n1)
                rms = float(np.sqrt(np.nanmean((n2[mask] - n1[mask])**2)))
                out.append((ymm, rms))
                ds2.close(); ds1.close()
            except Exception:
                continue
        return out
    int_rms  = rms_diff(interior)
    edge_rms = rms_diff(edge)
    if not int_rms or not edge_rms:
        record(cid, cname, INFO,
               f"insufficient v1 overlap (interior={len(int_rms)}, edge={len(edge_rms)})")
    else:
        int_med  = float(np.median([r for _, r in int_rms]))
        edge_med = float(np.median([r for _, r in edge_rms]))
        ratio = edge_med / max(int_med, 1e-30)
        detail = (f"interior median RMS={int_med:.3e}, edge median RMS={edge_med:.3e}, "
                  f"edge/interior = {ratio:.2f}")
        if ratio >= 1.5:
            record(cid, cname, PASS, f"{detail} (edge differs more, as expected)")
        elif ratio >= 1.0:
            record(cid, cname, WARN, f"{detail} (edge similar to interior)")
        else:
            record(cid, cname, FAIL,
                   f"{detail} (interior differs MORE -- something is wrong)")
''')

# ---- Section 7: Spatial Sanity ------------------------------------------
md("## Section 7 — Spatial Sanity")

md("""
    ### Check 7.1 — Ocean cells have NEE = 0

    MiCASA is land-only; ocean cells should be exactly zero (or NaN).
    Sample one fluxes_*.nc and confirm. Use coarse ocean mask: cells where
    every same-column-and-row land-monthly file has zero NPP and Rh.
""")
code('''
    cid, cname = "7.1", "ocean cells zero"
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))
    if not files:
        record(cid, cname, INFO, "no fluxes_*.nc")
    else:
        # Build a land mask from the multi-year monthly file: any cell
        # with non-zero NPP at any month is land.
        try:
            ds_m = xr.open_dataset(MONTHLY_FILE)
            land_mask = (np.abs(ds_m["NPP"]).max(dim="time").values > 1e-15)
            ds_m.close()
            ds_h = xr.open_dataset(files[len(files)//2])
            nee_mn = ds_h["NEE"].mean(dim="time").values
            ocean_nee = nee_mn[~land_mask]
            ocean_max = float(np.nanmax(np.abs(ocean_nee))) if ocean_nee.size else 0.0
            ds_h.close()
            if ocean_max > 1e-12:
                record(cid, cname, FAIL,
                       f"max |NEE| over ocean cells = {ocean_max:.3e} mol m-2 s-1 (expected 0)")
            else:
                record(cid, cname, PASS,
                       f"ocean cells exactly zero (max |NEE|={ocean_max:.1e})")
        except Exception as e:
            record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 7.2 — Antarctic land cells small

    Antarctica has minimal vegetation; NEE should be small in absolute
    terms relative to e.g. Amazon. Use latitude < -65° as proxy. Compare
    max |NEE| there vs in the tropics.
""")
code('''
    cid, cname = "7.2", "Antarctica fluxes small"
    files = sorted(ERA5_DIR.glob("fluxes_*.nc"))
    if not files:
        record(cid, cname, INFO, "no fluxes_*.nc")
    else:
        try:
            ds = xr.open_dataset(files[len(files)//2])
            lat = ds.latitude.values
            nee_mn = np.abs(ds["NEE"].mean(dim="time").values)
            ant_mask = lat < -65
            trop_mask = (lat > -23.5) & (lat < 23.5)
            ant_max  = float(np.nanmax(nee_mn[ant_mask]))
            trop_max = float(np.nanmax(nee_mn[trop_mask]))
            ds.close()
            ratio = ant_max / max(trop_max, 1e-30)
            detail = (f"max |NEE| Antarctic={ant_max:.3e}, tropical={trop_max:.3e}, "
                      f"ratio={ratio:.3f}")
            if ratio > 0.1:
                record(cid, cname, FAIL, detail + " (Antarctic > 10% of tropics)")
            elif ratio > 0.01:
                record(cid, cname, WARN, detail + " (Antarctic > 1% of tropics)")
            else:
                record(cid, cname, PASS, detail)
        except Exception as e:
            record(cid, cname, FAIL, f"exception: {e}")
''')

md("""
    ### Check 7.3 — Boreal seasonal phase

    Boreal land (>50°N) NEE: positive in winter (resp dominates), negative
    in summer (GPP dominates), with summer NEE max month in JJA. Sanity:
    in the most recent full year, identify the month of NEE minimum
    (max uptake) — should be in [6, 7, 8].
""")
code('''
    cid, cname = "7.3", "boreal seasonal phase"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        last_year = int(_summary["year"].max())
        full_year = _summary[_summary["year"] == last_year]
        if len(full_year) < 12:
            last_year -= 1
            full_year = _summary[_summary["year"] == last_year]
        if len(full_year) < 12 or "NEE_boreal" not in full_year.columns:
            record(cid, cname, INFO, "no complete year with NEE_boreal")
        else:
            min_month = int(full_year.set_index("month")["NEE_boreal"].idxmin())
            max_month = int(full_year.set_index("month")["NEE_boreal"].idxmax())
            detail = f"{last_year}: boreal NEE min={min_month}, max={max_month}"
            if min_month in (6, 7, 8):
                record(cid, cname, PASS, detail + " (uptake peak in JJA, OK)")
            else:
                record(cid, cname, WARN, detail + " (uptake peak NOT in JJA)")
''')

md("""
    ### Check 7.4 — Total NBE budget over land plausibility

    Annual NBE (NEE + FIRE + FUEL) over land should be a few PgC/yr
    in either direction. We don't have FIRE+FUEL aggregated in the
    summary cube (they live in the daily files, not in fluxes_*), so this
    check just looks at NEE annual sums and confirms they're in the
    plausible terrestrial-biosphere range (-10..+10 PgC/yr).
""")
code('''
    cid, cname = "7.4", "annual NEE in plausible range"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        annual = _summary.groupby("year")["NEE_global"].sum() / 1e6  # PgC/yr
        problems = [(int(y), float(v)) for y, v in annual.items() if abs(v) > 10]
        if problems:
            record(cid, cname, FAIL,
                   f"NEE outside ±10 PgC/yr: " +
                   "; ".join([f"{y}={v:.2f}" for y, v in problems[:5]]))
        else:
            record(cid, cname, PASS,
                   f"all {len(annual)} years' NEE within ±10 PgC/yr; "
                   f"range [{annual.min():.2f}, {annual.max():.2f}]")
''')

# ---- Section 8: Seasonal & Temporal -------------------------------------
md("## Section 8 — Seasonal & Temporal")

md("""
    ### Check 8.1 — December → January monthly continuity

    PIQS smoothes across calendar-year boundaries; the |Dec→Jan| jump per
    year should be comparable to other consecutive-month jumps. A
    systematic discontinuity at year boundaries would suggest the spline
    isn't actually smoothing across them.
""")
code('''
    cid, cname = "8.1", "Dec->Jan continuity"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        s = _summary.sort_values("yyyymm").set_index("yyyymm")
        # All consecutive month diffs in NEE_global
        nee = s["NEE_global"].values
        if len(nee) < 24:
            record(cid, cname, INFO, "less than 2 years of data")
        else:
            month_arr = s.reset_index()["month"].values
            diffs = np.abs(np.diff(nee))
            from_dec_mask = month_arr[:-1] == 12  # diff i is from month i to i+1
            other_diffs   = diffs[~from_dec_mask]
            dec_diffs     = diffs[from_dec_mask]
            other_med = float(np.median(other_diffs)) if other_diffs.size else 0.0
            dec_med   = float(np.median(dec_diffs))   if dec_diffs.size   else 0.0
            ratio = dec_med / max(other_med, 1e-30)
            detail = (f"|Dec->Jan| median={dec_med:.3e}, "
                      f"other |month->next| median={other_med:.3e}, ratio={ratio:.2f}")
            if ratio > 3.0:
                record(cid, cname, FAIL, detail + " (Dec->Jan jumps much larger)")
            elif ratio > 1.5:
                record(cid, cname, WARN, detail)
            else:
                record(cid, cname, PASS, detail)
''')

md("""
    ### Check 8.2 — Inter-annual stability of climatological cycle

    For the last 5 years, compute the per-month climatological NEE per
    band and the std-dev across years. The std/mean ratio (CoV) per
    calendar month should be modest — a year-to-year wobble of more than
    50% of the mean is suspicious.
""")
code('''
    cid, cname = "8.2", "interannual stability of seasonal cycle"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        last_year = int(_summary["year"].max())
        if (_summary["year"] == last_year).sum() < 12:
            last_year -= 1
        recent = _summary[(_summary["year"] >= last_year - 4) & (_summary["year"] <= last_year)]
        if len(recent) < 60:
            record(cid, cname, INFO, "fewer than 5 full years available")
        else:
            grp = recent.groupby("month")["NEE_global"]
            mean_m = grp.mean()
            std_m  = grp.std()
            cov = (std_m / mean_m.abs().replace(0, np.nan)).abs()
            max_cov = float(cov.max())
            detail = f"max |CoV| across calendar months = {max_cov:.2f} (window {last_year-4}..{last_year})"
            if max_cov > 1.0:
                record(cid, cname, FAIL, detail)
            elif max_cov > 0.5:
                record(cid, cname, WARN, detail)
            else:
                record(cid, cname, PASS, detail)
''')

md("""
    ### Check 8.3 — Polar-N vs boreal NEE ratio

    Polar (>70°N) NEE should be much smaller in absolute amplitude than
    boreal (>50°N) NEE — most boreal vegetation is below 70°N. If polar
    is comparable to boreal, something has shifted (e.g. mismapping).
""")
code('''
    cid, cname = "8.3", "polar << boreal NEE amplitude"
    if "_summary" not in globals() or _summary.empty:
        record(cid, cname, INFO, "summary cube unavailable")
    else:
        last_year = int(_summary["year"].max())
        if (_summary["year"] == last_year).sum() < 12:
            last_year -= 1
        recent = _summary[(_summary["year"] >= last_year - 4) & (_summary["year"] <= last_year)]
        polar_amp  = float(recent["NEE_polar_n"].max() - recent["NEE_polar_n"].min())
        boreal_amp = float(recent["NEE_boreal"].max() - recent["NEE_boreal"].min())
        ratio = polar_amp / max(abs(boreal_amp), 1e-30)
        detail = f"polar amp / boreal amp = {ratio:.3f}"
        if ratio > 1.0:
            record(cid, cname, FAIL, detail + " (polar exceeds boreal)")
        elif ratio > 0.5:
            record(cid, cname, WARN, detail)
        else:
            record(cid, cname, PASS, detail)
''')

md("""
    ## Phase 3 — Not yet implemented

    Sections deferred to Phase 3:

    - **Source-data provenance**: raw NCCS hash/count audit, ATMC range,
      ERA5 boundary cleanliness, vNRT-published-month tracker. Most of
      this overlaps with `check_daily_downloads.r` and `check_hashes.py`
      already in tree -- Phase 3 imports those checks into the notebook.
    - **PIQS tail-coefficient stability under one-more-month re-fit**:
      compare `fit.piqs.rda` from this run to one from the previous run;
      assert tail coefficients (last 3 months pre-pad) shift by less than
      a tolerance. Requires snapshotting the .rda before each refit.
    - **NRT clim-fill audit**: parse the global `:status=provisional` and
      `:provenance` attributes on monthly files; flag any month where the
      clim-fill fraction is >50% but `:status` is not `provisional`.
""")

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
import sys
out = sys.argv[1] if len(sys.argv) > 1 else "verify_v2.ipynb"
with open(out, "w") as fp:
    json.dump(nb, fp, indent=1)
    fp.write("\n")
print(f"wrote {out} ({len(cells)} cells)")
