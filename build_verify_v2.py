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

# ---- Phase 2 stubs ------------------------------------------------------
md("""
    ## Phase 2 / 3 — Not yet implemented

    Sections to add in subsequent phases:

    - **Section 5: Comparison vs prior products** — spatial correlation of v2
      `ERA5/fluxes_*.nc` against v1's `ERA5/fluxes_*.nc` (where they overlap),
      time series of global NEE/GPP/Rh, COVID-2020 anomaly visibility.
    - **Section 6: Spatial sanity** — Antarctica fluxes ~0; Amazon, boreal,
      monsoon hotspots seasonal pattern; Antarctic mask consistency.
    - **Section 7: Seasonal sanity** — NH/SH out-of-phase, latitude
      amplitude profile, year-to-year stability.
    - **Section 8: Budget identities** — NEE = -NPP + Rh - ATMC where ATMC
      flows through; NBE = NEE + FIRE + FUEL.
    - **Section 9: NRT-specific** — clim-fill fraction per provisional month;
      vNRT->v1 symlink integrity per day; PIQS tail-coefficient stability
      under one-more-month re-fit (compare two consecutive runs).

    Phase 1 covers structural invariants and schema. Phase 2 covers
    "looks reasonable". Phase 3 covers source-data provenance.
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
