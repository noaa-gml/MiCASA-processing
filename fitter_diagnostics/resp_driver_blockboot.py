#!/usr/bin/env python3
"""Spatial block-bootstrap CIs for the soil-temp respiration driver.

Replaces the i.i.d.-cell resample in resp_driver_ci_plots.py (which mislabeled
itself a 'block bootstrap' and, by treating ~15k spatially autocorrelated land
cells as independent, produced CIs far too tight). Here the resampling unit is a
BxB-degree spatial block, so the CI reflects the field's spatial decorrelation.

Data: full-year 2019 matched PCHIP pair (air vs soil respiration driver),
  air  = ERA5_bp_B_pchip_air   (respiration_temperature_driver = airtemp)
  soil = ERA5_bp_C_pchip_soil  (respiration_temperature_driver = soiltemp)
All 12 months -> addresses both the CI-method bug and the 'July+Jan only' /
under-powering critique.

For each month: per-cell diurnal amplitude (max-min of the 24-hour climatology),
ratio soil/air for resp and NEE on land cells. Then:
  - naive i.i.d.-cell bootstrap   (the OLD method, for contrast)
  - spatial block bootstrap        (resample BxB-deg blocks with replacement)
at block sizes 5/10/20 deg, global + latitude bands + pooled-annual.
"""
import numpy as np, netCDF4 as nc, datetime as dt
from dateutil.parser import parse as dparse

RNG = np.random.default_rng(7)
AIR  = "ERA5_bp_B_pchip_air"
SOIL = "ERA5_bp_C_pchip_soil"
MONTHS = [f"2019{m:02d}" for m in range(1, 13)]
NLAT, NLON = 180, 360

def hod_of(d):
    t = d.variables["time"]; vals = np.asarray(t[:])
    base = dparse(t.units.split("since")[1].strip()).replace(tzinfo=None)
    scale = {"days":86400,"day":86400,"hours":3600,"seconds":1}[t.units.split()[0]]
    secs = (vals*scale).astype("int64")
    return ((secs//3600) % 24).astype(int)

def amp(path, var):
    """Per-cell diurnal amplitude (max-min over hour-of-day climatology)."""
    d = nc.Dataset(path); hod = hod_of(d)
    x = np.asarray(d.variables[var][:]); lat = np.asarray(d.variables["latitude"][:])
    d.close()
    cyc = np.stack([x[hod==h].mean(0) for h in range(24)])   # (24,lat,lon)
    return cyc.max(0)-cyc.min(0), x.mean(0), lat

def wmedian(x, w):
    o = np.argsort(x); cw = np.cumsum(w[o])
    return x[o][np.searchsorted(cw, 0.5*cw[-1])]

# ---- accumulate small per-month ratio fields (free the hourly arrays) -------
lat = np.arange(NLAT)  # placeholder; set from data
fields = {"resp": [], "NEE": []}
mnmask = []
for ym in MONTHS:
    rec = {}
    for comp in ["resp", "NEE"]:
        a_air, mn_air, lat = amp(f"{AIR}/fluxes_{ym}.nc", comp)
        a_soil, _, _       = amp(f"{SOIL}/fluxes_{ym}.nc", comp)
        ratio = a_soil/np.clip(a_air, 1e-30, None)
        rec[comp] = (ratio, mn_air)
    # land mask from resp monthly mean
    land = np.isfinite(rec["resp"][1]) & (np.abs(rec["resp"][1]) > 1e-9)
    mnmask.append(land)
    fields["resp"].append(rec["resp"][0])
    fields["NEE"].append(rec["NEE"][0])
    print(f"loaded {ym}: {land.sum()} land cells")

W2 = np.clip(np.cos(np.deg2rad(lat)), 0, None)[:, None]*np.ones((1, NLON))
LATG = lat[:, None]*np.ones((1, NLON))

def block_ids(bs):
    """Integer block id per (lat,lon) cell for a bs-degree tiling."""
    li = (np.arange(NLAT)//bs)[:, None]*np.ones((1, NLON), int)
    lj = (np.arange(NLON)//bs)[None, :]*np.ones((NLAT, 1), int)
    return li*1000 + lj

def boot_iid(r, w, B=2000):
    n = r.size; out = np.empty(B)
    for b in range(B):
        s = RNG.integers(0, n, n); out[b] = wmedian(r[s], w[s])
    return np.percentile(out, [2.5, 97.5])

def boot_block(r, w, blk, B=2000):
    """Resample whole blocks with replacement; pool their cells."""
    uniq = np.unique(blk)
    members = {u: np.where(blk == u)[0] for u in uniq}
    nb = uniq.size; out = np.empty(B)
    for b in range(B):
        pick = uniq[RNG.integers(0, nb, nb)]
        idx = np.concatenate([members[u] for u in pick])
        out[b] = wmedian(r[idx], w[idx])
    return np.percentile(out, [2.5, 97.5])

def pool(comp, latrange=None):
    """Pool all 12 months: flat arrays of ratio, weight, blockid(per bs), lat."""
    rs, ws, ls = [], [], []
    blkmaps = {bs: [] for bs in (5, 10, 20)}
    for k in range(12):
        land = mnmask[k].copy()
        if latrange is not None:
            a, b = latrange; land &= (LATG >= a) & (LATG <= b)
        rs.append(fields[comp][k][land]); ws.append(W2[land]); ls.append(LATG[land])
        for bs in (5, 10, 20):
            blkmaps[bs].append(block_ids(bs)[land])
    r = np.concatenate(rs); w = np.concatenate(ws)
    blk = {bs: np.concatenate(blkmaps[bs]) for bs in (5, 10, 20)}
    return r, w, blk

print("\n================ SPATIAL BLOCK BOOTSTRAP (full-year 2019) ================")
for comp in ["resp", "NEE"]:
    r, w, blk = pool(comp)
    med = wmedian(r, w)
    iid = boot_iid(r, w)
    print(f"\n[{comp}] area-weighted median soil/air amplitude ratio = {med:.4f}")
    print(f"   naive i.i.d.-cell  95% CI [{iid[0]:.4f}, {iid[1]:.4f}]  (width {iid[1]-iid[0]:.4f})  <- OLD METHOD")
    for bs in (5, 10, 20):
        ci = boot_block(r, w, blk[bs])
        print(f"   block {bs:2d}deg       95% CI [{ci[0]:.4f}, {ci[1]:.4f}]  (width {ci[1]-ci[0]:.4f})  blocks~{np.unique(blk[bs]).size}")

print("\n---- latitude bands (resp, block 10deg, annual) ----")
BANDS = {"boreal 50-70N": (50, 70), "NH-temp 25-50N": (25, 50),
         "tropics 25S-25N": (-25, 25), "SH-temp 50-25S": (-50, -25)}
for nm, lr in BANDS.items():
    r, w, blk = pool("resp", lr)
    if r.size < 20: continue
    med = wmedian(r, w); ci = boot_block(r, w, blk[10])
    print(f"   {nm:16s} {med:.4f}  95% CI [{ci[0]:.4f}, {ci[1]:.4f}]")

print("\n---- per-month NEE ratio (block 10deg) : seasonal coverage ----")
for k, ym in enumerate(MONTHS):
    land = mnmask[k]
    r = fields["NEE"][k][land]; w = W2[land]; blk = block_ids(10)[land]
    med = wmedian(r, w); ci = boot_block(r, w, blk, B=800)
    print(f"   {ym}: NEE ratio {med:.4f}  95% CI [{ci[0]:.4f}, {ci[1]:.4f}]")
print("\nDONE")
