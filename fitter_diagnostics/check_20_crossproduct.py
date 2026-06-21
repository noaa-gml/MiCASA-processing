#!/usr/bin/env python3
"""Section-20 cross-product comparison — the real computation behind verify_v2
Checks 20.1 and 20.2 (previously deferred stubs).

20.1  v2-vs-v1 lat-band annual NEE: do the V2 changes (corrected 0.1->1deg
      aggregation §3.1, polar clip) shift mass between latitude bands? Compare
      per-band area-weighted annual mean NEE (= Rh - NPP) from the V1 monthly
      product vs the V2-pipeline monthly product over the overlap 2001..2024.
      Diurnalize preserves monthly means, so the monthly-product per-band annual
      NEE equals the shipped (diurnalized) per-band annual NEE to the polar-clip
      residual — so this is the cheap, exact proxy the original check deferred.

20.2  Global carbon-budget context: compute MiCASA global annual NBE
      (= Rh - NPP + FIRE + FUEL) in PgC/yr and place it against the Global Carbon
      Budget land-sink range (Friedlingstein et al. 2023). CASA-only (no ATMC) is
      NOT expected to close the growth-rate budget — the offset is the ATMC-type
      term the inversion supplies (§5.2) — so this is reported as context (INFO).
"""
import numpy as np, netCDF4 as nc, glob, os

V1 = "../MiCASA_v1/monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc"
V2GLOB = "monthly_1x1/MiCASA_v1_flux_x360_y180_monthly_2*.nc"
R = 6.371e6; D2R = np.pi/180.0; SPY = 365.25*86400.0

def band_masks(lat):
    return {"nh_mid":(lat>=30)&(lat<60), "trop":(lat>=-30)&(lat<30),
            "sh_mid":(lat>=-60)&(lat<-30), "boreal":(lat>=60)}

def cell_area(lat, nlon):
    # 1deg grid: area(phi) = R^2 * dlon * (sin(phi+dphi/2)-sin(phi-dphi/2))
    dphi = D2R*1.0; dlon = D2R*1.0
    a = R*R*dlon*(np.sin(lat*D2R+dphi/2)-np.sin(lat*D2R-dphi/2))  # per cell, per lon
    return np.repeat(a[:,None], nlon, axis=1)                     # (lat,lon) m^2

# ---- 20.1: per-band annual NEE, V1 vs V2 -----------------------------------
d1 = nc.Dataset(V1); lat = d1.variables["latitude"][:].data
nlon = d1.dimensions["longitude"].size
w = np.cos(np.radians(lat))
masks = band_masks(lat)
# v1 time -> years
tv = d1.variables["time"]; import datetime as dt
from dateutil.parser import parse as dparse
base = dparse(tv.units.split("since")[1].strip()).replace(tzinfo=None)
scale = {"days":1,"day":1,"hours":1/24,"seconds":1/86400}[tv.units.split()[0]]
yrs1 = np.array([(base+dt.timedelta(days=float(v)*scale)).year for v in tv[:]])
nee1 = (-d1.variables["NPP"][:].data + d1.variables["Rh"][:].data)   # (t,lat,lon)
d1.close()

def perband_annual_from_cube(nee, yrs):
    out = {}
    for bk,m in masks.items():
        sub = nee[:, m, :].mean(axis=2)                 # (t, lat-in-band)
        wt  = (sub * w[m]).sum(axis=1)/w[m].sum()        # (t,)
        ann = {}
        for y in np.unique(yrs):
            sel = yrs==y
            if sel.any(): ann[y] = wt[sel].mean()
        out[bk] = ann
    return out

v1b = perband_annual_from_cube(nee1, yrs1)

# v2 from per-month files
v2_files = sorted(glob.glob(V2GLOB))
v2_acc = {bk:{} for bk in masks}
for vf in v2_files:
    ym = os.path.basename(vf)[-9:-3]; yr=int(ym[:4])
    if yr<2001 or yr>2024: continue
    d = nc.Dataset(vf)
    nee = np.squeeze(-d.variables["NPP"][:].data + d.variables["Rh"][:].data)  # (lat,lon)
    d.close()
    for bk,m in masks.items():
        val = ((nee[m,:].mean(axis=1))*w[m]).sum()/w[m].sum()
        v2_acc[bk].setdefault(yr,[]).append(val)
v2b = {bk:{y:np.mean(v) for y,v in d.items()} for bk,d in v2_acc.items()}

print("=== 20.1  v2-vs-v1 lat-band annual NEE (gC m-2 s-1, area-wt mean, 2001-2024) ===")
worst=0.0
for bk in masks:
    ys=[y for y in range(2001,2025) if y in v1b[bk] and y in v2b[bk]]
    m1=np.mean([v1b[bk][y] for y in ys]); m2=np.mean([v2b[bk][y] for y in ys])
    rel=abs(m2-m1)/max(abs(m1),1e-30)*100; worst=max(worst,rel)
    print(f"  {bk:7s} v1={m1:+.3e}  v2={m2:+.3e}  rel diff {rel:.2f}%")
print(f"  => max per-band rel diff = {worst:.2f}%  ({'PASS <5%' if worst<5 else 'WARN' if worst<15 else 'FAIL'})")

# ---- 20.2: global annual NBE budget context --------------------------------
area = cell_area(lat, nlon)                   # (lat,lon) m^2
# recompute global NBE per year from v2 monthly files (Rh-NPP+FIRE+FUEL)
nbe_yr = {}
for vf in v2_files:
    yr=int(os.path.basename(vf)[-9:-3][:4])
    if yr<2001 or yr>2024: continue
    d=nc.Dataset(vf)
    f = np.squeeze(-d.variables["NPP"][:].data + d.variables["Rh"][:].data
                   + d.variables["FIRE"][:].data + d.variables["FUEL"][:].data)  # gC/m2/s
    d.close()
    pgC = np.nansum(f*area)*SPY/1e15          # gC/s -> PgC/yr
    nbe_yr.setdefault(yr,[]).append(pgC)
nbe = {y:np.mean(v) for y,v in nbe_yr.items()}
vals=np.array([nbe[y] for y in sorted(nbe)])
print("\n=== 20.2  MiCASA global annual NBE (= Rh-NPP+FIRE+FUEL), PgC/yr ===")
print(f"  2001-2024 mean = {vals.mean():+.2f} PgC/yr ; range [{vals.min():+.2f}, {vals.max():+.2f}]")
print( "  GCB2023 net land sink (incl. fire, excl. LUC): ~ -1.9 to -3.5 PgC/yr decadal")
print(f"  => CASA-only NBE offset vs observed land sink ~ {abs(vals.mean()-(-2.6)):.1f} PgC/yr")
print( "     (of order the ~3 PgC/yr ATMC term the inversion supplies — §5.2; INFO, not a closure)")
print("\nDONE")
