#!/usr/bin/env python3
"""Direct PIQS vs PCHIP budget diff — the measurement behind the §0 'budget
invariance' claim that was previously argued-not-diffed.

Computes global annual NEE (PgC/yr) from the two FULL diurnalized products:
  PCHIP (shipped V2) : MiCASA_v2/ERA5/fluxes_YYYYMM.nc
  PIQS  (V1-style)   : MiCASA_v1_piqs/ERA5/fluxes_YYYYMM.nc
then the long-term trend, the 2015-16 El Nino anomaly, and the 2020 COVID anomaly,
for BOTH, and the PIQS-PCHIP difference. If the fitter switch cannot move the
science signal, the differences are ~0 (to the polar-clip residual / FP).
"""
import numpy as np, netCDF4 as nc, calendar, glob, os

PROD = {"PCHIP": "ERA5", "PIQS": "../MiCASA_v1_piqs/ERA5"}
YEARS = list(range(2001, 2025))            # complete years, matches §20 window
R = 6.371e6; D2R = np.pi/180.0; GC_PER_MOL = 12.011

# cell area (1deg), from any file's latitude
d0 = nc.Dataset(f"{PROD['PCHIP']}/fluxes_200101.nc"); lat = d0.variables["latitude"][:]
nlon = d0.dimensions["longitude"].size; d0.close()
acell = R*R*D2R*(np.sin(lat*D2R+D2R/2)-np.sin(lat*D2R-D2R/2))
area = np.repeat(acell[:, None], nlon, axis=1)        # (lat,lon) m^2

def annual_series(era5dir):
    out = {}
    for y in YEARS:
        tot = 0.0
        for m in range(1, 13):
            f = f"{era5dir}/fluxes_{y}{m:02d}.nc"
            if not os.path.exists(f): continue
            ds = nc.Dataset(f)
            nee = ds.variables["NEE"][:]                # (time,lat,lon) mol m-2 s-1
            ds.close()
            secs = calendar.monthrange(y, m)[1] * 86400.0
            mean_flux = np.nanmean(nee, axis=0)         # (lat,lon) time-mean
            tot += np.nansum(mean_flux * area) * secs * GC_PER_MOL / 1e15  # PgC
        out[y] = tot
        print(f"  {os.path.basename(os.path.dirname(era5dir+'/x'))[:18]:18s} {y}: {tot:+.4f} PgC/yr", flush=True)
    return out

def trend(series):
    ys = np.array(sorted(series)); v = np.array([series[y] for y in ys])
    return np.polyfit(ys, v, 1)[0], v.mean()

print("=== global annual NEE (PgC/yr) ===")
S = {}
for name, d in PROD.items():
    print(f"-- {name} ({d}) --")
    S[name] = annual_series(d)

print("\n=== PIQS vs PCHIP (the budget-invariance measurement) ===")
for name in PROD:
    tr, mn = trend(S[name])
    e = np.mean([S[name][y] for y in (2015, 2016)]); c = S[name][2020]
    base = np.mean([S[name][y] for y in (2010,2011,2012,2013,2014,2017,2018,2019)])
    print(f"{name:6s}: mean {mn:+.4f} | trend {tr:+.5f} PgC/yr/yr | 2015-16 anom {e-base:+.4f} | 2020 anom {c-base:+.4f}")
trP, mnP = trend(S["PCHIP"]); trQ, mnQ = trend(S["PIQS"])
maxdiff = max(abs(S["PIQS"][y]-S["PCHIP"][y]) for y in YEARS)
maxrel  = max(abs(S["PIQS"][y]-S["PCHIP"][y])/max(abs(S["PCHIP"][y]),1e-9) for y in YEARS)
print(f"\nDelta (PIQS - PCHIP): trend {trQ-trP:+.2e} PgC/yr/yr ; mean {mnQ-mnP:+.2e} PgC/yr")
print(f"max |annual diff| over 2001-2024: {maxdiff:.2e} PgC/yr ({100*maxrel:.4f}% of annual)")
print("=> fitter switch moves the global annual budget by this much (≈0 confirms §0).")
print("DONE")
