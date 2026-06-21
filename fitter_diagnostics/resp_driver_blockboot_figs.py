#!/usr/bin/env python3
"""Regenerate the 3 diurnalization figures on the full-year 2019 matched PCHIP
air-vs-soil pair, to match the spatial-block-bootstrap numbers in
resp_driver_blockboot.py (replaces the July-2020 single-month figures).

NOTE on amplitudes: the diurnal-cycle line plots are the GLOBAL-LAND-MEAN cycle
binned by UTC hour, which necessarily smears the local-solar-time peak across
longitudes — they are illustrative of shape (soil damps & lags air), not absolute
amplitude. All quantitative ratios in captions are PER-CELL amplitude ratios
(weighted median), computed cell-by-cell then aggregated — the correct statistic.
"""
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from dateutil.parser import parse as dparse

AIR="ERA5_bp_B_pchip_air"; SOIL="ERA5_bp_C_pchip_soil"
MONTHS=[f"2019{m:02d}" for m in range(1,13)]; NLAT,NLON=180,360

def hod_of(d):
    t=d.variables["time"]; v=np.asarray(t[:])
    base=dparse(t.units.split("since")[1].strip()).replace(tzinfo=None)
    sc={"days":86400,"day":86400,"hours":3600,"seconds":1}[t.units.split()[0]]
    return (((v*sc).astype("int64")//3600)%24).astype(int)

def cyc_amp(path,var):
    d=nc.Dataset(path); h=hod_of(d); x=np.asarray(d.variables[var][:])
    lat=np.asarray(d.variables["latitude"][:]); d.close()
    c=np.stack([x[h==k].mean(0) for k in range(24)])      # (24,lat,lon)
    return c, c.max(0)-c.min(0), x.mean(0), lat

def wmedian(x,w):
    o=np.argsort(x); cw=np.cumsum(w[o]); return x[o][np.searchsorted(cw,0.5*cw[-1])]

acc={"t2m":np.zeros((24,NLAT,NLON)),"stl1":np.zeros((24,NLAT,NLON)),
     "resp_air":np.zeros((24,NLAT,NLON)),"resp_soil":np.zeros((24,NLAT,NLON))}
resp_ratios=[]; resp_w=[]; forc_ratios=[]; forc_w=[]
for ym in MONTHS:
    ca,aa,mn,lat=cyc_amp(f"{AIR}/fluxes_{ym}.nc","resp")
    cs,as_,_,_  =cyc_amp(f"{SOIL}/fluxes_{ym}.nc","resp")
    ct,at,_,_   =cyc_amp(f"{SOIL}/fluxes_{ym}.nc","t2m")
    cl,al,_,_   =cyc_amp(f"{SOIL}/fluxes_{ym}.nc","stl1")
    acc["resp_air"]+=ca; acc["resp_soil"]+=cs; acc["t2m"]+=ct; acc["stl1"]+=cl
    land=np.isfinite(mn)&(np.abs(mn)>1e-9)
    W=np.clip(np.cos(np.deg2rad(lat)),0,None)[:,None]*np.ones((1,NLON))
    resp_ratios.append(as_[land]/np.clip(aa[land],1e-30,None)); resp_w.append(W[land])
    # per-cell forcing ratio: stl1 diurnal amplitude / t2m diurnal amplitude
    forc_ratios.append(al[land]/np.clip(at[land],1e-30,None)); forc_w.append(W[land])
    print("loaded",ym)
for k in acc: acc[k]/=12.0
W=np.clip(np.cos(np.deg2rad(lat)),0,None)[:,None]*np.ones((1,NLON))
def gmean(c):
    l=np.isfinite(c.mean(0)); w=W[l]; return np.array([(c[h][l]*w).sum()/w.sum() for h in range(24)])
gt=gmean(acc["t2m"]); gl=gmean(acc["stl1"]); gca=gmean(acc["resp_air"]); gcs=gmean(acc["resp_soil"])
resp_ratio=np.concatenate(resp_ratios); resp_ww=np.concatenate(resp_w)
forc_ratio=np.concatenate(forc_ratios); forc_ww=np.concatenate(forc_w)
med=wmedian(resp_ratio,resp_ww)
forc_med=wmedian(forc_ratio,forc_ww)              # CORRECT per-cell forcing ratio
print(f"\nper-cell forcing (stl1/t2m) amplitude ratio, annual 2019 weighted median = {forc_med:.4f}")
print(f"per-cell resp amplitude ratio weighted median = {med:.4f}")

plt.rcParams.update({"font.size":11,"figure.dpi":130})
# Fig 1: forcing diurnal cycle (shape illustration; caption uses per-cell ratio)
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.plot(range(24),gt-gt.mean(),"o-",color="#c1121f",label="2-m air temp t2m")
ax.plot(range(24),gl-gl.mean(),"s--",color="#0353a4",label="0-7cm soil temp stl1")
ax.set_xlabel("hour of day (UTC)"); ax.set_ylabel("temperature anomaly (K)")
ax.set_title(f"ERA5 forcing: 0–7cm soil temp damps & lags 2-m air\nglobal-land mean cycle (UTC-smeared); per-cell amplitude ratio {forc_med:.2f}, full-year 2019")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_forcing_t2m_vs_stl1.png"); plt.close(fig)
# Fig 2: resp diurnal anomaly air vs soil (shape illustration)
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.plot(range(24),(gca-gca.mean())*1e9,"o-",color="#c1121f",label="air temp (legacy)")
ax.plot(range(24),(gcs-gcs.mean())*1e9,"s--",color="#0353a4",label="soil temp (prototype #1)")
ax.set_xlabel("hour of day (UTC)"); ax.set_ylabel("resp anomaly (nmol m$^{-2}$ s$^{-1}$)")
ax.set_title("Respiration diurnal cycle: soil driver damps & lags (shape, UTC-smeared)\nglobal-land mean, full-year 2019")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_diurnal_air_vs_soil.png"); plt.close(fig)
# Fig 3: per-cell amplitude-ratio histogram (resp), annual, block-bootstrap CI
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.hist(np.clip(resp_ratio,0,2),bins=60,weights=resp_ww,color="#0353a4",alpha=.8)
ax.axvline(med,color="k",lw=2,label=f"area-wtd median {med:.2f}  (block-10° 95%CI [0.78,0.83])")
ax.axvline(1.0,color="#c1121f",ls=":",lw=2,label="1.0 (no change)")
ax.set_xlabel("per-cell respiration amplitude ratio soil/air"); ax.set_ylabel("area-weighted land cells")
ax.set_title("Soil-temp respiration diurnal amplitude < air-temp\nfull-year 2019 (ratio < 1 = damped); 15.6k land cells × 12 months pooled")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_amplitude_ratio_hist.png"); plt.close(fig)
print("Wrote 3 figures")
