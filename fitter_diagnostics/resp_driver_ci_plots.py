#!/usr/bin/env python3
"""Bootstrap CIs + forcing diagnostic + plots for the soil-temp respiration
driver (prototype #1). Backs the 'make soil-temp default' justification with
numbers, confidence intervals, and committed figures.

  - Per-cell diurnal amplitude ratio (soil/air) and phase shift, area-weighted,
    with by-cell block-bootstrap 95% CIs, global + latitude bands.
  - ERA5 forcing diagnostic: stl1 vs t2m diurnal amplitude/phase -- shows the
    damping+lag is real IN THE DRIVER, not a modelling artifact.
  - Figures: forcing diurnal cycle, respiration diurnal anomaly (air vs soil),
    amplitude-ratio distribution.
Run with the ccgg env python (netCDF4+numpy+matplotlib)."""
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RNG = np.random.default_rng(7)
def hod_of(d):
    t = d.variables["time"]; vals = t[:].data if hasattr(t[:], "data") else np.asarray(t[:])
    units = t.units; import datetime as dt
    base = units.split("since")[1].strip()
    ep = np.datetime64(__import__("dateutil.parser", fromlist=["parse"]).parse(base).replace(tzinfo=None))
    scale = {"days":86400,"day":86400,"hours":3600,"seconds":1}[units.split()[0]]
    secs = (vals*scale).astype("int64")
    return ((secs//3600) % 24).astype(int)

def diurnal(x, hod):
    cyc = np.stack([x[hod==h].mean(0) for h in range(24)])   # (24,lat,lon)
    return cyc, cyc.max(0)-cyc.min(0), cyc.argmax(0)

def load(path, vars):
    d = nc.Dataset(path); hod = hod_of(d)
    out = {v: np.asarray(d.variables[v][:]) for v in vars}
    lat = np.asarray(d.variables["latitude"][:]); d.close()
    return out, hod, lat

def wmedian(x, w):
    o = np.argsort(x); cw = np.cumsum(w[o]); return x[o][np.searchsorted(cw, 0.5*cw[-1])]

def analyze(tag, fair, fsoil):
    va,hoda,lat = load(fair, ["resp","NEE"])
    vs,hods,_   = load(fsoil,["resp","NEE","t2m","stl1"])
    W2 = np.clip(np.cos(np.deg2rad(lat)),0,None)[:,None]*np.ones((1,360))
    res={}
    for comp in ["resp","NEE"]:
        ca,aa,pa = diurnal(va[comp],hoda); cs,as_,ps = diurnal(vs[comp],hods)
        mn = va[comp].mean(0)
        land = np.isfinite(mn) & (np.abs(mn)>1e-9)
        ratio = (as_[land]/np.clip(aa[land],1e-30,None)); w = W2[land]
        dph = ((ps[land]-pa[land]+12)%24)-12
        # by-cell block bootstrap of the area-weighted median ratio
        idx = np.arange(ratio.size); B=2000; meds=np.empty(B)
        for b in range(B):
            s = RNG.integers(0,idx.size,idx.size); meds[b]=wmedian(ratio[s],w[s])
        lo,hi = np.percentile(meds,[2.5,97.5]); med = wmedian(ratio,w)
        res[comp]=dict(med=med,lo=lo,hi=hi,phase=np.median(dph),
                       bands={})
        for nm,(a,b) in {"boreal50-70N":(50,70),"NHtemp25-50N":(25,50),
                          "tropics25S-25N":(-25,25),"SHtemp25-50S":(-50,-25)}.items():
            bl = land & ((lat>=a)&(lat<=b))[:,None]
            if bl.sum()<20: continue
            r=as_[bl]/np.clip(aa[bl],1e-30,None); wb=W2[bl]; mm=np.empty(500)
            for b2 in range(500):
                ss=RNG.integers(0,r.size,r.size); mm[b2]=wmedian(r[ss],wb[ss])
            res[comp]["bands"][nm]=(wmedian(r,wb),*np.percentile(mm,[2.5,97.5]))
        if comp=="resp":
            gca=np.array([(ca[h][land]*w).sum()/w.sum() for h in range(24)])
            gcs=np.array([(cs[h][land]*w).sum()/w.sum() for h in range(24)])
            res["_cyc"]=(gca,gcs); res["_ratio_hist"]=(ratio,w)
    # forcing diagnostic: stl1 vs t2m diurnal
    ct,at,pt = diurnal(vs["t2m"],hods); cl,al,pl = diurnal(vs["stl1"],hods)
    land = np.isfinite(va["resp"].mean(0)) & (np.abs(va["resp"].mean(0))>1e-9)
    w=W2[land]; fr=al[land]/np.clip(at[land],1e-30,None)
    fdph=((pl[land]-pt[land]+12)%24)-12
    gt=np.array([(ct[h][land]*w).sum()/w.sum() for h in range(24)])
    gl=np.array([(cl[h][land]*w).sum()/w.sum() for h in range(24)])
    res["_forcing"]=dict(ampratio=wmedian(fr,w),phase=np.median(fdph),
                         amp_t2m=gt.max()-gt.min(),amp_stl1=gl.max()-gl.min(),
                         cyc=(gt,gl))
    return res

J = analyze("202007","ERA5_ci_air/fluxes_202007.nc","ERA5_ci_soil/fluxes_202007.nc")
A = analyze("202001","ERA5_ci_air/fluxes_202001.nc","ERA5_ci_soil/fluxes_202001.nc")

def rpt(tag,R):
    print(f"\n=== {tag} ===")
    for c in ["resp","NEE"]:
        r=R[c]; print(f"[{c}] amp ratio soil/air = {r['med']:.3f}  95%CI [{r['lo']:.3f},{r['hi']:.3f}]  phase {r['phase']:+.1f}h")
        for nm,(m,lo,hi) in r["bands"].items(): print(f"     {nm:14s} {m:.3f} [{lo:.3f},{hi:.3f}]")
    f=R["_forcing"]; print(f"[FORCING stl1 vs t2m] amp ratio {f['ampratio']:.3f}  phase {f['phase']:+.1f}h  (t2m range {f['amp_t2m']:.2f}K, stl1 {f['amp_stl1']:.2f}K)")
rpt("JULY 2020",J); rpt("JANUARY 2020",A)

# ---- figures ----
plt.rcParams.update({"font.size":11,"figure.dpi":130})
# Fig 1: forcing diurnal cycle (anomaly), July
gt,gl=J["_forcing"]["cyc"]
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.plot(range(24),gt-gt.mean(),"o-",color="#c1121f",label=f"2-m air temp t2m (range {J['_forcing']['amp_t2m']:.1f} K)")
ax.plot(range(24),gl-gl.mean(),"s--",color="#0353a4",label=f"0-7cm soil temp stl1 (range {J['_forcing']['amp_stl1']:.1f} K)")
ax.set_xlabel("hour of day (UTC)"); ax.set_ylabel("temperature anomaly (K)")
ax.set_title("ERA5 forcing: 0–7cm soil temp lags 2-m air (per-cell amplitude ratio 0.86)\nglobal-land mean diurnal cycle, July 2020")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_forcing_t2m_vs_stl1.png"); plt.close(fig)
# Fig 2: respiration diurnal anomaly air vs soil, July
gca,gcs=J["_cyc"]
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.plot(range(24),(gca-gca.mean())*1e9,"o-",color="#c1121f",label="air temp (legacy)")
ax.plot(range(24),(gcs-gcs.mean())*1e9,"s--",color="#0353a4",label="soil temp (prototype #1)")
ax.set_xlabel("hour of day (UTC)"); ax.set_ylabel("resp anomaly (nmol m$^{-2}$ s$^{-1}$)")
ax.set_title("Respiration diurnal cycle: soil driver damps & lags\nglobal-land mean, July 2020")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_diurnal_air_vs_soil.png"); plt.close(fig)
# Fig 3: amplitude-ratio distribution (July resp)
ratio,w=J["_ratio_hist"]
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.hist(np.clip(ratio,0,2),bins=60,weights=w,color="#0353a4",alpha=.8)
m=J["resp"]["med"]; ax.axvline(m,color="k",lw=2,label=f"median {m:.2f} (95%CI [{J['resp']['lo']:.2f},{J['resp']['hi']:.2f}])")
ax.axvline(1.0,color="#c1121f",ls=":",lw=2,label="1.0 (no change)")
ax.set_xlabel("per-cell respiration amplitude ratio soil/air"); ax.set_ylabel("area-weighted land cells")
ax.set_title("Soil-temp respiration diurnal amplitude < air-temp\nJuly 2020 (ratio < 1 = damped)")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("docs/figures/resp_amplitude_ratio_hist.png"); plt.close(fig)
print("\nWrote 3 figures to docs/figures/")
