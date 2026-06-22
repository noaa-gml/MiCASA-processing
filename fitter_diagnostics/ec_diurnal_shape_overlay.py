#!/usr/bin/env python3
"""Mean nighttime respiration DIURNAL SHAPE: observed (AmeriFlux) vs Q10(air) vs Q10(soil).

Distinct from the within-night anomaly test (ec_resp_driver_validation.py Test B, which
tied at R2~=0.003 because night-to-night weather noise swamps the per-night signal).
Here we average each site's nighttime NEE (~=RECO) BY HOUR-OF-DAY over the whole record,
which averages out the weather noise and reveals the CLIMATOLOGICAL diurnal shape. We
then ask whether that observed mean shape matches Q10(air temp) or Q10(soil temp) better.

Physics: air temp swings ~5-10 C overnight (steep cooling); 0-7 cm soil temp is damped &
lagged (gentle). If observed RECO declines gently / is near-flat overnight, soil's shape
fits; if it tracks the steep air cooling, air's shape fits. Q10 ref T cancels under the
unit-mean normalization, so only the SHAPE is compared. Night = SW<20, u*>0.2, NEE>0;
raw (non-gap-filled) sensors only; >=50 samples/hour and >=6 night hours per site.
Caveat: only the NIGHT arc of the diurnal cycle is observable (daytime RECO hidden by GPP).
"""
import numpy as np, pandas as pd, glob, os, re, warnings
warnings.simplefilter("ignore")
ROOT="/work2/noaa/co2/kaushik/drought/ameriflux_US"
USTAR_MIN=0.2; SWNIGHT=20.0; Q10=1.5; MIN_PER_HR=50; MIN_HRS=6

igbp={}
try:
    meta=pd.read_csv(f"{ROOT}/NAm_US_sites.csv",encoding="latin-1",engine="python",on_bad_lines="skip")
    for _,r in meta.iterrows():
        try: igbp[str(r["Site Id"]).strip()]=str(r["Vegetation Abbreviation (IGBP)"]).strip()
        except Exception: pass
except Exception: pass

def _parse(c,base):
    m=re.match(rf"^{base}(_PI)?(_F)?_(\d+)_(\d+)_(\d+|[A-Z])$",c)
    if not m: return None
    is_pi,is_f=bool(m.group(1)),bool(m.group(2))
    return (0 if (not is_pi and not is_f) else (1 if (is_pi and not is_f) else 2)), int(m.group(4))
def pick(hdr,base,shallowest=False,allow_gapfill=False):
    c=[]
    for col in hdr:
        if "QC" in col or "SSITC" in col: continue
        p=_parse(col,base)
        if p is None: continue
        rr,vert=p
        if rr==2 and not allow_gapfill: continue
        c.append(((vert if shallowest else 0,rr),col))
    if not c: return None
    c.sort(); return c[0][1]
def pick_nee(hdr):
    for b in ("NEE","FC"):
        x=pick(hdr,b,allow_gapfill=False)
        if x: return x
    return None
def q10f(Tc): return Q10**(Tc/10.0)            # ref cancels under unit-mean norm
def norm(x): return x/np.mean(x)

rows=[]; cyc=[]
for f in sorted(glob.glob(f"{ROOT}/EC_sitedata/*/AMF_*_BASE_HH_*.csv")):
    site=os.path.basename(f).split("_")[1]
    try:
        hdr=pd.read_csv(f,skiprows=2,nrows=0).columns.tolist()
        nee=pick_nee(hdr); ta=pick(hdr,"TA",allow_gapfill=True)
        ts=pick(hdr,"TS",shallowest=True); sw=pick(hdr,"SW_IN"); ust=pick(hdr,"USTAR")
        if not (nee and ta and ts and sw and ust): continue
        df=pd.read_csv(f,skiprows=2,usecols=["TIMESTAMP_START",nee,ta,ts,sw,ust],na_values=[-9999,"-9999"])
        df=df.rename(columns={nee:"NEE",ta:"TA",ts:"TS",sw:"SW",ust:"USTAR"}).dropna()
        t=pd.to_datetime(df.TIMESTAMP_START.astype("int64").astype(str),format="%Y%m%d%H%M")
        df["hour"]=t.dt.hour
        n=df[(df.SW<SWNIGHT)&(df.USTAR>USTAR_MIN)&(df.NEE>0)]
        if len(n)<300: continue
        g=n.groupby("hour"); cnt=g.size(); hrs=cnt[cnt>=MIN_PER_HR].index.values
        if len(hrs)<MIN_HRS: continue
        obs=g["NEE"].mean().loc[hrs].values
        air=q10f(g["TA"].mean().loc[hrs].values); soil=q10f(g["TS"].mean().loc[hrs].values)
        if np.nanstd(g["TA"].mean().loc[hrs].values)<0.3 and np.nanstd(g["TS"].mean().loc[hrs].values)<0.3: continue
        o,a,s=norm(obs),norm(air),norm(soil)
        ra=float(np.sqrt(np.mean((o-a)**2))); rs=float(np.sqrt(np.mean((o-s)**2)))
        rows.append(dict(site=site,igbp=igbp.get(site,"?"),nhours=len(hrs),npts=len(n),
                         rmse_air=ra,rmse_soil=rs,winner=("soil" if rs<ra else "air"),
                         obs_amp=float(o.max()-o.min()),air_amp=float(a.max()-a.min()),
                         soil_amp=float(s.max()-s.min())))
        for h,oo,aa,ss in zip(hrs,o,a,s): cyc.append(dict(site=site,igbp=igbp.get(site,"?"),hour=int(h),obs=oo,air=aa,soil=ss))
    except Exception: continue

R=pd.DataFrame(rows); C=pd.DataFrame(cyc)
R.to_csv("fitter_diagnostics/ec_diurnal_shape_sites.csv",index=False)
C.to_csv("fitter_diagnostics/ec_diurnal_shape_cycles.csv",index=False)
print(f"=== Mean nighttime respiration diurnal-shape fit ({len(R)} sites) ===")
nsoil=(R.winner=="soil").sum()
from scipy.stats import binomtest
print(f"soil shape fits better (lower RMSE): {nsoil}/{len(R)} (binom p={binomtest(nsoil,len(R),0.5).pvalue:.4f})")
print(f"median RMSE  air {R.rmse_air.median():.3f} | soil {R.rmse_soil.median():.3f} | median (air-soil) {(R.rmse_air-R.rmse_soil).median():+.3f}")
print(f"median normalized night amplitude  observed {R.obs_amp.median():.3f} | air-model {R.air_amp.median():.3f} | soil-model {R.soil_amp.median():.3f}")
print("by biome:")
for b,gg in R.groupby("igbp"):
    if len(gg)>=3: print(f"  {b:5s} n={len(gg):2d}: soil-wins {(gg.winner=='soil').sum()}/{len(gg)}  medianAmp obs {gg.obs_amp.median():.3f} air {gg.air_amp.median():.3f} soil {gg.soil_amp.median():.3f}")
print("\nsites (sorted by RMSE_air-RMSE_soil, + = soil better):")
print(R.assign(d=R.rmse_air-R.rmse_soil).sort_values("d",ascending=False)[["site","igbp","nhours","rmse_air","rmse_soil","winner","obs_amp","air_amp","soil_amp"]].to_string(index=False))
print("\ncycles -> ec_diurnal_shape_cycles.csv ; sites -> ec_diurnal_shape_sites.csv")
print("DONE")
