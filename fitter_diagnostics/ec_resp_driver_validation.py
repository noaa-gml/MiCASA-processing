#!/usr/bin/env python3
"""Eddy-covariance validation of the respiration temperature driver (soil vs air),
deepened. Open gate for the soil-temp diurnalization flip (§2).

At night SW<20 -> NEE ≈ ecosystem respiration (no GPP, no partitioning model). For
each AmeriFlux site we test whether nighttime respiration is better explained by 2-m
air temperature or by SHALLOW soil temperature, three ways:
  1. R^2 of ln(resp) vs each temperature separately.
  2. COMPETITIVE: multiple regression ln(resp) ~ z(TA) + z(TS); the larger |standardized
     beta| is the temperature that wins when both compete (decisive, controls for the
     air-soil correlation).
  3. Q10 implied by each single fit (Q10 = exp(10*slope)); physical range ~1.5-3.
Stratified by IGBP biome. Prefers RAW (non-gap-filled) flux and the SHALLOWEST soil
temp (closest to ERA5 0-7 cm stl1).
"""
import numpy as np, pandas as pd, glob, os, re, warnings
warnings.simplefilter("ignore")
ROOT = "/work2/noaa/co2/kaushik/drought/ameriflux_US"

# IGBP biome per site
igbp = {}
try:
    meta = pd.read_csv(f"{ROOT}/NAm_US_sites.csv", encoding="latin-1", engine="python", on_bad_lines="skip")
    for _, r in meta.iterrows():
        try: igbp[str(r["Site Id"]).strip()] = str(r["Vegetation Abbreviation (IGBP)"]).strip()
        except Exception: pass
except Exception as e:
    print("biome metadata unavailable:", e)

def pick_nee(hdr):
    # prefer RAW turbulent NEE/FC over gap-filled (_F) to avoid temperature-based fill
    for pref in ("NEE_PI_1_1_1","FC_1_1_1","NEE_1_1_1","NEE_PI_F_1_1_1","FC_PI_F_1_1_1"):
        if pref in hdr: return pref
    for c in hdr:
        if re.match(r"^(NEE|FC)_(PI_)?\d", c) and "QC" not in c and "_F_" not in c: return c
    for c in hdr:
        if re.match(r"^(NEE|FC)_", c) and "QC" not in c: return c
    return None
def pick_ta(hdr):
    for pref in ("TA_1_1_1","TA_PI_F_1_1_1"):
        if pref in hdr: return pref
    return next((c for c in hdr if c.startswith("TA_") and "QC" not in c), None)
def pick_ts(hdr):
    # shallowest soil temp: TS[_PI_F]_H_V_R -> smallest V (vertical=depth index)
    cands=[]
    for c in hdr:
        m=re.match(r"^TS(_PI_F)?_(\d+)_(\d+)_(\d+)$", c)
        if m and "QC" not in c: cands.append((int(m.group(3)), 0 if m.group(1) else 1, c))  # (depth, prefer non-PI? no: prefer raw shallow)
    if not cands: return None
    cands.sort()  # shallowest first
    return cands[0][2]
def pick_sw(hdr): return next((c for c in hdr if c.startswith("SW_IN")), None)

rows=[]
for f in sorted(glob.glob(f"{ROOT}/EC_sitedata/*/AMF_*_BASE_HH_*.csv")):
    site=os.path.basename(f).split("_")[1]
    try:
        hdr=pd.read_csv(f,skiprows=2,nrows=0).columns.tolist()
        nee,ta,ts,sw=pick_nee(hdr),pick_ta(hdr),pick_ts(hdr),pick_sw(hdr)
        if not (nee and ta and ts and sw): continue
        df=pd.read_csv(f,skiprows=2,usecols=[nee,ta,ts,sw],na_values=[-9999,"-9999"])
        df=df.rename(columns={nee:"NEE",ta:"TA",ts:"TS",sw:"SW"}).dropna()
        n=df[(df.SW<20)&(df.NEE>0)]
        if len(n)<150 or n.TA.std()<0.5 or n.TS.std()<0.5: continue
        y=np.log(n.NEE.values); TA=n.TA.values; TS=n.TS.values
        def r2(T):
            b=np.polyfit(T,y,1); p=np.polyval(b,T)
            return 1-np.sum((y-p)**2)/np.sum((y-y.mean())**2), float(np.exp(10*b[0]))
        r2a,q10a=r2(TA); r2s,q10s=r2(TS)
        # competitive: standardized multiple regression
        za=(TA-TA.mean())/TA.std(); zs=(TS-TS.mean())/TS.std()
        X=np.column_stack([za,zs,np.ones_like(za)])
        beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
        rows.append(dict(site=site, igbp=igbp.get(site,"?"), n=len(n),
                         r2_air=r2a, r2_soil=r2s, q10_air=q10a, q10_soil=q10s,
                         beta_air=beta[0], beta_soil=beta[1]))
    except Exception:
        continue

R=pd.DataFrame(rows)
print(f"=== EC respiration-driver validation (AmeriFlux nighttime), {len(R)} sites ===\n")
def frac(mask): return f"{mask.sum()}/{len(R)} ({100*mask.mean():.0f}%)"
print("(1) SEPARATE R^2 of ln(nighttime resp) vs temperature:")
print(f"    median R^2 air {R.r2_air.median():.3f} | soil {R.r2_soil.median():.3f} | "
      f"median ΔR^2 (soil-air) {(R.r2_soil-R.r2_air).median():+.3f}")
print(f"    soil R^2 > air R^2 : {frac(R.r2_soil>R.r2_air)}")
print("(2) COMPETITIVE multiple regression (both temps, standardized) -- which wins:")
print(f"    |beta_soil| > |beta_air| : {frac(R.beta_soil.abs()>R.beta_air.abs())}")
print(f"    median |beta_soil| {R.beta_soil.abs().median():.3f} vs |beta_air| {R.beta_air.abs().median():.3f}")
print("(3) implied Q10 (physical range ~1.5-3):")
print(f"    median Q10 air {R.q10_air.median():.2f} | soil {R.q10_soil.median():.2f}")
print(f"    soil Q10 in [1.3,4] : {frac(R.q10_soil.between(1.3,4))} ; air in [1.3,4] : {frac(R.q10_air.between(1.3,4))}")
print("\n(4) by IGBP biome (soil R^2 > air, median ΔR^2):")
for b,g in R.groupby("igbp"):
    if len(g)>=3:
        print(f"    {b:5s} n={len(g):2d}: soil-wins {100*(g.r2_soil>g.r2_air).mean():3.0f}%  medianΔR² {(g.r2_soil-g.r2_air).median():+.3f}")
R.sort_values("r2_soil",ascending=False).to_csv("fitter_diagnostics/ec_resp_driver_validation_sites.csv",index=False)
print(f"\nper-site table ({len(R)} sites) -> fitter_diagnostics/ec_resp_driver_validation_sites.csv")
print("DONE")
