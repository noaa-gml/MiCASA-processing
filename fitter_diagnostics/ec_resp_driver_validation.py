#!/usr/bin/env python3
"""Eddy-covariance validation of the respiration temperature driver (soil vs air)
— the open gate for the soil-temp diurnalization flip (§2).

Does OBSERVED ecosystem respiration follow soil temperature or 2-m air temperature?

  TEST B (primary, non-circular, broad coverage): at night SW<20 -> NEE≈RECO (no GPP,
    no partitioning model). Fit ln(nighttime respiration) vs TA and vs TS; the
    temperature with the higher R^2 (and a physical Q10) is the better driver.
    Needs only NEE, TA, TS, SW -> available at most AmeriFlux BASE sites.
  TEST A (secondary): for the few sites that ship partitioned RECO_PI, compare the
    observed RECO diurnal-cycle amplitude to the air-Q10 vs soil-Q10 predictions.
    Flagged because RECO_PI partitioning may itself use a temperature (circular).
Data: /work2/noaa/co2/kaushik/drought/ameriflux_US/EC_sitedata/*/AMF_*_BASE_HH_*.csv
"""
import numpy as np, pandas as pd, glob, os, warnings
warnings.simplefilter("ignore")
ROOT = "/work2/noaa/co2/kaushik/drought/ameriflux_US/EC_sitedata"; Q10 = 1.5

def pick(cols, *cands):
    for c in cands:
        if c in cols: return c
    return None
def diurnal(df, col): return df.groupby("hour")[col].mean().reindex(range(24)).values
def amp(x): return np.nanmax(x)-np.nanmin(x)

night_rows, diur_rows = [], []
for f in sorted(glob.glob(f"{ROOT}/*/AMF_*_BASE_HH_*.csv")):
    site = os.path.basename(f).split("_")[1]
    try:
        hdr = pd.read_csv(f, skiprows=2, nrows=0).columns.tolist()
        ta = pick(hdr, "TA_PI_F_1_1_1", "TA_1_1_1")
        ts = pick(hdr, "TS_PI_F_1_1_1", "TS_1_1_1", "TS_PI_F_2_1_1", "TS_2_1_1") \
             or next((c for c in hdr if c.startswith("TS_") and "QC" not in c), None)
        nee= pick(hdr, "NEE_PI_F_1_1_1", "NEE_PI_1_1_1", "FC_PI_F_1_1_1", "FC_1_1_1")
        sw = pick(hdr, "SW_IN_PI_F_1_1_1", "SW_IN_1_1_1")
        reco = pick(hdr, "RECO_PI_F_1_1_1", "RECO_PI_1_1_1")
        if not (ta and ts and nee and sw): continue
        use = [c for c in {"TIMESTAMP_START", ta, ts, nee, sw, reco} if c]
        df = pd.read_csv(f, skiprows=2, usecols=use, na_values=[-9999, "-9999"])
        ren = {ta:"TA", ts:"TS", nee:"NEE", sw:"SW"}
        if reco: ren[reco] = "RECO"
        df = df.rename(columns=ren)
        df["hour"] = (df["TIMESTAMP_START"].astype("int64") // 100 % 100)
        # --- TEST B: nighttime respiration vs TA / TS ---
        n = df.dropna(subset=["NEE","TA","TS","SW"])
        n = n[(n["SW"] < 20) & (n["NEE"] > 0)]
        if len(n) > 200 and n["TS"].std() > 0.5 and n["TA"].std() > 0.5:
            y = np.log(n["NEE"].values); r2 = {}
            for T, k in ((n["TA"].values,"a"), (n["TS"].values,"s")):
                b = np.polyfit(T, y, 1); pr = np.polyval(b, T)
                r2[k] = 1 - np.sum((y-pr)**2)/np.sum((y-np.mean(y))**2)
            night_rows.append(dict(site=site, n=len(n), r2_air=r2["a"], r2_soil=r2["s"]))
        # --- TEST A: observed RECO diurnal amplitude (only if RECO present) ---
        if reco:
            d = df.dropna(subset=["RECO","TA","TS"])
            if len(d) >= 5000:
                rc, tac, tsc = diurnal(d,"RECO"), diurnal(d,"TA"), diurnal(d,"TS")
                if np.all(np.isfinite(rc)) and amp(rc) > 0:
                    rm = np.nanmean(rc)
                    pa = Q10**(tac/10); pa *= rm/np.nanmean(pa)
                    ps = Q10**(tsc/10); ps *= rm/np.nanmean(ps)
                    diur_rows.append(dict(site=site, obs_over_air=amp(rc)/amp(pa),
                                          soil_over_air=amp(ps)/amp(pa)))
    except Exception:
        continue

NB = pd.DataFrame(night_rows); DA = pd.DataFrame(diur_rows)
print(f"=== EC respiration-driver validation (AmeriFlux) ===\n")
print(f"TEST B — nighttime respiration vs temperature (non-circular): {len(NB)} sites")
if len(NB):
    print(f"  median R^2 vs AIR temp  : {NB.r2_air.median():.3f}")
    print(f"  median R^2 vs SOIL temp : {NB.r2_soil.median():.3f}")
    w = (NB.r2_soil > NB.r2_air).sum()
    print(f"  sites where SOIL explains nighttime respiration better : {w}/{len(NB)} ({100*w/len(NB):.0f}%)")
    print(f"  median ΔR^2 (soil − air): {(NB.r2_soil-NB.r2_air).median():+.3f}")
print(f"\nTEST A — observed RECO diurnal amplitude vs air/soil Q10: {len(DA)} sites (RECO_PI only; partitioning caveat)")
if len(DA):
    print(f"  median observed amp / air-Q10 amp : {DA.obs_over_air.median():.2f}  (1=tracks air, <1=damped toward soil)")
    print(f"  median soil-Q10 amp / air-Q10 amp : {DA.soil_over_air.median():.2f}")
NB.to_csv("fitter_diagnostics/ec_resp_driver_validation_sites.csv", index=False)
print("\nper-site night table -> fitter_diagnostics/ec_resp_driver_validation_sites.csv")
print("DONE")
