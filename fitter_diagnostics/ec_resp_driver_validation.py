#!/usr/bin/env python3
"""Eddy-covariance validation of the respiration temperature driver (soil vs air) -- v2,
hardened after adversarial review. Motivating analysis for the OPT-IN soil-temp
diurnalization driver (DIURNALIZATION_ALTERNATIVES.md); NOT part of the V1->V2 switch.

At night (SW<20 W/m^2) NEE ~= ecosystem respiration: no GPP, no partitioning model.
We test whether nighttime respiration is better explained by 2-m AIR temperature or by
SHALLOW SOIL temperature, separating two DISTINCT questions:

  TEST A -- SEASONAL driver. Whole-record nighttime ln(NEE) vs T. Dominated by the
    winter-summer temperature range. Answers: which temperature tracks the seasonal
    magnitude of respiration. (This is what the old script measured.)

  TEST B -- WITHIN-DAY driver (the question the diurnalization actually uses). Remove
    each NIGHT's mean from ln(NEE), TA, TS; regress the within-night anomalies. Answers:
    which temperature drives the SUB-DAILY shape of respiration. Soil temp has a small,
    lagged diurnal amplitude, so A and B can disagree -- and B is the relevant test for
    setting the diurnal cycle. CAVEAT: only the nighttime arc of the diurnal cycle is
    observable (daytime RECO is hidden under GPP); B tests that arc.

Hardening vs v1:
  * RAW sensors only. Gap-filled (_PI_F / _F) FORBIDDEN for NEE (nighttime gap-fill IS a
    temperature-driven RECO model -> circular), TS and SW (modeled fill contaminates the
    driver and the night mask). Robust naming: VAR[_PI][_F]_H_V_R, R digit or letter.
  * u* (USTAR) filter: drop low-turbulence nights (storage/drainage inflate apparent resp).
  * Shallowest soil temp by vertical index (closest to ERA5 0-7 cm stl1).
  * Block bootstrap by NIGHT for honest CIs (half-hourly data is autocorrelated ~0.98;
    i.i.d. N is fiction). Per-site win is "decisive" only if the 95% CI of dR2 excludes 0;
    otherwise "tie". Significance across sites = binomial on decisive sites only.
  * Degenerate sites (max R2 < 0.02 -- no temperature signal) dropped, not coin-flipped.
  * Report r(TA,TS) and VIF (collinearity inflates the competitive betas).
  * Report ALL biomes and per-site DROP REASON. Flag Q10>3.5 as unphysical.
"""
import numpy as np, pandas as pd, glob, os, re, warnings
warnings.simplefilter("ignore")
ROOT = "/work2/noaa/co2/kaushik/drought/ameriflux_US"
USTAR_MIN = 0.2          # m/s, standard nighttime turbulence threshold
SWNIGHT   = 20.0         # W/m^2, night definition
MIN_NIGHTS = 30          # need enough night-blocks to bootstrap
MIN_PTS    = 150
TIE_R2     = 0.01        # |dR2| below this with CI spanning 0 = tie
NBOOT      = 400
SEED       = 12345

# ---- IGBP biome per site -------------------------------------------------
igbp = {}
try:
    meta = pd.read_csv(f"{ROOT}/NAm_US_sites.csv", encoding="latin-1", engine="python", on_bad_lines="skip")
    for _, r in meta.iterrows():
        try: igbp[str(r["Site Id"]).strip()] = str(r["Vegetation Abbreviation (IGBP)"]).strip()
        except Exception: pass
except Exception as e:
    print("biome metadata unavailable:", e)

# ---- robust column selection (RAW preferred; gap-fill forbidden) ---------
def _parse(c, base):
    # match VAR[_PI][_F]_H_V_R with R a digit or single letter; return (rawrank, vert) or None
    m = re.match(rf"^{base}(_PI)?(_F)?_(\d+)_(\d+)_(\d+|[A-Z])$", c)
    if not m: return None
    is_pi, is_f = bool(m.group(1)), bool(m.group(2))
    rawrank = 0 if (not is_pi and not is_f) else (1 if (is_pi and not is_f) else 2)  # raw<pi<gapfilled
    return rawrank, int(m.group(4))           # group(4) = vertical (depth) index

def pick(hdr, base, shallowest=False, allow_gapfill=False):
    cands = []
    for c in hdr:
        if "QC" in c or "SSITC" in c: continue
        p = _parse(c, base)
        if p is None: continue
        rawrank, vert = p
        if rawrank == 2 and not allow_gapfill: continue       # forbid gap-filled
        key = (vert if shallowest else 0, rawrank)
        cands.append((key, c))
    if not cands: return None
    cands.sort()
    return cands[0][1]

def pick_nee(hdr):
    # measured turbulent nighttime flux ONLY -- never gap-filled (that is a modeled RECO)
    for base in ("NEE", "FC"):
        c = pick(hdr, base, allow_gapfill=False)
        if c: return c
    return None

# ---- per-site load + tests ----------------------------------------------
rng = np.random.default_rng(SEED)
rows, drops = [], []

def fit_r2_q10(T, y):
    b = np.polyfit(T, y, 1); p = np.polyval(b, T)
    ss = np.sum((y - y.mean())**2)
    r2 = 1.0 - np.sum((y - p)**2)/ss if ss > 0 else 0.0
    return r2, float(np.exp(10*b[0]))

def dR2(T_air, T_soil, y):
    return fit_r2_q10(T_soil, y)[0] - fit_r2_q10(T_air, y)[0]

for f in sorted(glob.glob(f"{ROOT}/EC_sitedata/*/AMF_*_BASE_HH_*.csv")):
    site = os.path.basename(f).split("_")[1]
    try:
        hdr = pd.read_csv(f, skiprows=2, nrows=0).columns.tolist()
        nee = pick_nee(hdr)
        ta  = pick(hdr, "TA", allow_gapfill=True)      # air temp: raw preferred, gapfill last resort
        ts  = pick(hdr, "TS", shallowest=True)          # soil temp: raw only, shallowest
        sw  = pick(hdr, "SW_IN")                          # raw only
        ust = pick(hdr, "USTAR")
        miss = [n for n,v in [("NEE",nee),("TA",ta),("TS",ts),("SW",sw),("USTAR",ust)] if not v]
        if miss:
            drops.append((site, "missing:"+",".join(miss))); continue
        use = ["TIMESTAMP_START", nee, ta, ts, sw, ust]
        df = pd.read_csv(f, skiprows=2, usecols=use, na_values=[-9999, "-9999"])
        df = df.rename(columns={nee:"NEE", ta:"TA", ts:"TS", sw:"SW", ust:"USTAR"})
        df = df.dropna(subset=["TIMESTAMP_START","NEE","TA","TS","SW","USTAR"])
        # night id: shift back 12 h so a full night groups under one date
        t = pd.to_datetime(df.TIMESTAMP_START.astype("int64").astype(str), format="%Y%m%d%H%M")
        df["night"] = (t - pd.Timedelta(hours=12)).dt.strftime("%Y%m%d")
        df["season"] = t.dt.month.map({12:"DJF",1:"DJF",2:"DJF",3:"MAM",4:"MAM",5:"MAM",
                                       6:"JJA",7:"JJA",8:"JJA",9:"SON",10:"SON",11:"SON"})
        n = df[(df.SW < SWNIGHT) & (df.USTAR > USTAR_MIN) & (df.NEE > 0)].copy()
        nnights = n.night.nunique()
        if len(n) < MIN_PTS or nnights < MIN_NIGHTS:
            drops.append((site, f"few_pts(n={len(n)},nights={nnights})")); continue
        if n.TA.std() < 0.5 or n.TS.std() < 0.5:
            drops.append((site, "low_temp_std")); continue

        y  = np.log(n.NEE.values); TA = n.TA.values; TS = n.TS.values
        # TEST A -- seasonal (whole-record)
        r2a, q10a = fit_r2_q10(TA, y)
        r2s, q10s = fit_r2_q10(TS, y)
        if max(r2a, r2s) < 0.02:
            drops.append((site, f"degenerate(maxR2={max(r2a,r2s):.3f})")); continue
        # competitive standardized regression + collinearity
        za = (TA-TA.mean())/TA.std(); zs = (TS-TS.mean())/TS.std()
        rcorr = float(np.corrcoef(TA, TS)[0,1]); vif = 1.0/max(1e-6, 1.0 - rcorr**2)
        beta = np.linalg.lstsq(np.column_stack([za, zs, np.ones_like(za)]), y, rcond=None)[0]

        # TEST B -- within-night anomalies (remove each night's mean)
        g = n.groupby("night")
        ya = y - g["NEE"].transform(lambda v: np.log(v).mean()).values
        TAa = TA - g["TA"].transform("mean").values
        TSa = TS - g["TS"].transform("mean").values
        # keep only nights with within-night temperature spread (else anomalies are noise)
        spread = (g["TA"].transform("std").values > 0.2) | (g["TS"].transform("std").values > 0.2)
        m = spread & np.isfinite(ya) & np.isfinite(TAa) & np.isfinite(TSa)
        r2a_w = fit_r2_q10(TAa[m], ya[m])[0] if m.sum() > MIN_PTS else np.nan
        r2s_w = fit_r2_q10(TSa[m], ya[m])[0] if m.sum() > MIN_PTS else np.nan

        # block bootstrap by NIGHT -> CI on dR2 for A and B
        nights = n.night.values
        uniq = np.unique(nights)
        idx_by_night = {u: np.where(nights == u)[0] for u in uniq}
        bootA, bootB = [], []
        for _ in range(NBOOT):
            samp = rng.choice(uniq, size=len(uniq), replace=True)
            ii = np.concatenate([idx_by_night[u] for u in samp])
            bootA.append(dR2(TA[ii], TS[ii], y[ii]))
            mb = m[ii]
            bootB.append(dR2(TAa[ii][mb], TSa[ii][mb], ya[ii][mb]) if mb.sum() > MIN_PTS else np.nan)
        loA, hiA = np.percentile(bootA, [2.5, 97.5])
        bB = np.array(bootB); bB = bB[np.isfinite(bB)]
        loB, hiB = (np.percentile(bB, [2.5, 97.5]) if len(bB) > 20 else (np.nan, np.nan))

        def verdict(d, lo, hi):
            if not np.isfinite(lo): return "na"
            if lo > 0: return "soil"
            if hi < 0: return "air"
            return "tie"

        rows.append(dict(site=site, igbp=igbp.get(site, "?"), n=len(n), nights=nnights,
                         r2_air=r2a, r2_soil=r2s, dR2_seasonal=r2s-r2a,
                         seas_lo=loA, seas_hi=hiA, seasonal=verdict(r2s-r2a, loA, hiA),
                         r2_air_wd=r2a_w, r2_soil_wd=r2s_w, dR2_withinday=r2s_w-r2a_w,
                         wd_lo=loB, wd_hi=hiB, withinday=verdict(r2s_w-r2a_w, loB, hiB),
                         q10_air=q10a, q10_soil=q10s, beta_air=beta[0], beta_soil=beta[1],
                         r_ta_ts=rcorr, vif=vif))
    except Exception as e:
        drops.append((site, f"error:{type(e).__name__}")); continue

R = pd.DataFrame(rows)
print(f"=== EC respiration-driver validation v2 (AmeriFlux nighttime, u*>{USTAR_MIN}) ===")
print(f"{len(R)} sites pass; {len(drops)} dropped\n")

def binom(k, nn):
    from scipy.stats import binomtest
    return binomtest(k, nn, 0.5).pvalue

def summarize(col, lab):
    v = R[col]
    soil = (v == "soil").sum(); air = (v == "air").sum(); tie = (v == "tie").sum()
    dec = soil + air
    p = binom(soil, dec) if dec > 0 else np.nan
    print(f"{lab}: decisive soil {soil} | air {air} | tie {tie}  "
          f"(binomial soil-vs-air on {dec} decisive sites p={p:.4f})")

print("--- TEST A: SEASONAL driver (whole-record nighttime ln NEE vs T) ---")
print(f"  median R2 air {R.r2_air.median():.3f} | soil {R.r2_soil.median():.3f} | "
      f"median dR2(soil-air) {R.dR2_seasonal.median():+.3f}")
summarize("seasonal", "  block-bootstrap verdict")
print("--- TEST B: WITHIN-DAY driver (within-night anomalies -- the diurnalization question) ---")
print(f"  median R2 air {R.r2_air_wd.median():.3f} | soil {R.r2_soil_wd.median():.3f} | "
      f"median dR2(soil-air) {R.dR2_withinday.median():+.3f}")
summarize("withinday", "  block-bootstrap verdict")
print("--- collinearity / competitive regression ---")
print(f"  median r(TA,TS) {R.r_ta_ts.median():.2f} | median VIF {R.vif.median():.1f} "
      f"(high => competitive betas variance-inflated; treat as ONE evidence line w/ Test A)")
print(f"  median |beta_soil| {R.beta_soil.abs().median():.3f} vs |beta_air| {R.beta_air.abs().median():.3f}")
print("--- Q10 sanity (physical ~1.5-3; >3.5 flagged) ---")
print(f"  median Q10 air {R.q10_air.median():.2f} | soil {R.q10_soil.median():.2f}")
print(f"  Q10_soil>3.5 (unphysical): {(R.q10_soil>3.5).sum()}/{len(R)}  "
      f"sites={list(R[R.q10_soil>3.5].site)}")
print("--- by IGBP biome (ALL biomes; seasonal soil-win fraction) ---")
for b, gg in R.groupby("igbp"):
    note = "" if len(gg) >= 3 else "  (n<3, noisy)"
    sw_win = (gg.seasonal == "soil").sum(); wd_win = (gg.withinday == "soil").sum()
    print(f"  {b:5s} n={len(gg):2d}: seasonal soil {sw_win}/{len(gg)} | within-day soil {wd_win}/{len(gg)}{note}")
print("--- DROP REASONS ---")
dd = pd.Series([r for _, r in drops]).str.replace(r"\(.*\)", "", regex=True).str.split(":").str[0]
for reason, cnt in dd.value_counts().items():
    print(f"  {reason:18s}: {cnt}")
ndrop_ts = sum(1 for _, r in drops if "missing" in r and "TS" in r)
print(f"  (of 'missing', {ndrop_ts} lack a usable raw TS sensor -- selection caveat)")

R.sort_values("dR2_seasonal", ascending=False).to_csv(
    "fitter_diagnostics/ec_resp_driver_validation_sites.csv", index=False)
pd.DataFrame(drops, columns=["site","drop_reason"]).to_csv(
    "fitter_diagnostics/ec_resp_driver_validation_drops.csv", index=False)
print(f"\nper-site -> ec_resp_driver_validation_sites.csv ; drops -> ec_resp_driver_validation_drops.csv")
print("DONE")
