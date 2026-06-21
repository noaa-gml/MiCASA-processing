#!/usr/bin/env python3
"""Figures for the PIQS+linear-fallback rebuttal in V1_TO_V2_JUSTIFICATION.md."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"figure.dpi":130})

# Fig 1: hybrid patch-discontinuity distribution (jump / local flux envelope)
v = pd.read_csv("fitter_diagnostics/linear_fallback_discontinuity.csv")["jump_over_env"].values
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.hist(np.clip(v,0,12),bins=60,color="#9d0208",alpha=.85)
ax.axvline(1,color="k",ls=":",lw=2,label="1× env (jump = entire local flux)")
ax.axvline(3,color="k",ls="--",lw=2,label="3× env")
ax.set_xlabel("patch discontinuity / local flux envelope")
ax.set_ylabel("patched cell-edges")
ax.set_title("PIQS+linear-fallback injects large flux discontinuities\n52% of patched edges exceed 3× the local flux; PCHIP = 0 (C⁰)")
ax.legend(); ax.grid(alpha=.3); fig.tight_layout()
fig.savefig("docs/figures/linear_fallback_discontinuity.png"); plt.close(fig)

# Fig 2: sign-safety vs continuity tradeoff -- PCHIP dominates the good corner
# x = sub-monthly wrong-sign rate (%), y = month-edge discontinuity (median /env)
methods = {
  "PIQS (V1)":        dict(x=29.3, y=0.0,  c="#c1121f", note="overshoot pieces; but 302-mo NRT rewrite"),
  "continuous linear":dict(x=36.9, y=0.0,  c="#e85d04", note="knot-flip; unbounded resonance"),
  "PIQS+linear\nfallback (Andy's)": dict(x=0.0, y=5.0, c="#9d0208", note="0 sign-flip but 5×-env jumps; 302-mo NRT"),
  "PPM":              dict(x=0.0,  y=0.018,c="#6a4c93", note="small edge jumps"),
  "PCHIP (V2)":       dict(x=0.65, y=0.0,  c="#0353a4", note="0.65% segs, C⁰, local"),
}
fig,ax=plt.subplots(figsize=(7.0,4.6))
off={"PIQS (V1)":(8,-4),"continuous linear":(8,-4),"PIQS+linear\nfallback (Andy's)":(10,-2),
     "PPM":(-8,16),"PCHIP (V2)":(8,-16)}
ha={"PPM":"right"}
for m,d in methods.items():
    ax.scatter(d["x"],d["y"],s=170,color=d["c"],zorder=3,edgecolor="k",linewidth=.6)
    ax.annotate(m,(d["x"],d["y"]),xytext=off.get(m,(6,6)),textcoords="offset points",
                fontsize=10,weight="bold",ha=ha.get(m,"left"))
ax.set_xlabel("sub-monthly wrong-sign rate (%)  →  worse")
ax.set_ylabel("month-edge discontinuity (median / env)  →  worse")
ax.set_title("PCHIP is sign-safe AND continuous; each alternative sacrifices one",fontsize=11,pad=10)
ax.set_xlim(-2,40); ax.set_ylim(-0.6,6); ax.grid(alpha=.3)
ax.axhspan(-0.4,0.3,xmin=0,xmax=0.1,color="#0353a4",alpha=.07)
fig.tight_layout(); fig.savefig("docs/figures/fitter_tradeoff_scatter.png"); plt.close(fig)
print("Wrote docs/figures/linear_fallback_discontinuity.png and fitter_tradeoff_scatter.png")
