#!/usr/bin/env python3
"""Fig. nova do paper 2: alpha_U inferido contra a linha de base de RG.
Le output/campaign/full/summary.csv. Sem solver, sem cadeia --- so replotagem.
R-dimensional: R14 em km, alpha_U adimensional, inclinacao em km^-1."""
import csv, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from math import factorial
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISP={"QHC19A":"QHC19-A","SLy":"SLy","SPGM4":"SPG(M4)","QHC21D":"QHC21-D",
      "BSK24":"BSk24","DD2":"DD2","NL3":"NL3"}
rows=sorted(csv.DictReader(open(ROOT+"/output/rev2/summary.csv")),
            key=lambda r: float(r["R14_GR"]))
R=np.array([float(r["R14_GR"]) for r in rows])
a=np.array([float(r["alpha_U_q50"]) for r in rows])
lo=a-np.array([float(r["alpha_U_q16"]) for r in rows])
hi=np.array([float(r["alpha_U_q84"]) for r in rows])-a
nm=[DISP.get(r["eos"],r["eos"]) for r in rows]
sl,ic=np.polyfit(R,a,1); cross=-ic/sl
pr=pearsonr(R,a); sp=spearmanr(R,a)
pexact=2.0/factorial(len(R))
plt.rcParams.update({"font.size":11,"axes.labelsize":12,"mathtext.fontset":"cm"})
fig,ax=plt.subplots(figsize=(6.4,4.6))
xs=np.linspace(R.min()-0.35,R.max()+0.35,60)
ax.plot(xs,sl*xs+ic,"-",color="#b23a48",lw=1.5,zorder=1)
ax.axhline(0,color="0.35",ls="--",lw=1.2,zorder=0)
ax.errorbar(R,a,yerr=[lo,hi],fmt="o",ms=6,color="#1f5a94",capsize=3,lw=1.3,zorder=3)
ax.plot([cross],[0.0],"s",color="#b23a48",ms=8,mfc="white",mew=1.8,zorder=4)
off={"QHC19-A":(-7,-15),"SLy":(11,-3),"SPG(M4)":(-7,9),"QHC21-D":(5,-17),
     "BSk24":(8,7),"DD2":(9,6),"NL3":(-11,7)}
for k,x,y in zip(nm,R,a):
    d=off.get(k,(6,6))
    ax.annotate(k,(x,y),textcoords="offset points",xytext=d,fontsize=9,
                ha="right" if d[0]<0 else "left")
ax.annotate(r"$\alpha_{\mathcal{U}}=0$ at $R_{1.4}^{\rm GR}=%.2f$ km"%cross,
            xy=(cross,0.0), xytext=(13.75,-0.155), fontsize=9.5, color="#b23a48",
            ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color="#b23a48", lw=1.1,
                            connectionstyle="arc3,rad=-0.2"))
ax.set_xlabel(r"$R_{1.4}^{\rm GR}$ [km]"); ax.set_ylabel(r"$\alpha_{\mathcal{U}}$")
ax.grid(alpha=0.2)
for s_ in ("top","right"): ax.spines[s_].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT+"/paper2/alpha_vs_R14GR.png",dpi=220)
out=dict(slope_per_km=float(sl),intercept=float(ic),zero_crossing_km=float(cross),
         pearson_r=float(pr[0]),pearson_p=float(pr[1]),
         spearman_rho=float(sp.statistic),exact_permutation_p=float(pexact),n=len(R))
json.dump(out,open(ROOT+"/output/rev2/alpha_vs_R14.json","w"),indent=2)
print("  inclinacao %+.4f /km   intercepto %+.3f   cruza zero em %.2f km"%(sl,ic,cross))
print("  Pearson r = %+.4f (p=%.2e)   Spearman rho = %+.3f   p exato de permutacao = %.1e"
      %(pr[0],pr[1],sp.statistic,pexact))
