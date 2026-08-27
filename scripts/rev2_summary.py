#!/usr/bin/env python3
"""Constroi output/rev2/summary.csv a partir das cadeias refeitas (init largo, 1600
passos, tratamento identico entre as sete EoS). Mesmo esquema do summary original.
Colunas da cadeia: 0=ens 1=step 2=walker 3=log10L 4=L 5=alpha 6=w 7=logpost
                   8=Mmax 9=Rmax 10=R14 11=Cmax 12=rhoc 13=accepted"""
import csv, os, sys
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BURN=400   # 25% de 1600
NAMES=["SLy","BSK24","DD2","NL3","QHC19A","QHC21D","SPGM4"]
gr={r["eos"]:(float(r["Mmax_GR"]),float(r["R14_GR"]))
    for r in csv.DictReader(open(ROOT+"/output/campaign/full/summary.csv"))}
cols=["eos","samples","Mmax_GR","R14_GR","P_alpha_lt_0","P_alpha_gt_0_w_lt_0"]
for p in ("log10_lambda","alpha_U","w_U","Mmax","R14"):
    cols += [f"{p}_q16",f"{p}_q50",f"{p}_q84"]
out=[]; missing=[]
for nm in NAMES:
    f=f"{ROOT}/output/rev2/{nm}/mcmc_chain.dat"
    if not os.path.exists(f): missing.append(nm); continue
    d=np.loadtxt(f); d=d[d[:,1]>=BURN]
    a,w=d[:,5],d[:,6]
    row={"eos":nm,"samples":len(d),"Mmax_GR":gr[nm][0],"R14_GR":gr[nm][1],
         "P_alpha_lt_0":float(np.mean(a<0)),
         "P_alpha_gt_0_w_lt_0":float(np.mean((a>0)&(w<0)))}
    for p,idx in (("log10_lambda",3),("alpha_U",5),("w_U",6),("Mmax",8),("R14",10)):
        q=np.percentile(d[:,idx],[16,50,84])
        row[f"{p}_q16"],row[f"{p}_q50"],row[f"{p}_q84"]=q
    out.append(row)
if missing: print("  FALTAM (nao gravadas):", ", ".join(missing), file=sys.stderr)
os.makedirs(ROOT+"/output/rev2",exist_ok=True)
with open(ROOT+"/output/rev2/summary.csv","w",newline="") as fh:
    wtr=csv.DictWriter(fh,fieldnames=cols); wtr.writeheader(); wtr.writerows(out)
print("  summary.csv com %d EoS"%len(out))
for r in sorted(out,key=lambda x: x["R14_GR"]):
    print("   %-8s R14_GR=%5.2f  alpha=%+.3f [%+.3f,%+.3f]  P(a<0)=%4.1f%%"
          %(r["eos"],r["R14_GR"],r["alpha_U_q50"],r["alpha_U_q16"],r["alpha_U_q84"],
            100*r["P_alpha_lt_0"]))
