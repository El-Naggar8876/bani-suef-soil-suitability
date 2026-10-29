"""Recompute every quantitative claim in the submitted manuscript."""
import pathlib, numpy as np, pandas as pd, geopandas as gpd, rasterio
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

REPO = pathlib.Path(__file__).resolve().parents[1].parent / "Bani_Suief_Paper_FF"
S1,S4,S5,S6,S7,S7B,S3B = (REPO/"outputs"/s for s in ("stage1","stage4","stage5","stage6","stage7","stage7b","stage3b"))
def rd(p):
    with rasterio.open(p) as s:
        a=s.read(1).astype('float64')
        if s.nodata is not None and np.isfinite(s.nodata): a[a==s.nodata]=np.nan
        return a
out=[]
def chk(claim, manuscript, computed, ok=None, note=""):
    if ok is None:
        try: ok = abs(float(manuscript)-float(computed))<=max(0.01,abs(float(computed))*0.015)
        except Exception: ok = str(manuscript)==str(computed)
    out.append(dict(claim=claim, in_manuscript=manuscript, recomputed=computed,
                    match="OK" if ok else "MISMATCH", note=note))

# --- study area ------------------------------------------------------------
sa = gpd.read_file(REPO/"Layers/Study_Area_Last.shp")
if sa.crs is None or sa.crs.to_epsg()!=32636: sa=sa.to_crs(32636)
chk("study-area polygon area (km2)", 641.98, round(float(sa.area.sum())/1e6,2))
idx=rd(S5/"ales_index_q50_0_100cm.tif"); valid=np.isfinite(idx)
chk("valid pixels", 715784, int(valid.sum()))
chk("valid mask area (km2)", 641.98, round(int(valid.sum())*0.0009,2), note="valid pixels x 900 m2")
with rasterio.open(S5/"ales_index_q50_0_100cm.tif") as s: H,W=s.height,s.width
chk("raster dimensions", "1107 x 1379", f"{W} x {H}")
chk("pct pixels valid", 47, round(100*valid.sum()/(H*W),1))

# --- baseline ALES ---------------------------------------------------------
for dep, mq50, mp in [("0_30cm",17.7,0.08), ("0_100cm",29.9,0.19)]:
    a=rd(S5/f"ales_index_q50_{dep}.tif"); v=a[np.isfinite(a)]
    p=rd(S5/f"p_at_least_S2_{dep}.tif"); pv=p[np.isfinite(p)]
    chk(f"mean ALES q50 {dep}", mq50, round(float(v.mean()),2))
    chk(f"mean P(>=S2) {dep}", mp, round(float(pv.mean()),3))
    q05=rd(S5/f"ales_index_q05_{dep}.tif"); q95=rd(S5/f"ales_index_q95_{dep}.tif")
    w=(q95-q05); w=w[np.isfinite(w)]
    chk(f"mean PI width {dep}", 30.2 if dep=="0_30cm" else 33.0, round(float(w.mean()),2))
    cls=rd(S5/f"ales_class_q50_{dep}.tif")
    tot=np.isfinite(a).sum()
    for c,nm,man in [(3,"S3",4.2 if dep=="0_30cm" else 55.5),(4,"N1",95.7 if dep=="0_30cm" else 44.5)]:
        chk(f"{nm} share {dep} (%)", man, round(100*float(((cls==c)&valid).sum())/tot,1))

# --- counterfactuals -------------------------------------------------------
for dep in ["0_30cm","0_100cm"]:
    for sc,man30,man100 in [("CF_EC",1.5,2.1),("CF_Ks",31.0,21.1),("CF_SOM",0.1,0.3),("CF_full",34.9,24.8)]:
        a=rd(S6/f"cf_{sc}_{dep}_delta_q50.tif"); v=a[np.isfinite(a)]
        chk(f"{sc} mean delta {dep}", man30 if dep=="0_30cm" else man100, round(float(v.mean()),2))
# Ks share of closable gap
for dep,man in [("0_30cm",89),("0_100cm",85)]:
    k=rd(S6/f"cf_CF_Ks_{dep}_delta_q50.tif"); f=rd(S6/f"cf_CF_full_{dep}_delta_q50.tif")
    m=np.isfinite(k)&np.isfinite(f)
    chk(f"Ks share of CF_full {dep} (%)", man, round(100*float(k[m].mean()/f[m].mean()),1))

# --- land use --------------------------------------------------------------
lu=rd(S7/"landuse_class.tif"); lv=lu[np.isfinite(lu)]
for code,nm,man in [(2,"double crop",39.5),(3,"single crop",9.6),(4,"bare desert",14.4),(5,"sparse veg",7.8),(1,"water",2.1)]:
    chk(f"land use {nm} (%)", man, round(100*float((lv==code).mean()),1), note="denominator = all raster pixels incl. outside mask")
tot_all=lu.size
for code,nm,man in [(2,"double crop",39.5),(3,"single crop",9.6)]:
    chk(f"land use {nm} (% of full grid)", man, round(100*float((lu==code).sum())/tot_all,1))

# --- policy priority -------------------------------------------------------
pri=rd(S7/"policy_priority_0_100cm.tif")
nh=int((pri==1).sum())
chk("HIGH-priority pixels", 177052, nh)
chk("HIGH-priority area (km2)", 159, round(nh*0.0009,1))
cropped=np.isin(lu,[2,3])&valid
chk("HIGH as % of cropped land", 32.7, round(100*nh/int(cropped.sum()),1))
chk("HIGH as % of mapped area", 24.7, round(100*nh/int(valid.sum()),1))
cfd=rd(S6/"cf_CF_Ks_0_100cm_delta_q50.tif")
n1=(rd(S5/"ales_class_q50_0_100cm.tif")==4)
thr=np.nanpercentile(cfd[cropped&n1],75)
chk("top-quartile CF_Ks threshold", 27.8, round(float(thr),1))

# --- cropped fraction by class --------------------------------------------
cls=rd(S5/"ales_class_q50_0_100cm.tif")
for c,nm,man in [(4,"N1",0.999),(3,"S3",0.56)]:
    m=(cls==c)&valid
    chk(f"cropped fraction of {nm} 0-100cm", man, round(float(np.isin(lu[m],[2,3]).mean()),3))
chk("N1 pixel count 0-100cm", 318336, int(((cls==4)&valid).sum()))
cls30=rd(S5/"ales_class_q50_0_30cm.tif")
m=(cls30==4)&valid
chk("cropped fraction of N1 0-30cm", 0.78, round(float(np.isin(lu[m],[2,3]).mean()),3))

# --- IWQI ------------------------------------------------------------------
w=pd.read_csv(S1/"water_samples_iwqi.csv")
col=[c for c in w.columns if "iwqi" in c.lower()][0]
chk("IWQI sample min", 63.07, round(float(w[col].min()),2))
chk("IWQI sample max", 86.38, round(float(w[col].max()),2))
chk("IWQI sample mean", 82.01, round(float(w[col].mean()),2))
iq=rd(S4/"iwqi_q50.tif"); iv=iq[np.isfinite(iq)]
chk("kriged IWQI min", 70.1, round(float(iv.min()),1))
chk("kriged IWQI max", 84.7, round(float(iv.max()),1))
chk("kriged IWQI mean", 80.9, round(float(iv.mean()),1))
i5=rd(S4/"iwqi_q05.tif"); i9=rd(S4/"iwqi_q95.tif"); ww=(i9-i5); ww=ww[np.isfinite(ww)]
chk("IWQI mean interval width", 28.6, round(float(ww.mean()),1))

# --- Table 2 source data ---------------------------------------------------
pf=pd.read_csv(S1/"soil_per_profile_0_100cm.csv"); pf["Sand_total"]=pf.Coarse_sand+pf.Fine_sand
for prop,mn,md,me,mx in [("EC",1.07,3.0,3.5,8.43),("ESP",3.9,6.4,6.5,9.4),
                         ("CaCO3",1.7,7.0,7.78,21.9),("SOM",0.65,0.85,0.86,1.10),
                         ("CEC",9.08,32.80,29.16,38.33),("Clay",11.42,50.43,43.46,53.75),
                         ("Sand_total",24.58,33.33,37.48,75.32),("Ks",0.06,0.31,0.85,6.06)]:
    s=pf[prop].dropna()
    chk(f"Table 2 {prop} min", mn, round(float(s.min()),2))
    chk(f"Table 2 {prop} median", md, round(float(s.median()),2))
    chk(f"Table 2 {prop} mean", me, round(float(s.mean()),2))
    chk(f"Table 2 {prop} max", mx, round(float(s.max()),2))
chk("Ks 75th percentile", 0.94, round(float(pf.Ks.quantile(0.75)),2))
chk("Ks IQR lower", 0.18, round(float(pf.Ks.quantile(0.25)),2))
chk("pct profiles Ks<1.5", 75, round(100*float((pf.Ks<1.5).mean()),1), note="manuscript says 'over 75%'")

# --- crop survey -----------------------------------------------------------
cp=pd.read_csv(S7B/"crop_points_extracted.csv")
chk("survey points", 993, len(cp))
chk("points on cropland", 981, int(cp.landuse.isin([2,3]).sum()))
chk("user's accuracy", 0.988, round(float(cp.landuse.isin([2,3]).mean()),3))
chk("pct on double crop", 89.7, round(100*float((cp.landuse==2).mean()),1))
chk("pct on single crop", 9.1, round(100*float((cp.landuse==3).mean()),1))
k=pd.read_csv(S7B/"crop_kpis.csv").set_index("crop")
for crop,n,pn1,mdk in [("Wheat",322,52.2,23.24),("Clover",523,38.4,21.88),
                       ("Citrus",30,43.3,23.17),("Banana",3,66.7,25.01),
                       ("Sugar Cane",15,33.3,20.76),("Medical Plants",21,66.7,23.93)]:
    chk(f"{crop} n", n, int(k.loc[crop,"n"]))
    chk(f"{crop} %N1+N2", pn1, round(float(k.loc[crop,"pct_N1_or_N2_%"]),1))
    chk(f"{crop} mean CF_Ks delta", mdk, round(float(k.loc[crop,"mean_delta_Ks"]),2))
# the %HIGH denominator problem
xt=pd.read_csv(S7B/f"crop_x_priority_0_100cm.csv").set_index("CROP")
for crop,man in [("Wheat",52.9),("Citrus",53.8),("Clover",36.6),("Medical Plants",14.3)]:
    r=xt.loc[crop]; among_pri=100*r.HIGH/(r.HIGH+r.MEDIUM+r.WATCH)
    among_all=100*r.HIGH/float(k.loc[crop,"n"])
    chk(f"{crop} %HIGH (pipeline denominator = priority-classified points)", man, round(float(among_pri),1))
    chk(f"{crop} %HIGH (share of ALL points of that crop)", man, round(float(among_all),1),
        note="manuscript wording implies this denominator")

df=pd.DataFrame(out)
df.to_csv(pathlib.Path(__file__).parent/"outputs_rev"/"manuscript_number_check.csv",index=False)
bad=df[df.match=="MISMATCH"]
print(f"checked {len(df)} claims | OK {len(df)-len(bad)} | MISMATCH {len(bad)}\n")
print("================ MISMATCHES ================")
print(bad.to_string(index=False))
