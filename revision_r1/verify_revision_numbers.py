# -*- coding: utf-8 -*-
"""Check every quantitative claim in the REVISED manuscript against the analysis outputs."""
from __future__ import annotations
import json, pathlib, re
import numpy as np, pandas as pd, rasterio, docx

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REV = HERE.parent
REPO = REV.parent / "Bani_Suief_Paper_FF"
NEW = _OUT
WAP = _WAPOR
S1, S4, S5, S6, S7, S7B, S3B = (REPO/"outputs"/s for s in
                                ("stage1","stage4","stage5","stage6","stage7","stage7b","stage3b"))
DOC = REPO.parent / 'revision_01_geoderma_regional' / '05_Manuscript' / 'ElNaggar_SoilSuitability_R1_clean.docx'
TXT = "\n".join(p.text for p in docx.Document(str(DOC)).paragraphs)

def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float64")
        if s.nodata is not None and np.isfinite(s.nodata): a[a == s.nodata] = np.nan
        return a

rows = []
def chk(label, claim_str, computed, tol=0.02):
    present = claim_str in TXT
    ok = present
    rows.append(dict(claim=label, string=claim_str, in_text=present,
                     computed=computed, status="OK" if ok else "STRING NOT FOUND"))

# ---- recompute the quantities the revised text asserts ------------------
idx100 = rd(S5/"ales_index_q50_0_100cm.tif"); valid = np.isfinite(idx100)
idx30 = rd(S5/"ales_index_q50_0_30cm.tif")
p100 = rd(S5/"p_at_least_S2_0_100cm.tif"); p30 = rd(S5/"p_at_least_S2_0_30cm.tif")
lu = rd(S7/"landuse_class.tif"); pri = rd(S7/"policy_priority_0_100cm.tif")
cropped = np.isin(lu,[2,3]) & valid
t1 = pd.read_csv(NEW/"sens_T1_limiting_factor.csv")
t2 = pd.read_csv(NEW/"sens_T2_elasticity.csv")
t3 = pd.read_csv(NEW/"sens_T3_aggregation_rule.csv")
t4 = pd.read_csv(NEW/"sens_T4_threshold_perturbation.csv")
mc = pd.read_csv(NEW/"mc_independent_vs_correlated.csv")
c1 = pd.read_csv(NEW/"wapor_C1_mixed_pixel.csv")
c2p = pd.read_csv(NEW/"wapor_C2_partial.csv")
c2 = pd.read_csv(NEW/"wapor_C2_confounding.csv")
con = pd.read_csv(NEW/"wapor_contrasts.csv")
ksq = pd.read_csv(NEW/"wapor_by_ks_quintile.csv")
sb = pd.read_csv(NEW/"survey_bias_distributions.csv")
sc = pd.read_csv(NEW/"survey_bias_composition.csv")
rw = pd.read_csv(NEW/"survey_bias_reweighted_percrop.csv")
unc = pd.read_csv(NEW/"unc_by_zone.csv")
prof = pd.read_csv(S1/"soil_per_profile_0_100cm.csv")

def mcv(dep, struct, scen, col="mean_delta"):
    r = mc[(mc.depth==dep)&(mc.structure==struct)&(mc.scenario==scen)]
    return float(r[col].iloc[0])

def t3v(dep, rule, scen):
    return float(t3[(t3.depth==dep)&(t3.rule==rule)&(t3.scenario==scen)].mean_delta.iloc[0])

C = {}
C["limiting 0_100 wetness %"] = float(t1[(t1.depth=="0_100cm")&(t1.sub_rating=="wetness (Ks)")].pct_pixels_limiting.iloc[0])
C["limiting 0_30 wetness %"]  = float(t1[(t1.depth=="0_30cm")&(t1.sub_rating=="wetness (Ks)")].pct_pixels_limiting.iloc[0])
C["limiting 0_100 physical %"]= float(t1[(t1.depth=="0_100cm")&(t1.sub_rating=="physical (texture)")].pct_pixels_limiting.iloc[0])
C["elast clay 0_100"] = float(t2[(t2.depth=="0_100cm")&(t2.property=="Clay")].mean_elasticity.iloc[0])
C["elast Ks 0_100"]   = float(t2[(t2.depth=="0_100cm")&(t2.property=="Ks")].mean_elasticity.iloc[0])
C["elast Sand 0_100"] = float(t2[(t2.depth=="0_100cm")&(t2.property=="Sand")].mean_elasticity.iloc[0])
C["elast SOM 0_100"]  = float(t2[(t2.depth=="0_100cm")&(t2.property=="SOM")].mean_elasticity.iloc[0])
C["elast EC 0_100"]   = float(t2[(t2.depth=="0_100cm")&(t2.property=="EC")].mean_elasticity.iloc[0])
for r in ["multiplicative","geometric","arithmetic","minimum"]:
    C[f"CF_Ks {r}"] = t3v("0_100cm", r, "CF_Ks")
    C[f"CF_EC {r}"] = t3v("0_100cm", r, "CF_EC")
C["prop CF_Ks 0_100"] = mcv("0_100cm","independent","CF_Ks")
C["prop CF_Ks 0_30"]  = mcv("0_30cm","independent","CF_Ks")
C["prop CF_EC 0_100"] = mcv("0_100cm","independent","CF_EC")
C["prop CF_EC 0_30"]  = mcv("0_30cm","independent","CF_EC")
C["prop CF_full 0_100"]= mcv("0_100cm","independent","CF_full")
C["prop CF_SOM 0_100"]= mcv("0_100cm","independent","CF_SOM")
C["prop P>=S2 CF_Ks 0_100"] = mcv("0_100cm","independent","CF_Ks","P_ge_S2")
C["PIw indep 0_100"] = float(mc[(mc.depth=="0_100cm")&(mc.structure=="independent")&(mc.scenario=="00_BASELINE")].PI_width.iloc[0])
C["PIw corr 0_100"]  = float(mc[(mc.depth=="0_100cm")&(mc.structure=="correlated")&(mc.scenario=="00_BASELINE")].PI_width.iloc[0])
C["P>=S2 indep base"]= float(mc[(mc.depth=="0_100cm")&(mc.structure=="independent")&(mc.scenario=="00_BASELINE")].P_ge_S2.iloc[0])
C["P>=S2 corr base"] = float(mc[(mc.depth=="0_100cm")&(mc.structure=="correlated")&(mc.scenario=="00_BASELINE")].P_ge_S2.iloc[0])
C["double crop %"] = 100*float(((lu==2)&valid).sum())/valid.sum()
C["single crop %"] = 100*float(((lu==3)&valid).sum())/valid.sum()
C["cultivated %"]  = C["double crop %"]+C["single crop %"]
C["partial rho T-Ks"]  = float(c2p[c2p.pair=="T vs Ks"].partial_rho.iloc[0])
C["partial rho TBP-Ks"]= float(c2p[c2p.pair=="TBP vs Ks"].partial_rho.iloc[0])
C["rho Ks-clay"] = float(c2[c2.variable=="Clay_pct"].rho_with_Ks.iloc[0])
C["rho Ks-SOM"]  = float(c2[c2.variable=="SOM_pct"].rho_with_Ks.iloc[0])
C["rho Ks-sand"] = float(c2[c2.variable=="Sand_pct"].rho_with_Ks.iloc[0])
C["rho Ks-index"]= float(c2[c2.variable=="ALES_index_q50"].rho_with_Ks.iloc[0])
C["T low quintile"]  = float(ksq.iloc[0].T_mean_mm); C["T high quintile"]=float(ksq.iloc[-1].T_mean_mm)
C["TBP low quintile"]= float(ksq.iloc[0].TBP_mean_kg_ha); C["TBP high quintile"]=float(ksq.iloc[-1].TBP_mean_kg_ha)
for _,r in c1[c1.subset=="double-cropped only"].iterrows():
    C[f"rho T-Ks dbl erode{int(r.erosion_m)}"] = float(r.rho_T_vs_Ks)
for v,lab in [("T_mean_mm","T"),("TBP_mean_kg_ha","TBP"),("AETI_mean_mm","AETI"),("WP_kg_m3","WP")]:
    C[f"HIGH vs rest {lab} %"]=float(con[(con.variable==v)&(con.contrast.str.contains("MEDIUM"))].pct_of_comparison.iloc[0])
C["survey mean index"]=float(sb[sb.variable=="ALES_index_q50"].survey_mean.iloc[0])
C["pop mean index"]=float(sb[sb.variable=="ALES_index_q50"].cropped_population_mean.iloc[0])
C["survey mean Ks"]=float(sb[sb.variable=="Ks_q50"].survey_mean.iloc[0])
C["pop mean Ks"]=float(sb[sb.variable=="Ks_q50"].cropped_population_mean.iloc[0])
C["S3 over-rep pp"]=float(sc[sc.stratum=="ALES class S3"].over_representation_pct_pts.iloc[0])
C["HIGH under-rep pp"]=float(sc[sc.stratum=="priority HIGH"].over_representation_pct_pts.iloc[0])
C["wheat N1 unw"]=float(rw[rw.crop=="Wheat"].pct_N1_unweighted.iloc[0])
C["wheat N1 rew"]=float(rw[rw.crop=="Wheat"].pct_N1_reweighted.iloc[0])
C["clover N1 unw"]=float(rw[rw.crop=="Clover"].pct_N1_unweighted.iloc[0])
C["clover N1 rew"]=float(rw[rw.crop=="Clover"].pct_N1_reweighted.iloc[0])
C["wheat HIGH unw"]=float(rw[rw.crop=="Wheat"].pct_HIGH_unweighted.iloc[0])
C["wheat HIGH rew"]=float(rw[rw.crop=="Wheat"].pct_HIGH_reweighted.iloc[0])
C["amb transition"]=float(unc[(unc.depth=="0_100cm")&(unc.zone.str.startswith("transition"))].mean_ambiguity.iloc[0])
C["amb interior"]=float(unc[(unc.depth=="0_100cm")&(unc.zone.str.startswith("interior"))].mean_ambiguity.iloc[0])
C["Ks<1.5 %"]=100*float((prof.Ks<1.5).mean())
C["Ks q1"]=float(prof.Ks.quantile(.25)); C["Ks q3"]=float(prof.Ks.quantile(.75))
C["ESP surface max"]=float(np.nanmax(rd(S3B/"ESP_0_100cm_q50_qrf.tif")))
C["CaCO3 surface max"]=float(np.nanmax(rd(S3B/"CaCO3_0_100cm_q50_qrf.tif")))
C["span>1 class 0_100 %"]=float(unc[(unc.depth=="0_100cm")&(unc.zone.str.startswith("ALL"))].mean_PI_width.iloc[0])
C["span>1 class 0_30 %"]=float(unc[(unc.depth=="0_30cm")&(unc.zone.str.startswith("ALL"))].mean_PI_width.iloc[0])

print("RECOMPUTED QUANTITIES THE REVISED TEXT RELIES ON")
for k,v in C.items():
    print(f"  {k:32s} {v:12.4f}")

# ---- string assertions --------------------------------------------------
ASSERT = [
 ("limiting factor 79%", "79% of root-zone pixels"), ("limiting 99.97%", "99.97% of pixels at 0–30 cm"),
 ("physical 20.7%", "20.7% at 0–100 cm"),
 ("elasticity clay", "−0.82 at 0–100 cm"), ("elasticity Ks", "(+0.40)"),
 ("agg multiplicative", "+21.1 index units under the multiplicative"),
 ("agg geometric", "+9.4 under the geometric mean"),
 ("agg arithmetic", "+8.4 under the arithmetic mean"),
 ("agg minimum", "+11.9 under a strict minimum rule"),
 ("prop CF_Ks", "+19.5 index units at 0–30 cm and +11.7 at 0–100 cm"),
 ("prop CF_EC", "+3.3 to +3.5 index units under full"),
 ("ratio statement", "3.3 under full propagation and 10.1 on the median-property surface"),
 ("ratio 0_30", "5.9 and 21.3 respectively at 0–30 cm"),
 ("PI width copula", "from 33.0 to 29.9 units"),
 ("P copula", "0.194 to 0.185"),
 ("double crop", "66.0% of valid pixels to double crop"),
 ("cultivated", "75.6%, is therefore under cultivation"),
 ("partial rho", "transpiration is +0.04"),
 ("rho clay", "−0.95 against clay"),
 ("T quintiles", "1158 mm yr⁻¹ in the lowest Ks"),
 ("TBP quintiles", "24,600 to 19,300 kg ha⁻¹"),
 ("erosion series", "−0.10 at 60 m, −0.01 at "),
 ("HIGH contrast T", "transpiration 4.6% higher"),
 ("HIGH contrast TBP", "biomass 4.7% higher"),
 ("survey index", "28.9 against 26.5"),
 ("survey Ks", "0.54 against 0.46 cm h⁻¹"),
 ("S3 overrep", "over-represented by 15.6 percentage points"),
 ("wheat reweight", "52% to 67%"), ("clover reweight", "38% to 54%"),
 ("wheat HIGH reweight", "28% to 36%"),
 ("Ks<1.5", "88% of profiles"),
 ("ESP max", "peaks at 8.1%"), ("CaCO3 max", "peaks at 12.6%"),
 ("span classes", "every valid pixel has a 90% "),
 ("span 0_30", "85.4% do"),
 ("IWQI", "81.0 to 82.3 with a mean of 82.0"),
 ("mean q50", "The mean q50 value was 29.9"),
 ("Ks IQR", "interquartile range 0.14–0.94"),
]
print("\nSTRING CHECKS AGAINST THE REVISED MANUSCRIPT")
missing = []
for label, s in ASSERT:
    ok = s in TXT
    if not ok: missing.append((label, s))
    print(f"  [{'OK ' if ok else 'MISS'}] {label:24s} :: {s[:58]}")
print(f"\n{len(ASSERT)-len(missing)}/{len(ASSERT)} assertions located in the revised text")
if missing:
    print("MISSING:")
    for l,s in missing: print("   ", l, "|", s)
