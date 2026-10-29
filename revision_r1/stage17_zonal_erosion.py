"""Stage 17 - are the zonal WaPOR contrasts the same edge artefact as the gradient?

Stage 13 reports a zonal comparison between the intervention-priority zones and
the rest of the cropped land, and stage 13b shows that the raw correlation
between predicted Ks and transpiration is produced by 100 m WaPOR cells that
straddle the cultivated margin, by eroding the cropped mask inward until the
correlation disappears.

The zonal contrasts were never put through that same correction, and they are
exposed to it: cropped S3 land lies disproportionately along the margin, since
almost all bare-desert pixels classify as S3. This recomputes every zonal
contrast under the same erosion widths stage 13b uses, so that the raw and the
corrected figures can be compared directly. The result is Table S8 of the
revised manuscript.
"""
from __future__ import annotations
import sys, pathlib
import numpy as np, pandas as pd, rasterio
from scipy import ndimage

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from revpaths import REPO, OUT, WAPOR

S5, S7 = REPO / "outputs" / "stage5", REPO / "outputs" / "stage7"
S3B = REPO / "outputs" / "stage3b"
DEPTH = "0_100cm"
BLOCK = 33                       # 1 km blocks on the 30 m grid
NBOOT = 200
RNG = np.random.default_rng(20261028)


def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            a[a == s.nodata] = np.nan
        return a


def block_ci(a_vals, a_blk, b_vals, b_blk, nboot=NBOOT):
    """Moving-block bootstrap CI on mean(a) - mean(b), resampling 1 km blocks."""
    ga = {b: np.where(a_blk == b)[0] for b in np.unique(a_blk)}
    gb = {b: np.where(b_blk == b)[0] for b in np.unique(b_blk)}
    ka, kb = np.array(list(ga)), np.array(list(gb))
    out = np.empty(nboot)
    for i in range(nboot):
        sa = np.concatenate([ga[b] for b in RNG.choice(ka, ka.size, replace=True)])
        sb = np.concatenate([gb[b] for b in RNG.choice(kb, kb.size, replace=True)])
        out[i] = a_vals[sa].mean() - b_vals[sb].mean()
    return np.percentile(out, (2.5, 97.5))


def main():
    T = rd(WAPOR / "wapor_T_mean_2021_2025_30m.tif")
    AETI = rd(WAPOR / "wapor_AETI_mean_2021_2025_30m.tif")
    TBP = rd(WAPOR / "wapor_TBP_mean_2021_2025_30m.tif")
    pri = rd(S7 / f"policy_priority_{DEPTH}.tif")
    lu = rd(S7 / "landuse_class.tif")
    cls = rd(S5 / f"ales_class_q50_{DEPTH}.tif")
    idx = rd(S5 / f"ales_index_q50_{DEPTH}.tif")

    H, W = pri.shape
    rr, cc = np.mgrid[0:H, 0:W]
    blk = ((rr // BLOCK) * 10000 + (cc // BLOCK)).ravel()

    fin = np.isfinite(idx) & np.isfinite(T) & np.isfinite(AETI) & np.isfinite(TBP)
    cropped, dbl = np.isin(lu, [2, 3]), (lu == 2)
    lays = {"T": T, "AETI": AETI, "TBP": TBP}

    rows = []
    for sub_name, m0 in [("all cropped", cropped), ("double-cropped", dbl)]:
        for e_px in (0, 2, 4, 7):
            m = ndimage.binary_erosion(m0, np.ones((2 * e_px + 1,) * 2)) if e_px else m0
            m = m & fin
            if m.sum() < 2000:
                continue
            groups = {"HIGH": m & (pri == 1),
                      "MED+WATCH": m & np.isin(pri, [2, 3]),
                      "S3": m & (cls == 3),
                      "N1": m & (cls == 4)}
            for a, b in [("HIGH", "MED+WATCH"), ("HIGH", "S3"), ("N1", "S3")]:
                ma, mb = groups[a], groups[b]
                if ma.sum() < 500 or mb.sum() < 500:
                    continue
                r = dict(subset=sub_name, erosion_m=e_px * 30,
                         n_a=int(ma.sum()), n_b=int(mb.sum()),
                         contrast=f"{a} - {b}")
                for vn, v in lays.items():
                    av, bv = v[ma], v[mb]
                    r[f"{vn}_pct"] = round(100 * (av.mean() - bv.mean()) / bv.mean(), 1)
                    if vn == "T":
                        lo, hi = block_ci(av, blk[ma.ravel()], bv, blk[mb.ravel()])
                        r["T_pct_lo"] = round(100 * lo / bv.mean(), 1)
                        r["T_pct_hi"] = round(100 * hi / bv.mean(), 1)
                rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "wapor_zonal_eroded.csv", index=False)
    print("===== zonal WaPOR contrasts under progressive erosion =====")
    print(df.to_string(index=False))
    print("\nwrote", OUT / "wapor_zonal_eroded.csv")


if __name__ == "__main__":
    main()
