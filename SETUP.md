# How to run this analysis on your own computer

This is the **complete, no-experience-assumed** guide. If you have never used
Python, a terminal, or Git, that is fine — every command you need is written out
below and you can copy and paste it.

Read the whole of **Step 0** before you start. It saves hours.

---

## Step 0 — What you must have before you begin

This repository contains the **code only**. It does not contain the field data,
and it does not contain the satellite imagery. This is deliberate: the soil and
crop data belong to a manuscript that is still under review, and the satellite
stack is a 3 GB file that GitHub will not accept.

So the code alone **cannot run**. You need two things from the corresponding
author (Dr Ahmed El-Naggar, `a.elnaggar@un-ihe.org`):

**Bundle A — the field data** (small, a few MB)

```
Analysis/
    Analysis_Banisuef_New_FF__March_2025.xls
    water_analyses Benisueif_FFF.xls
Layers/
    Study_Area_Last.shp   (plus its .shx .dbf .prj companions)
    Soil_Profiles.shp     (plus its .shx .dbf .prj companions)
    Water_Samples.shp     (plus its .shx .dbf .prj companions)
Crop_March2025_11/
    Crop_March2025_11.shp (plus its .shx .dbf .prj companions)
```

> A "shapefile" is not one file — it is a set of 4 to 6 files that share a name
> and differ only in extension. They must all travel together or nothing works.

**Bundle B — the pre-computed satellite stack** (large, ~1–3 GB)

```
outputs/stage2/covariates_at_profiles.csv
outputs/stage2/covariates_at_water_samples.csv
outputs/stage2b_local_stack/covariate_stack_30m.tif
```

Bundle B is what Stage 2 would have produced. Stage 2 talks to Google Earth
Engine, needs a Google Cloud service-account key that cannot be shared, and takes
hours. **Accepting Bundle B lets you skip Stage 2 entirely.** Unless you
specifically want to rebuild the satellite covariates yourself, take Bundle B.

You do **not** need a Google account, an Earth Engine account, or any API key to
reproduce every result in the manuscript from Bundle A + Bundle B.

---

## Step 1 — Install the two programs you need

You need **Git** (to download the code) and **Miniforge** (to run Python).

### Windows

1. Git — go to <https://git-scm.com/download/win>. The download starts on its
   own. Run the installer and click **Next** on every screen. The defaults are
   correct.
2. Miniforge — go to <https://conda-forge.org/download/> and download the
   Windows x86_64 installer. Run it. When it asks, choose **"Just Me"**. Leave
   everything else at the default.

After both finish, open the **Start menu**, type `Miniforge Prompt`, and open it.
A black window appears with a line ending in `(base) C:\Users\yourname>`.
**Every command in this guide is typed into that window.** Not into PowerShell,
not into Command Prompt — into the Miniforge Prompt.

### macOS / Linux

```bash
# Git is usually already installed. Check with:
git --version

# Install Miniforge:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Close and reopen your terminal afterwards. You should see `(base)` at the start
of the prompt line.

---

## Step 2 — Download the code

In the Miniforge Prompt, type these lines one at a time, pressing Enter after
each:

```
cd %USERPROFILE%\Documents
git clone https://github.com/El-Naggar8876/bani-suef-soil-suitability.git
cd bani-suef-soil-suitability
```

On macOS/Linux the first line is `cd ~/Documents` instead.

You now have a folder called `bani-suef-soil-suitability` inside Documents. This
folder is called **the project folder** for the rest of this guide. Everything
happens inside it.

---

## Step 3 — Build the Python environment

Still in the same window, still inside the project folder:

```
conda env create -f environment.yml
conda activate bani-suef-suitability
```

The first command downloads roughly 1.5 GB of scientific libraries and takes
**10 to 30 minutes**. It will look frozen at "Solving environment". It is not
frozen. Let it work.

When it is done, your prompt changes from `(base)` to
`(bani-suef-suitability)`. That prefix is how you know the environment is
active.

> **Important:** every time you close the window and come back later, you must
> run `conda activate bani-suef-suitability` again before running any script.
> If you forget, you will get `ModuleNotFoundError` and think something is
> broken. Nothing is broken — you just forgot this line.

<details>
<summary>Alternative: pip instead of conda (not recommended on Windows)</summary>

```
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # macOS/Linux
pip install -r requirements.txt
```

On Windows, `geopandas`, `rasterio` and `fiona` frequently fail to install this
way because they need compiled geospatial libraries. Use conda unless you know
you have GDAL working.
</details>

---

## Step 4 — Put the data in place

Unzip **Bundle A** and **Bundle B** so that the project folder looks exactly
like this:

```
bani-suef-soil-suitability/
├── Analysis/                    <- from Bundle A
├── Layers/                      <- from Bundle A
├── Crop_March2025_11/           <- from Bundle A
├── data/                        <- already there from GitHub
├── notebooks/                   <- already there from GitHub
├── outputs/
│   ├── stage2/                  <- from Bundle B
│   └── stage2b_local_stack/     <- from Bundle B
├── environment.yml
└── README.md
```

Two mistakes people make here:

- **Double folders.** Unzipping sometimes creates
  `Analysis/Analysis/the-file.xls`. The `.xls` files must sit *directly* inside
  `Analysis/`, with nothing in between.
- **Wrong level.** The folders go beside `notebooks/`, not inside it.

These folders are listed in `.gitignore`, so Git will correctly ignore them and
you cannot accidentally publish the unreleased data.

---

## Step 5 — Check everything before you run anything

This repository includes a checker that inspects your setup and tells you in
plain English what is missing. Run it now:

```
python notebooks/check_setup.py
```

It prints a list with `OK` or `MISSING` beside every required file. **Do not
move on until every item under "Bundle A" and "Bundle B" says OK.** Fixing a
missing file now takes one minute; discovering it three hours into Stage 5 does
not.

---

## Step 6 — Run the analysis

Move into the scripts folder and run the stages **in this order**. Each one
prints its progress and writes its results into `outputs/`.

```
cd notebooks

python stage1_data_audit.py
python stage3a_qrf_cv.py
python stage3b_predict_maps.py
python stage4_iwqi_surface.py
python stage5_ales_montecarlo.py
python stage6_counterfactual.py
python stage7_cropping_mismatch.py
python stage7b_crop_validation.py
python stage8_figures.py
python stage9_graphical_abstract.py
```

Notice that **Stage 2 is not in the list** — Bundle B replaced it. If you did
decide to rebuild the satellite stack yourself, see the appendix at the bottom.

Expect **2 to 4 hours** in total on a machine with 16 GB of RAM. Stages 5 and 6
account for most of it because they each run hundreds of Monte Carlo
simulations. The window will look idle for long stretches. That is normal.

You can stop after any stage and resume later — each stage reads the previous
stage's saved files from disk, so nothing is lost when you close the window.
Just remember `conda activate bani-suef-suitability` when you come back.

### What each stage does

| Stage | Script | What it produces |
|-------|--------|------------------|
| 1 | `stage1_data_audit.py` | Cleans the 60 soil profiles and 20 water samples; computes the irrigation water quality index |
| 3a | `stage3a_qrf_cv.py` | Compares four mapping methods using spatial cross-validation |
| 3b | `stage3b_predict_maps.py` | Maps 8 soil properties across the whole area, with uncertainty |
| 4 | `stage4_iwqi_surface.py` | Turns 20 water samples into a continuous water-quality surface |
| 5 | `stage5_ales_montecarlo.py` | 200 simulations of the land-suitability rating — the core result |
| 6 | `stage6_counterfactual.py` | Asks "what if we fixed drainage / salinity / …?" for 7 scenarios |
| 7 | `stage7_cropping_mismatch.py` | Compares suitability against what is actually being farmed |
| 7b | `stage7b_crop_validation.py` | Validates against the 993-point March 2025 field survey |
| 8 | `stage8_figures.py` | Manuscript Figures 1–6 |
| 9 | `stage9_graphical_abstract.py` | The graphical abstract |

---

## Step 7 — Look at the results

Everything lands in the `outputs/` folder inside the project folder.

- **Figures** — `outputs/figures/`. These open in any image viewer or PDF reader.
- **Tables** — the `.csv` files scattered through `outputs/stage1/`,
  `outputs/stage3a/`, `outputs/stage6/`, `outputs/stage7b/`. These open in Excel.
- **Maps** — the `.tif` files. These are geospatial rasters, not ordinary
  pictures. Open them with [QGIS](https://qgis.org) (free) rather than an image
  viewer, or they will look like meaningless grey squares.

To confirm you reproduced the published numbers, check these three:

| Published value | Exact file | Where in it |
|-----------------|------------|-------------|
| PICP₉₀ = 0.85 | `outputs/stage3a/cv_metrics.csv` | mean of the `PICP90` column for `model == QRF` |
| User's accuracy = 0.988 | `outputs/stage7b/crop_kpis.csv` | the accuracy row |
| Drainage is the binding constraint | `outputs/stage6/gap_decomposition_0_100cm.csv` | the `CF_Ks` scenario has the largest mean Δ |

All random seeds are fixed at `RANDOM_STATE = 42`, so your numbers should match
the manuscript exactly, not just approximately. If they differ, something in the
input data differs — say so before assuming the code is wrong.

---

## When something goes wrong

| What you see | What it means | What to do |
|---|---|---|
| `ModuleNotFoundError: No module named 'geopandas'` | The environment is not active | Run `conda activate bani-suef-suitability` |
| `FileNotFoundError: ...Analysis_Banisuef_New_FF__March_2025.xls` | Bundle A is missing or in the wrong place | Re-read Step 4, then run `check_setup.py` |
| `FileNotFoundError: ...covariate_stack_30m.tif` | Bundle B is missing | Ask for Bundle B; it is too big for email — expect a cloud link |
| `FileNotFoundError: ...covariates_at_profiles.csv` | You have Bundle B's `.tif` but not its two `.csv` files | Both parts of Bundle B are needed |
| `DriverError` or `.shx` errors | A shapefile arrived incomplete | Ask for the folder again as a zip, not as individual files |
| `MemoryError`, or the computer freezes | Stages 5–6 need real RAM | Close other applications; 16 GB is the practical minimum |
| `'python' is not recognized` | Wrong window | Use the Miniforge Prompt, not PowerShell |

If you are stuck, open an issue at
<https://github.com/El-Naggar8876/bani-suef-soil-suitability/issues> and paste
**the last 20 lines** of the error. The bottom line of a Python error is the one
that matters, and it is the line people most often leave out.

---

## Appendix A — Rebuilding the satellite stack yourself (Stage 2)

Skip this unless you have a reason to distrust Bundle B.

Stages 2 and 2b download Sentinel-1, Sentinel-2, SRTM, CHIRPS, TerraClimate and
ERA5-Land through Google Earth Engine. They authenticate with a **Google Cloud
service-account key**, not with the ordinary `earthengine authenticate` browser
login. You must create your own key:

1. Register for Earth Engine at <https://earthengine.google.com/signup/> and
   wait for approval (typically 1–2 days).
2. In the Google Cloud console, create a project, enable the Earth Engine API,
   create a **service account**, and download its **JSON key**.
3. Register that service account at
   <https://signup.earthengine.google.com/#!/service_accounts>.
4. Save the JSON file in the project folder as:

   ```
   .secrets/gee_service_account.json
   ```

   The `.secrets/` folder is git-ignored. **Never commit this file, and never
   send it to anyone** — it grants access to your Google Cloud billing account.

Then run, before Stage 3a:

```
python stage2_covariate_stack.py     # writes outputs/stage2/*.csv
python stage2b_download_stack.py     # writes outputs/stage2b_local_stack/*.tif
```

Both are required. They are not alternatives to each other:
`stage2` produces the covariate values sampled at the 60 profile locations
(needed by Stages 3a and 3b), while `stage2b` produces the wall-to-wall raster
stack (needed by Stages 3b, 4, 5, 6, 7 and 8). Stage 2 also starts an
asynchronous export to Google Drive, which is a leftover from development and is
not used by any later stage — you can ignore it.

Combined runtime depends on the Earth Engine queue and is typically several
hours.

---

## Appendix B — For the corresponding author: preparing the bundles

Run these from wherever your working copy with the real data lives.

**Bundle A** — zip these three folders together, preserving folder names:

```
Analysis/  Layers/  Crop_March2025_11/
```

Check before sending: each `.shp` must be accompanied by its `.shx`, `.dbf` and
`.prj`. Zipping the folders (rather than selecting files) handles this
automatically. Typical size: a few MB, so email or OneDrive both work.

**Bundle B** — zip these three files, preserving the `outputs/` folder
structure:

```
outputs/stage2/covariates_at_profiles.csv
outputs/stage2/covariates_at_water_samples.csv
outputs/stage2b_local_stack/covariate_stack_30m.tif
```

Typical size: 1–3 GB. Too large for email — use OneDrive, Google Drive or
WeTransfer.

**Never include** `.secrets/gee_service_account.json` in either bundle. It is a
credential, not data.

At manuscript acceptance, both bundles should be deposited on Zenodo under the
companion dataset DOI, at which point this whole appendix can be replaced by a
download link and the analysis becomes reproducible by anyone.
