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

**Bundle B — the satellite covariates.** Here you have a choice of two routes.

```
outputs/stage2/covariates_at_profiles.csv
outputs/stage2/covariates_at_water_samples.csv
outputs/stage2b_local_stack/covariate_stack_30m.tif
```

These three files are what Stages 2 and 2b produce. You can either **receive**
them or **rebuild** them:

| | Route 1 — receive the archived stack | Route 2 — rebuild it yourself |
|---|---|---|
| **You get** | Byte-for-byte the inputs behind the published figures | Your own freshly downloaded satellite composites |
| **Answers** | "Do I reproduce the published numbers *exactly*?" | "Does this method work *independently*?" |
| **You need** | A 1–3 GB download link | An Earth Engine account and a Google Cloud project |
| **Setup time** | Minutes | 1–2 days waiting for Google's approval |
| **Run time** | None | Several hours |
| **Do this if** | You are a coauthor, reviewer, or checking the results | You are extending the work, or the archive is unavailable |

**If you are a coauthor or reviewer, take Route 1.** It is faster and it is the
only route that verifies the published numbers. Ask the corresponding author for
the three files above.

For Route 2, see **Appendix A**. You will edit one file, `config.py`, to add your
own Google Cloud project ID.

> **An honest caveat about Route 2.** Google periodically reprocesses the
> Sentinel archives. A stack rebuilt in 2027 will be very close to, but not
> identical with, the one used for the manuscript, so your final numbers will
> land near the published ones rather than exactly on them. That is a property
> of satellite archives, not a bug in this code, and no amount of seed-fixing
> can remove it. Route 1 exists precisely because of this.

With Route 1, you need **no** Google account, Earth Engine account, or API key
to reproduce every result in the manuscript.

---

## Step 1 — Install the two programs you need

You need **Git** (to download the code) and **conda** (to run Python).

> ### ⚠️ Read this before installing anything
>
> **Commands in this guide are labelled Windows or macOS/Linux. Run only the
> ones for your machine.** Copying a macOS command into a Windows window
> produces confusing errors like
> `curl: (3) URL using bad/illegal format or missing URL` — that is Windows not
> understanding Mac syntax, not a broken download link.
>
> **You may already have both programs.** Check first, and skip whatever you
> already have.

### Step 1a — Check what you already have

Open a terminal:

- **Windows** — Start menu → type `Miniforge Prompt` and open it. If there is no
  Miniforge Prompt, try `Anaconda Prompt`. If neither exists, open
  `Command Prompt`.
- **macOS** — open `Terminal` from Applications → Utilities.
- **Linux** — open your usual terminal.

Type these two lines, pressing Enter after each:

```
git --version
conda --version
```

Now read the results:

| What you see | What it means |
|---|---|
| Both print a version number | **You are done with Step 1. Go to Step 2.** |
| `conda` prints a version | conda is installed — skip the Miniforge install below |
| `git` prints a version | Git is installed — skip the Git install below |
| `'...' is not recognized` / `command not found` | That one is missing — install it below |

There is one more shortcut. If the start of your prompt line already reads
`(base)`, like this:

```
(base) C:\Users\yourname>
```

then conda is installed **and** already switched on. That is all Step 1 was ever
trying to achieve.

### Step 1b — Install only what was missing

<details open>
<summary><b>Windows</b></summary>

**Git** — go to <https://git-scm.com/download/win>. The download starts on its
own. Run the installer and click **Next** on every screen; the defaults are
correct.

**conda** — go to <https://conda-forge.org/download/> and download the
**Windows x86_64** installer (a `.exe` file). Run it. When it asks, choose
**"Just Me"**. Leave everything else at the default.

Afterwards, open the **Start menu**, type `Miniforge Prompt`, and open that. Your
prompt should end with `(base) C:\Users\yourname>`.

Anaconda works just as well as Miniforge — if you already have it, use its
`Anaconda Prompt` and do not install Miniforge on top.

</details>

<details>
<summary><b>macOS / Linux</b> — do not run these on Windows</summary>

Git usually comes pre-installed. To install conda:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Close and reopen your terminal afterwards. You should see `(base)` at the start
of the prompt line.

</details>

Whichever terminal shows `(base)` is the window you will use for the whole rest
of this guide. On Windows that is normally the Miniforge Prompt or Anaconda
Prompt — not PowerShell.

---

## Step 2 — Download the code

> **Already have the folder?** If you were sent the project folder directly, or
> you already cloned it, skip the `git clone` line. Just `cd` into the folder you
> have and go to Step 3.

Type these lines one at a time, pressing Enter after each.

**Windows:**

```
cd %USERPROFILE%\Documents
git clone https://github.com/El-Naggar8876/bani-suef-soil-suitability.git
cd bani-suef-soil-suitability
```

**macOS / Linux:**

```
cd ~/Documents
git clone https://github.com/El-Naggar8876/bani-suef-soil-suitability.git
cd bani-suef-soil-suitability
```

The only difference is the first line. Everything after Step 2 is identical on
all three systems.

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

Notice that **Stage 2 is not in the list.** On Route 1 the archived stack
replaced it. On Route 2 you already ran Stages 2 and 2b during Appendix A, so by
this point they are done either way.

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

## Appendix A — Route 2: rebuilding the satellite stack yourself

This is Route 2 from Step 0. Take it if you are extending the work, or if the
archived stack is unavailable. Skip it if you were given Bundle B.

Stages 2 and 2b download Sentinel-1, Sentinel-2, SRTM, CHIRPS, TerraClimate and
ERA5-Land through Google Earth Engine, then assemble them into the 53-band 30 m
stack. Everything specific to *you* lives in one file, `config.py`, in the
project folder. **That is the only file you edit.**

### A1 — Get an Earth Engine account

1. Register at <https://earthengine.google.com/signup/>. Approval typically
   takes 1–2 days, so start here.
2. In the Google Cloud console <https://console.cloud.google.com/>, create a
   project (or use an existing one) and enable the **Earth Engine API** on it.
3. Note the project **ID**. It is shown in the project selector at the top of
   the console and looks like `ee-yourname` or `my-project-464812`. The ID is
   often *not* the same as the display name — copy the ID.

### A2 — Tell this repository who you are

Open `config.py` in the project folder with any text editor — Notepad is fine.
Near the top you will see:

```python
GEE_PROJECT: str | None = None      # e.g. "ee-yourname" or "my-project-464812"
```

Change it to your project ID, keeping the quotation marks:

```python
GEE_PROJECT: str | None = "ee-yourname"
```

Save the file. That is the only change most people need to make. The rest of
`config.py` is documentation, plus a section marked **FROZEN** listing the
values that must stay as they are for the run to count as a reproduction — the
March 2025 date lock, the 30 m grid, EPSG:32636, and the random seed.

### A3 — Log in

Run this once. It opens a browser window; approve the access request.

```
earthengine authenticate
```

Earth Engine remembers you on this machine afterwards.

### A4 — Test before committing to a long download

```
python config.py
```

This does a one-second round-trip to Earth Engine and prints either `OK` or a
specific explanation of what is wrong. Run it now rather than discovering a
credential problem forty minutes into a download.

### A5 — Run both stages

```
cd notebooks
python stage2_covariate_stack.py     # writes outputs/stage2/*.csv
python stage2b_download_stack.py     # writes outputs/stage2b_local_stack/*.tif
```

**Both are required — they are not alternatives.** `stage2` produces the
covariate values sampled at the 60 profile locations, which Stages 3a and 3b
need. `stage2b` produces the wall-to-wall raster, which Stages 3b, 4, 5, 6, 7
and 8 need. Running only one leaves later stages unable to start.

Stage 2 also kicks off an asynchronous export to Google Drive. That is a
leftover from development, no later stage reads it, and you can ignore it.

Combined runtime depends on the Earth Engine queue and is typically several
hours. Afterwards, continue from Stage 3a in Step 6.

### Using a service-account key instead

The original authors ran these stages unattended with a Google Cloud
service-account key rather than a browser login. If you have such a key, place
it at `.secrets/gee_service_account.json` and `config.py` will find and prefer
it automatically — no project ID needed, since the key carries its own.

**A service-account key is a credential, not data.** It is git-ignored for a
reason. Never commit it, never email it, and never include it in a Zenodo
archive. Anyone holding that file can spend money on your Google Cloud account.

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
