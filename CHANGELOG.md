# Changelog

All notable changes to this repository are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-10-27

Stages added for the first revision of the Geoderma Regional manuscript
(GEODRS-D-26-00837). All live in `revision_r1/` and write to
`outputs/stage_rev/`; they locate the repository root themselves, so they run
from either the repository or the manuscript revision folder.

### Added
- `revision_r1/stage10_sensitivity.py` - structural sensitivity of the ALES-Arid
  index: limiting-factor census, normalised elasticities, the decomposition
  repeated under four aggregation rules, and perturbation of the r_wetness
  lower optimum. Establishes that the ranking of constraints does not depend on
  the multiplicative aggregation rule.
- `revision_r1/stage11_copula_mc.py` - Monte Carlo propagation with a Gaussian
  copula reproducing the measured rank correlations among the eight soil
  properties, compared against the independent draws used in v1.0.0.
- `revision_r1/stage12_spatial_uncertainty.py` - classification-ambiguity
  surface and its drivers.
- `revision_r1/stage13_wapor_check.py` and `stage13b_wapor_robustness.py` -
  independent comparison against FAO WaPOR v3 transpiration, evapotranspiration
  and biomass, with mixed-pixel, confounding and crop-composition controls.
- `revision_r1/stage14_survey_bias.py` - quantification and post-stratified
  correction of the corridor-sampling bias in the 993-point crop survey.
- `revision_r1/stage15_cv_sampling.py` - full cross-validation metrics and
  sampling-density statistics.
- `revision_r1/wapor_to_gtiff.py` - WaPOR windows reprojected onto the 30 m
  analysis grid; outputs in `outputs/wapor/`.
- `revision_r1/verify_manuscript_numbers.py` and `verify_revision_numbers.py` -
  recompute every quantitative claim in the submitted and revised manuscripts
  against the pipeline outputs.

### Fixed
- `CITATION.cff` and `.zenodo.json` carried placeholders for the repository URL,
  the DOI and the author list. All five authors are now listed with
  affiliations; ORCIDs remain to be supplied.

### Note on the published description
The sub-rating curves described in Section 2.7 of the submitted manuscript did
not match those implemented in `notebooks/stage5_ales_montecarlo.py`. The code
was correct and the published results reproduce exactly from it; the manuscript
description has been corrected to match, and the revised Table 1 states every
threshold as implemented.


## [1.0.0] — 2026-05-11

### Added
- Initial public release of the seven-stage analysis pipeline accompanying the manuscript submission to *Geoderma Regional*.
- Stages 1–7 covering data audit, covariate-stack assembly, QRF spatial cross-validation and prediction, IWQI kriging, Monte Carlo Sys–Verheye ALES-Arid, counterfactual gap decomposition, cropping confrontation, and 993-point ground-truth validation.
- Auxiliary stages 8 and 9 for publication figures and the graphical abstract.
- `requirements.txt` and `environment.yml` for pip and conda installation.
- `CITATION.cff` and `.zenodo.json` for machine-readable metadata.
- MIT licence.
