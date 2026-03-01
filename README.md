# LCM-Spatial-MR (Oregon Marine Reserves)
Factor analysis + spatial preference patterns (Moran’s I / GAM) + latent-class covariates.

This repository contains analysis code and supporting files for my research on **preference segmentation for Oregon Marine Reserves**, including:
- **Factor analysis** (attitudinal constructs)
- **Spatial diagnostics** (e.g., Moran’s I) and **spatial GAM surfaces**
- **Latent class model (LCM) covariates / class-membership specifications** (Latent GOLD inputs)

> **Data are not included in this GitHub repository.** See **Data access** below.
>
> ## Repository structure

```text
LCM-Spatial-MR/
  README.md
  .gitignore
  .gitattributes
  code/
    Factor_Analysis/
      Factor_analysis2.R
    Spatial_Analysis/
      spatial_Moran_GAM_bi_fin1.R
      spatial_Moran_GAM_classmemb_fin3.R
    LCM_Covariates/
      Best_cov_syntax.txt
      Best_cov.lgf
  paper/
    Chap_3_24_paper.pdf
  data/        (NOT included)
  output/      (NOT included)


---

## Requirements

Latent GOLD

R (recommended: R >= 4.2) and packages commonly used in the scripts, e.g.:
- `sf`, `dplyr`, `tidyr`, `ggplot2`, `spdep`, `mgcv`, `readr`, `stringr`
---

## How to run (high level)

1. Place required input files in `data/` (not provided here).
2. Run scripts from the `code/` folder in this order (typical workflow):
   - `code/Factor_Analysis/Factor_analysis2.R`
   - `code/Spatial_Analysis/spatial_Moran_GAM_bi_fin!.R`
   - `code/Spatial_Analysis/spatial_Moran_GAM_classmemb_fin3.R`
   - LCM covariate model: `code/LCM_Covariates/`

Outputs should be saved to `output/` (not tracked).

---

## Data access

- If the dataset is restricted/non-public: email me to request access.

Contact: **mpignatari@usu.edu** (or **mpignatari@hotmail.com**)

---

## Paper / preprint

A PDF version is provided in `paper/`:
- `paper/Chap_3_24_paper.pdf`

---

## Citation

If you use this repository, please cite:

Pignatari, M. (2026). GitHub repository:  
https://github.com/mpignatari/LCM-Spatial-MR

---

## License

No license specified yet. If you want others to reuse code easily, consider adding an MIT License.



