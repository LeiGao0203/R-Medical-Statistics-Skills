# NHANES 2017–2018 LPA example

This example applies a reproducible latent profile analysis to public NHANES 2017–2018 data. It uses six continuous cardiometabolic indicators: BMI, waist circumference, mean systolic blood pressure, HbA1c, total cholesterol, and reverse-coded HDL-C.

## Reproduce

From the repository root:

```bash
Rscript advanced-stats/latent-profile-analysis/scripts/download_nhanes_2017_2018.R
Rscript example/nhanes-lpa/analysis/01_nhanes_lpa.R
```

Required R packages are `haven`, `dplyr`, `tidyr`, `ggplot2`, `readr`, `tibble`, and `mclust`.

## Directory structure

- `data/raw/`: downloaded CDC XPT files; these are public-use source files and should be retained with their official documentation and terms.
- `data/derived/`: merged, adult, complete-case analysis data.
- `analysis/`: executable R analysis script.
- `results/`: model-selection table, posterior classification, class sizes, original-scale profile means, standardization parameters, and session information.
- `figures/`: standardized profile plot.
- `metadata/`: data URLs and variable mapping.

## Interpretation boundary

The example fits an unweighted `mclust` model. NHANES is a complex survey, so this is a sample-internal exploratory LPA, not a national population estimate. `WTMEC2YR`, `SDMVPSU`, and `SDMVSTRA` are retained for a future survey-aware sensitivity analysis. The resulting profiles are empirical patterns, not validated disease subtypes or causal exposures.

## Data source

CDC/NCHS, [NHANES 2017–2018 data portal](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017). Review the [NHANES overview](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overview.aspx?BeginYear=2017), [laboratory overview](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overviewlab.aspx?BeginYear=2017), and analytic guidelines before using the data for publication.
