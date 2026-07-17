# R Medical Statistics Skills

<p align="center">
  <img src="assets/github-promo.png" alt="R Medical Statistics Skills promotional banner">
</p>

[中文版](README.md)

An AI coding agent skills collection for medical statistics and R-based statistical analysis. The project uses a generic `SKILL.md` structure that can be used by Codex and by other coding agents that support skills, rules, or knowledge-base directories. You can also send the relevant directories to an agent and let it install or read them for its own runtime. The project includes basic statistics, advanced statistics, commonly used methods in medical literature, and helper skills for creating plain R scripts, statistical reports, and Jupyter notebooks.

## Contents

- `basic-stats/`: Basic medical statistics skills, including t-tests, ANOVA, chi-square tests, correlation analysis, ROC analysis, sample size estimation, and statistical plotting.
- `advanced-stats/`: Advanced statistics skills, including ANCOVA, multiple regression, logistic regression, survival analysis, PCA, structural equation modeling, multilevel models, and latent profile analysis.
- `literature-stats/`: Skills for methods commonly seen in medical literature, including propensity score methods, Fine-Gray models, restricted cubic splines, subgroup analysis, and trend tests.
- `r-script/`: Original R script skill for generating reproducible `.R` analysis scripts for RStudio, command-line R, and non-notebook users.
- `quarto-report/`: Original Quarto/R Markdown report skill for generating medical statistics reports exportable to HTML, Word, or PDF.
- `jupyter-notebook/`: Original Jupyter Notebook skill for creating, organizing, and validating reproducible notebooks.
- `example/`: Reproducible examples based on public datasets, including demo data, notebooks, R scripts, statistical tables, and generated figures.

## Examples

The `example/` directory provides one complete representative medical statistics case study showing how these skills turn public data into a reproducible R analysis workflow. The example includes raw data in `data/`, analysis scripts and result tables in `analysis/`, and report-ready figures in `analysis/figures/`.

### Global Bone Marrow Cancer Dataset

- **Dataset**: [Global Bone Marrow Cancer Dataset](https://www.kaggle.com/datasets/zkskhurram/global-bone-marrow-cancer-dataset)
- **Example directory**: `example/global-bone-marrow-cancer-dataset/`
- **Size**: Country-level and trend datasets covering myeloma/leukemia incidence, survival, bone marrow transplant access, hematologist availability, treatment patterns, and 2000-2026 trends.

This example demonstrates a public-health and health-services research workflow: whether cancer burden, treatment access, and survival outcomes differ systematically across regions. The README now uses `analysis/advanced_bone_marrow_analysis.ipynb` as the primary showcase. The notebook is written with the R kernel, saves figures to `analysis/advanced_figures/`, and saves model result tables to `analysis/advanced_results/`. The original topic-based R scripts are still kept in `analysis/`, making it easy to compare the basic scripted workflow with the advanced notebook workflow.

Main outputs in `analysis/`:

- `table1_by_continent.csv`: Descriptive statistics by continent. Europe and Oceania show higher myeloma 5-year survival, BMT access scores, and hematologists per million than regions such as Africa.
- `01_descriptive_stats_plots.R` to `06_trend_analysis.R`: Topic-based R scripts covering descriptive statistics, correlation tests, ANOVA/chi-square tests, regression, PCA/clustering, and trend analysis.
- `advanced_bone_marrow_analysis.ipynb`: Advanced R notebook workflow covering bootstrap regression stability, nonlinear GAM modeling, PCA-derived country phenotypes, k-means clustering, clustered heatmaps, therapy-era trend analysis, and income-region survival trajectories.
- `advanced_results/`: Advanced model tables, including standardized bootstrap coefficients, GAM smooth terms, PCA country phenotypes, cluster profiles, and interrupted trend models.
- `advanced_figures/`: Report-ready PNG figures for README, manuscripts, reports, or presentations.

Representative advanced figures are shown below.

**Bootstrap forest plot for standardized coefficient stability**

![Bootstrap coefficient stability](example/global-bone-marrow-cancer-dataset/analysis/advanced_figures/advanced_bootstrap_forest.png)

**GAM-based nonlinear BMT access-survival relationship**

![GAM BMT access and survival](example/global-bone-marrow-cancer-dataset/analysis/advanced_figures/advanced_gam_bmt_survival.png)

**PCA-derived country phenotypes with cluster interpretation**

![PCA country phenotypes](example/global-bone-marrow-cancer-dataset/analysis/advanced_figures/advanced_pca_country_phenotypes.png)

**Clustered heatmap of country-level burden, access, and outcome indicators**

![Clustered indicator heatmap](example/global-bone-marrow-cancer-dataset/analysis/advanced_figures/advanced_clustered_indicator_heatmap.png)

**Therapy-era trend panel**

![Therapy-era trend panel](example/global-bone-marrow-cancer-dataset/analysis/advanced_figures/advanced_therapy_era_trends.png)

Example plotting code:

```r
library(mgcv)
library(ggplot2)

gam_fit <- gam(
  Myeloma_5Y_Survival_Pct ~
    s(BMT_Access_Score, k = 5) +
    s(Hematologists_Per_Million, k = 5) +
    s(Myeloma_Incidence_Per_100K, k = 5),
  data = country,
  method = "REML"
)

gam_pred <- predict(gam_fit, newdata = gam_grid, se.fit = TRUE)
```

## Recommended workflows

- **Plain R scripts**: Use `r-script/` for users who do not use Jupyter. It produces `analysis.R` files that can be sourced in RStudio or run with `Rscript analysis.R`.
- **Jupyter Notebook**: Use for interactive exploration, teaching walkthroughs, step-by-step explanation, and `.ipynb` deliverables.
- **Quarto / R Markdown reports**: Use `quarto-report/` for formal reports, manuscript appendices, project summaries, and exportable HTML, Word, or PDF outputs.

## Latent profile analysis example

`advanced-stats/latent-profile-analysis/` provides a medical LPA workflow covering continuous-indicator selection, candidate Gaussian mixture models, class-number selection, posterior probabilities, assignment uncertainty, stability, and boundaries for external-variable inference.

### NHANES cardiometabolic profile example

In [`example/nhanes-lpa-subagent/`](example/nhanes-lpa-subagent/), an independent sub-agent used NHANES 2017–2018 adult BMI, waist circumference, and mean systolic blood pressure to identify empirical cardiometabolic profiles. Of 5,265 MEC-examined adults, 4,754 were complete cases; the selected candidate was `G=4, EEE`, with a smallest class of 5.5% and mean maximum posterior probability of 0.856. This is an unweighted, sample-internal exploration rather than a national estimate or validated clinical subtype. The full process is available in [`nhanes_lpa_report.Rmd`](example/nhanes-lpa-subagent/analysis/nhanes_lpa_report.Rmd) and [`nhanes_lpa_report.html`](example/nhanes-lpa-subagent/results/nhanes_lpa_report.html).

![Candidate model comparison](example/nhanes-lpa-subagent/figures/nhanes_lpa_model_selection.png)

![Standardized latent profiles](example/nhanes-lpa-subagent/figures/nhanes_lpa_profile.png)

![Profile sizes and classification results](example/nhanes-lpa-subagent/figures/nhanes_lpa_class_sizes.png)

## Install

Default Codex install:

```bash
curl -fsSL https://raw.githubusercontent.com/LeiGao0203/R-Medical-Statistics-Skills/main/install.sh | bash
```

For other coding agents, point the installer at that agent's skills, rules, or knowledge-base directory:

```bash
curl -fsSL https://raw.githubusercontent.com/LeiGao0203/R-Medical-Statistics-Skills/main/install.sh | AGENT_SKILLS_DIR=/path/to/agent/skills bash
```

You can also send this repository URL or the command above to an agent and let it complete the installation with its own tools. Restart or refresh the agent after installation. The skills will then be available for relevant medical statistics, R scripting, statistical report, and notebook tasks.

## License

This repository uses multiple licenses:

- Content in `basic-stats/`, `advanced-stats/`, and `literature-stats/` that is related to or adapted from *R语言实战医学统计* is adapted from [R_medical_stat](https://github.com/ayueme/R_medical_stat) by 阿越就是我 and is released under CC BY-SA 4.0.
- `r-script/`, `quarto-report/`, and `jupyter-notebook/` are original content and are released under Apache License 2.0.
- Data under `example/` comes from the corresponding public data sources; reuse should follow the terms listed on the original dataset pages.

See [LICENSE](LICENSE) for details.

## Attribution

Part of the medical statistics skill content is organized and adapted from *R语言实战医学统计*. When redistributing, modifying, or adapting related content, please preserve the original attribution and CC BY-SA 4.0 license notice.

## Contributing

Contributions are welcome, including new statistical methods, corrections to R examples, improvements to method-selection logic, and reproducible examples. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before contributing.
