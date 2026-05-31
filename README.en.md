# R Medical Statistics Skills

[中文版](README.md)

A Codex skills collection for medical statistics and R-based statistical analysis. The project includes basic statistics, advanced statistics, commonly used methods in medical literature, and helper skills for creating plain R scripts, statistical reports, and Jupyter notebooks.

## Contents

- `basic-stats/`: Basic medical statistics skills, including t-tests, ANOVA, chi-square tests, correlation analysis, ROC analysis, sample size estimation, and statistical plotting.
- `advanced-stats/`: Advanced statistics skills, including ANCOVA, multiple regression, logistic regression, survival analysis, PCA, structural equation modeling, and multilevel models.
- `literature-stats/`: Skills for methods commonly seen in medical literature, including propensity score methods, Fine-Gray models, restricted cubic splines, subgroup analysis, and trend tests.
- `r-script/`: Original R script skill for generating reproducible `.R` analysis scripts for RStudio, command-line R, and non-notebook users.
- `quarto-report/`: Original Quarto/R Markdown report skill for generating medical statistics reports exportable to HTML, Word, or PDF.
- `jupyter-notebook/`: Original Jupyter Notebook skill for creating, organizing, and validating reproducible notebooks.
- `example/`: Reproducible examples based on public Kaggle datasets, including demo data, notebooks, R scripts, statistical tables, and generated figures.

## Examples

The `example/` directory currently includes three medical and health statistics examples that demonstrate typical workflows supported by these skills on real datasets. Each example usually contains raw data in `data/`, analysis code and result tables in `analysis/`, and visualization outputs in `analysis/figures/`.

### 1. Teen Mental Health

- **Dataset**: [Teen Mental Health](https://www.kaggle.com/datasets/argonnxx/teen-mental-health)
- **Example directory**: `example/teen-mental-health/`
- **Size**: 1200 records covering social media use, sleep, stress, anxiety, addiction level, mental health risk score, and depression labels among teenagers.

This example links adolescent digital behavior with mental health risk. It is useful for demonstrating a complete workflow that starts with exploratory analysis, moves through group comparisons, and ends with predictive modeling. The analysis first builds a Table 1 stratified by `digital_wellbeing_flag`, compares Healthy, Moderate, and At Risk groups, and then fits a logistic regression model for depression labels.

Main outputs in `analysis/`:

- `table1_by_wellbeing.csv`: Table 1 stratified by digital wellbeing group. The At Risk group has mean daily social media use of about 7.09 hours, compared with 4.84 hours in the Moderate group and 2.55 hours in the Healthy group.
- `multivariate_logistic_depression.csv`: Multivariable logistic regression results. Daily social media use and stress level are positively associated with depression labels, while longer sleep duration shows a protective association.
- `roc_depression_analysis.csv`: ROC summary for the depression prediction model, with AUC about 0.991, sensitivity about 0.968, and specificity about 0.953.
- `analysis/figures/`: Correlation heatmap, risk-score histogram, PCA plots, forest plot, and ROC curve.

Representative figures:

- `example/teen-mental-health/analysis/figures/correlation_heatmap.pdf`
- `example/teen-mental-health/analysis/figures/forest_plot_depression.pdf`
- `example/teen-mental-health/analysis/figures/roc_depression.pdf`

Example plotting code:

```r
library(ggplot2)
library(pROC)

roc_obj <- roc(df$depression_label, df$pred_prob)

ggroc(roc_obj, linewidth = 1.1, color = "#2C7FB8") +
  geom_abline(linetype = "dashed", color = "grey60") +
  theme_bw() +
  labs(
    title = "ROC curve for depression prediction",
    subtitle = paste0("AUC = ", round(auc(roc_obj), 3)),
    x = "1 - Specificity",
    y = "Sensitivity"
  )
```

### 2. Lung Cancer

- **Dataset**: [Lung Cancer](https://www.kaggle.com/datasets/mysarahmadbhat/lung-cancer)
- **Example directory**: `example/lung-cancer/`
- **Size**: 309 questionnaire records covering age, sex, smoking, anxiety, chronic disease, fatigue, coughing, shortness of breath, swallowing difficulty, chest pain, and lung cancer labels.

This example is closer to a medical questionnaire or case-control analysis. It demonstrates categorical-variable cleaning, baseline characteristic tables, group comparisons, univariable screening, multivariable logistic regression, and ROC evaluation. The goal is to identify symptoms and exposures associated with the lung cancer label and to show the resulting model's discrimination.

Main outputs in `analysis/`:

- `table1_baseline_characteristics.csv`: Baseline table comparing lung cancer and non-lung cancer groups. The lung cancer group shows higher proportions of yellow fingers, anxiety, peer pressure, fatigue, allergy, wheezing, and related symptoms.
- `univariate_logistic_regression.csv` and `multivariate_logistic_regression.csv`: Univariable and multivariable logistic regression results. In the multivariable model, coughing, swallowing difficulty, fatigue, chronic disease, and peer pressure show high odds ratios.
- `roc_analysis.csv`: ROC summary for the multivariable model, with AUC about 0.965, sensitivity about 0.863, and specificity about 0.974.
- `pca_loadings.csv`: PCA loadings for exploring the structure of symptom variables.

Representative figures:

- `example/lung-cancer/analysis/figures/bar_symptoms_prevalence.pdf`
- `example/lung-cancer/analysis/figures/forest_plot_univariate.pdf`
- `example/lung-cancer/analysis/figures/roc_curve_multivariable.pdf`
- `example/lung-cancer/analysis/figures/pca_biplot.pdf`

Example plotting code:

```r
library(ggplot2)

ggplot(symptom_summary, aes(x = reorder(symptom, prevalence), y = prevalence, fill = lung_cancer)) +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  labs(
    title = "Symptom prevalence by lung cancer status",
    x = NULL,
    y = "Prevalence",
    fill = "Lung cancer"
  )
```

### 3. Global Bone Marrow Cancer Dataset

- **Dataset**: [Global Bone Marrow Cancer Dataset](https://www.kaggle.com/datasets/zkskhurram/global-bone-marrow-cancer-dataset)
- **Example directory**: `example/global-bone-marrow-cancer-dataset/`
- **Size**: Country-level and trend datasets covering myeloma/leukemia incidence, survival, bone marrow transplant access, hematologist availability, treatment patterns, and 2000-2026 trends.

This example demonstrates a public-health and health-services research workflow: whether cancer burden, treatment access, and survival outcomes differ systematically across regions. Unlike the two notebook-based examples above, this example is organized as a set of R scripts, making it a useful showcase for the `r-script/` skill and reproducible `.R` analysis pipelines.

Main outputs in `analysis/`:

- `table1_by_continent.csv`: Descriptive statistics by continent. Europe and Oceania show higher myeloma 5-year survival, BMT access scores, and hematologists per million than regions such as Africa.
- `01_descriptive_stats_plots.R` to `06_trend_analysis.R`: Topic-based R scripts covering descriptive statistics, correlation tests, ANOVA/chi-square tests, regression, PCA/clustering, and trend analysis.
- `04_regression.R`: Multivariable linear regression for myeloma 5-year survival and logistic regression for high-survival countries.
- `05_pca_cluster_survival.R`: PCA and clustering of country-level indicators to explore combined patterns in incidence, survival, and health-care resources.

Representative figures:

- `example/global-bone-marrow-cancer-dataset/analysis/figures/boxplot_survival_by_continent.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/scatter_bmt_vs_survival.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/regression_bmt_survival.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/pca_biplot.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/kmeans_cluster.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/trend_survival.pdf`

Example plotting code:

```r
library(ggplot2)

ggplot(country, aes(x = BMT_Access_Score, y = Myeloma_5Y_Survival_Pct, color = Continent)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw() +
  labs(
    title = "BMT access and 5-year myeloma survival",
    x = "BMT access score",
    y = "Myeloma 5-year survival (%)",
    color = "Continent"
  )
```

## Recommended workflows

- **Plain R scripts**: Use `r-script/` for users who do not use Jupyter. It produces `analysis.R` files that can be sourced in RStudio or run with `Rscript analysis.R`.
- **Jupyter Notebook**: Use for interactive exploration, teaching walkthroughs, step-by-step explanation, and `.ipynb` deliverables.
- **Quarto / R Markdown reports**: Use `quarto-report/` for formal reports, manuscript appendices, project summaries, and exportable HTML, Word, or PDF outputs.

## Install

Copy the skill directories into your local Codex skills directory:

```bash
mkdir -p ~/.codex/skills
find advanced-stats basic-stats literature-stats -mindepth 1 -maxdepth 1 -type d -exec cp -R {} ~/.codex/skills/ \;
cp -R r-script ~/.codex/skills/
cp -R quarto-report ~/.codex/skills/
cp -R jupyter-notebook ~/.codex/skills/
```

Restart Codex after installation. The skills will then be available for relevant medical statistics, R scripting, statistical report, and notebook tasks.

## License

This repository uses multiple licenses:

- Content in `basic-stats/`, `advanced-stats/`, and `literature-stats/` that is related to or adapted from *R语言实战医学统计* is adapted from [R_medical_stat](https://github.com/ayueme/R_medical_stat) by 阿越就是我 and is released under CC BY-SA 4.0.
- `r-script/`, `quarto-report/`, and `jupyter-notebook/` are original content and are released under Apache License 2.0.
- Data under `example/` comes from the corresponding Kaggle datasets; reuse should follow the license terms listed on the Kaggle dataset pages.

See [LICENSE](LICENSE) for details.

## Attribution

Part of the medical statistics skill content is organized and adapted from *R语言实战医学统计*. When redistributing, modifying, or adapting related content, please preserve the original attribution and CC BY-SA 4.0 license notice.

## Contributing

Contributions are welcome, including new statistical methods, corrections to R examples, improvements to method-selection logic, and reproducible examples. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before contributing.
