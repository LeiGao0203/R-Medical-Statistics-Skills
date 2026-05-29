---
name: "r-script"
description: "Use when the user asks for R scripts (`.R`), RStudio-ready analysis, command-line R/Rscript workflows, or explicitly says they do not want Jupyter/notebooks. Create reproducible plain R analysis scripts for medical statistics tasks."
---

# R Script Skill

Create reproducible plain R analysis scripts for users who prefer RStudio, command-line R, or simple `.R` files instead of Jupyter notebooks.

## When to use

- The user asks for an `.R` script, RStudio workflow, command-line R, or `Rscript`.
- The user says they do not use Jupyter, notebooks, or `.ipynb`.
- The task is a medical statistics analysis that should be runnable as a script.
- The user needs code they can paste into RStudio or run from a terminal.

## Default output shape

For small tasks, create one script:

```text
analysis.R
```

For larger or reusable analyses, use a compact project layout:

```text
analysis/
├── analysis.R
├── data/
├── outputs/
│   ├── figures/
│   └── tables/
└── README.md
```

Only create the larger layout when the user needs multiple outputs, repeated execution, or handoff to another person.

## Workflow

1. Clarify the analysis target.
Identify the statistical method, outcome variable, grouping/exposure variables, covariates, data source, and intended report format.

2. Choose script scope.
Use a single `analysis.R` for focused tasks. Use the compact project layout for analyses with imported data, saved tables, saved figures, or multiple steps.

3. Write a top-to-bottom R script.
Organize the script with short section comments:

```r
# 1. Packages
# 2. Data import
# 3. Data checks
# 4. Descriptive statistics
# 5. Assumption checks
# 6. Main analysis
# 7. Tables and figures
# 8. Report-ready interpretation
```

4. Keep the statistical method environment neutral.
Use the relevant medical statistics skill for method choice and R code. Do not introduce notebook-only assumptions such as cells, markdown displays, or inline notebook output.

5. Make dependencies explicit.
Prefer base R when it is enough. When packages are needed, list them at the top:

```r
packages <- c("tidyverse", "broom")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))
```

Use only packages that the analysis actually needs.

6. Save outputs deliberately.
For scripts that generate artifacts, write tables to `outputs/tables/` and figures to `outputs/figures/`. Use stable filenames and avoid overwriting the input data.

7. Validate when possible.
Run the script with:

```bash
Rscript analysis.R
```

If R is unavailable, say so explicitly and provide the exact command the user can run locally.

## Script style

- Use clear object names such as `df`, `fit`, `summary_table`, and `plot_roc`.
- Keep transformations readable; avoid dense one-liners for clinical or teaching analyses.
- Print concise intermediate checks: `str(df)`, `summary(df)`, missingness counts, and group sizes.
- Include assumption checks before inferential tests when the method requires them.
- End with a short report-ready interpretation in comments or `cat()` output.

## RStudio notes

When the user mentions RStudio:

- Keep paths relative to the project folder.
- Avoid requiring terminal-only setup unless necessary.
- Mention that the user can open `analysis.R` and run it top-to-bottom with Source.

## Relationship to notebooks

Use this skill instead of `jupyter-notebook` when the user wants plain R. Use `jupyter-notebook` only when the user asks for `.ipynb`, notebooks, interactive tutorials, or notebook-style exploratory work.
