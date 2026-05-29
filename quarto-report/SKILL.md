---
name: "quarto-report"
description: "Use when the user asks for Quarto, R Markdown, `.qmd`, `.Rmd`, statistical analysis reports, report-ready medical statistics results, or exports to HTML, Word, or PDF. Create reproducible R-based reports with narrative, code chunks, tables, figures, and interpretation."
---

# Quarto Report Skill

Create reproducible statistical analysis reports using Quarto (`.qmd`) or R Markdown (`.Rmd`). Prefer Quarto unless the user specifically asks for R Markdown or an existing `.Rmd` workflow.

## When to use

- The user asks for a report, statistical analysis report, research report, manuscript appendix, or project summary.
- The user asks for Quarto, `.qmd`, R Markdown, `.Rmd`, HTML, Word, or PDF output.
- The task needs narrative explanation plus executable R code.
- The user wants tables, figures, and interpretation in one deliverable.

## When not to use

- Use `r-script` when the user only needs a plain `.R` script.
- Use `jupyter-notebook` when the user asks for `.ipynb`, interactive notebooks, or notebook-style teaching material.
- Use the relevant medical statistics skill for method selection and statistical code before shaping the report.

## Output choice

Prefer `.qmd`:

```text
report.qmd
outputs/
├── figures/
└── tables/
```

Use `.Rmd` only when the user asks for R Markdown, has an existing R Markdown project, or needs compatibility with an older workflow.

## Recommended report skeleton

```markdown
---
title: "Statistical Analysis Report"
format:
  html: default
  docx: default
execute:
  echo: true
  warning: false
  message: false
---

# Research Objective

# Data and Variables

# Statistical Methods

# Descriptive Analysis

# Main Results

# Sensitivity or Subgroup Analysis

# Conclusion
```

For Chinese reports, use Chinese headings:

```markdown
# 研究目的

# 数据与变量

# 统计方法

# 描述性分析

# 主要结果

# 敏感性分析或亚组分析

# 结论
```

## Workflow

1. Lock the report purpose.
Identify audience, output format, statistical method, source data, outcome variables, grouping/exposure variables, covariates, and required tables or figures.

2. Choose `.qmd` or `.Rmd`.
Default to `report.qmd`. Use `report.Rmd` only when requested.

3. Build the report top-to-bottom.
Write short narrative sections around focused R chunks. Keep each code chunk responsible for one step: packages, data import, cleaning, descriptive statistics, model fitting, diagnostics, tables, figures, and interpretation.

4. Keep analysis code reproducible.
Use relative paths, stable output filenames, explicit package loading, and deterministic settings such as `set.seed()` when simulation, bootstrap, matching, or resampling is used.

5. Present medical statistics clearly.
Include:
- Variables and coding rules
- Missing-data handling
- Statistical assumptions or diagnostics
- Effect estimates with confidence intervals
- p values when appropriate
- Clinically meaningful interpretation

6. Save outputs deliberately.
When generating files, write tables to `outputs/tables/` and figures to `outputs/figures/`. The rendered report should still be readable without inspecting those folders.

7. Validate when possible.
For Quarto:

```bash
quarto render report.qmd
```

For R Markdown:

```bash
Rscript -e "rmarkdown::render('report.Rmd')"
```

If Quarto, R, or required packages are unavailable, say so explicitly and provide the exact command the user can run locally.

## Chunk style

Use named chunks when helpful:

````markdown
```{r}
#| label: setup
#| include: false
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

```{r}
#| label: load-packages
library(tidyverse)
library(broom)
```
````

For reports intended for non-technical readers, set `echo: false` at the document or chunk level when the user wants a polished output without visible code.

## Report writing style

- Keep method descriptions specific enough to reproduce the analysis.
- Do not overstate statistical significance as clinical importance.
- Report uncertainty, not only p values.
- Prefer concise tables and figures over raw console dumps.
- End with a short conclusion that matches the evidence from the analysis.

## Relationship to other skills

The statistical skills decide what analysis to run. This skill decides how to package the analysis as a report. Use `r-script` for runnable code handoff and `jupyter-notebook` for interactive exploration.
