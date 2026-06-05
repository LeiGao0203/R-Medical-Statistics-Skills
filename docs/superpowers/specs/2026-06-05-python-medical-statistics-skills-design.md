# Python Medical Statistics Skills Design

Date: 2026-06-05

## Objective

Create a new, independently named, independently published open source project:

```text
Python-Medical-Statistics-Skills
```

The project will provide medical statistics skills and reproducible Python analysis workflows for AI coding agents. It should be understandable as a companion in spirit to the existing R-focused medical statistics skills project, but it must stand on its own as a Python-native repository with its own name, README, installer, examples, git history, and GitHub publishing path.

The local repository should be created as a sibling of the current R project:

```text
/Users/leigao/Documents/
├── R-Medical-Statistics-Skills/
└── Python-Medical-Statistics-Skills/
```

## Scope

The first release will be an MVP that is small enough to publish quickly and complete enough to demonstrate real value.

It will include:

- A generic `SKILL.md` based directory structure for coding agents.
- Python-oriented medical statistics method guidance.
- Reproducible `.py` analysis workflow guidance.
- Jupyter notebook scaffolding guidance and helper script.
- One complete example analysis.
- A one-line installer similar in spirit to the R project installer.
- Basic install tests.
- Bilingual README files.
- Open source license and attribution files.

It will not include a Python package in the first release. A `pyproject.toml`, reusable Python library, or CLI package can be added later if the project evolves beyond skills and examples.

## Recommended Repository Layout

```text
Python-Medical-Statistics-Skills/
├── README.md
├── README.en.md
├── LICENSE
├── NOTICE
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
├── install.sh
├── tests/
│   └── install_test.sh
├── basic-stats/
│   ├── ttest/
│   │   └── SKILL.md
│   ├── anova/
│   │   └── SKILL.md
│   ├── chisq/
│   │   └── SKILL.md
│   ├── correlation/
│   │   └── SKILL.md
│   ├── nonparametric/
│   │   └── SKILL.md
│   ├── roc/
│   │   └── SKILL.md
│   └── sample-size/
│       └── SKILL.md
├── advanced-stats/
│   ├── logistic-reg/
│   │   └── SKILL.md
│   ├── multi-reg/
│   │   └── SKILL.md
│   ├── survival/
│   │   └── SKILL.md
│   └── pca/
│       └── SKILL.md
├── literature-stats/
│   ├── ps-matching/
│   │   └── SKILL.md
│   └── subgroup-analysis/
│       └── SKILL.md
├── python-script/
│   └── SKILL.md
├── jupyter-notebook/
│   ├── SKILL.md
│   ├── scripts/
│   │   └── new_notebook.py
│   ├── assets/
│   │   ├── experiment-template.ipynb
│   │   └── tutorial-template.ipynb
│   └── references/
│       ├── experiment-patterns.md
│       ├── notebook-structure.md
│       ├── quality-checklist.md
│       └── tutorial-patterns.md
└── example/
    └── lung-cancer/
        ├── data/
        └── analysis/
```

## MVP Skill Set

The first release will include 13 method skills:

- `basic-stats/ttest`
- `basic-stats/anova`
- `basic-stats/chisq`
- `basic-stats/correlation`
- `basic-stats/nonparametric`
- `basic-stats/roc`
- `basic-stats/sample-size`
- `advanced-stats/logistic-reg`
- `advanced-stats/multi-reg`
- `advanced-stats/survival`
- `advanced-stats/pca`
- `literature-stats/ps-matching`
- `literature-stats/subgroup-analysis`

Each method skill should include:

- Trigger-oriented frontmatter.
- When to use the method.
- When not to use the method.
- Data format expectations.
- Assumptions and diagnostic checks.
- Python package dependencies.
- Method selection guidance.
- A standard workflow.
- Minimal runnable Python examples.
- Report-ready interpretation guidance.

## Python Analysis Stack

The repository should be Python-native. Preferred packages:

- `pandas` and `numpy` for data manipulation.
- `scipy` for classical statistical tests.
- `statsmodels` for regression models, ANOVA, confidence intervals, and inference.
- `scikit-learn` for ROC, prediction metrics, PCA, preprocessing, and model validation.
- `lifelines` for survival analysis.
- `matplotlib` and `seaborn` for figures.
- `pingouin` may be recommended for convenient tests where it improves readability, but core examples should avoid unnecessary dependencies.

Each skill should avoid pretending that Python APIs are identical to R APIs. Code examples should use idiomatic Python data frames, explicit imports, and clear object names.

## Workflow Skills

### `python-script`

The Python script skill should generate reproducible `.py` workflows for users who want a plain script rather than a notebook. Its default output should be:

```text
analysis.py
```

For larger analyses, it should recommend:

```text
analysis/
├── analysis.py
├── data/
├── outputs/
│   ├── figures/
│   └── tables/
└── README.md
```

The default script structure should be:

```text
1. Imports
2. Configuration and paths
3. Data import
4. Data checks
5. Descriptive statistics
6. Assumption checks
7. Main analysis
8. Tables and figures
9. Report-ready interpretation
```

### `jupyter-notebook`

The notebook skill should reuse the strong pattern from the current project: template notebooks, a standard-library helper script, and references for experiment/tutorial structure and quality checks.

The notebook helper should not require third-party packages just to scaffold a notebook.

## Installer

The installer should mirror the proven behavior of the R project installer:

- Default target: `${CODEX_HOME:-$HOME/.codex}/skills`
- Override target with `AGENT_SKILLS_DIR`
- Support local execution from a checked-out repository.
- Support explicit source with `PY_MED_STATS_SOURCE`
- Support archive install with `PY_MED_STATS_ARCHIVE_URL`
- Copy only skill directories, not examples or repository metadata, into the target skills directory.

The expected one-line install command for Codex should be:

```bash
curl -fsSL https://raw.githubusercontent.com/LeiGao0203/Python-Medical-Statistics-Skills/main/install.sh | bash
```

For other coding agents:

```bash
curl -fsSL https://raw.githubusercontent.com/LeiGao0203/Python-Medical-Statistics-Skills/main/install.sh | AGENT_SKILLS_DIR=/path/to/agent/skills bash
```

## Example Analysis

The MVP should include one complete example under `example/lung-cancer/`.

The example should show:

- Data import with `pandas`.
- Missingness and variable type checks.
- Descriptive statistics and a Table 1 style output.
- Basic hypothesis tests.
- Logistic regression.
- ROC analysis.
- At least two report-ready figures.
- Saved CSV outputs and figures.

The example can reuse the existing public lung cancer dataset structure if license terms are acceptable, but the Python repository should make its own analysis outputs and README description.

## Documentation

`README.md` should be the primary Chinese README. `README.en.md` should be the English version.

The README should include:

- Project identity and positioning.
- Contents.
- Recommended workflows.
- Example showcase.
- One-line install.
- Other-agent install.
- License and attribution.
- Contributing notes.

The English README should mirror the Chinese README closely enough that the project is understandable to international users.

## Licensing And Attribution

The repository should use Apache License 2.0 for original workflow skills, notebook tooling, installer, and examples created for this project.

If method skill text is adapted from the existing R medical statistics skills or from `R_medical_stat`, the repository must clearly preserve attribution and license compatibility. The safest first-release approach is:

- Treat method knowledge adapted from the R project as CC BY-SA 4.0 derived content.
- Keep a `NOTICE` file describing the relationship to the source material.
- State the mixed-license arrangement clearly in both README files.
- Avoid copying long passages verbatim where a Python-native rewrite is feasible.

## Publishing Plan

First release milestones:

1. Create `/Users/leigao/Documents/Python-Medical-Statistics-Skills`.
2. Initialize git.
3. Add the repository skeleton.
4. Add README, README.en, license, notice, contributing, and code of conduct.
5. Add 13 MVP method skills.
6. Add `python-script` and `jupyter-notebook` workflow skills.
7. Add notebook scaffold script and templates.
8. Add one lung cancer example.
9. Add install script and install tests.
10. Run local validation.
11. Commit locally.
12. Create the GitHub repository.
13. Push `main`.
14. Optionally tag `v0.1.0`.

## Acceptance Criteria

The MVP is ready for GitHub when:

- The repository exists at `/Users/leigao/Documents/Python-Medical-Statistics-Skills`.
- `README.md` and `README.en.md` explain the project clearly.
- `install.sh` installs all MVP skills into a target skills directory.
- `tests/install_test.sh` passes.
- The 13 method skills exist and use Python-oriented dependencies and examples.
- `python-script/SKILL.md` and `jupyter-notebook/SKILL.md` exist.
- The notebook scaffold script creates valid notebook JSON.
- The example analysis contains runnable Python code and saved outputs, or clearly documents any unavailable runtime dependency.
- `git status --short` is clean after commit.

## Open Decisions

The MVP design is approved with these decisions already fixed:

- Independent name.
- Independent local directory.
- Independent GitHub repository.
- MVP scope rather than full migration.
- Skills plus reproducible examples, not a Python package for the first release.

The implementation plan should still decide exact wording, example dataset handling, and how much of each method skill is rewritten from scratch versus adapted with attribution.
