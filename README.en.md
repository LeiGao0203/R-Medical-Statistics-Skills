# R Medical Statistics Skills

[中文版](README.md)

A Codex skills collection for medical statistics and R-based statistical analysis. The project includes basic statistics, advanced statistics, commonly used methods in medical literature, and a Jupyter Notebook skill for creating and maintaining reproducible notebooks.

## Contents

- `basic-stats/`: Basic medical statistics skills, including t-tests, ANOVA, chi-square tests, correlation analysis, ROC analysis, sample size estimation, and statistical plotting.
- `advanced-stats/`: Advanced statistics skills, including ANCOVA, multiple regression, logistic regression, survival analysis, PCA, structural equation modeling, and multilevel models.
- `literature-stats/`: Skills for methods commonly seen in medical literature, including propensity score methods, Fine-Gray models, restricted cubic splines, subgroup analysis, and trend tests.
- `jupyter-notebook/`: Original Jupyter Notebook skill for creating, organizing, and validating reproducible notebooks.
- `example/`: Reserved for future demos, notebooks, datasets, and usage examples.

## Install

Copy the skill directories into your local Codex skills directory:

```bash
mkdir -p ~/.codex/skills
find advanced-stats basic-stats literature-stats -mindepth 1 -maxdepth 1 -type d -exec cp -R {} ~/.codex/skills/ \;
cp -R jupyter-notebook ~/.codex/skills/
```

Restart Codex after installation. The skills will then be available for relevant medical statistics and notebook tasks.

## License

This repository uses multiple licenses:

- Content in `basic-stats/`, `advanced-stats/`, and `literature-stats/` that is related to or adapted from *R语言实战医学统计* is adapted from [R_medical_stat](https://github.com/ayueme/R_medical_stat) by 阿越就是我 and is released under CC BY-SA 4.0.
- `jupyter-notebook/` is original content and is released under Apache License 2.0.
- Future content under `example/` should declare its own license in the relevant file or directory.

See [LICENSE](LICENSE) for details.

## Attribution

Part of the medical statistics skill content is organized and adapted from *R语言实战医学统计*. When redistributing, modifying, or adapting related content, please preserve the original attribution and CC BY-SA 4.0 license notice.

## Contributing

Contributions are welcome, including new statistical methods, corrections to R examples, improvements to method-selection logic, and reproducible examples. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before contributing.
