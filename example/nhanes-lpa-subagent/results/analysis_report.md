# NHANES 2017–2018 minimal exploratory LPA

Run time: 2026-07-17 00:14:14
Analysis N: 4754; random seed: 20260716.

## Skill and analysis specification

This run followed `advanced-stats/latent-profile-analysis/SKILL.md`, `references/model-selection.md`, and `references/datasets.md`. The three pre-specified continuous indicators were BMI (BMXBMI), waist circumference (BMXWAIST), and the mean of available BPXSY1–BPXSY4 systolic readings (SBP_MEAN).

The local raw directory contained: BMX_J.XPT, BPX_J.XPT, DEMO_J.XPT, GHB_J.XPT, HDL_J.XPT, TCHOL_J.XPT. Only DEMO_J.XPT, BMX_J.XPT, and BPX_J.XPT were needed for these three indicators; the other local files were not used.

## Sample selection and missingness

Files were merged by SEQN. Participants were restricted to MEC-examined adults (RIDSTATR == 2 and RIDAGEYR >= 20). Ordinary `mclust` requires observed indicators, so complete cases on all three indicators were analyzed; no imputation was performed.

| stage | n | rule |
| --- | --- | --- |
| DEMO rows | 9254 | DEMO_J.XPT |
| Inner merge of DEMO/BMX/BPX | 8704 | Merge by SEQN |
| MEC examined | 8704 | RIDSTATR == 2 |
| Adults age >=20 | 5265 | RIDAGEYR >= 20 |
| Complete cases | 4754 | BMI, waist, and mean SBP all observed |

| indicator | label | unit | adult_n | missing_n | missing_pct |
| --- | --- | --- | --- | --- | --- |
| BMXBMI | Body mass index | kg/m^2 | 5265.00 | 90.00 | 1.71 |
| BMXWAIST | Waist circumference | cm | 5265.00 | 328.00 | 6.23 |
| SBP_MEAN | Mean systolic blood pressure | mmHg | 5265.00 | 266.00 | 5.05 |

## Standardization and model selection

The three indicators were standardized using the complete-case mean and SD. The derived file keeps the original-scale variables, z-score variables, MEC weight (WTMEC2YR), PSU (SDMVPSU), and strata (SDMVSTRA). Candidate models were G=1–4 with `EEI`, `VVI`, and `EEE` covariance structures. `mclust` BIC/AIC/SABIC/ICL were calculated with the larger BIC/ICL values indicating better relative fit.

For a non-singleton candidate to be considered eligible, the minimum class proportion had to be at least 5%, mean maximum posterior probability at least 0.80, and the proportion with maximum posterior below 0.70 no more than 20%. The selected model was the eligible candidate with the largest BIC; ICL and classification diagnostics are retained for joint interpretation.

Selected model: G=4 EEE.

| G | model | status | BIC | SABIC | ICL | min_class_n | min_class_prop | avg_max_posterior | uncertain_lt_0_70 | relative_entropy | quality_ok | selected |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4.000 | EEE | ok | -30696.923 | -30630.193 | -34119.352 | 262.000 | 0.055 | 0.856 | 0.180 | 0.740 | TRUE | TRUE |
| 3.000 | EEE | ok | -30752.848 | -30698.828 | -32082.673 | 265.000 | 0.056 | 0.945 | 0.065 | 0.873 | TRUE | FALSE |
| 2.000 | EEE | ok | -31649.269 | -31607.960 | -32674.973 | 482.000 | 0.101 | 0.956 | 0.050 | 0.844 | TRUE | FALSE |
| 1.000 | EEE | ok | -32264.236 | -32235.637 | -32264.236 | 4754.000 | 1.000 | 1.000 | 0.000 | 1.000 | FALSE | FALSE |
| 4.000 | VVI | ok | -32626.068 | -32540.271 | -34846.442 | 689.000 | 0.145 | 0.904 | 0.114 | 0.832 | TRUE | FALSE |
| 3.000 | VVI | ok | -34186.726 | -34123.173 | -36209.159 | 1111.000 | 0.234 | 0.912 | 0.104 | 0.806 | TRUE | FALSE |
| 4.000 | EEI | ok | -34603.227 | -34546.030 | -37010.825 | 286.000 | 0.060 | 0.897 | 0.134 | 0.817 | TRUE | FALSE |
| 3.000 | EEI | ok | -34966.308 | -34921.821 | -36841.122 | 501.000 | 0.105 | 0.916 | 0.099 | 0.821 | TRUE | FALSE |
| 2.000 | VVI | ok | -36619.894 | -36578.585 | -38232.930 | 1865.000 | 0.392 | 0.929 | 0.085 | 0.755 | TRUE | FALSE |
| 2.000 | EEI | ok | -37104.818 | -37073.042 | -38377.699 | 1261.000 | 0.265 | 0.943 | 0.070 | 0.807 | TRUE | FALSE |
| 1.000 | EEI | ok | -40521.603 | -40502.537 | -40521.603 | 4754.000 | 1.000 | 1.000 | 0.000 | 1.000 | FALSE | FALSE |
| 1.000 | VVI | ok | -40521.603 | -40502.537 | -40521.603 | 4754.000 | 1.000 | 1.000 | 0.000 | 1.000 | FALSE | FALSE |

## Classification quality and key results

| profile | label | n | proportion | average_assigned_posterior | n_below_0_70 |
| --- | --- | --- | --- | --- | --- |
| 1.000 | Profile 1 | 430.000 | 0.090 | 0.714 | 210.000 |
| 2.000 | Profile 2 | 400.000 | 0.084 | 0.836 | 102.000 |
| 3.000 | Lower cardiometabolic burden | 3662.000 | 0.770 | 0.875 | 484.000 |
| 4.000 | Higher cardiometabolic burden | 262.000 | 0.055 | 0.856 | 60.000 |

| assigned_class | n_assigned | posterior_class_1 | posterior_class_2 | posterior_class_3 | posterior_class_4 |
| --- | --- | --- | --- | --- | --- |
| 1.000 | 430.000 | 0.714 | 0.003 | 0.256 | 0.026 |
| 2.000 | 400.000 | 0.008 | 0.836 | 0.153 | 0.002 |
| 3.000 | 3662.000 | 0.089 | 0.029 | 0.875 | 0.007 |
| 4.000 | 262.000 | 0.061 | 0.007 | 0.076 | 0.856 |

Overall mean maximum posterior probability was 0.856; 856 of 4754 participants (18%) had maximum posterior probability below 0.70 and were retained.

Original-scale hard-class profile means:

| profile | label | indicator_label | unit | mean | n_hard_class |
| --- | --- | --- | --- | --- | --- |
| 1.000 | Profile 1 | Body mass index | kg/m^2 | 37.988 | 430.000 |
| 1.000 | Profile 1 | Waist circumference | cm | 128.018 | 430.000 |
| 1.000 | Profile 1 | Mean systolic blood pressure | mmHg | 122.641 | 430.000 |
| 2.000 | Profile 2 | Body mass index | kg/m^2 | 27.451 | 400.000 |
| 2.000 | Profile 2 | Waist circumference | cm | 97.255 | 400.000 |
| 2.000 | Profile 2 | Mean systolic blood pressure | mmHg | 169.251 | 400.000 |
| 3.000 | Lower cardiometabolic burden | Body mass index | kg/m^2 | 27.774 | 3662.000 |
| 3.000 | Lower cardiometabolic burden | Waist circumference | cm | 96.055 | 3662.000 |
| 3.000 | Lower cardiometabolic burden | Mean systolic blood pressure | mmHg | 122.460 | 3662.000 |
| 4.000 | Higher cardiometabolic burden | Body mass index | kg/m^2 | 46.690 | 262.000 |
| 4.000 | Higher cardiometabolic burden | Waist circumference | cm | 126.956 | 262.000 |
| 4.000 | Higher cardiometabolic burden | Mean systolic blood pressure | mmHg | 128.431 | 262.000 |

The accompanying figure plots standardized profile means. For orientation, Profile 1 has: Body mass index=37.99; Waist circumference=128.02; Mean systolic blood pressure=122.64.

## Limitations

- Complete-case analysis can change the target population and may be biased if missingness is informative.
- Ordinary `mclust::Mclust()` does not use NHANES WTMEC2YR, SDMVPSU, or SDMVSTRA. This is an unweighted, sample-internal exploratory LPA, not a nationally representative estimate; a formal study needs a survey-aware or sensitivity-analysis strategy.
- Gaussian mixtures are sensitive to skewness, outliers, covariance restrictions, initialization, and local solutions. The candidate comparison is not a substitute for replication or substantive validation.
- The three indicators define the profiles and should not be reused as independent external validation variables. Downstream analyses should account for classification uncertainty.
- The data are cross-sectional; profiles are descriptive and do not establish causality or clinical subtypes.

## Output files

The script is `analysis/agent_lpa_analysis.R`. Outputs are `data/derived/analysis_data.csv`, `results/model_selection.csv`, `results/class_sizes.csv`, `results/profile_means_original_scale.csv`, `results/posterior_classification.csv`, `results/analysis_report.md`, `results/session_info.txt`, and `figures/nhanes_lpa_profile.png`.
