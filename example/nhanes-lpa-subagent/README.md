# NHANES LPA：子 agent 独立分析示例

这是由独立子 agent 读取 `advanced-stats/latent-profile-analysis/` 后完成的可复现 NHANES 2017–2018 LPA 示例。它没有复用原有 `example/nhanes-lpa/` 的分析脚本，而是重新选择了三个连续体检指标：BMI、腰围和可用收缩压读数的均值。

## 运行

从仓库根目录执行：

```bash
Rscript example/nhanes-lpa-subagent/analysis/agent_lpa_analysis.R
```

所需 R 包：`haven`、`mclust` 和 `ggplot2`。

## 目录结构

- `data/raw/`：CDC/NCHS NHANES 2017–2018 公开 XPT 文件。
- `data/derived/analysis_data.csv`：按 `SEQN` 合并、限制 MEC 成人并完成病例筛选后的分析数据。
- `analysis/agent_lpa_analysis.R`：子 agent 实际运行的独立 R 脚本。
- `analysis/nhanes_lpa_report.Rmd`：本次交付的 R Markdown 过程报告源文件。
- `results/analysis_report.md`：从读取 skill、样本筛选到模型解释的过程与结果记录。
- `results/nhanes_lpa_report.html`：已渲染的 HTML 报告。
- `results/model_selection.csv`：`G=1–4`、`EEI/VVI/EEE` 候选模型的 BIC、SABIC、ICL 和分类质量。
- `results/class_sizes.csv`、`profile_means_original_scale.csv`、`posterior_classification.csv`：类别规模、原始量表剖面均值和个体后验概率。
- `figures/nhanes_lpa_profile.png`：标准化剖面图。

## 子 agent 结果摘要

- MEC 成人：5,265 人；三个指标完整病例：4,754 人。
- 缺失率：BMI 1.71%、腰围 6.23%、平均收缩压 5.05%。
- 最终探索性模型：`G=4, EEE`；最小类别 262 人（5.5%）。
- 平均最大后验概率 0.856；856 人最大后验概率低于 0.70，但均被保留。

该分析使用普通 `mclust`，没有使用 NHANES 的权重、PSU 和分层变量进行复杂抽样推断；因此只能作为样本内探索性 LPA，不能直接解释为美国总体估计、临床亚型或因果结论。

数据来源：[NHANES 2017–2018 官方数据入口](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017)。

## 渲染报告

从仓库根目录执行：

```r
rmarkdown::render(
  "example/nhanes-lpa-subagent/analysis/nhanes_lpa_report.Rmd",
  output_dir = "example/nhanes-lpa-subagent/results",
  knit_root_dir = getwd()
)
```
