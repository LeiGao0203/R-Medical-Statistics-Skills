# R Medical Statistics Skills

[English version](README.en.md)

一套面向医学统计和 R 语言分析场景的 Codex skills。项目包含基础统计、高级统计、医学文献常用统计方法，以及用于生成普通 R 脚本、统计报告和维护 Jupyter Notebook 的辅助 skill。

## Contents

- `basic-stats/`: t 检验、方差分析、卡方检验、相关分析、ROC、样本量、统计绘图等基础医学统计 skill。
- `advanced-stats/`: 协方差分析、多元回归、Logistic 回归、生存分析、PCA、结构方程、多水平模型等进阶统计 skill。
- `literature-stats/`: 倾向评分、Fine-Gray、限制性立方样条、亚组分析、趋势检验等医学文献常见方法 skill。
- `r-script/`: 原创 R 脚本 skill，用于为 RStudio、命令行 R 或非 notebook 用户生成可复现 `.R` 分析脚本。
- `quarto-report/`: 原创 Quarto/R Markdown 报告 skill，用于生成可导出 HTML、Word 或 PDF 的医学统计分析报告。
- `jupyter-notebook/`: 原创 Jupyter Notebook skill，用于创建、整理和验证可复现 notebook。
- `example/`: 基于公开 Kaggle 数据集整理的可复现示例，包含演示数据、notebook、R 脚本、统计表和图形结果。

## Examples

`example/` 目录目前包含 3 个医学与健康统计分析示例，可用于展示本项目 skills 在真实数据上的典型工作流。每个示例通常包含 `data/` 原始数据、`analysis/` 分析代码与结果表，以及 `analysis/figures/` 中的可视化结果。

### 1. Teen Mental Health

- **数据来源**：[Teen Mental Health](https://www.kaggle.com/datasets/argonnxx/teen-mental-health)
- **示例目录**：`example/teen-mental-health/`
- **数据规模**：1200 条青少年社交媒体使用、睡眠、压力、焦虑、成瘾程度、心理健康风险评分和抑郁标签记录。

这个示例把青少年数字行为与心理健康风险联系起来，适合演示“探索性分析 + 组间比较 + 预测建模”的完整流程。分析从 `digital_wellbeing_flag` 分层的基线表开始，比较 Healthy、Moderate 和 At Risk 三组在社交媒体时长、睡眠、压力、焦虑和学业表现上的差异，再进一步构建抑郁标签的 Logistic 回归模型。

主要分析结果保存在 `analysis/`：

- `table1_by_wellbeing.csv`：按数字健康分组的 Table 1。At Risk 组平均每日社交媒体使用约 7.09 小时，高于 Moderate 组的 4.84 小时和 Healthy 组的 2.55 小时。
- `multivariate_logistic_depression.csv`：多因素 Logistic 回归结果。每日社交媒体使用、压力水平与抑郁标签呈正相关，睡眠时长呈保护性关联。
- `roc_depression_analysis.csv`：抑郁预测模型 ROC 结果，AUC 约 0.991，灵敏度约 0.968，特异度约 0.953。
- `analysis/figures/`：包含相关热图、风险评分直方图、PCA 图、森林图和 ROC 曲线。

代表性图表：

- `example/teen-mental-health/analysis/figures/correlation_heatmap.pdf`
- `example/teen-mental-health/analysis/figures/forest_plot_depression.pdf`
- `example/teen-mental-health/analysis/figures/roc_depression.pdf`

示例图表代码片段：

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

- **数据来源**：[Lung Cancer](https://www.kaggle.com/datasets/mysarahmadbhat/lung-cancer)
- **示例目录**：`example/lung-cancer/`
- **数据规模**：309 条肺癌风险因素问卷记录，包含年龄、性别、吸烟、焦虑、慢性病、疲劳、咳嗽、气促、吞咽困难、胸痛和肺癌标签。

这个示例更接近医学问卷或病例-对照数据分析场景，适合演示分类变量整理、基线特征表、组间检验、单因素筛选、多因素 Logistic 回归和 ROC 评价。分析目标是识别与肺癌标签相关的症状和暴露因素，并展示模型的区分能力。

主要分析结果保存在 `analysis/`：

- `table1_baseline_characteristics.csv`：肺癌组与非肺癌组基线表。肺癌组在黄手指、焦虑、同伴压力、疲劳、过敏、喘息等变量上的比例更高。
- `univariate_logistic_regression.csv` 和 `multivariate_logistic_regression.csv`：单因素和多因素 Logistic 回归结果。多因素模型中，咳嗽、吞咽困难、疲劳、慢性病、同伴压力等变量显示较高 OR。
- `roc_analysis.csv`：多因素模型 ROC 结果，AUC 约 0.965，灵敏度约 0.863，特异度约 0.974。
- `pca_loadings.csv`：PCA 载荷结果，用于观察症状变量之间的综合结构。

代表性图表：

- `example/lung-cancer/analysis/figures/bar_symptoms_prevalence.pdf`
- `example/lung-cancer/analysis/figures/forest_plot_univariate.pdf`
- `example/lung-cancer/analysis/figures/roc_curve_multivariable.pdf`
- `example/lung-cancer/analysis/figures/pca_biplot.pdf`

示例图表代码片段：

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

- **数据来源**：[Global Bone Marrow Cancer Dataset](https://www.kaggle.com/datasets/zkskhurram/global-bone-marrow-cancer-dataset)
- **示例目录**：`example/global-bone-marrow-cancer-dataset/`
- **数据规模**：包含国家层面的骨髓瘤/白血病发病率、生存率、骨髓移植可及性、血液科医生数量、治疗方式和 2000-2026 年趋势数据。

这个示例用于展示更偏公共卫生和卫生服务研究的问题：不同地区的血液肿瘤负担、治疗资源可及性与生存结局之间是否存在系统差异。与前两个 notebook 示例不同，该示例主要使用拆分的 R 脚本组织分析流程，适合展示 `r-script/` skill 生成可复现 `.R` 分析管线的能力。

主要分析结果保存在 `analysis/`：

- `table1_by_continent.csv`：按大陆汇总的描述性统计。欧洲和大洋洲的骨髓瘤 5 年生存率、BMT 可及性评分和每百万人血液科医生数量整体高于非洲等地区。
- `01_descriptive_stats_plots.R` 至 `06_trend_analysis.R`：按分析主题拆分的 R 脚本，覆盖描述性统计、相关检验、ANOVA/卡方检验、回归、PCA/聚类和趋势分析。
- `04_regression.R`：以骨髓瘤 5 年生存率为结局构建多元线性回归，并用 Logistic 回归分析高生存率国家的相关因素。
- `05_pca_cluster_survival.R`：对国家层面指标进行 PCA 和聚类，探索国家在发病率、生存率和医疗资源上的综合分组。

代表性图表：

- `example/global-bone-marrow-cancer-dataset/analysis/figures/boxplot_survival_by_continent.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/scatter_bmt_vs_survival.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/regression_bmt_survival.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/pca_biplot.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/kmeans_cluster.pdf`
- `example/global-bone-marrow-cancer-dataset/analysis/figures/trend_survival.pdf`

示例图表代码片段：

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

## 推荐使用方式

- **普通 R 脚本**：不使用 Jupyter 的用户优先使用 `r-script/`，输出 `analysis.R`，可在 RStudio 中 Source，或用 `Rscript analysis.R` 在命令行运行。
- **Jupyter Notebook**：适合交互式探索、教学演示、逐步解释和需要 `.ipynb` 交付的任务。
- **Quarto / R Markdown 报告**：使用 `quarto-report/`，适合正式报告、论文附录、课题汇报，以及可导出的 HTML、Word、PDF。

## Install

将各个 skill 目录复制到本机 Codex skills 目录：

```bash
mkdir -p ~/.codex/skills
find advanced-stats basic-stats literature-stats -mindepth 1 -maxdepth 1 -type d -exec cp -R {} ~/.codex/skills/ \;
cp -R r-script ~/.codex/skills/
cp -R quarto-report ~/.codex/skills/
cp -R jupyter-notebook ~/.codex/skills/
```

重启 Codex 后，相关 skill 会在对应医学统计、R 脚本、统计报告或 notebook 任务中触发。

## License

本项目采用混合许可证：

- `basic-stats/`、`advanced-stats/`、`literature-stats/` 中与《R语言实战医学统计》相关的内容，改编自阿越就是我的开源项目 [R_medical_stat](https://github.com/ayueme/R_medical_stat)，按 CC BY-SA 4.0 发布。
- `r-script/`、`quarto-report/` 和 `jupyter-notebook/` 为原创内容，按 Apache License 2.0 发布。
- `example/` 中的数据来自对应 Kaggle 数据集；再次使用时请遵守 Kaggle 页面及原数据集的授权条款。

详见 [LICENSE](LICENSE)。

## Attribution

统计类 skill 的部分内容基于《R语言实战医学统计》整理和改写。再次分发、修改或演绎相关内容时，请保留原作者署名和 CC BY-SA 4.0 授权信息。

## Contributing

欢迎补充新的统计方法、修正代码示例、改进方法选择逻辑或增加可复现实例。提交前请阅读 [CONTRIBUTING.md](CONTRIBUTING.md)。
