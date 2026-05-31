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

`example/` 目录提供一个较完整、代表性更强的医学统计案例，用于展示本项目 skills 如何把公开数据整理成可复现的 R 分析流程。示例包含 `data/` 原始数据、`analysis/` 分析脚本与结果表，以及 `analysis/figures/` 中可直接用于报告的图表。

### Global Bone Marrow Cancer Dataset

- **数据来源**：[Global Bone Marrow Cancer Dataset](https://www.kaggle.com/datasets/zkskhurram/global-bone-marrow-cancer-dataset)
- **示例目录**：`example/global-bone-marrow-cancer-dataset/`
- **数据规模**：包含国家层面的骨髓瘤/白血病发病率、生存率、骨髓移植可及性、血液科医生数量、治疗方式和 2000-2026 年趋势数据。

这个示例用于展示更偏公共卫生和卫生服务研究的问题：不同地区的血液肿瘤负担、治疗资源可及性与生存结局之间是否存在系统差异。与前两个 notebook 示例不同，该示例主要使用拆分的 R 脚本组织分析流程，适合展示 `r-script/` skill 生成可复现 `.R` 分析管线的能力。

主要分析结果保存在 `analysis/`：

- `table1_by_continent.csv`：按大陆汇总的描述性统计。欧洲和大洋洲的骨髓瘤 5 年生存率、BMT 可及性评分和每百万人血液科医生数量整体高于非洲等地区。
- `01_descriptive_stats_plots.R` 至 `06_trend_analysis.R`：按分析主题拆分的 R 脚本，覆盖描述性统计、相关检验、ANOVA/卡方检验、回归、PCA/聚类和趋势分析。
- `04_regression.R`：以骨髓瘤 5 年生存率为结局构建多元线性回归，并用 Logistic 回归分析高生存率国家的相关因素。
- `05_pca_cluster_survival.R`：对国家层面指标进行 PCA 和聚类，探索国家在发病率、生存率和医疗资源上的综合分组。

代表性图表如下，README 中使用 PNG 版本直接展示；对应 PDF 版本仍保留在同一目录，便于报告排版或论文附录使用。

**不同大陆的骨髓瘤 5 年生存率**

![Myeloma 5-year survival by continent](example/global-bone-marrow-cancer-dataset/analysis/figures/boxplot_survival_by_continent.png)

**BMT 可及性与骨髓瘤 5 年生存率**

![BMT access and myeloma 5-year survival](example/global-bone-marrow-cancer-dataset/analysis/figures/scatter_bmt_vs_survival.png)

**国家层面指标 PCA 双标图**

![PCA biplot for country-level indicators](example/global-bone-marrow-cancer-dataset/analysis/figures/pca_biplot.png)

**2000-2026 年全球生存率趋势**

![Global survival trend](example/global-bone-marrow-cancer-dataset/analysis/figures/trend_survival.png)

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
