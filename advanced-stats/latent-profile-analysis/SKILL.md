---
name: latent-profile-analysis
description: "医学统计中的潜在剖面分析（LPA）工作流：用连续临床、心理或生物标志物指标识别潜在人群剖面，比较高斯混合模型的类别数与方差—协方差结构，评估分类不确定性和稳定性，并分析协变量或远端结局。Use when a user asks for latent profile analysis, LPA, Gaussian mixture profiles, patient subtypes based on continuous indicators, tidyLPA, mclust, BCH, R3STEP, or a reproducible medical LPA report. Do not use ordinary clustering when probabilistic profile modeling and classification uncertainty are required; use latent class analysis for mainly categorical indicators."
---

# 潜在剖面分析（LPA）

## 核心定位

把 LPA 作为“测量模型—分类诊断—外部变量推断”三个环节处理。LPA 用同一基线时间窗内的多个连续指标估计高斯有限混合模型，输出类别比例、类别特异的均值/方差、个体后验概率和最可能类别。类别是模型对样本异质性的经验性表示，不应未经外部验证写成天然疾病亚型或因果暴露。

详细模型映射、数据集获取和报告清单见：

- [references/model-selection.md](references/model-selection.md)：`mclust` 协方差模型、类别数选择、分类质量和敏感性分析。
- [references/datasets.md](references/datasets.md)：对话中提到的数据集的公开性、访问限制、推荐指标和数据治理注意事项。

## 何时使用与何时换方法

使用 LPA 的典型条件：

- 每行是一名受试者，指标在明确的同一时间窗内测量；
- 剖面指标主要是连续变量，例如症状量表、生命体征、实验室指标或功能量表总分；
- 研究问题是“是否存在不同指标组合的人群”，而不是单一结局预测。

改用其他方法的情形：

- 指标主要为二分类/有序分类变量：使用 LCA 或广义混合模型；
- 同一指标跨多个时间点的轨迹：考虑 LCGA、GMM 或 LTA；
- 已有明确标签、目标是预测标签：使用监督分类；
- 只需按距离分组且不需要后验概率：使用聚类分析。

## 标准医学分析流程

1. **锁定研究人群和时间零点。** 规定纳入标准、重复就诊处理、基线时间窗和结局起算点；避免把治疗后变量用于基线剖面。
2. **预先指定指标。** 指标应有共同临床主题，避免把人口学变量、治疗方案、结局和高度重复的量表条目混进剖面定义。
3. **审查数据质量。** 检查变量类型、量表方向、范围、偏态、异常值、相关性、重复患者和缺失模式。高斯混合模型对异常值和重尾分布敏感。
4. **处理缺失。** 先报告每个指标的缺失率和完整病例选择。`mclust` 不接受缺失值；正式研究不能用未经论证的均值插补替代缺失机制分析。多重插补后要处理类别标签置换和类别数不一致。
5. **决定尺度。** 量纲不同或关注相对模式时，用基线分析样本的均值和标准差标准化；原始量表均值仍需保留用于临床解释。标准化不能修复偏态或异常值。
6. **比较候选模型。** 预先设定类别数范围（例如 1–6）和至少两种合理的协方差结构。使用 `mclust` 的 `EEI`、`VVI`、`EEE`、`VVV` 可近似覆盖常见的 LPA Model 1、2、3、6；不要只拟合最灵活的模型。
7. **综合选择类别数。** 同时看 BIC、SABIC/AIC、ICL、BLRT（若可用）、最小类别绝对人数、平均后验概率、收敛和临床可解释性。`mclust` 的 `bic` 是“越大越好”的定义，不能直接按传统软件的“越小越好”排序。
8. **诊断分类质量。** 输出每一类人数、平均后验概率、后验概率矩阵、分类不确定个体比例和剖面图。不要为了提高熵值而删除边界个体。
9. **命名并解释剖面。** 先按估计均值和原始量表临床阈值命名；使用“低症状负担型”“疼痛—疲劳突出型”等描述性名称，避免未经验证的生物学机制命名。
10. **检查稳定性。** 对最终模型做不同初始值/重复拟合、bootstrap 或发现集—验证集复制，并比较类别比例、均值和形状是否稳定。
11. **分析外部变量。** 预测因素用 R3STEP 或偏差校正方法；连续远端结局通常考虑 BCH；分类结局可考虑 DCAT。最可能类别的普通回归只能作为明确标注的敏感性分析，不能忽略分类误差。
12. **报告限制。** 说明抽样权重、复杂抽样设计、多中心聚集、测量不变性、缺失、类别稳定性和观察性研究的非因果性质。

## R 实现选择

优先用 `mclust` 做透明、可脚本化的候选模型网格；可用 `tidyLPA` 简化剖面图、`get_fit()` 和 `get_data()` 工作流。示例模块使用 `haven::read_xpt()` 读取 NHANES XPT 文件，并用 `mclust::Mclust()` 显式循环模型和类别数，以便保留每个候选解、分类结果和原始量表解释。

```r
library(mclust)

set.seed(20260716)
fit <- Mclust(
  data = scale(indicator_data),
  G = 1:6,
  modelNames = c("EEI", "VVI", "EEE", "VVV"),
  verbose = FALSE
)

# mclust::Mclust 返回：
# fit$classification：最可能类别
# fit$z：每个人对各类别的后验概率
# fit$bic：选中模型的 BIC（越大越好）
# fit$parameters：类别均值、方差/协方差和混合比例
```

## 分析过程的交付方式

统计方法与交付格式分开选择：

- 用户只要可运行代码或 RStudio/命令行执行：调用 `r-script`，输出 `.R`。
- 用户要带方法叙述、表格、图形和解释的医学分析报告：调用 `quarto-report`，优先 `.qmd`；已有 R Markdown 工作流时使用 `.Rmd`。
- 用户明确要交互式、分步探索：调用 `jupyter-notebook`，使用 R kernel 生成 `.ipynb`。

无论选择哪一种格式，都应复用同一套数据筛选、随机种子、候选模型、结果表和图形命名。不要因为换成 R Markdown 或 Notebook 而重新选择类别数或隐藏分类不确定性。本仓库的子 agent 示例选择 R Markdown 作为最终过程报告，同时保留可独立运行的 `.R` 源脚本。

## 仓库内可运行案例

完整 NHANES 2017–2018 案例位于 [example/nhanes-lpa](../../example/nhanes-lpa/)，保留如下结构：

```text
example/nhanes-lpa/
├── data/raw/       # CDC 官方 XPT 原始文件；不在 skill 安装时复制
├── data/derived/   # 合并、筛选和标准化后的分析数据
├── analysis/       # 下载和分析脚本
├── results/        # 模型选择、剖面均值、分类诊断、session 信息
├── figures/        # 剖面图
└── metadata/       # 数据源、变量和下载清单
```

从仓库根目录运行：

```bash
Rscript example/nhanes-lpa/analysis/01_nhanes_lpa.R
```

案例使用 BMI、腰围、平均收缩压、HbA1c、总胆固醇和 HDL-C 作为心代谢连续指标；由于 `mclust` 标准接口不直接处理 NHANES 复杂抽样权重，案例把剖面识别标记为未加权样本内探索，并保留 `WTMEC2YR`、`SDMVPSU`、`SDMVSTRA` 供后续调查设计敏感性研究使用。不要将该示例直接解释为美国总体的复杂抽样推断。

## 医学论文最小报告集

报告研究对象和时间窗、指标定义/方向、缺失处理、标准化方式、候选类别数、协方差结构、软件和版本、随机种子、所有合理候选模型的拟合指标、最终类别人数与原始量表均值、平均后验概率、分类不确定性、稳定性分析，以及外部变量方法如何传播分类误差。

> 许可证：本模块作为 `advanced-stats/` 医学统计内容发布，遵循仓库中的 CC BY-SA 4.0 许可证说明；NHANES 数据、R 软件包和第三方资料遵循各自的使用条款。
