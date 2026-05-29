---
name: medical-stat-factor-analysis
description: "R语言医学统计：探索性因子分析（EFA）。使用psych和FactoMineR包进行因子提取（主成分法/主轴因子法/极大似然法）、因子旋转（varimax/oblimin）、因子解释及信度分析（Cronbach's α）。TRIGGER when user mentions 因子分析、EFA、因子载荷、Cronbach's α、信度分析、量表结构，or asks about exploring latent factor structure of a scale. SKIP for 主成分分析、结构方程模型（验证性因子分析）。"
---

# 探索性因子分析 (Exploratory Factor Analysis, EFA)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**什么时候用：**

- 存在无法直接观测的潜在变量（不可测现象），需要通过多个可观测指标间接反映（如用语言能力、记忆能力、思维能力反映意识清醒状态）
- 需要对量表/问卷进行结构探索，确定观测变量背后隐含的因子个数及各变量在因子上的载荷
- 需要简化数据结构，用少数几个公共因子解释多个变量之间的相关性
- 量表编制或文化适应阶段，不清楚条目归属于几个维度时进行初步探索
- 医学综合评价体系中需要提取潜在的"综合因子""水平因子""数量因子"等专业意义明确的维度

**什么时候不用：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 重点在综合原始变量信息而非解释变量间关系 | 主成分分析（PCA）（`prcomp()`） |
| 已有理论假设，需要验证特定因子结构与数据的拟合程度 | 验证性因子分析（CFA）（`lavaan::cfa()`） |
| 需要降维后做主成分回归预测因变量 | 主成分回归（PCR）（`pls::pcr()`） |
| 需要将样本分成若干类别 | 聚类分析（`stats::hclust()`） |
| 各变量之间相关性很低（KMO < 0.6） | 不适合因子分析，数据不存在提取公共因子的基础 |

**医学研究常见应用：**

- 医院医疗质量综合评价：门诊人次、出院人数、病床利用率、病床周转次数、平均住院天数、治愈好转率、病死率、诊断符合率、抢救成功率等指标提取综合因子
- 生命质量量表（SF-36）的结构探索：将36个条目归为生理健康、心理健康等维度
- 疼痛评估工具编制：将多个疼痛描述词归纳为感觉疼痛、情感疼痛等少数维度
- 中医证候规范研究：数十个症状指标降维为少数证候因子

## 前置条件

**R包：**

```r
# 核心因子分析包
install.packages("psych")
# 多方法确定因子个数
install.packages("parameters")
# 综合因子分析（备选方案，提供更多旋转方法和结果可视化）
install.packages("FactoMineR")
# KMO/Bartlett检验及模型诊断
install.packages("performance")
# 可视化辅助
install.packages(c("see", "pheatmap"))
# 读入SPSS数据
install.packages("foreign")
```

**数据格式：**

- 每一行代表一个观测样本，每一列代表一个观测变量
- 所有变量必须是数值型连续变量
- 因子分析使用相关矩阵，变量量纲不一致时需先标准化（可设 `cor = "cor"` 或在分析前调用 `scale()`）
- 建议样本量与变量数之比至少为5:1，一般要求样本量 ≥100 或变量数的5-10倍

**统计前提：**

- KMO取样充分性检验（Kaiser-Meyer-Olkin）：KMO > 0.6 勉强可接受，> 0.7 一般，> 0.8 良好，> 0.9 极佳
- Bartlett球形检验：p < 0.05 拒绝变量间相互独立的零假设，表示数据适合因子分析
- 变量之间应有足够的相关性（相关矩阵中多数相关系数绝对值 > 0.3）
- 不存在严重的多重共线性（行列式值不宜过小）

## 方法选择决策树

```
你的分析目标 →
├── 确定因子个数 →
│   ├── 初步探索 → 碎石图（fa.parallel） + 平行分析
│   ├── 多方法综合 → parameters::n_factors()（含18种方法）
│   ├── 特征值准则 → 保留特征值 > 1 的因子（Kaiser准则）
│   └── 累积方差贡献率 → 保留至累积贡献率 > 60%-70% 的因子数
├── 因子提取方法（fm参数） →
│   ├── 数据满足多元正态分布 → 最大似然法（fm="ml"），可进行假设检验和拟合指数
│   ├── 数据偏离正态或样本量较小 → 主轴因子法（fm="pa"），对分布要求宽松
│   ├── 一般情况下的默认选择 → 最小残差法（fm="minres"）
│   └── 特殊情况 → 加权最小二乘（fm="wls"）或广义加权最小二乘（fm="gls"）
├── 因子旋转方法（rotate参数） →
│   ├── 假设因子之间相互独立 → 正交旋转（如 rotate="varimax" 最大方差法）
│   ├── 理论或经验表明因子间可能存在相关 → 斜交旋转（rotate="oblimin"）
│   └── 初步探索、尚未决定 → 先不旋转（rotate="none"）观察初始载荷
└── 信度分析（补充） →
    ├── 已确定量表各条目归属哪个因子 → Cronbach's α（psych::alpha()）
    └── α < 0.7 → 考虑删除某些条目或增加条目数量
```

## 标准工作流

### 步骤1：数据准备与相关性检验

- 读入数据，提取数值型变量用于分析
- 计算相关矩阵，初步观察变量间相关性
- 执行 KMO 检验和 Bartlett 球形检验（用 `performance::check_factorstructure()` 或 `psych::KMO()` + `psych::cortest.bartlett()`）

### 步骤2：确定最佳因子个数

- 绘制碎石图并进行平行分析（`psych::fa.parallel()`）
- 使用 `parameters::n_factors()` 综合18种方法的结果
- 结合特征值 > 1 准则和累积方差贡献率（一般 > 60%）
- 结合专业知识判断因子个数的合理性

### 步骤3：执行因子分析

- 选择因子提取方法（ml / pa / minres / wls / gls）
- 先不旋转，观察初始因子载荷矩阵，评估各因子的意义
- 尝试不同的提取方法和因子个数，比较累积方差贡献率和载荷模式

### 步骤4：因子旋转

- 若初始载荷矩阵各因子的专业意义模糊，进行因子旋转
- 正交旋转（varimax）使各因子彼此独立，简化因子解释
- 斜交旋转（oblimin）允许因子间相关，在心理量表等场景中更符合实际
- 旋转后重新审视各因子对应的高载荷变量（一般 > 0.4 或 > 0.5），为因子命名

### 步骤5：结果报告（如何写论文中的统计描述）

论文中因子分析部分的报告要点：
- 报告 KMO 值和 Bartlett 球形检验结果：「KMO 值为 0.XX，Bartlett 球形检验 χ² = XX.XX，P < 0.001，表明数据适合进行因子分析」
- 报告因子提取方法、旋转方法、提取的因子个数及依据
- 报告累积方差贡献率：「前 X 个因子累积方差贡献率为 68.XX%」
- 列出旋转后的因子载荷矩阵（表格形式），标注每个因子的高载荷变量
- 结合专业背景为各因子命名并解释其意义
- 若涉及量表信度分析，报告各维度的 Cronbach's α 系数

## 代码示例

```r
# 加载包
library(psych)
library(parameters)
library(performance)
library(see)
library(pheatmap)

# === 数据准备 ===
df <- foreign::read.spss("datasets/例22-02.sav", to.data.frame = TRUE, reencode = "utf-8")
names(df) <- c("年", "月", "门诊人次", "出院人数", "病床利用率",
               "病床周转次数", "平均住院天数", "治愈好转率",
               "病死率", "诊断符合率", "抢救成功率")
df.use <- df[, -c(1, 2)]  # 去掉年、月两列

# === 前提条件检验：KMO 和 Bartlett ===
performance::check_factorstructure(df.use)
# 或分别调用：
# psych::KMO(df.use)            # KMO取样充分性
# psych::cortest.bartlett(df.use) # Bartlett球形检验

# === 确定最佳因子个数 ===

# 方法1：碎石图 + 平行分析
set.seed(1)
fa.parallel(df.use, fa = "fa", fm = "ml")
## Parallel analysis suggests that the number of factors = 3

# 方法2：多方法综合判断
n <- parameters::n_factors(df.use)
n
## The choice of 3 dimensions is supported by 6 (31.58%) methods out of 19

# 查看各方法详情
as.data.frame(n)
##    n_Factors              Method              Family
## 2          2             Bentler             Bentler
## 5          2       Velicer's MAP        Velicers_MAP
## 10         3   Parallel analysis               Scree
## 11         3    Kaiser criterion               Scree
## ...（共18种方法）

summary(n)
##   n_Factors n_Methods Variance_Cumulative
## 1         2         5           0.4936863
## 2         3         6           0.6385558
## 3         4         1           0.7072393

# 可视化因子数选择
plot(n) + theme_modern()

# === 因子分析（4因子，最大似然法，最大方差旋转） ===
fa.res <- fa(df.use, nfactors = 4, rotate = "varimax", fm = "ml")
fa.res
## Factor Analysis using method =  ml
## Standardized loadings (pattern matrix):
##                ML3   ML1   ML2   ML4   h2    u2 com
## 门诊人次     -0.31  0.23 -0.03  0.92 1.00 0.005 1.4
## 出院人数      0.75  0.16  0.24  0.27 0.72 0.276 1.6
## 病床利用率   -0.10  0.83  0.03  0.07 0.71 0.289 1.0
## 病床周转次数  0.46  0.84  0.09  0.26 1.00 0.005 1.8
## 平均住院天数 -0.64 -0.23  0.24  0.21 0.57 0.435 1.8
## 治愈好转率   -0.09 -0.09  0.98 -0.10 1.00 0.005 1.1
## 病死率       -0.20 -0.18 -0.42 -0.06 0.26 0.743 1.9
## 诊断符合率   -0.56  0.02 -0.10  0.18 0.36 0.642 1.3
## 抢救成功率    0.70 -0.04  0.04 -0.21 0.55 0.455 1.2
##
##                        ML3  ML1  ML2  ML4
## SS loadings           2.15 1.58 1.29 1.12
## Proportion Var        0.24 0.18 0.14 0.12
## Cumulative Var        0.24 0.41 0.56 0.68
## Proportion Explained  0.35 0.26 0.21 0.18
## Cumulative Proportion 0.35 0.61 0.82 1.00
##
## Tucker Lewis Index of factoring reliability =  0.931
## RMSEA index =  0.055  and the 90 % confidence intervals are  0 0.235
## BIC =  -14.67

# === 提取各组件 ===
fa.res$loadings        # 因子载荷矩阵
fa.res$communality     # 公因子方差 (h2)
fa.res$uniquenesses    # 独特性 (u2)
fa.res$complexity      # 复杂性 (com)
fa.res$Vaccounted      # 载荷平方和、方差解释比例等
fa.res$STATISTIC       # 似然比卡方
fa.res$PVAL            # P值

# === 尝试不同的提取/旋转方法 ===

# 最小残差法，不旋转
fa.res_minres <- fa(df.use, nfactors = 4, rotate = "none", fm = "minres")
## 4因子累积方差解释: 0.68

# 主轴因子法，不旋转
fa.res_pa <- fa(df.use, nfactors = 4, rotate = "none", fm = "pa")
## 4因子累积方差解释: 0.69

# 斜交旋转 (oblimin)
fa.res_oblimin <- fa(df.use, nfactors = 4, rotate = "oblimin", fm = "ml")

# === 结果可视化 ===

# 因子路径图
fa.diagram(fa.res)

# 因子载荷热图（更直观）
pheatmap::pheatmap(t(fa.res$loadings), cluster_rows = FALSE, cluster_cols = FALSE)

# 因子载荷散点图矩阵
factor.plot(fa.res, labels = colnames(df.use), show.points = TRUE)

# === 信度分析（补充：Cronbach's α） ===
# 假定量表各条目已分组，对各维度分别计算信度
# 以因子1（病床利用因子）为例，其高载荷变量为病床利用率、病床周转次数
# alpha_res <- psych::alpha(df.use[, c("病床利用率", "病床周转次数")])
# alpha_res$total$raw_alpha  # 提取 Cronbach's α 系数
```

## 结果解读指南

**因子载荷矩阵（Standardized loadings / Pattern Matrix）：**

- 每一行是一个原始变量，每一列是一个公共因子
- 载荷绝对值越大，说明该变量与该因子的相关性越强（通常 > 0.4 或 > 0.5 视为显著载荷）
- 正载荷 = 正相关，负载荷 = 负相关
- 一个变量在某个因子上载荷高、在其他因子上载荷低 → 该变量"纯净"，因子结构清晰

**公因子方差（h2, communality）：**

- 每个原始变量的方差中被公共因子解释的比例，值域 [0, 1]
- h2 越接近 1，说明该变量被公共因子解释得越充分
- 例：病床周转次数 h2 = 1.00，表示其方差几乎完全被4个公共因子解释
- 病死率 h2 = 0.26，表示公共因子只能解释其 26% 的方差，大部分信息被独特因子解释

**独特性（u2, uniqueness）：**

- u2 = 1 − h2，不能由公共因子解释的方差比例
- u2 越高，该变量与公共因子的关联越弱，在因子分析中的贡献越小

**复杂性（com, complexity）：**

- 反映该变量在多个因子上载荷的分散程度
- com 越接近 1：变量主要与一个因子相关，结构简单
- com 越大：变量与多个因子都有一定关联，交叉载荷明显
- 例：病床利用率 com = 1.0，仅与一个因子相关；平均住院天数 com = 1.8，与多个因子相关

**方差解释表：**

- `SS loadings`：每个因子的载荷平方和（特征值），反映该因子解释的总方差
- `Proportion Var`：该因子解释的方差占总方差的比例
- `Cumulative Var`：前 k 个因子累积解释的方差比例
  - 例：Cumulative Var = 0.68，表示4个因子解释了总方差的 68%
- 累积方差贡献率一般希望 > 60%，> 70% 较理想

**模型拟合指标：**

- `Tucker Lewis Index (TLI)`：一般 > 0.9 表示拟合良好。例中 0.931，拟合可接受
- `RMSEA`：一般 < 0.05 表示接近拟合，< 0.08 表示合理拟合。例中 0.055，拟合尚可
- `BIC`：比较不同模型的指标，值越小越好
- `RMSR`（残差均方根）：越接近 0 拟合越好
- 似然比卡方（Likelihood Chi Square）：p > 0.05 表示模型拟合良好（不拒绝模型）

**因子命名：**

- 因子3（ML3）：出院人数(0.75)、平均住院天数(−0.64)、诊断符合率(−0.56)、抢救成功率(0.70) → **综合因子**
- 因子1（ML1）：病床利用率(0.83)、病床周转次数(0.84) → **病床利用因子**
- 因子2（ML2）：治愈好转率(0.98)、病死率(−0.42) → **医疗水平因子**
- 因子4（ML4）：门诊人次(0.92) → **数量因子**

## 常见问题与注意事项

**Q1：因子分析和主成分分析有什么区别？**

PCA 重在"综合"原始变量信息，提取的主成分是原始变量的线性组合，目标是解释尽可能多的方差。EFA 重在"解释"原始变量间的相关关系，假设存在不可观测的潜在因子支配了变量间的相关性。如果目标是降维 + 计算综合得分，用 PCA；如果目标是探索潜变量结构 + 命名因子，用 EFA。

**Q2：提取几个因子有没有金标准？**

没有绝对标准。一般综合以下依据：
- 碎石图和平行分析的建议（最常用）
- `parameters::n_factors()` 多方法投票结果
- 特征值 > 1 的因子数（Kaiser准则）
- 累积方差贡献率 > 60%-70%
- 因子的专业可解释性（最重要的一条：提取的因子必须有合理的专业意义）

**Q3：正交旋转和斜交旋转怎么选？**

正交旋转（varimax）保持因子间相互独立，结果简洁，解释便利，适合大部分探索性场景。斜交旋转（oblimin）允许因子间相关，在心理学量表、态度问卷等场景中更贴近实际（因为心理特质往往存在自然相关）。如果斜交旋转后因子间相关系数较小（< 0.3），可考虑退回到正交旋转。

**Q4：出现 Heywood 案例（h2 > 1 或 u2 < 0）怎么办？**

这表明因子解异常，常见原因：提取的因子过多、样本量不足、变量间相关性过强。应对措施：减少因子个数、增加样本量、更换提取方法（如改用主轴因子法 pa）、检查变量间是否存在共线性。

**Q5：SPSS 中因子分析默认用的是主成分法，和 R 的 psych 有什么区别？**

SPSS 的"因子分析"模块默认使用主成分提取法（实际上做的是 PCA），而 R 的 `psych::fa()` 默认使用最小残差法（minres）。这是概念上的区别——SPSS 的默认操作究其原理更接近 PCA 而非 EFA。需要用真正的 EFA 时，在 SPSS 中选择"主轴因子法"（principal axis factoring）；在 R 中直接使用 `fm="pa"` 或 `fm="ml"`。

**Q6：如何报告因子分析结果？**

论文中典型报告模式：
> 采用探索性因子分析（最大似然法 + 最大方差旋转）对9个医疗质量指标进行分析。KMO 值为 0.69，Bartlett 球形检验 χ² = 119.03，P < 0.001，说明数据适合因子分析。根据碎石图及平行分析，结合累计方差贡献率，提取4个因子，累计方差贡献率为 68.23%。旋转后的因子载荷矩阵见表X。根据各因子包含的指标内容，将4个因子分别命名为：综合因子、病床利用因子、医疗水平因子和数量因子。

**Q7：因子分析需要多少样本量？**

经验规则：样本量与变量数之比 ≥ 5:1，最少也需要 ≥ 2:1。样本量绝对数一般建议 ≥ 100。但在探索性工作或预实验中，样本量偏少也可进行，只是结果的稳定性较差，结论需要谨慎外推。
