---
name: medical-stat-pca
description: "R语言医学统计：主成分分析（PCA）。用于数据降维、提取主成分、计算主成分得分，通过prcomp()和FactoMineR实现。TRIGGER when user mentions 主成分分析、PCA、数据降维、KMO、Bartlett球形检验、主成分得分、prcomp、特征值、方差贡献率，or asks about dimensionality reduction. SKIP for 因子分析（探索性因子分析）、主成分回归、聚类分析、判别分析。"
---

# 主成分分析 (Principal Component Analysis, PCA)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**什么时候用：**

- 多指标数据中存在较多指标且指标间有较强相关性，希望用少数几个综合指标概括原始指标的主要信息
- 需要消除原始变量之间的多重共线性，提取出相互正交的综合指标
- 需要将高维数据降维到2-3维以便可视化，直观展示样本间的关系
- 医学研究中需要计算综合评分（如生长发育评分、生活质量评分、疾病严重度评分）
- 机器学习中的数据预处理步骤，减少特征数量、加速模型训练

**什么时候不用：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 需要识别不可直接观测的潜变量（潜在因子） | 探索性因子分析（EFA）（`psych::fa()`） |
| 降维后需要用主成分做回归预测因变量 | 主成分回归（PCR）（`pls::pcr()`） |
| 需要将样本分组而不是降维 | 聚类分析（`stats::hclust()` 或 `stats::kmeans()`） |
| 需要对新样本进行分类判别 | 判别分析（`MASS::lda()`） |
| 各变量之间相关性很低，不适合提取共同信息 | 不存在降维的基础，考虑其他分析方法 |

**医学研究常见应用：**

- 儿童生长发育综合评价：身高、体重、胸围、头围、坐高、肺活量等多个指标提取综合发育因子
- 中医证候学：数＋甚至数十个症状指标降维为少数几个证候主成分
- 代谢组学/蛋白质组学：数百个代谢物/蛋白浓度降维为主成分，用于区分疾病状态

## 前置条件

**R包：**

```r
# 基础安装已包含 prcomp()、biplot()、screeplot()
# KMO 和 Bartlett 球形检验 / 综合 PCA 分析
install.packages("psych")
# 数据合成检验
install.packages("performance")
# 主成分个数选择
install.packages("parameters")
# 可视化增强
install.packages(c("see", "ggplot2"))
```

```r
library(psych)
library(performance)
library(parameters)
library(see)
library(ggplot2)
```

**数据格式要求：**

- 所有变量必须是连续型数值变量（不可包含分类/因子变量，`Species` 等分组变量需先排除）
- 数据框格式，每行一个观测，无缺失值
- 变量量纲不同时建议标准化（`scale. = TRUE`）

**统计假设与前提检验：**

- 变量间存在足够的相关性：KMO 检验统计量（Overall MSA）> 0.6 较理想，> 0.5 勉强可用
- Bartlett 球形检验 p < 0.05：拒绝"相关矩阵为单位矩阵"的原假设，表明变量间存在显著相关性
- 各变量最好近似服从多元正态分布（PCA 对正态性不严格，但严重偏态会影响结果）

## 方法选择决策树

```
你的数据情况 →
├── 变量间相关性较强（KMO > 0.6, Bartlett p < 0.05） → 使用 prcomp() 进行 PCA
│   ├── 需要快速实现，变量数适中 → stats::prcomp(scale. = TRUE)
│   └── 需要更丰富的输出和配套绘图函数 → FactoMineR::PCA()
├── 变量间相关性较弱（KMO < 0.5, Bartlett p > 0.05） →
│   └── 数据不适合 PCA，检查变量选择或考虑其他方法（如聚类分析）
├── 含分类变量/分组变量 →
│   └── 先将分类变量排除，仅对连续变量做 PCA；分组变量可用于后续着色可视化
└── 需要进一步用主成分做预测 →
    └── 主成分回归（PCR）：pls::pcr() 或 stats::lm(y ~ ., data = scores)
```

## 标准工作流

### 步骤1：数据准备与探索

```r
# 查看数据结构
str(iris)

# PCA 仅使用连续变量，排除分类变量（Species）
iris_num <- iris[, -5]

# 查看基本描述统计
summary(iris_num)
```

### 步骤2：前提条件检验

```r
# 1. 相关性矩阵：检查变量间是否存在足够的相关性
cor(iris_num)

# 2. KMO 检验：衡量偏相关性是否适合做 PCA
psych::KMO(iris_num)
# 看 Overall MSA：> 0.6 较理想，> 0.5 勉强可用

# 3. Bartlett 球形检验
psych::cortest.bartlett(iris_num)
# p < 0.05 表示适合 PCA

# 以上两步可合并为一行
performance::check_factorstructure(iris_num)
```

### 步骤3：执行 PCA

```r
pca.res <- prcomp(iris_num, scale. = TRUE, center = TRUE)

# 查看完整摘要（标准差、方差贡献率、累积贡献率）
summary(pca.res)
```

### 步骤4：确定保留的主成分个数

```r
# 方法1：碎石图
screeplot(pca.res, type = "lines")

# 方法2：Kaiser 准则（特征值 > 1）
pca.res$sdev^2  # 特征值

# 方法3：parameters 包多方法综合选择
n <- parameters::n_components(iris_num)
n
summary(n)
plot(n) + theme_modern()
```

保留规则：累积方差贡献率达到 70%~90%；或特征值 > 1。二者结合判断，无绝对标准。

### 步骤5：结果报告（论文书写）

在论文中应报告：KMO 值、Bartlett 球形检验 χ² 和 p 值、保留的主成分个数及其理由（特征值 > 1 或累积贡献率 ≥ 某个阈值）、各主成分的特征值与方差贡献率、主成分载荷矩阵、各主成分的命名与解释。

示例："KMO 检验值为 0.54，Bartlett 球形检验 χ² = 706.96，p < 0.001，表明数据适合进行主成分分析。根据 Kaiser 准则（特征值 > 1）和碎石图，保留前 2 个主成分，累计解释总方差的 95.81%。"

## 代码示例

### 完整 PCA 分析流程（以 iris 数据集为例）

```r
library(psych)
library(performance)
library(parameters)
library(see)
library(ggplot2)

# ============================================
# 1. 数据准备
# ============================================
str(iris)
## 'data.frame':    150 obs. of  5 variables:
##  $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
##  $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
##  $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
##  $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
##  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...

iris_num <- iris[, -5]  # 排除分类变量 Species

# ============================================
# 2. 相关性检验
# ============================================
cor(iris_num)
##              Sepal.Length Sepal.Width Petal.Length Petal.Width
## Sepal.Length    1.0000000  -0.1175698    0.8717538   0.8179411
## Sepal.Width    -0.1175698   1.0000000   -0.4284401  -0.3661259
## Petal.Length    0.8717538  -0.4284401    1.0000000   0.9628654
## Petal.Width     0.8179411  -0.3661259    0.9628654   1.0000000

# ============================================
# 3. KMO 和 Bartlett 球形检验
# ============================================
psych::KMO(iris_num)
## Kaiser-Meyer-Olkin factor adequacy
## Call: psych::KMO(r = iris_num)
## Overall MSA =  0.54
## MSA for each item =
## Sepal.Length  Sepal.Width Petal.Length  Petal.Width
##         0.58         0.27         0.53         0.63

psych::cortest.bartlett(iris_num)
## $chisq
## [1] 706.9592
##
## $p.value
## [1] 1.92268e-149
##
## $df
## [1] 6

# 或者二者合并检查
performance::check_factorstructure(iris_num)
## - Sphericity: Bartlett's test suggests sufficient significant correlation
##   (Chisq(6) = 706.96, p < .001).
## - KMO: overall MSA = 0.54.

# ============================================
# 4. 执行 PCA
# ============================================
pca.res <- prcomp(iris_num, scale. = TRUE, center = TRUE)

# ---- 载荷矩阵 ----
pca.res$rotation
##                     PC1         PC2        PC3        PC4
## Sepal.Length  0.5210659 -0.37741762  0.7195664  0.2612863
## Sepal.Width  -0.2693474 -0.92329566 -0.2443818 -0.1235096
## Petal.Length  0.5804131 -0.02449161 -0.1421264 -0.8014492
## Petal.Width   0.5648565 -0.06694199 -0.6342727  0.5235971

# ---- 特征值 ----
pca.res$sdev^2
## [1] 2.91849782 0.91403047 0.14675688 0.02071484

# ---- 主成分得分（前6个样本） ----
head(pca.res$x)
##            PC1        PC2         PC3          PC4
## [1,] -2.257141 -0.4784238  0.12727962  0.024087508
## [2,] -2.074013  0.6718827  0.23382552  0.102662845
## [3,] -2.356335  0.3407664 -0.04405390  0.028282305
## [4,] -2.291707  0.5953999 -0.09098530 -0.065735340
## [5,] -2.381863 -0.6446757 -0.01568565 -0.035802870
## [6,] -2.068701 -1.4842053 -0.02687825  0.006586116

# ---- 方差贡献率与累积贡献率 ----
summary(pca.res)
## Importance of components:
##                           PC1    PC2     PC3     PC4
## Standard deviation     1.7084 0.9560 0.38309 0.14393
## Proportion of Variance 0.7296 0.2285 0.03669 0.00518
## Cumulative Proportion  0.7296 0.9581 0.99482 1.00000

# ============================================
# 5. 结果可视化
# ============================================

# 双标图：同时展示样本得分和变量载荷
biplot(pca.res)

# 碎石图：帮助判断保留几个主成分
screeplot(pca.res, type = "lines")

# ============================================
# 6. 确定最佳主成分个数
# ============================================
n <- parameters::n_components(iris_num)
n
## # Method Agreement Procedure:
##
## The choice of 1 dimensions is supported by 6 (46.15%) methods out of 13
## (Bentler, Optimal coordinates, Acceleration factor, Parallel analysis,
##  Kaiser criterion, Velicer's MAP).

summary(n)
##   n_Factors n_Methods Variance_Cumulative
## 1         1         6           0.7213402
## 2         2         3           0.8667486
## 3         3         3           0.8912973

# 可视化各方法的选择结果
plot(n) + theme_modern()
```

## 结果解读指南

### 1. KMO 与 Bartlett 球形检验

| 输出项 | 含义 | 解读 |
|--------|------|------|
| `Overall MSA` | 总体抽样充分性度量 | 越接近 1 越好，> 0.6 理想；0.5-0.6 勉强可用；< 0.5 不宜 PCA |
| `MSA for each item` | 每个变量的 MSA | 低于 0.5 的变量可考虑剔除 |
| Bartlett `chisq` | 卡方值 | 值越大，变量间相关性越强 |
| Bartlett `p.value` | p 值 | p < 0.05 拒绝单位矩阵原假设，适合 PCA |

### 2. 载荷矩阵 (`$rotation`)

- **载荷的绝对值**：表示原始变量对主成分的贡献程度。绝对值越大，该变量对该主成分越重要
- **载荷的符号**：正号表示正相关，负号表示负相关
- **主成分计算公式**：例如第一主成分 PC1 = 0.521 × Sepal.Length + (-0.269) × Sepal.Width + 0.580 × Petal.Length + 0.565 × Petal.Width

### 3. 特征值与方差贡献率

- **特征值** (`$sdev^2`)：反映每个主成分包含的信息量。特征值越大说明该主成分越重要
- **Proportion of Variance**：方差贡献率，等于该主成分特征值 / 所有特征值之和。表示该主成分解释原始数据总方差的比例
- **Cumulative Proportion**：累积方差贡献率。前 k 个主成分的方差贡献率之和，用于确定保留多少个主成分

### 4. 主成分得分 (`$x`)

- 主成分得分是原始数据投影到主成分方向上的坐标值
- 可用于后续排序、聚类分析、可视化，或作为新变量输入回归模型
- PC1 得分越高或越低，表示该样本在第一主成分代表的特征方向上越极端

### 5. 双标图解读

- **点（样本）**：距离近的点表示在主成分空间特征相似
- **箭头（变量）**：箭头方向表示变量与主成分的相关方向；箭头长度表示变量对主成分的重要性
- 两个箭头夹角接近 0° 表示正相关，接近 180° 表示负相关，接近 90° 表示相关性较弱

## 常见问题与注意事项

**Q1：Standard deviation（标准差）和特征值是什么关系？**

特征值 = 标准差的平方，即 `pca.res$sdev^2`。`summary()` 输出中的 `Standard deviation` 即 `pca.res$sdev`。

**Q2：prcomp() 中 `scale. = TRUE` 什么时候必须用？**

当各变量量纲不同时（如身高 cm vs 体重 kg vs 肺活量 mL），必须标准化使各变量贡献均等。如果各变量量纲相同或认为变量本身变异性大小有意义时，可不标准化。

**Q3：KMO < 0.5 怎么办？**

说明数据不适合 PCA。可检查是否混入了分类变量，或考虑剔除 MSA 值特别低的变量后重试。若仍不理想，数据本身变量间缺乏相关性，不适合降维。

**Q4：主成分分析 vs 因子分析有什么区别？**

PCA 是对原始变量方差的线性变换，不需要假设潜在因子模型；因子分析假设存在不可直接观测的潜在公因子。PCA 目标是最大化方差解释（降维），因子分析目标是揭示数据背后的潜在结构。在医学论文中，二者常被混用，注意区分。

**Q5：R 的 prcomp() 和 SPSS 的 PCA 结果一样吗？**

核心数学相同，但部分软件输出的载荷可能是"载荷平方"（特征向量 × √特征值），R 的 `$rotation` 输出的是特征向量。若需与 SPSS 对比，可用 `pca.res$rotation %*% diag(pca.res$sdev)` 获得与 SPSS 一致的载荷矩阵。

**Q6：`parameters::n_components()` 输出了多种方法的不同结论，听谁的？**

以多数方法一致的选择为参考（`n_Methods` 最大值对应行）。同时结合：
- 累积方差贡献率 ≥ 70%~90%
- Kaiser 准则（特征值 > 1）
- 专业领域可解释性

三种原则综合判断，而非仅看统计量。

**Q7：PCA 需要正态分布假设吗？**

PCA 对变量分布的正态性假设并不严格，因为它是基于协方差（或相关）矩阵的计算。但严重偏态分布可能导致少数异常值主导主成分方向。建议先检查数据分布，必要时做变量变换（如对数变换）。

**Q8：样本量多少合适？**

一般建议样本量至少是变量数的 5~10 倍（如 4 个变量至少需要 20~40 个样本），样本量过小会导致 PCA 结果不稳定。
