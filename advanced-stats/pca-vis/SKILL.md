---
name: medical-stat-pca-vis
description: "R语言医学统计：主成分分析结果可视化。涵盖碎石图（scree plot）、变量载荷图、样本得分图、双标图（biplot）的绘制，使用factoextra和ggplot2包。TRIGGER when user mentions PCA图、主成分图、碎石图、biplot、载荷图、得分图、fviz、3D PCA，or asks about visualizing PCA results. SKIP for 主成分分析计算（仅需prcomp/PCA统计结果）、因子分析可视化。"
---

# 主成分分析可视化 (PCA Visualization)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用PCA可视化的典型场景：**

- 已完成PCA分析，需直观呈现结果：碎石图判断主成分个数、载荷图解释变量-成分关系、得分图展示样本分布
- 探索不同组别样本在主成分空间上的聚类趋势
- 需要同时展示变量载荷和样本分布的双标图（biplot）用于论文发表
- 多变量降维后需三维可视化展示前3个主成分间的关系
- 需自定义美化PCA图以符合期刊配色或样式要求

**不使用PCA可视化的场景：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 未执行PCA分析，需要计算主成分 | `prcomp()` / `FactoMineR::PCA()` |
| 需确定保留主成分个数的统计准则 | 平行分析（`psych::fa.parallel()`）、Kaiser准则 |
| 需降维后用主成分做回归 | 主成分回归（PCR，`pls::pcr()`） |
| 需识别潜在因子结构而非降维展示 | 探索性因子分析（另见因子分析技能） |
| 仅需变量间相关关系热图 | `corrplot` / `pheatmap` |
| 样本数极大（>10万） | 增量PCA或稀疏PCA |

## 前置条件

**R 包安装：**

```r
install.packages(c("factoextra", "FactoMineR", "ggplot2", "corrplot", 
                   "ggsci", "scatterplot3d", "ggpubr"))
```

**数据格式要求：**

- 需提供已完成PCA分析的结果对象，支持两种来源：
  - `FactoMineR::PCA()` 返回的对象（推荐，包含 cos2/contrib 等额外统计结果）
  - `prcomp()` / `princomp()` 返回的对象（R基础PCA）
- 变量应为连续型数值变量，不含缺失值
- 如需按组上色或分面，需额外提供分组变量（因子型）
- 数据建议标准化（`scale.unit = TRUE` 或 `scale. = TRUE`）

## 方法选择决策树

```
你的可视化需求 →
├── 选择主成分个数 →
│   ├── 方差解释度视图 → 碎石图 fviz_eig(pca.res, addlabels = TRUE)
│   └── 查看数值结果 → get_eigenvalue(pca.res)
│
├── 展示变量与主成分的关系 →
│   ├── 变量投影方向图 → fviz_pca_var(pca.res)
│   ├── 变量相关性强弱 → fviz_cos2(pca.res, choice = "var", axes = 1:2)
│   ├── 变量贡献大小 → fviz_contrib(pca.res, choice = "var", axes = 1)
│   └── 矩阵形式展示 → corrplot(res.var$cos2, is.corr = FALSE)
│
├── 展示样本在主成分空间的分布 →
│   ├── 不分组的散点图 → fviz_pca_ind(pca.res)
│   ├── 按分组上色 + 置信椭圆 → fviz_pca_ind(pca.res, col.ind = group, addEllipses = TRUE)
│   ├── 按cos2映射颜色/大小 → fviz_pca_ind(pca.res, col.ind = "cos2", pointsize = "cos2")
│   └── 按contrib映射 → fviz_pca_ind(pca.res, col.ind = "contrib", pointsize = "contrib")
│
├── 同时展示变量和样本 →
│   ├── 基础双标图 → fviz_pca_biplot(pca.res, col.ind = group)
│   ├── 自定义变量颜色按分组 → fviz_pca_biplot(pca.res, col.var = factor(var_group))
│   └── 精细化控制（填充/边框/透明度） → fviz_pca_biplot() + ggpubr::fill_palette()
│
├── 需要ggplot2底层自由控制 →
│   ├── 提取PC得分 → pca.res$x（prcomp）或 get_pca_ind(pca.res)$coord
│   ├── 提取变量载荷 → pca.res$rotation（prcomp）或 get_pca_var(pca.res)$coord
│   └── ggplot2自由绘制 → ggplot(df, aes(PC1, PC2)) + geom_point() + stat_ellipse()
│
└── 需要三维可视化 →
    └── scatterplot3d(tmp[, 1:3], color = group_color, pch = 15)
```

## 标准工作流

### 步骤 1：执行 PCA 分析

使用 `FactoMineR::PCA()` 完成分析，`graph = FALSE` 关闭默认图形，由 factoextra 接管可视化。

```r
rm(list = ls())
library(factoextra)
library(FactoMineR)

pca.res <- PCA(iris[, -5], graph = FALSE, scale.unit = TRUE)
pca.res
## **Results for the Principal Component Analysis (PCA)**
## The analysis was performed on 150 individuals, described by 4 variables
```

PCA结果对象包含 `$eig`（特征值）、`$var`（变量结果：coord/cor/cos2/contrib）、`$ind`（样本结果：coord/cos2/contrib）。

### 步骤 2：特征值提取与碎石图

```r
get_eigenvalue(pca.res)
##       eigenvalue variance.percent cumulative.variance.percent
## Dim.1 2.91849782       72.9624454                    72.96245
## Dim.2 0.91403047       22.8507618                    95.81321
## Dim.3 0.14675688        3.6689219                    99.48213
## Dim.4 0.02071484        0.5178709                   100.00000

fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 100))
```

前2个主成分累积贡献率达 95.81%。通常保留特征值 > 1 或累积贡献率 > 80% 的主成分。

### 步骤 3：变量结果提取

```r
res.var <- get_pca_var(pca.res)

res.var$cor        # 变量与主成分的相关系数（标准化后等于coord）
## Petal.Length 与 Dim.1 相关系数 0.992，Sepal.Width 与 Dim.2 相关系数 0.883

res.var$cos2       # 坐标平方（表示质量），同变量所有主成分之和为1
## Petal.Length/Dim.1 = 0.983, Sepal.Width/Dim.2 = 0.779

res.var$contrib    # 贡献度（%），同主成分上所有变量之和为100
## Dim.1: Petal.Length 33.7%, Petal.Width 31.9%, Sepal.Length 27.2%, Sepal.Width 7.3%
## Dim.2: Sepal.Width 85.2%, Sepal.Length 14.2%, 其余 < 1%
```

`coord` 是变量在主成分上的投影坐标。如 `Sepal.Width` 在 Dim.1 = -0.460、Dim.2 = 0.883，由此坐标确定箭头方向。

### 步骤 4：变量结果可视化

```r
# 变量载荷箭头图（基础）
fviz_pca_var(pca.res)

# 按cos2着色
fviz_pca_var(pca.res, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# 按contrib着色
fviz_pca_var(pca.res, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# 黑白版：透明度映射cos2
fviz_pca_var(pca.res, alpha.var = "cos2")

# 变量cos2条形图
fviz_cos2(pca.res, choice = "var", axes = 1:2)

# 变量contrib条形图（指定主成分）
fviz_contrib(pca.res, choice = "var", axes = 1)
fviz_contrib(pca.res, choice = "var", axes = 1:2)

# corrplot矩阵形式
library("corrplot")
corrplot(res.var$cos2, is.corr = FALSE)
corrplot(res.var$contrib, is.corr = FALSE)
```

解读：Dim.1 主要代表花瓣指标（Petal.Length + Petal.Width），Dim.2 主要代表花萼宽度（Sepal.Width）。

### 步骤 5：维度描述

```r
res.desc <- dimdesc(pca.res, axes = c(1, 2), proba = 0.05)
res.desc$Dim.1
##              correlation       p.value
## Petal.Length   0.9915552 3.369916e-133
## Petal.Width    0.9649790  6.609632e-88
## Sepal.Length   0.8901688  2.190813e-52
## Sepal.Width   -0.4601427  3.139724e-09
```

`dimdesc()` 展示每个主成分维度与原始变量的相关性及显著性检验，用于给主成分赋予生物学/医学命名。

### 步骤 6：样本得分图

```r
res.ind <- get_pca_ind(pca.res)
head(res.ind$coord)  # 样本在各主成分的得分（score）

# 基本分布图
fviz_pca_ind(pca.res)

# 经典样式：按组着色 + 置信椭圆
fviz_pca_ind(pca.res,
             geom.ind = "point",
             col.ind = iris$Species,
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE,
             legend.title = "Groups")

# cos2映射颜色
fviz_pca_ind(pca.res, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# cos2映射点大小
fviz_pca_ind(pca.res, pointsize = "cos2", pointshape = 21, fill = "#E7B800", repel = TRUE)

# cos2颜色 + contrib大小
fviz_pca_ind(pca.res, col.ind = "cos2", pointsize = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# 样本的cos2和contrib条形图
fviz_cos2(pca.res, choice = "ind", axes = 1:2)
fviz_contrib(pca.res, choice = "ind", axes = 1:2)
```

### 步骤 7：双标图 (Biplot)

同时展示变量箭头和样本分布，是论文中最常见的主成分可视化形式。

```r
# 基础双标图
fviz_pca_biplot(pca.res, 
                col.ind = iris$Species, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species")

# 高级：分别控制样本填充色和变量颜色
fviz_pca_biplot(pca.res, 
                geom.ind = "point", pointshape = 21, pointsize = 2.5,
                fill.ind = iris$Species, col.ind = "black",
                col.var = factor(c("sepal", "sepal", "petal", "petal")),
                legend.title = list(fill = "Species", color = "Clusters"),
                repel = TRUE) +
  ggpubr::fill_palette("jco") + ggpubr::color_palette("npg")

# 样本按物种着色，变量按贡献着色+透明度
fviz_pca_biplot(pca.res, 
                geom.ind = "point", fill.ind = iris$Species,
                col.ind = "black", pointshape = 21, pointsize = 2,
                palette = "jco", addEllipses = TRUE,
                alpha.var = "contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                legend.title = list(fill = "Species", color = "Contrib", alpha = "Contrib"))
```

### 步骤 8：ggplot2 自定义绘制

当 factoextra 默认样式不满足需求时，可手动提取数据用 ggplot2 绘制：

```r
# prcomp 计算 PCA
pca.res <- prcomp(iris[, -5], scale. = TRUE, center = TRUE)

# 提取PC得分并与分组信息合并
tmp <- as.data.frame(pca.res$x)
tmp$species <- iris$Species

# ggplot2 自定义绘图
library(ggplot2)
library(ggsci)

ggplot(tmp, aes(PC1, PC2)) +
  geom_point(aes(color = species)) +
  stat_ellipse(aes(fill = species), alpha = 0.2, geom = "polygon", type = "norm") +
  scale_fill_aaas() + scale_color_aaas() + theme_bw()
```

### 步骤 9：3D PCA 图

```r
library(scatterplot3d)

scatterplot3d(tmp[, 1:3],
              color = rep(c("#00AFBB", "#E7B800", "#FC4E07"), each = 50),
              pch = 15, lty.hide = 2)
legend("topleft", c('Setosa', 'Versicolor', 'Virginica'),
       fill = c("#00AFBB", "#E7B800", "#FC4E07"), box.col = NA)
```

## 代码示例

完整可复现的PCA可视化流程：

```r
library(factoextra)
library(FactoMineR)
library(corrplot)
library(ggplot2)
library(ggsci)
library(scatterplot3d)

pca.res <- PCA(iris[, -5], graph = FALSE, scale.unit = TRUE)

fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_var(pca.res, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_contrib(pca.res, choice = "var", axes = 1:2)

fviz_pca_ind(pca.res, geom.ind = "point", col.ind = iris$Species,
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, legend.title = "Groups")

fviz_pca_biplot(pca.res, col.ind = iris$Species, palette = "jco",
                addEllipses = TRUE, label = "var", col.var = "black",
                repel = TRUE, legend.title = "Species")

# ggplot2版本
tmp <- as.data.frame(prcomp(iris[, -5], scale. = TRUE, center = TRUE)$x)
tmp$species <- iris$Species
ggplot(tmp, aes(PC1, PC2)) +
  geom_point(aes(color = species)) +
  stat_ellipse(aes(fill = species), alpha = 0.2, geom = "polygon", type = "norm") +
  scale_fill_aaas() + scale_color_aaas() + theme_bw()

# 3D版本
scatterplot3d(tmp[, 1:3], pch = 15, lty.hide = 2,
              color = rep(c("#00AFBB", "#E7B800", "#FC4E07"), each = 50))
legend("topleft", c('Setosa', 'Versicolor', 'Virginica'),
       fill = c("#00AFBB", "#E7B800", "#FC4E07"), box.col = NA)
```

## 结果解读指南

**碎石图**：纵轴为方差百分比，横轴为主成分序号。曲线趋于平缓的"拐点"即建议保留的主成分个数。常配合特征值 = 1 的 Kaiser 准则线。

**变量载荷图**：箭头方向代表原始变量投影方向，箭头越长表示该变量对两个主成分的解释度越高。同向为正相关，正交近似独立，反向为负相关。

**变量贡献图**：条形图展示各变量对指定主成分的贡献百分比。红色虚线标识期望平均贡献（100 / 变量数 %），高于该线的变量是该主成分的主导变量。

**样本得分图**：每个点为一样本，距离越近表示在多维特征空间中越相似。置信椭圆（默认 95% 水平）可判断组间分离趋势，但不等于统计检验。

**双标图**：同时展示样本点和变量箭头。样本相对变量箭头的投影位置解释该样本在该变量方向上的相对大小。

**维度描述**：列出与各主成分显著相关的原始变量，用于给主成分赋予含义命名。

## 常见问题与注意事项

**cos2 和 contrib 的区别？**

`cos2`（余弦平方）表示变量被该主成分捕获的变异比例，同一变量所有主成分 cos2 之和为 1。`contrib`（贡献度）表示变量对该主成分构建的贡献百分比，同主成分所有变量 contrib 之和为 100。

**变量载荷图坐标和相关系数为什么一样？**

使用 `scale.unit = TRUE` 标准化后，`res.var$coord` 等于 `res.var$cor`。两者数值相等但解释侧重不同：coord 用于构图，cor 用于统计解读。

**PCA() 和 prcomp() 结果在 factoextra 中有什么区别？**

`PCA()` 返回的对象更丰富（含 cos2、contrib），fviz 函数自动提取。`prcomp()` 也可直接使用，但 `addEllipses` 等高级参数需配合额外变量。推荐统一使用 `PCA()`。

**fviz 系列函数支持哪些自定义？**

底层基于 ggpubr 的 `ggscatter`，支持所有 ggplot2 图层叠加（`+ theme_bw()` 等）和 ggpubr 调色板（`fill_palette("jco")`），图表可通过 `ggsave()` 导出。

**分组变量没有在 PCA 计算中纳入，可视化时如何加上？**

分组变量不应纳入 PCA 计算。计算时仅纳入数值变量（`iris[, -5]`），可视化时通过 `col.ind = iris$Species` 或 `fill.ind = group_var` 映射给颜色/形状。

**置信椭圆代表什么？可以解读为组间差异显著吗？**

椭圆形基于多元正态假设的 95% 置信区域。椭圆不重叠不等于组间差异不显著（PCA 不能直接做显著性检验），统计推断需用 MANOVA 或 PERMANOVA。

**参考资料**：http://www.sthda.com/ 提供 factoextra 包的完整文档和图库。
