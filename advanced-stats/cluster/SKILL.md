---
name: medical-stat-cluster
description: "R语言医学统计：聚类分析。涵盖系统聚类（层次聚类）、K-means聚类、变量聚类，用于发现数据中的自然分组或对变量进行降维分组。TRIGGER when user mentions 聚类分析、层次聚类、K-means、系统聚类、hclust、kmeans、PAM、pam，or asks about grouping subjects or variables. SKIP for 判别分析（有分类标签）、主成分分析。"
---

# 聚类分析 (Cluster Analysis)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

聚类分析用于将观察对象或变量按照其相似性进行分组，使得同一组内的对象尽可能相似，不同组间的对象尽可能不同。与判别分析不同，聚类分析不需要预先知道分类标签，属于无监督学习方法。

**适用条件：**
- 需要对样本进行自然分组，如根据患者的临床指标将患者分为不同亚型
- 探索性数据分析，发现数据中可能存在的内在结构
- 变量降维分组，将多个相关变量聚合成若干变量簇，每簇选择一个代表性变量

**不适用场景：**
- 已有明确的分类标签，需要建立分类规则 → 使用判别分析（`lda()` / `qda()`）
- 需要对变量进行线性组合降维 → 使用主成分分析（`prcomp()` / `princomp()`）
- 数据量极大（>10万条记录）→ K-means 可能收敛缓慢，考虑 MiniBatchKMeans 或抽样后聚类

**医学研究常见应用：**
- 根据基因表达谱将肿瘤患者分为不同分子亚型
- 根据实验室检查指标将疾病表型分类
- 药物化合物结构分类
- 居民膳食模式分类

## 前置条件

**R包安装：**

```r
install.packages(c("flexclust", "NbClust", "factoextra", "cluster", "rattle", "psych"))
```

**数据格式要求：**

- 数据必须为数值型矩阵或数据框（因子变量需先转换为哑变量）
- 行代表观测对象，列代表聚类变量
- 聚类前需进行标准化（`scale()`），消除量纲影响。不同变量量纲差异大时，未标准化的距离完全由数值最大的变量主导
- 缺失值需提前处理（删除或插补），`hclust()` 和 `kmeans()` 不接受含 `NA` 的数据

**统计假设：**

- 层次聚类：无严格分布假设，但对异常值敏感，建议标准化后使用
- K-means 聚类：假设各类呈球形分布，各组方差相近；对异常值敏感
- PAM 聚类：无分布假设，对异常值较稳健，可接受混合数据类型

**距离与连接方法选择：**

距离计算（`dist()` 的 `method` 参数）：
- `"euclidean"`：欧几里得距离，最常用，适用于连续型变量
- `"manhattan"`：曼哈顿距离，对异常值相对不敏感
- `"maximum"`：切比雪夫距离

层次聚类的连接方法（`hclust()` 的 `method` 参数）：
- `"average"`：平均连接法，对异常值相对稳健，推荐
- `"complete"`：完全连接法，倾向于产生紧凑的类
- `"single"`：单一连接法，易产生链状效应，不推荐
- `"ward.D2"`：离差平方和法，倾向于产生大小均匀的类

## 方法选择决策树

```
你的数据情况 →
├── 样本量小（n < 150），希望可视化聚类过程、不确定聚类个数
│   └── 使用系统聚类（层次聚类）hclust() + NbClust确定个数 + cutree()切分
│
├── 样本量大（n ≥ 150），变量为连续型，假设各类呈球形
│   ├── 能事先确定聚类个数 → 直接使用 kmeans()
│   └── 不能确定聚类个数 → NbClust(method="kmeans") 或 fviz_nbclust() 确定个数后使用 kmeans()
│
├── 数据可能存在异常值，或包含混合数据类型，或对稳健性要求高
│   └── 使用 PAM 聚类 cluster::pam()
│
└── 需要对变量而非样本进行聚类（变量降维）
    └── 使用系统聚类（对转置后的数据矩阵）或 varclus()（Hmisc包）
```

> K-means 需要预设聚类个数 K。如无法确定 K，使用 `NbClust()` 或 `factoextra::fviz_nbclust()` 通过肘部法、CH 指数等多种准则判断最佳 K 值。

## 标准工作流

### 步骤1：数据准备与探索

- 去除 ID 列、分类标签列等非分析变量
- 处理缺失值（`na.omit()` 删除或用均值/中位数插补）
- 检查各变量分布，识别潜在的异常值
- 对数据进行标准化：`scale(data)`

### 步骤2：确定最佳聚类个数

- 使用 `NbClust::NbClust()` 综合多种评判准则
- 层次聚类可结合树状图肉眼判断切割高度
- K-means 聚类使用肘部法（`factoextra::fviz_nbclust()`）

### 步骤3：执行聚类分析

- 层次聚类：`hclust(dist(std_data), method = "average")` → `cutree(hc, k = best_k)`
- K-means：`kmeans(std_data, centers = best_k, nstart = 25)`（`nstart` 建议 ≥ 25，多次随机初始值取最优解）
- PAM：`cluster::pam(data, k = best_k, stand = TRUE)`

### 步骤4：结果可视化

- 层次聚类：`plot(hc)` + `rect.hclust(hc, k = best_k)`
- 通用：`factoextra::fviz_cluster(fit, data = std_data, ellipse = TRUE, ellipse.type = "t")`

### 步骤5：结果报告（论文写作）

- 报告聚类方法（层次/K-means/PAM）、所用距离/连接方式、标准化方式
- 报告各聚类个数及每类样本量
- 报告各类的变量均值（原始尺度）以描述各类特征
- 层次聚类可附聚类树状图；K-means/PAM 可附二维聚类投影图

## 代码示例

### 示例1：层次聚类（系统聚类）

数据来自 `flexclust` 包中的 `nutrient` 数据集，包含 27 种食物（牛肉、鱼肉等）的 5 个营养成分指标。

```r
library(flexclust)
library(NbClust)

data(nutrient, package = "flexclust")
row.names(nutrient) <- tolower(row.names(nutrient))

dim(nutrient)
## [1] 27  5

psych::headTail(nutrient)
##                 energy protein fat calcium iron
## beef braised       340      20  28       9  2.6
## hamburger          245      21  17       9  2.7
## beef roast         420      15  39       7    2
## beef steak         375      19  32       9  2.6
## ...                ...     ... ...     ...  ...
## salmon canned      120      17   5     159  0.7
## sardines canned    180      22   9     367  2.5
## tuna canned        170      25   7       7  1.2
## shrimp canned      110      23   1      98  2.6

# 标准化后聚类
nutrient.scaled <- scale(nutrient)

# 欧几里得距离 + 平均连接法
h.clust <- hclust(dist(nutrient.scaled, method = "euclidean"),
                  method = "average")

plot(h.clust, hang = -1, main = "层次聚类", sub = "", xlab = "",
     cex.lab = 1.0, cex.axis = 1.0, cex.main = 2)

# 确定最佳聚类个数
nc <- NbClust(nutrient.scaled, distance = "euclidean",
              min.nc = 2, max.nc = 10,
              method = "average")
## *** : 根据多数原则，最佳聚类数目为 2

# 查看各准则投票情况
barplot(table(nc$Best.nc[1, ]),
        xlab = "聚类数目",
        ylab = "评判准则个数")

# 划分为 5 类
cluster <- cutree(h.clust, k = 5)
table(cluster)
## cluster
##  1  2  3  4  5
##  7 16  1  2  1

# 绘制带分类标识的聚类树
plot(h.clust, hang = -1, main = "", xlab = "")
rect.hclust(h.clust, k = 5)
```

### 示例2：结合 NbClust 和 fviz_nbclust 确定 K 值（K-means）

数据来自 `rattle` 包中的 `wine` 数据集，178 种葡萄酒，13 个化学成分变量。

```r
library(rattle)
library(NbClust)
library(factoextra)

data(wine, package = "rattle")
df <- scale(wine[, -1])  # 去掉第一列 Type 标签

# 方法A：NbClust 综合评判
set.seed(123)
nc <- NbClust(df, min.nc = 2, max.nc = 15, method = "kmeans")
## * Among all indices:
## * 2 proposed 2 as the best number of clusters
## * 19 proposed 3 as the best number of clusters
## * According to the majority rule, the best number of clusters is 3

barplot(table(nc$Best.nc[1, ]),
        xlab = "聚类数目", ylab = "评判准则个数")

# 方法B：factoextra 肘部法
set.seed(123)
fviz_nbclust(df, kmeans, k.max = 15)  # CH指数＋肘部法，建议K=3
```

### 示例3：执行 K-means 聚类

```r
# K-means 聚类，nstart=25 多次随机初始值取最优
set.seed(123)
fit.km <- kmeans(df, centers = 3, nstart = 25)
fit.km
## K-means clustering with 3 clusters of sizes 51, 62, 65
##
## Cluster means:
##      Alcohol      Malic        Ash Alcalinity   Magnesium     Phenols
## 1  0.1644436  0.8690954  0.1863726  0.5228924 -0.07526047 -0.97657548
## 2  0.8328826 -0.3029551  0.3636801 -0.6084749  0.57596208  0.88274724
## 3 -0.9234669 -0.3929331 -0.4931257  0.1701220 -0.49032869 -0.07576891
##
## Within cluster sum of squares by cluster:
## [1] 326.3537 385.6983 558.6971
##  (between_SS / total_SS =  44.8 %)

# 聚类可视化
fviz_cluster(fit.km, data = df,
             ellipse = TRUE,
             ellipse.type = "t",
             geom = "point",
             palette = "lancet",
             ggtheme = theme_bw())
```

### 示例4：PAM 聚类（围绕中心点的划分）

```r
library(cluster)

set.seed(123)
fit.pam <- pam(wine[-1], k = 3, stand = TRUE)
fit.pam
## Medoids:
##      ID Type Alcohol Malic  Ash Alcalinity Magnesium Phenols Flavanoids ...
## 36   35    1   13.48  1.81 2.41       20.5       100    2.70       2.98 ...
## 107 106    2   12.25  1.73 2.12       19.0        80    1.65       2.03 ...
## 149 148    3   13.32  3.24 2.38       21.5        92    1.93       0.76 ...

# 聚类可视化
clusplot(fit.pam, main = "PAM cluster")

# 或使用 factoextra 美化
fviz_cluster(fit.pam,
             ellipse = TRUE,
             ellipse.type = "t",
             geom = "point",
             palette = "aaas",
             ggtheme = theme_bw())
```

### 示例5：层次聚类树状图精细定制

```r
dhc <- as.dendrogram(h.clust)

# 给标签按分类添加颜色
clusMember <- cutree(h.clust, k = 5)
labelColors <- c("#CDB380", "#036564", "#EB6841", "#EDC951", "#487AA1")

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar,
                            list(cex = 1.5, pch = 20,
                                 col = labCol,
                                 lab.col = labCol,
                                 lab.font = 2, lab.cex = 1))
  }
  n
}

diyDendro <- dendrapply(dhc, colLab)
plot(diyDendro, main = "DIY Dendrogram")
legend("topright",
       legend = paste("Cluster", 1:5),
       col = labelColors,
       pch = 20, bty = "n", pt.cex = 2, cex = 1,
       text.col = "black")
```

## 结果解读指南

**层次聚类树状图：**
- 纵轴（Height）表示合并时的距离，高度越高说明被合并的两类差异越大
- 横轴为各观测的标签排列，`hang = -1` 将标签对齐在底部方便读取
- `rect.hclust()` 添加的红色矩形框标示了各类的范围
- 切割高度通过 `k` 参数（指定聚类数）隐式确定，也可用 `h` 参数直接指定高度

**K-means 结果组件：**
- `centers` / `Cluster means`：各类在各变量上的均值（标准化后）。正值表示该变量高于总体均值，负值表示低于总体均值。通过比较各类的聚类中心，可描述各类的特征分布
- `size`：各类的样本量
- `withinss`：各类的组内平方和，值越小同类对象越紧凑
- `totss`：总平方和；`betweenss / totss` 表示组间差异占总差异的比例，通常 40%~80% 较为理想。过低说明聚类效果差，过高可能过拟合
- `tot.withinss`：总组内平方和，越小越好

**PAM 结果组件：**
- `medoids`：各类的中心点（medoid），指的是最能代表该类的一个实际观测，而非虚拟的均值点。其各变量值为原始尺度
- `silinfo`：轮廓系数（silhouette），用于评价聚类质量，取值 [-1, 1]。接近 1 表示聚类效果好（同类紧凑、异类分离），接近 0 表示处于类边界，负值表示可能错分

**论文中统计描述示例：**
> 采用系统聚类法（欧几里得距离、平均连接法）对 27 种食物的营养成分（能量、蛋白质、脂肪、钙、铁）进行聚类分析。经 NbClust 包多准则评判，最优聚类数为 5 类。第 1 类（n=7）为高钙低能量食物，第 2 类（n=16）为中等能量中等蛋白质食物……聚类树状图见图 1。

## 常见问题与注意事项

**Q1：为什么聚类前要标准化？**
不同变量量纲差异大时，未标准化的距离完全由数值最大的变量主导。例如能量（范围 100-400）和铁含量（范围 0.5-3.0），不做标准化时能量将主导聚类结果。使用 `scale()` 后各变量均值为 0、标准差为 1，处于同等权重。

**Q2：层次聚类 vs K-means 如何选择？**
- 不确定聚类个数、想要可视化聚类过程 → 层次聚类
- 数据量大（>150 条记录）、追求速度 → K-means
- 实际工作中两者可以结合：先用层次聚类确定聚类个数，再用 K-means 重新聚类取长补短

**Q3：K-means 中 `nstart` 参数的作用？**
K-means 的结果依赖于初始聚类中心的选择。`nstart` 越大（建议 ≥ 25），算法尝试的初始中心数越多，越容易找到全局最优解。`set.seed()` 确保结果可复现。

**Q4：K-means 和 PAM 的区别？**
- K-means 使用欧几里得距离，基于各类均值；PAM 可以使用任意距离，基于实际观测（medoid）
- K-means 对异常值敏感（一个极端值会大幅移动均值）；PAM 更稳健
- PAM 可处理混合数据类型，但计算量更大

**Q5：聚类结果稳定吗？如何验证？**
聚类结果受随机初始值影响。建议：
- 设 `set.seed()` 保证复现
- 多次运行比较，看分组是否一致
- 使用 `NbClust` 多准则交叉验证最佳聚类数
- 对 K-means 加大 `nstart`（≥ 25）
- 用轮廓系数（silhouette）评估聚类质量

**Q6：异常值如何处理？**
- 层次聚类对异常值敏感：聚类前用箱线图检查，考虑 Winsorize 或移除
- K-means 对异常值敏感：改用 PAM 聚类
- 异常值可能自身构成一类，不一定需要移除

**Q7：如何报告聚类结果？**
论文中建议报告以下内容：
- 聚类目的、所用方法（距离、连接法/算法）、标准化方式
- 确定聚类个数的方法与依据（如 NbClust 的多数原则、肘部法）
- 各类的样本量及在各变量上的均数 ± 标准差（便于读者理解类特征）
- 聚类质量评价指标（组间平方和比例、轮廓系数等）
- 可视化结果（树状图或聚类投影图）
