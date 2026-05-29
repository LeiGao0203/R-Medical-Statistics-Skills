---
name: medical-stat-nonparametric
description: "R语言医学统计：非参数检验。涵盖Wilcoxon符号秩检验（配对）、Mann-Whitney U检验（两独立样本）、Kruskal-Wallis H检验（多组）、Friedman检验（随机区组）。TRIGGER when user mentions 非参数检验、Wilcoxon、Mann-Whitney、Kruskal-Wallis、秩和检验，or when data violates normality assumptions. SKIP for 正态分布数据（用t检验/ANOVA）。"
---

# 秩转换的非参数检验 (Nonparametric Tests)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用非参数检验：**

- 数据不满足正态性假设（经 Shapiro-Wilk 检验 p ≤ 0.05）且无法通过变量转换改善
- 结局变量为等级资料（如疗效：治愈/显效/有效/无效）或开口资料（如 <0.01）
- 样本量较小（n < 30）且总体分布未知
- 数据存在明显的极端值/离群值，参数方法不稳健
- 配对设计差值不满足正态分布 → Wilcoxon符号秩检验
- 两组独立样本数据不满足正态或方差不齐 → Mann-Whitney U检验
- 多组独立样本（完全随机设计）数据不满足正态或方差齐 → Kruskal-Wallis H检验
- 随机区组设计（重复测量）数据不满足球对称假设 → Friedman M检验

**何时不使用非参数检验：**

| 情况 | 应使用的替代方法 |
|------|----------------|
| 数据满足正态性和方差齐性 | 配对t检验 / 独立样本t检验 / 单因素方差分析 |
| 两分类变量关联性分析 | 卡方检验 / Fisher精确检验 |
| 需要控制协变量 | 秩转换后的协方差分析（Quade's ANCOVA或art包） |
| 大样本且正态性轻微偏离 | 参数检验具有一定稳健性，可酌情使用 |

医学研究中的常见应用：两种治疗方法疗效的等级比较、不同药物组实验室指标（偏态数据）的比较、不同时间点重复测量等级资料的比较。

## 前置条件

**R 包依赖：**

```r
# 基础包（无需安装）
# stats —— wilcox.test(), kruskal.test(), friedman.test(), pairwise.wilcox.test()

# 第三方包（多重比较需要）
install.packages("PMCMRplus")  # Nemenyi检验（K-W后多重比较）、Quade检验
```

**非参数检验的前提条件：**

非参数检验对总体分布没有要求（分布自由，distribution-free），但仍有一些基本要求：

1. **独立性**：除配对/区组设计外，各观测值相互独立
2. **随机抽样**：样本来自研究总体的随机抽样
3. **分布形状**：Wilcoxon秩和检验（两独立样本）假设两组分布形状相似，仅在位置上可能有差异；若分布形状不同，检验的实际上是随机优势（stochastic superiority）而非中位数移位
4. **数据类型**：至少为顺序尺度（ordinal）——能排出大小顺序即可

**数据格式要求：**

- **配对设计**：两个等长数值向量，每个受试者两种条件下的测量值一一对应
- **两独立样本**：格式一——两个数值向量分别存放两组数据；格式二——一列数值结局变量 + 一列分组变量（长格式）
- **多组独立样本**（Kruskal-Wallis）：一列数值结局变量 + 一列多水平分组变量
- **随机区组设计**（Friedman）：矩阵格式，每行是一个区组，每列是一种处理

## 方法选择决策树

```
你的数据和设计情况 →
│
├── 配对设计（同一受试者前后测量、自身配对）
│   ├── 差值满足正态性 → 配对t检验（t.test(x1, x2, paired = TRUE)）
│   └── 差值不满足正态性 → Wilcoxon符号秩检验 → wilcox.test(x1, x2, paired = TRUE)
│
├── 两独立样本（完全随机分组）
│   ├── 满足正态性+方差齐性 → 独立样本t检验（t.test(x ~ group, var.equal = TRUE)）
│   └── 不满足正态性 → Wilcoxon秩和检验（Mann-Whitney U） → wilcox.test(x ~ group)
│
├── 多组独立样本（完全随机设计，≥3组）
│   ├── 满足正态性+方差齐性 → 单因素方差分析（aov(x ~ group)）
│   └── 不满足前提条件 → Kruskal-Wallis H检验 → kruskal.test(x ~ group)
│       └── 若整体有统计学意义 → Nemenyi检验（事后多重比较）
│
└── 多组相关样本（随机区组设计/重复测量，≥3个时间点）
    ├── 满足球对称假设 → 重复测量方差分析
    └── 不满足球对称 → Friedman M检验 → friedman.test(矩阵)
        └── 若整体有统计学意义 → Quade检验（事后两两比较）
```

**配对设计、独立样本与重复测量设计的快速区分：**

- 配对设计 = 每个受试者贡献2个测量值（治疗前后、左右侧等）
- 两独立样本 = 每个受试者只属于一个组，贡献1个值
- 随机区组设计 = 每个区组接受所有处理（≥3种），每种处理1个测量值

## 标准工作流

### 步骤1：数据准备与探索

读入数据，查看数据结构，绘制箱线图初步观察各组分布情况。数据格式见「前置条件」。

### 步骤2：前提条件检验

对原始数据或差值（配对设计）进行正态性检验（`shapiro.test()`）。若 p > 0.05 满足正态性，应优先使用参数检验（t检验/ANOVA）因为参数检验效能更高。若 p ≤ 0.05 或数据本身为等级资料，使用本章的非参数方法。

### 步骤3：执行统计分析

根据研究设计选择对应的非参数检验（见决策树）。注意 `wilcox.test()` 中 `paired = TRUE` 用于配对设计，`paired = FALSE`（默认）用于独立样本。对于存在相同秩次（ties）的数据，R会自动进行校正。

### 步骤4：多重比较（若多组检验显著）

Kruskal-Wallis H检验或Friedman M检验整体显著后，需要进行事后两两比较以确定具体哪些组间存在差异。使用 `PMCMRplus` 包中的 Nemenyi检验（独立样本）或 Quade检验（区组设计）。

### 步骤5：结果报告

论文中报告格式示例："三组药物死亡率比较，Kruskal-Wallis H检验结果显示差异有统计学意义（χ² = 9.74, df = 2, p = 0.008）。Nemenyi检验两两比较结果显示，Drug_C组死亡率低于Drug_A组（q = 3.998, p = 0.013）。"

## 代码示例

### 6.1 配对样本比较的Wilcoxon符号秩检验

Wilcoxon符号秩检验（signed-rank test）是对配对差值进行秩次分析，是配对t检验的非参数替代。

使用课本**例8-1**的数据——12名受试者两种处理的测量值：

```r
test1 <- c(60, 142, 195, 80, 242, 220, 190, 25, 198, 38, 236, 95)
test2 <- c(76, 152, 243, 82, 240, 220, 205, 38, 243, 44, 190, 100)

boxplot(test1, test2)
```

```r
wilcox.test(test1, test2, paired = TRUE, alternative = "two.sided", exact = FALSE)
## 
##  Wilcoxon signed rank test with continuity correction
## 
## data:  test1 and test2
## V = 11.5, p-value = 0.06175
## alternative hypothesis: true location shift is not equal to 0
```

结果与课本一致，p = 0.062 > 0.05，差异无统计学意义。

> **R 4.4.0+ 注意**：`wilcox.test()` 在 R 4.4.0 及以后的版本中，formula 方法不再支持 `paired` 参数。若使用 `wilcox.test(x ~ group, paired = TRUE)` 会报错 `cannot use 'paired' in formula method`，请改用两向量形式 `wilcox.test(x1, x2, paired = TRUE)`。

### 6.2 两独立样本比较的Wilcoxon秩和检验

Wilcoxon秩和检验（rank sum test），也称 Mann-Whitney U 检验，是两独立样本t检验的非参数替代。

使用课本**例8-3**的数据——两组独立样本的指标值：

```r
RD1 <- c(2.78, 3.23, 4.20, 4.87, 5.12, 6.21, 7.18, 8.05, 8.56, 9.60)
RD2 <- c(3.23, 3.50, 4.04, 4.15, 4.28, 4.34, 4.47, 4.64, 4.75, 4.82, 4.95, 5.10)
```

**方法一：两向量形式**

```r
wilcox.test(RD1, RD2, paired = FALSE, correct = FALSE, exact = FALSE)
## 
##  Wilcoxon rank sum test
## 
## data:  RD1 and RD2
## W = 86.5, p-value = 0.08049
## alternative hypothesis: true location shift is not equal to 0
```

**方法二：formula 形式（长数据格式）**

```r
value <- c(RD1, RD2)
group <- rep(c("A", "B"), c(length(RD1), length(RD2)))
wilcox.test(value ~ group)
## 
##  Wilcoxon rank sum test with continuity correction
## 
## data:  value by group
## W = 86.5, p-value = 0.08134
## alternative hypothesis: true location shift is not equal to 0
```

参数说明：`correct = FALSE` 取消连续性校正，`exact = FALSE` 使用正态近似计算p值（样本量较大或有相同秩次时需要）。

### 6.3 完全随机设计多个样本比较的 Kruskal-Wallis H 检验

Kruskal-Wallis H检验是单因素方差分析的非参数替代，用于比较三组及以上独立样本的分布差异。

使用课本**例8-5**的数据——三种药物处理后的死亡率：

```r
death_rate <- c(32.5, 35.5, 40.5, 46, 49, 16, 20.5, 22.5, 29, 36,
                6.5, 9.0, 12.5, 18, 24)
drug <- rep(c("Drug_A", "Drug_B", "Drug_C"), each = 5)
mydata <- data.frame(death_rate, drug)

str(mydata)
## 'data.frame':    15 obs. of  2 variables:
##  $ death_rate: num  32.5 35.5 40.5 46 49 16 20.5 22.5 29 36 ...
##  $ drug      : chr  "Drug_A" "Drug_A" "Drug_A" "Drug_A" ...
```

```r
kruskal.test(death_rate ~ drug, data = mydata)
## 
##  Kruskal-Wallis rank sum test
## 
## data:  death_rate by drug
## Kruskal-Wallis chi-squared = 9.74, df = 2, p-value = 0.007673
```

p = 0.008 < 0.05，三组间差异有统计学意义。继续使用例8-6数据演示另一个案例：

```r
data8_6 <- data.frame(
  days = c(2, 2, 2, 3, 4, 4, 4, 5, 7, 7,
           5, 5, 6, 6, 6, 7, 8, 10, 12,
           3, 5, 6, 6, 6, 7, 7, 9, 10, 11, 11),
  type = c(rep("9D", 10), rep("11C", 9), rep("DSC", 11))
)

kruskal.test(days ~ type, data = data8_6)
## 
##  Kruskal-Wallis rank sum test
## 
## data:  days by type
## Kruskal-Wallis chi-squared = 9.9405, df = 2, p-value = 0.006941
```

#### Kruskal-Wallis H检验后的多重比较（Nemenyi检验）

Kruskal-Wallis检验显著后，需要事后两两比较。使用 `PMCMRplus` 包的 Nemenyi检验：

```r
library(PMCMRplus)

data8_6$type <- factor(data8_6$type)

res <- kwAllPairsNemenyiTest(days ~ type, data = data8_6)
## Warning in kwAllPairsNemenyiTest.default(...): Ties are present, p-values are not corrected.

summary(res)
##                q value Pr(>|q|)  
## 9D - 11C == 0    3.628 0.027794 *
## DSC - 11C == 0   0.177 0.991411  
## DSC - 9D == 0    3.998 0.013057 *
```

结果显示 DSC 与 9D、9D 与 11C 之间差异有统计学意义。

### 6.4 随机区组设计多个样本比较的Friedman M检验

Friedman M检验是重复测量方差分析（或随机化区组设计的方差分析）的非参数替代。

使用课本**例8-9**的数据——8名受试者在4种频率下的反应率：

```r
# 读入SPSS数据
df <- foreign::read.spss("datasets/例08-09.sav", to.data.frame = TRUE)

str(df)
## 'data.frame':    8 obs. of  4 variables:
##  $ a: num  8.4 11.6 9.4 9.8 8.3 8.6 8.9 7.8
##  $ b: num  9.6 12.7 9.1 8.7 8 9.8 9 8.2
##  $ c: num  9.8 11.8 10.4 9.9 8.6 9.6 10.6 8.5
##  $ d: num  11.7 12 9.8 12 8.6 10.6 11.4 10.8
```

```r
M <- as.matrix(df)  # friedman.test() 需要矩阵格式

friedman.test(M)
## 
##  Friedman rank sum test
## 
## data:  M
## Friedman chi-squared = 15.152, df = 3, p-value = 0.001691
```

p = 0.002 < 0.01，四种频率下的反应率差异有统计学意义。

#### Friedman M检验后的两两比较（Quade检验）

```r
# 构建矩阵（课本例8-9数据）
df <- matrix(
  c(8.4, 11.6, 9.4, 9.8, 8.3, 8.6, 8.9, 7.8,
    9.6, 12.7, 9.1, 8.7, 8, 9.8, 9, 8.2,
    9.8, 11.8, 10.4, 9.9, 8.6, 9.6, 10.6, 8.5,
    11.7, 12, 9.8, 12, 8.6, 10.6, 11.4, 10.8),
  byrow = FALSE, nrow = 8,
  dimnames = list(1:8, LETTERS[1:4])
)

library(PMCMRplus)
quadeAllPairsTest(df, dist = "Normal")
##   A       B       C     
## B 0.2200  -       -     
## C 0.0017  0.0644  -     
## D 1.7e-07 7.7e-05 0.0860

# 获取详细结果表格
res <- quadeAllPairsTest(df, dist = "Normal")
toTidy(res)
##   group1 group2 statistic      p.value alternative
## 1      B      A  1.226488 2.200150e-01   two.sided
## 2      C      A  3.526154 1.686568e-03   two.sided
## 3      C      B  2.299666 6.440153e-02   two.sided
## 4      D      A  5.549859 1.715396e-07   two.sided
## 5      D      B  4.323371 7.683144e-05   two.sided
## 6      D      C  2.023706 8.600089e-02   two.sided
##                                           method distribution p.adjust.method
## 1 Quade's testwith standard-normal approximation            z            holm
## ...
```

结果采用 Holm 校正，D 与 A、D 与 B、C 与 A 之间的差异有统计学意义（p < 0.05）。

## 结果解读指南

**wilcox.test()（配对）输出解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|---------|
| `V = ` | Wilcoxon符号秩检验统计量 | 基于正秩和计算，样本量较小时查Wilcoxon符号秩检验界值表 |
| `p-value` | p值 | p < 0.05 表示配对差值的中位数（或分布位置）与0的差异有统计学意义 |

**wilcox.test()（独立样本）输出解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|---------|
| `W = ` | Wilcoxon秩和检验统计量（Mann-Whitney U统计量） | 反映两组秩次的差异程度 |
| `p-value` | p值 | p < 0.05 表示两组分布位置差异有统计学意义 |
| 连续性校正（continuity correction） | 对相同秩次（ties）的调整 | 存在相同观察值时自动开启 |

**kruskal.test() 输出解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|---------|
| `Kruskal-Wallis chi-squared` | H统计量，近似服从χ²分布 | 值越大，组间差异越大 |
| `df` | 自由度 = 组数 - 1 | 用于χ²检验 |
| `p-value` | p值 | p < 0.05 表示至少有两组之间的分布差异有统计学意义 |

**friedman.test() 输出解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|---------|
| `Friedman chi-squared` | M统计量 | 值越大，处理间差异越大 |
| `df` | 自由度 = 处理数 - 1 | |
| `p-value` | p值 | p < 0.05 表示至少有两种处理间差异有统计学意义 |

**关于效应量：**

非参数检验不直接提供均值差或Cohen's d等效应量指标。常用的非参数效应量包括：
- 秩双列相关系数（rank-biserial correlation）：用于两独立样本
- Kendall's W（协调系数）：用于Friedman检验，值在0~1之间
- η²_H（基于H统计量的效应量）：`η²_H = H / (N - 1)`

## 常见问题与注意事项

**Q1：非参数检验和参数检验的根本区别是什么？**

参数检验（t检验、方差分析）假设数据来自特定分布（如正态分布），检验均数差异。非参数检验不假设总体分布，基于数据的**秩次**（排序位置）进行推断，因此适用于等级资料、偏态分布资料和有离群值的数据。代价是：当数据满足参数检验条件时，非参数检验的检验效能（power）低于参数检验。

**Q2：Wilcoxon秩和检验 和 Mann-Whitney U检验 是同一个检验吗？**

在R中，`wilcox.test(..., paired = FALSE)` 同时实现了Wilcoxon秩和检验和Mann-Whitney U检验，两者在数学上是等价的（W值和U值可以互相转换）。但严格来说，统计学历史上两者提出时略有差异——Wilcoxon秩和检验基于两组混合秩次，Mann-Whitney U检验基于两两比较的U统计量。

**Q3：存在较多相同秩次（ties）时怎么办？**

相同观察值会被赋予平均秩次。`wilcox.test()` 和 `kruskal.test()` 会自动进行ties校正。当ties较多时，`exact = FALSE`（正态近似）通常比精确检验更合适；`correct = TRUE`（连续性校正，默认）会使p值略微保守。

**Q4：多组非参数检验显著后，为什么不能直接用 pairwise.wilcox.test 做两两比较？**

`pairwise.wilcox.test()` 虽然能做两两比较并校正p值，但它每次只使用两组数据——这会导致秩次在每次比较中重新计算，改变了检验的基础。Nemenyi检验、Dunn检验等方法基于所有组的整体秩次进行两两比较，更合适。建议Kruskal-Wallis后用 `kwAllPairsNemenyiTest()`（PMCMRplus包），Friedman后用 `quadeAllPairsTest()` 或 `frdAllPairsNemenyiTest()`。

**Q5：为什么Wilcoxon符号秩检验的p值有时和SPSS不一致？**

可能原因：
- R默认使用连续性校正（`correct = TRUE`），SPSS在某些版本不使用
- 存在差值为0的情况处理方式不同（R扣除此类观测并警告，SPSS可能直接排除）
- 精确p值的计算方法不同。如果样本量允许（n < 50且无ties），可以设置 `exact = TRUE` 获得精确p值

**Q6：什么时候用 Nemenyi检验，什么时候用 Quade检验？**

- Nemenyi检验和Dunn检验：用于**完全随机设计的独立样本**（Kruskal-Wallis后的多重比较）
- Quade检验：用于**随机区组设计的相关样本**（Friedman后的多重比较），考虑了区组效应
- 原则上选择与前面主检验设计一致的多重比较方法

**Q7：单侧检验如何操作？**

将 `alternative` 参数设为 `"less"` 或 `"greater"`（默认为 `"two.sided"`）。医学研究绝大多数情况使用双侧检验。

**常见错误提醒：**

1. 配对设计数据错误使用了两独立样本的秩和检验——配对数据两组不独立，必须用 `paired = TRUE`
2. 多组比较后未做事后两两比较就直接报告各组差异——整体检验显著仅表示"至少两组不同"，不能指明具体哪两组
3. 在多组比较中直接做多次Wilcoxon检验代替Kruskal-Wallis + 事后检验——多次比较会显著增加I类错误
4. 忽略了R 4.4.0+中`wilcox.test()` formula方法不支持`paired`参数的更新——请改用两向量形式
5. 误以为非参数检验没有前提条件——虽然不要求正态性，但独立样本秩和检验要求两组分布形状相似
