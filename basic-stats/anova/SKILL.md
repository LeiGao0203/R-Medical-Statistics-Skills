---
name: medical-stat-anova
description: "R语言医学统计：方差分析。涵盖单因素方差分析、随机区组设计方差分析、两两比较（LSD/TukeyHSD/Dunnett/SNK）、Welch方差分析、Bartlett/Levene方差齐性检验。TRIGGER when user mentions 方差分析、ANOVA、多组比较、组间差异、单因素、one-way ANOVA，or asks about comparing means across 3+ groups with continuous data. SKIP for 两组比较（用t检验）、分类数据（用卡方检验）、重复测量（用重复测量方差分析）、需要控制协变量（用协方差分析）。"
---

# 多样本均数比较的方差分析 (ANOVA: Analysis of Variance)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用方差分析：**

- 比较3组及以上独立样本的连续型变量均值差异（完全随机设计的单因素方差分析，如：4种降脂药对低密度脂蛋白的影响比较）
- 配对区组设计中多组均数的比较（随机区组设计方差分析，如：不同药物在同一批动物上的疗效比较，设区组控制个体差异）
- 拉丁方设计中三个因素（处理、行区组、列区组）的方差分析（如：不同药物在不同皮肤部位、不同受试者上的疗效比较）
- 两阶段交叉设计中处理因素、阶段、个体的综合分析
- 结局变量为连续型数值变量，分组变量为多分类变量（≥3个水平）

**何时不用方差分析：**

| 情况 | 替代方法 |
|------|----------|
| 1~2组均数比较 | t检验 |
| 分类数据（率/构成比）的比较 | 卡方检验 / Fisher精确检验 |
| 数据严重偏离正态或方差不齐且无法处理 | Kruskal-Wallis秩和检验（非参数检验） |
| 需要控制协变量 | 协方差分析（ANCOVA） |
| 重复测量设计（同一受试者多次测量） | 重复测量方差分析 |
| 多个因变量同时分析 | 多变量方差分析（MANOVA） |

## 前置条件

**R包依赖：**

```r
# 基础包（无需安装）
# stats —— aov(), summary(), TukeyHSD(), bartlett.test(), oneway.test(), pairwise.t.test()

# 需要安装的包
install.packages("car")         # leveneTest() 方差齐性检验
install.packages("gplots")      # plotmeans() 组均值可视化
install.packages("PMCMRplus")   # lsdTest(), dunnettTest(), snkTest() 多重比较
```

**统计假设（方差分析的前提条件）：**

1. **独立性**（independence）：各观测值相互独立，由研究设计保证
2. **正态性**（normality）：各组的因变量应来自正态分布总体。可使用Shapiro-Wilk检验各组的正态性；方差分析对轻中度偏离正态有一定的稳健性
3. **方差齐性**（homogeneity of variance）：各组总体方差相等。使用Bartlett检验或Levene检验判断；Bartlett检验对正态性敏感，Levene检验更稳健
4. **数据类型**：因变量为连续型数值变量，自变量为多分类因子型变量

**数据格式要求：**

- 长格式（long format）：一列数值型因变量 + 一列多分类分组变量（factor型）
- 随机区组设计：额外包含一列区组变量（factor型）
- 拉丁方设计：包含处理因素、行区组、列区组各一列
- 两阶段交叉设计：包含受试者ID、阶段、处理类型、测量值各一列

## 方法选择决策树

```
你的数据情况 →
├── 研究设计是什么？
│   ├── 完全随机设计（各组独立，无配对关系）
│   │   ├── 方差齐性（Bartlett/Levene p > 0.05）
│   │   │   └── 使用 标准单因素方差分析 → aov(y ~ group)
│   │   └── 方差不齐（Bartlett/Levene p ≤ 0.05）
│   │       └── 使用 Welch方差分析 → oneway.test(y ~ group, var.equal = FALSE)
│   │
│   └── 随机区组设计（同一批次/配对个体在不同处理下测量）
│       └── 使用 随机区组设计方差分析 → aov(y ~ block + group)
│
├── 是否需要两两比较？（ANOVA显著后）
│   ├── 所有组之间两两比较 → TukeyHSD(fit) 或 LSD (宽松)/Bonferroni (保守)
│   ├── 仅与一个对照组比较 → dunnettTest(fit)
│   └── 均数排序/同质子集 → snkTest(fit)（SNK-q检验）
│
└── 特殊设计？
    ├── 拉丁方设计 → aov(y ~ treatment + row_block + col_block)
    └── 两阶段交叉设计 → aov(y ~ phase + treatment + subject_id)
```

**多重比较方法选择指南：**

| 方法 | R函数 | 适用场景 | 特点 |
|------|-------|----------|------|
| LSD-t | `PMCMRplus::lsdTest()` | 探索性分析 | 不校正多重比较，最宽松，易出假阳性 |
| Tukey HSD | `TukeyHSD()` | 所有组之间两两比较（最常用）| 控制Family-wise error rate，较均衡 |
| Bonferroni | `pairwise.t.test(..., p.adjust.method = "bonferroni")` | 严格控制I类错误 | 最保守，组数多时效能很低 |
| Dunnett-t | `PMCMRplus::dunnettTest()` | 多个实验组 vs 一个对照组 | 仅作g-1次比较，效能高于Tukey |
| SNK-q | `PMCMRplus::snkTest()` | 探索均数同质子集 | 控制每个step的错误率 |

## 标准工作流

### 步骤1：数据准备与探索

读入数据或手动构造数据框。检查数据结构（`str()`），计算各组描述性统计量。可视化探查数据分布。

```r
# 数据格式：一列分组，一列测量值
str(data1)
## 'data.frame':  120 obs. of  2 variables:
##  $ trt   : chr  "group1" "group1" "group1" ...
##  $ weight: num  3.53 4.59 4.34 2.66 ...

# 分组描述统计
tapply(data1$weight, data1$trt, mean)
tapply(data1$weight, data1$trt, sd)
tapply(data1$weight, data1$trt, length)

# 箱线图探索
boxplot(weight ~ trt, data = data1, xlab = "组别", ylab = "测量值")
```

### 步骤2：前提条件检验

检验方差齐性（ANOVA最重要的前提条件）。

```r
# Bartlett检验（数据需近似正态）
bartlett.test(weight ~ trt, data = data1)

# Levene检验（对非正态更稳健，推荐）
library(car)
leveneTest(weight ~ trt, data = data1)
```

### 步骤3：执行统计分析

根据方差齐性检验结果选择合适方法，执行ANOVA及多重比较。

```r
# 方差齐 → 标准ANOVA
fit <- aov(weight ~ trt, data = data1)
summary(fit)

# 方差不齐 → Welch ANOVA
oneway.test(weight ~ trt, data = data1, var.equal = FALSE)

# ANOVA显著后 → 多重比较
TukeyHSD(fit)  # 或 lsdTest(fit) / dunnettTest(fit)
```

### 步骤4：结果解读

查看F值、组间/组内自由度、p值。p < 0.05表示至少有一组与其他组不同；多重比较结果明确具体哪些组之间有差异。

### 步骤5：结果报告

论文中报告格式示例："四种药物的低密度脂蛋白降低值分别为 group1 3.43±0.58、group2 2.72±0.65、group3 2.70±0.72、group4 1.97±0.52。方差分析显示四组间差异有统计学意义（F = 24.88, df1 = 3, df2 = 116, p < 0.001）。Tukey HSD两两比较显示，group1与其他三组均有显著差异，group2与group4、group3与group4之间差异有统计学意义，group2与group3之间差异无统计学意义。"

## 代码示例

### 完全随机设计单因素方差分析（课本例4-2）

```r
trt <- c(rep("group1", 30), rep("group2", 30), rep("group3", 30), rep("group4", 30))

weight <- c(3.53, 4.59, 4.34, 2.66, 3.59, 3.13, 3.30, 4.04, 3.53, 3.56, 3.85, 4.07, 1.37,
            3.93, 2.33, 2.98, 4.00, 3.55, 2.64, 2.56, 3.50, 3.25, 2.96, 4.30, 3.52, 3.93,
            4.19, 2.96, 4.16, 2.59, 2.42, 3.36, 4.32, 2.34, 2.68, 2.95, 2.36, 2.56, 2.52,
            2.27, 2.98, 3.72, 2.65, 2.22, 2.90, 1.98, 2.63, 2.86, 2.93, 2.17, 2.72, 1.56,
            3.11, 1.81, 1.77, 2.80, 3.57, 2.97, 4.02, 2.31, 2.86, 2.28, 2.39, 2.28, 2.48,
            2.28, 3.48, 2.42, 2.41, 2.66, 3.29, 2.70, 2.66, 3.68, 2.65, 2.66, 2.32, 2.61,
            3.64, 2.58, 3.65, 3.21, 2.23, 2.32, 2.68, 3.04, 2.81, 3.02, 1.97, 1.68, 0.89,
            1.06, 1.08, 1.27, 1.63, 1.89, 1.31, 2.51, 1.88, 1.41, 3.19, 1.92, 0.94, 2.11,
            2.81, 1.98, 1.74, 2.16, 3.37, 2.97, 1.69, 1.19, 2.17, 2.28, 1.72, 2.47, 1.02,
            2.52, 2.10, 3.71)

data1 <- data.frame(trt, weight)
head(data1)
##      trt weight
## 1 group1   3.53
## 2 group1   4.59
## 3 group1   4.34
## 4 group1   2.66
## 5 group1   3.59
## 6 group1   3.13

boxplot(weight ~ trt, data = data1)

fit <- aov(weight ~ trt, data = data1)
summary(fit)
##              Df Sum Sq Mean Sq F value   Pr(>F)
## trt           3  32.16  10.719   24.88 1.67e-12 ***
## Residuals   116  49.97   0.431
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### 随机区组设计方差分析（课本例4-4）

```r
weight <- c(0.82, 0.65, 0.51, 0.73, 0.54, 0.23, 0.43, 0.34, 0.28, 0.41, 0.21,
            0.31, 0.68, 0.43, 0.24)
block <- c(rep(c("1", "2", "3", "4", "5"), each = 3))
group <- c(rep(c("A", "B", "C"), 5))

data4_4 <- data.frame(weight, block, group)
head(data4_4)
##   weight block group
## 1   0.82     1     A
## 2   0.65     1     B
## 3   0.51     1     C
## 4   0.73     2     A
## 5   0.54     2     B
## 6   0.23     2     C

fit <- aov(weight ~ block + group, data = data4_4)
summary(fit)
##             Df Sum Sq Mean Sq F value  Pr(>F)
## block        4 0.2284 0.05709   5.978 0.01579 *
## group        2 0.2280 0.11400  11.937 0.00397 **
## Residuals    8 0.0764 0.00955
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

> **注意**：随机区组设计中 formula 的书写顺序为 `block + group`，block在前可先从总变异中扣除区组间变异。

### Welch方差分析（方差不齐时使用）

```r
# Welch方差分析，不要求方差齐性
oneway.test(weight ~ trt, data = data1, var.equal = FALSE)
## 
##  One-way analysis of means (not assuming equal variances)
## 
## data:  weight and trt
## F = 24.003, num df = 3.000, denom df = 60.096, p-value = 2.047e-10
```

### 多重比较

**LSD-t检验（课本例4-7）：**

```r
library(PMCMRplus)

fit <- aov(weight ~ trt, data = data1)
res <- lsdTest(fit)
# 也可以：lsdTest(weight ~ trt, data = data1)

summary(res)
##                      t value   Pr(>|t|)
## group2 - group1 == 0  -4.219 4.8872e-05 ***
## group3 - group1 == 0  -4.322 3.2889e-05 ***
## group4 - group1 == 0  -8.639 3.5772e-14 ***
## group3 - group2 == 0  -0.102    0.91871
## group4 - group2 == 0  -4.420 2.2345e-05 ***
## group4 - group3 == 0  -4.318 3.3397e-05 ***

# 可视化（相同字母表示差异不显著）
plot(res)
```

**TukeyHSD检验：**

```r
TukeyHSD(fit)
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = weight ~ trt, data = data1)
## 
## $trt
##                      diff        lwr        upr     p adj
## group2-group1 -0.71500000 -1.1567253 -0.2732747 0.0002825
## group3-group1 -0.73233333 -1.1740587 -0.2906080 0.0001909
## group4-group1 -1.46400000 -1.9057253 -1.0222747 0.0000000
## group3-group2 -0.01733333 -0.4590587  0.4243920 0.9996147
## group4-group2 -0.74900000 -1.1907253 -0.3072747 0.0001302
## group4-group3 -0.73166667 -1.1733920 -0.2899413 0.0001938

# 可视化（置信区间包含0表示差异不显著）
par(las = 2)
par(mar = c(5, 8, 4, 2))
plot(TukeyHSD(fit))
```

**Dunnett-t检验（实验组与对照组比较）：**

```r
library(PMCMRplus)

res <- dunnettTest(fit)
# 或者 dunnettTest(weight ~ trt, data = data1)

summary(res)
##                      t value   Pr(>|t|)
## group2 - group1 == 0  -4.219 0.00012148 ***
## group3 - group1 == 0  -4.322 0.00010083 ***
## group4 - group1 == 0  -8.639 1.4655e-14 ***

plot(res)
```

> Dunnett-t检验仅做实验组与对照组的比较（g-1次），不做实验组之间的两两比较。

**Bonferroni校正的两两比较：**

```r
# pairwise.t.test 是基础R中最灵活的多重比较函数
# 支持多种p值校正方法：bonferroni, holm, BH, BY, fdr 等
pairwise.t.test(data1$weight, data1$trt, p.adjust.method = "bonferroni")
## 
##  Pairwise comparisons using t tests with pooled SD 
## 
## data:  data1$weight and data1$trt 
## 
##        group1  group2  group3 
## group2 8.8e-05 -       -      
## group3 5.9e-05 1.0000  -      
## group4 6.4e-13 4.0e-05 6.0e-05
## 
## P value adjustment method: bonferroni
```

**SNK-q检验（课本例4-9）：**

```r
library(PMCMRplus)

data4_4$group <- factor(data4_4$group)

fit <- aov(weight ~ group, data = data4_4)
res <- snkTest(fit)

summary(res)
##            q value Pr(>|q|)
## B - A == 0  -2.526 0.099390 .
## C - A == 0  -4.209 0.028913 *
## C - B == 0  -1.684 0.256834

plot(res)
```

### 方差齐性检验

**Bartlett检验（课本例4-10）：**

```r
bartlett.test(weight ~ trt, data = data1)
## 
##  Bartlett test of homogeneity of variances
## 
## data:  weight by trt
## Bartlett's K-squared = 5.2192, df = 3, p-value = 0.1564
# p > 0.05 表示满足方差齐性，可进行标准ANOVA
```

**Levene检验：**

```r
library(car)

leveneTest(weight ~ trt, data = data1)
## Levene's Test for Homogeneity of Variance (center = median)
##        Df F value Pr(>F)
## group   3   1.493 0.2201
##       116
# p > 0.05 表示满足方差齐性
```

> Levene检验以中位数为中心（默认），对非正态数据更稳健，是SPSS中默认的方差齐性检验方法。

### 拉丁方设计方差分析（课本例4-5）

```r
psize <- c(87, 75, 81, 75, 84, 66, 73, 81, 87, 85, 64, 79, 73, 73, 74, 78, 73, 77,
           77, 68, 69, 74, 76, 73, 64, 64, 72, 76, 70, 81, 75, 77, 82, 61, 82, 61)
drug <- c("C", "B", "E", "D", "A", "F", "B", "A", "D", "C", "F", "E",
          "F", "E", "B", "A", "D", "C", "A", "F", "C", "B", "E", "D",
          "D", "C", "F", "E", "B", "A", "E", "D", "A", "F", "C", "B")
col_block <- c(rep(1:6, 6))
row_block <- c(rep(1:6, each = 6))

mydata <- data.frame(psize, drug, col_block, row_block)
mydata$col_block <- factor(mydata$col_block)
mydata$row_block <- factor(mydata$row_block)

fit <- aov(psize ~ drug + row_block + col_block, data = mydata)
summary(fit)
##             Df Sum Sq Mean Sq F value Pr(>F)
## drug         5  667.1  133.43   3.906 0.0124 *
## row_block    5  250.5   50.09   1.466 0.2447
## col_block    5   85.5   17.09   0.500 0.7723
## Residuals   20  683.2   34.16
```

### 两阶段交叉设计方差分析（课本例4-6）

```r
contain <- c(760, 770, 860, 855, 568, 602, 780, 800, 960, 958, 940, 952,
             635, 650, 440, 450, 528, 530, 800, 803)
phase <- rep(c("phase_1", "phase_2"), 10)
type <- c("A", "B", "B", "A", "A", "B", "A", "B", "B", "A",
          "B", "A", "A", "B", "B", "A", "A", "B", "B", "A")
testid <- rep(1:10, each = 2)

mydata <- data.frame(testid, phase, type, contain)
mydata$testid <- factor(mydata$testid)

fit <- aov(contain ~ phase + type + testid, data = mydata)
summary(fit)
##             Df Sum Sq Mean Sq  F value   Pr(>F)
## phase        1    490     490    9.925   0.0136 *
## type         1    198     198    4.019   0.0799 .
## testid       9 551111   61235 1240.195 1.32e-11 ***
## Residuals    8    395      49
```

### 组均值可视化（带95%置信区间）

```r
library(gplots)

plotmeans(weight ~ trt, data = data1,
          xlab = "组别", ylab = "测量值",
          main = "各组件本均值及95%置信区间")
```

## 结果解读指南

**summary(aov()) 输出各组成部分的解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `Df` | 自由度 | 组间df = k-1（k为组数）；组内df = N-k；区组df = b-1 |
| `Sum Sq` | 离均差平方和（SS） | 组间SS：组均值与总均值差异；组内SS：个体值与组均值差异（误差） |
| `Mean Sq` | 均方（MS = SS/df） | 组间MS除以组内MS得到F值 |
| `F value` | F统计量 | F = MS_between / MS_within，越大组间差异越大 |
| `Pr(>F)` | p值 | p < 0.05 → 至少有一组与其他组不同；p ≥ 0.05 → 各组间差异无统计学意义 |

**TukeyHSD 输出解读：**

| 输出项 | 含义 |
|--------|------|
| `diff` | 两组均数之差 |
| `lwr` / `upr` | 差值的95%置信区间下/上限 |
| `p adj` | 校正后的p值（控制family-wise error rate） |

p adj < 0.05 表示该对比较差异有统计学意义。可视化图中置信区间包含0的pair表示差异不显著。

**Levene检验/Bartlett检验解读：**

- H₀：各组方差相等（方差齐性）
- p > 0.05 → 不拒绝H₀，可认为方差齐性满足，使用标准ANOVA
- p ≤ 0.05 → 拒绝H₀，方差不齐，使用Welch方差分析（`oneway.test()`）或考虑Kruskal-Wallis秩和检验

**多重比较方法选择依据：**

- TukeyHSD：最常用，适合所有组两两比较，控制总体I类错误
- LSD：不校正，易出假阳性，仅适合探索性分析或3组时使用
- Dunnett：专门用于实验组与对照组的比较（不关心实验组之间比较）
- Bonferroni：最保守，适用于严格控制的验证性研究
- SNK：给出同质子集分组，适合探索均数排序模式

## 常见问题与注意事项

**Q1：方差分析显著（p < 0.05），是否一定要做多重比较？**

是的。方差分析只能告诉你"至少有一组与其他组不同"，不能告诉你具体哪两组有差异。必须通过多重比较（TukeyHSD、LSD等）才能明确具体差异来源。

**Q2：随机区组设计中 formula 顺序有讲究吗？**

有。`aov(y ~ block + group)` 中 block 放在前面可以优先从总变异中扣除区组间变异，从而更准确地评估处理因素（group）的效应。如果顺序反过来（`aov(y ~ group + block)`），在非平衡设计中SS分配会不同。

**Q3：Bartlett检验和Levene检验应该选哪个？**

- Bartlett检验：要求数据服从正态分布，对正态性偏离敏感；数据确实正态时效能较高
- Levene检验：对非正态数据更稳健，是SPSS的默认方差齐性检验方法
- 建议：常规使用Levene检验（`car::leveneTest()`），仅在确认数据正态时考虑Bartlett检验

**Q4：方差不齐怎么办？**

三个层次的处理方案：
1. 使用Welch方差分析（`oneway.test(y ~ group, var.equal = FALSE)`），这是最直接的替代方法
2. 尝试数据转换（对数转换、平方根转换等）使之满足方差齐性
3. 改用非参数检验：Kruskal-Wallis秩和检验（`kruskal.test(y ~ group)`）

**Q5：`pairwise.t.test()` 和 `TukeyHSD()` 的区别？**

- `pairwise.t.test()`：调用的是两两t检验，通过 `p.adjust.method` 参数支持多种p值校正方法（bonferroni、holm、BH、fdr等），更灵活
- `TukeyHSD()`：基于学生化范围分布（studentized range distribution），专为ANOVA设计，同时给出差值置信区间
- 对于标准两两比较，推荐使用 `TukeyHSD()`，因为它天然适配ANOVA框架且结果更完整

**Q6：SPSS和R的ANOVA结果有差异吗？**

基本一致。差别主要在于：
- SPSS使用Type III SS（适用于非平衡设计），R的 `aov()` 默认使用Type I SS（顺序型）
- 对于平衡设计的随机区组ANOVA，两者结果完全相同
- R中如需Type III SS，可使用 `car::Anova(fit, type = "III")`
- SPSS的LSD-t与R的 `lsdTest()` 的t值可能会有微小差异（约±0.01），不影响结论

**Q7：LSD、Tukey、Bonferroni的多重比较P值为什么不同？**

因为三者对多重比较的校正方式不同：
- LSD：不校正，直接用α=0.05做每次比较 → 最松
- Tukey HSD：基于q分布，控制所有比较的整体错误率（family-wise error rate）→ 适中
- Bonferroni：α' = α/m（m为比较次数）→ 最严，当组数较多时过于保守
- 在3组时三者结论通常一致；组数越多，Bonferroni越保守，LSD假阳性越多

**常见错误提醒：**

1. 多组比较错误地用多次两两t检验代替ANOVA——会放大I类错误（多重比较问题），应先用ANOVA判断总体差异，再做事后比较
2. 随机区组设计按完全随机设计分析——忽略了区组因素，损失统计效能，且可能得出错误结论
3. 组变量未设为factor型——R中字符型和因子型在 `aov()` 中行为不同，务必 `data$group <- factor(data$group)`
4. 忽略方差齐性检验直接进行标准ANOVA——当方差不齐且各组样本量不等时影响较大
5. 将p > 0.05解释为"各组相等"——ANOVA不显著仅说明未发现组间差异有统计学意义，不等同于各组均数相等
