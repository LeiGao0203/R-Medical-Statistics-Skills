---
name: medical-stat-ttest
description: "R语言医学统计：t检验。涵盖单样本t检验、配对t检验、两独立样本t检验、正态性检验和方差齐性检验。TRIGGER when user mentions t检验、均值比较、两组差异分析、配对设计，or asks about comparing means between 1-2 groups with continuous data. SKIP for 多组比较（用ANOVA）、分类数据（用卡方检验）、非正态数据（用非参数检验）。"
---

# t检验 (t-test)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用 t 检验：**

- 比较一组样本的均值与已知总体均值（单样本t检验，如：某地区成年男性血红蛋白含量与全国标准140 g/L的比较）
- 比较同一组受试者在两种条件下的测量值（配对t检验，如：治疗前后、左右侧配对、同体不同部位）
- 比较两个独立组的连续型变量均值差异（两独立样本t检验，如：两种降糖药物的疗效比较）
- 结局变量为连续型数值变量（如血糖、血压、血红蛋白），且只有1~2个组

**何时不用 t 检验：**

| 情况 | 替代方法 |
|------|----------|
| 3组及以上均数比较 | 单因素方差分析（ANOVA）或 Kruskal-Wallis 检验 |
| 分类数据（率/构成比） | 卡方检验 / Fisher 精确检验 |
| 数据严重偏离正态且无法转换 | Wilcoxon 符号秩检验（配对）或 Mann-Whitney U 检验（独立样本） |
| 需要控制协变量 | 协方差分析（ANCOVA）或多重线性回归 |
| 重复测量设计（≥3个时间点） | 重复测量方差分析 |

## 前置条件

**R 包依赖：**

```r
# 基础包（无需安装）
# stats    —— t.test(), shapiro.test(), var.test()
# foreign  —— read.spss() 读取SPSS数据

# 可选安装
install.packages("moments")  # 偏度/峰度计算
install.packages("tidyverse")# 数据整理（pivot_longer等）
```

**统计假设（t 检验的前提条件）：**

1. **正态性**（normality）：数据应来自正态分布总体。使用 Shapiro-Wilk 检验判断；小样本时正态性尤为重要
2. **方差齐性**（homogeneity of variance）：两独立样本t检验要求两组方差相等。使用 F 检验（`var.test()`）或 Levene 检验判断
3. **独立性**（independence）：各观测值之间相互独立（由研究设计保证）
4. **数据类型**：结局变量为连续型数值变量

**数据格式要求：**

- 单样本：一列数值型变量
- 配对样本：每个受试者占一行，待比较的两列数值变量（宽格式），或多行长格式
- 两独立样本：一列数值型结局变量 + 一列二分类分组变量（factor/character）

## 方法选择决策树

```
你的数据情况 →
├── 1组样本，与已知常模/标准值比较
│   └── 使用 单样本t检验 → t.test(x, mu = 总体均值)
│
├── 2组数据，配对设计（同一受试者、自身前后、左右配对）
│   ├── 差值满足正态性 → 使用 配对t检验 → t.test(x1, x2, paired = TRUE)
│   └── 差值不满足正态性 → 使用 Wilcoxon符号秩检验 → wilcox.test(x1, x2, paired = TRUE)
│
└── 2组数据，独立样本设计（完全随机分组）
    ├── 第一步：正态性检验（shapiro.test()）
    │   ├── 两组均满足正态性 → 继续第二步
    │   └── 任一组不满足 → 使用 Mann-Whitney U检验 或数据转换
    │
    └── 第二步：方差齐性检验（var.test()）
        ├── 方差齐（p > 0.05）→ 使用 Student's t检验 → t.test(x ~ group, var.equal = TRUE)
        └── 方差不齐（p ≤ 0.05）→ 使用 Welch近似t检验 → t.test(x ~ group, var.equal = FALSE)
```

## 标准工作流

### 步骤1：数据准备与探索

读入数据，检查结构（`str()`），计算各组描述性统计量（均值、标准差、样本量）。

```r
library(foreign)
df <- read.spss('datasets/例03-07.sav', to.data.frame = TRUE)
df$group <- c(rep('阿卡波糖', 20), rep('拜糖平', 20))
attributes(df)[3] <- NULL

# 分组描述统计
tapply(df$x, df$group, summary)
tapply(df$x, df$group, sd)
```

### 步骤2：前提条件检验

先检验正态性，再检验方差齐性（仅两独立样本需要）。

```r
# 正态性检验（Shapiro-Wilk）
shapiro.test(df$x[df$group == '阿卡波糖'])
shapiro.test(df$x[df$group == '拜糖平'])

# 方差齐性检验（F检验）
var.test(x ~ group, data = df)
```

### 步骤3：执行统计分析

根据步骤2的结果选择对应的t检验方法。所有t检验都通过 `t.test()` 函数完成，仅参数不同。

### 步骤4：结果解读

查看 t 值、自由度、p 值和95%置信区间。p < 0.05 表示差异有统计学意义，临床意义需结合效应量（均差及置信区间）。

### 步骤5：结果报告

论文中报告格式示例："阿卡波糖组血糖下降值为 2.07 ± 3.06，拜糖平组为 2.63 ± 2.41，两组比较差异无统计学意义（t = -0.64, df = 38, p = 0.525）。"

## 代码示例

### 单样本t检验（课本例3-5）

```r
library(foreign)
df <- read.spss('datasets/例03-05.sav', to.data.frame = T)
head(df)
##   no  hb
## 1  1 112
## 2  2 137
## 3  3 129
## 4  4 126
## 5  5  88
## 6  6  90

st <- t.test(df$hb, mu = 140, alternative = 'two.sided')
st
## 
##  One Sample t-test
## 
## data:  df$hb
## t = -2.1367, df = 35, p-value = 0.03969
## alternative hypothesis: true mean is not equal to 140
## 95 percent confidence interval:
##  122.1238 139.5428
## sample estimates:
## mean of x 
##  130.8333
```

### 配对样本t检验（课本例3-6）

```r
library(foreign)
df <- read.spss('datasets/例03-06.sav', to.data.frame = T)
head(df)
##   no    x1    x2
## 1  1 0.840 0.580
## 2  2 0.591 0.509
## 3  3 0.674 0.500
## 4  4 0.632 0.316
## 5  5 0.687 0.337
## 6  6 0.978 0.517

pt <- t.test(df$x1, df$x2, paired = TRUE, var.equal = TRUE)
pt
## 
##  Paired t-test
## 
## data:  df$x1 and df$x2
## t = 7.926, df = 9, p-value = 2.384e-05
## alternative hypothesis: true mean difference is not equal to 0
## 95 percent confidence interval:
##  0.1946542 0.3501458
## sample estimates:
## mean difference 
##          0.2724
```

**使用 formula 形式（需长数据格式）：**

```r
library(tidyverse)
df.l <- df |> pivot_longer(2:3, names_to = "group", values_to = "x")
t.test(x ~ group, data = df.l, paired = TRUE, var.equal = TRUE)
## 
##  Paired t-test
## 
## data:  x by group
## t = 7.926, df = 9, p-value = 2.384e-05
## ...
```

> **注意**：R 4.4.0 起，`t.test()` 的 formula 方法不再支持 `paired` 参数，建议使用两向量形式（`t.test(x1, x2, paired = TRUE)`）。

### 两独立样本t检验（课本例3-7）

```r
library(foreign)
df <- read.spss('datasets/例03-07.sav', to.data.frame = T)
df$group <- c(rep('阿卡波糖', 20), rep('拜糖平', 20))
attributes(df)[3] <- NULL
head(df)
##   no    x    group
## 1  1 -0.7 阿卡波糖
## 2  2 -5.6 阿卡波糖
## 3  3  2.0 阿卡波糖
## 4  4  2.8 阿卡波糖
## 5  5  0.7 阿卡波糖
## 6  6  3.5 阿卡波糖

# 等方差t检验（Student's t-test）
tt <- t.test(x ~ group, data = df, var.equal = TRUE)
tt
## 
##  Two Sample t-test
## 
## data:  x by group
## t = -0.64187, df = 38, p-value = 0.5248
## alternative hypothesis: true difference in means between group 阿卡波糖 and group 拜糖平 is not equal to 0
## 95 percent confidence interval:
##  -2.326179  1.206179
## sample estimates:
## mean in group 阿卡波糖   mean in group 拜糖平 
##                  2.065                  2.625

# Welch近似t检验（方差不齐时使用）
tw <- t.test(x ~ group, data = df, var.equal = FALSE)
tw
## 
##  Welch Two Sample t-test
## 
## data:  x by group
## t = -0.64187, df = 36.086, p-value = 0.525
## alternative hypothesis: true difference in means between group 阿卡波糖 and group 拜糖平 is not equal to 0
## 95 percent confidence interval:
##  -2.32926  1.20926
## sample estimates:
## mean in group 阿卡波糖   mean in group 拜糖平 
##                  2.065                  2.625
```

### 正态性检验与方差齐性检验

```r
library(moments)

# Shapiro-Wilk正态性检验
shapiro.test(df$mean)
## 
##  Shapiro-Wilk normality test
## 
## data:  df$mean
## W = 0.99409, p-value = 0.9444

# 偏度与峰度
skewness(df$mean)
## [1] 0.1423707
kurtosis(df$mean)
## [1] 3.045566

# D'Agostino偏度检验
agostino.test(df$mean)
##  D'Agostino skewness test
## skew = 0.14237, z = 0.61614, p-value = 0.5378

# Anscombe-Glynn峰度检验
anscombe.test(df$mean)
##  Anscombe-Glynn kurtosis test
## kurt = 3.04557, z = 0.41992, p-value = 0.6745

# F检验（两样本方差齐性检验）
var.test(x ~ group, data = df)
## 
##  F test to compare two variances
## 
## data:  x by group
## F = 1.5984, num df = 19, denom df = 19, p-value = 0.3153
## alternative hypothesis: true ratio of variances is not equal to 1
## 95 percent confidence interval:
##  0.6326505 4.0381795
## sample estimates:
## ratio of variances 
##           1.598361
```

## 结果解读指南

**t.test() 输出各组成部分的解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `t = ` | t 统计量的值 | 绝对值越大，两组差异越大；正负号表示方向 |
| `df = ` | 自由度（degrees of freedom） | 配对：n-1；等方差两样本：n1+n2-2；Welch校正：Satterthwaite公式 |
| `p-value = ` | p 值 | p < 0.05 拒绝 H₀（差异有统计学意义） |
| `95 percent confidence interval` | 均值差/均值的95%置信区间 | 不包含0（单样本不包含 μ₀）等价于 p < 0.05 |
| `sample estimates` | 样本均数估计值 | 描述性数值，用于报告 |

**常见 p 值表述对照：**

- p > 0.05：差异无统计学意义 → 不拒绝 H₀
- p < 0.05：差异有统计学意义 → 拒绝 H₀
- p < 0.01：差异有高度统计学意义
- p < 0.001：差异有极高度统计学意义

**shapiro.test() 解读：**

- H₀：数据来自正态分布总体
- p > 0.05 → 不拒绝 H₀，可认为数据满足正态性
- p ≤ 0.05 → 拒绝 H₀，数据不满足正态性，考虑非参数检验或数据转换
- Shapiro-Wilk 是医学统计中最常用的正态性检验

**var.test() 解读：**

- H₀：两总体方差相等（方差比 = 1）
- p > 0.05 → 方差齐，使用 Student's t 检验（`var.equal = TRUE`）
- p ≤ 0.05 → 方差不齐，使用 Welch 近似 t 检验（`var.equal = FALSE`）

**偏度与峰度判断：**

- 偏度（skewness）：表示分布对称性，0为完全对称，正值右偏，负值左偏
- 峰度（kurtosis）：表示分布陡峭程度，正态分布峰度 ≈ 3，>3 峰更尖，<3 峰更平
- 偏度、峰度同时不显著（对应检验 p > 0.05）可作为正态性的辅助证据

## 常见问题与注意事项

**Q1：配对t检验和两独立样本t检验的本质区别是什么？**

配对t检验检验的是**差值的均值是否为0**，其本质是对配对差值进行单样本t检验（H₀: μ_d = 0）。配对设计能消除个体间变异，提高检验效能。两独立样本t检验则直接比较两组均值差异。

**Q2：大样本（n > 50）时还需要正态性检验吗？**

根据中心极限定理，大样本下样本均值近似正态分布，t检验对正态性偏离有一定稳健性。但严重偏态的数据仍建议使用非参数检验或进行变量转换（如对数转换）。实际应用中，多数医学研究仍会报告正态性检验结果。

**Q3：t.test() 中 var.equal 到底应该设为 TRUE 还是 FALSE？**

如果事前经过 F 检验（`var.test()`）确认方差齐性，设置 `var.equal = TRUE`（经典的 Student's t 检验）。如果没有检验或方差不齐，设置 `var.equal = FALSE`（默认值，Welch 近似 t 检验）。某些统计教材建议默认使用 Welch t 检验，因为它在方差不等时更稳健，在方差相等时效率损失也很小。

**Q4：SPSS 和 R 的 t 检验结果有差异吗？**

基本一致。差别主要在于：
- SPSS 默认同时输出等方差和不等方差的 t 检验结果，并附 Levene 方差齐性检验
- R 的 `t.test()` 一次只输出一种（由 `var.equal` 参数控制）
- SPSS 使用 Levene 检验判断方差齐性，R 的 `var.test()` 是 F 检验。两者的结论通常一致

**Q5：单侧检验（one-sided test）如何实现？**

将 `alternative` 参数设为 `'less'`（H₁: μ < μ₀）或 `'greater'`（H₁: μ > μ₀），默认为 `'two.sided'`（双侧检验）。医学研究通常使用双侧检验，除非有充分先验依据支持单侧。

**Q6：正态性检验 p 值略小于 0.05，还能用 t 检验吗？**

若样本量较大（n > 30）且偏度/峰度不严重，t检验仍有一定稳健性。也可以尝试：
- 非参数检验：Wilcoxon 符号秩检验（`wilcox.test(..., paired = TRUE)`）或 Mann-Whitney U 检验（`wilcox.test(x ~ group)`）
- 数据转换：对数转换 `log(x)`、平方根转换 `sqrt(x)` 等

**常见错误提醒：**

1. 配对设计数据错误使用了两独立样本t检验——配对设计必须使用配对的 t 检验，否则检验效能严重降低
2. 多组比较错误地使用多次两两 t 检验——会增大 I 类错误，应使用 ANOVA + 多重比较校正
3. 忽略方差齐性检验直接使用等方差假设——当样本量不等时影响更大
4. 将 p > 0.05 解释为"两组相等"——不拒绝 H₀ ≠ 证明 H₀ 成立，只能说"未发现差异有统计学意义"
