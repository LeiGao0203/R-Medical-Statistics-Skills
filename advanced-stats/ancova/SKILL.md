---
name: medical-stat-ancova
description: "R语言医学统计：协方差分析。在比较组间差异时控制一个或多个连续型协变量的影响，结合ANOVA和线性回归的思想。TRIGGER when user mentions 协方差分析、ANCOVA、控制协变量、校正基线、控制混杂因素后的组间比较，or asks about adjusting for covariates in group comparisons. SKIP for 简单方差分析（直接比较组间差异不含协变量）、多元线性回归（多自变量预测因变量，不含分组因子）。"
---

# 协方差分析 (ANCOVA — Analysis of Covariance)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用 ANCOVA 的典型场景：**

- 临床试验中需要比较不同治疗组的疗效，但基线指标（如治疗前血糖、血压）在组间不均衡，需要校正基线值
- 比较各组终点指标时，存在已知的连续型混杂因素（如年龄、病程），需要扣除其影响
- 完全随机设计或随机区组设计中，除了分组因素和区组因素外，还有需要控制的定量协变量

**不使用 ANCOVA 的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 仅比较组间均数差异，没有协变量 | 单因素方差分析（`aov(y ~ group)`） |
| 多个自变量共同预测一个连续型因变量，不含分组因素 | 多元线性回归（`lm(y ~ x1 + x2 + x3)`） |
| 因变量为分类变量 | Logistic 回归 |
| 协变量与分组因素存在交互作用（回归线不平行） | 不可直接使用 ANCOVA，需考虑分层分析或纳入交互项 |

**医学研究常见应用：**

- 比较三种降糖药物的降糖效果，校正治疗前的糖化血红蛋白（HbA1c）水平（课本例13-1）
- 比较不同饲料喂养下猪的增重，校正初始体重并控制窝别（区组）效应（课本例13-2）

## 前置条件

**R 包安装：**

```r
install.packages(c("tidyverse", "HH", "rstatix", "car", "foreign", "patchwork"))
```

`HH` 包提供 `ancovaplot()` 函数，需从 CRAN 安装。

**数据格式要求：**

- 因变量 `y`：连续型数值变量
- 协变量 `x`：连续型数值变量（基线指标或混杂因素）
- 分组变量 `group`：因子型（factor），如 "A 组"、"B 组"、"C 组"
- 区组变量 `block`（随机区组设计时）：因子型

数据应为长格式（每个观测一行），若为宽格式（如 x1,y1,x2,y2,x3,y3），需先用 `pivot_longer()` 转换。

**统计假设（必须检验）：**

1. **正态性**：各组的因变量及协变量服从正态分布（或残差正态）
2. **独立性**：各观测相互独立
3. **方差齐性**：各样本总体方差齐（可事后检验残差）
4. **线性关系**：各组内因变量 y 与协变量 x 之间存在线性回归关系
5. **回归线平行（斜率齐性，最重要）**：各组回归斜率相同，即 group:x 交互项不显著

> 若回归线不平行（斜率齐性检验 p < 0.05），则 ANCOVA 结果不可靠，应分层分析或考虑其他方法。

## 方法选择决策树

```
你的数据情况 →
├── 完全随机设计，有一个协变量
│   └── 使用 Type I SS 的 ANCOVA：aov(y ~ x + group) 或 anova_test(type = 1)
│
├── 完全随机设计，想用 tidy 风格输出
│   └── 使用 rstatix::anova_test(y ~ x + group, type = 1)
│
├── 随机区组设计，有协变量 + 区组
│   ├── 先验假设分组最重要 → aov(y ~ x + block + group)
│   ├── 需要 Type II 检验 → car::Anova(fit)
│   └── 需要 Type III 检验 → anova_test(y ~ x + block + group, type = 3)
│
├── 需可视化 ANCOVA 结果
│   ├── 快速出图 → HH::ancovaplot(y ~ x + group)
│   └── 定制美化 → ggplot2 + geom_smooth(method = "lm")
│
└── 回归线不平行（交互项显著）
    └── 分层分析 / 纳入交互项后解释 / 倾向性评分方法
```

## 标准工作流

### 步骤 1：数据准备与探索

若数据为宽格式（每人一行，多列存储），使用 `pivot_longer()` 转换为长格式：

```r
library(tidyverse)

df_long <- df_wide %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "group"),
               names_pattern = "(.)(.)") %>%
  mutate(group = as.factor(group))
```

用 `str()` 或 `glimpse()` 确认数据包含 `x`（协变量）、`y`（因变量）、`group`（分组，factor 型）。

### 步骤 2：前提条件检验

**回归线平行假设（最关键）：**

```r
# 完全随机设计：检验 group:x 交互项
fit_parallel <- aov(y ~ x + group + x:group, data = df)
summary(fit_parallel)
# 若 x:group 的 p > 0.05，满足平行假设
```

**正态性检验：**

```r
# 对 ANCOVA 模型的残差进行检验
fit <- aov(y ~ x + group, data = df)
shapiro.test(residuals(fit))
```

### 步骤 3：执行统计分析

完全随机设计：

```r
# 基础 R 方法 — 协变量必须放在分组变量之前
fit <- aov(y ~ x + group, data = df)
summary(fit)

# rstatix tidy 方法
library(rstatix)
res <- anova_test(y ~ x + group, data = df, type = 1)
get_anova_table(res)
```

随机区组设计：

```r
# 基础 R
fit <- aov(y ~ x + block + group, data = df)
summary(fit)

# Type II SS（不假设因素顺序）
car::Anova(fit)

# rstatix Type III SS
res <- anova_test(y ~ x + block + group, data = df, type = 3)
get_anova_table(res)
```

### 步骤 4：结果解读

- 首先看协变量 `x` 的 `Pr(>F)`：若 p < 0.05，说明协变量对因变量有显著影响，纳入协变量是必要的
- 再看分组变量 `group` 的 `Pr(>F)`：若 p < 0.05，说明在控制协变量后，各组间的修正均数存在统计学差异
- 均方（Mean Sq）= 平方和（Sum Sq）/ 自由度（Df）；F = 组间均方 / 残差均方

以课本例 13-1 为例：F = 58.48, p < 0.001，可认为在扣除初始糖化血红蛋白的影响后，三组患者的降糖效果不同。

### 步骤 5：结果报告（论文写法）

> 采用协方差分析（ANCOVA），以治疗前糖化血红蛋白为协变量，比较三组治疗后的降糖效果。结果显示，在控制基线糖化血红蛋白后，三组间的降糖效果差异具有统计学意义（F(2, 86) = 58.48，p < 0.001）。经校正后三组均数 ± 标准误分别为……

**论文中需报告的内容：**
1. 组别 × 协变量的交互项检验结果（验证回归线平行）
2. ANCOVA 主效应结果（F 值、自由度、p 值）
3. 各组的修正均数（adjusted means）和标准误

## 代码示例

### 示例 1：完全随机设计（课本例 13-1）

研究三种降糖药的降糖效果，以治疗前糖化血红蛋白 `x` 为协变量，治疗后 `y` 为结局。

```r
library(tidyverse)
library(HH)
library(rstatix)

# 数据录入（宽格式）
df13_1 <- data.frame(
  x1 = c(10.8, 11.6, 10.6, 9.0, 11.2, 9.9, 10.6, 10.4, 9.6, 10.5,
         10.6, 9.9, 9.5, 9.7, 10.7, 9.2, 10.5, 11.0, 10.1, 10.7, 8.5,
         10.0, 10.4, 9.7, 9.4, 9.2, 10.5, 11.2, 9.6, 8.0),
  y1 = c(9.4, 9.7, 8.7, 7.2, 10.0, 8.5, 8.3, 8.1, 8.5, 9.1, 9.2, 8.4,
         7.6, 7.9, 8.8, 7.4, 8.6, 9.2, 8.0, 8.5, 7.3, 8.3, 8.6, 8.7,
         7.6, 8.0, 8.8, 9.5, 8.2, 7.2),
  x2 = c(10.4, 9.7, 9.9, 9.8, 11.1, 8.2, 8.8, 10.0, 9.0, 9.4, 8.9,
         10.3, 9.3, 9.2, 10.9, 9.2, 9.2, 10.4, 11.2, 11.1, 11.0, 8.6,
         9.3, 10.3, 10.3, 9.8, 10.5, 10.7, 10.4, 9.4),
  y2 = c(9.2, 9.1, 8.9, 8.6, 9.9, 7.1, 7.8, 7.9, 8.0, 9.0, 7.9, 8.9,
         8.9, 8.1, 10.2, 8.5, 9.0, 8.9, 9.8, 10.1, 8.5, 8.1, 8.6, 8.9,
         9.6, 8.1, 9.9, 9.3, 8.7, 8.7),
  x3 = c(9.8, 11.2, 10.7, 9.6, 10.1, 9.8, 10.1, 10.3, 11.0, 10.5,
         9.2, 10.1, 10.4, 10.0, 8.4, 10.1, 9.3, 10.5, 11.1, 10.5, 9.7,
         9.2, 9.3, 10.4, 10.0, 10.3, 9.9, 9.4, 8.3, 9.2),
  y3 = c(7.6, 7.9, 9.0, 7.8, 8.5, 7.5, 8.3, 8.2, 8.4, 8.1, 7.0, 7.7,
         8.0, 6.6, 6.1, 8.1, 7.8, 8.4, 8.2, 8.0, 7.6, 6.9, 6.7, 8.1,
         7.4, 8.2, 7.6, 7.8, 6.6, 7.2)
)

# 转换为长格式
df13_11 <- df13_1 %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "group"),
               names_pattern = "(.)(.)") %>%
  mutate(group = as.factor(group))

glimpse(df13_11)
## Rows: 90
## Columns: 3
## $ group <fct> 1, 2, 3, 1, 2, 3, ...
## $ x     <dbl> 10.8, 10.4, 9.8, 11.6, 9.7, 11.2, ...
## $ y     <dbl> 9.4, 9.2, 7.6, 9.7, 9.1, 7.9, ...

# 单因素协方差分析 — 协变量必须在分组变量之前
fit <- aov(y ~ x + group, data = df13_11)
summary(fit)
##             Df Sum Sq Mean Sq F value Pr(>F)
## x            1  29.06  29.057  171.20 <2e-16 ***
## group        2  19.85   9.925   58.48 <2e-16 ***
## Residuals   86  14.60   0.170

# rstatix 版本（输出更整洁，含效应量 ges）
res <- anova_test(y ~ x + group, data = df13_11, type = 1)
get_anova_table(res)
## ANOVA Table (type I tests)
##   Effect DFn DFd       F        p p<.05   ges
## 1      x   1  86 171.199 3.64e-22     * 0.666
## 2  group   2  86  58.480 9.22e-17     * 0.576

# 可视化
ancovaplot(y ~ x + group, data = df13_11)
```

### 示例 2：随机区组设计（课本例 13-2）

比较三种饲料的增重效果，以初始体重 `x` 为协变量，12 个窝（block）为区组。

```r
library(foreign)
library(car)
library(rstatix)

# 读取 SPSS 数据
df <- foreign::read.spss("datasets/例13-02.sav",
                         to.data.frame = TRUE,
                         reencode = "utf-8")
df$block <- factor(df$block)

str(df)
## 'data.frame':   36 obs. of  4 variables:
##  $ x    : num  257 272 210 300 262 ...
##  $ y    : num  27 41.7 25 52 14.5 ...
##  $ group: Factor w/ 3 levels "A....","B....",..
##  $ block: Factor w/ 12 levels "1","2","3","4",..

# 基础 R — 注意顺序：协变量 > 区组 > 分组
fit <- aov(y ~ x + block + group, data = df)
summary(fit)
##             Df Sum Sq Mean Sq F value  Pr(>F)
## x            1  69073   69073 651.823 < 2e-16 ***
## block       11   4024     366   3.452 0.00711 **
## group        2    464     232   2.189 0.13692
## Residuals   21   2225     106

# Type II SS（不依赖因素顺序）
car::Anova(fit)
## Anova Table (Type II tests)
## Response: y
##           Sum Sq Df F value    Pr(>F)
## x         6174.2  1 58.2643 1.733e-07 ***
## block     3765.3 11  3.2302   0.01009 *
## group      463.9  2  2.1891   0.13692
## Residuals 2225.4 21

# rstatix Type III
res <- anova_test(y ~ x + block + group, data = df, type = 3)
get_anova_table(res)
## ANOVA Table (type III tests)
##   Effect DFn DFd      F        p p<.05   ges
## 1      x   1  21 58.264 1.73e-07     * 0.735
## 2  block  11  21  3.230 1.00e-02     * 0.629
## 3  group   2  21  2.189 1.37e-01       0.173
```

## 结果解读指南

| 输出项 | 含义 | 解读要点 |
|--------|------|---------|
| `x` 行 | 协变量效应 | p < 0.05 说明协变量与因变量线性关系显著，纳入正确 |
| `group` 行 | 分组主效应（修正后） | p < 0.05 说明校正协变量后组间差异有统计学意义 |
| `Df` | 自由度 | group 自由度为 k-1（k=组数）；残差自由度为 N-k-协变量数 |
| `Sum Sq` | 平方和 | 反映对应因素的变异大小，可比较不同因素相对重要性 |
| `Mean Sq` | 均方 = Sum Sq / Df | 标准化后的变异量度 |
| `F value` | F 统计量 = group_MS / residual_MS | 值越大，组间差异相对于组内变异越大 |
| `Pr(>F)` | p 值 | 在原假设（各组修正均数相等）成立时观察到该 F 值的概率 |
| `ges`（rstatix） | 广义 eta 平方（generalized eta-squared） | 效应量指标，值越大效应越大 |

**常用表述模板：**

- 组间差异显著：「在控制协变量 XXX 后，各组间的 YYY 差异具有统计学意义（F(df1, df2) = F值，p < 0.05）」
- 组间差异不显著：「在控制协变量 XXX 后，尚不能认为各组间的 YYY 存在差异（F(df1, df2) = F值，p = 0.XX）」

## 常见问题与注意事项

**Q1: 为什么 aov() 公式中协变量要放在分组变量前面？**

A: R 的 `aov()` 默认使用 Type I（序贯型）平方和计算，结果依赖于因素在公式中的顺序。将协变量放在前面意味着先扣除协变量解释的变异，再检验分组效应。这与 ANCOVA 的逻辑一致：先校正基线差异，再比较组别。

**Q2: Type I、II、III 平方和有什么区别？何时用哪种？**

A:
- **Type I（序贯型）**：因素顺序影响结果，适合完全随机设计 ANCOVA（协变量优先），使用 `aov()` 或 `anova_test(type = 1)`
- **Type II**：检验每个主效应时控制其它同阶项（不含交互），适合无交互项的模型，使用 `car::Anova()`
- **Type III**：检验每个效应时控制所有其他效应，适用于不平衡设计或有交互的模型，使用 `anova_test(type = 3)`

在随机区组设计中，三种 SS 结果可能不同，Type III 通常仅依赖 `rstatix::anova_test()` 获得。

**Q3: 协变量只能是连续变量吗？可以有多个协变量吗？**

A: ANCOVA 的协变量原则上是连续型变量。可以有多个协变量，公式写为 `aov(y ~ x1 + x2 + group)`，此时也称为多重协方差分析。但协变量越多，对样本量的要求越高，且多重共线性问题需注意。

**Q4: SPSS 结果和 R 结果不一致怎么办？**

A: 最常见的原因是平方和类型不同。SPSS 默认使用 Type III SS，与 R 的 `aov()` 默认 Type I 不同。使用 `rstatix::anova_test(type = 3)` 或 `car::Anova(..., type = 3)` 可获得与 SPSS 一致的结果。

**Q5: 如果回归线不平行怎么办？**

A: 若 `x:group` 交互项 p < 0.05，拒绝斜率齐性假设，说明协变量与因变量的关系因组而异。此时：
- 不可直接报告常规 ANCOVA 结果
- 可分别对每组做回归分析，定性比较
- 若分组因素本身关注度不高，可改为多元回归加交互项探索
- 可考虑使用倾向性评分匹配（PSM）等其他校正方法
