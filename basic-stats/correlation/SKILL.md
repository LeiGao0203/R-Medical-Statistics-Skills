---
name: medical-stat-correlation
description: "R语言医学统计：双变量相关与简单线性回归。涵盖Pearson相关系数、Spearman秩相关系数、简单线性回归、决定系数R²、曲线拟合。TRIGGER when user mentions 相关分析、相关系数、Pearson、Spearman、两变量关系、散点图、线性回归、lm()，or asks about the relationship between two continuous variables. SKIP for 多变量回归（多元线性回归）、分类结局（Logistic回归）、偏相关与典型相关分析。"
---

# 双变量回归与相关 (Correlation & Simple Linear Regression)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用双变量回归与相关：**

- 探索两个连续型变量之间是否存在关联及其方向和密切程度（相关分析）
- 以一个变量预测另一个变量的取值（简单线性回归）
- 两变量均为连续型数值变量（如体重与肾体积、药物浓度与效应值）
- 散点图呈现线性趋势时进行直线相关/回归分析
- 散点图呈现曲线趋势时进行曲线拟合（对数转换、指数转换等）
- 比较两条回归直线是否平行（如正常组与疾病组的回归关系比较）

**何时不使用本章方法：**

| 情况 | 替代方法 |
|------|----------|
| 3个及以上自变量的回归 | 多元线性回归（`lm()` 加多个自变量） |
| 二分类/多分类结局变量 | Logistic 回归 |
| 需要控制第三变量的混杂作用 | 偏相关分析（`ppcor::pcor()`） |
| 两变量为分类资料 | 卡方检验及其关联性指标（Cramer's V, phi 系数） |
| 多变量之间的相关矩阵 | 典型相关分析（CCA） |

## 前置条件

**R 包依赖：**

```r
# 基础包（无需安装）
# stats     —— lm(), cor(), cor.test(), predict()
# graphics  —— plot(), abline()

# 可选安装
install.packages("ggplot2")   # 散点图 + 回归线
install.packages("broom")     # 将模型输出整理为tidy格式
install.packages("foreign")   # 读取SPSS数据
```

**统计假设（简单线性回归的前提条件）：**

1. **线性**（linearity）：因变量 Y 与自变量 X 之间存在线性关系。通过散点图判断
2. **独立性**（independence）：各观测值之间相互独立（由研究设计保证）
3. **方差齐性**（homoscedasticity）：对于每个 X 值，Y 的方差应大致相等。通过残差图判断
4. **正态性**（normality of residuals）：对于每个 X 值，Y 值呈正态分布。大样本下对偏离不敏感

**Pearson 相关的附加条件：**
- 两变量均服从正态分布（双正态），或近似正态
- 当任一变量不满足正态或为等级资料时，应使用 Spearman 秩相关

**数据格式要求：**
- 两列数值型变量，每行代表一个观测对象
- 用于回归时明确区分自变量 X（预测因子）和因变量 Y（结局变量）

## 方法选择决策树

```
你的数据情况 →
├── 只想看两变量的关联性（无因果方向要求）
│   ├── 两变量均满足正态性 → 使用 Pearson相关 → cor.test(x, y, method = "pearson")
│   ├── 任一变量不满足正态性 → 使用 Spearman秩相关 → cor.test(x, y, method = "spearman")
│   └── 变量为等级资料/有序分类 → 使用 Spearman秩相关 → cor.test(x, y, method = "spearman")
│
├── 有明确的因果方向（X 影响 Y），或想做预测
│   ├── 散点图呈线性趋势 → 使用 简单线性回归 → lm(y ~ x)
│   │   ├── 需要回归系数检验 → summary(lm(...))
│   │   ├── 需要系数置信区间 → broom::tidy(lm(...), conf.int = TRUE)
│   │   ├── 需要估计总体均数CI → predict(lm(...), interval = "confidence")
│   │   └── 需要预测个体Y值PI → predict(lm(...), interval = "prediction")
│   │
│   └── 散点图呈曲线趋势 → 使用 曲线拟合
│       ├── 趋势像对数函数 → lm(y ~ log10(x))
│       ├── 趋势像指数函数 → lm(log(y) ~ x)
│       └── 更复杂曲线 → 多项式回归 or 样条回归（见高级章节）
│
└── 需要比较两组回归直线是否平行/重合
    └── 使用 交互项模型比较 → lm(y ~ x * group) vs lm(y ~ x + group)
        └── 交互项P值判断是否平行，组别系数P值判断截距是否相等
```

## 标准工作流

### 步骤1：数据准备与探索

读入数据，检查结构，做散点图以初步判断两变量关系趋势。

```r
# 构造数据（课本例9-5：体重与双肾体积的关系）
df <- data.frame(
  weight = c(43, 74, 51, 58, 50, 65, 54, 57, 67, 69, 80, 48, 38, 85, 54),
  kv = c(217.22, 316.18, 231.11, 220.96, 254.70, 293.84, 263.28,
         271.73, 263.46, 276.53, 341.15, 261.00, 213.20, 315.12, 252.08)
)

str(df)
## 'data.frame':    15 obs. of  2 variables:
##  $ weight: num  43 74 51 58 50 65 54 57 67 69 ...
##  $ kv    : num  217 316 231 221 255 ...

# 散点图初探
library(ggplot2)
ggplot(df, aes(weight, kv)) +
  geom_point(size = 4) +
  labs(x = "体重(kg) X", y = "双肾体积(ml) Y") +
  theme_classic()
```

### 步骤2：前提条件检验

正态性检验（仅 Pearson 相关和回归残差需要）。小样本时查看散点图趋势是否线性。

```r
# 正态性检验
shapiro.test(df$weight)
shapiro.test(df$kv)

# 对于回归分析，拟合后检查残差正态性
fit <- lm(kv ~ weight, data = df)
shapiro.test(residuals(fit))
```

### 步骤3：执行统计分析

根据决策树选择方法：想预测/有因果关系用 `lm()`，纯关联用 `cor.test()`，非正态用 `method = "spearman"`。

### 步骤4：结果解读

- Pearson/Spearman：看 r/rho 的大小和正负方向，置信区间是否包含0，p值是否小于0.05
- 简单线性回归：看 R²（模型解释度）、回归系数（效应大小及方向）、F检验（整体模型显著性）、t检验（单个系数显著性）
- 曲线拟合：看转换后的 R² 和回归系数，判断拟合优劣

### 步骤5：结果报告

论文中报告格式示例："体重与双肾体积呈显著正相关，Pearson 相关系数为0.875（95% CI: 0.658–0.958, P < 0.001）。" 或 "以体重预测双肾体积的线性回归方程为 Ŷ = 1.662 + 0.139X，R² = 0.778，方程具有统计学意义（F = 20.97, P = 0.004）。"

## 代码示例

### 简单线性回归（课本例9-1）

```r
df9_1 <- data.frame(
  x = c(13, 11, 9, 6, 8, 10, 12, 7),
  y = c(3.54, 3.01, 3.09, 2.48, 2.56, 3.36, 3.18, 2.65)
)

fit <- lm(y ~ x, data = df9_1)
summary(fit)
## 
## Call:
## lm(formula = y ~ x, data = df9_1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.21500 -0.15937 -0.00125  0.09583  0.30667 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  1.66167    0.29700   5.595  0.00139 **
## x            0.13917    0.03039   4.579  0.00377 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.197 on 6 degrees of freedom
## Multiple R-squared:  0.7775, Adjusted R-squared:  0.7404 
## F-statistic: 20.97 on 1 and 6 DF,  p-value: 0.003774
```

回归方程：截距 = 1.662，斜率 = 0.139（即 Ŷ = 1.662 + 0.139X）。F 检验 P = 0.004，方程整体有统计学意义。`Multiple R-squared: 0.7775` 表示模型解释了 Y 总变异的 77.75%。

### 回归系数的置信区间（例9-3 + 例9-4）

```r
library(broom)

# 回归系数及其95%CI
broom::tidy(fit, conf.int = TRUE)
## # A tibble: 2 x 7
##   term        estimate std.error statistic p.value conf.low conf.high
##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
## 1 (Intercept)    1.66     0.297       5.59 0.00139   0.935      2.39 
## 2 x              0.139    0.0304      4.58 0.00377   0.0648     0.214

# 当 X=12 时，总体均数 μ_Y|X 的95%可信区间
new_x <- data.frame(x = 12)
predict(fit, newdata = new_x, interval = "confidence", level = 0.95)
##        fit      lwr      upr
## 1 3.331667 3.079481 3.583852

# 当 X=12 时，个体 Y 值的95%预测区间
predict(fit, newdata = new_x, interval = "prediction", level = 0.95)
##        fit      lwr      upr
## 1 3.331667 2.787731 3.875602
```

> 可信区间（confidence interval）是总体均数 μ_Y|X 的范围，预测区间（prediction interval）是个体值 Y 的范围。预测区间宽度 > 可信区间。

### Pearson 相关分析（课本例9-5）

```r
df <- data.frame(
  weight = c(43, 74, 51, 58, 50, 65, 54, 57, 67, 69, 80, 48, 38, 85, 54),
  kv = c(217.22, 316.18, 231.11, 220.96, 254.70, 293.84, 263.28,
         271.73, 263.46, 276.53, 341.15, 261.00, 213.20, 315.12, 252.08)
)

# 仅计算相关系数 r
cor(df$weight, df$kv)
## [1] 0.8754315

# 计算相关系数 + 假设检验 + 置信区间
cor.test(~ weight + kv, data = df)
## 
##  Pearson's product-moment correlation
## 
## data:  weight and kv
## t = 6.5304, df = 13, p-value = 1.911e-05
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.6584522 0.9580540
## sample estimates:
##       cor 
## 0.8754315
```

结论：r = 0.875，95% CI (0.658, 0.958)，p < 0.001，具有统计学意义。

### 散点图 + 回归线可视化

```r
library(ggplot2)

ggplot(df, aes(weight, kv)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_vline(xintercept = mean(df$weight), linetype = 2) +
  geom_hline(yintercept = mean(df$kv), linetype = 2) +
  labs(x = "体重(kg) X", y = "双肾体积(ml) Y") +
  theme_classic()
```

> 图中虚线为各自均数，辅助判断协变方向。Mean-corrected 散点中，点在一/三象限为正向协变。

### R² 的计算

```r
# R² = r²（对于简单线性回归而言）
summary(lm(weight ~ kv, data = df))
## Multiple R-squared:  0.7664, Adjusted R-squared:  0.7484
## F-statistic: 42.65 on 1 and 13 DF,  p-value: 1.911e-05

# 验证：R² = r²
0.8754315^2
## [1] 0.76638
```

### Spearman 秩相关（课本例9-8）

```r
library(foreign)
df9_8 <- foreign::read.spss("datasets/例09-08.sav", to.data.frame = TRUE)

# 计算 Spearman 秩相关系数
cor(df9_8$x, df9_8$y, method = "spearman")
## [1] 0.9050568

# 同时获得 P 值
cor.test(df9_8$x, df9_8$y, method = "spearman")
## 
##  Spearman's rank correlation rho
## 
## data:  df9_8$x and df9_8$y
## S = 92, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9050568

# 带连续性校正（cor.test 默认在计算P值时进行校正）
cor.test(df9_8$x, df9_8$y, method = "spearman", continuity = TRUE)
## 
##  Spearman's rank correlation rho
## 
## data:  df9_8$x and df9_8$y
## S = 92, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9050568
```

### 两条回归直线的比较（例9-1 vs 例9-9）

```r
df9_1 <- data.frame(
  x = c(13, 11, 9, 6, 8, 10, 12, 7),
  y = c(3.54, 3.01, 3.09, 2.48, 2.56, 3.36, 3.18, 2.65)
)

df9_9 <- foreign::read.spss("datasets/例09-09.sav", to.data.frame = TRUE)

# 直接 anova(lm1, lm2) 会报错（样本量不同）
# 正确方法：合并数据 + 交互项检验

df9_1$group <- "group9_1"
df9_9$group <- "group9_9"

df9 <- rbind(df9_1, df9_9[, -1])
df9$group <- factor(df9$group)

# 不包含交互项模型
model_no_interaction <- lm(y ~ x + group, data = df9)
# 包含交互项模型
model_interaction <- lm(y ~ x * group, data = df9)

# 比较交互项显著性 → 判断两直线是否平行
anova(model_no_interaction, model_interaction)
## Analysis of Variance Table
## 
## Model 1: y ~ x + group
## Model 2: y ~ x * group
##   Res.Df     RSS Df Sum of Sq      F Pr(>F)
## 1     15 0.62211                           
## 2     14 0.59101  1  0.031103 0.7368 0.4052

# 交互项 P = 0.405 > 0.05，还不能认为两条回归直线不平行
```

### 曲线拟合

**例9-11：对数转换自变量 X**

```r
df9_11 <- foreign::read.spss("datasets/例09-11.sav", to.data.frame = TRUE)

# 散点图
library(ggplot2)
ggplot(df9_11, aes(x, y)) + geom_point(size = 4)

# 对数转换后的散点图
ggplot(df9_11, aes(log10(x), y)) + geom_point(size = 4)

# 对数转换后的直线回归
f9_11 <- lm(y ~ log10(x), data = df9_11)
summary(f9_11)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  110.060      4.095   26.88 0.000113 ***
## log10(x)      36.115      2.968   12.17 0.001195 ** 
## ---
## Multiple R-squared:  0.9801, Adjusted R-squared:  0.9735
```

**例9-12：对数转换因变量 Y**

```r
df9_12 <- foreign::read.spss("datasets/例09-12.sav", to.data.frame = TRUE)

# 对数转换后的散点图
ggplot(df9_12, aes(x, log(y))) + geom_point(size = 4)

f9_12 <- lm(log(y) ~ x, data = df9_12)
summary(f9_12)
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  4.037159   0.084103   48.00 5.08e-16 ***
## x           -0.037974   0.002284  -16.62 3.86e-10 ***
## ---
## Multiple R-squared:  0.9551, Adjusted R-squared:  0.9516
## F-statistic: 276.4 on 1 and 13 DF,  p-value: 3.858e-10
```

## 结果解读指南

**lm() summary 输出各组成部分的解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `(Intercept) Estimate` | 截距 b₀ | X=0 时 Y 的预测值，需结合专业意义判断是否合理 |
| `x Estimate` | 回归系数 b₁（斜率） | X 每增加1单位，Y 平均改变 b₁ 单位；符号表方向 |
| `Std. Error` | 回归系数的标准误 | 估计的精度，越小越精确 |
| `t value` (对系数) | 回归系数=0 的 t 检验 | t = Estimate / Std.Error |
| `Pr(>|t|)` (对系数) | 回归系数=0 的 p 值 | p < 0.05 → 该变量对 Y 有统计学显著影响 |
| `Multiple R-squared` | 决定系数 R² | 模型解释 Y 总变异的百分比（0~1，越大越好） |
| `Adjusted R-squared` | 调整 R² | 校正了自变量个数，简单回归中略小于 R² |
| `F-statistic` | 整体模型检验 F 值 | 检验所有自变量系数是否同时为0 |
| `F-statistic p-value` | 模型整体 p 值 | p < 0.05 → 回归方程整体有统计学意义 |

**cor.test() 输出各组成部分的解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `cor` / `rho` | 相关系数 | 取值 -1 ~ 1，0 为无线性相关；符号表示方向 |
| `t` | 相关系数的 t 检验统计量 | 用于计算 p 值 |
| `df` | 自由度 | n - 2 |
| `p-value` | 相关系数是否 ≠ 0 的检验 | p < 0.05 → 相关系数有统计学意义 |
| `95% confidence interval` | r 的95%置信区间 | 不包含0 等价于 p < 0.05；描述 r 的可能范围 |

**相关系数 r 的常用判断标准：**

| |r| 范围 | 相关程度 |
|--------|--------|
| 0.0 ~ 0.3 | 弱相关（微弱/无） |
| 0.3 ~ 0.5 | 低度相关 |
| 0.5 ~ 0.7 | 中度相关 |
| 0.7 ~ 0.9 | 高度相关 |
| 0.9 ~ 1.0 | 极高度相关 |

**r vs R² vs 回归系数 b 的区别：**

- **r（Pearson 相关系数）**：衡量两变量线性关联度，无因果方向。r 为正表示正相关，r 为负表示负相关
- **R²（决定系数）**：回归模型中 Y 被 X 解释的变异比例。对于简单线性回归，R² = r²
- **b₁（回归系数 / 斜率）**：有方向含义——X 每变化1单位，Y 期望改变 b₁ 单位。b₁ 和 r 符号一致，但 b₁ 的大小取决于变量单位

## 常见问题与注意事项

**Q1：Pearson 相关和 Spearman 秩相关怎么选？**

Pearson 相关要求两变量均服从正态分布（双正态），反映的是两者间的线性关联程度。Spearman 秩相关无需正态性假设，计算的是等级相关性，适用于：
- 任一变量不满足正态分布
- 变量为等级资料/有序分类
- 变量间存在单调但非线性的关系
- 数据存在离群值时更稳健

实际应用中若不确定分布情况，可两种方法都算，若结果方向一致则结论更可靠。

**Q2：r 值很大但 p > 0.05，怎么回事？**

可能原因：样本量极小。相关系数的假设检验效力高度依赖样本量，即使 r = 0.8，n = 5 时也可能 p > 0.05。p > 0.05 意味着不能排除总体中 ρ = 0 的可能性。建议结合置信区间一起报告。

**Q3：相关分析和回归分析先做哪个？**

通常先做散点图 → 相关分析（是否存在关联）→ 若存在关联且有因果关系假设 → 回归分析（量化因果关系+预测）。如果纯粹是探索性研究，两者可同时报告。

**Q4：R² = 0.78 算好模型吗？**

对于医学研究中的简单线性回归（仅1个自变量），R² > 0.7 已属较高水平的解释度。但模型好坏还需结合：
- 残差诊断（独立性、正态性、等方差性）
- 预测区间的实际宽度是否满足临床决策需要
- 研究领域的惯例（部分领域 R² > 0.3 即可接受）

**Q5：比较两条回归直线时直接 anova(fit1, fit2) 报错怎么办？**

错误 `models were not all fitted to the same size of dataset` 说明两个模型基于不同样本量的数据。正确做法：合并两个数据集，添加分组变量，用 `lm(y ~ x + group)` 和 `lm(y ~ x * group)` 两个模型，通过交互项的 F 检验（`anova(model_no_interaction, model_interaction)`）判断两直线是否平行。组别系数 P 值判断截距是否相等。

**Q6：曲线拟合 vs 直线回归，如何决策？**

- 先画散点图看趋势，不要盲目做直线回归
- 若趋势明显非直线（如对数曲线、指数曲线），直接做直线回归会低估关联
- 常用策略：假设数据生成机制（如化学反应、药代动力学模型）指导转换形式
- 简单线性回归的 R² 和曲线拟合的 R² 不能直接比较（因变量不同），只能在同一种因变量形式下比较

**常见错误提醒：**

1. 将相关分析结果直接解释为因果关系——相关 ≠ 因果，因果关系需要研究设计（如随机对照试验）保证
2. 对非线性关系强行使用 Pearson 相关——应先用散点图检查，必要时使用 Spearman 或曲线拟合
3. 简单线性回归忽略残差诊断——即使 R² 很高，残差存在模式也意味着模型设定有误
4. 混淆总体均数的可信区间和个体值的预测区间——前者描述 μ_Y|X 的不确定性，后者描述个体预测值 Y 的散布范围，后者一定比前者宽
5. 相关系数假显著——当散点图呈现 U 型或其他非线性形态时，即使实际存在强关系，Pearson r 仍可能接近于 0
