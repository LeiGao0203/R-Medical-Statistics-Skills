---
name: medical-stat-sphericity
description: "R语言医学统计：球对称检验。使用Mauchly检验判断重复测量方差分析的球对称假设是否成立，若违反则进行Greenhouse-Geisser或Huynh-Feldt校正。TRIGGER when user mentions 球对称、Mauchly检验、Greenhouse-Geisser、Huynh-Feldt、重复测量前提检验，or asks about sphericity assumption for repeated measures ANOVA. SKIP for 一般方差齐性检验（Bartlett/Levene）、t检验的前提检验。"
---

# 球对称检验 (Sphericity Test / Mauchly's Test)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

- **必须使用**：在进行重复测量方差分析（repeated measures ANOVA）之前，必须先行球对称检验，判断重复测量误差的协方差矩阵经正交对比变换后是否与单位矩阵成比例。
- **组内因素≥3个水平**：当重复测量因素有3个及以上水平时，球对称检验才有意义（2水平时自动满足）。
- **何时不用**：单次测量的组间比较（如独立样本t检验、单因素ANOVA）不需要球对称检验，此时应使用方差齐性检验（Levene检验/Bartlett检验）。

典型医学场景：
- 同一批患者在治疗前、治疗中期、治疗结束时的某指标测量值（纵向重复测量）
- 交叉设计试验中不同阶段的数据
- 同一受试对象在不同条件下的测量值

## 前置条件

```r
# 内置函数，无需额外安装包
# 如需读取SPSS数据：
install.packages("foreign")
```

**数据格式**：宽格式（wide format），每一行是一个受试对象，每一列是一个时间点/条件下的测量值。示例：

| 受试者 | t0   | t45  | t90  | t135 |
|--------|------|------|------|------|
| 1      | 5.32 | 5.32 | 4.98 | 4.65 |
| 2      | 5.32 | 5.26 | 4.93 | 4.70 |

**统计假设**：
- H0：数据满足球对称假设（各时间点间差值的方差相等）
- H1：数据不满足球对称假设
- 检验水准：常取 α = 0.05（也可取 α = 0.10，因为Mauchly检验在小样本时不够敏感）

## 方法选择决策树

```
你的数据情况 →
├── 单组（无分组因素），多次重复测量 → mauchly.test(lm(data ~ 1), X = ~ 1)
├── 多组（有分组因素），多次重复测量 → mauchly.test(lm(data ~ 1), M = ~ group + times, X = ~ times)
└── 数据为长格式 → 先用 tidyr::pivot_wider() 转为宽格式，或直接使用 ezANOVA() / anova_test() 同时输出球对称检验结果
```

**违反后的处理**：
```
Mauchly检验 p < 0.05 →
├── 首选 Greenhouse-Geisser 校正（ε < 0.75 时更保守，推荐默认使用）
└── 备选 Huynh-Feldt 校正（ε > 0.75 时可用，更宽松）
```

## 标准工作流

### 步骤1：数据准备与探索

将数据读入R，确认数据的行数（受试对象数）和列数（时间点数），并检查数据是否有缺失值。

### 步骤2：数据格式转换

`mauchly.test()` 要求输入为矩阵（matrix）格式。使用 `as.matrix()` 将数据框转为矩阵。如果数据有分组变量需要单独设置。

### 步骤3：执行球对称检验

单组数据使用 `mauchly.test(lm(data ~ 1), X = ~ 1)`。
多组数据使用 `mauchly.test(lm(data ~ 1), M = ~ group + times, X = ~ times)`。

### 步骤4：结果解读

查看W统计量和p值。若 p < 0.05，拒绝球对称假设，后续重复测量方差分析需使用校正后的自由度。

### 步骤5：结果报告（论文写法）

示例：「采用Mauchly球对称检验判断数据是否满足球对称假设，结果显示 W = 0.063, p = 0.008，拒绝球对称假设。因此，组内效应的F检验采用Greenhouse-Geisser法进行自由度校正。」

## 代码示例

### 示例1：单组数据（课本表12-3，8名受试者4个时间点的血糖值）

```r
library(foreign)

df <- foreign::read.spss("datasets/表12-3重复测量ANOVA.sav",
                         to.data.frame = TRUE, reencode = "utf-8")
str(df)
## 'data.frame':    8 obs. of  4 variables:
##  $ t0  : num  5.32 5.32 5.94 5.49 5.71 6.27 5.88 5.32
##  $ t45 : num  5.32 5.26 5.88 5.43 5.49 6.27 5.77 5.15
##  $ t90 : num  4.98 4.93 5.43 5.32 5.43 5.66 5.43 5.04
##  $ t135: num  4.65 4.7 5.04 5.04 4.93 5.26 4.93 4.48

df <- as.matrix(df)

mauchly.test(lm(df ~ 1), X = ~ 1)
## 
##  Mauchly's test of sphericity
##  Contrasts orthogonal to
##  ~1
## 
## data:  SSD matrix from lm(formula = df ~ 1)
## W = 0.06273, p-value = 0.008207
```

### 示例2：多组数据（课本例12-3，3种诱导方法，5个时间点血压值）

```r
library(foreign)

df1 <- foreign::read.spss("datasets/例12-03.sav", to.data.frame = TRUE)
str(df1)
## 'data.frame':    15 obs. of  7 variables:
##  $ No   : num  1 2 3 4 5 6 7 8 9 10 ...
##  $ group: Factor w/ 3 levels "A","B","C": 1 1 1 1 1 2 2 2 2 2 ...
##  $ t0   : num  120 118 119 121 127 121 122 128 117 118 ...
##  $ t1   : num  108 109 112 112 121 120 121 129 115 114 ...
##  $ t2   : num  112 115 119 119 127 118 119 126 111 116 ...
##  $ t3   : num  120 126 124 126 133 131 129 135 123 123 ...
##  $ t4   : num  117 123 118 120 126 137 133 142 131 133 ...

# 按组别合并测量值，转为矩阵
df2 <- as.matrix(cbind(df1[1:5, 3:7], df1[6:10, 3:7], df1[11:15, 3:7]))

# 创建分组和时间因子（顺序必须与数据列一致）
times <- ordered(rep(1:5, 3))
group <- factor(rep(c("A", "B", "C"), each = 5))

mauchly.test(lm(df2 ~ 1), M = ~ group + times, X = ~ times)
## 
##  Mauchly's test of sphericity
##  Contrasts orthogonal to
##  ~times
## 
##  Contrasts spanned by
##  ~group + times
## 
## data:  SSD matrix from lm(formula = df2 ~ 1)
## W = 0.427, p-value = 0.279
```

### 示例3：便捷方法——使用 ezANOVA 自动输出球对称检验

```r
install.packages("ez")
library(ez)

# 需先将数据转为长格式
df_long <- tidyr::pivot_longer(df, cols = starts_with("t"),
                               names_to = "time", values_to = "value")

ezANOVA(data = df_long, dv = value, wid = id,
        within = time, detailed = TRUE)
# 输出中自动包含 Mauchly's Test 及 Greenhouse-Geisser / Huynh-Feldt 校正结果
```

## 结果解读指南

| 输出项 | 含义 | 判读标准 |
|--------|------|----------|
| W (Mauchly's W) | 球对称检验统计量，取值范围 (0, 1]。W = 1 表示完全满足球对称；W 越小，偏离越严重。 | 配合 p 值使用 |
| p-value | 球对称假设的显著性检验 | p < 0.05 → 拒绝球对称假设，需校正 |
| Greenhouse-Geisser ε | 校正系数 epsilon，取值范围 [1/(k-1), 1]，k 为重复测量水平数。ε 越小违反越严重。 | 乘以原自由度得到校正后的自由度 |
| Huynh-Feldt ε | 校正系数，通常比 G-G 的 ε 略大（更不保守）。当真实 ε 接近 1 时 H-F 偏差更小。 | ε > 0.75 时可考虑使用 H-F 校正 |

**论文报告示例**：
- 满足球对称：「Mauchly球对称检验显示数据满足球对称假设（W = 0.427, p = 0.279），因此采用未经校正的一元方差分析结果。」
- 违反球对称：「Mauchly球对称检验显示数据不满足球对称假设（W = 0.063, p = 0.008），故采用Greenhouse-Geisser法对组内效应的自由度进行校正（ε = 0.512）。」

## 常见问题与注意事项

**Q1：Mauchly检验的样本量很小时可靠吗？**
样本量较小时检验功效不足，容易得出「满足球对称」的假阴性结论。此时可适当放宽检验水准至 α = 0.10，或直接选择对球对称假设稳健的多变量方法（MANOVA），或直接使用校正结果。

**Q2：2个时间点的重复测量需要做球对称检验吗？**
不需要。当重复测量因素仅有2个水平时，球对称假设自动满足，因为只有1个差值方差需要估计。

**Q3：Greenhouse-Geisser 和 Huynh-Feldt 校正怎么选？**
- 当 ε < 0.75 时，建议使用 Greenhouse-Geisser 校正，更为保守可靠
- 当 ε > 0.75 时，Huynh-Feldt 校正的 I 类错误控制更接近名义水平
- 实际应用中直接使用 G-G 校正是稳妥选择，大多数统计软件默认采用此法

**Q4：SPSS和R的球对称检验结果一致吗？**
一致。两者的检验方法相同（Mauchly, 1940），输出W和p值理论上应完全相同。SPSS在「Repeated Measures」对话框中自动输出，R需手动调用 `mauchly.test()` 或通过 `ezANOVA()` / `anova_test()` 获取。

**Q5：mauchly.test() 的函数签名让我困惑，怎么用？**
核心理解：第一个参数是 `lm()` 拟合的模型（一般为 `lm(data ~ 1)`），`X = ~ 1` 用于单组情景，多组时用 `X = ~ times` 指定组内因素、`M = ~ group + times` 指定完整模型公式。

**注意事项**：
- `mauchly.test()` 要求输入为数值**矩阵**（matrix），不是数据框。使用 `as.matrix()` 前请确认所有列均为数值型。
- 数据中不应有缺失值。如有缺失，应先进行缺失值处理或删除不完整观测。
- 对于大数据或复杂设计，建议直接使用 `ez::ezANOVA()` 或 `rstatix::anova_test()`，它们会自动计算并输出Mauchly检验结果和校正值。
