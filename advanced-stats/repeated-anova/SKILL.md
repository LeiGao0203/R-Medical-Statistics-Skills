---
name: medical-stat-repeated-anova
description: "R语言医学统计：重复测量方差分析。涵盖单因素重复测量、两因素重复测量（含一个组间因素）、球对称检验及Greenhouse-Geisser/Huynh-Feldt校正，以及组间多重比较、时间趋势正交多项式分析和时间点事前检验。TRIGGER when user mentions 重复测量、重复测量方差分析、repeated measures ANOVA、within-subject design、不同时间点比较、前后测量比较、受试者内设计、aov Error，or asks about analyzing repeated measurements over time. SKIP for 独立样本方差分析、广义估计方程(GEE)、混合效应模型/多水平模型、协方差分析(ANCOVA)、配对t检验（仅两个时间点无交互效应）。"
---

# 重复测量方差分析 (Repeated Measures ANOVA)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用重复测量方差分析：**

- 同一组受试者在多个时间点重复测量同一指标（如：治疗前/治疗后血压变化，麻醉诱导后t0-t4五个时间点的血压监测）
- 含有一个组间因素（如不同治疗方案）和一个组内因素（如测量时间）的两因素设计
- 仅含有组内因素（单因素重复测量）的设计，如多个不同条件/处理对同一批受试者的效应
- 结局变量为连续型数值变量，重复测量次数≥2
- 需同时评估组间效应、时间效应及其交互效应

**何时不用重复测量方差分析：**

| 情况 | 替代方法 |
|------|----------|
| 仅两个时间点、无分组因素，只需前后比较 | 配对t检验 |
| 数据严重违反球对称假设且校正后仍不满足 | 混合效应模型（lme4）或多水平模型 |
| 结局变量为非正态分布的分类变量或计数资料 | 广义估计方程（GEE） |
| 需要控制协变量 | 协方差分析（ANCOVA） |
| 各组为独立样本，非同一受试者重复测量 | 独立样本方差分析 / 多因素方差分析 |

**常见医学研究场景：**

- 临床试验中不同治疗方案在多个随访时间点的疗效比较
- 麻醉诱导研究中不同诱导方法在各时点对生理指标的影响
- 药物浓度在不同时间点对不同组别受试者的动态变化
- 康复治疗中不同干预手段前后多个时间点的功能指标评估

## 前置条件

**R包依赖：**

```r
# 基础包（无需安装）
# stats —— aov(), summary(), contr.poly(), interaction.plot()
# graphics —— boxplot()

# 需要安装的包
install.packages("tidyverse")   # dplyr + tidyr 数据清洗与转换
install.packages("reshape2")    # melt() 数据转换，长宽格式切换
install.packages("rstatix")     # anova_test() 含球对称检验及校正、t_test() 事前检验
install.packages("PMCMRplus")   # lsdTest() 组间多重比较
install.packages("ggplot2")     # 高级绘图
install.packages("foreign")     # read.spss() 读取SPSS格式数据
```

**统计假设：**

1. **正态性**（normality）：各处理水平下的因变量应近似服从正态分布；当样本量较大时有一定稳健性
2. **球对称性**（sphericity）：不同时间点测量值之差分的方差相等。使用Mauchly球对称检验判断；不满足时使用Greenhouse-Geisser（GGe）或Huynh-Feldt（HFe）ε校正
3. **独立性**（independence）：不同受试者之间的观测值相互独立（由研究设计保证）
4. **数据类型**：因变量为连续型数值变量，组内因素（时间）为多水平因子，组间因素为分类因子

**数据格式要求：**

- **长格式**（long format）：每行代表一个受试者在某个时间点的一次测量
  - 受试者ID列（factor型）
  - 组间因素列（factor型）：如治疗分组
  - 组内因素/时间列（factor型）：如t0、t1、t2……
  - 因变量列（numeric型）：测量值

```
##   No group times  hp
## 1  1     A    t0 120
## 2  1     A    t1 108
## 3  2     A    t0 118
## 4  2     A    t1 109
```

> 原始数据通常为宽格式（每个时间点一列），分析前必须用 `pivot_longer()` 或 `melt()` 转为长格式。

## 方法选择决策树

```
你的数据情况 →
├── 仅组内因素（单因素重复测量）
│   ├── 2个水平 → 配对t检验 或 aov(y ~ time + Error(ID/time))
│   └── ≥3个水平 → aov(y ~ time + Error(ID/time))，需球对称检验，不满足时用GG/HF校正
│
├── 组内因素 + 组间因素（两因素重复测量）
│   ├── 两因素两水平（如：治疗前后×两组）→ aov(y ~ time * group + Error(ID/time))
│   └── 两因素多水平（如：5个时点×3组）→ 同上，需球对称检验及多重比较
│
├── 需多重比较
│   ├── 组间比较 → PMCMRplus::lsdTest() / TukeyHSD()（与其他方差分析一致）
│   ├── 时间趋势探索 → 正交多项式变换 + aov(..., split=...)
│   └── 时间点两两比较 → rstatix::t_test(..., paired=T, ref.group="t0")
│
└── 球对称不满足（Mauchly p < 0.05）
    ├── ε > 0.75 → 使用Huynh-Feldt校正
    ├── ε < 0.75 或不知 → 使用Greenhouse-Geisser校正（更保守）
    └── 严重偏离 → 改用混合效应模型（lme4::lmer()）或多水平模型
```

## 标准工作流

### 步骤1：数据准备与探索

```
1. 读取数据 → foreign::read.spss() 或 read.csv()
2. 宽格式转长格式 → tidyverse::pivot_longer() 或 reshape2::melt()
3. 因子化关键变量 → ID、分组、时间均转为 factor 类型
4. 描述性探索 → 分组计算各时点均值，绘制交互效应图或箱线图
```

### 步骤2：球对称检验

```
1. 使用 rstatix::anova_test() 自动输出 Mauchly 球对称检验结果
2. 若 p > 0.05 → 球对称满足，直接使用未校正的 F 值和 p 值
3. 若 p < 0.05 → 球对称不满足，查阅校正结果（GGe 和 HFe）
```

### 步骤3：执行重复测量方差分析

```
1. 基础方法：aov(dv ~ within * between + Error(ID/within), data = ...)
2. 推荐方法：rstatix::anova_test(dv, wid, within, between)
   - 自动输出两因素（组内+组间）的ANOVA表
   - 同时输出Mauchly球对称检验及GG/HF校正
3. summary() 查看结果
```

### 步骤4：多重比较（若主效应或交互效应显著）

```
1. 组间比较：PMCMRplus::lsdTest() 或 TukeyHSD()
2. 时间趋势：正交多项式 contr.poly() + aov(..., split=...)
3. 时间点两两比较：rstatix::t_test(..., paired=T, ref.group="baseline")
```

### 步骤5：结果可视化与报告

```
1. interaction.plot() 绘制交互效应轮廓图
2. boxplot() 绘制分组×时间箱线图
3. ggplot2 绘制分组均值趋势折线图
4. 论文中报告：各时点均数±标准差、F值、p值、球对称检验W及p值，校正方法及调整后自由度与p值
```

## 代码示例

### 示例1：两因素两水平（治疗前后×两组，课本例12-1）

```r
library(tidyverse)
library(foreign)

# 读取数据
df12_1 <- foreign::read.spss("datasets/12-1.sav", to.data.frame = T)

str(df12_1)
## 'data.frame':    20 obs. of  5 variables:
##  $ n    : num  1 2 3 4 5 6 7 8 9 10 ...
##  $ x1   : num  130 124 136 128 122 118 116 138 126 124 ...
##  $ x2   : num  114 110 126 116 102 100 98 122 108 106 ...
##  $ group: Factor w/ 2 levels "处理组","对照组": 1 1 1 1 1 1 1 1 1 1 ...

# 宽格式转长格式
df12_11 <- df12_1[,1:4] %>%
  pivot_longer(cols = 2:3, names_to = "time", values_to = "hp") %>%
  mutate_if(is.character, as.factor)

df12_11$n <- factor(df12_11$n)

head(df12_11)
## # A tibble: 6 × 4
##   n     group  time     hp
##   <fct> <fct>  <fct> <dbl>
## 1 1     处理组 x1      130
## 2 1     处理组 x2      114
## 3 2     处理组 x1      124
## 4 2     处理组 x2      110
## 5 3     处理组 x1      136
## 6 3     处理组 x2      126

# 重复测量方差分析
f1 <- aov(hp ~ time * group + Error(n/time), data = df12_11)
summary(f1)
##
## Error: n
##           Df Sum Sq Mean Sq F value Pr(>F)
## group      1  202.5   202.5   1.574  0.226
## Residuals 18 2315.4   128.6
##
## Error: n:time
##            Df Sum Sq Mean Sq F value   Pr(>F)
## time        1 1020.1  1020.1   55.01 7.08e-07 ***
## time:group  1  348.1   348.1   18.77 0.000401 ***
## Residuals  18  333.8    18.5

# 交互效应图
with(df12_11,
     interaction.plot(time, group, hp, type = "b", col = c("red","blue"),
                      pch = c(12,16), main = "两因素两水平重复测量方差分析"))
```

### 示例2：两因素多水平（5个时点×3组，课本例12-3）

```r
library(tidyverse)
library(foreign)

df12_3 <- foreign::read.spss("datasets/例12-03.sav", to.data.frame = T,
                             reencode = "utf-8")

str(df12_3)
## 'data.frame':    15 obs. of  7 variables:
##  $ No   : num  1 2 3 4 5 6 7 8 9 10 ...
##  $ group: Factor w/ 3 levels "A","B","C": 1 1 1 1 1 2 2 2 2 2 ...
##  $ t0   : num  120 118 119 121 127 121 122 128 117 118 ...
##  $ t1   : num  108 109 112 112 121 120 121 129 115 114 ...
##  $ t2   : num  112 115 119 119 127 118 119 126 111 116 ...
##  $ t3   : num  120 126 124 126 133 131 129 135 123 123 ...
##  $ t4   : num  117 123 118 120 126 137 133 142 131 133 ...

# 转为长数据
df12_31 <- df12_3 %>%
  pivot_longer(cols = 3:7, names_to = "times", values_to = "hp")

df12_31$No <- factor(df12_31$No)
df12_31$times <- factor(df12_31$times)

# 基础方法
f2 <- aov(hp ~ times * group + Error(No/(times)), data = df12_31)
summary(f2)
##
## Error: No
##           Df Sum Sq Mean Sq F value Pr(>F)
## group      2  912.2   456.1   5.783 0.0174 *
## Residuals 12  946.5    78.9
##
## Error: No:times
##             Df Sum Sq Mean Sq F value   Pr(>F)
## times        4 2336.5   584.1   106.6  < 2e-16 ***
## times:group  8  837.6   104.7    19.1 1.62e-12 ***
## Residuals   48  263.1     5.5

# rstatix 方法（推荐，含球对称检验及校正）
library(rstatix)

anova_test(data = df12_31,
           dv = hp,
           wid = No,
           within = times,
           between = group)
## ANOVA Table (type II tests)
##
##        Effect DFn DFd       F        p p<.05   ges
## 1       group   2  12   5.783 1.70e-02     * 0.430
## 2       times   4  48 106.558 3.02e-23     * 0.659
## 3 group:times   8  48  19.101 1.62e-12     * 0.409
##
## $`Mauchly's Test for Sphericity`
##        Effect     W     p p<.05
## 1       times 0.293 0.178
## 2 group:times 0.293 0.178
##
## $`Sphericity Corrections`
##        Effect   GGe      DF[GG]    p[GG] p[GG]<.05   HFe      DF[HF]    p[HF] p[HF]<.05
## 1       times 0.679 2.71, 32.58 1.87e-16         * 0.896 3.59, 43.03 4.65e-21         *
## 2 group:times 0.679 5.43, 32.58 4.26e-09         * 0.896 7.17, 43.03 2.04e-11         *

# 可视化
with(df12_31,
     interaction.plot(times, group, hp, type = "b",
                      col = c("red","blue","green"),
                      pch = c(12,16,20),
                      main = "两因素多水平重复测量方差分析"))
```

### 示例3：组间多重比较（LSD法，基于例12-3）

```r
library(reshape2)
library(PMCMRplus)

# 数据转换（使用reshape2::melt）
df.l <- melt(df12_3, id.vars = c("No","group"),
             variable.name = "times", value.name = "hp")
df.l$No <- factor(df.l$No)

# 查看各组均值
df.l |> group_by(group) |> summarise(mm = mean(hp))
##   group    mm
## 1 A      120.
## 2 B      124.
## 3 C      128.

# LSD组间两两比较
summary(lsdTest(hp ~ group, data = df.l))
##            t value  Pr(>|t|)
## B - A == 0   2.175 0.0329218   *
## C - A == 0   3.860 0.0002446 ***
## C - B == 0   1.686 0.0962097   .
```

### 示例4：时间趋势正交多项式分析

```r
# 正交多项式对比矩阵
contrasts(df.l$times) <- contr.poly(5)

# A组时间趋势
f1 <- aov(hp ~ times, data = df.l[df.l$group == "A",])
summary(f1,
        split = list(times = list(liner = 1, quadratic = 2,
                                   cubic = 3, biquadrate = 4)))
##                     Df Sum Sq Mean Sq F value   Pr(>F)
## times                4  475.4   118.9   5.580 0.003486 **
##   times: liner       1   84.5    84.5   3.967 0.060229 .
##   times: quadratic   1   26.4    26.4   1.240 0.278655
##   times: cubic       1  364.5   364.5  17.113 0.000511 ***
##   times: biquadrate  1    0.0     0.0   0.001 0.972627

# B组时间趋势
f2 <- aov(hp ~ times, data = df.l[df.l$group == "B",])
summary(f2, split = list(times = list(liner = 1, quadratic = 2,
                                       cubic = 3, biquadrate = 4)))
##                     Df Sum Sq Mean Sq F value   Pr(>F)
## times                4 1017.0   254.3   9.757 0.000152 ***
##   times: liner       1  662.5   662.5  25.421 6.24e-05 ***
##   times: quadratic   1  296.2   296.2  11.367 0.003034 **
##   times: cubic       1    3.9     3.9   0.150 0.702229
##   times: biquadrate  1   54.4    54.4   2.088 0.163954

# C组时间趋势
f3 <- aov(hp ~ times + Error(No/times),
          data = df.l[df.l$group == "C",])
summary(f3, split = list(times = list(liner = 1, quadratic = 2,
                                       cubic = 3, biquadrate = 4)))
```

### 示例5：时间点事前检验（配对t检验，与t0基值比较）

```r
library(rstatix)

df.l |>
  group_by(group) |>
  t_test(hp ~ times, ref.group = "t0", paired = TRUE)
## # A tibble: 12 × 11
##    group .y.   group1 group2    n1    n2 statistic    df         p    p.adj p.adj.signif
##  1 A     hp    t0     t1         5     5     8.35      4 0.001     0.004   **
##  2 A     hp    t0     t2         5     5     1.77      4 0.152     0.304   ns
##  3 A     hp    t0     t3         5     5    -3.64      4 0.022     0.066   ns
##  4 A     hp    t0     t4         5     5     0.147     4 0.89      0.89    ns
##  5 B     hp    t0     t1         5     5     1.72      4 0.16      0.16    ns
##  6 B     hp    t0     t2         5     5     4.35      4 0.012     0.024   *
##  7 B     hp    t0     t3         5     5    -8.37      4 0.001     0.003   **
##  8 B     hp    t0     t4         5     5   -16.7       4 0.0000747 0.000299 ***
##  9 C     hp    t0     t1         5     5     1.44      4 0.223     0.292   ns
## 10 C     hp    t0     t2         5     5     4.75      4 0.009     0.028   *
## 11 C     hp    t0     t3         5     5    -5.12      4 0.007     0.028   *
## 12 C     hp    t0     t4         5     5    -1.80      4 0.146     0.292   ns
```

## 结果解读指南

### aov() 输出解读

`summary()` 输出包含**两个误差分层**的方差分析表：

**Error: ID（受试者间变异）—— 组间效应：**
| 列名 | 含义 |
|------|------|
| Df | 自由度（组间 + 残差） |
| Sum Sq | 离均差平方和 |
| Mean Sq | 均方（组间 / 组内） |
| F value | F统计量 = 组间均方 / 残差均方 |
| Pr(>F) | p值：p < 0.05 表示不同组间的平均差异具有统计学意义 |

**Error: ID:time（受试者内变异）—— 组内效应：**
| 行 | 含义 |
|----|------|
| time | 时间主效应：不同时间点测量值是否存在差异 |
| time:group | 交互效应：不同组别在时间上的变化趋势是否不同 |
| Residuals | 受试者内残差 |

**解读逻辑：**
1. 先看交互效应 `time:group`：若交互效应显著（p < 0.05），说明不同组随时间变化的趋势不同，需进一步做简单效应分析
2. 若交互效应不显著，再看主效应（time 和 group）

### rstatix::anova_test() 输出解读

与基础 `aov()` 一致，额外提供：

- **ges**：广义eta方（generalized eta-squared），效应量指标，值越大说明该因素的解释力越强
- **Mauchly球对称检验**：W值（0~1之间），p > 0.05 表示球对称假设成立；p < 0.05 说明球对称不满足，需参考校正结果
- **Sphericity Corrections**：
  - **GGe**（Greenhouse-Geisser ε）：更保守，适用于ε较小时
  - **HFe**（Huynh-Feldt ε）：适用于ε较大时（> 0.75）
  - 校正后的自由度含小数点，F值不变，p值按校正后的分母自由度重新计算
  - `p[GG]<.05` 列和 `p[HF]<.05` 列：校正后仍然显著时标记 `*`

### 正交多项式时间趋势解读

- `liner`（一次方/线性趋势）：p < 0.05 表示测量值随时间呈线性变化
- `quadratic`（二次方）：p < 0.05 表示时间趋势呈U型或倒U型
- `cubic`（三次方）：p < 0.05 表示趋势有更复杂的弯曲形态
- `biquadrate`（四次方）：一般较少关注，高阶项显著意味着趋势形态非常复杂

### 论文报告示例

> 采用重复测量方差分析比较三种诱导方法在不同时点（T0~T4）的血压变化。球对称检验结果显示，Mauchly W = 0.293, p = 0.178，符合球对称假设。结果表明，时间主效应具有统计学意义（F(4,48) = 106.6, p < 0.001），组别×时间的交互效应也具有统计学意义（F(8,48) = 19.1, p < 0.001）。事后比较显示……

## 常见问题与注意事项

**Q1: 重复测量方差分析与普通方差分析的区别是什么？**

A: 重复测量方差分析的误差项分为两层（受试者间 + 受试者内），而普通方差分析只有一个误差项。`aov()` 中的 `Error(ID/time)` 即指定了误差分层结构——ID为受试者间误差，ID:time为受试者内误差。如果用普通方差分析处理重复测量数据，会错误地将受试者内变异归入总残差，导致效应被低估或高估。

**Q2: 球对称假设不满足时怎么办？**

A: 查阅 `anova_test()` 输出的 Sphericity Corrections 部分：
- ε > 0.75 时报告 HFe 的校正结果
- ε < 0.75 时报告 GGe 的校正结果（更保守）
- 若校正后 p 值仍 < 0.05，结论不变
- 若校正前后 p 值结论反转（一个显著一个不显著），应报告校正结果
- 球对称严重不满足时考虑改用混合效应模型（`lme4::lmer()`）

**Q3: 数据必须是长格式吗？**

A: 是的。`aov()` 和 `anova_test()` 都需要长格式（每行一次测量）。原始数据通常是宽格式（每个时间点一列），用 `tidyr::pivot_longer()` 或 `reshape2::melt()` 转换。

**Q4: 多重比较时如何选择？**

A: 分三种情况：
1. 组间比较（哪个处理组更好）→ 与其他方差分析相同，LSD、Tukey、Dunnett等方法均可
2. 时间趋势探索（随时间如何变化）→ 用正交多项式分析各次方趋势
3. 时间点两两比较（哪个时间点出现显著变化）→ 用配对t检验，推荐以基线为参考组

**Q5: rstatix::anova_test() 和 aov() 如何选择？**

A: 建议优先使用 `rstatix::anova_test()`：
- 自动输出 Mauchly 球对称检验（`aov()` 不输出）
- 自动输出 GG/HF 校正（`aov()` 不输出）
- 输出广义eta方（ges）作为效应量指标
- 输出格式整洁，便于提取与报告

**Q6: SPSS中重复测量方差分析与R结果不一致怎么办？**

A: SPSS默认输出球对称检验及校正，且多重比较的校正方法（如Bonferroni、Sidak）与R默认方法有时不同，两者结果可能略有差异。确保关键结论（交互效应和主效应的显著/不显著）一致即可。R中需手动选择校正方法（如 `p.adjust.method = "bonferroni"`），与SPSS对齐。
