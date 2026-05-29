---
name: medical-stat-multi-anova
description: "R语言医学统计：多因素方差分析。涵盖析因设计、正交设计、嵌套设计、裂区设计的方差分析，分析多个因素及其交互作用对连续型因变量的影响。TRIGGER when user mentions 多因素方差分析、析因设计、交互作用、正交设计、嵌套设计、裂区设计、主效应，or asks about analyzing two or more factors simultaneously using ANOVA. SKIP for 单因素方差分析（单因素ANOVA）、重复测量设计（repeated measures）、协方差分析（ANCOVA）。"
---

# 多因素方差分析 (Multi-factor ANOVA)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**什么时候用：**

- 研究中有两个或两个以上分类自变量（因素），每个因素有多个水平，因变量为连续型变量
- 需要同时分析多个因素的**主效应**（main effect）以及因素间的**交互作用**（interaction effect）
- 实验设计包括：析因设计（完全交叉设计）、正交设计、嵌套设计、裂区设计
- 典型医学场景：同时比较不同药物和不同剂量对疗效的影响、不同手术方式和不同术后时间对愈合效果的影响、不同培养条件下多个菌株的生长差异

**什么时候不用：**

- 只有一个因素 → 使用单因素方差分析（one-way ANOVA）
- 数据不符合正态性和方差齐性 → 使用非参数检验（如 Friedman 检验）或数据变换
- 同一受试对象在不同时间点重复测量 → 使用重复测量方差分析（repeated measures ANOVA），需要指定 `Error()` 误差项
- 需要控制一个或多个连续型协变量 → 使用协方差分析（ANCOVA）
- 因变量是分类变量 → 使用卡方检验或 logistic 回归

## 前置条件

**R包：**

```r
# 基础安装已包含 aov()，无需额外安装用于基础方差分析
# 可视化可选择性安装：
install.packages("gplots")  # plotmeans()
install.packages("HH")      # interaction2wt()
```

**数据格式要求：**

数据必须为长格式（long format），每行一个观测值：
- 因素列（自变量）：应为 `factor` 类型
- 因变量列：连续型数值型 `numeric`

**统计假设（针对每组组合）：**

1. **独立性**：各观测值相互独立
2. **正态性**：各单元格残差近似服从正态分布。可使用 `shapiro.test(residuals(model))` 检验
3. **方差齐性**：各单元格方差相等。可使用 `bartlett.test(y ~ interaction(x1, x2), data = df)` 检验
4. **无显著离群值**：检查箱线图或标准化残差

## 方法选择决策树

```
你的实验设计和研究目的 →
├── 多个因素，所有因素水平组合均有观测，想分析所有主效应和交互作用
│   ├── 2个因素 → 两因素析因设计方差分析（14.2 节）
│   ├── 3个或更多因素 → 三因素析因设计方差分析（14.3 节）
│   └── 公式: aov(y ~ A * B * C, data = df)
├── 因素较多但只关心部分主效应和交互，采用正交表安排实验
│   └── 正交设计方差分析（14.4 节）
│       公式: aov(y ~ A + B + C + D + A:B, data = df)
├── 因素存在层级关系，次级因素嵌套在上级因素之内
│   └── 嵌套设计方差分析（14.5 节）
│       公式: aov(y ~ factor1 / factor2, data = df)
│       等价: aov(y ~ factor1 + factor1:factor2, data = df)
└── 实验单位分两级，一级因素只作用于一级单位，二级因素作用于二级单位
    └── 裂区设计方差分析（14.6 节）
        公式: aov(y ~ factorA * factorB + Error(id/factorB), data = df)
```

## 标准工作流

### 步骤1：数据准备与探索

将数据整理为长格式，确保因素变量已转换为 `factor` 类型。使用 `str()` 检查数据结构，使用 `interaction.plot()` 或 `plotmeans()` 初步探索因素间关系。

### 步骤2：前提条件检验

- 使用 `shapiro.test(residuals(model))` 检验正态性
- 使用 `bartlett.test(y ~ interaction(A, B), data = df)` 检验方差齐性
- 使用箱线图检查离群值：`boxplot(y ~ A * B, data = df)`

### 步骤3：执行统计分析

根据研究设计选择对应公式调用 `aov()`，对结果调用 `summary()` 获取方差分析表。模型公式中：
- `A * B` 展开为 `A + B + A:B`
- `A / B` 展开为 `A + A:B`
- 因素间**交叉**（crossed）用 `*`，因素间**嵌套**（nested）用 `/`

### 步骤4：结果解读

查看 `summary()` 输出的方差分析表。重点关注：
- 各因素的 F 值和 P 值
- 交互作用项是否显著（若显著，主效应的解释需谨慎）
- 残差的自由度、平方和和均方误差

### 步骤5：结果报告

在论文中报告：各因素的自由度（df）、均方（MS）、F 值、P 值。当交互作用显著时，需进一步做简单效应分析或两两比较。

## 代码示例

### 两因素析因设计（课本例11-1）

```r
# 2×2 factorial: suture method × time → axonal passage rate
df11_1 <- data.frame(
  x1 = rep(c("外膜缝合","束膜缝合"), each = 10),
  x2 = rep(c("缝合1个月","缝合2个月"), each = 5),
  y = c(10,10,40,50,10,30,30,70,60,30,10,20,30,50,30,50,50,70,60,30)
)
str(df11_1)

f1 <- aov(y ~ x1 * x2, data = df11_1)
summary(f1)
##             Df Sum Sq Mean Sq F value Pr(>F)  
## x1           1    180     180   0.600 0.4499  
## x2           1   2420    2420   8.067 0.0118 *
## x1:x2        1     20      20   0.067 0.7995  
## Residuals   16   4800     300

# 可视化交互效应
interaction.plot(df11_1$x2, df11_1$x1, df11_1$y,
                 type = "b", col = c("red","blue"), pch = c(12,15),
                 xlab = "缝合时间", ylab = "轴突通过率")
```

### I×J 析因设计（课本例11-2）

```r
# 3×3 factorial: drug A dose × drug B dose → analgesia time
df11_2 <- data.frame(
  druga = rep(c("1mg","2.5mg","5mg"), each = 3),
  drugb = rep(c("5微克","15微克","30微克"), each = 9),
  y = c(105,80,65,75,115,80,85,120,125,115,105,80,125,130,90,65,
        120,100,75,95,85,135,120,150,180,190,160)
)

f2 <- aov(y ~ druga * drugb, data = df11_2)
summary(f2)
##             Df Sum Sq Mean Sq F value  Pr(>F)   
## druga        2   6572    3286   8.470 0.00256 **
## drugb        2   7022    3511   9.050 0.00190 **
## druga:drugb  4   7872    1968   5.073 0.00647 **
## Residuals   18   6983     388
```

### 三因素析因设计（课本例11-3）

```r
# 5×2×2 factorial: uniform type × humidity × activity → thermal sensation
df11_3 <- foreign::read.spss("datasets/例11-03-5种军装热感觉5-2-2.sav",
                             to.data.frame = T, reencode = "UTF-8")
df11_3$a <- factor(df11_3$a)

f3 <- aov(x ~ a * b * c, data = df11_3)
summary(f3)
##             Df Sum Sq Mean Sq F value   Pr(>F)    
## a            4   5.20    1.30   3.024   0.0224 *  
## b            1   9.94    9.94  23.138 6.98e-06 ***
## c            1 283.35  283.35 659.485  < 2e-16 ***
## a:b          4   1.94    0.48   1.128   0.3491    
## a:c          4   1.48    0.37   0.862   0.4905    
## b:c          1  12.68   12.68  29.514 5.82e-07 ***
## a:b:c        4   1.61    0.40   0.937   0.4472    
## Residuals   80  34.37    0.43
```

### 正交设计（课本例11-4）

```r
# Orthogonal design: only main effects + A:B interaction
df11_4 <- data.frame(
  a = rep(c("5度","25度"), each = 4),
  b = rep(c(0.5, 5.0), each = 2),
  c = c(10, 30),
  d = c(6.0, 8.0, 8.0, 6.0, 8.0, 6.0, 6.0, 8.0),
  x = c(86,95,91,94,91,96,83,88)
)
df11_4$a <- factor(df11_4$a)
df11_4$b <- factor(df11_4$b)
df11_4$c <- factor(df11_4$c)
df11_4$d <- factor(df11_4$d)

f4 <- aov(x ~ a + b + c + d + a:b, data = df11_4)
summary(f4)
##             Df Sum Sq Mean Sq F value Pr(>F)  
## a            1    8.0     8.0     3.2 0.2155  
## b            1   18.0    18.0     7.2 0.1153  
## c            1   60.5    60.5    24.2 0.0389 *
## d            1    4.5     4.5     1.8 0.3118  
## a:b          1   50.0    50.0    20.0 0.0465 *
## Residuals    2    5.0     2.5
```

### 嵌套设计（课本例11-6）

```r
# Nested design: temperature nested within catalyst type
df11_6 <- data.frame(
  factor1 = factor(rep(c("A","B","C"), each = 6)),
  factor2 = factor(rep(c(70,80,90,55,65,75,90,95,100), each = 2)),
  y = c(82,84,91,88,85,83,65,61,62,59,56,60,71,67,75,78,85,89)
)

# "/" in formula means factor2 is nested in factor1
f <- aov(y ~ factor1 / factor2, data = df11_6)
# 等价写法: aov(y ~ factor1 + factor1:factor2, data = df11_6)
summary(f)
##                 Df Sum Sq Mean Sq F value   Pr(>F)    
## factor1          2 1956.0   978.0  177.82 5.83e-08 ***
## factor1:factor2  6  401.0    66.8   12.15 0.000716 ***
## Residuals        9   49.5     5.5
```

### 裂区设计（课本例11-7）

```r
# Split-plot: drug (factorA) on rabbit level, toxin (factorB) on injection site
df11_7 <- data.frame(
  factorA = factor(rep(c("a1","a2"), each = 10)),
  factorB = factor(rep(c("b1","b2"), 10)),
  id = factor(rep(c(1:10), each = 2)),
  y = c(15.75,19.00,15.50,20.75,15.50,18.50,17.00,20.50,16.50,20.00,
        18.25,22.25,18.50,21.50,19.75,23.50,21.50,24.75,20.75,23.75)
)

# Error(id/factorB): factorB nested within each rabbit (id)
f <- aov(y ~ factorA * factorB + Error(id/factorB), data = df11_7)
summary(f)
## Error: id
##           Df Sum Sq Mean Sq F value   Pr(>F)    
## factorA    1  63.01   63.01   28.01 0.000735 ***
## Residuals  8  18.00    2.25
##
## Error: id:factorB
##                 Df Sum Sq Mean Sq F value   Pr(>F)    
## factorB          1  63.01   63.01  252.05 2.48e-07 ***
## factorA:factorB  1   0.11    0.11    0.45    0.521    
## Residuals        8   2.00    0.25
```

## 结果解读指南

**方差分析表（ANOVA table）各列含义：**

- `Df`：自由度。因素自由度 = 水平数 - 1，残差自由度 = 总观测数 - 所有参数个数
- `Sum Sq`：离均差平方和，表示该因素能解释的变异大小
- `Mean Sq`：均方 = Sum Sq / Df，即单位自由度的平均变异
- `F value`：F 统计量 = Mean Sq(因素) / Mean Sq(残差)。F 越大，该因素效应越可能真实存在
- `Pr(>F)`：P 值。通常以 0.05 为界

**解读优先级：**

1. 先看**交互作用**：若交互作用 P < 0.05，说明一个因素的效应依赖于另一个因素的水平，此时主效应的解释需谨慎，应考虑做简单效应分析
2. 再看**主效应**：若交互作用不显著，单独解读各因素的主效应
3. 最后关注**残差**：残差均方是计算 F 值的分母，过大可能说明模型拟合不佳或数据变异大

**论文报告示例：**

> 经析因设计方差分析，A 因素主效应具有统计学意义（F(1, 16) = 0.60, P = 0.4499），B 因素主效应具有统计学意义（F(1, 16) = 8.07, P = 0.0118），A 与 B 的交互作用无统计学意义（F(1, 16) = 0.067, P = 0.7995）。

## 常见问题与注意事项

**Q: `aov(y ~ A * B)` 和 `aov(y ~ A + B)` 有什么区别？**

`A * B` 包含 A 主效应、B 主效应、A:B 交互作用。`A + B` 只包含两个主效应，不含交互项。若研究关心交互作用，必须使用 `*`。

**Q: 嵌套设计的 `factor1 / factor2` 公式的含义？**

`factor1 / factor2` 在 R 公式中等价于 `factor1 + factor1:factor2`。表示 factor2 是嵌套在 factor1 之内的——factor2 的各个水平在不同 factor1 水平下没有对应关系。例如，催化剂 A、B、C 各自设定了不同温度进行实验，温度就嵌套在催化剂内。

**Q: 裂区设计中 `Error(id/factorB)` 的作用？**

裂区设计的误差项分为两部分：一级单位（id）的误差用于检验 A 因素主效应；二级单位（id:factorB）的误差用于检验 B 因素及其与 A 的交互作用。`Error(id/factorB)` 告诉 R 为模型指定这两层误差项。缺少 `Error()` 会错误地使用同一误差项检验所有效应，得到错误的 P 值。

**Q: R 结果和 SPSS 结果一致吗？**

对于普通析因设计，两者完全一致。对于裂区设计，SPSS 手工操作较繁琐且容易指定错误差项，而 R 的 `Error()` 参数使指定更直观。注意查看 R 输出的两个 `Error:` 分层结果——分别对应一级和二级单位层面的检验。

**Q: 交互作用显著后怎么办？**

交互作用显著意味着一个因素在不同水平下的效应不同。此时可以：
1. 固定一个因素水平，对另一个因素做简单效应分析
2. 使用 `emmeans` 包进行边际均值对比和事后检验
3. 用 `interaction.plot()` 绘图直观展示交互模式

**Q: 方差分析显著后需要两两比较吗？**

如果某个因素有三个及以上水平且主效应显著，通常需要事后两两比较（post-hoc test），使用 `TukeyHSD()` 函数或 `multcomp` 包。注意：当交互作用也显著时，简单效应分析的解读优先于主效应的两两比较。

**Q: 多因素方差分析和多因素回归分析的关系？**

实质上，多因素方差分析的 `aov()` 等价于线性模型 `lm()`，只是输出格式不同。`summary.aov()` 给出方差分析表，而 `summary.lm()` 给出回归系数。可根据报告需要选择使用。
