---
name: medical-stat-chisq
description: "R语言医学统计：卡方检验。涵盖独立样本四格表卡方检验（Pearson/连续性校正/Fisher确切概率法）、配对四格表McNemar检验、R×C列联表卡方检验、Cochran-Mantel-Haenszel分层检验及频数拟合优度检验。TRIGGER when user mentions 卡方检验、四格表、列联表、率或构成比的比较、分类变量关联性分析、McNemar检验、Fisher精确检验，or asks about analyzing categorical data. SKIP for 连续型变量比较（t检验/ANOVA）、ROC分析、回归建模。"
---

# 卡方检验 (Chi-squared Test)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用卡方检验：**

- 比较两个或多个样本率/构成比的差异（如：两种疗法的有效率比较、三组治疗方案的治愈率比较）
- 检验两个分类变量之间是否存在关联性（如：血型与MN血型系统之间的关联）
- 配对设计的分类资料比较（如：两种检测方法对同一样本的阳性/阴性检出率比较）
- 分层列联表分析——控制混杂因素后检验行变量与列变量的关联（CMH检验）
- 检验样本频数分布是否符合某一理论分布（拟合优度检验）

**何时不用卡方检验：**

| 情况 | 替代方法 |
|------|----------|
| 连续型变量的均值比较 | t检验 / 方差分析（ANOVA） |
| 单向有序R×C表（分组无序、结局有序） | 秩转换的非参数检验（Wilcoxon/Mann-Whitney） |
| 双向有序属性相同的R×C表（一致性评价） | Kappa检验（`vcd::Kappa()`） |
| 需要控制多个混杂因素或建模 | Logistic回归 / 对数线性模型 |
| 结局变量为连续型，需要分析趋势 | 线性回归 / Spearman等级相关 |

## 前置条件

**R 包依赖：**

```r
# 基础包（无需安装）
# stats    —— chisq.test(), fisher.test(), mcnemar.test(), mantelhaen.test(), prop.test()

# 可选安装
install.packages("gmodels")    # CrossTable() — 四格表综合输出，类似SPSS
install.packages("vcd")        # assocstats() — 列联系数 / Kappa() — 一致性检验
install.packages("DescTools")  # MHChisqTest() — 线性趋势 / BreslowDayTest() — 齐性检验 / CohenKappa()
install.packages("rcompanion") # pairwiseNominalIndependence() — 多重比较
install.packages("vcdExtra")   # CMHtest() — 两变量CMH检验
```

**统计假设与应用条件：**

1. **独立性**：各观测值相互独立（由研究设计保证）
2. **期望频数**：Pearson卡方要求理论频数 T ≥ 5 的格子比例不低于80%，且任意格子 T ≥ 1。不满足时需用校正公式或Fisher精确检验
3. **样本量**：四格表卡方的样本量 n ≥ 40 是使用Pearson卡方的基本条件之一
4. **数据类型**：分组变量和结局变量均为分类变量（二分类或多分类）

**数据格式要求：**

- **原始数据**（data.frame）：每行一个观测，包含分组变量和结局变量列
- **汇总数据**（matrix/table）：行×列列联表，由 `table()` 或 `matrix()` 构造
- **频数数据**（需转换）：三列数据（分组、结局、频数），可通过 `matrix()` 或 `xtabs()` 转换为列联表
- **配对四格表**：2×2 matrix，行和列分别代表两种方法/时间点的检测结果
- **CMH检验**：3维 array，第三维为分层变量，可通过 `array()` 或 `xtabs()` 构造

## 方法选择决策树

```
你的数据情况 →
├── 四格表（2×2）资料
│   ├── 独立样本设计（完全随机分组）
│   │   ├── n ≥ 40 且所有 T ≥ 5 → Pearson卡方检验（chisq.test(..., correct = FALSE)）
│   │   ├── n ≥ 40 但有 1 ≤ T < 5 → 连续性校正卡方（chisq.test(..., correct = TRUE)）
│   │   └── n < 40 或任意 T < 1 → Fisher确切概率法（fisher.test()）
│   │
│   └── 配对设计（自身前后、两种方法检测同一样本）
│       └── McNemar检验（mcnemar.test(..., correct = TRUE)）
│
├── R×C（行×列）列联表资料
│   ├── 双向无序（两变量均为无序分类）
│   │   ├── 目的：多个样本率/构成比比较 → Pearson卡方（chisq.test()）
│   │   └── 目的：分析关联性 → Pearson卡方 + 列联系数（vcd::assocstats()）
│   │
│   ├── 单向有序（分组无序，结局有序）
│   │   └── 秩转换的非参数检验（Wilcoxon/Mann-Whitney），不属于卡方检验
│   │
│   ├── 双向有序属性相同（如两种方法测同一等级指标）
│   │   └── Kappa一致性检验（vcd::Kappa() 或 DescTools::CohenKappa()）
│   │
│   └── 双向有序属性不同（如年龄组×疗效等级）
│       ├── 目的：检验组间疗效差异 → 秩转换的非参数检验
│       ├── 目的：检验线性趋势 → MHChisqTest() 或 Cochran-Armitage趋势检验
│       └── 目的：等级相关 → Spearman等级相关
│
├── 分层列联表（控制混杂因素）
│   └── Cochran-Mantel-Haenszel检验（mantelhaen.test()）
│       ├── 需检验层间效应值齐性 → BreslowDayTest() 或 woolf检验
│       └── X和Y均为有序 → 使用相关统计量（MH检验自动给出）
│
└── 频数分布拟合优度
    └── chisq.test(x, p = 理论概率向量)
```

## 标准工作流

### 步骤1：数据准备与探索

原始个体数据用 `table()` 生成列联表；汇总频数数据用 `matrix()` 手工构造列联表。用 `str()` 检查数据结构，确认分组变量和结局变量均为因子型。

```r
# 原始数据→列联表
mytable <- table(df$group, df$outcome)
mytable

# 或手工构造矩阵
my_mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                 dimnames = list(Group = c("A", "B"),
                                 Outcome = c("有效", "无效")))
```

### 步骤2：检查理论频数

对列联表计算期望频数，判断应采用Pearson卡方、校正卡方还是Fisher精确检验。

```r
chisq.test(mytable)$expected  # 查看期望频数
# 判断规则：n≥40且所有T≥5→Pearson; n≥40且1≤T<5→校正; n<40或T<1→Fisher
```

### 步骤3：执行统计分析

根据数据类型和研究目的选择对应方法执行检验。对R×C表需要注意区分双向无序、单向有序、双向有序等不同情况。

### 步骤4：结果解读

查看χ²值、自由度、p值。若p < 0.05，说明组间率/构成比差异有统计学意义。对多重比较，需校正检验水准 α' = α / (比较次数)。

### 步骤5：结果报告

论文中报告格式示例："两种疗法有效率比较，Pearson χ² = 12.857, df = 1, p = 0.0003，差异有统计学意义，治疗组有效率（95.2%）高于安慰剂组（78.1%）。"

## 代码示例

### 四格表卡方检验（课本例7-1）— Pearson卡方

```r
# 构造原始数据
ID <- seq(1, 200)
treat <- factor(c(rep("treated", 104), rep("placebo", 96)))
impro <- factor(c(rep("marked", 99), rep("none", 5),
                  rep("marked", 75), rep("none", 21)))
data1 <- data.frame(ID, treat, impro)

mytable <- table(data1$treat, data1$impro)
mytable
##           marked none
##   placebo     75   21
##   treated     99    5

# Pearson卡方（服从n≥40且T≥5条件）
chisq.test(mytable, correct = FALSE)
##
##  Pearson's Chi-squared test
##
## data:  mytable
## X-squared = 12.857, df = 1, p-value = 0.0003362
```

### 四格表卡方检验（课本例7-2）— 连续性校正

```r
per <- matrix(c(46, 6, 18, 8),
              nrow = 2, byrow = TRUE,
              dimnames = list(group = c("胞磷胆碱", "神经节苷脂"),
                              effect = c("有效", "无效")))
per
##             effect
## group        有效 无效
##   胞磷胆碱     46    6
##   神经节苷脂   18    8

# 连续性校正卡方（本例有1≤T<5，需校正）
chisq.test(per, correct = TRUE)
## Warning: Chi-squared approximation may be incorrect
##
##  Pearson's Chi-squared test with Yates' continuity correction
##
## data:  per
## X-squared = 3.1448, df = 1, p-value = 0.07617
```

### gmodels::CrossTable() — 一站式输出（类似SPSS）

```r
library(gmodels)
CrossTable(data1$treat, data1$impro, digits = 4,
           expected = TRUE, chisq = TRUE, fisher = TRUE, mcnemar = TRUE,
           format = "SPSS")
# 输出同时包含：四格表（含期望频数）、Pearson卡方、Yates校正、Fisher精确检验、
# McNemar检验，一次性提供所有常用选项
```

### 配对四格表 McNemar 检验（课本例7-3）

```r
ana <- matrix(c(11, 12, 2, 33), nrow = 2, byrow = TRUE,
              dimnames = list(免疫荧光 = c("阳性", "阴性"),
                              乳胶凝集 = c("阳性", "阴性")))
ana
##         乳胶凝集
## 免疫荧光 阳性 阴性
##     阳性   11   12
##     阴性    2   33

mcnemar.test(ana, correct = TRUE)
##
##  McNemar's Chi-squared test with continuity correction
##
## data:  ana
## McNemar's chi-squared = 5.7857, df = 1, p-value = 0.01616
```

### Fisher 确切概率法（课本例7-4）

```r
hbv <- matrix(c(4, 18, 5, 6), nrow = 2, byrow = TRUE,
              dimnames = list(组别 = c("预防注射组", "非预防组"),
                              效果 = c("阳性", "阴性")))
hbv
##             效果
## 组别       阳性 阴性
##   预防注射组    4   18
##   非预防组      5    6

fisher.test(hbv)
##
##  Fisher's Exact Test for Count Data
##
## data:  hbv
## p-value = 0.121
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  0.03974151 1.76726409
## sample estimates:
## odds ratio
##  0.2791061
```

### R×C列联表 — 多个样本率比较（课本例7-6）

```r
df <- read.csv("datasets/例07-06.csv", header = TRUE)
# 频数数据转换为矩阵
M <- matrix(df$f, nrow = 3, byrow = TRUE,
            dimnames = list(trt = c("物理", "药物", "外用"),
                            effect = c("有效", "无效")))
M
##       effect
## trt    有效 无效
##   物理  199    7
##   药物  164   18
##   外用  118   26

chisq.test(M, correct = FALSE)
##
##  Pearson's Chi-squared test
##
## data:  M
## X-squared = 21.038, df = 2, p-value = 2.702e-05

# prop.test() 也适用于两列率的比较
prop.test(M, correct = TRUE)
## X-squared = 21.038, df = 2, p-value = 2.702e-05
```

### R×C列联表 — 双向无序关联性检验（课本例7-8）+ 列联系数

```r
blood <- matrix(c(431, 490, 902, 388, 410, 800,
                  495, 587, 950, 137, 179, 32),
                nrow = 4, byrow = TRUE,
                dimnames = list(abo = c("o", "a", "b", "ab"),
                                mn = c("m", "n", "mn")))
chisq.test(blood, correct = FALSE)
## X-squared = 213.16, df = 6, p-value < 2.2e-16

library(vcd)
assocstats(blood)
##                     X^2 df P(> X^2)
## Likelihood Ratio 248.14  6        0
## Pearson          213.16  6        0
## Contingency Coeff.: 0.188
## Cramer's V        : 0.136
```

### R×C列联表 — 双向有序线性趋势检验（课本例7-9）

```r
ather <- matrix(c(70, 22, 4, 2, 27, 24, 9, 3,
                  16, 23, 13, 7, 9, 20, 15, 14),
                nrow = 4, byrow = TRUE,
                dimnames = list(age = c("20~", "30~", "40~", "≥50"),
                                level = c("-", "+", "++", "+++")))

# 总卡方
chisq.test(ather)
## X-squared = 71.432, df = 9, p-value = 7.97e-12

# 线性回归分量卡方
library(DescTools)
MHChisqTest(ather)
## X-squared = 63.389, df = 1, p-value = 1.696e-15

# 非线性分量卡方 = 71.432 - 63.389 = 8.043，df = 9 - 1 = 8
pchisq(q = 8.043, df = 8, lower.tail = FALSE)
## [1] 0.4292811
```

### 多个样本率间的多重比较（卡方分割法）

```r
# 两两比较，检验水准 α' = 0.05 / (K×(K-1)/2 + 1)
# K=3时 α' = 0.05 / (3×2/2+1) = 0.0125

chisq.test(M[1:2, ], correct = FALSE)       # 物理 vs 药物
## X-squared = 6.756, df = 1, p-value = 0.009343  # p < 0.0125，有统计学意义

chisq.test(M[c(1, 3), ], correct = FALSE)   # 物理 vs 外用
## X-squared = 21.323, df = 1, p-value = 3.881e-06 # p < 0.0125，有统计学意义

chisq.test(M[2:3, ], correct = FALSE)       # 药物 vs 外用
## X-squared = 4.591, df = 1, p-value = 0.03214   # p > 0.0125，无统计学意义

# 或使用 rcompanion 包一步完成
library(rcompanion)
pairwiseNominalIndependence(M)
##    Comparison  p.Chisq p.adj.Chisq
## 1 物理 : 药物 9.34e-03    1.40e-02
## 2 物理 : 外用 3.88e-06    1.16e-05
## 3 药物 : 外用 3.21e-02    3.21e-02
```

### Cochran-Mantel-Haenszel 分层检验（课本例7-12）

```r
# 构造三维数组：行=心肌梗死，列=口服避孕药，层=年龄分层
myo <- array(c(17, 47, 121, 944, 12, 158, 14, 663),
             dim = c(2, 2, 2),
             dimnames = list(心肌梗死 = c("病例", "对照"),
                             口服避孕药 = c("是", "否"),
                             年龄分层 = c("<40岁", "≥40岁")))

mantelhaen.test(myo, correct = FALSE)
##
##  Mantel-Haenszel chi-squared test without continuity correction
##
## data:  myo
## Mantel-Haenszel X-squared = 24.184, df = 1, p-value = 8.755e-07
## alternative hypothesis: true common odds ratio is not equal to 1
## 95 percent confidence interval:
##  1.930775 4.933840
## sample estimates:
## common odds ratio
##          3.086444

# Breslow-Day 层间效应值齐性检验
library(DescTools)
BreslowDayTest(myo)
## X-squared = 0.23409, df = 1, p-value = 0.6285

# Woolf法检验层间效应值齐性
woolf <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  1 - pchisq(sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2), k - 1)
}
woolf(myo)
## [1] 0.6400154
```

### 原始数据形式进行CMH检验

```r
myoo.tab <- xtabs(~口服避孕药 + 心肌梗死 + 年龄分层, data = myoo)
mantelhaen.test(myoo.tab, correct = FALSE)
# 结果与 array 形式相同

# 两变量CMH检验（vcdExtra::CMHtest）
vcdExtra::CMHtest(my_matrix, types = "cor")
```

### 频数分布拟合优度卡方检验（课本例7-13）

```r
x <- c(26, 51, 75, 63, 38, 17, 9)
p <- c(0.0854, 0.2102, 0.2585, 0.2120, 0.1304, 0.0641, 0.0394)

chisq.test(x = x, p = p)
##
##  Chi-squared test for given probabilities
##
## data:  x
## X-squared = 2.0377, df = 6, p-value = 0.9162
```

## 结果解读指南

**chisq.test() 输出各组成部分的解读：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `X-squared` | χ² 统计量 | 值越大，观测频数与期望频数差异越大 |
| `df` | 自由度 | 四格表：1；R×C表：(行数-1)×(列数-1) |
| `p-value` | p值 | p < 0.05 拒绝H₀（组间率/构成比差异有统计学意义） |

**fisher.test() 输出额外信息：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `odds ratio` | 优势比 (OR) | >1 表示暴露组事件风险更高，<1 表示保护因素 |
| `95% confidence interval` | OR的95% CI | 不包含1等价于 p < 0.05 |

**mantelhaen.test() 输出：**

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| `Mantel-Haenszel X-squared` | 控制分层变量后的关联检验 | p < 0.05 表示在控制混杂后仍有统计学关联 |
| `common odds ratio` | 公共OR（各层加权平均） | Breslow-Day p > 0.05 则该OR有效代表各层 |
| `95% confidence interval` | 公共OR的95% CI | 不包含1等价于统计学意义 |

**mcnemar.test() 解读：**

- H₀：两种方法的阳性率相同（b = c）
- 只基于不一致的对子（b和c）进行比较
- p < 0.05 表示两种方法的检出率差异有统计学意义

**四格表方法选择的临界值记忆口诀：**

```
n≥40且T≥5 → Pearson卡方，不校正
n≥40有T在1-5间 → Yates连续性校正 / Fisher都可以
n<40 或 T<1 → 只能用Fisher确切概率法
配对设计 → McNemar检验
```

**assocstats() 列联系数解读：**

- **Contingency Coefficient（列联系数）**：0 ~ √((k-1)/k)，k为行列数中较小值
- **Cramer's V**：0 ~ 1，值越大关联性越强
- 两者均表示关联强度，非因果关系

## 常见问题与注意事项

**Q1：R的 `chisq.test()` 默认 `correct = TRUE` 和 SPSS 默认行为不同怎么办？**

R的 `chisq.test()` 默认对四格表进行Yates连续性校正（`correct = TRUE`）。SPSS默认不做校正。若要与SPSS结果一致，设置 `correct = FALSE`。对于符合Pearson卡方条件（n≥40且T≥5）的数据，推荐 `correct = FALSE`。

**Q2：R×C表卡方检验有统计学意义后，如何做两两比较？**

使用卡方分割法：将多组拆分为多个四格表分别检验，但检验水准需校正为 α' = α / [K×(K-1)/2 + 1]，其中K为组数。也可使用 `rcompanion::pairwiseNominalIndependence()` 直接获得校正后p值。

**Q3：什么情况下卡方检验会给出警告"Chi-squared approximation may be incorrect"？**

当期望频数T < 5的格子数超过20%或有格子T < 1时，卡方近似不可靠。此时应：
- 四格表：改用Fisher精确检验
- R×C表：合并相邻的行/列，或使用Fisher精确检验的扩展

**Q4：Fisher精确检验和Pearson卡方检验结果不同时，以哪个为准？**

以Fisher精确检验为准。Fisher精确检验基于超几何分布直接计算精确概率，不依赖大样本近似。当两种方法结果不一致时，说明卡方近似的条件可能不满足，应优先使用Fisher结果。

**Q5：配对四格表为什么不能直接用Pearson卡方检验？**

配对四格表的行和列来自同一组受试者的两次测量，数据不独立。Pearson卡方的独立性假设不成立。配对设计应使用McNemar检验，它只分析不一致的观测对子（即b≠c的格子）。

**Q6：CMH检验和直接分层两两卡方检验的区别？**

- CMH检验给出控制分层变量后的综合关联检验和一个公共OR估计
- 直接分层检验：各层分别做卡方检验，无法给出整体结论
- CMH的前提是Breslow-Day检验p > 0.05（各层OR一致），若不满足需分报告各层结果

**Q7：`prop.test()` 和 `chisq.test()` 对于多个样本率比较有区别吗？**

对于两列（如有效/无效）的R×C表，两者结果完全等价，χ²值和p值相同。`prop.test()` 额外提供各组的样本率估计值。对于超过两列的多分类列联表，只能用 `chisq.test()`。

**常见错误提醒：**

1. 配对设计错误使用独立样本卡方——配对四格表必须用 McNemar 检验，否则结论完全错误
2. 单向有序R×C表（疗效等级数据）错误用卡方——应使用秩转换的非参数检验
3. R×C表卡方显著后不做多重比较直接下"各组间有差异"的笼统结论——应用卡方分割法明确哪些组间有差异
4. 忽略理论频数条件盲目使用Pearson卡方——需要先检查 `chisq.test(x)$expected` 再决定方法
5. 拟合优度检验的概率向量 `p` 之和不为1——`chisq.test()` 要求p各元素之和必须为1
