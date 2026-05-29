---
name: medical-stat-ps-matching
description: "R语言医学统计：倾向性评分匹配（PSM）。通过MatchIt包进行1:1或1:n匹配，平衡两组的基线特征，减少观察性研究中的选择偏倚。TRIGGER when user mentions 倾向性评分匹配、PSM、MatchIt、匹配后均衡性检验、SMD、标准化差异、love plot、cobalt，or asks about matching in observational studies. SKIP for 倾向性评分加权(IPTW)、倾向性评分分层、倾向性评分回归调整、倾向性评分DID."
---

# 倾向性评分：匹配 (Propensity Score Matching / PSM)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用 PSM 的典型场景：**

- 观察性研究（队列研究、病例对照研究、横断面研究）中，需要比较两组处理因素的效果，但基线不均衡
- 协变量较多时，倾向性评分将多个混杂因素的影响综合为一个值（PS 值），降低维度
- 处理因素为二分类变量（如用药 vs 未用药、手术 vs 保守治疗），结果变量可为连续型、分类型或生存数据
- 需要产生匹配后数据集，用于后续任何统计建模（t 检验、回归、生存分析等）

**不使用 PSM 的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 协变量很少（≤2 个），数据量充足 | 可直接用协方差分析或多因素回归调整 |
| 不希望损失样本量（PSM 会丢弃未匹配者） | 倾向性评分加权（IPTW），保留全部样本 |
| 处理因素为连续型变量或多分类变量 | 需特殊方法（`WeightIt` 包）或转换为二分类 |
| 仅需分层比较，不要求个体配对 | 倾向性评分分层（stratification） |
| 已有匹配后数据，想进一步校正残余混杂 | 倾向性评分回归调整（covariate adjustment） |

**PSM 的通用四步流程：**

1. 估计倾向性评分值（PS）
2. 利用 PS 值进行匹配
3. 均衡性检验及模型评价
4. 在匹配后数据上进行处理效应估计

## 前置条件

**R 包安装：**

```r
install.packages(c("MatchIt", "cobalt", "tableone", "ggplot2", "cowplot"))
```

如需使用机器学习方法估计 PS：

```r
install.packages(c("randomForest", "nnet", "rpart"))
```

**核心数据要求：**

- `treatment`（处理因素）：二分类变量（0/1 或 factor），代表是否接受处理
- `covariates`（协变量/混杂）：连续型或分类型变量，这些变量是两组间不平衡的来源
- 数据应为**完整数据**，匹配前需处理缺失值（删除、均值填补、多重插补等），`na.omit()` 是最快但不一定最优的方式
- 处理组（如干预组）一般应少于对照组（`matchit()` 默认以 0 组为对照，1 组为处理；可通过 `estimand` 参数改变）

**核心概念：**

- **倾向性评分（Propensity Score, PS）**：以处理因素为因变量、混杂因素为自变量，通过 logistic 回归（或其他模型）估计的每个研究对象接受处理的概率。即 `P(Treatment = 1 | covariates)`
- **SMD（Standardized Mean Difference，标准化均数差异）**：衡量协变量在组间差异的标准化指标，SMD < 0.1 通常被认为已均衡
- **VR（Variance Ratio，方差比）**：衡量连续变量在组间方差的比值，VR 越接近 1 越均衡，VR > 2.0 或 VR < 0.5 提示不均衡

## 方法选择决策树

```
你的需求 →
├── PS 值估计方法
│   ├── 传统医学统计，需要可解释性 → logistic 回归（`distance = "logit"`，默认）
│   ├── 协变量与处理因素关系非线性 → 在公式中增加二次项 `I(x^2)` 或交互项 `x1:x2`
│   ├── 数据量大、协变量维度高 → 广义可加模型（`distance = "GAMlogit"`）
│   ├── 追求预测精度，可接受黑箱模型 → 随机森林、神经网络（`distance = "nnet"`）或自定义 PS 值传入
│   └── 多个协变量、需扩展公式 → 先用 `glm()` 手动计算 PS，再传给 `distance` 参数
│
├── 匹配方法（`method` 参数）
│   ├── 最常用、简单 → 最近邻匹配（`method = "nearest"`，默认）
│   ├── 类别变量需精确控制 → 精确匹配（`exact = c("gender", "site")`），但样本损失大
│   ├── 协变量维数高、需在匹配中直接优化平衡 → 遗传匹配（`method = "genetic"`）
│   ├── 需要利用全部样本生成权重 → 全匹配（`method = "full"`），每个处理都有控制组的加权组合
│   └── 协变量相关性复杂 → 最优匹配（`method = "optimal"`）
│
├── 匹配策略
│   ├── 样本充足、保守方案 → 1:1 无放回（`ratio = 1`, `replace = FALSE`）
│   ├── 对照组远多于处理组、想保留信息 → 1:n 有放回（`ratio = 3`, `replace = TRUE`），一般不超过 1:4
│   ├── 需要提高匹配精度 → 设置卡钳值（`caliper = 0.2`），即 PS 差距 > 0.2 个标准差不配对
│   └── 丢弃共同支持域之外的个体 → `discard = "both"`（丢弃两组极端样本）
│
├── 均衡性判断标准
│   ├── 连续型协变量 → SMD < 0.1 且 VR ∈ [0.5, 2.0] → 均衡
│   ├── 分类型协变量 → SMD < 0.1 → 均衡
│   └── 也可辅以假设检验（t 检验 p > 0.05、卡方检验 p > 0.05），但学界更推荐 SMD
│
└── 均衡后 → 用匹配后数据（`match.data()` 提取）进行后续分析（t 检验、回归、生存分析等）
```

## 标准工作流

### 步骤 1：数据准备与缺失值处理

倾向性评分对缺失值敏感。匹配前必须处理缺失值。可选方案：
- 直接删除缺失行（`na.omit()`）——最简单但可能带来偏倚
- 多重插补后逐套匹配再合并
- KNN、随机森林等方法插补

务必在匹配前使用 `tableone` 查看原始数据的基线情况和 SMD，作为匹配前后对比的基线。

### 步骤 2：执行匹配

使用 `MatchIt::matchit()` 指定公式、数据、PS 估计方法和匹配策略。PS 估计结果通过 `$distance` 获取，匹配矩阵通过 `$match.matrix` 查看。

### 步骤 3：均衡性检验

使用 `cobalt::bal.tab()` 计算匹配前后的 SMD 和 VR，以 SMD < 0.1 和 VR 在 [0.5, 2.0] 区间内为均衡标准。`summary(m.out, standardize = TRUE)` 也会给出结果，但其对分类变量的 SMD 计算不够准确，推荐优先使用 `cobalt`。

### 步骤 4：提取匹配后数据

`match.data(m.out)` 返回匹配后的数据集，其中包含 `distance`（PS 值）、`weights`（匹配权重）和 `subclass`（匹配对编号）三列。后续分析直接使用此数据集。

### 步骤 5：在匹配后数据上分析处理效应

用匹配后数据执行你要做的统计分析——t 检验、卡方检验、线性回归、logistic 回归、Cox 回归等，直接比较两组结果即可。

### 步骤 6：可视化

- `cobalt::love.plot()`：绘制匹配前后各协变量 SMD 变化的 Love Plot
- `cobalt::bal.plot()`：绘制匹配前后各协变量分布的密度图、柱状图或 ECDF 图
- `plot(m.out, type = "hist")` / `plot(m.out, type = "jitter")`：MatchIt 自带图形，但较简陋

## 代码示例

### 1. 数据准备

```r
set.seed(2020)

# 构造模拟数据：吸烟(Smoke)对心血管疾病(CVD)的影响，年龄和性别为混杂
x.Gender <- rep(0:1, c(400, 600))
x.Age <- round(abs(rnorm(1000, mean = 45, sd = 15)))

z <- (x.Age - 45) / 15 - (x.Age - 45)^2 / 100 + 2 * x.Gender
tps <- exp(z) / (1 + exp(z))
Smoke <- as.numeric(runif(1000) < tps)
z.y <- x.Gender + 0.3 * x.Age + 5 * Smoke - 20
y <- exp(z.y) / (1 + exp(z.y))
CVD <- as.numeric(runif(1000) < y)

# 随机添加缺失值
x.Age.mask <- rbinom(1000, 1, 0.2)
x.Age <- ifelse(x.Age.mask == 1, NA, x.Age)

data <- data.frame(x.Age, x.Gender, Smoke, CVD)
head(data)
##   x.Age x.Gender Smoke CVD
## 1    51        0     1   0
## 2    50        0     0   0
## 3    29        0     0   0
## 4    28        0     0   0
## 5     3        0     0   0
## 6    56        0     1   1
```

### 2. 匹配前基线资料表

```r
library(tableone)

table2 <- CreateTableOne(
  vars = c("x.Age", "x.Gender", "CVD"),
  data = data,
  factorVars = c("x.Gender", "CVD"),
  strata = "Smoke",
  smd = TRUE
)
table2 <- print(table2, smd = TRUE, showAllLevels = TRUE,
                noSpaces = TRUE, printToggle = FALSE)
table2
##                    Stratified by Smoke
##                     level 0               1              p        test SMD    
##   n                 ""    "549"           "451"          ""       ""   ""     
##   x.Age (mean (SD)) ""    "42.76 (19.69)" "47.04 (8.14)" "<0.001" ""   "0.284"
##   x.Gender (%)      "0"   "299 (54.5)"    "101 (22.4)"   "<0.001" ""   "0.698"
##                     "1"   "250 (45.5)"    "350 (77.6)"   ""       ""   ""     
##   CVD (%)           "0"   "452 (82.3)"    "230 (51.0)"   "<0.001" ""   "0.705"
##                     "1"   "97 (17.7)"     "221 (49.0)"   ""       ""   ""
```

匹配前 `x.Age` 和 `x.Gender` 的 p 值均 < 0.001，SMD 分别为 0.284 和 0.698，表明两组基线严重不平衡。

### 3. 标准 PSM（logistic 回归 + 1:1 最近邻匹配）

```r
library(MatchIt)

data.complete <- na.omit(data)  # 先处理缺失

m.out <- matchit(Smoke ~ x.Age + x.Gender,
                 data = data.complete,
                 distance = "logit",
                 method = "nearest",
                 replace = FALSE,
                 ratio = 1)
m.out
## A matchit object
##  - method: 1:1 nearest neighbor matching without replacement
##  - distance: Propensity score
##              - estimated with logistic regression
##  - number of obs.: 831 (original), 738 (matched)
##  - target estimand: ATT
##  - covariates: x.Age, x.Gender
```

### 4. 查看匹配矩阵和丢弃情况

```r
head(m.out$match.matrix)
##    [,1] 
## 1  "204"
## 6  "283"
## 10 "56" 
## 12 "41" 
## 20 "79" 
## 26 "84"

table(m.out$discarded)
## FALSE 
##   831
```

### 5. 均衡性检验（cobalt 包）

```r
library(cobalt)

# SMD 阈值设为 0.1
bal.tab(m.out, m.threshold = 0.1, un = TRUE)
## Balance Measures
##              Type Diff.Un Diff.Adj        M.Threshold
## distance Distance  0.9050   0.5448                   
## x.Age     Contin.  0.5247  -0.0329     Balanced, <0.1
## x.Gender   Binary  0.3330   0.2222 Not Balanced, >0.1
##
## Balance tally for mean differences
##                    count
## Balanced, <0.1         1
## Not Balanced, >0.1     1
##
## Sample sizes
##           Control Treated
## All           462     369
## Matched       369     369
## Unmatched      93       0
```

匹配后 `x.Age` 的 SMD = -0.0329 < 0.1，已均衡；`x.Gender` 的 SMD = 0.2222 > 0.1，仍不均衡。

### 6. VR 检验

```r
bal.tab(m.out, v.threshold = 2)
## Balance Measures
##              Type Diff.Adj V.Ratio.Adj      V.Threshold
## distance Distance   0.5448      0.7605     Balanced, <2
## x.Age     Contin.  -0.0329      0.1908 Not Balanced, >2
## x.Gender   Binary   0.2222           .                 
##
## Balance tally for variance ratios
##                  count
## Balanced, <2         1
## Not Balanced, >2     1
```

`x.Age` 按 VR 标准也没有均衡（0.1908，远低于 1 且 < 0.5 提示方差不均）。

### 7. 提取匹配后数据并作统计检验

```r
mdata <- match.data(m.out)
head(mdata)
##   x.Age x.Gender Smoke CVD  distance weights subclass
## 1    51        0     1   0 0.2583040       1        1
## 2    50        0     0   0 0.2545807       1       16
## 6    56        0     1   1 0.2774451       1        2
## 9    71        0     0   1 0.3397803       1      237
## 10   47        0     1   1 0.2436248       1        3
## 12   59        0     1   1 0.2893402       1        4

# t 检验：x.Age 已无差异
t.test(x.Age ~ Smoke, data = mdata)
## t = 0.25327, df = 503.47, p-value = 0.8002

# 卡方检验：x.Gender 仍有差异
chisq.test(mdata$x.Gender, mdata$Smoke, correct = FALSE)
## X-squared = 41.342, df = 1, p-value = 1.278e-10
```

### 8. 换公式增加交互项 / 二次项

```r
m.out2 <- matchit(Smoke ~ I(x.Age^2) + x.Age + x.Gender,
                  data = data.complete)
```

### 9. 使用随机森林计算 PS

```r
library(randomForest)

data.complete$Smoke <- factor(data.complete$Smoke)
rf.out <- randomForest(Smoke ~ x.Age + x.Gender, data = data.complete)
eps <- rf.out$votes[, 2]  # 提取预测为 1 的概率

m.out3 <- matchit(formula = Smoke ~ x.Age + x.Gender,
                  data = data.complete,
                  distance = eps,
                  method = "nearest",
                  replace = TRUE,
                  discard = "both",
                  ratio = 2)
```

### 10. 精确匹配（处理性别不均衡）

```r
m.out4 <- matchit(Smoke ~ I(x.Age^2) + x.Age + x.Gender,
                  data = data.complete,
                  distance = "logit",
                  method = "nearest",
                  exact = c("x.Gender"),
                  replace = FALSE,
                  ratio = 1)
## Warning: Fewer control units than treated units in some `exact` strata;
## not all treated units will get a match.

bal.tab(m.out4, m.threshold = 0.1)
## Balance tally for mean differences
##                    count
## Balanced, <0.1         4
## Not Balanced, >0.1     0
##
## Sample sizes
##           Control Treated
## Matched       147     147
## Unmatched     315     222
```

精确匹配确实让所有协变量均衡了，但仅保留 294/831 = 35% 的样本。

### 11. 结果可视化（cobalt 包）

```r
library(cowplot)

plot_grid(
  bal.plot(m.out, var.name = "x.Age", which = "both", grid = TRUE),
  bal.plot(m.out, var.name = "x.Gender", which = "both", grid = TRUE),
  bal.plot(m.out, var.name = "x.Age", which = "both", grid = TRUE, type = "ecdf"),
  love.plot(bal.tab(m.out, m.threshold = 0.1),
            stat = "mean.diffs", grid = TRUE, stars = "raw", abs = FALSE)
)
```

- **密度图/柱状图**：展示匹配前后协变量在两组间的分布对比
- **ECDF 图**：累计分布函数对比，匹配后两条线越重合越均衡
- **Love Plot**：以竖线标出 SMD = ±0.1 阈值，匹配后的点落在竖线之间即为均衡

## 结果解读指南

### matchit() 输出解读

| 输出项 | 含义 |
|--------|------|
| `method: 1:1 nearest neighbor matching` | 匹配方法为 1:1 最近邻 |
| `distance: Propensity score` | PS 值由 logistic 回归估计 |
| `number of obs.: 831 (original), 738 (matched)` | 原始 831 人，成功匹配 738 人 |
| `target estimand: ATT` | 估计目标为处理组的平均处理效应（ATT），即仅在处理组中可找到匹配对照的个体 |

### bal.tab() 输出解读

| 指标 | 含义 | 均衡标准 |
|------|------|----------|
| `Diff.Un`（匹配前 SMD） | 匹配前组间标准化差异 | 越大越不平衡 |
| `Diff.Adj`（匹配后 SMD） | 匹配后组间标准化差异 | < 0.1 均衡 |
| `V.Ratio.Adj`（匹配后 VR） | 匹配后方差比 | [0.5, 2.0] 均衡 |
| `M.Threshold` / `V.Threshold` | 自动判断均衡与否 | — |

### 论文中如何报告

- "本研究采用倾向性评分匹配以控制组间混杂因素的不平衡。以吸烟状态为因变量，年龄和性别为自变量，通过 logistic 回归估计倾向性评分。采用 1:1 最近邻匹配，卡钳值设为 0.2。匹配后，各协变量的标准化均数差异（SMD）均 < 0.1，表明平衡良好。"
- 统计论文中通常会在表 1 报告匹配前后的 SMD 以及各协变量的组间比较结果（p 值）。

## 常见问题与注意事项

**Q1: 为什么我的匹配后 SMD 还 > 0.1？**

五种试法：① 换 PS 计算公式，增加二次项 `I(x^2)` 或交互项 `x1:x2`；② 换匹配方法（`method = "genetic"` 等）；③ 使用精确匹配（`exact = c("变量名")`），但会丢弃大量样本；④ 增加样本量（最根本的方法）；⑤ 匹配后再用回归或分层进一步调整。

**Q2: PS 就是 p 值吗？**

倾向性评分的确等价于 logistic 回归的预测概率值（p），但严格来说是参数 `π`（事件的概率），它是多协变量综合的结果。用 `glm(..., family = binomial())$fitted.values` 就可以自己计算，跟 `matchit()` 默认的 `distance = "logit"` 结果完全一样。

**Q3: cobalt 和 tableone 的 SMD 不一致怎么办？**

优先采用 `cobalt` 的结果。`summary(m.out)` 默认会把分类变量当连续变量计算 SMD，不够准确；`tableone` 也有 subtle issues。`cobalt::bal.tab()` 对分类变量（Binary）做了正确区分和处理，参考文献：Zhang Z et al., *Balance diagnostics after propensity score matching*, Ann Transl Med 2019;7:16。

**Q4: SMD 的阈值必须是 0.1 吗？**

0.1 是常用的参考阈值（约对应 Cohen's d 的小效应量），部分文献接受 0.2 甚至 0.25 作为粗平衡的标准。严格的临床试验倾向使用 0.1，探索性研究可适当放宽到 0.2。SPSS 中没有 SMD 的标准选项，需要手动计算或麻烦地二次编程。

**Q5: 匹配后怎么后续分析？**

用 `match.data()` 提取匹配后数据，然后在新的 `data.frame` 上做任何分析（t 检验、线性回归、logistic 回归、Cox 回归等）即可。注意：如果使用的是有放回匹配（`replace = TRUE`），部分对照个体会多次出现（被重复匹配），此时需要利用 `weights` 列进行加权分析。

**Q6: 缺失值如何处理？**

倾向性评分对缺失值很敏感。本章演示用 `na.omit()` 简单删除，但更好的做法包括：多重插补（`mice` 包）、KNN 插补、随机森林插补（`missForest` 包）。更高级的做法是在插补后的每套数据上分别做匹配再合并结果，但这在临床上并不总是强制执行。

**Q7: 倾向性评分匹配 vs. 倾向性评分加权 vs. 倾向性评分分层，怎么选？**

| 方法 | 优势 | 劣势 |
|------|------|------|
| 匹配（PSM） | 生成可直接分析的配对数据，结果直观 | 丢弃未匹配者，样本损失 |
| 加权（IPTW） | 保留全部样本，可推广到全人群 | 极端权重可能放大偏倚 |
| 分层 | 简单易理解 | 层数多时某些层样本量太小 |
| 回归调整 | 灵活可结合多种模型 | 依赖模型正确性 |

**Q8: 处理因素是多分类或连续型怎么办？**

`MatchIt` 主要面向二分类处理因素。多分类处理可使用 `WeightIt` 包的广义倾向性评分（generalized propensity score）功能。连续型处理也可用 `WeightIt` 或 `cbps` 等方法。

**Q9: 倾向性评分只能平衡测量到的协变量吗？**

是的。PSM 只能平衡纳入模型的、已测量的协变量，无法消除未测量混杂因素带来的偏倚。因此，研究设计和数据收集阶段对关键混杂的控制至关重要。可以通过敏感性分析（如 E-value）来评估未测量混杂的潜在影响。

---

**参考资料：**

1. Zhang Z, Kim HJ, Lonjon G, et al. Balance diagnostics after propensity score matching. *Ann Transl Med* 2019;7:16. DOI: 10.21037/atm-20-3998
2. https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html
3. https://ngreifer.github.io/cobalt/
