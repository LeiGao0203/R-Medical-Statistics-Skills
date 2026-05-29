---
name: medical-stat-ps-weighting
description: "R语言医学统计：倾向性评分逆概率加权（IPTW）。通过加权创建伪总体，使两组协变量分布均衡，包括ATE和ATT权重。TRIGGER when user mentions IPTW、逆概率加权、倾向性评分加权、stabilized weight、ATE、ATT、overlap weighting、重叠加权，or asks about weighting methods for causal inference. SKIP for PSM匹配、倾向性评分分层、倾向性评分回归。"
---

# 倾向性评分：加权 (Propensity Score Weighting)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用倾向性评分加权的典型场景：**

- 观察性研究中需要控制多个混杂因素，但不想像匹配那样损失样本量
- 需要估计总体平均处理效应（ATE），即推广到全人群的因果效应
- 处理组和对照组协变量分布严重不均衡，需通过权重构造伪总体使其均衡
- SMD 或 p 值显示基线不均衡的多协变量数据集
- 研究目标人群是两组协变量分布重叠部分（ATO），使用重叠加权

**不使用的场景：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 需要直接一对一匹配研究对象 | 倾向性评分匹配（PSM） |
| 希望将 PS 作为协变量放入回归模型 | 倾向性评分回归 |
| 希望按 PS 五分位数分层分析 | 倾向性评分分层 |
| 极端 PS 值较多（接近 0 或 1），权重方差过大 | 考虑使用重叠加权或截断权重（trimming） |
| 治疗因素有多水平（≥3 组） | 多分类倾向性评分加权（`PSweight` 包支持） |

**医学研究常见应用：**

- 基于电子健康档案或注册数据库评估药物/手术的因果效应
- 并存多个协变量不均衡的回顾性队列研究
- 敏感性分析中替代匹配或回归方法验证结果稳健性

## 前置条件

**R 包安装：**

```r
# 核心包
install.packages(c("twang", "tableone", "cobalt", "survey"))

# 高级/便捷包
install.packages("PSweight")
```

**核心数据要求：**

- `treatment`：二分类处理因素变量（0/1，factor 类型）
- `outcome`：结局变量（连续型使用线性回归/gaussian family，二分类使用 logistic 回归/binomial family）
- `covariates`：需要均衡的混杂变量（连续型或分类型均可）
- 数据不允许含有缺失值（缺失需先进行插补或删除）

**科学假设：**

- **无不可测量混杂（no unmeasured confounding）**：所有影响处理和结局的混杂因素都被测量并纳入 PS 模型
- **正性假设（positivity）**：每个个体接受不同处理水平的概率均 > 0，即 PS 值不能为 0 或 1
- **正确指定模型**：PS 模型（通常为 logistic 回归）正确地包含和指定了与处理和结局相关的协变量
- **一致性（consistency）**：观测到的处理对每个个体潜在结局的映射是明确的

## 方法选择决策树

```
你的分析目标 →
├── 估计总体平均处理效应（ATE），推广到全部研究对象
│   ├── 手动实现 → 权重 = 1/PS（处理组）、1/(1-PS)（对照组）
│   ├── 使用 R 包 → twang::ps()、PSweight::PSweight(weight="IPW")
│   └── 希望权重方差小，避免极端值影响 → 稳定化权重（stabilized weight）
│
├── 估计处理组的平均处理效应（ATT，即干预组的处理效应）
│   ├── 手动实现 → 权重 = 1（处理组）、PS/(1-PS)（对照组）
│   └── 使用 PSweight 包 → weight="treated"
│
├── 目标人群为两组协变量分布重叠部分（ATO）
│   ├── 特点：两组权重之和最小、伪总体方差最小、SMD 几乎为 0
│   └── 使用 PSweight 包 → weight="overlap"，权重 = 1-PS（处理组）、PS（对照组）
│
├── 需要处理极端权重（某些个体权重过大，导致估计不稳定）
│   ├── 在 1st 和 99th 百分位数截断 → 手动 clip 权重
│   └── 改用重叠加权（overlap weighting）替代 IPTW
│
└── 加权后进行效应估计
    ├── 连续型结局 → survey::svyglm()（Gaussian family）
    ├── 二分类结局 → survey::svyglm()（Binomial family）
    └── 生存结局 → survey::svycoxph()
```

## 标准工作流

### 步骤 1：数据准备与探索

```r
library(twang)
library(tableone)
library(cobalt)

data(lindner, package = "twang")

# 将分类变量转换为 factor 类型
lindner[, c(3, 4, 6, 7, 8, 10)] <- lapply(lindner[, c(3, 4, 6, 7, 8, 10)], factor)
str(lindner)
## 'data.frame':  996 obs. of  11 variables:
##  $ lifepres       : num
##  $ cardbill       : int
##  $ abcix          : Factor w/ 2 levels "0","1"
##  $ stent          : Factor w/ 2 levels "0","1"
##  $ height         : int
##  $ female         : Factor w/ 2 levels "0","1"
##  $ diabetic       : Factor w/ 2 levels "0","1"
##  $ acutemi        : Factor w/ 2 levels "0","1"
##  $ ejecfrac       : int
##  $ ves1proc       : Factor w/ 6 levels "0","1","2","3","4","5"
##  $ sixMonthSurvive: logi
```

其中 `abcix` 是处理因素变量，`sixMonthSurvive` 是二分类结局变量，`cardbill` 是连续型结局变量，其余变量是协变量。

### 步骤 2：加权前基线均衡性评估

```r
# 选择协变量列
covs <- colnames(lindner)[c(1, 4:10)]

# 使用 tableone 查看加权前数据情况
tab <- CreateTableOne(
  vars = covs,
  strata = "abcix",
  data = lindner
)
print(tab, showAllLevels = TRUE, smd = TRUE)
```

使用 `cobalt` 包更便捷地检查各协变量均衡性：

```r
covariates <- subset(lindner, select = c(1, 4:10))
bal.tab(covariates, treat = lindner$abcix, s.d.denom = "pooled",
        m.threshold = 0.1, un = TRUE, v.threshold = 2)
```

### 步骤 3：计算倾向性评分

倾向性评分是给定协变量条件下接受处理的概率，通常使用 logistic 回归估计：

```r
psfit <- glm(abcix ~ stent + height + female + diabetic + acutemi +
               ejecfrac + ves1proc,
             data = lindner, family = binomial())
ps <- psfit$fitted.values
```

也可使用随机森林、GBM 等方法估计 PS（twang 包支持）。

### 步骤 4：计算权重并分配给数据集

**ATE 权重（IPW / IPTW）：**

- 处理组权重 = 1 / PS
- 对照组权重 = 1 / (1 - PS)

```r
iptw <- ifelse(lindner$abcix == 1, 1/ps, 1/(1-ps))
lindner$iptw <- iptw
```

**ATT 权重：**

- 处理组权重 = 1
- 对照组权重 = PS / (1 - PS)

**重叠加权（overlap weighting）：**

- 处理组权重 = 1 - PS
- 对照组权重 = PS

### 步骤 5：加权后均衡性验证

```r
bal.tab(covariates, treat = lindner$abcix, s.d.denom = "pooled",
        weights = lindner$iptw,
        m.threshold = 0.1, un = TRUE, v.threshold = 2)
## Balance tally for mean differences
##                    count
## Balanced, <0.1        13
## Not Balanced, >0.1     0

## Effective sample sizes
##            Control Treated
## Unadjusted  298.    698.
## Adjusted    202.27  671.09
```

`Diff.Adj` 即加权后的 SMD，所有协变量 SMD < 0.1 表明均衡性良好。加权后有效样本量发生变化。

### 步骤 6：加权后分析

**使用 survey 包进行加权回归分析：**

```r
library(survey)

# 创建 survey 设计对象
df <- svydesign(ids = ~1, data = lindner, weights = ~ iptw)

# 加权后基线表
tab_IPTW <- svyCreateTableOne(vars = covs, strata = "abcix", data = df, test = TRUE)
print(tab_IPTW, showAllLevels = TRUE, smd = TRUE)

# 加权 logistic 回归（二分类结局）
f <- glm(sixMonthSurvive ~ abcix + stent + height + female + diabetic +
           acutemi + ejecfrac + ves1proc,
         data = lindner, family = binomial(),
         weights = iptw)
summary(f)
```

**使用 PSweight 包进行重叠加权（推荐方式）：**

```r
library(PSweight)

formula.ps <- abcix ~ stent + height + female + diabetic + acutemi + ejecfrac + ves1proc

PSweight <- PSweight(
  ps.formula = formula.ps,
  data = lindner,
  weight = "overlap",      # 重叠加权；可选 "IPW"、"treated"
  yname = "cardbill",      # 结局变量
  family = "gaussian",     # 连续型用 gaussian，二分类用 binomial
  ps.method = "glm",
  out.method = "glm"
)

summary(PSweight)
## Closed-form inference:
## Contrast:
##             0 1
## Contrast 1 -1 1
##             Estimate Std.Error       lwr    upr  Pr(>|z|)
## Contrast 1 1134.9079    1.3253 1132.3103 1137.5 < 2.2e-16 ***

# 均衡性检验
SumStat <- SumStat(ps.formula = formula.ps, data = lindner, weight = "overlap")
SumStat[["ess"]]   # 有效样本量
##   unweighted  overlap
## 0        298 287.4367
## 1        698 569.5826

summary(SumStat)   # 加权前后 SMD 对比
## unweighted result
##           Mean 0  Mean 1   SMD
## stent      0.584   0.705 0.254
## ...
## overlap result
##           Mean 0  Mean 1 SMD
## stent      0.633   0.633   0
## height   171.464 171.464   0
## ...

plot(SumStat)      # 均衡性图形展示
```

### 步骤 7：结果报告（论文写法）

在论文中报告加权分析结果时，应包含以下要素：

1. 倾向性评分模型的构建方法（如 logistic 回归）和纳入的协变量列表
2. 加权方法（ATE / ATT / ATO）和权重公式
3. 加权前后各协变量的 SMD 对比表（或 Love plot）
4. 加权后的有效样本量
5. 加权模型中处理因素的效应估计值（OR / HR / β）、95% CI 和 p 值
6. 作为敏感性分析，可报告未加权分析的结果以作对比

示例描述："采用基于 logistic 回归的倾向性评分逆概率加权（IPTW）控制组间协变量不均衡，加权后所有协变量标准化均数差（SMD）均 < 0.1，表明组间均衡性良好。加权 logistic 回归结果显示，abcix 与 6 个月生存率的关联具有统计学意义（OR = 5.96, 95% CI: 3.19–11.13, p < 0.001）。"

## 代码示例

```r
# ==========================================
# 完整 IPTW 分析流程
# ==========================================
library(twang)
library(tableone)
library(cobalt)
library(survey)

# -- 1. 数据加载与预处理 --
data(lindner, package = "twang")
lindner[, c(3, 4, 6, 7, 8, 10)] <- lapply(
  lindner[, c(3, 4, 6, 7, 8, 10)], factor)

# -- 2. 协变量定义 --
covs <- colnames(lindner)[c(1, 4:10)]
covariates <- subset(lindner, select = c(1, 4:10))

# -- 3. 加权前均衡性检查 --
tab <- CreateTableOne(vars = covs, strata = "abcix", data = lindner)
print(tab, showAllLevels = TRUE, smd = TRUE)
## 观察各协变量 p 值和 SMD 判断不均衡程度

bal.tab(covariates, treat = lindner$abcix, s.d.denom = "pooled",
        m.threshold = 0.1, un = TRUE, v.threshold = 2)
##     ves1proc_1     Binary -0.1894 Not Balanced, >0.1
##     ves1proc_2     Binary  0.1360 Not Balanced, >0.1

# -- 4. 计算倾向性评分 --
psfit <- glm(abcix ~ stent + height + female + diabetic + acutemi +
               ejecfrac + ves1proc,
             data = lindner, family = binomial())
ps <- psfit$fitted.values

# -- 5. 计算 IPTW 权重 --
iptw <- ifelse(lindner$abcix == 1, 1/ps, 1/(1-ps))
lindner$iptw <- iptw

# -- 6. 加权后均衡性验证 --
bal.tab(covariates, treat = lindner$abcix, s.d.denom = "pooled",
        weights = lindner$iptw,
        m.threshold = 0.1, un = TRUE, v.threshold = 2)
## Balance tally for mean differences
##     Balanced, <0.1        13
##     Not Balanced, >0.1     0

# -- 7. 加权后回归分析 --
f <- glm(sixMonthSurvive ~ abcix + stent + height + female +
           diabetic + acutemi + ejecfrac + ves1proc,
         data = lindner, family = binomial(), weights = iptw)
summary(f)
## abcix1       1.785e+00  3.215e-01   5.551 2.84e-08 ***
## exp(1.785) ≈ 5.96 → abcix 使 6 个月死亡风险降低（OR 约 5.96）

# ==========================================
# 重叠加权（PSweight 包）
# ==========================================
library(PSweight)

rm(list = ls())
data(lindner, package = "twang")

formula.ps <- abcix ~ stent + height + female + diabetic + acutemi +
  ejecfrac + ves1proc

PSweight <- PSweight(
  ps.formula = formula.ps,
  data = lindner,
  weight = "overlap",
  yname = "cardbill",
  family = "gaussian",
  ps.method = "glm",
  out.method = "glm"
)
summary(PSweight)
##             Estimate Std.Error       lwr    upr  Pr(>|z|)
## Contrast 1 1134.9079    1.3253 1132.3103 1137.5 < 2.2e-16 ***

SumStat <- SumStat(ps.formula = formula.ps, data = lindner, weight = "overlap")
summary(SumStat)
## overlap result
##           Mean 0  Mean 1 SMD
## stent      0.633   0.633   0
## ...所有协变量 SMD = 0
```

## 结果解读指南

**倾向性评分（PS）的实质：**
PS 是给定协变量条件下接受处理的概率（0–1）。并非疗效指标，而是均衡工具。

**权重含义：**
- ATE 权重下，PS 极高（≈1）的对照对象或 PS 极低（≈0）的处理对象会获得极大权重。这是 IPTW 的固有局限。
- 重叠加权自动对分布重叠区域的个体给予更高权重，两端极端值权重趋近于 0，因此 SMD 几乎完美均衡且方差更小。

**SMD（标准化均数差）解读：**
- 通常 SMD < 0.1 视为组间均衡
- 加权后 SMD 应较未加权 SMD 显著减小
- 若加权后仍有 SMD ≥ 0.1，需重新考虑 PS 模型设定（交互项、非线性项、替代估计方法）

**有效样本量（ESS）：**
- 加权后 ESS 小于原始样本量，反映权重分布不均带来的信息损失
- 重叠加权的 ESS 损失通常小于标准 IPTW

**处理效应估计：**
- IPTW logistic 回归中，abcix 的 OR = exp(1.785) ≈ 5.96，表明使用 abcix 的患者 6 个月存活率更高
- 需报告 OR 的 95% CI 和 p 值
- `PSweight::summary()` 直接输出效应估计、标准误、置信区间和 p 值

**R 自带 glm 中 weights 参数的注意事项：**
R 的 `glm()` 和 `lm()` 的 `weights` 参数并不等于样本的抽样权重（sampling weight）。建议使用 `survey::svyglm()` 进行加权回归，以获得正确的标准误估计。

## 常见问题与注意事项

**Q1: IPTW 和 ATT 的权重公式各自适用什么场景？**

| 权重类型 | 目标人群 | 处理组权重 | 对照组权重 | 适用场景 |
|----------|---------|-----------|-----------|---------|
| ATE | 所有研究对象 | 1/PS | 1/(1-PS) | 希望推广到全人群 |
| ATT | 处理组对象 | 1 | PS/(1-PS) | 仅关心已接受治疗者的效应 |
| ATO (overlap) | 重叠人群 | 1-PS | PS | 希望减少极端值影响 |

**Q2: 稳定化权重（stabilized weight）是什么？**

稳定化 IPTW = 处理率/PS（处理组）或 (1-处理率)/(1-PS)（对照组），其中处理率为样本总体的处理比例。稳定化权重可减少极端值导致的权重方差膨胀。

**Q3: 权重过大（>10）怎么办？**

方案一：改用重叠加权（overlap weighting）
方案二：截断权重（weight trimming），如将权重限制在 1st–99th 百分位
方案三：检查 PS 模型设定，纳入更多重要协变量或协变量交互项

**Q4: PSweight 包的 `weight` 参数有哪些选项？**

- `"IPW"`：逆概率加权（ATE）
- `"treated"`：仅针对处理组的加权（ATT）
- `"overlap"`：重叠加权（ATO）
- `"matching"`：与匹配对应的权重
- `"entropy"`：熵权重

**Q5: 加权后如何选择回归模型？**

- 连续型结局 → `family = "gaussian"`（线性模型）
- 二分类结局 → `family = "binomial"`（logistic 模型）
- 推荐使用 `survey::svyglm()` 而非 `glm(weights=...)`，因为后者不能正确处理加权数据的标准误

**Q6: PS 的估计方法有哪些选择？**

`PSweight` 的 `ps.method` 参数支持：
- `"glm"`：logistic 回归（最简单、最常用）
- `"gbm"`：广义提升模型（自动选择非线性关系）
- `"superlearner"`：超级学习者（集成多种方法）

**Q7: twang 包 vs PSweight 包 vs PSW 包如何选？**

- `PSweight`：功能最全面，支持多分类处理、多种权重类型、均衡性检验一站式
- `twang`：擅长使用 GBM 估计 PS，适合复杂非线性关系
- `PSW`：更轻量，侧重于二分类的加权分析和均衡性评估

**Q8: 加权分析的优缺点总结？**

优点：
- 保留了全部样本，不丢失信息
- 可以估计 ATE 和 ATT 不同的因果效应
- 适用于大多数回归模型（logistic、Cox、线性）
- 重叠加权对极端值不敏感，均衡性极佳

缺点：
- 极端 PS 值会导致权重方差巨大
- 依赖 PS 模型正确设定
- R 基础函数对加权数据处理不够完善，需借助 survey 包
- 对无不可测量混杂的假设敏感性高

---

**参考文献：**

涂博祥, 秦婴逸, 吴骋, 等. 倾向性评分加权方法介绍及R软件实现[J]. 中国循证医学杂志, 2022, 22(3): 365–372.
