---
name: medical-stat-sample-size
description: "R语言医学统计：样本量计算。使用pwr包进行常见研究设计的样本量估计，涵盖t检验（单样本/两样本/配对）、方差分析、卡方检验、相关分析、比例比较等。TRIGGER when user mentions 样本量计算、样本量估计、power analysis、statistical power、pwr包，or asks about how many subjects needed for a study. SKIP for 事后功效分析（post-hoc）、统计检验本身。"
---

# 样本量计算 (Sample Size Calculation)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用样本量计算：**

- 研究设计阶段，确定需要多少受试对象以保证足够的统计功效
- 临床试验方案中估算样本量，满足伦理委员会和研究注册要求
- 基金申请书/毕业论文中陈述样本量依据
- 已知预期效应值、显著性水平和统计功效，反推所需样本量
- 已知样本量、显著性水平和统计功效，反推可检测到的最小效应值

**常见研究设计对应的pwr函数：**

| 研究设计 | pwr函数 | 效应值类型 |
|----------|---------|-----------|
| 单样本t检验（均数与已知总体比较） | `pwr.t.test(type="one.sample")` | d |
| 两独立样本t检验 | `pwr.t.test(type="two.sample")` | d |
| 配对t检验 | `pwr.t.test(type="paired")` | d |
| 多样本均数比较（单因素ANOVA） | `pwr.anova.test()` | f |
| 单样本率与已知总体率比较 | `pwr.p.test()` + `ES.h()` | h |
| 两独立样本率比较 | `pwr.2p.test()` + `ES.h()` | h |
| 多样本率比较（卡方检验） | `pwr.chisq.test()` + `ES.w2()` | w |
| 直线相关分析 | `pwr.r.test()` | r |

**何时不用样本量计算：**

| 情况 | 替代方法 |
|------|----------|
| 事后功效分析（post-hoc power） | 可算但意义有限；优先报告置信区间宽度 |
| 已完成研究，直接做统计检验 | t检验/方差分析/卡方检验等 |
| 复杂的存活分析、多因素回归设计 | 需要使用PASS软件或专门的公式 |
| 需要调整脱落率 | 在pwr结果上除以(1-脱落率)即可 |

## 前置条件

**R 包依赖：**

```r
install.packages("pwr")  # 功效分析核心包
```

**核心概念：功效分析的四个要素**

- **样本量（n）**：每组所需受试对象数（两样本t检验的n为每组的例数）
- **显著性水平（α / sig.level）**：一类错误概率，通常设为0.05
- **统计功效（power / 1-β）**：二类错误概率的补集，通常设为0.8或0.9
- **效应值（effect size）**：标准化的组间差异，是样本量计算的关鍵参数

计算样本量就是已知其中三个参数，求第四个参数的过程。不同研究设计的效应值类型不同。

**效应值 Cohen 约定标准：**

| 检验类型 | 效应值参数 | 小 | 中 | 大 |
|----------|-----------|----|----|-----|
| t检验 | d = (μ1 - μ2) / σ | 0.2 | 0.5 | 0.8 |
| 方差分析 | f | 0.1 | 0.25 | 0.4 |
| 卡方检验 | w | 0.1 | 0.3 | 0.5 |
| 比例比较 | h | 0.2 | 0.5 | 0.8 |
| 相关分析 | r | 0.1 | 0.3 | 0.5 |

## 方法选择决策树

```
你的研究设计 →
├── 均数比较
│   ├── 1组样本，与已知总体均数比较
│   │   └── pwr.t.test(type="one.sample", d = (样本均数-总体均数)/σ)
│   ├── 2组，配对设计（同体前后/左右）
│   │   └── pwr.t.test(type="paired", d = 差值均数/差值标准差)
│   ├── 2组，独立样本设计
│   │   └── pwr.t.test(type="two.sample", d = (μ1 - μ2) / σ)
│   └── ≥3组，独立样本设计（单因素ANOVA）
│       └── pwr.anova.test(k=组数, f=组间SD/组内SD)
│           ⚠ R计算f值困难，推荐使用PASS软件
│
├── 率比较
│   ├── 1个率，与已知总体率比较
│   │   └── pwr.p.test(h = ES.h(样本率, 总体率), alternative=...)
│   ├── 2个率，独立样本
│   │   └── pwr.2p.test(h = ES.h(率1, 率2), alternative=...)
│   │   └── 或 power.prop.test(p1=率1, p2=率2) 【base R】
│   └── ≥3个率（行×列表/卡方检验）
│       └── pwr.chisq.test(w = ES.w2(概率矩阵/组数), df = (行-1)×(列-1))
│
└── 相关分析
    └── 两连续变量直线相关
        └── pwr.r.test(r = 预期相关系数, alternative="two.sided")
```

## 标准工作流

### 步骤1：确定研究设计类型

明确比较类型（均数/率/相关）、组数、配对与否，确定使用哪个pwr函数。

### 步骤2：确定已知参数

- **sig.level**：通常0.05（双侧）或0.025（单侧）
- **power**：通常0.8（β=0.2）或0.9（β=0.1）
- **effect size**：从文献综述、预实验或Cohen约定标准估计

### 步骤3：计算效应值

```r
# t检验——直接算d
d <- (mu1 - mu2) / sigma

# 率比较——用ES.h()算h
h <- ES.h(p1, p2)

# 卡方检验——用ES.w2()算w
prob <- rbind(c(有效率), c(无效率))  # 2行×k列矩阵
w <- ES.w2(prob / k)

# 相关——直接用相关系数r
r <- 0.8
```

### 步骤4：执行样本量计算

```r
library(pwr)

# 示例：两样本t检验，α=0.05，power=0.9，d=0.6
result <- pwr.t.test(d = 0.6, sig.level = 0.05, power = 0.9,
                     type = "two.sample", alternative = "two.sided")
```

### 步骤5：结果解读与报告

查看输出中的 n。对于两样本设计，`n` 是 **每组** 所需例数。总样本量为 n×2。向上取整：n=59.35 → 每组至少60例。

论文报告格式示例："以α=0.05（双侧），β=0.1（把握度90%），预期两组间标准化均数差d=0.6，采用两独立样本t检验，经pwr包计算每组需60例，共需120例研究对象。"

## 代码示例

```r
library(pwr)

# ============================================================
# 10.1 t检验的样本量计算
# ============================================================

# --- 单样本t检验（课本例36-3）---
# 用药治疗矽肺患者，σ=25，α=0.05，β=0.1，辨别出尿矽排出量平均增加10mg/L
# d = 10/25 = 0.4, 单侧检验
pwr.t.test(d = 10/25, sig.level = 0.05, power = 1 - 0.1,
           type = "one.sample", alternative = "greater")
##      One-sample t test power calculation
##               n = 54.90553
##               d = 0.4
##       sig.level = 0.05
##           power = 0.9
##     alternative = greater

# base R 等价写法
power.t.test(delta = 10, sd = 25, sig.level = 0.05, power = 1 - 0.1,
             type = "one.sample", alternative = "one.sided")
##               n = 54.90553

# --- 两样本t检验（课本例36-4）---
# 比较A、B两种处理的平均血流量增加，(μ1-μ2)/σ = 0.6, α=0.05, β=0.1
pwr.t.test(d = 0.6, sig.level = 0.05, power = 1 - 0.1,
           type = "two.sample", alternative = "two.sided")
##      Two-sample t test power calculation
##               n = 59.35155
##               d = 0.6
##       sig.level = 0.05
##           power = 0.9
##     alternative = two.sided
## NOTE: n is number in *each* group

# ============================================================
# 10.2 多样本均数比较（ANOVA，课本例36-5）
# ============================================================
# ⚠ 效应量f计算困难，R中难以直接计算。推荐使用PASS软件。
# pwr.anova.test(k = 4, f = ???, sig.level = 0.05, power = 1 - 0.1)

# ============================================================
# 10.3 样本率和已知总体率的比较（课本例36-6）
# ============================================================
# 常规方法有效率80%, 新方法预期90%, α=0.05, β=0.1
ES.h(0.9, 0.8)
## [1] 0.2837941

pwr.p.test(h = ES.h(0.9, 0.8), sig.level = 0.05, power = 1 - 0.1,
           alternative = "greater")
##      proportion power calculation for binomial distribution (arcsine transformation)
##               h = 0.2837941
##               n = 106.3315
##       sig.level = 0.05
##           power = 0.9
##     alternative = greater
# 注意：与课本结果(137.1)有差异，因计算方法不同，建议用PASS或套课本公式

# ============================================================
# 10.4 两独立样本率的比较（课本例36-7）
# ============================================================
# 甲药有效率60%, 乙药85%, α=0.05, 1-β=0.9
ES.h(0.85, 0.60)
## [1] 0.5740396

pwr.2p.test(h = ES.h(0.85, 0.60), sig.level = 0.05, power = 1 - 0.1,
            alternative = "two.sided")
##      Difference of proportion power calculation for binomial distribution
##               h = 0.5740396
##               n = 63.77382
##       sig.level = 0.05
##           power = 0.9
##     alternative = two.sided
## NOTE: same sample sizes

# base R 等价写法
power.prop.test(p1 = 0.85, p2 = 0.6, sig.level = 0.05, power = 1 - 0.1,
                alternative = "two.sided")
##               n = 64.93465
##              p1 = 0.85
##              p2 = 0.6
## NOTE: n is number in *each* group

# ============================================================
# 10.5 多样本率的比较（卡方检验，课本例36-8）
# ============================================================
# 3种方法治疗消化性溃疡，甲40%, 乙50%, 丙65%, α=0.05, β=0.1
prob <- rbind(c(0.4, 0.5, 0.65),  # 有效率
              c(0.6, 0.5, 0.35))  # 无效率

ES.w2(prob / 3)  # 有几组就除以几
## [1] 0.2055947

pwr.chisq.test(w = ES.w2(prob / 3), df = 2,      # df = (3-1)*(2-1)
               sig.level = 0.05, power = 1 - 0.1)
##      Chi squared power calculation
##               w = 0.2055947
##               N = 299.3655
##              df = 2
##       sig.level = 0.05
##           power = 0.9
## NOTE: N is the number of observations

# ============================================================
# 10.6 直线相关分析（课本例36-9）
# ============================================================
# 血硒与发硒相关系数0.8, α=0.05, β=0.1
pwr.r.test(r = 0.8, sig.level = 0.05, power = 1 - 0.1,
           alternative = "two.sided")
##      approximate correlation power calculation (arctangh transformation)
##               n = 11.16238
##               r = 0.8
##       sig.level = 0.05
##           power = 0.9
##     alternative = two.sided
```

## 结果解读指南

**pwr.t.test() 输出解读：**

- `n`：每组所需样本量。单样本和配对设计是总样本量；两样本设计是每组样本量（注意 `NOTE: n is number in *each* group`）
- `d`：效应值，即标准化的均值之差 d = (μ1 - μ2) / σ。0.2=小，0.5=中，0.8=大
- `sig.level`：设定的显著性水平（α），默认0.05
- `power`：统计功效（1-β），反映探测到真实差异的概率

**pwr.2p.test() / pwr.p.test() 输出解读：**

- `h`：经arcsine变换后的效应值，由 `ES.h(p1, p2)` 计算
- `n`：每组样本量；单样本率检验(`pwr.p.test`)中为总样本量

**pwr.chisq.test() 输出解读：**

- `w`：由 `ES.w2()` 计算的行×列表效应值
- `N`：总样本量（非每组样本量），注意 `NOTE: N is the number of observations`
- `df`：自由度 = (行数-1) × (列数-1)

**pwr.r.test() 输出解读：**

- `r`：预期相关系数
- `n`：所需总样本量

**通用原则：** 计算出的 `n` 一律向上取整（如 59.35 → 60）。考虑脱落率时，n_最终 = n / (1 - 脱落率)，例如脱落率20% → n / 0.8。

## 常见问题与注意事项

**Q1：pwr包和PASS软件，哪个更好？**

PASS更专业、操作更简单（点选式），支持更复杂的设计类型，但没有Mac版且收费。pwr包免费跨平台，适合常见的研究设计。对于复杂设计（如多组ANOVA的f值计算、等效性/非劣效设计），推荐PASS。

**Q2：效应值从哪里来？**

三个来源：（1）文献综述——查相似研究的效应值；（2）预实验——用小样本数据估计；（3）Cohen约定标准——无先验信息时，按小/中/大选择。小效应需要更大样本量。

**Q3：pwr.t.test 中的 n 是每组还是总共？**

- `type="two.sample"` 时，`n` 是 **每组** 例数（注意输出中的NOTE）
- `type="one.sample"` 或 `"paired"` 时，`n` 是总例数

**Q4：pwr.anova.test 为什么我的结果和课本不一样？**

效应量 f 的计算在R中比较困难，f = sqrt(Σ(μi - μgrand)²/k) / σ_within，需要组间方差和组内方差（pooled SD）。pwr包本身不提供便捷的计算函数，因此ANOVA样本量计算推荐使用PASS。

**Q5：单侧还是双侧检验？**

默认使用双侧（`alternative="two.sided"`），除非有明确的研究假设表明只能朝一个方向变化。单侧检验所需样本量更少，但审稿人可能会质疑。

**Q6：power 设多少合适？**

通常0.80（β=0.2）或0.90（β=0.1）。power越高需要的样本量越大。基金申请和新药临床试验通常要求≥0.80，确证性临床试验建议≥0.90。

**Q7：base R 的 power.t.test() 和 pwr.t.test() 有什么区别？**

- `power.t.test()` 使用原始尺度参数（delta, sd），不需要手动算d
- `pwr.t.test()` 使用标准化效应值 d，更统一
- 两者在t检验场景下结果等价；pwr包覆盖的研究设计类型更多

**Q8：样本量计算结果不取整怎么办？**

一律向上取整。如59.35 → 每组至少60例。取整后可以反算实际power，确保不低于设定值。
