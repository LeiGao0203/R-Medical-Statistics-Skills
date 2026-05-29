---
name: medical-stat-discrete-distrib
description: "R语言医学统计：离散型变量的概率分布。涵盖二项分布、Poisson分布及其在医学中的应用，包括概率计算、拟合优度检验。TRIGGER when user mentions 二项分布、Poisson分布、离散分布、率的置信区间、概率计算，or asks about rare event probability or binomial probability. SKIP for 连续型分布（t分布、正态分布）。"
---

# 几种离散型变量的分布及其应用 (Discrete Distributions)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用二项分布 (binomial distribution)**：
- 每次试验只有两个可能结果（成功/失败、患病/未患病、有效/无效）
- 已知固定的试验次数 n，每个观测独立，成功概率 π 恒定
- 医学应用：率的置信区间估计、单组率的假设检验、两独立样本率的比较

**使用 Poisson 分布 (Poisson distribution)**：
- 描述单位时间/面积/体积内稀有事件的发生次数
- 事件发生独立，总体均数 λ 恒定
- 医学应用：罕见病发病率、空气中颗粒物计数、单位体积细菌数、放射物计数
- 当二项分布的 n 很大而 π 很小时（n≥100, π≤0.01），二项分布近似 Poisson 分布

**不适用此方法的情况**：
- 连续型数据的比较 → 使用 t 检验或方差分析
- 分类数据的假设检验（多分类/列联表） → 使用卡方检验
- 存在混杂因素的率比较 → 使用 logistic 回归或 Poisson 回归
- 事件间不独立或有过离散现象 → 使用负二项回归

## 前置条件

```r
# 本章所有函数均来自 R 内置 stats 包，无需额外安装
library(stats)  # 内置，提供 binom.test, poisson.test, prop.test, chisq.test
```

**数据格式**：
- `binom.test`: 提供成功次数 x 和试验总次数 n
- `poisson.test`: 提供事件计数 x，可选时间基数 T 和参考率 r
- `prop.test`: 提供成功次数 x 和总数 n，或提供 2×2 矩阵/table
- `chisq.test`（拟合优度）: 提供观测频数向量 x 和理论概率向量 p

**统计假设**：
- 二项分布：n 次独立试验，每次试验概率 π 相同
- Poisson 分布：事件发生相互独立，在给定区间内发生率恒定
- 正态近似法（大样本）：要求 n≥50 且 np 和 n(1-p) 均≥5（二项）；λ≥20（Poisson）

## 方法选择决策树

```
数据类型 →
├── 已知试验次数 n，二元结果（成功/失败）
│   ├── 估计率 → 总体率的区间估计（3.1.1）
│   │   ├── 小样本 → binom.test() 精确法
│   │   └── 大样本 → prop.test() 或正态近似公式
│   ├── 单样本率与已知总体率比较 → 样本率和总体率的比较（3.1.2）
│   │   ├── 小样本 → binom.test(x, n, p=π₀, alternative=...)
│   │   └── 大样本 → prop.test(x, n, p=π₀, alternative=..., correct=FALSE)
│   └── 两独立样本率的比较 → 两样本率的比较（3.1.3）
│       ├── 原始数据 → prop.test(矩阵, correct=FALSE)
│       └── 拟合优度 → chisq.test(x=观测, p=理论概率)
│
├── 稀有事件计数，无固定 n（单位时间/空间内的发生数）
│   ├── 估计发生率 → 总体均数的区间估计（3.2.1）
│   │   ├── 小λ → poisson.test() 精确法
│   │   └── 大λ（λ≥20） → 正态近似法 x ± zα/₂·√x
│   ├── 单样本与已知总体率比较 → 样本均数和总体均数的比较（3.2.2）
│   │   ├── 小λ → poisson.test(x, T=时间基数, r=参考率, alternative=...)
│   │   │   或手动: 1 - ppois(q=x-1, lambda=T*r)
│   │   └── 大样本小率 → prop.test(x, n, p=π₀, ...) 正态近似
│   └── 两独立样本率的比较 → 两个样本均数的比较（3.2.3）
│       └── poisson.test(c(x1,x2), c(T1,T2))
│
└── 事件有过离散（方差 > 均数） → 负二项分布 / 负二项回归（略）
```

## 标准工作流

### 步骤1：确认分布类型
- 判断数据是二元计数（有固定 n）还是稀有事件计数（无固定 n）
- 对于二项问题，明确 n（试验总数）和 x（成功数），计算样本率 p = x/n
- 对于 Poisson 问题，明确 x（事件计数）和 T（时间/空间基数），计算单位率 r = x/T

### 步骤2：选择精确法或近似法
- 小样本或低计数：使用精确法（`binom.test`, `poisson.test`）
- 大样本：可使用正态近似法（`prop.test` 或手动公式）
- 正态近似条件：n·p 和 n·(1-p) ≥ 5（二项）；λ ≥ 20（Poisson）

### 步骤3：执行检验
- 确定检验方向：`"two.sided"`（双侧）、`"greater"`（大于）、`"less"`（小于）
- 设定置信水平：默认 `conf.level = 0.95`
- 连续性校正：正态近似时可设 `correct = FALSE`（更接近课本公式）

### 步骤4：结果解读
- 查看 p-value：p < 0.05 拒绝 H₀，差异有统计学意义
- 查看置信区间：区间是否包含参考值/零值
- 查看样本估计值：反映效应的方向和大小

### 步骤5：结果报告
- 报告样本率/发生率及 95% CI
- 报告检验统计量、p 值（精确到小数点后 3-4 位）
- 用专业中文表述结论（见各示例）

## 代码示例

以下代码可直接复制到 R 中运行。所有函数均为 base R 内置，无需额外安装。

```r
library(stats)

# ========== 二项分布 ==========

# --- 3.1.1 总体率的区间估计 ---

# 例6-2: 13名术后妇女6人受孕，95%CI (小样本 - 精确法)
binom.test(x = 6, n = 13, conf.level = 0.95)
## 95 percent confidence interval:
##  0.1922324 0.7486545
## sample estimates:
## probability of success
##              0.4615385

# 例6-3: 100人有55人有效，95%CI (大样本 - 正态近似法)
Sp <- sqrt((0.55 * (1-0.55)) / 100)
0.55 + 1.96 * Sp
## [1] 0.6475088
0.55 - 1.96 * Sp
## [1] 0.4524912

# 或者直接使用 prop.test()
prop.test(x = 55, n = 100, correct = FALSE)
## 95 percent confidence interval:
##  0.4524460 0.6438546

# --- 3.1.2 样本率和总体率的比较 ---

# 例6-4: 已知受孕率0.55，10人中9人受孕，A方法是否优于B？(单侧 - 精确法)
binom.test(x = 9, n = 10, p = 0.55, alternative = "greater")
## p-value = 0.02326
## 95 percent confidence interval:
##  0.6058367 1.0000000
## 结论：拒绝H0，A方法受孕率高于B方法

# 例6-5: 已知治愈率0.6，10人中9人有效，两种药物是否不同？(双侧 - 精确法)
binom.test(x = 9, n = 10, p = 0.6, alternative = "two.sided")
## p-value = 0.05865
## 结论：不拒绝H0，尚不能认为两种药物的疗效不同

# 例6-6: 180名患者治愈117人，新方法是否优于常规？(大样本 - 正态近似法)
prop.test(x = 117, n = 180, p = 0.45, alternative = "greater", correct = FALSE)
## X-squared = 29.091, p-value = 3.453e-08
## 95 percent confidence interval:
##  0.5896943 1.0000000
## 结论：拒绝H0，新疗法效果优于常规疗法

# --- 3.1.3 两样本率的比较 ---

# 例6-7: 颈椎病发病的性别差异 (两独立样本率比较)
t67 <- matrix(c(36, 84, 22, 88), ncol = 2, byrow = TRUE,
              dimnames = list(c("male", "female"), c("success", "failure")))
t67
##        success failure
## male        36      84
## female      22      88

prop.test(t67, correct = FALSE)
## X-squared = 3.0433, df = 1, p-value = 0.08107
## 95 percent confidence interval:
##  -0.01095102  0.21095102
## sample estimates: prop 1 = 0.3, prop 2 = 0.2
## 结论：不拒绝H0，尚不能认为该职业人群颈椎病的发病有性别差异

# 例6-8: 家族集聚性分析 (拟合优度卡方检验)
x <- c(26, 10, 28, 18)
p <- c(0.13265, 0.38235, 0.36735, 0.11765)
chisq.test(x = x, p = p)
## X-squared = 42.949, df = 3, p-value = 2.523e-09
## 结论：拒绝H0，认为此种疾病存在家族集聚性

# ========== Poisson分布 ==========

# --- 3.2.1 总体均数的区间估计 ---

# 例6-10: 粉尘21个/升，95%CI 和 99%CI (精确法)
poisson.test(x = 21)
## 95 percent confidence interval:
##  12.99933 32.10073

poisson.test(x = 21, conf.level = 0.99)
## 99 percent confidence interval:
##  11.06923 35.94628

# 例6-11: 大样本正态近似法 (λ=68, λ≥20 可用)
68 - 1.96 * sqrt(68)
## [1] 51.83743
68 + 1.96 * sqrt(68)
## [1] 84.16257

# --- 3.2.2 样本均数和总体均数的比较 ---

# 例6-12: 母亲吸烟与先心病 (单侧 - 精确法)
# 120名吸烟母亲子女中4例先心病，一般人群发病率 0.008/人
poisson.test(x = 4, T = 120, r = 0.008, alternative = "greater")
## p-value = 0.01663
## 95 percent confidence interval:
##  0.01138599        Inf
## sample estimates: event rate = 0.03333333

# 手动计算：P(X ≥ 4 | λ = 120*0.008 = 0.96)
1 - ppois(q = 4 - 1, lambda = 120 * 0.008)
## [1] 0.01663305
## 结论：拒绝H0，母亲吸烟会增加小孩先心病的发病风险

# 例6-13: 亲缘血统与精神发育不全 (大样本正态近似)
prop.test(x = 123, n = 25000, p = 0.003, alternative = "greater", correct = FALSE)
## X-squared = 30.812, p-value = 1.421e-08
## 95 percent confidence interval:
##  0.004243748 1.000000000
## 结论：拒绝H0，有亲缘血统婚配关系的后代精神发育不全发生率高于一般人群

# --- 3.2.3 两个样本均数的比较 ---

# 例6-14: 两种纯净水大肠杆菌数比较 (两Poisson率比较)
poisson.test(c(4, 7), c(1, 1))
## count1 = 4, expected count1 = 5.5, p-value = 0.5488
## 95 percent confidence interval:
##  0.1226664 2.2477580
## sample estimates: rate ratio = 0.5714286
## 结论：不拒绝H0，尚不能认为有差别

# 例6-15: 罕见病发病的地域差异 (不同时间基数的两样本率比较)
poisson.test(c(32, 12), c(4, 3))
## count1 = 32, expected count1 = 25.143, p-value = 0.04653
## 95 percent confidence interval:
##  1.002761 4.264145
## sample estimates: rate ratio = 2
## 结论：拒绝H0，认为该疾病的发病存在地域性差异
```

## 结果解读指南

**binom.test 输出解读**：
- `number of successes = x, number of trials = n`: 确认输入的观测数据
- `p-value`: 双侧检验默认检验 H₀: π = 0.5（未指定 p 参数时）；若指定 p 参数则检验 H₀: π = p
- `X percent confidence interval`: 此处为 Clopper-Pearson 精确置信区间，非近似区间
- `sample estimates: probability of success`: 样本率 p̂ = x/n

**prop.test 输出解读**：
- `X-squared`: 卡方统计量（基于正态近似的 z²）
- `p-value`: p < 0.05 拒绝 H₀，差异有统计学意义
- `X percent confidence interval`: 基于 Wilson 评分法的率置信区间
- 两样本比较时：`prop 1` 和 `prop 2` 分别为两组样本率

**poisson.test 输出解读**：
- `number of events = x, time base = T`: 事件计数和时间/空间基数
- `event rate = x/T`: 单位发生率估计
- `rate ratio`: 两样本比较时的率比（rate ratio），RR = 1 表示无差异
- 置信区间不包含 1 → p < 0.05 → 差异有统计学意义

**chisq.test（拟合优度）输出解读**：
- `X-squared`: 卡方统计量 = Σ[(O-E)²/E]
- `df`: 自由度 = 类别数 - 1
- p < 0.05 表示观测分布与理论分布差异有统计学意义

**结论表述模板**：
- p < 0.05: "拒绝 H₀，接受 H₁，可认为……差异有统计学意义"
- p ≥ 0.05: "不拒绝 H₀，尚不能认为……存在差异"（注意：不能说"无差异"，只能说证据不足）

## 常见问题与注意事项

**Q1: binom.test 和 prop.test 有什么区别？**
- `binom.test`: 基于二项分布的精确检验，任何样本量均可用。小样本时推荐使用。
- `prop.test`: 基于正态近似的卡方检验（与 chisq.test 等价）。大样本时与精确法结果接近。支持连续性校正参数 `correct`。

**Q2: 什么时候用 Poisson 分布而不是二项分布？**
- 二项分布：已知试验总次数 n，关注的是"n 次中有多少次成功"
- Poisson 分布：无固定试验次数，关注的是"一定时间/空间内发生多少次"
- 当 n ≥ 100 且 π ≤ 0.01 时，二项分布可近似为 Poisson 分布（λ = nπ）

**Q3: 两样本率比较时 prop.test 的正确用法？**
- 传入 2×2 矩阵时自动进行两样本率比较（如例6-7）
- 矩阵的行=分组，列=发生/未发生（或成功/失败）
- `correct = FALSE` 结果更接近课本按公式手算的值

**Q4: poisson.test 中 T 参数是什么含义？**
- T 是 time base（时间/空间基数），默认为 1
- 如每 120 人中有 4 例，T=120 表示观测人年数，rate = 4/120 = 0.0333/人
- 两样本比较时：`poisson.test(c(x1,x2), c(T1,T2))` 比较两种率

**Q5: ppois 与 poisson.test 计算 p 值有何不同？**
- `poisson.test(x=4, T=120, r=0.008, alternative="greater")` 是精确 Poisson 检验
- `1 - ppois(q=4-1, lambda=120*0.008)` 计算 P(X ≥ 4)，λ = 120 × 0.008 = 0.96
- 注意 `ppois(q, lambda)` 计算 P(X ≤ q)，故用 `q = x-1` 求上侧概率

**Q6: chisq.test 用于拟合优度时 p 参数如何设定？**
- p 是理论概率向量，必须为非负数且和为 1
- `chisq.test(x, p)` 检验观测频数 x 是否符合理论分布 p
- 若 x 中有 0 值（即某类别观测频数为 0），可能需要合并类别或使用精确方法
