---
name: medical-stat-gee
description: "R语言医学统计：广义估计方程（GEE）。用于分析纵向数据/重复测量数据的边际模型（population-averaged），处理分类结局和连续结局，使用geepack包，可指定多种作业相关矩阵。TRIGGER when user mentions GEE、广义估计方程、作业相关矩阵、work correlation、population-averaged model、纵向数据分析，or asks about marginal models for clustered data. SKIP for 混合效应模型（个体水平）、重复测量方差分析。"
---

# 广义估计方程 (Generalized Estimating Equations, GEE)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用 GEE：**

- 纵向数据/重复测量数据：对同一组受试者在多个时间点进行重复测量
- 结局变量可以是连续型（用 `family = gaussian`）或二分类/计数（用 `family = binomial` / `family = poisson`）
- 数据存在缺失值（GEE 对随机缺失 MAR 较为稳健）
- 不满足重复测量方差分析的前提条件（残差正态性、球对称假设等）
- 研究目标是人群平均效应（population-averaged effect），而非个体水平效应
- 存在需要调整的混杂因素，且允许加入交互项

**何时不用 GEE，改用其他方法：**

- 需要个体水平的随机效应推断 → 使用混合效应模型（`lme4::glmer`、`nlme::lme`）
- 因变量是连续正态且满足球对称假设，无缺失值 → 重复测量方差分析（`anova_test`）
- 仅有一个时间点或数据独立 → 广义线性模型（`glm`）
- 需要完整的似然函数进行模型比较（AIC/BIC）→ 混合效应模型（GEE 并非基于似然，AIC 不可用，只能用 QIC）

**常见医学应用：**

- 临床试验中多次随访的疗效评估（如抗抑郁药物在不同时间点的有效率）
- 慢性病管理研究中多次测量的生理指标变化
- 流行病学队列研究中的纵向结局分析

## 前置条件

**必需的 R 包：**

```r
install.packages("geepack")     # GEE 核心包
install.packages("ggplot2")     # 可视化
install.packages("dplyr")       # 数据操作
install.packages("broom")       # 结果整理（OR 和 95%CI 提取）
```

**数据格式要求：**

- 必须是长格式（long format）：每行代表一个受试者在一个时间点的一次测量
- 必须包含唯一标识个体的 `id` 变量（因子型或字符型均可）
- 必须包含时间变量（如 `time`，可以是数值型或因子型）
- 数据集按 `id` 排序更为安全，虽然 `geeglm` 会尝试自动处理

**统计假设：**

- 独立性：不同个体之间相互独立，但同一个体内的重复测量可相关
- 边际模型假设：结局的期望值与线性预测项之间的关系通过指定的连接函数（link function）正确设定
- 作业相关矩阵不必完全正确：GEE 使用稳健标准误差（sandwich estimator），即使作业相关矩阵指定错误，回归系数的估计仍然一致
- 缺失机制：假设为 MCAR 或 MAR 时结果相对稳健

## 方法选择决策树

```
你的纵向/重复测量数据 →
├── 结局变量为连续型（正态） →
│   ├── 球对称假设成立 + 无缺失值 → 重复测量方差分析
│   └── 不满足上述条件 → GEE (family = gaussian)
├── 结局变量为二分类（0/1） →
│   └── GEE (family = binomial)
├── 结局变量为计数（次数/事件数） →
│   └── GEE (family = poisson)
└── 需要个体水平推断 / 需要 AIC 比较模型 →
    └── 混合效应模型（而非 GEE）

作业相关矩阵 (corstr) 的选择 →
├── 不确定 → 先使用 "independence"（稳健标准误差会校正）
├── 测量次数少 + 间隔相近 → "exchangeable"（等相关结构）
├── 有自然时间顺序 + 间隔不等 → "ar1"（一阶自相关）
├── 测量次数少 + 无结构规律 → "unstructured"（无结构相关）
└── 不确定最优 → 遍历所有 corstr，选 QIC 最小者

时间变量的编码 →
├── 时间点等距 + 假设线性趋势 → 数值型 time
└── 时间点不规则 / 想捕捉非线性变化 → factor(time)
```

## 标准工作流

### 步骤1：数据准备与探索

- 将数据读入为 long format，确认包含 `id` 列和 `time` 列
- 将分类变量设为因子型（`factor()`），指定参考水平（`relevel()`）
- 使用 `tapply()` 或 `group_by() %>% summarise()` 计算各分组在各时间点的结局汇总统计量
- 绘制随时间变化的分组折线图（含置信区间），初步判断交互作用是否存在

### 步骤2：模型构建

- 使用 `geepack::geeglm()` 建立 GEE 模型
- 公式中必须指定 `id = id`（指定聚类变量）
- 指定 `family`（binomial / gaussian / poisson）
- 从图中判断是否有交互作用，如有则加入 `drug*time` 交互项
- 选择合适的 `corstr`（可从 "independence" 开始作为基准）

### 步骤3：作业相关矩阵选择（可选但推荐）

- 遍历 `c("independence", "exchangeable", "ar1", "unstructured")` 
- 对每个模型计算 `QIC()`（类似 AIC 的准则，用于 GEE 模型比较）
- 选择 QIC 最小的相关结构作为最终模型

### 步骤4：结果解读

- `summary()` 查看系数估计值、标准误差、Wald 检验 P 值
- 使用 `broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)` 提取 OR 及其 95%CI
- 若有交互项，计算不同条件下的组合效应

### 步骤5：论文中如何报告

示例：本研究采用广义估计方程（GEE）分析两种抗抑郁药物在三个随访时间点对治疗有效率的影响。模型使用 logit 链接函数和二项分布族，作业相关矩阵采用独立结构。结果显示，在调整诊断严重程度后，新药与时间存在显著交互作用（OR = 2.77, 95% CI: 1.92–3.98, p < 0.001），表明新药的抗抑郁疗效随治疗时间显著增强。重度诊断患者的治疗有效率为轻度患者的 27%（OR = 0.27, 95% CI: 0.20–0.36, p < 0.001）。

## 代码示例

### 完整可运行代码

```r
# 加载包
library(geepack)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)

# 读取数据
dat <- read.csv("datasets/depression.csv", stringsAsFactors = TRUE)
dat$id <- factor(dat$id)
dat$drug <- relevel(dat$drug, ref = "standard")
head(dat, n = 3)
##   diagnose     drug id time depression
## 1     mild standard  1    0          1
## 2     mild standard  1    1          1
## 3     mild standard  1    2          1

# 变量说明：
# diagnose: 抑郁症严重程度 (mild/severe)
# drug: 药物类型 (standard/new)
# id: 受试者编号
# time: 随访时间点 (0, 1, 2)
# depression: 疗效 (1=有效, 0=无效)

# ----- 数据探索 -----
# 各分组各时间点的有效率
with(dat, tapply(depression, list(diagnose, drug, time), mean)) |>
  ftable() |>
  round(2)
##                     0    1    2
##
## mild   standard  0.51 0.59 0.68
##        new       0.53 0.79 0.97
## severe standard  0.21 0.28 0.46
##        new       0.18 0.50 0.83

# 计算汇总数据用于画图
summary_dat <- dat %>%
  group_by(diagnose, drug, time) %>%
  summarise(n = n(),
            normal_rate = mean(depression),
            .groups = "drop") %>%
  mutate(se = sqrt(normal_rate * (1 - normal_rate) / n),
         lower = pmax(0, normal_rate - 1.96 * se),
         upper = pmin(1, normal_rate + 1.96 * se))

# 绘制有效率随时间变化图
ggplot(summary_dat, aes(time, normal_rate, color = drug)) +
  geom_line(linewidth = 1.3, aes(linetype = drug)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 0.9) +
  facet_wrap(vars(diagnose),
             labeller = labeller(diagnose = c(mild = "轻度抑郁", severe = "重度抑郁"))) +
  labs(
    title = "图1：两种药物治疗轻重度抑郁症各时间点有效率比较（95%置信区间）",
    x = "记录时间",
    y = "有效率",
    color = "药物类型",
    linetype = "药物类型"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = c("new" = "#2E8B57", "standard" = "#DC143C"))

# ----- 建立 GEE 模型 -----
# 图中可见药物和时间存在交互作用，在模型中纳入交互项
dep_gee <- geeglm(depression ~ diagnose + drug * time,
                  data = dat,
                  id = id,
                  family = binomial,         # 二分类结局
                  corstr = "independence")    # 作业相关矩阵

# ----- 模型结果 -----
summary(dep_gee)
##
## Call:
## geeglm(formula = depression ~ diagnose + drug * time, family = binomial,
##     data = dat, id = id, corstr = "independence")
##
##  Coefficients:
##                Estimate  Std.err   Wald Pr(>|W|)
## (Intercept)    -0.02799  0.17419  0.026    0.872
## diagnosesevere -1.31391  0.14598 81.006  < 2e-16 ***
## drugnew        -0.05960  0.22854  0.068    0.794
## time            0.48241  0.11994 16.179 5.76e-05 ***
## drugnew:time    1.01744  0.18769 29.385 5.93e-08 ***
## ---
## Correlation structure = independence
## Number of clusters:   340  Maximum cluster size: 3

# 提取 OR 值和 95%CI（对数尺度结果指数化）
broom::tidy(dep_gee, exponentiate = TRUE, conf.int = TRUE)
## # A tibble: 5 × 7
##   term           estimate std.error statistic  p.value conf.low conf.high
##   <chr>             <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
## 1 (Intercept)       0.972     0.174    0.0258 8.72e-1    0.691     1.37
## 2 diagnosesevere    0.269     0.146   81.0    0          0.202     0.358
## 3 drugnew           0.942     0.229    0.0680 7.94e-1    0.602     1.47
## 4 time              1.62      0.120   16.2    5.76e-5    1.28      2.05
## 5 drugnew:time      2.77      0.188   29.4    5.93e-8    1.91      4.00

# ----- QIC：选择最优作业相关矩阵 -----
corstrs <- c("independence", "exchangeable", "ar1", "unstructured")
qics <- list()
for (i in 1:length(corstrs)) {
  dep_gee1 <- geeglm(depression ~ diagnose + drug * time,
                     data = dat, id = id, family = binomial,
                     corstr = corstrs[i])
  qics[[i]] <- QIC(dep_gee1)
}
do.call(rbind, qics)
##       QIC QICu Quasi Lik   CIC params QICC
## [1,] 1172 1172      -581 5.140      5 1172
## [2,] 1172 1172      -581 5.139      5 1172
## [3,] 1172 1172      -581 5.139      5 1172
## [4,] 1172 1172      -581 5.087      5 1173
```

## 结果解读指南

**模型输出各部分含义：**

| 输出项 | 含义 |
|--------|------|
| `Estimate` | 回归系数 β，在 logit（对数优势）尺度上 |
| `Std.err` | 稳健标准误差（sandwich estimator），已校正聚类内的相关性 |
| `Wald` | Wald 检验统计量，用于检验系数是否为 0 |
| `Pr(>|W|)` | P 值，< 0.05 表示该变量与结局的关联有统计学意义 |
| `Correlation structure` | 模型使用的作业相关矩阵类型 |
| `Number of clusters` | 独立个体数（受试者数量） |
| `Maximum cluster size` | 单个个体的最大重复测量次数 |

**如何用中文解释结果（以 depression 数据为例）：**

- `diagnosesevere` 的系数为 -1.31，OR = exp(-1.31) = 0.27。解释：在相同药物和随访时间下，重度抑郁患者的治疗有效率为轻度患者的 27%，差异有统计学意义（P < 0.001）。**基线诊断越严重，治疗有效的可能性显著越低。**

- `drugnew` 的系数为 -0.06，OR = 0.94，P = 0.79。解释：在治疗开始时（time = 0），新药组与标准药组的有效率无统计学差异，**两组基线均衡**。

- `time` 的系数为 0.48，OR = exp(0.48) = 1.62，P < 0.001。解释：使用标准药时，每增加一个时间单位，治疗有效率约增加 62%。

- `drugnew:time` 交互项系数为 1.02，OR = 2.77，P < 0.001。解释：新药组随时间改善的速度是标准药组的 2.77 倍。在 time = 2 时，新药组有效率是标准药组的 exp(-0.0596 + 1.017 × 2) ≈ 7.2 倍。

- QIC 值：四个相关矩阵的 QIC 非常接近（均为 1172），表示本数据中相关结构的选择对模型影响很小。选择 QIC 最小的即可。

## 常见问题与注意事项

**Q: GEE 和混合效应模型的区别是什么？**

A: GEE 是边际模型（population-averaged），回答“总体平均效果”的问题；混合效应模型是条件模型（subject-specific），回答“对特定个体效果”的问题。如果研究目的是做个体预测或需要随机效应估计，用混合效应模型。如果是检验组间平均差异（如临床试验疗效比较），GEE 更加常用。

**Q: 作业相关矩阵如何选择？**

A: GEE 的回归系数估计对作业相关矩阵的错误指定具有稳健性（sandwich estimator）。一般策略：（1）先用独立结构作基准;（2）遍历多个结构计算 QIC，选择 QIC 最小的;（3）如果 QIC 相差很小（如本例），选最简单的（independence 或 exchangeable）。

**Q: 时间变量用数值型还是因子型？**

A: 时间点等距且假设线性变化趋势时，使用数值型 `time`（节约自由度）;如果时间点不规则或想捕捉非线性趋势（如先升后降），使用 `factor(time)` 并在模型中纳入与分组变量的交互项。

**Q: GEE 可以处理缺失值吗？**

A: GEE 在随机缺失（MAR）假设下结果较为稳健，但缺失值会降低统计效能。若缺失为非随机缺失（MNAR），需要考虑敏感性分析。

**Q: geepack 和 gee 包有什么区别？**

A: `geepack` 更现代化，支持 `QIC()` 和多种连接函数，是首选。`gee` 包较为老旧但也可使用。`glmtoolbox` 提供类似功能。三者用法基本一致。

**Q: SPSS 中如何做 GEE？**

A: SPSS 中通过「分析 → 广义线性模型 → 广义估计方程」实现，操作路径与 R 中的 `geeglm` 逻辑一致。SPSS 输出以 QIC 替代 AIC，Wald 卡方替代 F 检验。
