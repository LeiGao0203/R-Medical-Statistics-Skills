---
name: medical-stat-subgroup-forest
description: "R语言医学统计：亚组分析与多因素回归森林图的区别及实际应用。厘清两种常被混淆的森林图类型——亚组分析的森林图（每亚组内单独回归）vs 多因素回归的森林图（同一模型中哑变量各水平与参考组比较）。TRIGGER when user mentions 亚组分析森林图 vs 多因素森林图、两种森林图的区别、森林图类型选择、亚组分析还是多因素回归，or asks about distinguishing subgroup forest plots from multivariable regression forest plots. SKIP for 单纯亚组分析绘制（见subgroup-analysis技能）、单纯多因素回归建模（见对应回归技能）、Meta分析森林图。"
---

# 亚组分析和多因素回归的森林图比较 (Subgroup Analysis vs Multivariable Forest Plot)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

你是否在论文中见过这两种森林图——它们都包含 HR/OR 值、95% 可信区间、P 值，看起来几乎一模一样？实际上，它们表达的含义截然不同：

**多因素回归的森林图：**
- 将多个变量同时纳入一个回归模型，展示每个变量（经哑变量编码后）相对于参考水平的效应量
- 本质是「一个回归模型的结果可视化」
- 图中每个分类变量有一个 reference 行（HR=1，无 P 值、无 CI）

**亚组分析的森林图：**
- 在数据的不同子集（如男/女、老/少）中分别做单因素或多因素回归，汇总各子集的结果
- 本质是「多个独立分析的结果汇总」
- 每行都是一个完整分析的结果，均有自己的 HR/CI/P 值

**不使用本方法的情况：**
| 你的需求 | 应使用的方法 |
|----------|-------------|
| 仅需要做亚组分析绘图 | 参见 subgroup-analysis / subgroup-onecode 技能 |
| 仅需要做多因素回归建模 | 参见 logistic-reg / survival / multi-reg 技能 |
| Meta 分析中的合并效应量森林图 | Meta 分析专用方法（meta 包） |
| 展示单因素回归中每个变量的效应 | 多因素回归森林图（本技能 Section 1），仅用单个变量即可 |

## 前置条件

**R 包安装：**

```r
install.packages(c("tidyverse", "survival", "survminer", "broom", "forestploter"))
```

**数据格式要求：**
- 生存数据：需含时间变量和状态变量（如 `time` 和 `status`，0=删失/1=事件）
- 暴露/处理变量：至少一个二分类或多分类分组变量
- 协变量：所有需纳入回归的变量需为数值型或 factor 型
- 亚组分析需有可用于分层的分组变量（如性别、年龄组）

**统计前提：**
- 多因素回归需满足所用回归模型的基本假定（Cox 需满足 PH 假定，逻辑回归需满足独立性等）
- 亚组分析不改变原始回归模型的假定，但每个亚组内样本量需足够——某亚组样本量过少会导致估计不稳定甚至无法收敛

## 方法选择决策树

```
你的分析目标 →
├── 想在同一张图中展示一个回归模型中各分类变量的效应量
│   ├── 分类变量有 ≥2 个水平 → 多因素回归森林图（使用 ggforest / forestploter）
│   │   └── 需为每个分类变量设置参考水平，参考行 HR=1, CI=NA, P=NA
│   └── 变量为连续型 → 无需哑变量编码，直接报告 per unit 或 per SD 效应
├── 想知道某一暴露因素在不同人群亚组中效应是否一致
│   ├── 暴露 × 分组的交互项 P<0.05 → 效应修饰存在，报告亚组分析森林图
│   └── 暴露 × 分组的交互项 P≥0.05 → 可能不需要展示亚组森林图（视期刊要求）
├── 审稿人要求「多因素回归 + 亚组分析的森林图」
│   └── 在每个亚组内分别做多因素回归 → 每个亚组内展示多个 HR/CI → 需自行构建数据框
└── 仅需可视化单因素回归中某变量的各分类水平效应 → 两者均可，推荐多因素回归森林图方式
```

## 标准工作流

### 步骤1：数据准备
- 加载 `survival::colon` 数据集（或你自己的数据）
- 选择分析人群，将分类变量全部转为 factor 类型
- 明确因变量（时间+状态）和自变量（暴露+协变量）

### 步骤2：多因素回归及森林图
- 构建含全部协变量的回归模型（Cox / Logistic / 线性）
- 用 `survminer::ggforest()` 快速出图（Cox 回归专用）
- 或用 `forestploter` 包自定义数据框绘制（通用方法）

### 步骤3：亚组分析及森林图
- 确定分层变量（如性别、年龄分组）
- 在每个亚组内分别做回归分析
- 提取每个亚组的 HR/OR、95%CI、P 值
- 用 `forestploter::forest()` 汇总绘制

### 步骤4：对比两种森林图的结果解读
- 多因素回归森林图中：同一变量不同水平的 HR 反映的是哑变量编码后与参考水平的比较
- 亚组分析森林图中：同一暴露变量在不同亚组中的 HR，反映的是效应修饰（交互作用）
- 论文中必须按「你所采用的方法对应的含义」去解读，不可混用

### 步骤5：论文报告格式
- 多因素回归森林图："图 X 展示了多因素 Cox 回归中各自变量与生存结局的关联。以 Obs 组为参考，Lev+5FU 组的 HR=1.68 (95%CI: 1.33-2.13, P<0.001)……"
- 亚组分析森林图："图 X 展示了 xxx 暴露与结局的亚组分析结果。P for interaction = 0.xx，提示效应修饰不（具）有统计学意义……"

## 代码示例

### 多因素回归森林图（完整流程）

```r
rm(list = ls())
library(tidyverse)
library(survival)
library(broom)
library(forestploter)

# 1. 数据准备
# survival::colon 数据集：1858行×16列，结肠癌生存数据
# 变量：rx(治疗), sex, age, obstruct, perfor, adhere, differ, extent, surg, node4
# 因变量：time(生存时间), status(0=删失,1=事件)

df <- colon %>%
  mutate(rx = as.numeric(rx)) %>%
  filter(etype == 1, !rx == 2) %>% # 仅保留Obs组和Lev+5FU组
  select(time, status, rx, sex, age, obstruct, perfor, adhere,
         differ, extent, surg, node4) %>%
  mutate(
    sex     = factor(sex, levels = c(0, 1), labels = c("female", "male")),
    age     = factor(ifelse(age > 65, ">65", "<=65"), levels = c(">65", "<=65")),
    obstruct = factor(obstruct, levels = c(0, 1), labels = c("No", "Yes")),
    perfor  = factor(perfor, levels = c(0, 1), labels = c("No", "Yes")),
    adhere  = factor(adhere, levels = c(0, 1), labels = c("No", "Yes")),
    differ  = factor(differ, levels = c(1, 2, 3),
                     labels = c("well", "moderate", "poor")),
    extent  = factor(extent, levels = c(1, 2, 3, 4),
                     labels = c("submucosa", "muscle", "serosa", "contiguous")),
    surg    = factor(surg, levels = c(0, 1), labels = c("short", "long")),
    node4   = factor(node4, levels = c(0, 1), labels = c("No", "Yes")),
    rx      = factor(ifelse(rx == 3, 0, 1), levels = c(0, 1))
  )
str(df)
## 'data.frame': 619 obs. of 12 variables

# 2. 多因素Cox回归
fit_multi <- coxph(Surv(time, status) ~ ., data = df)
summary(fit_multi)
## n= 606, number of events= 292
##                       coef exp(coef) se(coef)     z Pr(>|z|)
## rx1               0.521198  1.684043 0.120261 4.334  1.47e-05 ***
## sexmale          -0.125724  0.881858 0.118615 -1.060   0.2892
## node4Yes          0.811284  2.250796 0.123699  6.559  5.43e-11 ***

# 3. 快速出图：survminer::ggforest (Cox专用)
library(survminer)
ggforest(fit_multi, fontsize = 1)
# 注意：ggforest 自动为每个分类变量添加 reference 行

# 4. 自定义森林图：forestploter (通用，支持任意回归模型)
multidf <- broom::tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2)),
    `HR(95%CI)` = paste0(estimate, "(", conf.low, "-", conf.high, ")")
  ) %>%
  select(term, estimate, p.value, conf.low, conf.high, `HR(95%CI)`)
multidf
## # A tibble: 13 × 6
##    term          estimate p.value conf.low conf.high `HR(95%CI)`
##    <chr>            <dbl>   <dbl>    <dbl>     <dbl> <chr>
##  1 rx1               1.68    0        1.33      2.13 1.68(1.33-2.13)
##  2 sexmale           0.88    0.29     0.70      1.11 0.88(0.7-1.11)
##  3 node4Yes          2.25    0        1.77      2.87 2.25(1.77-2.87)
##  ...

# 保存后手动添加 reference 行（每个分类变量的参考水平需要补一行HR=1）
# 此处直接使用已整理好的数据（含变量名分组行、reference行）
plot_df <- read.csv(file = "./datasets/multidf.csv", check.names = FALSE)

# 将P值的NA替换为空格（参考行不显示P值）
plot_df[, c(3)][is.na(plot_df[, c(3)])] <- " "
# 添加空列用于展示可信区间
plot_df$` ` <- paste(rep(" ", nrow(plot_df)), collapse = " ")

# 绘制
p <- forest(
  data    = plot_df[, c(1, 6, 7, 3)],
  lower   = plot_df$conf.low,
  upper   = plot_df$conf.high,
  est     = plot_df$estimate,
  ci_column = 3,
  sizes   = 1,
  ref_line = 1,
  xlim    = c(0.1, 4)
)
print(p)
```

### 亚组分析的数据框构建思路（伪代码）

```r
# 亚组分析：在每个亚组内分别跑Cox回归，提取结果后合并
# 详细实现见 subgroup-analysis 技能

# 核心思路示意：
# 1. 按亚组变量拆分数据
# 2. 每个亚组内拟合 coxph(Surv(time,status) ~ exposure, data = subset)
# 3. 从每个拟合中提取 HR、95%CI、P 值
# 4. 将所有亚组结果用 rbind() 合并为数据框
# 5. 用 forestploter::forest() 绘制
```

## 结果解读指南

**多因素回归森林图中各列的含义：**

| 列名 | 含义 | 关键注意点 |
|------|------|-----------|
| term（或 subgroup） | 变量名和水平 | 格式为"变量名+水平名"（如 rx1、sexmale），某些行仅为变量名作标题 |
| estimate | HR/OR（效应量） | 哑变量编码后的特定水平 vs 参考水平的比；reference 行始终为 1 |
| conf.low/conf.high | 95% 可信区间下限/上限 | reference 行 CI 为 [1, 1]（实际表示方法上通常空白或写1） |
| p.value | Wald 检验 P 值 | reference 行无 P 值（NA 或空白）：参考组无"自己比自己"的检验 |

**两类森林图的核心区别：**

| 对比维度 | 多因素回归森林图 | 亚组分析森林图 |
|----------|-----------------|---------------|
| 分析次数 | 1 次（一个模型） | 每个亚组 1 次（多个模型） |
| 每行的含义 | 同一模型中某变量某水平 vs 参考水平 | 不同亚组中暴露因素的效应 |
| 是否有 reference 行 | 是（每个多水平分类变量各一个） | 否 |
| 回答的问题 | "哪些因素是独立的危险因素？" | "暴露效应在不同人群中是否一致？" |
| 解读重点 | 每个变量的 HR/OR 大小和 P 值 | P for interaction 和效应方向的一致性 |

**论文中解读示例（对照上表）：**

多因素回归森林图："多因素分析显示，治疗方式（Lev+5FU vs Obs：HR=1.68, 95%CI 1.33-2.13, P<0.001）和阳性淋巴结>4个（HR=2.25, 95%CI 1.77-2.87, P<0.001）是独立的预后因素。"

亚组分析森林图："亚组分析显示，治疗获益在各亚组间方向一致：男性组 HR=1.55 (95%CI 1.10-2.18)，女性组 HR=1.82 (95%CI 1.28-2.59)。P for interaction = 0.45，提示性别对治疗效果无显著效应修饰。"

## 常见问题与注意事项

**Q1: 两幅森林图长得几乎一模一样，怎么一眼区分？**

看是否有 reference 行。多因素回归森林图中每个多水平分类变量必定有一行标为 reference（或有 HR=1, CI=1-1, 无 P 值），亚组分析森林图中没有这一行。另外，多因素回归森林图的标题通常是"Multivariable model"，亚组分析通常明确标注"Subgroup analysis"。

**Q2: 什么时候用多因素回归森林图，什么时候用亚组分析森林图？**

研究目的是评估多个候选因素与结局的独立关联时，用多因素回归森林图。研究目的是检查某特定暴露/干预在不同亚组中效应是否一致时，用亚组分析森林图。两者回答完全不同的临床问题。

**Q3: 能否在亚组分析的每个亚组中都做多因素回归，再画森林图？**

可以，这是"多因素回归的亚组分析"。做法是在每个亚组内都拟合多因素回归模型（含多个协变量），每个亚组中每个变量都产出自己的 HR/CI。这样画出的森林图中每个亚组下会有多行（对应多个变量）。这在某些研究中用于展示校正混杂后的亚组效应，但实现较为复杂。

**Q4: 为什么多因素回归森林图中 reference 的 HR 是 1？**

分类变量在回归中会自动变为哑变量，以某一水平为参照（HR=1），其余水平与该参照比较得出各自的 HR。因此 reference 的 HR 定义为 1，严格来说它不是一个"估计值"，而是"定义的基准"。

**Q5: 多因素回归森林图中能不能不显示 reference 行？**

可以。部分期刊不要求显示 reference 行，此时只需在构建数据框时删除 reference 对应的行即可。但保留 reference 行有助于读者理解每个变量的参照组是什么。

**Q6: 亚组分析的 P for interaction 不显著，还应该画森林图吗？**

若 P for interaction ≥ 0.05，说明各亚组效应方向/大小一致，此时亚组分析森林图不是必需的。但许多期刊仍要求展示亚组分析森林图（尤其是大型 RCT 的事后分析），即使交互不显著。建议遵从期刊要求。

**Q7: 可以用 ggplot2 代替 forestploter 绘制森林图吗？**

可以。ggplot2 更灵活，但需手动完成坐标轴、误差线、参考线的绘制。forestploter 专为森林图设计，出图更快且自带 publication-ready 的默认样式。详见"其他资源"中的相关链接。

**Q8: SPSS 中能实现这两种森林图吗？**

SPSS 原生的回归输出不直接生成森林图。多因素回归森林图在 SPSS 中需手动整理输出表格再在其他工具中绘图；亚组分析则需在 SPSS 中按分组变量执行多次回归再汇总。R 在这两类森林图的绘制上远比 SPSS 方便。
