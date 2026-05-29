---
name: medical-stat-ps-stratification
description: "R语言医学统计：倾向性评分回归调整与分层分析。将倾向性评分作为协变量纳入回归模型，或按评分分层分析，是PSM的替代方法。TRIGGER when user mentions 倾向性评分回归、倾向性评分分层、PS as covariate、ps regression adjustment，or asks about non-matching propensity score methods. SKIP for PSM匹配、IPTW加权。"
---

# 倾向性评分：回归和分层 (Propensity Score: Regression & Stratification)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

倾向性评分回归和分层适用于：

- 观察性研究中存在多个混杂因素，需要通过倾向性评分进行控制
- 两组基线资料不均衡，但数据量不足以做匹配，或匹配会损失大量样本
- 不需要像 PSM 那样严格 1:1 匹配的灵活场景
- 倾向性评分回归：最简单直接的方法，将 PS 作为协变量纳入回归模型，适合快速校正
- 倾向性评分分层：按 PS 分层后在每层内分析，再汇总效应，适合处理因素和对照组 PS 分布重叠较好的情况

不适用场景：
- 处理组和对照组 PS 范围几乎完全不重叠 → 分层无意义，考虑 PSM 或 IPTW
- 需要估计 ATT（处理组平均处理效应）→ 倾向性评分匹配更合适
- 需要模拟随机对照试验的平衡效果 → PSM 更直接
- 数据量极大且对因果推断要求高 → IPTW 可能是更优选择

## 前置条件

**R 包依赖**：

```r
install.packages(c("tidyverse", "broom"))
```

**数据要求**：
- 二分类处理因素变量（如 0/1 编码）
- 连续型或可适当编码为数值的结局变量
- 多个混杂因素变量（连续型或分类型，分类型需转为 factor）
- 数据中无缺失值（或已做插补处理）

**方法假设**：
- 无不可观测混淆（strongly ignorable treatment assignment）：所有影响处理分配和结局的混杂因素都已被测量并纳入
- 倾向性评分模型（通常是逻辑回归）正确指定
- 分层分析要求各层内两组的 PS 范围有足够重叠（common support）

## 方法选择决策树

```
你的数据情况 →
├── 只需快速校正混杂，不关心样本损失 →
│   └── 使用倾向性评分回归（PS 作为协变量直接纳入回归）
├── 希望分亚组分析，检查各层效果是否一致 →
│   └── 使用倾向性评分分层（分 5~10 层）
├── 两组 PS 范围大量不重叠 →
│   └── 不建议分层和回归，换用 PSM 或 IPTW
└── 处理因素为多分类 →
    └── 使用多分类逻辑回归计算 GPS（广义倾向性评分），再回归或分层
```

## 标准工作流

### 步骤 1：数据准备与探索

检查数据维度、变量类型，明确处理因素、结局变量、混杂因素。查看两组的基线情况，确认存在不均衡，需要校正。

```r
library(tidyverse)

ecls <- read.csv("datasets/ecls.csv") %>%
  dplyr::select(c5r2mtsc_std, catholic, race_white, w3momed_hsb, p5hmage,
                w3momscr, w3dadscr) %>%
  na.omit()

dim(ecls)
# [1] 5548    7

glimpse(ecls)
```

检查两组间结局变量均值和混杂因素分布差异：

```r
ecls %>%
  group_by(catholic) %>%
  summarise(
    n_students = n(),
    mean_math = mean(c5r2mtsc_std),
    std_error = sd(c5r2mtsc_std) / sqrt(n_students)
  )
# # A tibble: 2 x 4
#   catholic n_students mean_math std_error
#      <int>      <int>     <dbl>     <dbl>
# 1        0       4597     0.156    0.0144
# 2        1        951     0.221    0.0277
```

```r
# 连续型混杂因素的组间差异
ecls %>%
  group_by(catholic) %>%
  select(p5hmage, w3momscr, w3dadscr) %>%
  summarise_all(list(~mean(., na.rm = T)))
# # A tibble: 2 x 4
#   catholic p5hmage w3momscr w3dadscr
#      <int>   <dbl>    <dbl>    <dbl>
# 1        0    37.8     43.8     42.6
# 2        1    39.8     47.5     45.8
```

```r
# 分类型混杂因素的组间差异
tab <- xtabs(~race_white + catholic, data = ecls)
chisq.test(tab, correct = F)
# X-squared = 48.596, df = 1, p-value = 3.145e-12
```

### 步骤 2：计算倾向性评分

以处理因素为因变量，所有混杂因素为自变量，构建逻辑回归模型。模型的预测概率即为倾向性评分。

```r
m_ps <- glm(catholic ~ race_white + w3momed_hsb + p5hmage + w3momscr + w3dadscr,
            family = binomial(), data = ecls)

pr_score <- predict(m_ps, type = "response")

ecls_ps <- ecls %>%
  mutate(ps = pr_score)
```

可查看 PS 在两组间的分布：

```r
prs_df <- data.frame(pr_score = pr_score, catholic = m_ps$model$catholic)
labs <- paste("Actual school type attended:", c("Catholic", "Public"))
prs_df %>%
  mutate(catholic = ifelse(catholic == 1, labs[1], labs[2])) %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~catholic) +
  xlab("Probability of going to Catholic school") +
  theme_bw()
```

PS 分布偏态时，可对 PS 做 log 变换后再进行后续分析。

### 步骤 3：倾向性评分回归

将 PS 和处理因素一同纳入结局变量的线性回归（或逻辑回归、Cox 回归等，视结局类型而定）：

```r
psl <- lm(c5r2mtsc_std ~ catholic + ps, data = ecls_ps)
summary(psl)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.58249    0.02929 -19.885  < 2e-16 ***
# catholic    -0.10772    0.03241  -3.324 0.000893 ***
# ps           4.48236    0.15873  28.239  < 2e-16 ***
# ---
# Residual standard error: 0.8934 on 5545 degrees of freedom
# Multiple R-squared:  0.1263, Adjusted R-squared:  0.126
```

关键结果：`catholic` 的系数 -0.10772（p = 0.000893），校正 PS 后处理因素仍有统计学意义。

### 步骤 4：倾向性评分分层

首先检查两组 PS 范围是否重叠：

```r
ecls_ps %>% group_by(catholic) %>%
  summarise(range = range(ps))
# # A tibble: 4 x 2
# # Groups:   catholic [2]
#   catholic  range
#      <int>  <dbl>
# 1        0 0.0370
# 2        0 0.477
# 3        1 0.0492
# 4        1 0.404
```

两组范围基本一致（0.037~0.477 vs 0.049~0.404），可以分层。若两组 PS 范围差异大，应按交集部分分层。

一般按 PS 分 5~10 层，可用等分法或百分位数法。示例按 0.1、0.2、0.3 为切点分 4 层：

```r
ecls_pslevel <- ecls_ps %>%
  mutate(
    ps_level = case_when(
      ps <= 0.1           ~ "level_1",
      ps > 0.1 & ps <= 0.2 ~ "level_2",
      ps > 0.2 & ps <= 0.3 ~ "level_3",
      TRUE                ~ "level_4"
    ),
    p5hmage = as.double(p5hmage),
    across(where(is.integer), as.factor)
  )
```

### 步骤 5：分层后均衡性检验与效应估计

在各层内检验混杂因素的均衡性和结局变量的差异。均衡性好的层：混杂因素在两组间均无显著差异。效应一致的层：结局变量在两组间均存在显著差异。

**连续型变量**（t 检验）：

```r
ecls_pslevel %>%
  pivot_longer(cols = c(1, 5:7), names_to = "variates", values_to = "values") %>%
  group_nest(ps_level, variates) %>%
  dplyr::mutate(
    tt  = map(data, ~ t.test(values ~ catholic, data = .x)),
    res = map_dfr(tt, broom::tidy)
  ) %>%
  unnest(res)
```

**分类型变量**（卡方检验）：

```r
ecls_pslevel %>%
  group_split(ps_level) %>%
  map(~ chisq.test(.$race_white, .$catholic, correct = F)) %>%
  map_dbl("p.value")
# [1] 0.4755703 0.8423902 0.5696924 0.2667193
```

### 步骤 6：结果报告

论文中可这样描述：

> 采用倾向性评分回归调整混杂因素，将倾向性评分作为协变量纳入多元线性回归模型。结果表明，校正混杂后两组学生成绩差异具有统计学意义（β = -0.108，p = 0.001）。进一步按倾向性评分四分位数分层分析，各层内混杂因素分布均衡（P 均 > 0.05）。

## 代码示例

```r
# 完整工作流：倾向性评分回归 + 分层
library(tidyverse)
library(broom)

# 1. 数据加载
ecls <- read.csv("datasets/ecls.csv") %>%
  dplyr::select(c5r2mtsc_std, catholic, race_white, w3momed_hsb,
                p5hmage, w3momscr, w3dadscr) %>%
  na.omit()

# 2. 基线查看
ecls %>%
  group_by(catholic) %>%
  summarise(
    n = n(),
    mean_outcome = mean(c5r2mtsc_std),
    se = sd(c5r2mtsc_std) / sqrt(n)
  )
# # A tibble: 2 x 4
#   catholic     n mean_outcome     se
#      <int> <int>        <dbl>  <dbl>
# 1        0  4597        0.156 0.0144
# 2        1   951        0.221 0.0277

with(ecls, t.test(c5r2mtsc_std ~ catholic))
# Welch Two Sample t-test
# t = -2.0757, p-value = 0.03809

# 3. 计算倾向性评分
m_ps <- glm(catholic ~ race_white + w3momed_hsb + p5hmage + w3momscr + w3dadscr,
            family = binomial(), data = ecls)

ecls_ps <- ecls %>%
  mutate(ps = predict(m_ps, type = "response"))

# 4. 倾向性评分回归
model_ps <- lm(c5r2mtsc_std ~ catholic + ps, data = ecls_ps)
summary(model_ps)

# 5. 倾向性评分分层
ecls_pslevel <- ecls_ps %>%
  mutate(
    ps_level = case_when(
      ps <= 0.1            ~ "Q1",
      ps > 0.1 & ps <= 0.2 ~ "Q2",
      ps > 0.2 & ps <= 0.3 ~ "Q3",
      TRUE                 ~ "Q4"
    ),
    across(where(is.integer), as.factor)
  )

# 6. 分层后均衡性检验 —— 连续型变量
ecls_pslevel %>%
  pivot_longer(cols = c(c5r2mtsc_std, p5hmage, w3momscr, w3dadscr),
               names_to = "variate", values_to = "value") %>%
  group_nest(ps_level, variate) %>%
  mutate(
    ttest = map(data, ~ t.test(value ~ catholic, data = .x)),
    res   = map_dfr(ttest, broom::tidy)
  ) %>%
  unnest(res) %>%
  dplyr::select(ps_level, variate, estimate1, estimate2, p.value)

# 7. 分层后均衡性检验 —— 分类变量
ecls_pslevel %>%
  group_split(ps_level) %>%
  map(~ chisq.test(.$race_white, .$catholic, correct = F)) %>%
  map_dbl("p.value")
```

## 结果解读指南

**倾向性评分回归输出解读**：

- `catholic` 的 Estimate = -0.108：在控制倾向性评分后，天主教学校学生的标准化成绩比公立学校学生平均低 0.108 个标准差
- p = 0.000893 < 0.05：学校类型的效应有统计学意义
- `ps` 的 Estimate = 4.48：PS 对成绩有显著预测作用，验证了纳入 PS 的合理性
- 注意：本示例中符号由正变负，说明未校正前可能因混杂因素（如家庭背景更好）掩盖了真实效应方向

**倾向性评分分层输出解读**：

- 查看输出中 `p.value` 列：若各层内混杂因素的 p > 0.05，说明分层后该混杂因素在两组间均衡
- 理想模式：混杂因素 p 均 > 0.05（均衡），结局变量 p 均 < 0.05（有效应）
- 若分层后结局变量不再显著，可能的原因：分层过细导致每层样本量不足；PS 并未充分捕获所有混杂信息；数据本身不适合分层分析

**与 PSM 对比**：

| 特征 | PS 回归 | PS 分层 | PSM |
|------|---------|---------|-----|
| 样本保留 | 全部 | 全部 | 匹配后减少 |
| 操作复杂度 | 最低 | 中等 | 较高 |
| 混杂控制 | 依赖模型假设 | 依赖分层合理性 | 通过匹配直观平衡 |
| 常用场景 | 快速校正 | 敏感性分析 | 基线严重不均衡 |

## 常见问题与注意事项

**Q1：为什么分层后混杂因素还是不平衡？**

分层效果取决于 PS 模型的质量和分层策略。可能原因：PS 模型遗漏了重要交互项或非线性项；分层切点不合理；样本量不足以支撑高层数。可尝试增加 PS 模型的复杂度（加入交互项、平方项）、调整分层切点（如用五分位数而非等距切点）。

**Q2：应该分几层合适？**

文献一般建议 5~10 层。层数太少则层内仍有残余混杂；层数太多则每层样本量过少、估计不稳定。可用敏感性分析：分别分 5 层、10 层看结论是否一致。

**Q3：PS 回归和 PS 分层哪个更好？**

PS 回归最简单，但依赖模型线性假设。PS 分层对模型假设要求更弱，但分层策略的选择有一定主观性。实践中两者常并用：PS 回归作为主要分析，PS 分层作为敏感性分析。如果两者结论一致，结果更可信。

**Q4：PS 分布偏态怎么办？**

可对 PS 做 logit 变换（`log(ps/(1-ps))`），然后用变换后的值进行回归或分层。logit 变换后的 PS 通常更接近正态分布，有利于线性模型假设。代码示例：

```r
ecls_ps <- ecls_ps %>%
  mutate(ps_logit = log(ps / (1 - ps)))
```

**Q5：这个例子中分层结果不理想，说明什么？**

原示例中分层后仅有 level_3 层内结局变量有差异，其余层均无差异，提示该数据可能不适合分层分析。实践中应尝试不同方法（PS 回归、PSM、IPTW）并比较结果，选择最适合当前数据的方法。

**Q6：SPSS 用户迁移到 R 的注意事项**：

SPSS 中倾向性评分的实现需要使用 `PROPENSITY SCORE` 插件或手动逻辑回归取预测概率。R 中直接用 `glm()` + `predict(type = "response")` 即可，更灵活。分层后的循环检验在 SPSS 中需要重复操作或用 macro，R 中用 `purrr::map` + `group_split` 可一行代码完成。
