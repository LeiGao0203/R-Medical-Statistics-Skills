---
name: medical-stat-subgroup-analysis
description: "R语言医学统计：亚组分析及森林图绘制。在回归模型中按不同亚组（如性别、年龄层）分层展示效应量，使用forestploter包绘制亚组分析和多因素回归模型的森林图。TRIGGER when user mentions 亚组分析、森林图、森林图绘制、分层效应、交互作用图、subgroup analysis、forest plot，or asks about stratified effect estimates. SKIP for Meta分析的森林图、简单回归结果展示、单因素生存曲线图。"
---

# 亚组分析及森林图绘制 (Subgroup Analysis and Forest Plot)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用亚组分析的典型场景：**

- 在回归模型（Cox 回归、Logistic 回归、线性回归等）中，按不同分类变量（如性别、年龄段、肿瘤分期）分层展示干预/暴露的效应量
- 探索某一治疗方法在不同特征人群间的效应是否存在差异（效应修饰 / 交互作用初筛）
- 控制混杂因素的方法之一：通过分层分析观察各亚组效应量与总效应量的差异（辛普森悖论问题）
- 高分 SCI 论文常见：用森林图直观展示每个亚组的 HR/OR 及其 95% CI、p 值、各组样本量

**不使用亚组分析 / 森林图的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 仅展示单一模型的回归系数汇总 | `broom::tidy()` 表格即可，无需亚组分析 |
| 数据来自多个研究的汇总（Meta 分析） | 使用 `meta` / `metafor` 包做 Meta 分析森林图 |
| 需要统计检验交互作用是否显著 | 在模型中添加交互项（`treatment * subgroup`），用似然比检验或 Wald 检验 |
| 仅展示单因素模型各变量的 HR | 单因素多变量表格即可，无需分层 |
| 仅展示 K-M 生存曲线 | 使用 `ggsurvplot()` 或 `survminer::ggsurvplot()` |

**医学研究常见应用：**

- 新药临床试验中，按性别、年龄（>65 / ≤65）、BMI 分层展示主要终点
- 观察性研究中，按合并症（有 / 无）分层展示暴露与结局的关联
- 肿瘤预后研究中，按肿瘤分期、淋巴结转移状态分层展示治疗效应

## 前置条件

**R 包安装：**

```r
install.packages(c("survival", "tidyverse", "broom", "forestploter"))
```

**核心数据要求：**

- 生存数据（Cox 回归亚组分析）：`time`（生存时间，≥0）、`status`（结局事件：0 = 删失，1 = 事件）
- 暴露/干预变量：需为二分类因子型，表示两组对照关系（如治疗组 vs 对照组）
- 亚组变量：需为因子型（factor），每个取值代表不同的亚组水平（如 "male" / "female"）
- 其他协变量（选填）：多因素回归亚组分析时使用

**关键注意：**
- 分类变量必须转为因子型，这样 R 回归函数会自动进行哑变量编码
- 连续型变量若要作为亚组变量，需先按临界值做二分类（如年龄 > 65 / ≤ 65）
- 亚组变量的缺失值需事先清理（`drop_na()`），否则会导致某一层结果缺失

## 方法选择决策树

```
你的研究目标 →
├── 仅需快速实现分层效应量计算
│   ├── 仅一个亚组变量、层数少 → 手动子集拟合（逐步在每个亚组内单独回归）
│   └── 多个亚组变量、层数多 → tidyverse + purrr 批量亚组分析
│
├── 需要绘制森林图
│   ├── 简单快速 → forestploter::forest() 默认参数
│   ├── 需要 NEJM 风格美化 → forestploter::forest() + forest_theme() 自定义主题
│   └── 需要用 ggplot2 生态 → ggforestplot / ggplot2 自建森林图
│
├── 数据为 Cox 回归 → coxph(Surv(time, status) ~ treatment, data = subset)
├── 数据为 Logistic 回归 → glm(outcome ~ treatment, family = binomial, data = subset)
├── 数据为线性回归 → lm(outcome ~ treatment, data = subset)
│
└── 需要展示不同模型/协变量调整程度的森林图
    └── 分层构建多个模型结果，合并数据后统一用 forestploter 绘制
```

## 标准工作流

### 步骤 1：数据准备与探索

使用 `survival` 包中的 `colon` 数据集演示。该数据集包含 1858 行结肠癌患者生存数据，16 列变量。

```r
rm(list = ls())
library(survival)

str(colon)  ## 'data.frame': 1858 obs. of 16 variables

# 数据预处理：筛选 Obs 组和 Lev+5FU 组，变量转为因子
suppressMessages(library(tidyverse))

df <- colon %>%
  mutate(rx = as.numeric(rx)) %>%
  filter(etype == 1, !rx == 2) %>%
  select(time, status, rx, sex, age, obstruct, perfor, adhere,
         differ, extent, surg, node4) %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), labels = c("female", "male")),
    age = ifelse(age > 65, ">65", "<=65"),
    age = factor(age, levels = c(">65", "<=65")),
    obstruct = factor(obstruct, levels = c(0, 1), labels = c("No", "Yes")),
    perfor   = factor(perfor,   levels = c(0, 1), labels = c("No", "Yes")),
    adhere   = factor(adhere,   levels = c(0, 1), labels = c("No", "Yes")),
    differ   = factor(differ,   levels = c(1, 2, 3),
                      labels = c("well", "moderate", "poor")),
    extent   = factor(extent,   levels = c(1, 2, 3, 4),
                      labels = c("submucosa", "muscle", "serosa", "contiguous")),
    surg     = factor(surg,     levels = c(0, 1), labels = c("short", "long")),
    node4    = factor(node4,    levels = c(0, 1), labels = c("No", "Yes")),
    rx       = ifelse(rx == 3, 0, 1),
    rx       = factor(rx, levels = c(0, 1))
  )

str(df)  ## 'data.frame': 619 obs. of 12 variables
```

### 步骤 2：不分亚组的整体分析（总效应量）

先对所有数据进行单因素 Cox 回归，获得总体效应量，作为后续亚组分析的参照基准：

```r
fit <- coxph(Surv(time, status) ~ rx, data = df)
broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
## # A tibble: 1 × 7 — estimate = 1.67, p < 0.001, 95% CI: 1.32–2.11
```

结果解读：整体人群中 Lev+5FU 组的复发风险优于 Obs 组（HR = 1.67, p < 0.001）。但需注意**辛普森悖论**：总体有效的方法在不同亚组（如单独男性或女性）中可能表现为无效，这正是亚组分析的意义。

### 步骤 3：手动单亚组分析

以 `sex` 为亚组变量，分别在男性、女性中拟合 Cox 回归：

```r
# 男性亚组
fit_male <- coxph(Surv(time, status) ~ rx, data = df[df$sex == "male", ])
broom::tidy(fit_male, exponentiate = TRUE, conf.int = TRUE)
## HR = 2.29, 95% CI: 1.60–3.26, p < 0.001

# 女性亚组
fit_female <- coxph(Surv(time, status) ~ rx, data = df[df$sex == "female", ])
broom::tidy(fit_female, exponentiate = TRUE, conf.int = TRUE)
## HR = 1.32, 95% CI: 0.96–1.80, p = 0.088
```

结果发现：男性亚组 HR = 2.29（p < 0.001）治疗效应显著，女性亚组 HR = 1.32（p = 0.088）治疗效应不显著。两个亚组的效应方向一致但强度不同，提示存在性别交互作用的可能。手动方式适合亚组数量少的情况，亚组多时应使用批量分析。

### 步骤 4：tidyverse + purrr 批量亚组分析

核心思路：将宽数据转换为长数据 → 按亚组变量和水平分组 → 在每个组内拟合模型 → 提取结果。

```r
# 宽转长
dfl <- df %>%
  pivot_longer(cols = 4:ncol(.), names_to = "var", values_to = "value") %>%
  arrange(var)

# 批量对每个亚组拟合 Cox 回归
ress <- dfl %>%
  group_nest(var, value) %>%
  drop_na(value) %>%
  mutate(
    model = map(data, ~ coxph(Surv(time, status) ~ rx, data = .x)),
    res   = map(model, broom::tidy, conf.int = TRUE, exponentiate = TRUE)
  ) %>%
  select(var, value, res)

# 计算每个亚组中两种治疗方式的人数
ss <- dfl %>%
  group_by(var, value, rx) %>%
  drop_na(value) %>%
  summarise(sample_size = n(), .groups = "drop") %>%
  select(var, value, rx, sample_size)

# 合并亚组结果与样本量
resss <- ress %>%
  left_join(ss, by = c("var", "value")) %>%
  unnest(c(res, rx, sample_size)) %>%
  pivot_wider(names_from = "rx", values_from = "sample_size",
              names_prefix = "rx_") %>%
  select(-c(term, std.error, statistic)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>%
  mutate(`HR(95%CI)` = paste0(estimate, "(", conf.low, "-", conf.high, ")"))

```

### 步骤 5：合并总体效应量与亚组效应量

```r
# 提取总体效应量
fit <- coxph(Surv(time, status) ~ rx, data = df)
res_all <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

df %>% count(rx)

res_all <- res_all %>%
  mutate(
    var = "All people", value = " ",
    rx_0 = 304, rx_1 = 315,
    across(where(is.numeric), ~ round(.x, digits = 2))
  ) %>%
  mutate(`HR(95%CI)` = paste0(estimate, "(", conf.low, "-", conf.high, ")")) %>%
  select(var, value, estimate, p.value, conf.low, conf.high, rx_0, rx_1, `HR(95%CI)`)

# 合并所有结果
resss <- bind_rows(res_all, resss)
```

### 步骤 6：森林图绘制

使用 `forestploter` 包绘制 NEJM 风格的森林图。需要先准备一个包含分组标签和空占位行的数据框：

```r
library(forestploter)
library(grid)

# 读入整理好格式的数据（可直接从 R 对象构建，也可导出 CSV 在 Excel 中调整后读回）
# 核心：subgroup 列为标签列，estimate / conf.low / conf.high 为效应量和 CI，
# HR(95%CI) 为文字列，p.value 为 p 值列，rx_0 / rx_1 为两个治疗组的人数列
# 亚组变量名的行（如 "adhere"、"age"）对应的 estimate / CI / p 值均留空（NA）

# 数据准备（将分组标签行的 NA 替换为空格以隐藏图形成分）
plot_df <- read.csv("datasets/resss.csv", check.names = FALSE)
# 格式：Subgroup, estimate, p.value, conf.low, conf.high, rx_0, rx_1, HR(95%CI)
# 亚组标签行（如 "adhere", "age"）的 estimate 等列为 NA

plot_df[, c(3, 6, 7)][is.na(plot_df[, c(3, 6, 7)])] <- " "
plot_df$` ` <- paste(rep(" ", nrow(plot_df)), collapse = " ")

# 基础森林图
p <- forest(
  data     = plot_df[, c(1, 6, 7, 9, 8, 3)],
  lower    = plot_df$conf.low,
  upper    = plot_df$conf.high,
  est      = plot_df$estimate,
  ci_column = 4,
  sizes    = (plot_df$estimate + 0.001) * 0.3,
  ref_line = 1,
  xlim     = c(0.1, 4)
)
print(p)
```

### 步骤 7：森林图美化（NEJM 风格）

```r
# 自定义主题
tm <- forest_theme(
  base_size        = 12,
  ci_lwd           = 1.5,
  refline_lwd      = 1.5,
  refline_lty      = "dashed",
  refline_col      = "grey20",
  footnote_cex     = 0.8,
  footnote_fontface = "italic",
  footnote_col     = "grey30",
  core = list(bg_params = list(fill = c("#FFFFFF", "#f5f7f6"), col = NA))
)

p <- forest(
  data      = plot_df[, c(1, 6, 7, 9, 8, 3)],
  lower     = plot_df$conf.low,
  upper     = plot_df$conf.high,
  est       = plot_df$estimate,
  ci_column = 4,
  sizes     = (plot_df$estimate + 0.001) * 0.3,
  ref_line  = 1,
  xlim      = c(0.1, 4),
  arrow_lab = c("Obs better", "Lev+5-FU better"),
  theme     = tm
)
print(p)
```

## 代码示例

### 示例 1：Cox 回归亚组分析精简版

```r
rm(list = ls())
library(survival)
suppressMessages(library(tidyverse))
library(broom)

# 数据准备（同标准工作流步骤1）
df <- colon %>%
  mutate(rx = as.numeric(rx)) %>%
  filter(etype == 1, !rx == 2) %>%
  select(time, status, rx, sex, age, obstruct, perfor, adhere,
         differ, extent, surg, node4) %>%
  mutate(across(c(sex, obstruct, perfor, adhere, differ, extent, surg, node4), factor),
         age = factor(ifelse(age > 65, ">65", "<=65"), levels = c(">65", "<=65")),
         rx  = factor(ifelse(rx == 3, 0, 1), levels = c(0, 1)))

# 批量亚组分析
dfl <- df %>%
  pivot_longer(cols = 4:ncol(.), names_to = "var", values_to = "value") %>%
  arrange(var)

ress <- dfl %>%
  group_nest(var, value) %>%
  drop_na(value) %>%
  mutate(
    res = map(data, ~ broom::tidy(
      coxph(Surv(time, status) ~ rx, data = .x),
      conf.int = TRUE, exponentiate = TRUE))
  ) %>%
  select(var, value, res) %>%
  unnest(res) %>%
  select(var, value, estimate, p.value, conf.low, conf.high)
```

### 示例 2：Logistic 回归 / 线性回归亚组分析

```r
# Logistic 回归：将 coxph 替换为 glm
glm_result <- dfl %>%
  group_nest(var, value) %>%
  drop_na(value) %>%
  mutate(
    res = map(data, ~ broom::tidy(
      glm(outcome ~ exposure, data = .x, family = binomial),
      conf.int = TRUE, exponentiate = TRUE))
  )

# 线性回归：连续结局的亚组分析
lm_result <- dfl %>%
  group_nest(var, value) %>%
  drop_na(value) %>%
  mutate(
    res = map(data, ~ broom::tidy(
      lm(outcome ~ exposure, data = .x),
      conf.int = TRUE))
  )
```

### 示例 3：森林图完整绘制

```r
library(forestploter)

# 假设 plot_df 已经整理好（Subgroup, estimate, p.value, conf.low, conf.high,
#                                    rx_0, rx_1, `HR(95%CI)`）
# 注意：亚组变量名行（如 "sex"、"age"）的 estimate = NA

plot_df[, c("estimate", "p.value", "conf.low", "conf.high")][
  is.na(plot_df[, c("estimate", "p.value", "conf.low", "conf.high")])] <- " "
plot_df$` ` <- paste(rep(" ", nrow(plot_df)), collapse = " ")

tm <- forest_theme(
  base_size = 12, ci_lwd = 1.5,
  refline_lwd = 1.5, refline_lty = "dashed", refline_col = "grey20",
  core = list(bg_params = list(fill = c("#FFFFFF", "#f5f7f6"), col = NA))
)

p <- forest(
  data      = plot_df[, c("Subgroup", "rx_0", "rx_1", "HR(95%CI)", " ", "p.value")],
  lower     = as.numeric(plot_df$conf.low),
  upper     = as.numeric(plot_df$conf.high),
  est       = as.numeric(plot_df$estimate),
  ci_column = 4,
  sizes     = 0.3,
  ref_line  = 1,
  xlim      = c(0.1, 4),
  arrow_lab = c("Favors Control", "Favors Treatment"),
  theme     = tm
)
print(p)
```

## 结果解读指南

### 亚组分析结果解读

| 输出项 | 说明 |
|--------|------|
| `estimate` | 效应量（HR / OR / β 系数），取决于模型类型 |
| `p.value` | Wald 检验 p 值，p < 0.05 表示该亚组内暴露与结局存在统计学关联 |
| `conf.low / conf.high` | 效应量的 95% 置信区间上下界 |
| `rx_0 / rx_1` | 该亚组中对照组和治疗组的样本量 |
| `HR(95%CI)` | 格式化的文字列，森林图表格中直接展示用 |

**解读要点：**

1. **总效应量**（All people 行）：反映在所有患者中的暴露效应
2. **亚组间比较**：比较各亚组的 HR 方向（同向 / 反向）和强度，但不能仅凭 p 值判断差异是否显著 —— p 值仅反映该亚组内部效应是否 ≠ 1，不反映亚组间差异是否显著
3. **交互作用**：要正式检验亚组间差异，需在模型中引入交互项（`treatment * subgroup`），而非仅比对各亚组的 p 值
4. **可信区间宽度**：小样本亚组的 CI 宽、精度低，结论需谨慎

**论文报告示例**：Lev+5FU 治疗组在总体人群中降低了结肠癌复发风险（HR = 1.67, 95% CI: 1.32–2.11, p < 0.001）。按性别分层后，男性患者治疗获益显著（HR = 2.29, 95% CI: 1.60–3.26, p < 0.001），女性患者治疗获益无统计学显著性（HR = 1.32, 95% CI: 0.96–1.80, p = 0.088）。

### 森林图解读

- **横轴**：效应量（HR / OR / β）的值，ref_line = 1 表示无效应
- **方块**：点估计（方块大小 ≈ 样本量的平方根或效应量的绝对值）
- **水平线**：95% 置信区间，线越长精度越低
- **虚线**：ref_line，CI 跨过此线表示 p > 0.05（无统计学意义）
- **箭头标记**：xlim 范围外区域的文字标签，如 "Obs better" / "Lev+5-FU better"

## 常见问题与注意事项

### Q1：亚组分析与交互作用检验有什么区别？

亚组分析是**描述性**的，仅展示不同亚组层内的效应量，不能断言"效应量因亚组不同而有差异"。交互作用检验是**推断性**的 —— 在模型中加入 `treatment * subgroup` 交互项，用似然比检验或 Wald 检验判断交互项系数是否 ≠ 0。论文中应同时报告亚组分析森林图（可视化）和交互作用 p 值（统计学检验）。

### Q2：多分类亚组变量如何处理？

`group_nest(var, value)` 会自动将每个水平作为独立亚组。例如 `differ` 有 well / moderate / poor 三个水平，会分别生成 3 个结果行。在森林图中，三个子行共用父行标签 "differ"。需注意当多分类时，参考水平是谁取决于因子编码顺序。

### Q3：如何把 p = 0.000 的结果在森林图中显示为 "< 0.001"？

在数据准备阶段，将 p.value 列转为字符型再替换：`p_display = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))`，森林图表格列使用 `p_display` 而非原始 `p.value`。

### Q4：连续型变量可以做亚组分析吗？

需先二分类化。按中位数或临床临界值拆分（如 BMI ≥ 25 / < 25、年龄 > 65 / ≤ 65）。不建议随意切分（会损失信息），更好的做法是使用交互项或样条回归探索非线性交互效应。

### Q5：亚组分析中需要校正协变量吗？

协变量调整的亚组分析中，在 `coxph()` 公式中加入协变量：`coxph(Surv(time, status) ~ rx + age + sex, data = .x)`。区分：不校正展示原始效应 vs 全模型校正展示独立效应。森林图中应标注是否校正及校正了哪些变量。

### Q6：为什么 forestploter 的 estimate 列不能有 0 或 NA？

`sizes` 参数使用 `(estimate + 0.001) * 0.3` 计算方块大小，不能为 0 或负值，否则绘图时方块消失。对于亚组标签行（estimate 为 NA），将对应行的 estimate 列设为空格字符串 " " 可隐藏该行的图形元素。

### Q7：`Publish` 包的亚组分析有什么问题？

`Publish::subgroupAnalysis()` 在进行亚组分析时拟合的是错误的模型（未正确处理各亚组模型中的参考水平），可能导致结果不准确。建议使用本技能中展示的 tidyverse + purrr 方案，逐层拟合独立的单因素回归，确保每个亚组的模型使用正确的子集数据和唯一的参考水平。

### Q8：SPSS 与 R 的亚组分析有何差异？

SPSS 中可通过"拆分文件"（Split File）实现分层分析，但每次只能按一个变量拆分，需手动记录输出结果。R 的 tidyverse + purrr 可一次性批量处理所有亚组变量并输出结构化结果。SPSS 的森林图需额外模块，R 的 forestploter 直接输出出版级图形。

### Q9：forestploter 如何导出图形？

```r
pdf("forest_plot.pdf", width = 8, height = 10)
print(p)
dev.off()

png("forest_plot.png", width = 8, height = 10, units = "in", res = 300)
print(p)
dev.off()
```

### 关键提醒

- 亚组分析是探索性的，不能替代交互作用检验 —— 论文中应同时报告两者的结果
- `group_nest()` 分组后每个水平样本量独立，小样本亚组结果不稳定，解释需谨慎
- 批量拟合模型时出现 "Loglik converged before variable" 警告是正常的，表示该亚组数据完美预测（分离），可跳过或标记该结果
- forestploter 的 `xlim` 需足够宽以容纳所有 CI，超出 `xlim` 范围的 CI 会自动显示为箭头
- 分组的顺序影响森林图的外观，建议在 Excel 中手动排列好顺序再读入 R 绘图
