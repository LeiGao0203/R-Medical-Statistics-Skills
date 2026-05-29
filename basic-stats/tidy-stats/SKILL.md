---
name: medical-stat-tidy-stats
description: "R语言医学统计：tidy风格批量统计分析。使用rstatix包进行管道友好的批量t检验、Wilcoxon检验、方差分析、Kruskal-Wallis检验、卡方检验和Fisher检验，自动输出整洁数据框格式结果。TRIGGER when user mentions 批量检验、同时对多个变量做检验、tidy统计、rstatix、管道操作统计，or asks about running multiple statistical tests efficiently across many variables. SKIP for 单个统计检验、复杂模型构建。"
---

# tidy流统计分析 (Tidy-style Statistical Analysis)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用 tidy 流统计分析：**

- 需要同时对多个结局变量执行相同类型的统计检验（如同时对4项量表评分做t检验）
- 需要批量进行正态性检验、方差齐性检验等前提条件检查
- 希望结果以整洁数据框（tibble）形式输出，方便后续可视化或导出
- 使用管道操作（`%>%`）串联数据处理和统计分析全流程
- SPSS用户迁移到R，希望实现"一键多变量同时检验"的操作体验
- 临床试验基线表多变量组间比较（Table 1）

**何时不使用 tidy 流统计分析：**

| 情况 | 替代方法 |
|------|----------|
| 单个变量的假设检验 | 直接使用 `t.test()` / `wilcox.test()` / `chisq.test()` 等基础R函数 |
| 复杂模型（多元回归、生存分析等） | 直接使用 `lm()` / `glm()` / `coxph()` 等建模函数 |
| 需要自定义检验逻辑（条件分支） | 结合 `purrr::map()` 或 for 循环自行编写 |
| 计量/纯统计学习题 | 使用 `pwr` 包或 `stats` 基础函数手动计算 |

## 前置条件

**R 包依赖：**

```r
install.packages("rstatix")
install.packages("tidyverse")   # 数据整形与管道
```

**rstatix 核心函数与 base R 等价对照：**

| rstatix 函数 | Base R 等价 | 适用场景 |
|---|---|---|
| `shapiro_test()` | `shapiro.test()` | 正态性检验 |
| `levene_test()` | `car::leveneTest()` | 方差齐性检验 |
| `t_test()` | `t.test()` | 两组均值比较（自动 Welch 校正） |
| `wilcox_test()` | `wilcox.test()` | 两/多组秩和检验 |
| `anova_test()` | `aov()` + `summary()` | 方差分析（含效应量 ges） |
| `kruskal_test()` | `kruskal.test()` | Kruskal-Wallis 秩和检验 |
| `chisq_test()` | `chisq.test()` | 分类变量关联性检验 |
| `fisher_test()` | `fisher.test()` | Fisher 精确检验 |

**数据格式要求：**

- **必须使用长数据格式**（每行一个观测，含分组变量列、变量名列、数值列）
- 宽数据需先用 `pivot_longer()` 转换
- 分组变量应为 factor 类型
- 所有 rstatix 函数均兼容 `%>%` 管道和 `group_by()` 分组

**统计假设：**

- t检验：正态性、方差齐性（Welch t检验对方差不齐稳健）
- ANOVA：正态性、方差齐性（`levene_test()`）
- 卡方检验：期望频数不宜过小（<5的格数不超过20%）

## 方法选择决策树

```
你的数据分析任务 →
├── 需同时对多个连续变量做检验
│   ├── 两组比较
│   │   ├── 数据满足正态性 → rstatix::t_test()
│   │   └── 数据不满足正态性 → rstatix::wilcox_test()
│   ├── 多组（≥3）比较
│   │   ├── 数据满足正态性 + 方差齐性 → rstatix::anova_test()
│   │   └── 数据不满足条件 → rstatix::kruskal_test()
│   └── 需先检验前提条件
│       ├── 正态性 → group_by + shapiro_test()
│       └── 方差齐性 → group_by + levene_test()
│
├── 需同时对多个分类变量做检验
│   ├── 样本量大、期望频数充足 → rstatix::chisq_test()
│   └── 小样本或期望频数<5 → rstatix::fisher_test()
│
└── 仅单个变量检验 → 使用 base R 函数（t.test / wilcox.test / chisq.test）
```

## 标准工作流

### 步骤1：数据准备——宽数据转长数据

```r
library(tidyverse)
df <- read.csv("datasets/20210801.csv", header = TRUE)
df_l <- df %>%
  pivot_longer(
    cols = 2:5,           # 待检验的变量列
    names_to = "变量",     # 变量名列
    values_to = "积分"     # 数值列
  ) %>%
  dplyr::mutate_if(is.character, as.factor)
```

将宽格式的"一个变量一列"转为长格式的"变量名-值"键值对结构，这是 rstatix 管道批量操作的核心前提。

### 步骤2：批量正态性检验

按变量和分组进行 Shapiro-Wilk 检验：

```r
library(rstatix)
df_l %>% group_by(变量, 组别) %>% shapiro_test(积分)
```

输出为整洁数据框，每行对应一个变量-组别的正态性检验结果（statistic + p值）。

### 步骤3：批量方差齐性检验

按变量进行 Levene 检验（仅按变量分组，不按组别）：

```r
df_l %>% group_by(变量) %>% levene_test(积分 ~ 组别)
```

### 步骤4：执行批量统计分析

根据前提条件检验结果，选择合适检验：

```r
# 批量t检验（Welch校正，默认不等方差）
df_l %>% group_by(变量) %>% t_test(积分 ~ 组别)

# 批量Wilcoxon秩和检验
df_l %>% group_by(变量) %>% wilcox_test(积分 ~ 组别)

# 批量方差分析
df_l %>% group_by(变量) %>% anova_test(积分 ~ 组别)

# 批量Kruskal-Wallis检验
df_l %>% group_by(变量) %>% kruskal_test(积分 ~ 组别)
```

结果直接为 tibble 格式，可用 `filter(p < 0.05)` 快速筛选有统计学意义的变量。

### 步骤5：结果报告

论文中可表述为："对4项评分指标分别进行独立样本t检验，结果显示两组在排便困难（t=0.089, p=0.929）、生活质量（t=-0.101, p=0.920）、粪便性状（t=-0.356, p=0.723）和排便时间（t=-0.273, p=0.786）方面的差异均无统计学意义。"

## 代码示例

### 完整示例：多变量批量t检验

```r
library(tidyverse)
library(rstatix)

# 1. 读入宽数据
df <- read.csv("datasets/20210801.csv", header = TRUE)
str(df)
## 'data.frame':    60 obs. of  5 variables:
##  $ 组别    : chr  "实验组" "实验组" "实验组" "实验组" ...
##  $ 排便困难: int  12 11 15 14 11 13 11 14 13 15 ...
##  $ 生活质量: int  87 94 95 85 101 91 84 89 84 92 ...
##  $ 粪便性状: int  2 3 2 2 3 3 3 3 3 3 ...
##  $ 排便时间: int  2 2 2 2 2 2 2 3 2 2 ...

# 2. 宽→长转换
df_l <- df %>%
  pivot_longer(cols = 2:5, names_to = "变量", values_to = "积分") %>%
  dplyr::mutate_if(is.character, as.factor)

head(df_l)
## # A tibble: 6 × 3
##   组别   变量      积分
##   <fct>  <fct>    <int>
## 1 实验组 排便困难    12
## 2 实验组 生活质量    87
## 3 实验组 粪便性状     2
## 4 实验组 排便时间     2
## 5 实验组 排便困难    11
## 6 实验组 生活质量    94

# 3. 同时正态性检验
df_l %>% group_by(变量, 组别) %>% shapiro_test(积分)
## # A tibble: 8 × 5
##   组别   变量     variable statistic        p
##   <fct>  <fct>    <chr>        <dbl>    <dbl>
## 1 对照组 粪便性状 积分         0.452 1.73e- 9
## 2 实验组 粪便性状 积分         0.404 5.98e-10
## 3 对照组 排便困难 积分         0.871 1.78e- 3
## 4 实验组 排便困难 积分         0.915 1.95e- 2
## 5 对照组 排便时间 积分         0.577 3.91e- 8
## 6 实验组 排便时间 积分         0.597 6.64e- 8
## 7 对照组 生活质量 积分         0.974 6.47e- 1
## 8 实验组 生活质量 积分         0.962 3.46e- 1

# 4. 同时方差齐性检验
df_l %>% group_by(变量) %>% levene_test(积分 ~ 组别)
## # A tibble: 4 × 5
##   变量       df1   df2 statistic      p
##   <fct>    <int> <int>     <dbl>  <dbl>
## 1 粪便性状     1    58    0.127  0.723
## 2 排便困难     1    58    0.226  0.636
## 3 排便时间     1    58    0.0746 0.786
## 4 生活质量     1    58    3.96   0.0514

# 5. 同时t检验（Welch校正，默认var.equal = FALSE）
df_l %>% group_by(变量) %>% t_test(积分 ~ 组别)
## # A tibble: 4 × 9
##   变量     .y.   group1 group2    n1    n2 statistic    df     p
## * <fct>    <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>
## 1 粪便性状 积分  对照组 实验组    30    30   -0.356   57.5 0.723
## 2 排便困难 积分  对照组 实验组    30    30    0.0890  57.4 0.929
## 3 排便时间 积分  对照组 实验组    30    30   -0.273   58.0 0.786
## 4 生活质量 积分  对照组 实验组    30    30   -0.101   52.4 0.92
```

### 加入效应量（Cohen's d）

```r
df_l %>% group_by(变量) %>% cohens_d(积分 ~ 组别)
```

### 批量卡方检验示例

```r
# 假设有多组分类变量需要比较
df_l %>%
  group_by(variable_name) %>%
  chisq_test(outcome ~ group)
```

### 结果筛选与导出

```r
# 筛选有统计学意义的变量
df_l %>%
  group_by(变量) %>%
  t_test(积分 ~ 组别) %>%
  filter(p < 0.05)

# 导出为CSV
df_l %>%
  group_by(变量) %>%
  t_test(积分 ~ 组别) %>%
  write.csv("t_test_results.csv")
```

## 结果解读指南

**rstatix 输出对照表（以 t_test 为例）：**

| 列名 | 含义 | 解读要点 |
|------|------|----------|
| `.y.` | 因变量名 | 即为 `values_to` 指定的列名 |
| `group1`, `group2` | 比较的两组 | 对照/实验或任意两组 |
| `n1`, `n2` | 各组样本量 | 确认分组均衡性 |
| `statistic` | 检验统计量 | t值、F值、卡方值等 |
| `df` | 自由度 | t检验用 Welch-Satterthwaite 校正自由度 |
| `p` | p 值 | p < 0.05 表示差异有统计学意义 |

**所有 rstatix 函数输出均为 tibble**，可直接用 `filter()`、`select()`、`mutate()` 等 dplyr 函数进一步处理。

**SPSS vs rstatix：**
- SPSS 中选中多个变量一键输出多个表格
- rstatix 中 `group_by(变量) %>% t_test(...)` 一次性输出所有变量的检验结果于同一个整洁表格
- rstatix 输出可直接用于 ggplot2 可视化，SPSS 需要手动整理

## 常见问题与注意事项

**Q1：宽数据 vs 长数据，哪种格式适合 rstatix？**

rstatix 要求长数据。宽数据中每个变量占一列，需要 `pivot_longer()` 转换后使用。长数据的核心是"变量名列 + 数值列 + 分组列"的三列结构。

**Q2：t_test() 默认是 Student t 还是 Welch t？**

rstatix 的 `t_test()` 默认 `var.equal = FALSE`，即使用 Welch 校正 t 检验。这与 `t.test()` 的默认行为一致。如需等方差假定，显式设置 `var.equal = TRUE`。

**Q3：如何实现多组组间两两比较（post-hoc）？**

```r
# Dunnett多重比较（与参照组比较）
df_l %>%
  group_by(变量) %>%
  t_test(积分 ~ 组别, ref.group = "对照组")

# Tukey HSD事后检验（ANOVA后）
df_l %>%
  group_by(变量) %>%
  tukey_hsd(积分 ~ 组别)
```

**Q4：一次检验多个变量是否增加I类错误风险？**

是的。同时对多个变量进行相同检验属于多重比较，会增加假阳性概率。解决方案：
- 使用 Bonferroni 校正：`p.adjust(results$p, method = "bonferroni")`
- 使用 FDR 校正（更宽松）：`p.adjust(results$p, method = "BH")`
- 在结果中将原始 p 值和校正后 p 值并排呈现

**Q5：shapiro_test() 返回8行结果，但 p 值小于0.05的不止一个，该怎么办？**

这是正态性检验的典型结果。对于不满足正态性的变量，考虑：
- 换用 `wilcox_test()` 进行非参数检验
- 数据转换（对数转换等）后重新检验
- 注意：Shapiro-Wilk 在小样本时检验功效较低，大样本时又过于敏感

**常见错误提醒：**

1. 宽数据直接传入 rstatix 函数而不先 `pivot_longer()`——rstatix 的 formula 接口只接受"数值 ~ 分组"的单变量公式
2. 分组变量为 character 而非 factor——影响排序和参考组选择
3. 忘记 `group_by(变量)` 导致对全部数据做一次检验——结果只有一行而非每变量一行
4. 对同一变量在不同数据集上分别运行检验后再合并——这正是 rstatix 希望消除的重复代码模式
