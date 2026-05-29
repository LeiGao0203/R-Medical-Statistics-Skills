---
name: medical-stat-subgroup-onecode
description: "R语言医学统计：一行代码实现亚组分析。使用jstable包一次性完成多个亚组的效应估计和森林图，极大简化工作流。支持Cox回归、logistic回归、svyglm和svycoxph。TRIGGER when user mentions 一键亚组分析、批量亚组、一行代码森林图、jstable、高效亚组分析、TableSubgroupMultiCox，or asks about streamlining subgroup analysis. SKIP for 手动亚组分析（purrr/dplyr手写循环）、Meta分析森林图、多因素回归森林图比较。"
---

# 亚组分析1行代码实现 (Subgroup Analysis in One Line)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用一行代码亚组分析的典型场景：**

- 需要在多个分层变量（如性别、年龄组、肿瘤分期）下同时评估处理效应（HR 或 OR），快速产出表格
- 需要森林图展示各亚组的效应大小及交互作用 P 值，但不想手写 `purrr` 循环拼接
- 处理因素为二分类（如治疗 vs 对照），结果变量为生存数据或二分类结局
- 使用了复杂调查设计（`svyglm` / `svycoxph`），jstable 同样支持

**不使用的场景：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 需要完全控制每个亚组的模型细节（协变量调整、分层变量不同） | 手动 `purrr` 循环 + `coxph` / `glm`（见第44章） |
| 需要比较亚组分析和多因素回归的森林图差异 | 见第46章《亚组分析和多因素回归的森林图比较》 |
| 属于 Meta 分析的亚组森林图 | 使用 `meta` 包或 `metafor` 包 |
| 需要自定义森林图的每一个元素（字体、颜色、间距） | jstable 输出可直接传入 `forestploter` 进一步美化 |

## 前置条件

**R 包安装：**

```r
install.packages(c("jstable", "survival", "forestploter", "tidyverse"))
# 如需 GitHub 最新版：
# remotes::install_github("jinseob2kim/jstable")
```

**核心数据要求：**

- 处理变量：二分类（0/1 或 factor）
- 亚组变量：分类变量（factor 型），每个变量会被自动识别所有水平
- 生存数据：需 `time`（生存时间，数值型）和 `status`（事件状态，0/1，1 表示发生事件）
- 分类变量必须为 **factor 类型**，否则 jstable 不会自动进行哑变量编码
- 缺失值需提前处理（`na.omit()` 或插补），jstable 自身不处理缺失

**jstable 支持的模型类型：**

| 函数 | 模型 | 结局类型 |
|------|------|----------|
| `TableSubgroupMultiCox()` | Cox 比例风险模型 | 生存数据（time + status） |
| `TableSubgroupMultiGLM()` | 广义线性模型 | 二分类/连续型 |
| 对 `svyglm` 和 `svycoxph` 对象 | 复杂调查设计 | 同上 |

## 方法选择决策树

```
你的需求 →
├── 模型类型
│   ├── 生存结局（time + status） → TableSubgroupMultiCox()
│   ├── 二分类或连续型结局 → TableSubgroupMultiGLM()
│   └── 复杂调查设计数据 → 先构建 svyglm/svycoxph 对象，再传入对应函数
│
├── 森林图需求
│   ├── 快速出图，用 jstable 输出直接画 → 用 forestploter::forest()，数据稍作清洗即可
│   ├── 需要 HR(95%CI) 格式的列 → 在 plot_df 中手动 paste() 拼接
│   │    如：plot_df$HR_CI <- paste0(sprintf("%.2f", plot_df[["Point Estimate"]]),
│   │        " (", sprintf("%.2f", plot_df$Lower), "-", sprintf("%.2f", plot_df$Upper), ")")
│   └── 需要高度定制的森林图 → 将 jstable 输出作为数据源，用 ggplot2 手动绘制
│
├── 亚组变量选择
│   ├── 所有二分类/多分类变量都参与亚组分析 → var_subgroups = names(df)[c(…)]
│   ├── 仅展示关键临床变量 → 手选 3–5 个最有临床意义的变量
│   └── 变量太多导致图太拥挤 → 优先选取交互作用 p < 0.1 的变量展示
│
└── P for interaction 解读
    ├── p < 0.05 → 该亚组变量与处理因素存在显著的交互作用，效应在亚组间有差异
    └── p ≥ 0.05 → 无显著交互作用，各亚组效应方向一致
```

## 标准工作流

### 步骤 1：数据准备

分类变量全部转为 factor，处理变量为二分类（factor 或 0/1）。确保无缺失值。以 `survival::colon` 数据集为例，筛选 etype == 1（复发事件），仅保留 Obs 组和 Lev+5FU 组。

### 步骤 2：一行代码执行亚组分析

调用 `TableSubgroupMultiCox(formula, var_subgroups, data)`，指定生存公式、亚组变量字符向量和数据框。结果直接为 data.frame，包含各亚组的 HR、95%CI、两组事件率、各亚组内的 P 值，以及最右侧的 P for interaction。

### 步骤 3：清洗结果用于森林图

将结果中的 NA 替换为空格，数值列转为 numeric 类型，添加空白列用于 `forestploter` 的置信区间显示区域。

### 步骤 4：绘制森林图

使用 `forestploter::forest()` 指定展示列、`est`、`lower`、`upper`、`ci_column`、`ref_line` 和 `xlim`。可以将输出保存为 pdf/png。

## 代码示例

```r
library(jstable)
library(survival)
library(tidyverse)
library(forestploter)

rm(list = ls())

df <- colon %>%
  mutate(rx = as.numeric(rx)) %>%
  filter(etype == 1, !rx == 2) %>%
  select(time, status, rx, sex, age, obstruct, perfor, adhere,
         differ, extent, surg, node4) %>%
  mutate(
    sex     = factor(sex, levels = c(0, 1), labels = c("female", "male")),
    age     = ifelse(age > 65, ">65", "<=65"),
    age     = factor(age, levels = c(">65", "<=65")),
    obstruct = factor(obstruct, levels = c(0, 1), labels = c("No", "Yes")),
    perfor  = factor(perfor, levels = c(0, 1), labels = c("No", "Yes")),
    adhere  = factor(adhere, levels = c(0, 1), labels = c("No", "Yes")),
    differ  = factor(differ, levels = c(1, 2, 3),
                     labels = c("well", "moderate", "poor")),
    extent  = factor(extent, levels = c(1, 2, 3, 4),
                     labels = c("submucosa", "muscle", "serosa", "contiguous")),
    surg    = factor(surg, levels = c(0, 1), labels = c("short", "long")),
    node4   = factor(node4, levels = c(0, 1), labels = c("No", "Yes")),
    rx      = ifelse(rx == 3, 0, 1),
    rx      = factor(rx, levels = c(0, 1))
  )

str(df)
## 'data.frame':  619 obs. of  12 variables

res <- TableSubgroupMultiCox(
  formula = Surv(time, status) ~ rx,
  var_subgroups = c("sex", "age", "obstruct", "perfor", "adhere",
                    "differ", "extent", "surg", "node4"),
  data = df
)

res
##         Variable Count Percent Point Estimate Lower Upper rx=0 rx=1 P value
## rx       Overall   619     100           1.67  1.32  2.11 34.4 48.9  <0.001
## 1            sex  <NA>    <NA>           <NA>  <NA>  <NA> <NA> <NA>    <NA>
## 2         female   312    50.4           1.32  0.96   1.8 41.1 47.8   0.088
## 3           male   307    49.6           2.29   1.6  3.26 26.6 50.1  <0.001
## ... (共9个亚组变量，含P for interaction列)

plot_df <- res
plot_df[, c(2, 3, 9)][is.na(plot_df[, c(2, 3, 9)])] <- " "
plot_df$` ` <- paste(rep(" ", nrow(plot_df)), collapse = " ")
plot_df[, 4:6] <- apply(plot_df[, 4:6], 2, as.numeric)

# 如需 HR(95%CI) 格式
plot_df$HR_CI <- paste0(
  sprintf("%.2f", plot_df[["Point Estimate"]]),
  " (", sprintf("%.2f", plot_df$Lower),
  "-", sprintf("%.2f", plot_df$Upper), ")"
)

p <- forest(
  data     = plot_df[, c(1, 2, 3, 11, 9)],
  lower    = plot_df$Lower,
  upper    = plot_df$Upper,
  est      = plot_df[["Point Estimate"]],
  ci_column = 4,
  ref_line = 1,
  xlim     = c(0.1, 4)
)
print(p)
## 输出为 PDF: ggsave("forest.pdf", p, width = 10, height = 8)
```

## 结果解读指南

**输出表格各列含义：**

| 列名 | 含义 |
|------|------|
| `Variable` | 亚组变量名及其水平；`Overall` 行显示总人群效应 |
| `Count` | 该亚组内的样本量 |
| `Percent` | 该亚组占总样本百分比 |
| `Point Estimate` | 该亚组内处理因素的 HR（Cox）或 OR（GLM） |
| `Lower` / `Upper` | 效应估计的 95% 置信区间下限/上限 |
| `rx=0` / `rx=1` | 该亚组内对照组和处理组的事件率（%） |
| `P value` | 该亚组内处理效应的 P 值 |
| `P for interaction` | 亚组变量与处理因素交互作用的 P 值 |

**关键解读规则：**

- **P for interaction < 0.05**：效应在亚组间存在显著差异，处理效果在某一亚组中更强或更弱
- **P for interaction ≥ 0.05**：无显著交互作用，各亚组的效应方向一致，差异可能由随机误差导致
- 交叉线（ref_line = 1）穿过方块表示该亚组内处理效应不显著

**论文中如何报告：**

"采用 jstable 包的 TableSubgroupMultiCox 函数进行亚组分析，以性别、年龄、肠梗阻、肠穿孔、粘连、肿瘤分化程度、局部扩散、手术时长、淋巴结状态为分层变量。森林图展示了各亚组的 HR（95% CI）及交互作用 P 值。结果显示，Lev+5FU 治疗与观察组相比，总体 HR 为 1.67（95%CI: 1.32–2.11），性别亚组间的交互作用 P = 0.029，提示治疗效应在男性和女性中存在差异。"

## 常见问题与注意事项

**Q1: jstable 和手动 purrr 循环做亚组分析，有什么区别？**

| 特性 | jstable | purrr 手动循环 |
|------|---------|---------------|
| 代码量 | 1 行 | 10–20 行 |
| 输出格式 | 直接为可画图的 data.frame | 需自行拼接、规整 |
| 灵活性 | 较低，模型固定 | 较高，可自定义每个模型 |
| 交互作用 P | 自动计算 | 需手动添加交互项 |
| 适用场景 | 标准 Cox/GLM 亚组分析 | 复杂模型、自定义调整 |

**Q2: 分类变量必须是 factor 吗？**

是的。jstable 依赖 factor 的水平识别亚组内的各个类别。如果把分类变量保留为 numeric（0/1/2），各个数值会被当作连续变量处理，不会自动分组。

**Q3: 森林图中 NA 行（亚组标题行）怎么处理？**

jstable 的输出在每个亚组变量前插入 NA 行作为标题分隔。画图前把 `Count`、`Percent`、`P value` 列的 NA 替换为空格 `" "`，`Point Estimate`、`Lower`、`Upper` 保持 NA 即可（`forestploter` 会在该行留空）。

**Q4: 警告 "Loglik converged before variable" 是什么意思？**

在小样本亚组中，Cox 模型的对数似然在参数估计完成前就已收敛。常见原因是该亚组内某一组的全部对象均达标或均不达标。此时 HR 可能为 0 或 Inf，应谨慎解读。

**Q5: 可以自定义森林图中展示哪些列吗？**

`forestploter::forest()` 的 `data` 参数接受任意列组合。选择 `plot_df[, c(1, 2, 3, 11, 9)]` 即 Variable、Count、Percent、空白置信区间列、P value。调整列索引即可自定义。

**Q6: jstable 支持 logistic 回归吗？**

使用 `TableSubgroupMultiGLM()` 并在 formula 中指定二分类结局（`outcome ~ rx`）、`family = "binomial"`。输出为 OR 而非 HR。

**Q7: 森林图的 x 轴范围怎么调？**

使用 `forest()` 的 `xlim` 参数。如 `xlim = c(0.1, 6)`。若某个亚组 CI 超出 xlim，`forestploter` 会用箭头标示。

**Q8: 如何导出森林图为发表级质量？**

`forestploter` 输出的 `p` 对象可用 `ggsave("forest.pdf", p, width = 10, height = 8)` 导出为 PDF/SVG。

---

**参考资料：**

1. jstable 包文档：https://jinseob2kim.github.io/jstable/
2. forestploter 包文档：https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html
3. 手动亚组分析（第44章）：`1041-subgroupanalysis.html`
4. 亚组分析与多因素回归森林图比较（第46章）：`亚组分析和多因素回归的森林图.html`
