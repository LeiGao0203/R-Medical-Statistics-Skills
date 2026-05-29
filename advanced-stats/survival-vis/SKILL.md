---
name: medical-stat-survival-vis
description: "R语言医学统计：生存曲线可视化。使用survminer包绘制Kaplan-Meier生存曲线、累积风险曲线、添加risk table、标注中位生存时间和log-rank检验p值。TRIGGER when user mentions KM曲线图、生存曲线图、ggsurvplot、risk table、ncensor plot、中位生存时间标注、生存曲线加总曲线，or asks about plotting survival curves with survminer. SKIP for 生存分析的统计建模（Cox回归）、cut-point寻找、Fine-Gray竞争风险模型。"
---

# 生存曲线可视化 (Survival Curve Visualization)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用场景：**

- 绘制高质量的 Kaplan-Meier 生存曲线用于论文发表
- 比较两组或多组间生存曲线，同时标注 log-rank 检验 p 值
- 在生存曲线下方添加 risk table（at-risk 人数随时间变化表）
- 标注中位生存时间并绘制参考线
- 按多个分类变量进行分面（facet）展示
- 将多条生存函数（如 OS 和 PFS）合并在一张图上

**不使用本技能的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 需拟合 Cox 回归进行多因素分析 | `medical-stat-survival` 技能 |
| 需寻找连续变量最佳 cut-point | `surv_cutpoint()` 配合生存分析技能 |
| 需竞争风险模型的累积发生率曲线 | Fine-Gray 检验技能 |
| 需 Cox 模型诊断图 | `ggcoxzph()`，见生存分析技能 |

## 前置条件

**R 包安装：**

```r
install.packages(c("survival", "survminer"))
# ggsci 提供 NPG、Lancet、JCO 等期刊配色，可选安装
install.packages("ggsci")
```

**核心数据要求：**

- `time`：生存时间（连续型数值，≥0），单位统一
- `status`：结局事件指示变量，标准编码 1 = 事件，0 = 删失。注意 `survival::lung$status` 原始编码为 1 = 删失、2 = 死亡
- 分组变量：因子型，多分类需正确设置水平
- 分面变量：因子型，每个组合需有足够样本量

**前置计算：** 需先通过 `survfit(Surv(time, status) ~ group, data = data)` 构建生存函数对象。

## 方法选择决策树

```
你的可视化需求 →
├── 单条生存曲线（无分组）
│   └── ggsurvplot(fit, data = df, conf.int = TRUE)
│
├── 分组生存曲线（2 组或多组）
│   ├── 添加 log-rank p 值 → pval = TRUE
│   ├── 添加 risk table → risk.table = TRUE
│   ├── 标注中位生存时间 → surv.median.line = "hv"
│   ├── 添加删失计数图 → ncensor.plot = TRUE
│   └── 添加总体曲线 → add.all = TRUE
│
├── 多个分类变量，需分面展示
│   ├── 按一个变量分面 → ggsurvplot_facet(fit, data, facet.by = "rx")
│   ├── 按多个变量分面 → ggsurvplot_facet(fit, data, facet.by = c("rx", "adhere"))
│   └── 手动精细控制 → ggsurv$plot + facet_grid() / ggsurv$table + facet_grid()
│
├── 累积风险曲线（非生存概率）
│   └── fun = "cumhaz" 或 fun = "event"
│
├── 多个生存函数画在一张图上（如 OS + PFS）
│   └── ggsurvplot_combine(list(PFS = fit_pfs, OS = fit_os), data = df)
│
├── 根据某一变量分组绘制多张独立图
│   └── ggsurvplot_group_by(fit, data, group.by = c("rx", "adhere"))
│
├── 多个不同生存函数各自独立出图
│   └── ggsurvplot_list(list(fit1, fit2), data = df)
│
└── 最高度自定义（分别控制 plot / table / ncensor 样式）
    └── 保存 ggsurvplot 对象 → 修改 $plot / $table / $ncensor.plot
```

## 标准工作流

### 步骤 1：数据准备与生存函数构建

```r
library(survival)
library(survminer)

data("lung")
fit <- survfit(Surv(time, status) ~ sex, data = lung)
```

### 步骤 2：基础生存曲线

```r
# 最简调用
ggsurvplot(fit, data = lung)

# 更改删失标记形状和大小
ggsurvplot(fit, data = lung, censor.shape = "|", censor.size = 4)
```

### 步骤 3：发布级生存曲线（p 值 + risk table + 中位生存时间）

```r
ggsurvplot(
  fit, data = lung,
  size = 1, palette = "lancet",
  conf.int = TRUE, pval = TRUE,
  surv.median.line = "hv",
  risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
  legend.labs = c("Male", "Female"),
  ggtheme = theme_classic2()
)
```

### 步骤 4：累积风险曲线 / 累积事件曲线

```r
ggsurvplot(fit, fun = "cumhaz", conf.int = TRUE, palette = "lancet", ggtheme = theme_bw())
ggsurvplot(fit, fun = "event", conf.int = TRUE, palette = "grey", ggtheme = theme_pubclean())
```

### 步骤 5：多变量分面绘制

```r
data(colon)
fit3 <- survfit(Surv(time, status) ~ sex + rx + adhere, data = colon)

ggsurv <- ggsurvplot(fit3, data = colon,
  fun = "cumhaz", conf.int = TRUE,
  risk.table = TRUE, risk.table.col = "strata", ggtheme = theme_bw())

ggsurv$plot + facet_grid(rx ~ adhere)          # 分面生存曲线
ggsurv$table + facet_grid(rx ~ adhere, scales = "free")  # 分面 risk table
```

### 步骤 6：ggsurvplot_facet() 简化分面

```r
fit <- survfit(Surv(time, status) ~ sex, data = colon)
ggsurvplot_facet(fit, colon, facet.by = "rx", palette = "jco", pval = TRUE)
ggsurvplot_facet(fit, colon, facet.by = c("rx", "adhere"), palette = "jco", pval = TRUE)
```

### 步骤 7：根据变量分组独立出图

```r
fit <- survfit(Surv(time, status) ~ sex, data = colon)
ggsurv.list <- ggsurvplot_group_by(fit, colon, group.by = "rx",
  risk.table = TRUE, pval = TRUE, conf.int = TRUE, palette = "jco")
names(ggsurv.list)
## [1] "rx.Obs::sex"     "rx.Lev::sex"     "rx.Lev+5FU::sex"
```

## 代码示例

### 示例 1：完整发布级生存曲线（含所有组件）

```r
library(survival); library(survminer)
fit <- survfit(Surv(time, status) ~ sex, data = lung)

ggsurvplot(
  fit, data = lung,
  size = 1, palette = "lancet",
  conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  surv.median.line = "hv", conf.int.style = "step",
  risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
  ncensor.plot = TRUE, ncensor.plot.height = 0.25,
  xlab = "Time in days", xlim = c(0, 500), break.time.by = 100,
  legend.labs = c("Male", "Female"), ggtheme = theme_classic2(),
  title = "Survival curves", subtitle = "Based on Kaplan-Meier estimates",
  font.title = c(16, "bold", "darkblue"), font.x = c(14, "bold.italic", "red"),
  risk.table.title = "Note the risk set sizes",
  ncensor.plot.title = "Number of censorings"
)
```

### 示例 2：添加总体生存曲线

```r
fit <- surv_fit(Surv(time, status) ~ sex, data = lung)
ggsurvplot(fit, data = lung, risk.table = TRUE, pval = TRUE,
           surv.median.line = "hv", palette = "jco", add.all = TRUE)
```

### 示例 3：PFS 和 OS 画在同一张图上

```r
set.seed(123)
demo.data <- data.frame(
  os.time = colon$time, os.status = colon$status,
  pfs.time = sample(colon$time), pfs.status = colon$status,
  sex = colon$sex, rx = colon$rx, adhere = colon$adhere
)
pfs <- survfit(Surv(pfs.time, pfs.status) ~ 1, data = demo.data)
os  <- survfit(Surv(os.time, os.status) ~ 1, data = demo.data)
ggsurvplot_combine(list(PFS = pfs, OS = os), demo.data)
```

### 示例 4：同时绘制多个生存函数（独立出图）

```r
f1 <- survfit(Surv(time, status) ~ adhere, data = colon)
f2 <- survfit(Surv(time, status) ~ rx, data = colon)
ggsurvplot_list(list(sex = f1, rx = f2), colon, legend.title = list("sex", "rx"))
```

### 示例 5：精细化自定义三个组件

```r
ggsurv <- ggsurvplot(
  fit, data = lung,
  risk.table = TRUE, pval = TRUE, conf.int = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0, 500), xlab = "Time in days", break.time.by = 100,
  ggtheme = theme_light(), risk.table.y.text.col = TRUE,
  risk.table.height = 0.25, risk.table.y.text = FALSE,
  ncensor.plot = TRUE, ncensor.plot.height = 0.25,
  conf.int.style = "step", surv.median.line = "hv",
  legend.labs = c("Male", "Female")
)

ggsurv$plot <- ggsurv$plot + labs(
  title = "Survival curves", subtitle = "Based on Kaplan-Meier estimates")
ggsurv$table <- ggsurv$table + labs(
  title = "Note the risk set sizes")
ggsurv$ncensor.plot <- ggsurv$ncensor.plot + labs(
  title = "Number of censorings")

# 统一修改所有组件字体
customize_labels <- function (p, font.title = NULL, font.subtitle = NULL,
  font.x = NULL, font.y = NULL, font.xtickslab = NULL) {
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else list.plots <- original.p
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown(size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if(!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if(!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if(!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if(!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if(!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]] else list.plots
}

ggsurv <- customize_labels(ggsurv,
  font.title = c(16, "bold", "darkblue"), font.subtitle = c(15, "bold.italic", "purple"),
  font.x = c(14, "bold.italic", "red"), font.y = c(14, "bold.italic", "darkred"),
  font.xtickslab = c(12, "plain", "darkgreen"))
ggsurv
```

## 结果解读指南

### 生存曲线

- **阶梯状下降**：曲线从 1.0 开始，每次下降对应一个或多个终点事件
- **+ 标记**：删失数据位置，表示该患者在此时间点后失访或尚未到达终点
- **曲线间距**：曲线分开越远，两组生存差异越大
- **曲线交叉**：可能提示 PH 假设不成立

### 中位生存时间

- 生存概率首次降至 50% 时对应的随访时间，通过 `surv.median.line = "hv"` 标注
- 论文中常报告为 "Median OS: XX months for group A vs. XX months for group B"

### Log-rank 检验 p 值

- `pval = TRUE` 在图内显示 p 值，p < 0.05 表示组间生存曲线差异有统计学意义

### Risk Table

- 显示各随访时间点剩余 at-risk 人数，`risk.table = TRUE` 添加
- `risk.table.y.text.col = TRUE` 按分组颜色显示数字；`risk.table.y.text = FALSE` 显示为条形

### Ncensor Plot

- 显示各时间区间的删失事件计数，`ncensor.plot = TRUE` 在 risk table 下方添加

### ggsurvplot 返回对象结构

- `$plot`：生存曲线主图（ggplot2 对象）
- `$table`：risk table（ggplot2 对象）
- `$ncensor.plot`：删失计数图（ggplot2 对象）
- `$data.survplot`：绘图用原始数据

## 常见问题与注意事项

### Q1：ggsurvplot 报 "All aesthetics have length 1" 警告

由 `surv.median.line = "hv"` 触发，无害。可忽略或设置 `surv.median.line = "none"` 避免。

### Q2：Palette 配色有哪些可选方案？

支持 ggsci 期刊配色：`"npg"`、`"lancet"`、`"jco"`、`"nejm"`、`"aaas"`、`"jama"`；RColorBrewer 如 `"Dark2"`、`"Set1"`；也可直接传入颜色向量 `c("#E7B800", "#2E9FDF")`。

### Q3：ggsurvplot() / ggsurvplot_list() / ggsurvplot_combine() / ggsurvplot_group_by() 的区别？

- `ggsurvplot()`：一个 survfit 对象出一张图，分组以颜色区分
- `ggsurvplot_list()`：survfit 列表，每个元素出独立一张图
- `ggsurvplot_combine()`：survfit 列表，所有函数画在同一张图上（如 PFS + OS）
- `ggsurvplot_group_by()`：给定分组变量，按各水平输出多张独立图，效果等同 `ggsurvplot_facet()` 但不分面

### Q4：如何保存高质量图片用于论文？

```r
p <- ggsurvplot(fit, data = lung, ...)
ggsave("km_curve.pdf", plot = print(p), width = 8, height = 6, dpi = 300)
# or
ggsave("km_curve.png", plot = print(p), width = 8, height = 6, dpi = 600)
```

注意 `print(p)` 是必需的，因 `ggsurvplot` 返回的并非普通 ggplot 对象。

### Q5：如何把分面后的生存曲线和 risk table 用 gridExtra 组合？

```r
g2 <- ggplotGrob(curv_facet)
g3 <- ggplotGrob(tbl_facet)
min_ncol <- min(ncol(g2), ncol(g3))
g <- gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size = "last")
g$widths <- grid::unit.pmax(g2$widths, g3$widths)
grid::grid.newpage()
grid::grid.draw(g)
```

### Q6：lung 数据集 status 编码为 1=删失, 2=死亡，survfit 能正确处理吗？

`survival::lung` 被 survival 包内部标记为特殊类型，`Surv(time, status)` 能自动识别。但建议养成确认编码的习惯，使用 `?lung` 查阅帮助文档。

### 关键提醒

- `ggsurvplot()` 返回的并非普通 ggplot 对象，用 `ggsave()` 保存时需 `print(p)` 包裹
- `risk.table = TRUE` 会增加图像高度，根据分组数调整 `risk.table.height`（0.2-0.5）
- 多组（>4 组）时建议精简颜色或使用 `palette = "grey"` 避免视觉混乱
- 分面图中某组合样本量过小（< 10）时 KM 曲线不稳定，不建议展示
- `conf.int = TRUE` 绘制 95% 可信区间，样本量小时区间可能较宽，属正常现象
