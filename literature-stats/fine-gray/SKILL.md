---
name: medical-stat-fine-gray
description: "R语言医学统计：Fine-Gray竞争风险模型。处理存在竞争事件的生存数据，估计累积发生率函数（CIF），使用cmprsk包进行Fine-Gray检验和竞争风险回归。TRIGGER when user mentions 竞争风险、Fine-Gray、CIF、累积发生率、cmprsk、competing risks，or asks about survival analysis with competing events. SKIP for 标准Cox回归（无竞争风险）、Kaplan-Meier、生存曲线可视化。"
---

# Fine-Gray检验和竞争风险模型 (Fine-Gray Test & Competing Risks Model)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用竞争风险模型的典型场景：**

- 生存分析中存在多个可能的终点事件，其中某些事件会"阻碍"或"影响"感兴趣事件的发生
- 例如：研究肠癌复发（感兴趣事件），但患者可能因心梗死亡（竞争风险事件）而无法观察到复发
- 需估计某事件在竞争风险存在下的累积发生率（CIF, Cumulative Incidence Function）
- 需校正协变量后评估某因素对感兴趣事件风险的影响（Fine-Gray 多因素回归）

**不使用竞争风险模型的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 仅一个终点事件，无竞争风险 | 标准 Cox 回归 / log-rank 检验 |
| 仅需绘制生存曲线（K-M 曲线） | Kaplan-Meier + `ggsurvplot()` |
| 需要同时建模多个事件及其转换（如健康→复发→死亡） | 多状态模型（multi-state model） |
| 研究终点本身没有竞争关系 | Cox 回归即可，无需 Fine-Gray |

**医学研究常见应用：**

- 肿瘤研究中，感兴趣事件为肿瘤复发/进展，竞争事件为非肿瘤死亡
- 造血干细胞移植研究中，感兴趣事件为复发，竞争事件为移植相关死亡
- 心血管研究中，感兴趣事件为心肌梗死，竞争事件为卒中死亡

## 前置条件

**R 包安装：**

```r
install.packages(c("cmprsk", "casebase", "ggplot2"))
```

**核心数据要求：**

- `ftime`：生存时间/随访时间（连续型数值，≥0）
- `Status`（fstatus）：结局指示变量，三分类 —— 0 = 删失，1 = 感兴趣事件，2 = 竞争风险事件
  - 也可使用其他编码，但需在 `cuminc()` 中用 group 参数或 `crr()` 中用 failcode/cencode 显式指定
- 分组变量（单因素比较）：因子型（factor）
- 协变量（多因素回归）：需转为数值型（整数），分类变量需手工做哑变量编码

**关键注意：**
- `cmprsk::crr()` **不会自动对分类变量做哑变量编码**，必须手动 `model.matrix()` 或自行构建哑变量
- 与标准 Cox 回归不同，Fine-Gray 模型直接对 CIF 建模（非风险函数），回归结果为 subdistribution hazard ratio（SHR）

## 方法选择决策树

```
你的研究目标 →
├── 仅比较组间感兴趣的累积发生率（CIF）
│   └── 单因素 Fine-Gray 检验 → `cuminc(ftime, fstatus, group)`
│
├── 需要校正多个协变量，评估某因素的独立效应
│   ├── 协变量为连续/已编码数值型 → `crr(ftime, fstatus, cov, failcode=1, cencode=0)`
│   └── 协变量含多分类因子型 → 先用 `model.matrix()` 哑变量化，再传入 `crr()`
│
├── 需要可视化 CIF 曲线
│   ├── 快速查看 → `plot(cuminc_result)` 基础图形
│   ├── 分组较少（≤4组）→ 使用 ggplot2 长数据格式化绘制
│   └── 分组较多（≥5组）→ 分面（facet_wrap / facet_grid）
│
├── 同时关心多个结局事件的回归系数
│   └── 分别对每种事件（failcode=1, failcode=2）各拟合一个 `crr()` 模型
│
└── 需要提取模型结果用于制表
    └── 手动从 `summary(crr_result)` 提取 coef 表，或自行构建 data.frame
```

## 标准工作流

### 步骤 1：数据准备与探索

```r
library(cmprsk)

data("bmtcrr", package = "casebase")
str(bmtcrr)
## 'data.frame':    177 obs. of  7 variables:
##  $ Sex   : Factor w/ 2 levels "F","M"
##  $ D     : Factor w/ 2 levels "ALL","AML"
##  $ Phase : Factor w/ 4 levels "CR1","CR2","CR3","Relapse"
##  $ Age   : int  48 23 7 26 36 ...
##  $ Status: int  2 1 0 2 2 ...
##  $ Source: Factor w/ 2 levels "BM+PB","PB"
##  $ ftime : num  0.67 9.5 131.77 24.03 1.47 ...
```

变量说明：

| 变量 | 含义 | 类型 |
|------|------|------|
| `Sex` | 性别（F/M） | 因子 |
| `D` | 疾病类型（ALL / AML） | 因子 |
| `Phase` | 疾病阶段（CR1 / CR2 / CR3 / Relapse） | 因子 |
| `Age` | 年龄 | 整数 |
| `Status` | 结局：0 = 删失，1 = 复发，2 = 竞争风险 | 整数 |
| `Source` | 移植方式（BM+PB / PB） | 因子 |
| `ftime` | 生存时间（月） | 数值 |

研究背景：探讨骨髓移植+血液移植（BM+PB）与单纯血液移植（PB）治疗白血病的疗效。感兴趣事件为复发（Status=1），竞争风险事件为移植不良反应死亡（Status=2）。

### 步骤 2：单因素 Fine-Gray 检验（组间 CIF 比较）

```r
bmtcrr$Status <- factor(bmtcrr$Status)
f <- cuminc(bmtcrr$ftime, bmtcrr$Status, bmtcrr$D)
f
## Tests:
##        stat         pv df
## 1 2.8623325 0.09067592  1
## 2 0.4481279 0.50322531  1
```

### 步骤 3：绘制 CIF 曲线

```r
plot(f, xlab = 'Month', ylab = 'CIF', lwd = 2, lty = 1,
     col = c('red', 'blue', 'black', 'forestgreen'))
```

或用 ggplot2 绘制：

```r
library(ggplot2)

# 手工提取 cuminc 对象数据
ALL1 <- data.frame(ALL1_t = f[[1]][[1]], ALL1_C = f[[1]][[2]])
AML1 <- data.frame(AML1_t = f[[2]][[1]], AML1_C = f[[2]][[2]])
ALL2 <- data.frame(ALL2_t = f[[3]][[1]], ALL2_C = f[[3]][[2]])
AML2 <- data.frame(AML2_t = f[[4]][[1]], AML2_C = f[[4]][[2]])

tmp <- data.frame(
  month = c(ALL1$ALL1_t, AML1$AML1_t, ALL2$ALL2_t, AML2$AML2_t),
  cif   = c(ALL1$ALL1_C, AML1$AML1_C, ALL2$ALL2_C, AML2$AML2_C),
  type  = rep(c("ALL1", "AML1", "ALL2", "AML2"),
              c(nrow(ALL1), nrow(AML1), nrow(ALL2), nrow(AML2)))
)

ggplot(tmp, aes(month, cif)) +
  geom_line(aes(color = type, group = type), linewidth = 1.2) +
  theme_bw() +
  theme(legend.position = "top")
```

### 步骤 4：多因素 Fine-Gray 回归

```r
# 剔除结局变量，自变量转为整数
covs <- subset(bmtcrr, select = -c(ftime, Status))
covs[, c(1:3, 5)] <- lapply(covs[, c(1:3, 5)], as.integer)

# 拟合模型：failcode=1 指定复发为感兴趣事件，cencode=0 指定删失编码
f2 <- crr(bmtcrr$ftime, bmtcrr$Status, covs, failcode = 1, cencode = 0)
summary(f2)
## Competing Risks Regression
##
##           coef exp(coef) se(coef)      z p-value
## Sex     0.0494     1.051   0.2867  0.172 0.86000
## D      -0.4860     0.615   0.3040 -1.599 0.11000
## Phase   0.4144     1.514   0.1194  3.470 0.00052
## Age    -0.0174     0.983   0.0118 -1.465 0.14000
## Source  0.9526     2.592   0.5469  1.742 0.08200
##
##        exp(coef) exp(-coef)  2.5% 97.5%
## Sex        1.051      0.952 0.599  1.84
## D          0.615      1.626 0.339  1.12
## Phase      1.514      0.661 1.198  1.91
## Age        0.983      1.018 0.960  1.01
## Source     2.592      0.386 0.888  7.57
##
## Num. cases = 177
## Pseudo Log-likelihood = -267
## Pseudo likelihood ratio test = 23.6  on 5 df
```

### 步骤 5：结果报告

论文中可报告为：在控制了竞争风险事件（移植相关死亡）后，多因素 Fine-Gray 回归显示，疾病所处阶段（Phase）是患者复发的独立影响因素（SHR = 1.51, 95% CI: 1.20–1.91, p = 0.00052）。

## 代码示例

### 示例 1：完整单因素分析流程

```r
rm(list = ls())

library(cmprsk)

# 数据加载
data("bmtcrr", package = "casebase")

# 确保 Status 为因子型
bmtcrr$Status <- factor(bmtcrr$Status)

# Fine-Gray 检验：比较 ALL vs AML 的累计复发率
f <- cuminc(bmtcrr$ftime, bmtcrr$Status, bmtcrr$D)
f
## Tests:
##        stat         pv df
## 1 2.8623325 0.09067592  1
## 2 0.4481279 0.50322531  1
##
## Estimates and Variances:
## $est
##              20        40        60        80       100       120
## ALL 1 0.3713851 0.3875571 0.3875571 0.3875571 0.3875571 0.3875571
## AML 1 0.2414530 0.2663827 0.2810390 0.2810390 0.2810390        NA
## ALL 2 0.3698630 0.3860350 0.3860350 0.3860350 0.3860350 0.3860350
## AML 2 0.4439103 0.4551473 0.4551473 0.4551473 0.4551473        NA

# CIF 曲线绘制
plot(f, xlab = 'Month', ylab = 'CIF', lwd = 2, lty = 1,
     col = c('red', 'blue', 'black', 'forestgreen'))
```

### 示例 2：ggplot2 绘制 CIF 曲线（美化版）

```r
rm(list = ls())

library(cmprsk)
library(ggplot2)

data("bmtcrr", package = "casebase")
bmtcrr$Status <- factor(bmtcrr$Status)

f <- cuminc(bmtcrr$ftime, bmtcrr$Status, bmtcrr$D)

# 提取各组数据
extract_cif <- function(obj, index) {
  data.frame(time = obj[[index]][[1]], cdf = obj[[index]][[2]])
}

ALL1 <- extract_cif(f, 1)
AML1 <- extract_cif(f, 2)
ALL2 <- extract_cif(f, 3)
AML2 <- extract_cif(f, 4)

tmp <- data.frame(
  month = c(ALL1$time, AML1$time, ALL2$time, AML2$time),
  cif   = c(ALL1$cdf, AML1$cdf, ALL2$cdf, AML2$cdf),
  type  = rep(c("ALL复发", "AML复发", "ALL竞争事件", "AML竞争事件"),
              c(nrow(ALL1), nrow(AML1), nrow(ALL2), nrow(AML2)))
)

ggplot(tmp, aes(month, cif, color = type, linetype = type)) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = c("red", "blue", "black", "forestgreen")) +
  labs(x = "随访时间（月）", y = "累积发生率（CIF）") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank())
```

### 示例 3：多因素 Fine-Gray 回归（含哑变量编码）

```r
rm(list = ls())

library(cmprsk)

data("bmtcrr", package = "casebase")

# 分类变量转哑变量
X <- model.matrix(~ Sex + D + Phase + Age + Source - 1, data = bmtcrr)
X <- X[, -1, drop = FALSE]  # 去掉截距列的重复（crr 自带截距）

f_adj <- crr(bmtcrr$ftime, bmtcrr$Status, X, failcode = 1, cencode = 0)
summary(f_adj)
```

## 结果解读指南

### `cuminc()` 输出解读

- **Tests 表**：每一行对应一个事件类型（按编码升序排列，编码 1 = 感兴趣事件，编码 2 = 竞争风险事件）
  - `stat`：Fine-Gray 检验统计量
  - `pv`：p 值，p < 0.05 表示该事件在不同分组间的累积发生率差异有统计学意义
  - `df`：自由度（组数 - 1）
- **$est**：各时间点各组的 CIF 估计值（如 ALL 1 = 0.371 表示 20 个月时 ALL 组累计复发率为 37.1%）
- **$var**：对应的方差估计值

**解读范例**：Fine-Gray 检验显示，ALL 组与 AML 组的累计复发率无统计学差异（p = 0.091）；两组的累计竞争风险事件发生率也无统计学差异（p = 0.503）。

### `crr()` 输出解读

- **coef**：回归系数（log SHR），正值表示 subdistribution hazard 增加（危险因素），负值表示降低（保护因素）
- **exp(coef)**：SHR（subdistribution hazard ratio）。SHR > 1 危险因素，SHR < 1 保护因素
- **se(coef)**：回归系数的标准误
- **z**：Wald z 统计量（coef / se(coef)）
- **p-value**：Wald 检验 p 值，p < 0.05 表示在控制竞争风险后该变量对感兴趣事件有独立影响
- **2.5% / 97.5%**：SHR 的 95% 置信区间。区间不跨 1 → 有统计学意义
- **Pseudo likelihood ratio test**：模型整体似然比检验

**论文报告示例**：在控制了竞争风险事件（移植相关死亡）后，多因素 Fine-Gray 回归分析显示，Phase（疾病阶段）是复发的独立危险因素（SHR = 1.51, 95% CI: 1.20–1.91, p < 0.001）。

### CIF 曲线解读

- 纵轴：累积发生率（Cumulative Incidence Function, CIF）
- 横轴：时间（在对应单位的节点处 CIF 发生阶梯式变化）
- 实线/不同颜色区分不同组别和事件类型
- 曲线越高 → 该组该事件的累计发生率越高

**图例命名**：`ALL1` = ALL 组复发，`AML1` = AML 组复发，`ALL2` = ALL 组竞争事件，`AML2` = AML 组竞争事件。

## 常见问题与注意事项

### Q1：Fine-Gray 回归与标准 Cox 回归的根本区别是什么？

Cox 回归将竞争风险事件视为删失，假设删失与感兴趣事件独立。但当竞争风险事件的发生会"阻碍"感兴趣事件的观察时，这一假设不成立。Fine-Gray 模型直接对 CIF 建模，不对竞争风险事件做删失处理，而是将它们保留在 risk set 中，结果为 subdistribution hazard ratio（SHR），不可直接等同于标准 Cox 回归的 HR。

### Q2：`crr()` 中分类变量为什么不自动哑变量化？

`crr()` 的 `cov1` 参数要求传入**数值矩阵或数据框**，内部不做因子转换。若传入因子型数据框，函数可能报错或将其强制转为错误的连续数值。因此必须手动哑变量编码：单分类变量可用 `as.integer()`（二分类），多分类用 `model.matrix()`。

### Q3：什么时候应该用 Fine-Gray 而不是标准 Cox？

当研究中存在无论如何都"阻止"感兴趣事件发生的竞争事件时，应优先使用 Fine-Gray。例如肿瘤复发研究中，非肿瘤死亡使患者无法再观察到复发。相反，如果竞争事件不太可能影响感兴趣事件的发生（如失访、非相关性住院），标准 Cox 回归通常足够。

### Q4：数据中的 Status 编码不是 0/1/2 怎么办？

`cuminc()` 中只要分组变量的不同水平代表不同事件即可。`crr()` 中通过 `failcode` 指定感兴趣事件的编码，`cencode` 指定删失的编码，其他编码均被视为竞争风险事件。建议统一使用 0 = 删失，1 = 感兴趣事件，2 = 竞争风险事件的编码习惯。

### Q5：如何解释 SHR 与 HR 的不同？

SHR（subdistribution hazard ratio）的解读与 HR 类似（> 1 为危险，< 1 为保护），但其分母是那些"尚未发生感兴趣事件"的被研究者，其中仍包含"已经发生了竞争风险事件"的个体（在传统生存分析中这些个体已被剔除）。因此 SHR 反映了在竞争风险存在时，该因素对感兴趣事件累积发生率的直接效应。

### Q6：broom 包能提取 Fine-Gray 模型结果吗？

`broom` 暂不支持 `cmprsk::crr()` 对象。需手动从 `summary(f2)` 输出中提取系数表，或用 `f2$coef` 直接获取系数向量，自行构建 data.frame 用于制表。

### Q7：SPSS 与 R 的 Fine-Gray 分析有何差异？

SPSS 中 Fine-Gray 模型通常需额外安装扩展模块或在 STATA 中实现。R 的 `cmprsk` 包是 Fine-Gray 方法的原始作者之一发布的参考实现，结果更可靠且可编程化。SPSS 的默认输出中 SHR 较为直观，但手动编码分类变量的要求与 R 一致。

### 关键提醒

- `crr()` 不做哑变量编码——这是最常见的错误来源，务必先手工 encoding
- `cuminc()` 中的 `group` 变量只有两个水平时只输出两行 Tests 结果（一行对应一个事件）；超过两个水平时行数 = 事件类型数 × 水平数
- Pseudo likelihood ratio test 的 p 值报告的是模型整体显著性，不可替代单个变量的 Wald p 值
- 原始 `plot.cuminc()` 使用基础图形，y 轴标签 "Probability" 实际为 CIF 而非传统概率
