---
name: medical-stat-spline-reg
description: "R语言医学统计：限制性立方样条（RCS）回归。使用rms包在回归模型中拟合非线性关系，自动选择节点（knot），检验非线性p值，常用于Cox回归和Logistic回归的非线性建模。TRIGGER when user mentions 样条回归、RCS、限制性立方样条、spline、非线性建模、rcs()、节点选择，or asks about modeling non-linear dose-response relationships. SKIP for 多项式拟合、广义加性模型(GAM)。"
---

# 样条回归 (Spline Regression / Restricted Cubic Splines)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用限制性立方样条的典型场景：**

- 变量间存在明确的非线性关系，直线回归拟合效果差，需要捕捉复杂的非线性趋势
- 需要同时报告非线性效应的显著性检验（Nonlinear p 值），便于判断是否偏离线性假设
- 在 Cox 回归或 Logistic 回归中，希望可视化某连续变量与 HR/OR 之间的非线性剂量-反应关系
- 多项式回归拟合线在数据两端出现不自然的翘起（Runge 现象），需要更稳定的边缘行为
- 文献中常见的"RCS 平滑曲线 + 阴影置信区间 + HR=1 参考线"组合图

**不使用样条回归的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 非线性模式很简单（U型、倒U型），用二次项就能描述 | 多项式回归（`y ~ x + I(x^2)`） |
| 需要完全非参数、数据驱动的平滑 | 广义加性模型（GAM, `mgcv::gam()`） |
| 数据存在已知的生物断点/阈值 | 分段回归 |
| 变量间关系接近直线 | 简单线性回归 / Cox 回归 / Logistic 回归的线性项 |
| 自变量是分类变量 | 无需样条，直接用分类变量编码 |

**医学研究常见应用：**

- Cox 回归中年龄（连续变量）与死亡风险的 HR 曲线，确定保护因素和危险因素的年龄分界点
- Logistic 回归中 BMI 与患病 OR 之间的 J 型 / U 型关系
- 红细胞分布宽度（RDW）、中性粒细胞/淋巴细胞比值（NLR）等实验室指标与预后的非线性量效关系
- 生存分析中连续型生物标志物与结局事件的非线性风险曲线

## 前置条件

**R 包安装：**

```r
install.packages(c("rms", "ggplot2"))
```

**数据格式要求：**

- 自变量 `x`：连续型数值变量（如年龄、BMI、生化指标、随访时间）
- 因变量：
  - 线性回归：连续型数值变量 `y`
  - Logistic 回归：二分类变量（0/1 或 factor）
  - Cox 回归：`Surv(time, status)` 生存对象
- 数据框格式，每行一个观测

**统计假设：**

1. **线性回归中的样条**：残差仍需满足正态性、独立性和等方差性
2. **Logistic / Cox 回归中的样条**：满足对应回归模型的标准假定（对数线性假设、比例风险假定等）
3. **节点数选择**：通常 3–5 个节点即可，节点太少拟合不够灵活，节点太多容易过度拟合。Harrell 推荐默认 4 个节点，最多不超过 5 个
4. **样本量**：样条回归的灵活性随节点增加而提高，过小的样本量（< 100）不宜使用过多节点

## 方法选择决策树

```
你的数据情况 →
├── 因变量为连续型数值变量 → lm(y ~ rcs(x, 节数))
├── 因变量为二分类变量（0/1） →
│   ├── 使用 rms::lrm(binary ~ rcs(x, 节数) + covariates)
│   ├── 配合 anova() 检验 Nonlinear p 值
│   └── 配合 ggplot(Predict()) 绘制 OR 曲线
├── 因变量为生存数据（Surv 对象） →
│   ├── 使用 rms::cph(Surv(time, status) ~ rcs(x, 节数) + covariates)
│   ├── 配合 anova() 检验 Nonlinear p 值
│   └── 配合 ggplot(Predict(fit, age, fun = exp)) 绘制 HR 曲线
├── 不知道选几个节点 →
│   ├── 分别尝试 3、4、5 个节点
│   ├── 比较 AIC 值，选择 AIC 最小的模型
│   └── 实战中 4 个节点最常用
├── 需要在 HR/OR=1 处画竖线分割 →
│   ├── 使用 Predict() 查找 HR/OR 接近 1 的 x 值
│   └── 修改 dd$limits$age[2] 为该值后 update() 模型，重新绘图
└── 需要使用 splines 包的 ns() 而非 rms::rcs() → ns() 为自然三次样条，概念类似，但 rms 包的 Predict/ggplot 绘制最为便捷
```

## 标准工作流

### 步骤1：数据准备与包加载

加载 `rms` 和 `ggplot2`，如果是逻辑回归或生存分析，准备对应的数据和二分类结局/生存对象。

### 步骤2：打包数据（rms 特有步骤）

使用 `rms` 系列函数前，必须用 `datadist()` 打包数据并设置全局选项：

```r
dd <- datadist(mydata)
options(datadist = 'dd')
```

这一步为后续的 `Predict()` 绘图提供数据分布信息（参考值、分位数等）。**漏写这句会导致 Predict() 报错。**

### 步骤3：拟合样条模型

使用 `rcs()` 函数包裹自变量：

```r
# 线性回归
f <- lm(y ~ rcs(x, 5))

# Logistic 回归
f <- lrm(outcome ~ rcs(age, 5) + sex, data = titanic3)

# Cox 回归
f <- cph(Surv(time, death) ~ rcs(age, 4) + sex, data = data)
```

`rcs()` 中的第二个参数是节点数（knots），默认 4，可设置为 3、4、5 等。

### 步骤4：检验非线性显著性

```r
anova(f)
```

关注输出中 `Nonlinear` 行的 p 值：
- p < 0.05 → 该变量的效应存在显著的非线性成分，适合使用样条回归
- p ≥ 0.05 → 非线性不显著，可用简单的线性项代替

### 步骤5：可视化拟合曲线

**线性回归可视化：**

```r
plot(x, y)
lines(x, fitted(f), col = "red")

# 或使用 ggplot2
ggplot(df, aes(x, y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ rcs(x, 5), se = TRUE, color = "red") +
  theme_bw()
```

**Logistic 回归 / Cox 回归可视化：**

```r
# Logistic: 默认返回 log-odds，加 fun = plogis 返回概率
ggplot(Predict(f, age, sex)) +
  geom_hline(yintercept = 1, color = "grey20", linetype = 2) +
  theme_bw()

# Cox: fun = exp 返回 HR 值
ggplot(Predict(fit, age, fun = exp, ref.zero = TRUE)) +
  geom_hline(yintercept = 1, color = "grey20", linetype = 2) +
  theme_bw()
```

### 步骤6：论文结果报告

论文中可描述为："对年龄变量使用限制性立方样条（第 25、50、75 百分位数为节点）纳入 Cox 比例风险模型，非线性检验结果提示年龄与结局之间存在显著的非线性关联（Nonlinear P = 0.0055）。RCS 曲线显示，当年龄 < 48 岁时 HR < 1 属于保护因素，当年龄 ≥ 48 岁时 HR > 1 转为危险因素。"

## 代码示例

### 示例1：线性回归中的 RCS

```r
rm(list = ls())
x <- 1:100
k <- c(25, 50, 75)
u <- function(x) ifelse(x > 0, x, 0)
x2 <- u(x - k[1])
x3 <- u(x - k[2])
x4 <- u(x - k[3])
set.seed(1)
y <- 0.8 + 1*x + -1.2*x2 + 1.4*x3 + -1.6*x4 + rnorm(100, sd = 2.2)

# 加载R包
library(rms)
library(ggplot2)

# 普通直线回归作为对照
f_lm <- lm(y ~ x)
plot(x, y)
lines(x, fitted(f_lm), col = "red")

# 限制性立方样条回归（5个节点）
f <- lm(y ~ rcs(x, 5))
plot(x, y)
lines(x, fitted(f), col = "red")

# ggplot2 可视化
df.tmp <- data.frame(x = x, y = y)
ggplot(df.tmp, aes(x, y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",
              formula = y ~ rcs(x, 5),
              se = TRUE,
              color = "red") +
  theme_bw()
```

### 示例2：逻辑回归中的 RCS

```r
library(rms)
load(file = "datasets/titanic3.rdata")

# 打包数据
dd <- datadist(titanic3)
options(datadist = 'dd')

# 拟合逻辑回归，年龄用 sqrt 变换后纳入 RCS（5个节点）
f <- lrm(survived ~ rcs(sqrt(age), 5) + sex, data = titanic3)
f
## Logistic Regression Model
## lrm(formula = survived ~ rcs(sqrt(age), 5) + sex, data = titanic3)
## Obs          1046    LR chi2     328.06      R2       0.363    C       0.794
##  0            619    d.f.             5     R2(5,1046)0.266    Dxy     0.588
##  1            427    Pr(> chi2) <0.0001    R2(5,758.1)0.347    gamma   0.592
## max |deriv| 2e-07                            Brier    0.168    tau-a   0.284
##           Coef     S.E.    Wald Z Pr(>|Z|)
## Intercept   3.0936  0.5428   5.70 <0.0001
## age        -0.6383  0.1771  -3.60 0.0003
## age'        1.5544  0.6527   2.38 0.0172
## age''     -12.1583  8.8925  -1.37 0.1715
## age'''     15.8326 16.9397   0.93 0.3500
## sex=male   -2.4944  0.1549 -16.10 <0.0001

# 检验非线性
anova(f)
##                 Wald Statistics          Response: survived
##  Factor     Chi-Square d.f. P
##  age         14.97     4    0.0048
##   Nonlinear  12.65     3    0.0055
##  sex        259.17     1    <.0001
##  TOTAL      265.88     5    <.0001

# age 的 Nonlinear P = 0.0055 < 0.05，存在显著非线性

# 绘制年龄与 OR 的关系曲线
ggplot(Predict(f, age, sex)) +
  geom_hline(yintercept = 1, color = "grey20", linetype = 2) +
  theme_bw()
```

### 示例3：Cox 回归中的 RCS

```r
rm(list = ls())
library(rms)
library(ggplot2)

# 构造模拟生存数据
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male', 'Female'), n, rep = TRUE, prob = c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age - 50) + .8*(sex == 'Female'))
time <- -log(runif(n)) / h
label(time) <- 'Follow-up Time'
death <- ifelse(time <= cens, 1, 0)
time <- pmin(time, cens)
units(time) <- "Year"
data <- data.frame(age, sex, time, death)

# 打包数据
dd <- datadist(data)
options(datadist = 'dd')

# 拟合 Cox 模型（4个节点）
fit <- cph(Surv(time, death) ~ rcs(age, 4) + sex, data = data)
fit
## Cox Proportional Hazards Model
## cph(formula = Surv(time, death) ~ rcs(age, 4) + sex, data = data)
## Obs       1000    LR chi2     78.28     R2       0.083
## Events     183    d.f.            4    R2(4,1000)0.072
## Center -0.2861    Pr(> chi2) 0.0000     R2(4,183)0.334
##          Coef    S.E.   Wald Z Pr(>|Z|)
## age      -0.0173 0.0286 -0.61  0.5443
## age'      0.2040 0.0767  2.66  0.0079
## age''    -0.7500 0.2679 -2.80  0.0051
## sex=Male -0.6445 0.1488 -4.33  <0.0001

# 绘制年龄与 HR 的关系曲线
ggplot(Predict(fit, age, fun = exp, ref.zero = TRUE)) +
  geom_hline(yintercept = 1, color = "grey20", linetype = 2) +
  theme_bw()

# 查找 HR=1 时的年龄
# Predict(fit, age, fun = exp, ref.zero = T)
## 101 48.60445 Male 0.9818776
## 102 48.89330 Male 1.0087001  ← HR 跨越 1 的年龄约为 48 岁

# 重新设定 HR=1 的参考点为 48 岁，绘制带分界线的图
dd$limits$age[2] <- 48
fit <- update(fit)
ggplot(Predict(fit, age, fun = exp, ref.zero = TRUE)) +
  geom_hline(yintercept = 1, color = "steelblue", linetype = 2, linewidth = 1.2) +
  geom_vline(xintercept = 48, color = "red", linetype = 2, linewidth = 1.2) +
  theme_classic()
```

## 结果解读指南

**anova() 输出解读：**

| 输出项 | 含义 |
|--------|------|
| Factor | 模型中的自变量名称，如 age、sex |
| Chi-Square | Wald χ² 统计量，值越大事影响越大 |
| d.f. | 自由度。连续变量使用 rcs() 时，自由度 = 节点数 − 1 |
| P | 该变量的总体显著性 |
| Nonlinear 行 | 最重要的输出：检验该变量在扣除线性效应后是否仍存在显著的非线性成分 |

**解读要点：**

- `age` 行：年龄的总体效应检验，`d.f. = 4` 表示使用了 4 个相关的模型自由度（对应 5 个节点）
- `Nonlinear` 行（`d.f. = 3`）：年龄去掉 1 个线性自由度后的非线性效应检验
  - **Nonlinear P < 0.05**：该变量存在显著的非线性关系，支持使用样条回归
  - **Nonlinear P ≥ 0.05**：非线性不显著，使用普通线性项即可，样条回归无必要

**lrm() / cph() 系数输出解读：**

样条回归的输出中可以看到 `age'`, `age''`, `age'''` 等带撇号的系数项，这些都是样条基函数的系数。**不直接解读这些系数的数值**——它们只是构造曲线的"积木块"，单个系数没有临床意义。应通过 `Predict()` 绘图看整体曲线形状。

**论文报告要点：**

- 说明使用的节点数量和选取方法（如"取第 5、35、65、95 百分位数"或"取默认的 4 个节点"）
- 报告 Nonlinear P 值，判定是否存在非线性关系
- 描述曲线的走向：在哪些范围内 HR/OR > 1（危险因素），哪些范围内 < 1（保护因素）
- 附上带置信区间阴影的 RCS 拟合曲线图

## 常见问题与注意事项

**Q1：`rcs(x, k)` 中的 k 应该如何选择？**

Harrell 的建议：默认 4 个节点，在 3–5 个之间选择即可。一般可以分别尝试 k = 3, 4, 5，比较 AIC 值，取 AIC 最小者。如果有明确的先验知识（如已知某年龄是重要拐点），可将该值指定为节点之一。

可指定节点位置：

```r
rcs(x, c(20, 40, 60))  # 手动指定 3 个节点位置
rcs(x, 4)               # 4 个节点，位置自动按分位数选取
```

**Q2：为什么使用 rms 包前必须先执行 `datadist()` + `options()`？**

`rms` 包的 `Predict()` 函数需要知道每个变量的数据分布信息才能生成合理的预测值网格。如果你忘了执行这两行代码直接调用 `Predict()`，会收到错误提示。**务必记住：每次换数据集后都需要重新 `datadist()`。**

**Q3：`rcs()` 和 `splines::ns()` 有什么区别？**

| 方面 | `rms::rcs()` | `splines::ns()` |
|------|--------------|-----------------|
| 全称 | Restricted Cubic Splines | Natural Splines |
| 包 | rms（Regression Modeling Strategies） | splines（base R 自带） |
| 本质 | 两者都是自然三次样条（在两端约束为线性），概念上等价 |
| 优势 | 与 `lrm()` / `cph()` / `Predict()` 深度集成，绘图方便 | 无需 `datadist()` 打包，轻量级 |
| 选择 | 若使用 rms 系的建模函数（`lrm`, `cph`）同时需要绘制预测曲线，用 `rcs()` | 若仅需简单地在线性模型中加入样条项，用 `ns()` |

**Q4：多项式回归 vs 样条回归，什么时候选哪个？**

- 多项式回归：适用于简单且平滑的非线性模式（U型、温和的 S 型），公式直观可解释，但两端容易波动失真
- 样条回归：适用于更复杂、局部波动的非线性模式，两端被约束为线性所以更稳定，是文献中最常用的非线性建模方法

一般来说，样条回归比多项式回归更受推荐，尤其在医学文献中 RCS 已成为非线性建模的标准方法。

**Q5：Nonlinear P > 0.05，还可用样条回归吗？**

如果 Nonlinear P > 0.05，说明非线性成分不显著。此时：可以保留样条项用于图形展示，但更推荐简化为普通线性项（去掉 `rcs()` 函数包裹），因为线性模型更简洁。论文中可表述为"非线性检验结果提示不存在显著的非线性关系（Nonlinear P = x.xx），故使用线性项纳入模型。"

**Q6：如何找到 HR/OR 恰好等于 1 时对应的 x 值（保护/危险因素分界点）？**

使用 `Predict()` 将模型结果展开为数据框，找到 `yhat` 列（或 `exp(yhat)` 列）从 < 1 跳到 > 1 区间对应的 x 值，然后通过 `dd$limits$变量名[2] <- 该值` 重置参考点，`update()` 模型后重新绘图，此时 HR=1 的位置将自动移到该 x 值处。

**Q7：SPSS 中能做 RCS 吗？**

SPSS 标准界面不支持 RCS，需要通过安装扩展或编写 SPSS 语法实现。相比之下 R 的 `rms` 包是 RCS 分析的业界标准工具，因此文献中绝大多数 RCS 分析使用的都是 R 语言。

**Q8：节点数的自由度如何计算？**

节点数为 k 时，样条项的自由度（d.f.）= k − 1。其中线性分量占 1 个自由度，非线性分量占 k − 2 个自由度。这也是 `anova()` 中 `Nonlinear` 行显示 `d.f. = k − 2` 的原因（如 5 个节点→ d.f. 总量 = 4，非线性 = 3）。

**进一步阅读：**

- 聂博士的 RCS 系列合集（10+ 篇高质量教程）：微信搜索「RCS系列合集」
- F. Harrell. *Regression Modeling Strategies: With Applications to Linear Models, Logistic and Ordinal Regression, and Survival Analysis* (2015, Springer)
