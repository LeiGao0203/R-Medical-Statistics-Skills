---
name: medical-stat-polynomial
description: "R语言医学统计：多项式非线性拟合。使用多项式回归捕捉变量间的非线性关系，包括二次项、三次项拟合及模型比较。TRIGGER when user mentions 多项式回归、非线性关系、二次项、polynomial regression、曲线拟合，or asks about modeling non-linear relationships with polynomial terms. SKIP for 样条回归（RCS）、广义加性模型(GAM)。"
---

# 多项式拟合 (Polynomial Fitting)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用多项式拟合的典型场景：**

- 散点图显示变量间存在明显的曲线关系（U型、倒U型、S型等），用直线回归拟合效果差
- 需要捕捉自变量与因变量之间的非线性趋势，且对函数形式有简单的多项式假设
- 剂量-反应关系中存在先升后降或加速度变化的趋势
- 时间序列数据中的长期趋势拟合（如人口增长、疾病发病率变化）
- 作为基准模型，与更复杂的非线性方法（样条、GAM）进行对比

**不使用多项式拟合的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 非线性模式非常复杂，局部波动大 | 限制性立方样条（RCS）回归 |
| 需要灵活的非参数平滑 | 广义加性模型（GAM, `mgcv::gam()`） |
| 关系近似为直线 | 简单线性回归（`lm(y ~ x)`） |
| 变量间存在已知的拐点/阈值 | 分段回归 |
| 多项式拟合线在两端出现失真翘起 | 样条回归（自然样条/限制性立方样条） |

**医学研究常见应用：**

- 人口增长率随时间的变化（本文示例）
- 某生化指标与疾病风险的U型关系
- 药物剂量与疗效的非线性量效关系
- BMI与死亡率之间的J型/U型关联

## 前置条件

**R 包安装：**

```r
install.packages(c("car", "ggplot2"))
```

**数据格式要求：**

- 自变量 `x`：连续型数值变量（如年份、剂量、年龄）
- 因变量 `y`：连续型数值变量（如人口数、血压、疗效评分）
- 数据框格式，每行一个观测

**统计假设：**

1. **线性设定的可推广**：残差仍需满足正态性、独立性和等方差性（与线性回归一致的假定）
2. **无过度拟合**：多项式的次数不宜过高，通常 ≤ 3 次，过高会导致数据两端的拟合线异常波动
3. **自变量取值范围内拟合**：多项式外推预测非常不可靠——在原数据范围之外，拟合曲线可能迅速偏离实际

## 方法选择决策树

```
你的数据情况 →
├── 散点图显示简单U型或倒U型曲线 → 二次项多项式 (y ~ x + I(x^2))
├── 散点图显示S型或更复杂的单一方向弯曲 → 三次项多项式 (y ~ x + I(x^2) + I(x^3))
├── 不确定用几次项 →
│   ├── 依次拟合1次、2次、3次……用 anova() 做似然比检验
│   ├── 后一项的 p ≥ 0.05 → 选前一个模型即可
│   └── 后一项的 p < 0.05 → 说明需要更高次项
├── 自变量数量多或公式写法繁琐 → 使用 poly(x, degree) 简化写法
├── 拟合线在数据两端出现不自然的上翘/下坠 → 改用限制性立方样条（RCS）回归
└── 用于逻辑回归或Cox回归 → glm() 或 coxph() 中用法与 lm() 完全一致
```

## 标准工作流

### 步骤1：数据准备与探索

加载数据，绘制散点图判断是否存在非线性趋势：

```r
library(car)
data("USPop")
psych::headTail(USPop)

plot(population ~ year, data = USPop)
```

若散点图明显偏离直线，继续步骤2。

### 步骤2：拟合线性回归作为基准

```r
f <- lm(population ~ year, data = USPop)

plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f), col = "blue")
```

直观判断线性拟合是否欠佳，同时留下线性模型供后续比较。

### 步骤3：拟合多项式回归

逐步增加多项式次数：

```r
# 二次项
f1 <- lm(population ~ year + I(year^2), data = USPop)
plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f1))

# 三次项
f2 <- lm(population ~ year + I(year^2) + I(year^3), data = USPop)
plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f2))
```

### 步骤4：模型比较（似然比检验）

```r
# 线性 vs 二次项
anova(f, f1)

# 二次项 vs 三次项
anova(f1, f2)
```

比较依据：若 ANOVA 的 p 值 < 0.05，说明增加该次项后模型拟合显著改善；若 p ≥ 0.05，则当前次数已够，无需继续升高。

### 步骤5：选择最优模型并可视化

参见下方「代码示例」中的 ggplot2 部分，`geom_smooth()` 配合 `se = TRUE` 可同时绘制拟合曲线和置信区间。

### 步骤6：论文报告

论文中可描述为："通过散点图观察，变量间呈非线性趋势。采用多项式回归进行拟合，经似然比检验确认，二次项模型较线性模型显著改善（F = 1408.1, p < 0.001），而三次项较二次项无显著改善（F = 3.39, p = 0.082），故最终选用二次多项式回归模型。"

## 代码示例

完整可运行的代码，加载数据、拟合模型、模型比较、可视化：

```r
# 加载包和数据
library(car)
library(ggplot2)
data("USPop")

# 散点图探索
plot(population ~ year, data = USPop)

# 基准：线性回归
f <- lm(population ~ year, data = USPop)
plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f), col = "blue")

# 二次项多项式
f1 <- lm(population ~ year + I(year^2), data = USPop)
plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f1))

# 三次项多项式
f2 <- lm(population ~ year + I(year^2) + I(year^3), data = USPop)
plot(population ~ year, data = USPop)
lines(USPop$year, fitted(f2))

# 似然比检验
anova(f, f1)   # 线性 vs 二次项
## Model 1: population ~ year
## Model 2: population ~ year + I(year^2)
##   Res.Df     RSS Df Sum of Sq      F    Pr(>F)
## 1     20 12819.0
## 2     19   170.7  1     12648 1408.1 < 2.2e-16 ***

anova(f1, f2)  # 二次项 vs 三次项
## Model 1: population ~ year + I(year^2)
## Model 2: population ~ year + I(year^2) + I(year^3)
##   Res.Df    RSS Df Sum of Sq      F  Pr(>F)
## 1     19 170.66
## 2     18 143.64  1    27.027 3.3868 0.08227 .

# ggplot2 可视化
df.tmp <- data.frame(x = USPop$year, y = USPop$population)
ggplot(df.tmp, aes(x, y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),
              color = "red", linewidth = 1, se = TRUE) +
  theme_bw()
```

**简便写法 —— 使用 `poly()` 函数：**

对于高次项，`I(x^2) + I(x^3) + I(x^4) + ...` 写法繁琐。用 `poly()` 可极大简化：

```r
# 等价于 x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6)
f.6 <- lm(y ~ poly(x, 6))

# 可视化
plot(x, y)
lines(x, fitted(f.6))
```

注意：`poly()` 默认生成**正交多项式**，系数解释与 `I(x^2)` 不同，但拟合值和预测值完全一致。

**扩展至逻辑回归和 Cox 回归：**

```r
# 逻辑回归中的多项式项
glm(y ~ x + I(x^2), family = binomial, data = mydata)

# Cox 回归中的多项式项
library(survival)
coxph(Surv(time, status) ~ x + I(x^2), data = mydata)
```

公式写法与 `lm()` 完全相同，只需替换建模函数即可。

## 结果解读指南

**anova() 似然比检验输出解读：**

| 输出项 | 含义 |
|--------|------|
| Res.Df | 残差自由度，值越小说明模型使用的参数越多 |
| RSS | 残差平方和，越小说明模型拟合越好 |
| Df | 两模型间自由度的差值（增加的参数个数） |
| Sum of Sq | 增加该次项后 RSS 减少的量 |
| F | F 统计量，衡量模型改善是否显著（Sum of Sq 与 RSS 的比值经自由度校正） |
| Pr(>F) | p 值，< 0.05 说明增加该项有统计学意义 |

**如何确定最佳次数：**

1. 从低次到高次逐步使用 `anova(low, high)` 比较嵌套模型
2. 当 p ≥ 0.05 时停止，取前一个显著改善的模型
3. 多数医学应用中，2 次或 3 次已足够，一般不超过 3 次

**论文报告要点：**

- 描述判断非线性关系的依据（散点图趋势、似然比检验结果）
- 报告最终选用的多项式次数
- 给出模型比较的 F 值和 p 值
- 附上拟合曲线图（含原始数据散点和置信区间）

## 常见问题与注意事项

**Q1：`I(year^2)` 和 `poly(year, 2)` 有什么区别？**

| 写法 | 含义 | 系数解释 |
|------|------|----------|
| `y ~ x + I(x^2)` | 原始多项式：x 和 x² 作为两个普通变量纳入模型 | 系数直接对应原始 x 的单位变化效应 |
| `y ~ poly(x, 2)` | 正交多项式：将 x 的多项式转化为互不相关的正交基 | 系数不可直接解读为原始尺度的效应，但拟合值相同 |

建议：需要直观解释系数时用 `I(x^2)` 写法；仅需拟合线或预测时用 `poly()` 更便捷。

**Q2：高次多项式有什么风险？**

- **边缘效应**：x 取值两端的拟合线容易出现不自然的上翘或下坠（Runge 现象）
- **过度拟合**：次数过高时拟合曲线穿过每一个数据点，但泛化能力差
- **多重共线性**：`x`、`x²`、`x³` 之间高度相关，导致系数估计不稳定。`poly()` 的正交多项式可解决此问题

**Q3：多项式回归能用于 logistic 或 Cox 回归吗？**

可以。只需把 `lm()` 换成 `glm(family = binomial)` 或 `coxph()`，公式部分写法完全一致。

**Q4：如何判断是否真的存在非线性关系？**

- 先绘制散点图并用 `geom_smooth(method = "loess")` 观察局部趋势
- 对线性模型和二次项模型做 `anova()` 比较
- 多项式二次项系数的 t 检验 p < 0.05 也提示存在非线性成分

**Q5：多项式次数最多取到多少？**

没有硬性上限，但在医学统计实践中，通常 ≤ 3 次。超过 3 次时建议改用限制性立方样条（RCS）或 GAM，它们在保证灵活性的同时更稳定。

**Q6：SPSS 中也能做多项式回归，与 R 有什么差异？**

SPSS 的曲线估计中内置了线性、二次、三次等模型，操作更图形化，但不支持 R 中灵活的 `poly()` 正交多项式写法。R 的 `anova()` 嵌套模型比较也更直接。两者的拟合结果是等价的。
