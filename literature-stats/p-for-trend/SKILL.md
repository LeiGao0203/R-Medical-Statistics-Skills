---
name: medical-stat-p-for-trend
description: "R语言医学统计：趋势检验、交互作用检验和每标准差效应分析。涵盖有序变量趋势性检验（p for trend）、交互作用分析（p for interaction）及per 1-SD增量效应报告。TRIGGER when user mentions p for trend、p for interaction、per 1 SD、趋势P值、交互作用P值、per SD效应，or asks about trend tests and interaction in regression. SKIP for Cochran-Armitage趋势检验、亚组分析的森林图。"
---

# p-for-trend / p-for-interaction / per-1-SD 分析
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

本方法适用于回归分析（逻辑回归、Cox回归、线性回归等）中评估自变量与因变量的剂量反应关系及变量间交互作用。典型使用场景：

- **p for trend**：连续型变量按分位数/临床切点分箱后，检验各等级间因变量是否存在线性趋势。常用于暴露-反应关系研究、剂量-效应关系评估
- **p for interaction**：检验两个变量之间是否存在交互效应（效应修饰），常用于分层分析、敏感性分析
- **per 1 SD**：报告自变量每增加1个标准差的效应量，使不同量纲变量的效应可比

不使用本方法的场景：
- 分类变量（二分类/多分类）的趋势检验请使用 **Cochran-Armitage检验**
- 亚组分析的森林图绘制不属于本方法范畴
- 本方法假设线性趋势；若趋势为非线性，请使用多项式拟合或样条回归

## 前置条件

**R包依赖：**

```r
install.packages("foreign")   # 读取SPSS数据
install.packages("broom")     # 整理模型输出为tidy格式
install.packages("lmtest")    # 似然比检验
```

**数据要求：**

- 自变量：连续型变量或其分箱后的有序分类变量（作为数值型纳入即可计算p for trend）
- 因变量：二分类（逻辑回归）、生存时间+状态（Cox回归）或连续型（线性回归）
- 用于p for interaction时，两变量乘积项应非完全共线

**统计前提：**

- p for trend 假设因变量随自变量等级线性变化；若实际为U型等非线性关系，p for trend 可能不显著或产生误导
- p for interaction 需先确定主效应模型正确设定，交互项检验基于似然比检验或Wald检验

## 方法选择决策树

```
你的分析目标 →
├── 评估有序分类变量的线性趋势（各等级OR/HR并报告趋势P值）
│   ├── 自变量为连续型 → 先分箱为有序分类变量 → 作为因子纳入回归得各级OR/CI → 作为数值型纳入得p for trend
│   └── 自变量已是等级变量 → 作为因子纳入得各级OR/CI → 作为数值型纳入得p for trend
├── 检验两变量间的交互作用
│   ├── 两变量均为数值型/二分类 → 新建乘积项纳入模型，看交互项P值（方法1）
│   ├── 两变量均为数值型/二分类 → 构建含/不含交互项的两模型，似然比检验（方法2）
│   └── 多分类变量（哑变量化后）→ 必须用似然比检验比较两模型，不可直接看单个交互项P值
└── 报告每标准差效应（per 1 SD）→ scale()标准化自变量 → 纳入回归 → 报告效应量及95%CI
```

## 标准工作流

### 步骤1：数据准备与探索
- 读入数据（SPSS `.sav`、CSV等格式）
- 检查变量类型，确认自变量是连续型还是已分箱
- 若需分箱，按临床切点或分位数划分，并检查各箱样本量（避免某箱样本过少）

### 步骤2：p for trend — 连续型变量分箱后趋势检验
- 将连续型变量分箱后创建为因子变量
- 以因子变量纳入回归模型，得到每级相对于参考水平的OR/HR及95%CI
- 将有序等级变量以数值型纳入同一回归模型，该项的P值即为p for trend
- 注意：p for trend 检验的是线性趋势，数值型编码应当等距（如1,2,3,4）

### 步骤3：p for interaction — 交互作用检验
- 方法1（乘积项法）：新建两变量的乘积项，纳入回归模型，查看乘积项Wald检验P值
- 方法2（似然比检验法）：分别拟合无交互项和有交互项的两个模型，用 `lrtest()` 比较
- 多分类变量必须使用方法2

### 步骤4：per 1 SD — 标准化效应
- 使用 `scale()` 对自变量进行标准化（均数=0，标准差=1）
- 将标准化后的变量纳入回归模型，结果即为每增加1个SD的效应量

### 步骤5：结果报告（论文中的统计描述）
- p for trend：报告中写明趋势P值，附各级OR/CI表格
- p for interaction：报告交互P值，若显著需分层报告各亚组效应
- per 1 SD：报告"per 1-SD increase, OR = x.xx (95%CI: x.xx - x.xx), P = x.xxx"
- 所有方法均应说明校正了哪些协变量

## 代码示例

```r
# 加载包
library(foreign)
library(broom)
library(lmtest)

# 读取数据（孙振球版医学统计学第4版 例16-02）
df16_2 <- foreign::read.spss("datasets/例16-02.sav",
                             to.data.frame = TRUE,
                             use.value.labels = FALSE,
                             reencode = "utf-8")
# 变量说明：
# x1: 年龄（1=<45, 2=45-55, 3=55-65, 4=>65）
# x2: 高血压病史（1=有, 0=无）
# x7: BMI（1=<24, 2=24-26, 3=>26）
# y: 冠心病（1=是, 0=否）

# ======== p for trend ========
# 将x1作为数值型纳入逻辑回归，得到p for trend
f_trend <- glm(y ~ x1 + x2,
               data = df16_2,
               family = binomial())
broom::tidy(f_trend)
## # A tibble: 3 × 5
##   term        estimate std.error statistic p.value
##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
## 1 (Intercept)   -2.22      1.03      -2.15  0.0313
## 2 x1             0.712     0.423      1.68  0.0928  ← p for trend
## 3 x2             1.08      0.625      1.73  0.0840

# 将x1变为因子，得到各级OR和95%CI
df16_2$x1.f <- factor(df16_2$x1)
f_or <- glm(y ~ x1.f + x2,
            data = df16_2,
            family = binomial())
broom::tidy(f_or, conf.int = TRUE, exponentiate = TRUE)
## # A tibble: 5 × 7
##   term        estimate std.error statistic p.value conf.low conf.high
##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
## 1 (Intercept)    0.200     1.10     -1.47   0.142    0.0104      1.24
## 2 x1.f2          2.32      1.19      0.704  0.481    0.289      49.3
## 3 x1.f3          4.48      1.26      1.19   0.233    0.485     102.
## 4 x1.f4          9.42      1.63      1.38   0.169    0.508     438.
## 5 x2             2.94      0.639     1.69   0.0918   0.854      10.7

# ======== p for interaction ========
# 方法1：新建乘积项，看交互项P值
df16_2$x17 <- df16_2$x1 * df16_2$x7  # 年龄 × BMI交互项
f_int1 <- glm(y ~ x1 + x7 + x17,
              family = binomial(),
              data = df16_2)
summary(f_int1)
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## x17           0.9249     0.7489   1.235    0.217  ← p for interaction

# 方法2：似然比检验比较两个模型
f_noint <- glm(y ~ x1 + x7,
               family = binomial(),
               data = df16_2)
f_withint <- glm(y ~ x1 + x7 + x17,
                 family = binomial(),
                 data = df16_2)
lmtest::lrtest(f_noint, f_withint)
## Likelihood ratio test
## Model 1: y ~ x1 + x7
## Model 2: y ~ x1 + x7 + x17
##   #Df  LogLik Df  Chisq Pr(>Chisq)
## 1   3 -32.216
## 2   4 -31.254  1 1.9238     0.1654  ← p for interaction

# ======== per 1 SD ========
df16_2$weight <- rnorm(54, 70, 11)
df16_2$weight.scaled <- scale(df16_2$weight)  # 均值0，标准差1

f_sd <- glm(y ~ weight.scaled, data = df16_2)
broom::tidy(f_sd, conf.int = TRUE, exponentiate = TRUE)
## # A tibble: 2 × 7
##   term          estimate std.error statistic       p.value conf.low conf.high
##   <chr>            <dbl>     <dbl>     <dbl>         <dbl>    <dbl>     <dbl>
## 1 (Intercept)       1.62    0.0689     6.98  0.00000000525    1.41       1.85
## 2 weight.scaled     1.05    0.0696     0.727 0.470            0.918      1.21
##   ↑ 即 weight 每增加1个SD，冠心病OR = 1.05 (95%CI: 0.92-1.21)
```

## 结果解读指南

**p for trend 结果解读：**
- x1 的 P 值（0.0928）即为校正x2后的趋势检验P值
- P < 0.05 说明随着年龄等级升高，冠心病风险存在线性递增趋势（剂量-反应关系）
- 同时查看各级OR值：x1.f2→2.32, x1.f3→4.48, x1.f4→9.42，OR逐级增大也佐证趋势存在
- 注意：p for trend 不显著（P > 0.05）可能因为：(1) 确实无趋势；(2) 趋势非线性；(3) 样本量不足

**p for interaction 结果解读：**
- 方法1：x17（乘积项）P=0.217 > 0.05，年龄与BMI之间不存在统计学意义的交互作用
- 方法2：似然比检验 P=0.1654 > 0.05，加入交互项未能显著改善模型拟合，结论与方法1一致
- 两种方法通常结果相近；多分类变量仅能使用方法2
- P < 0.05 时说明交互作用存在，需分层报告各亚组的效应量，不可混用总效应

**per 1 SD 结果解读：**
- weight.scaled 系数的指数化值即 per 1 SD 的OR
- 上例：weight 每增加1个SD（约11 kg），冠心病OR = 1.05（95%CI: 0.92-1.21），P = 0.470
- 标准化后不同量纲变量的效应量可直接比较大小

**论文报告示例：**
> 表X显示各年龄组冠心病风险的OR值及p for trend。以<45岁组为参考，45-55岁组OR=2.32 (95%CI: 0.29-49.3)，55-65岁组OR=4.48 (95%CI: 0.48-102.3)，>65岁组OR=9.42 (95%CI: 0.51-437.5)。趋势检验结果p for trend = 0.093。

## 常见问题与注意事项

**Q1: p for trend 和 Cochran-Armitage 检验有何区别？**
Cochran-Armitage 专门针对分类变量（2×K列联表）的趋势检验，不需要回归模型。p for trend 基于回归模型，可校正协变量，适用于连续型自变量分箱后的趋势分析。两者应用场景不同，不可混用。

**Q2: 分箱时类别编码是否必须等距？**
是的。p for trend 假设线性趋势，因此数值型编码应当等距（如1,2,3,4）。若各组间实际间距不等（如0-10, 10-30, 30-100），建议使用各组中位数重新编码。

**Q3: 似然比检验（方法2）和乘积项直接看P值（方法1）结果不一致怎么办？**
两方法理论基础不同：方法1基于Wald检验，方法2基于似然比检验。对小样本或非线性情况，似然比检验更为稳健。多分类变量交互时仅能使用方法2。实际应用中建议两种方法都做，结果一致时更有说服力。

**Q4: per 1 SD 和直接使用原始单位（per unit）哪个更好？**
per 1 SD 使不同量纲变量的效应量可横向比较，在涉及多种暴露因素的论文中常用。但需要注意：SD受样本影响，不同研究SD不同导致效应量不可直接跨研究比较。建议同时报告原始单位效应和per SD效应。

**Q5: 交互作用P值不显著，是否还需做分层分析？**
交互P > 0.05 说明无统计学意义的效应修饰，此时各亚组效应方向/大小一致，通常不需要分层报告。但如果领域惯例或审稿人要求，可补充分层分析作为敏感性分析。

**Q6: 如何在Cox回归中实现这三类分析？**
代码结构完全相同，仅需将 `glm(..., family = binomial())` 替换为 `coxph()`（survival包）。p for trend、p for interaction、per 1 SD 均直接适用。
