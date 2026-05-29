---
name: medical-stat-logistic-reg
description: "R语言医学统计：Logistic回归。涵盖二分类Logistic回归、多重Logistic回归、有序Logistic回归、条件Logistic回归、OR值解读、逐步回归、模型评估（似然比检验、Hosmer-Lemeshow检验、C统计量、伪R²）。TRIGGER when user mentions Logistic回归、逻辑回归、logit、OR值、odds ratio、二分类结局、多分类结局、配对病例对照，or asks about modeling binary or categorical outcomes with predictors. SKIP for 多元线性回归（连续因变量）、Cox回归（生存数据）、Poisson回归（计数数据）、对数线性模型（列联表）。"
---

# Logistic回归 (Logistic Regression)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**什么时候用：**

- 因变量为**二分类变量**（如患病/未患病、死亡/存活），需要探索多个危险因素的影响 → 二项Logistic回归
- 因变量为**无序多分类变量**（如获取健康知识的途径：传统媒体/网络/社区宣传），需要分析影响因素 → 多项Logistic回归（multinomial）
- 因变量为**有序等级变量**（如疗效：无效/有效/痊愈），需要同时利用等级顺序信息 → 有序Logistic回归（ordinal）
- 数据来自**配对病例-对照研究**（如1:M匹配），需要控制匹配因素 → 条件Logistic回归
- 典型医学场景：筛选疾病危险因素、构建诊断预测模型、评估治疗方式的疗效差异、分析1:M匹配的病例-对照资料

**什么时候不用：**

- 因变量是连续型数值变量 → 使用多元线性回归
- 因变量是生存时间（含删失） → 使用Cox比例风险回归
- 因变量是计数数据（如住院次数、发病次数） → 使用Poisson回归或负二项回归
- 仅做分类资料的关联性检验（无多个自变量调整） → 使用卡方检验
- 需要分析名义因变量和多个分类自变量之间的复杂关联结构 → 考虑对数线性模型

## 前置条件

**R包：**

```r
# 基础安装已包含 glm()、step()
# 二项Logistic回归无需额外包
# 多项Logistic回归
install.packages("nnet")
# 有序Logistic回归
install.packages("MASS")
# 条件Logistic回归
install.packages("survival")
# 伪R²
install.packages("DescTools")
# 整洁输出
install.packages("broom")
# 平行线检验（有序回归）
install.packages("brant")
```

**数据格式要求：**

- 因变量：二分类变量用 `factor` 类型（两水平），多分类变量用 `factor` 类型（三水平及以上），有序变量用 `ordered = TRUE` 的 `factor` 类型
- 自变量：数值型变量可直接纳入，分类变量必须先转为 `factor` 类型（R会自动进行哑变量编码）
- 条件Logistic回归：数据需包含匹配对子的编号列，用于 `strata()`

**关键统计假设：**

1. **线性假设**（对连续自变量）：连续自变量与 logit(P) 之间呈线性关系。可先用 `boxplot(连续变量 ~ 因变量)` 初步判断
2. **无多重共线性**：自变量间不应高度相关。可使用 `vif()` 检验（来自 `car` 包）
3. **样本量充足**：一般要求事件数（较少的那一类）至少是自变量个数的 10-15 倍
4. **有序Logistic回归特需**：需满足**平行线假设**（proportional odds assumption），即自变量的回归系数在不同等级的切点处相等。使用 `brant::brant()` 检验
5. **条件Logistic回归**：匹配因素（如年龄、性别）已被对子内部自动控制，无需作为自变量纳入模型

## 方法选择决策树

```
你的因变量类型 →
├── 二分类（是/否、患病/未患病）
│   ├── 独立样本 → 二项Logistic回归（glm(y ~ ., family = binomial())）
│   └── 配对病例-对照（1:M匹配） → 条件Logistic回归（clogit(y ~ . + strata(i))）
├── 无序多分类（无等级顺序，如A型/B型/C型）
│   └── 多项Logistic回归（nnet::multinom(y ~ .)）
│       注意：需设置参考类别，R默认以因子第一水平为参考
└── 有序等级（疗效None/Mild/Complete、分期I/II/III）
    └── 有序Logistic回归（MASS::polr(Y ~ ., method = "logistic")）
        使用前需先做平行线检验（brant test），不满足时改用多项Logistic回归
```

## 标准工作流

### 步骤1：数据准备与探索

将分类自变量转为因子型，有序因变量转为有序因子型。使用 `str()` 检查变量类型，使用 `table(y)` 查看因变量分布。对于病例-对照数据，确认匹配对子的标识变量。

```r
# 分类变量转为因子
df[, vars] <- lapply(df[, vars], factor)
# 因变量转为有序因子
df$Y <- factor(df$Y, levels = c(1, 2, 3), labels = c("无效", "有效", "痊愈"), ordered = TRUE)
```

### 步骤2：前提条件检验

- 检查因变量各类别的频数：`table(df$y)` — 确保没有零频的类别以避免完全分离
- 对连续自变量检查与logit的线性关系：可用 `car::boxTidwell()` 或拟合样条后比较
- 检查多重共线性：`car::vif(glm(y ~ ., family = binomial(), data = df))`（VIF > 5 或 10 提示共线性）
- 有序Logistic回归：执行平行线检验 `brant::brant(fit)`

### 步骤3：执行统计分析

根据因变量类型选择对应的回归函数。二项Logistic回归使用 `glm(..., family = binomial())`，可通过 `step(model, direction = "backward")` 进行逐步回归筛选变量。

### 步骤4：结果解读

提取回归系数（β）、标准误、Wald统计量、P值、OR值及其95%可信区间。根据P值确定具有统计学意义的自变量，根据OR值大小和方向判断危险/保护因素。

### 步骤5：结果报告

论文中应报告：每个自变量的 β、SE、Wald χ²（或z值）、P值、OR 和 95% CI。通常制作一个多因素Logistic回归分析结果三线表。还需报告模型的整体显著性检验（似然比检验）和拟合优度指标（AIC、伪R²、Hosmer-Lemeshow检验结果）。

## 代码示例

### 二项Logistic回归（binomial logistic regression）

**示例数据：** 孙振球《医学统计学》第4版例16-2，26例冠心病患者和28例对照者，探讨冠心病危险因素。

变量说明：x1=年龄（1:<45, 2:45-55, 3:55-65, 4:>65），x2=高血压史（1有/0无），x3=高血压家族史，x4=吸烟，x5=高血脂，x6=动物脂肪摄入（0低/1高），x7=BMI（1:<24, 2:24-26, 3:>26），x8=A型性格，y=冠心病（1是/0否）。

```r
# 加载数据
df16_2 <- foreign::read.spss("datasets/例16-02.sav",
                             to.data.frame = TRUE,
                             use.value.labels = FALSE,
                             reencode = "utf-8")

str(df16_2)
## 'data.frame':  54 obs. of  11 variables:
##  $ x1: num  3 2 2 2 3 3 2 3 2 1 ...
##  $ x2: num  1 0 1 0 0 0 0 0 0 0 ...
##  $ x3: num  0 1 0 0 0 1 1 1 0 0 ...
##  ...  (x4~x8类似)
##  $ y : num  0 0 0 0 0 0 0 0 0 0 ...

# 将自变量和因变量转为因子型
# 注意：x1和x7是有序分类变量，其余为无序分类
df16_2[, c(2:10)] <- lapply(df16_2[, c(2:10)], factor)

# 拟合二项Logistic回归
f <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
         data = df16_2,
         family = binomial())

summary(f)
##
## Call:
## glm(formula = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
##     family = binomial(), data = df16_2)
##
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept) -5.46026    2.07370  -2.633  0.00846 **
## x12          0.85285    1.54399   0.552  0.58070
## x13          0.47754    1.59320   0.300  0.76438
## x14          3.44227    2.10985   1.632  0.10278
## x21          1.14905    0.93176   1.233  0.21750
## x31          1.66039    1.16857   1.421  0.15535
## x41          0.85994    1.32437   0.649  0.51613
## x51          0.73600    0.97088   0.758  0.44840
## x61          3.92067    1.57004   2.497  0.01252 *
## x72         -0.03467    1.13363  -0.031  0.97560
## x73         -0.38230    1.61710  -0.236  0.81311
## x81          2.46322    1.10484   2.229  0.02578 *
## ---
##     Null deviance: 74.786  on 53  degrees of freedom
## Residual deviance: 40.028  on 42  degrees of freedom
## AIC: 64.028

# 提取OR值：OR = exp(β)
exp(coef(f))
##  (Intercept)          x12          x13          x14          x21
##  0.004252469  2.346329320  1.612111759 31.257871683  3.155194147
##          x31          x41          x51          x61          x72
##  5.261381340  2.363023282  2.087573511 50.434470096  0.965919321
##          x73          x81
##  0.682290259 11.742555242

# OR值的95%可信区间
exp(confint(f))
##                    2.5 %       97.5 %
## (Intercept)  3.136876e-05    0.1376413
## x12          1.311093e-01   81.8646261
## x13          7.863610e-02   59.4639513
## ...
## x61          3.777465e+00 2159.5535363
## x81          1.666190e+00  148.0206875

# 提取Wald值：Wald = z^2 = (β/SE)^2
summary(f)$coefficients[, 3]^2
## (Intercept)         x12         x13         x14         x21
## 6.933188870 0.305111544 0.089843733 2.661883233 1.520790277
##         x31         x41         x51         x61         x72
## 2.018903576 0.421615676 0.574682148 6.235929079 0.000935592
##         x73         x81
## 0.055890396 4.970577395

# 模型整体检验：与空模型比较
f0 <- glm(y ~ 1, data = df16_2, family = binomial())
anova(f0, f, test = "Chisq")
## Analysis of Deviance Table
## Model 1: y ~ 1
## Model 2: y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
##   Resid. Df Resid. Dev Df Deviance  Pr(>Chi)
## 1        53     74.786
## 2        42     40.028 11   34.758 0.0002716 ***
# P < 0.001，模型整体有统计学意义

# 伪R²
DescTools::PseudoR2(f, which = c("McFadden", "CoxSnell", "Nagelkerke"))
##   McFadden   CoxSnell Nagelkerke
##  0.4647704  0.4746397  0.6331426

# 逐步回归（向后法）
f2 <- step(f, direction = "backward")
## 经过逐步筛选后最终纳入: y ~ x2 + x3 + x6 + x8
summary(f2)
##
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)  -3.0314     0.8965  -3.381 0.000722 ***
## x21           1.4715     0.7656   1.922 0.054617 .
## x31           1.2251     0.7543   1.624 0.104359
## x61           3.6124     1.3391   2.698 0.006985 **
## x81           1.8639     0.8045   2.317 0.020505 *
## ---
## AIC: 57.537
```

### 多项Logistic回归（multinomial logistic regression）

**示例数据：** 例16-5，2个社区314名成人，探索社区和性别对获取健康知识途径的影响。Y=1传统大众传媒，2网络，3社区宣传。

```r
library(nnet)

df <- read.csv("datasets/例16-05.csv", header = TRUE)

# 转为因子型，设置标签
df$X1 <- factor(df$X1, levels = c(0, 1), labels = c("社区1", "社区2"))
df$X2 <- factor(df$X2, levels = c(0, 1), labels = c("男", "女"))
df$Y  <- factor(df$Y,  levels = c(1, 2, 3),
                labels = c("传统大众传媒", "网络", "社区宣传"))

# 多项Logistic回归
fit <- multinom(Y ~ X1 + X2, data = df, model = TRUE)
## converged

summary(fit)
## Coefficients:
##          (Intercept)    X1社区2      X2女
## 网络       0.5484998 -1.3743147 0.4321069
## 社区宣传   0.3940422 -0.9933526 1.2266459
##
## Residual Deviance: 633.1508
## AIC: 645.1508

# 使用broom::tidy获得整洁结果（含P值）
broom::tidy(fit)
## # A tibble: 6 x 6
##   y.level  term        estimate std.error statistic   p.value
##   <chr>    <chr>          <dbl>     <dbl>     <dbl>     <dbl>
## 1 网络     (Intercept)    0.548     0.258      2.12 0.0337
## 2 网络     X1社区2       -1.37      0.320     -4.29 0.0000177
## 3 网络     X2女           0.432     0.327      1.32 0.186
## 4 社区宣传 (Intercept)    0.394     0.257      1.53 0.126
## 5 社区宣传 X1社区2       -0.993     0.295     -3.36 0.000766
## 6 社区宣传 X2女           1.23      0.299      4.10 0.0000413

# 手动计算P值的方式（备选）
z_stats <- summary(fit)$coefficients / summary(fit)$standard.errors
p_values <- (1 - pnorm(abs(z_stats))) * 2

# OR值
exp(coef(fit))
##          (Intercept)   X1社区2     X2女
## 网络        1.730655 0.2530129 1.540500
## 社区宣传    1.482963 0.3703330 3.409774

# OR值95% CI
exp(confint(fit))

# 模型整体检验
fit0 <- multinom(Y ~ 1, data = df, model = TRUE)
anova(fit0, fit)
## Likelihood ratio tests of Multinomial Models
##     Model Resid. df Resid. Dev   Test    Df LR stat.      Pr(Chi)
## 1       1       626   677.2069
## 2 X1 + X2       622   633.1508 1 vs 2     4  44.0561 6.245931e-09

# 预测类别和概率
pred  <- predict(fit, df, type = "class")
prob  <- predict(fit, df, type = "probs")  # 或 fitted(fit)

# 伪R²
DescTools::PseudoR2(fit, which = "all")
##        McFadden        CoxSnell      Nagelkerke
##      0.06505559       0.13090778       0.14803636
```

### 有序Logistic回归（ordinal logistic regression）

**示例数据：** 例16-4，84例患者临床试验，探索性别（X1，男=0/女=1）和治疗方法（X2，传统=0/新型=1）对疗效的影响。Y=1无效，2有效，3痊愈。

```r
library(MASS)

df <- read.csv("datasets/例16-04.csv", header = TRUE)

# 因变量设为有序因子
df$Y  <- factor(df$Y, levels = c(1, 2, 3),
                labels = c("无效", "有效", "痊愈"),
                ordered = TRUE)
df$X1 <- factor(df$X1, levels = c(0, 1), labels = c("男", "女"))
df$X2 <- factor(df$X2, levels = c(0, 1), labels = c("传统疗法", "新型疗法"))

# 有序Logistic回归
fit <- polr(Y ~ X1 + X2, data = df, Hess = TRUE, method = "logistic")
summary(fit)
##
## Coefficients:
##            Value Std. Error t value
## X1女       1.319     0.5381   2.451
## X2新型疗法 1.797     0.4718   3.809
##
## Intercepts:
##           Value  Std. Error t value
## 无效|有效 1.8128 0.5654     3.2061
## 有效|痊愈 2.6672 0.6065     4.3979
##
## Residual Deviance: 150.0294
## AIC: 158.0294

# 手动计算P值
p <- pnorm(abs(coef(summary(fit))[, "t value"]), lower.tail = FALSE) * 2
##         X1女   X2新型疗法    无效|有效    有效|痊愈
## 1.425572e-02 1.392807e-04 1.345300e-03 1.092866e-05

# OR值
exp(coef(fit))
##       X1女 X2新型疗法
##   3.738765   6.033338

# 平行线检验（有序回归的关键前提）
library(brant)
brant::brant(fit)
## --------------------------------------------
## Test for X2  df  probability
## --------------------------------------------
## Omnibus      1.83    2   0.4
## X1女         1.59    1   0.21
## X2新型疗法   0.01    1   0.94
## --------------------------------------------
## H0: Parallel Regression Assumption holds
# P > 0.05，平行线假设成立，可使用有序Logistic回归

# 模型整体检验
fit0 <- polr(Y ~ 1, data = df, Hess = TRUE, method = "logistic")
anova(fit0, fit)
## Likelihood ratio tests of ordinal regression models
##     Model Resid. df Resid. Dev   Test    Df LR stat.     Pr(Chi)
## 1       1        82   169.9159
## 2 X1 + X2        80   150.0294 1 vs 2     2  19.8865 4.80508e-05

# 预测类别和概率
pred <- predict(fit, df, type = "class")
prob <- predict(fit, df, type = "probs")

# 伪R²
DescTools::PseudoR2(fit, which = "all")
##        McFadden        CoxSnell      Nagelkerke
##       0.1170373       0.2108068       0.2429443
```

### 条件Logistic回归（conditional logistic regression）

**示例数据：** 例16-3，用1:2配对研究探讨喉癌发病危险因素，6个危险因素，25对数据（75条记录）。`i` 为对子编号。

```r
library(survival)

df <- foreign::read.spss("datasets/例16-03.sav", to.data.frame = TRUE)
str(df)
## 'data.frame':  75 obs. of  8 variables:
##  $ i : num  1 1 1 2 2 2 ...  # 对子编号
##  $ y : num  1 0 0 1 0 0 ...  # 1=病例, 0=对照
##  $ x1~x6: 6个危险因素（数值型，不需要转因子）

# 条件Logistic回归：strata(i)指定同一对子内的个体为一组
fit <- clogit(y ~ x1 + x2 + x3 + x4 + x5 + x6 + strata(i),
              data = df, method = "exact")

summary(fit)
##
##        coef exp(coef) se(coef)      z Pr(>|z|)
## x1  2.58880  13.31380  2.50172  1.035   0.3008
## x2  1.68796   5.40843  0.68545  2.463   0.0138 *
## x3  2.31944  10.16995  1.26096  1.839   0.0659 .
## x4 -3.88886   0.02047  1.90656 -2.040   0.0414 *
## x5 -0.49102   0.61200  1.19020 -0.413   0.6799
## x6  3.50899  33.41447  2.13723  1.642   0.1006
## ---
## Likelihood ratio test= 42.21  on 6 df,   p=2e-07
## Concordance= 0.91  (se = 0.064 )

# clogit的结果直接输出 coef, exp(coef), se(coef), z, Pr(>|z|)
# 不需要额外计算OR = exp(coef)，因为已包含在输出中
```

## 结果解读指南

**回归系数（β, Estimate）：**

- β 是 logit(P) 的回归系数。β > 0 表示该因素增加事件发生的概率（危险因素）；β < 0 表示降低概率（保护因素）
- Wald统计量 = (β / SE)² = z²，用于检验 β 是否等于 0
- P < 0.05 表示该自变量对因变量的影响具有统计学意义

**OR值（Odds Ratio）：**

- OR = exp(β)，表示该因素每增加一个单位（或与参照组相比）时，事件发生优势的倍数
- OR > 1：危险因素；OR < 1：保护因素；OR = 1：无关联
- 95% CI 跨过 1 表示无统计学意义（P > 0.05）
- 以「x12: OR = 2.35」为例，解读为：45~55岁人群患冠心病的风险是小于45岁人群的2.35倍，但P=0.58无统计学意义

**模型整体评价：**

- `Null deviance` 和 `Residual deviance`：两者差值越大越好，反映模型解释力。似然比检验 P < 0.05 表示模型有统计学意义
- AIC（赤池信息准则）：值越低越好，可用于比较不同模型的拟合效果
- 伪R²（McFadden / CoxSnell / Nagelkerke）：与线性回归R²不同，伪R²通常偏低，0.2~0.4已可接受。McFadden's R² > 0.2 已视为较好
- Concordance（C统计量，条件Logistic回归输出）：0.7-0.8可接受，0.8-0.9优秀，>0.9极好。反映模型区分能力，等价于ROC曲线下面积（AUC）

**哑变量编码的解读：**

R语言默认将因子的第一个水平作为参照组，输出中以 `变量名+水平` 形式出现（如 `x12` 表示 x1=2 组与 x1=1 组比较）。需要注意哪个是参照组。可通过 `relevel(df$x1, ref = "2")` 自定义参照组。

**有序Logistic回归特有输出解读：**

- `Intercepts`（截距/阈值）：`无效|有效` 和 `有效|痊愈` 是两个分割点的截距
- 平行线检验 `brant()`：H0 为回归系数在各级别间相同。P > 0.05 表示满足平行线假设，可以使用有序Logistic回归；P < 0.05 则应改用多项Logistic回归

## 常见问题与注意事项

**Q1: R和SPSS的参考水平默认设置不同，如何对齐？**

SPSS默认以**最后一个水平**为参考，R默认以**第一个水平**为参考。要让R与SPSS一致，使用 `relevel()` 设置参照：

```r
df$x1 <- relevel(df$x1, ref = "4")  # 以第4组为参照
```

或因变量使用 `df$y <- relevel(df$y, ref = "1")` 来确保1是参照。

**Q2: 自变量是连续型数值变量如何处理？**

连续型数值变量不应转为因子，直接纳入模型即可。但需检查线性假设（continuous predictor与logit(P)之间应为线性关系）。可在模型中加入平方项或使用限制性立方样条（RCS）处理非线性关系。

**Q3: 逐步回归应该用向前、向后还是步进法？**

向后法（`direction = "backward"`）优于向前法，因为向前法可能遗漏重要交互作用。步进法（`direction = "both"`）在医学研究中较常用。但现代医学统计越来越强调**基于专业知识的变量选择**而非纯数据驱动的自动筛选。逐步回归结果的OR值可能存在偏倚，应在论文中说明。

**Q4: 有序Logistic回归不满足平行线假设怎么办？**

改用**多项Logistic回归**（`nnet::multinom()`），放弃有序信息，将因变量视为无序多分类处理。或将因变量合并为二分（如无效 vs 有效+痊愈）改用二项Logistic回归。

**Q5: 出现完全分离（complete separation）或拟完全分离怎么办？**

当某个自变量的取值可以完美预测因变量时（如所有x=1的人y都为1），会导致β估计趋于无穷大、标准误异常大。解决方法：
- 使用 `brglm2::brglm()` 进行Firth's penalized likelihood估计
- 或使用 `glm(..., method = "brglmFit")` 作为替代
- 检查数据，考虑合并稀疏类别或增加样本量

**Q6: 如何报告Logistic回归结果？**

论文中最少需报告：OR值、95% CI、P值。规范做法为制作表格，包含：
- 每个自变量的 β（可省略）、SE、Wald χ²、P值、OR（95% CI）
- 表格脚注说明：模型整体检验的 χ² 值和 P 值、伪R²、模型纳入的变量
- 若做逐步回归，需说明纳入和排除标准（如 α_入=0.05, α_出=0.10）

**Q7: 条件Logistic回归与普通Logistic回归的核心区别？**

条件Logistic回归通过 `strata()` 将同一匹配组（对子）内的个体归为一层，每一层内单独估计条件似然，从而自动消除匹配因素（如年龄、性别）的混杂效应。这些匹配因素**不应**再作为自变量放入模型中。匹配比例可以是1:1、1:M或M:N。

**Q8: 自变量赋值的影响？**

二分类变量直接赋0/1或1/2均可，但0/1更直观（系数即为有/无的效果）。对于有序自变量（如病情严重程度），即使编码为1/2/3/4，如果作为数值变量纳入，模型将假定其效应为线性递增。如希望各等级独立估计效应，应将其转为因子型（设置哑变量）。
