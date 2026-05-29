---
name: medical-stat-poisson-nb-reg
description: "R语言医学统计：泊松回归与负二项回归。用于分析计数资料（发病率、发生次数）的影响因素，负二项回归适用于存在过离散的计数数据。TRIGGER when user mentions 泊松回归、负二项回归、Poisson回归、Negative Binomial回归、计数资料、count data、过离散（overdispersion）、发病密度、incidence density、偏移量（offset）、IRR，or asks about modeling count outcomes. SKIP for Logistic回归（二分类结局）、线性回归（连续型结局）、Cox回归（生存分析）。"
---

# 泊松回归和负二项回归 (Poisson and Negative Binomial Regression)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用场景：**
- 分析计数型结局变量（count data）的影响因素，如某疾病发病次数、住院次数、死亡人数等
- 罕见事件发生率研究（恶性肿瘤、先天性疾病等）
- 当因变量为发病密度（incidence density）时，需要设置偏移量（offset）调整观察人年
- 数据存在过度离散（overdispersion，方差 > 均值）时，需使用负二项回归替代泊松回归
- 队列研究中估计发病率比值（IRR / RR）

**不适合的场景：**
- 二分类结局变量（是否患病/死亡）→ 使用 Logistic 回归
- 连续型结局变量（血压值、BMI等）→ 使用线性回归
- 生存时间结局 → 使用 Cox 回归 / 生存分析
- 因变量服从二项分布且可计算比例（proportion）→ 考虑 Logistic 回归

## 前置条件

**R包安装：**

```r
install.packages("MASS")       # glm.nb() 负二项回归
install.packages("haven")      # 读取SPSS格式数据
install.packages("performance") # check_overdispersion() 过度离散检验
```

**数据格式要求：**
- 因变量：非负整数计数（0, 1, 2, …），代表事件发生次数
- 自变量：分类或连续变量均可，分类变量会被自动编码为哑变量
- 若需要有偏移量（如观察人年），需提供log变换后的偏移变量

**统计假设：**
- 泊松回归：计数资料的均值 ≈ 方差（等离散，equidispersion）
- 负二项回归：允许方差 > 均值，以离散参数 theta 衡量额外变异
- 各观测之间相互独立
- log 链接函数下自变量与因变量的对数期望值呈线性关系

## 方法选择决策树

```
你的计数资料 →
├── 均值 ≈ 方差（无过度离散）→ 泊松回归（glm, family=poisson）
├── 方差 > 均值（存在过度离散）→
│   ├── 只需调整标准误，不改变系数 → 拟泊松回归（family=quasipoisson）
│   └── 需要更精确建模 → 负二项回归（glm.nb, MASS包）
├── 不同观测有不同暴露量（观察人年/面积等）→ 加入 offset=log(暴露量)
└── 因变量为0/1二分类 → 改用 Logistic 回归
```

## 标准工作流

### 步骤1：数据准备与探索

读取数据，检查数据结构。对于频数表格式的数据，需先用 `rep()` 转换为原始数据格式（一行一个观测）。确认因变量为计数，检查各变量的因子水平和参照组设置。

```r
# 对于有偏移量的数据（如观察人年），确认偏移变量为正数
summary(data$N)
# 检查因变量分布
table(data$Y)
```

### 步骤2：前提条件检验

在建立泊松回归前，需检验数据是否存在过度离散。建立初步模型后用 `check_overdispersion()` 检验：

```r
library(performance)
f_pois <- glm(Y ~ X1 + X2, data = df, family = poisson())
check_overdispersion(f_pois)
# dispersion ratio > 1 且 p < 0.05 → 存在过度离散，改用负二项回归
```

也可通过残差偏差（residual deviance）和 Pearson 卡方检验评估模型拟合优度：

```r
# 残差偏差法
1 - pchisq(deviance(f_pois), df = f_pois$df.residual)
# Pearson卡方法
Pearson <- sum((df$Y - f_pois$fitted.values)^2 / f_pois$fitted.values)
1 - pchisq(Pearson, df = f_pois$df.residual)
# p < 0.05 提示拟合不佳，可能存在过度离散
```

### 步骤3：执行统计分析

根据步骤2结果选择模型：

```r
# 泊松回归（含偏移量）
f <- glm(Y ~ X1 + YEARGRP, data = data, family = poisson(), offset = log(N))

# 拟泊松回归（处理过度离散，系数不变，SE调整）
f_qp <- glm(Y ~ X1 + YEARGRP, data = data, family = quasipoisson(), offset = log(N))

# 负二项回归
library(MASS)
f_nb <- glm.nb(y ~ x, data = df)
```

### 步骤4：结果解读

查看回归系数和显著性：

```r
summary(f)
coef(f)          # 回归系数
confint(f)       # 系数的95%置信区间（profile likelihood方法）
exp(coef(f))     # RR/IRR值（发病率比值）
exp(confint(f))  # RR/IRR的95%置信区间
```

### 步骤5：结果报告（论文中的统计描述）

论文中报告格式示例：

> 在控制年龄因素后，砷暴露组因呼吸道疾病死亡的风险是非暴露组的 2.25 倍（IRR = 2.25，95%CI：1.77～2.85，p < 0.001），差异具有统计学意义。

或按要求的三线表格式列出各变量的 IRR、95%CI 和 p 值。

## 代码示例

### 示例1：泊松回归（含偏移量）

孙振球《医学统计学》第4版例18-1。研究砷暴露与因呼吸道疾病死亡之间的关系（回顾性队列研究），不同年龄组有不同的人年观察数，需使用偏移量。

```r
library(haven)

# 读取数据
data18_1 <- haven::read_sav("datasets/例18-01.sav", encoding = "GBK")
data18_1 <- haven::as_factor(data18_1)
str(data18_1)
## tibble [8 × 7] (S3: tbl_df/tbl/data.frame)
##  $ X1     : Factor w/ 2 levels "有暴露","无暴露"
##  $ YEARGRP: Factor w/ 4 levels "40-49岁组","50-59岁组",..
##  $ Y      : num 死亡数
##  $ N      : num 观察单位数（人年）

# 修改参照水平：无暴露作参照
levels(data18_1$X1)  ## [1] "有暴露" "无暴露"
data18_1$X1 <- factor(data18_1$X1, levels = c("无暴露", "有暴露"))

# 建立泊松回归模型，N为观察人年，取log后作为偏移量
f <- glm(Y ~ X1 + YEARGRP, data = data18_1, family = poisson(), offset = log(N))
summary(f)
## 
## Call:
## glm(formula = Y ~ X1 + YEARGRP, family = poisson(), data = data18_1, 
##     offset = log(N))
## 
## Coefficients:
##                  Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       -8.0086     0.2233 -35.859  < 2e-16 ***
## X1有暴露           0.8109     0.1210   6.699 2.09e-11 ***
## YEARGRP50-59岁组   1.4702     0.2453   5.994 2.04e-09 ***
## YEARGRP60-69岁组   2.3661     0.2372   9.976  < 2e-16 ***
## YEARGRP>=70岁组    2.6238     0.2548  10.297  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for poisson family taken to be 1)
##     Null deviance: 260.9304  on 7  degrees of freedom
## Residual deviance:   9.9303  on 3  degrees of freedom
## AIC: 61.342

# 提取结果
coef(f)
##      (Intercept)         X1有暴露 YEARGRP50-59岁组 YEARGRP60-69岁组 
##       -8.0086495        0.8108698        1.4701505        2.3661111 
##  YEARGRP>=70岁组 
##        2.6237532

confint(f)
##                       2.5 %    97.5 %
## (Intercept)      -8.4780952 -7.598323
## X1有暴露          0.5724818  1.047494
## YEARGRP50-59岁组  1.0094838  1.976026
## YEARGRP60-69岁组  1.9239329  2.858522
## YEARGRP>=70岁组   2.1411250  3.145495

exp(coef(f))     # IRR/RR值
##      (Intercept)         X1有暴露 YEARGRP50-59岁组 YEARGRP60-69岁组 
##     3.325736e-04     2.249864e+00     4.349890e+00     1.065587e+01 
##  YEARGRP>=70岁组 
##     1.378737e+01

exp(confint(f))  # IRR/RR的95%置信区间
##                         2.5 %       97.5 %
## (Intercept)      0.0002079745 5.012913e-04
## X1有暴露         1.7726609281 2.850500e+00
## YEARGRP50-59岁组 2.7441841472 7.214017e+00
## YEARGRP60-69岁组 6.8478373075 1.743573e+01
## YEARGRP>=70岁组  8.5090045641 2.323118e+01

# 拟合优度检验：残差偏差法
1 - pchisq(deviance(f), df = f$df.residual)
## [1] 0.01916785  # p < 0.05，拟合欠佳

# 拟合优度检验：Pearson卡方法
Pearson <- sum((data18_1$Y - f$fitted.values)^2 / f$fitted.values)
1 - pchisq(Pearson, df = f$df.residual)
## [1] 0.02137005  # p < 0.05，拟合欠佳

# 过度离散检验
library(performance)
check_overdispersion(f)
## # Overdispersion test
##        dispersion ratio = 3.231
##   Pearson's Chi-Squared = 9.692
##                 p-value = 0.021
## → dispersion ratio = 3.231，数据存在过度离散

# 改用拟泊松回归
f1 <- glm(Y ~ X1 + YEARGRP, data = data18_1, family = quasipoisson(), offset = log(N))
summary(f1)
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -8.0086     0.4014 -19.950 0.000275 ***
## X1有暴露           0.8109     0.2176   3.727 0.033640 *  
## YEARGRP50-59岁组   1.4702     0.4408   3.335 0.044555 *  
## YEARGRP60-69岁组   2.3661     0.4263   5.550 0.011534 *  
## YEARGRP>=70岁组    2.6238     0.4580   5.729 0.010558 *  
## ---
## (Dispersion parameter for quasipoisson family taken to be 3.23081)
## → 系数不变，标准误变大，p值变得更保守
```

### 示例2：负二项回归

孙振球《医学统计学》第4版例18-2。研究居住地类型与蚊虫幼虫滋生容器数的关系。频数表数据需先展开为原始格式。

```r
library(haven)
library(MASS)

# 读取数据（频数表格式）
data18_2 <- haven::read_sav("datasets/例18-02.sav", encoding = "GBK")
data18_2 <- haven::as_factor(data18_2)
str(data18_2)
## tibble [24 × 3] (S3: tbl_df/tbl/data.frame)
##  $ x    : Factor w/ 3 levels "农村","城市贫民区","城市"
##  $ y    : num  受蚊子幼虫滋生的容器数
##  $ count: num  频数

head(data18_2)
## # A tibble: 6 × 3
##   x         y count
##   <fct> <dbl> <dbl>
## 1 农村      0   136
## 2 农村      1    23
## 3 农村      2    10
## 4 农村      3     5
## 5 农村      4     2
## 6 农村      5     1

# 将频数表转换为原始数据格式
x1 <- rep(data18_2$x, data18_2$count)
y1 <- rep(data18_2$y, data18_2$count)
df <- data.frame(x1, y1)
str(df)  ## 'data.frame': 299 obs. of 2 variables

xtabs(~ x1 + y1, data = df)  # 确认与原始数据一致
##             y1
## x1             0   1   2   3   4   5   6  11
##   农村       136  23  10   5   2   1   1   1
##   城市贫民区  38   8   2   0   0   0   0   0
##   城市        67   5   0   0   0   0   0   0

# 负二项回归
f <- glm.nb(y1 ~ x1, data = df)
summary(f)
## 
## Call:
## glm.nb(formula = y1 ~ x1, data = df, init.theta = 0.3002652205, 
##     link = log)
## 
## Coefficients:
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -0.7100     0.1731  -4.102  4.1e-05 ***
## x1城市贫民区  -0.6762     0.4274  -1.582 0.113612    
## x1城市        -1.9572     0.5256  -3.724 0.000196 ***
## ---
## (Dispersion parameter for Negative Binomial(0.3003) family taken to be 1)
##     Null deviance: 174.95  on 298  degrees of freedom
## Residual deviance: 156.37  on 296  degrees of freedom
## AIC: 426.23
##               Theta:  0.3003 
##           Std. Err.:  0.0764 
##  2 x log-likelihood:  -418.2280

# 计算IRR的95%CI（以城市 vs 农村为例）
# 系数均值 ± 1.96 × 标准误 → 再取指数
cbind(exp(-1.9572 + 1.96 * 0.5256), exp(-1.9572 - 1.96 * 0.5256))
##          [,1]     [,2]
## [1,] 2.526917 19.83405
## → 农村家庭滋生蚊虫幼虫机会是城市家庭的 exp(1.9572) = 7.08 倍（95%CI：2.53～19.83）
```

## 结果解读指南

**模型输出核心要素：**

| 输出项 | 含义 | 解读示例 |
|--------|------|----------|
| `Estimate` | 回归系数 β | 正值→风险增加，负值→风险降低 |
| `exp(Estimate)` | IRR / RR（发病率比值） | IRR = 2.25 表示暴露组发生率是对照组的2.25倍 |
| `Std. Error` | 系数的标准误 | 用于计算置信区间 |
| `z value` / `t value` | 检验统计量 | 绝对值 > 1.96 通常 p < 0.05 |
| `Pr(>|z|)` | p 值 | p < 0.05 表示该变量有统计学意义 |
| `Null deviance` | 空模型偏差 | 自由度为 n-1 |
| `Residual deviance` | 残差偏差 | 越小越好，自由度越小，说明使用参数越多 |
| `AIC` | 赤池信息准则 | 用于模型比较，值越小拟合越好 |
| `Dispersion parameter` | 离散参数 | poisson默认为1；quasipoisson >1 表示过度离散 |
| `Theta` | 负二项离散参数 | 值越小过度离散越严重；接近1表示接近泊松分布 |

**结果解释要点：**
- 基准发病密度：exp(截距) = 所有自变量取参照水平时的事件发生率
- 各哑变量系数取 exp 后为相对于参照组的 IRR
- 若包含偏移量，结果自动调整为发病率密度（每单位暴露量的事件数）
- 论文表述模板："在控制XX因素后，XX组发生YY的风险是对照组的IRR倍（95%CI：XX～XX，p = XX）"

**poisson vs quasipoisson vs neg-binomial：**
- `poisson`：估计系数的点值和标准误，假设等离散；若数据存在过度离散，标准误被低估，p 值偏小（假阳性增加）
- `quasipoisson`：系数不变，标准误乘以 sqrt(dispersion)，p 值更保守
- `glm.nb`：系数和标准误均重新估计，模型更灵活，适合明显过度离散的数据

## 常见问题与注意事项

**Q1：什么时候需要设置偏移量（offset）？**
当不同观测的暴露量不同时。如不同年龄组的观察人年数不同，或不同区域面积不同。公式：`offset = log(暴露量)`。偏移量在模型中是系数固定为1的项，相当于在等式右侧减去 log(暴露量)，从而将计数转化为率。

**Q2：如何判断是否存在过度离散？**
三种方法：① 计算方差/均值比，>1提示过度离散；② 用 `check_overdispersion()` 检验，dispersion ratio > 1 且 p < 0.05；③ 残差偏差或 Pearson 卡方检验 p < 0.05 提示偏离泊松假设。

**Q3：拟泊松回归（quasipoisson）和负二项回归有什么区别？**
拟泊松回归只调整标准误，不改变系数估计值，相当于在标准误上乘以 sqrt(dispersion)。负二项回归重新估计所有参数，通过 theta 参数显式建模数据变异，灵活度更高。数据过度离散不严重时可先试 quasipoisson，离散严重时用负二项回归。

**Q4：SPSS与R在泊松回归输出上的差异？**
系数估计值（Estimate）通常一致，但标准误和置信区间的计算方法可能不同，导致略有差异。SPSS 默认使用 Wald 方法，R 的 `confint()` 默认使用 profile likelihood 方法（比 Wald 更准确）。可通过 `confint.default(f)` 获得与 SPSS 一致的 Wald 置信区间。

**Q5：分类变量的因子水平（参照组）如何设定？**
R 默认将因子中字母顺序最小的水平作为参照组。可通过 `factor(x, levels = c("对照组", "实验组"))` 将对照组设为第一水平（参照组）。`relevel()` 也可用于快速更改参照水平。

**Q6：计数数据中含有大量零怎么办？**
如果零值过多（zero-inflation），标准泊松和负二项可能拟合不佳。此时需考虑零膨胀泊松模型（zero-inflated Poisson, ZIP）或零膨胀负二项模型（ZINB），可使用 `pscl` 包中的 `zeroinfl()` 函数。

**注意事项：**
- 偏移量务必取自然对数 `log()`，否则模型系数解释不正确
- 过度离散检验应在泊松回归后进行，不要凭直觉跳过
- 拟泊松回归不提供 AIC（显示 NA），因为似然函数不是完全定义的
- 负二项回归的 theta 值越小说明过度离散越严重，theta 趋近无穷时等价于泊松回归
- 样本量很小时（观测数 < 20），过度离散检验的把握度不足，建议结合学科经验判断
- 因变量中的计数必须是非负整数，不应包含小数或负数
