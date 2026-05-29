---
name: medical-stat-multilevel
description: "R语言医学统计：多水平模型（混合效应模型）。用于处理层次结构数据（如患者嵌套于医院、重复测量嵌套于个体），涵盖随机截距模型、随机系数模型、交叉层交互效应和重复测量数据建模，使用lme4和nlme包。TRIGGER when user mentions 多水平模型、混合效应模型、随机效应、lme4、lmer、层次结构、多中心数据、重复测量，or asks about analyzing clustered/hierarchical data, multilevel models, mixed models, random intercept, random slope. SKIP for 广义估计方程(GEE)、固定效应单水平模型、普通线性回归（若数据无分层结构）。"
---
# 多水平模型 (Multilevel Models / Mixed Effects Models)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

### 适用情况
- **层次结构数据**：个体嵌套于群体（如患者嵌套于医院、学生嵌套于学校）
- **重复测量数据**：同一患者多个时间点的测量数据（时间点为水平1，患者为水平2）
- **多中心临床试验**：不同研究中心（水平2）的患者（水平1）
- **纵向追踪数据**：对同一批受试者进行多次随访
- **因变量为连续型**：使用 `lmer()` 拟合线性混合效应模型
- **因变量为离散型**：使用 `glmer()` 拟合广义线性混合效应模型（逻辑、泊松等）

### 不适用情况
- **数据无层次结构**：普通数据用 `lm()` 或 `glm()` 即可
- **仅关注群体平均效应**且数据独立性成立：使用广义估计方程 (GEE, `geepack`) 更合适
- **仅有固定效应**、无需分解方差成分：使用普通方差分析或线性回归
- **需要处理时间序列自相关结构**而非随机效应：考虑 `nlme::gls()` 或 `lme()` 中的相关结构设定

### 医学研究常见应用
- 多中心药物临床试验，分析中心间差异
- 学生成绩与学校特征的关系（教育数据，分析范式通用）
- 患者住院费用与医院特征的关联
- 纵向随访中血糖水平随时间的变化趋势

## 前置条件

### 必需 R 包

```r
install.packages(c("lme4", "lmerTest", "performance", "ggplot2", "dplyr"))
# lme4 用于拟合模型，lmerTest 补充提供 P 值，performance 计算 ICC 和模型比较
```

### 数据格式要求
- **长格式数据**：每行一个观测，包含分组变量（如 `schcode`）、因变量（如 `math`）、自变量（如 `ses`）
- **分组变量**须为因子或整数类型，用于指定随机效应的分组
- **水平2变量**：同一分组内所有个体共享的变量（如学校类型 `public`），在组内不变

### 统计假定
- **随机效应**服从正态分布，均值为0
- **残差**服从正态分布（线性混合模型），满足方差齐性
- **随机效应与残差**相互独立
- **ICC > 0**：若组内相关系数为0，说明无层次结构，不需要多水平模型

## 方法选择决策树

```
你的数据情况 →
├── 因变量连续，仅需考虑不同组截距不同
│   └── 使用随机截距模型（方差成分模型）
│       lmer(y ~ x1 + x2 + (1 | group), data = df)
│
├── 因变量连续，不同组间斜率和截距均不同
│   ├── 先画分组散点图 + 拟合线验证斜率变异
│   └── 使用随机系数模型
│       lmer(y ~ x + (x | group), data = df)
│
├── 含水平2变量，需跨层交互
│   └── 在随机系数模型中加入交互项
│       lmer(y ~ x * level2_var + (x | group), data = df)
│
├── 重复测量数据（每个患者多个时间点）
│   ├── 时间点为水平1，患者 ID 为水平2
│   ├── 仅需随机截距：lmer(y ~ time + (1 | id), data = df)
│   └── 需随机斜率：lmer(y ~ time + (time | id), data = df)
│
├── 因变量为二分类（如发病/未发病）
│   └── 广义线性混合模型
│       glmer(y ~ x + (1 | group), family = binomial, data = df)
│
├── 因变量为计数（如发病次数）
│   └── glmer(y ~ x + (1 | group), family = poisson, data = df)
│
└── 不确定是否需要多水平模型
    └── 先拟合空模型，计算 ICC
        null_model <- lmer(y ~ (1 | group), data = df)
        performance::icc(null_model)
        如果 ICC ≈ 0 → 用普通回归即可
```

## 标准工作流

### 步骤1：数据准备与探索

1. 加载数据，检查结构 `str(data)`
2. 确认分组变量和层次关系（学生→学校，患者→医院等）
3. 对连续变量作描述统计 `summary(data)`
4. 按分组绘制散点图 + 拟合线，直观观察截距/斜率差异

### 步骤2：前提条件检验

1. **拟合空模型**（仅含随机截距，无自变量）
2. **计算 ICC**：`performance::icc(null_model)`，判断数据是否存在层次结构
3. **残差诊断**：`plot(model)` 检查残差正态性和方差齐性
4. 若 ICC < 0.05 且不显著，可认为层次结构不明显，但若有先验理论支持，仍可使用

### 步骤3：执行统计分析

1. 从空模型开始，逐步加入固定效应（先水平1，后水平2）
2. 通过 `summary()` 查看固定效应和随机效应估计值
3. 使用似然比检验或 AIC/BIC 比较嵌套模型：
   ```r
   anova(model1, model2)  # 似然比检验（ML估计）
   ```
4. 确认是否需随机斜率：对比随机截距模型和随机系数模型的 AIC/BIC
5. 若需要，加入跨层交互项

### 步骤4：结果解读

- 固定效应解释方法同普通回归
- 关注随机效应的方差大小——方差越大说明群体差异越大
- 计算边际 R²（固定效应解释的变异）和条件 R²（固定+随机效应解释的变异）

### 步骤5：结果报告

论文中报告格式示例：
- "采用两水平线性混合效应模型分析数据，以学校为随机截距，社会经济地位为固定效应。结果显示，社会经济地位与数学成绩呈显著正相关（β = 3.96, SE = 0.14, p < 0.001），学校间差异的 ICC 为 0.138。模型的条件 R² 为 0.167。"

## 代码示例

### 示例1：空模型与ICC计算

```r
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(ggplot2)

data <- read.csv('datasets/heck2011.csv')
str(data)
## 'data.frame':   6871 obs. of  10 variables:
##  $ schcode : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ math    : num  47.1 63.6 57.7 53.9 58 ...
##  $ ses     : num  0.586 0.304 -0.544 -0.848 0.001 ...
##  $ female  : int  1 1 1 0 0 0 0 1 0 1 ...
##  $ public  : int  0 0 0 0 0 0 0 0 0 0 ...

# 空模型
null_model <- lmer(math ~ (1 | schcode), data = data)
summary(null_model)
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  schcode  (Intercept) 10.64    3.262
##  Residual             66.55    8.158
## Number of obs: 6871, groups:  schcode, 419
##
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)  57.6742     0.1883 416.0655   306.3   <2e-16 ***

# 计算ICC
performance::icc(null_model)
##     Adjusted ICC: 0.138
##   Unadjusted ICC: 0.138
```

### 示例2：添加水平1固定效应（随机截距模型）

```r
ses_l1 <- lmer(math ~ ses + (1 | schcode), data = data, REML = TRUE)
summary(ses_l1)
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  schcode  (Intercept)  3.469   1.863
##  Residual             62.807   7.925
##
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)   57.5960     0.1329  375.6989  433.36   <2e-16 ***
## ses            3.8739     0.1366 3914.6382   28.35   <2e-16 ***

# 与普通线性回归比较
f <- lm(math ~ ses, data = data)
compare_performance(f, ses_l1, metrics = "common")
## Name   |           Model |   AIC (weights) |   BIC (weights) |  RMSE |    R2
## ----------------------------------------------------------------------------
## f      |              lm | 48304.0 (<.001) | 48324.5 (<.001) | 8.131 | 0.143
## ses_l1 | lmerModLmerTest | 48219.1 (>.999) | 48246.4 (>.999) | 7.810 |

# 可信区间
confint(ses_l1)
##                 2.5 %    97.5 %
## (Intercept) 57.335234 57.856673
## ses          3.596455  4.152745
```

### 示例3：添加水平2固定效应

```r
ses_l1_public_l2 <- lmer(math ~ ses + public + (1 | schcode),
                         data = data, REML = TRUE)
summary(ses_l1_public_l2)
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)
## (Intercept)   57.63143    0.25535  381.81733 225.693   <2e-16 ***
## ses            3.87338    0.13673 3928.37427  28.329   <2e-16 ***
## public        -0.04859    0.29862  385.93649  -0.163    0.871

# 多模型比较
compare_performance(null_model, ses_l1, ses_l1_public_l2)
## Name             |   AIC (weights) | R2 (cond.) | R2 (marg.) |   ICC
## ---------------------------------------------------------------------
## null_model       | 48881.8 (<.001) |      0.138 |      0.000 | 0.138
## ses_l1           | 48219.1 (0.729) |      0.167 |      0.121 | 0.052
## ses_l1_public_l2 | 48221.1 (0.271) |      0.167 |      0.121 | 0.053
```

### 示例4：随机系数模型（含随机斜率）

```r
# 先按学校分组画图（取前10个学校）
data_sub <- data %>% filter(schcode <= 10)
ggplot(data_sub, aes(x = ses, y = math, colour = factor(schcode))) +
  geom_point() +
  geom_smooth(aes(group = schcode), method = "lm", se = FALSE,
              fullrange = TRUE) +
  labs(colour = "schcode")

# 拟合随机系数模型
ses_l1_random <- lmer(math ~ ses + (ses | schcode),
                      data = data, REML = TRUE)
summary(ses_l1_random)
## Random effects:
##  Groups   Name        Variance Std.Dev. Corr
##  schcode  (Intercept)  3.2042  1.7900
##           ses          0.7794  0.8828   -1.00
##  Residual             62.5855  7.9111
##
## Fixed effects:
##              Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)   57.6959     0.1315  378.6378  438.78   <2e-16 ***
## ses            3.9602     0.1408 1450.7730   28.12   <2e-16 ***

# 查看每个学校的随机效应
head(ranef(ses_l1_random)$schcode, 5)
##   (Intercept)        ses
## 1   0.9746643 -0.4806908
## 2   1.0450461 -0.5154022
## 3  -3.4842479  1.7183825
## 4   1.8810911 -0.9277279
## 5  -3.8147866  1.8813996
```

### 示例5：跨层交互效应

```r
crosslevel_model <- lmer(math ~ ses * public + (ses | schcode),
                         data = data, REML = TRUE)
summary(crosslevel_model)
## Fixed effects:
##               Estimate Std. Error         df t value Pr(>|t|)
## (Intercept)   57.72440    0.25183  382.39815 229.216   <2e-16 ***
## ses            4.42383    0.27427 1283.55623  16.130   <2e-16 ***
## public        -0.02632    0.29472  387.41741  -0.089   0.9289
## ses:public    -0.62520    0.31957 1363.95274  -1.956   0.0506 .
```

### 示例6：重复测量数据的两水平模型

```r
# 模拟降压药治疗数据
data12_1 <- data.frame(
  id = factor(c(1:10, 1:10)),
  stat = factor(rep(c("治疗前", "治疗后"), each = 10),
                levels = c("治疗前", "治疗后")),
  bp = c(130, 124, 136, 128, 122, 118, 116, 138, 126, 124,
         114, 110, 126, 116, 102, 100, 98, 122, 108, 106)
)

# 可视化
ggplot(data12_1, aes(stat, bp)) +
  geom_line(aes(color = id, group = id))

# 随机截距模型
f_rpt <- lmer(bp ~ stat + (1 | id), data = data12_1)
summary(f_rpt)
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)
## (Intercept) 126.2000     2.6153   9.6662   48.25 7.56e-13 ***
## stat治疗后  -16.0000     0.9888   9.0000  -16.18 5.83e-08 ***
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  id       (Intercept) 63.511   7.969
##  Residual              4.889   2.211
```

## 结果解读指南

### lmer 输出结构

输出分为两大部分：

**1. Random effects（随机效应）**
| 字段 | 含义 |
|------|------|
| Groups | 分组变量名称 |
| Name | 随机效应类型（Intercept 截距 / 自变量名 斜率） |
| Variance | 该随机效应的方差——越大表示群体间差异越大 |
| Std.Dev. | 标准差（方差的平方根） |
| Corr | 随机截距与随机斜率的相关系数 |
| Residual | 残差方差——模型无法解释的变异 |

**2. Fixed effects（固定效应）**

| 字段 | 含义 |
|------|------|
| Estimate | 回归系数（β）：自变量每变化1单位，因变量变化的量 |
| Std. Error | 标准误 |
| df | 自由度（Satterthwaite 近似） |
| t value | t 统计量 |
| Pr(>|t|) | P 值（lmerTest 提供） |

### ICC 解读

- ICC = 组间方差 / (组间方差 + 残差方差)
- ICC = 0：群体间无差异，无需多水平模型
- ICC = 0.138（本例）：13.8% 的成绩变异由学校差异造成，其余由学生个体差异造成
- 加入自变量后 ICC 下降，说明部分变异被自变量解释了

### R² 解读
- **Marginal R²**：仅固定效应解释的变异比例
- **Conditional R²**：固定效应 + 随机效应共同解释的变异比例

### 论文报告句式

- "多水平模型显示，[自变量] 与 [因变量] 显著相关（β = X.XX，95%CI: X.XX–X.XX，p = 0.XXX）。"
- "空模型的 ICC 为 0.XXX，表明 [比例]% 的总变异可归因于 [群体水平] 的差异。"
- "加入 [自变量] 后，组间方差由 X.XX 降至 X.XX，条件 R² 为 0.XXX。"

## 常见问题与注意事项

### lme4 vs lmerTest：P 值从哪里来？
- `lme4::lmer()` 默认**不输出 P 值**，因为混合模型中自由度的计算没有一致方法
- `lmerTest::lmer()` 使用 Satterthwaite 近似计算自由度，从而提供 P 值
- 加载 `library(lmerTest)` 后，`lmer()` 函数会被自动增强，`summary()` 会显示 P 值
- 加载顺序应为 `library(lmerTest)` 在 `library(lme4)` 之后

### REML vs ML：用哪个估计方法？
- **REML（限制性最大似然）**：默认方法，对随机效应方差估计更准确，推荐用于最终模型
- **ML（最大似然）**：用于模型比较（`anova()` 似然比检验），因为 REML 下比较固定效应不同的模型不可靠
- 规则：比较模型时用 `REML = FALSE`（ML），最终报告用 `REML = TRUE`

### "boundary (singular) fit" 警告
- 表示某些随机效应的方差接近0或相关系数为±1
- 常见原因：
  1. 随机效应太复杂而数据不足以估计（减少随机效应项）
  2. 分组数量太少（需要 ≥10-15 个组）
  3. 组内样本量太小
- 处理：简化模型（仅保留随机截距），或使用 `nlme::lme()` 调整优化器

### lmer 公式语法速查

| 公式写法 | 含义 |
|----------|------|
| `(1 | group)` | 仅随机截距 |
| `(x | group)` | 随机截距 + 随机斜率（等价 `(1 + x | group)`） |
| `(0 + x | group)` | 仅随机斜率，无随机截距 |
| `(1 | group1) + (1 | group2)` | 交叉随机效应 |
| `x * z` | 等价 `x + z + x:z`，包含交互项 |

### 与 SPSS 的结果差异
- SPSS MIXED 过程和 `lmer` 的默认估计方法可能有别（SPSS 默认 REML 时对自由度处理不同）
- P 值计算方式不同：SPSS 可能使用不同的 df 近似
- 系数估计通常一致，但 P 值和置信区间边界可能有微小差异

### 何时不需要多水平模型？
- 数据没有自然的分组/层次结构
- ICC 非常接近 0 且模型比较显示 AIC/BIC 没有改善
- 仅有少量分组（< 5 个组），随机效应方差难以准确估计——此时可将分组作为固定效应

### 广义混合效应模型（glmer）
- 当因变量为二分类时，用 `glmer(..., family = binomial)`
- 当因变量为计数时，用 `glmer(..., family = poisson)`
- 当数据存在过度离散时，考虑 `family = binomial` 后叠加观察水平随机效应 `(1 | obs)`，或使用 `glmmTMB` 包的负二项族

### 参考资料
- 冯国双：多水平模型系列推文（微信公众号，推荐入门阅读）
- Introduction to Multilevel Modelling: https://www.learn-mlms.com/
- CRAN Task View: Mixed Models: https://cran.r-project.org/web/views/MixedModels.html
- A Cheatsheet for Building Multilevel Models in R: https://paulrjohnson.net/blog/2022-11-01-multilevel-model-r-cheatsheet/
