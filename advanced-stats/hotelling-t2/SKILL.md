---
name: medical-stat-hotelling-t2
description: "R语言医学统计：多变量统计分析。涵盖多变量数据的描述性统计、Hotelling T²检验（多变量t检验），用于同时比较多个连续型结局变量在两组或多组间的差异。TRIGGER when user mentions 多变量分析、Hotelling T²、多个结局同时比较、multivariate t-test、多元正态检验、轮廓分析，or asks about comparing multiple outcomes simultaneously。SKIP for MANOVA多因素方差分析、多变量回归、主成分分析、聚类分析、判别分析。"
---

# 多变量数据的统计描述和统计推断 (Hotelling's T² / Multivariate Analysis)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用：**

- 每个观察对象记录多个连续型反应变量（结局指标），需同时比较这些变量在组间的差异
- 临床试验中同时记录收缩压、舒张压作为联合结局
- 多次重复测量数据，不考虑球对称假设时直接采用多变量方法
- 比较两组或多组轮廓（多个指标构成的曲线形状）是否平行、重合

**何时不应使用：**

- 只有一个结局变量 → 使用单变量t检验或方差分析
- 有多个分组因素但只有一个结局变量（多因素试验） → 使用多因素方差分析，这是单变量数据而非多变量数据
- 变量间存在自变量-因变量关系 → 使用多元回归
- 探索性降维 → 使用主成分分析、因子分析
- 样本间分类 → 使用判别分析、聚类分析

**医学研究常见应用场景：**

- 血脂指标（甘油三酯、总胆固醇、高密度脂蛋白）在患者组与对照组的差异
- 新生儿出生体重与身长的联合比较
- 减肥药物在服药后多个时间点体重的变化
- 健康问卷调查多个问题的回答在两个人群中的差异

## 前置条件

**Required R packages:**

```r
install.packages(c("haven", "ICSNP", "Hotelling", "MVN", "mvnormtest",
                   "rrcov", "ggplot2", "tidyr", "dplyr", "Hmisc", "profileR"))
# profileR可能需要从GitHub安装：
# devtools::install_github("cddesja/profileR", build_vignettes = TRUE)
```

**Required data format:**

- 宽格式（每行一个观察对象，每列一个变量）
- 所有反应变量必须为连续型数值（numeric/double）
- 分组变量需为factor或numeric
- 数据通常为`.sav`格式，使用`haven::read_sav()`读取

**Statistical assumptions:**

- 多元正态分布：多个反应变量的联合分布服从多元正态
- 组间协方差矩阵齐同（两组或多组时要求）
- 观察对象互相独立

## 方法选择决策树

```
你的数据情况 →
├── 单组资料（与已知总体均值向量比较）→ 单样本Hotelling T²检验
│     └── ICSNP::HotellingsT2(X, mu=c(…), test="f") 或 rrcov::T2.test()
├── 两组资料（比较两组的多个均值向量）→ 两样本Hotelling T²检验
│     ├── 要SPSS一致结果 → Hotelling::hotelling.test()
│     └── 只要P值正确 → ICSNP::HotellingsT2() 或 rrcov::T2.test()
├── 多组资料（≥3组，比较多个均值向量）→ MANOVA
│     └── manova(cbind(y1,y2,…) ~ group), summary(…, test="Wilks")
├── 重复测量数据（同一人多次测量）→ 单样本Hotelling T²检验（对改变值）
│     └── 不依赖球对称假设，作为重复测量ANOVA的替代
└── 比较两组轮廓 → 轮廓分析
      └── profileR::pbg(data, group), 按顺序：平行→相合→水平检验
```

## 标准工作流

### 步骤1：数据准备与探索

```r
library(haven)
data <- haven::read_sav("datasets/filename.sav")
str(data); head(data)        # 数据结构预览
colMeans(data)               # 均值向量
cov(data)                    # 协方差矩阵
cor(data)                    # 相关矩阵
library(Hmisc)
rcorr(as.matrix(data))$P     # 相关系数的P值矩阵
```

### 步骤2：前提条件检验

```r
# 多元正态性检验 — 两种方法，P>0.05满足
library(mvnormtest)
mshapiro.test(t(data))       # 注意转置

library(MVN)
mvn(data)$multivariateNormality  # 默认Henze-Zirkler法
```

### 步骤3：执行统计分析

根据决策树选择对应方法：

- 单样本：`ICSNP::HotellingsT2(X, mu=pop_mean, test="f")`
- 两样本：`Hotelling::hotelling.test(Y ~ group)`
- 多组：`manova(cbind(y1,y2,…) ~ group)` → `summary(…, test="Wilks")`
- 轮廓分析：`profileR::pbg(data=response_cols, group=group_var)`

### 步骤4：结果解读

- T²值/F值大小：反映组间多变量均数向量差异程度
- P值：`p < 0.05`表示组间多变量差异有统计学意义
- 结合各变量均值判断方向：如“治疗组甘油三酯和总胆固醇高于正常组，高密度脂蛋白低于正常组”

### 步骤5：结果报告

论文报告示例：
> 采用Hotelling T²检验同时比较两组新生儿的体重和身长，结果显示两组在体重和身长构成的二维平面分布上的差异有统计学意义（T²=9.87, F=4.58, P=0.031）。

> 多变量方差分析显示三组慢性胃炎儿童T细胞免疫功能（T3、T4、T8细胞百分比）间差异有统计学意义（Wilks Lambda=0.089, F=5.50, P=0.004）。

## 代码示例

### 示例1：多变量统计描述（课本例14-1）

15名正常成年男性血脂：x1=甘油三酯, x2=总胆固醇, x3=高密度脂蛋白。

```r
library(haven)
library(Hmisc)

data14_1 <- haven::read_sav("datasets/例14-01.sav")
str(data14_1)
## tibble [15 × 3] (S3: tbl_df/tbl/data.frame)
##  $ x1: num 1.06 0.98 0.85 0.96 0.98 0.99 1.01 1.02 1.02 1.1 ...
##  $ x2: num 2.56 2.42 2.35 2.55 2.65 2.6 2.35 2.89 2.54 2.64 ...
##  $ x3: num 1.93 1.8 1.68 1.34 2.55 2.33 1.93 1.8 1.68 1.34 ...

colMeans(data14_1)
##       x1       x2       x3
## 1.020000 2.728667 2.043333

cov(data14_1)
##             x1         x2          x3
## x1 0.005757143 0.01029286 0.009314286
## x2 0.010292857 0.08864095 0.080211905
## x3 0.009314286 0.08021190 0.186838095

cor(data14_1)
##           x1        x2        x3
## x1 1.0000000 0.4556331 0.2839967
## x2 0.4556331 1.0000000 0.6232882
## x3 0.2839967 0.6232882 1.0000000

rcorr(as.matrix(data14_1))$P
##            x1         x2        x3
## x1         NA 0.08785736 0.3049796
## x2 0.08785736         NA 0.0130474
## x3 0.30497964 0.01304740        NA
```

均值向量描述三个指标的平均水平，协方差矩阵描述变异程度，相关矩阵描述相关性。甘油三酯与总胆固醇r=0.456(P=0.088)，总胆固醇与高密度脂蛋白r=0.623(P=0.013)。

### 示例2：多元正态性检验

```r
library(mvnormtest)
library(MVN)

mshapiro.test(t(data14_1))
## W = 0.90331, p-value = 0.1069

mvn(data14_1)$multivariateNormality
##            Test        HZ   p value MVN
## 1 Henze-Zirkler 0.5079275 0.5592683 YES
## P > 0.05，满足多元正态分布
```

### 示例3：单样本Hotelling T²检验（课本例14-2）

5名怀疑冠心病男性 vs 已知正常人群总体均值（甘油三酯=1.02, 总胆固醇=2.73, 高密度脂蛋白=2.04）。

```r
library(ICSNP)
data14_2 <- haven::read_sav("datasets/例14-02.sav")

ICSNP::HotellingsT2(X = data14_2, mu = c(1.02, 2.73, 2.04), test = "f")
## T.2 = 2389.8, df1 = 3, df2 = 2, p-value = 0.0004183
## P < 0.01，冠心病男性血脂与正常男性有差异

rrcov::T2.test(x = data14_2, mu = c(1.02, 2.73, 2.04))
## T2 = 14338.9, F = 2389.8, df1 = 3, df2 = 2, p-value = 0.0004183
```

注意：`ICSNP`结果中`T.2`实际为F值（与`rrcov`的`F`一致），而非真正的T²值。

### 示例4：两样本Hotelling T²检验（课本例14-3）

孕期保健教育组（group=1, n=6）与对照组（group=2, n=7）婴儿体重和身长比较。

```r
library(haven)
library(Hotelling)
library(ICSNP)

data14_3 <- haven::read_sav("datasets/例14-03.sav")

tapply(data14_3[,-1], data14_3$group, FUN = colMeans)
## $`1`
##   weight   height
##  3.64750 51.66667
## $`2`
##    weight    height
##  3.148571 48.571429

# Hotelling包 — 结果与SPSS一致
with(data14_3, Hotelling::hotelling.test(cbind(weight, height) ~ group))
## Test stat:  9.4862, Numerator df: 2, Denominator df: 10, P-value: 0.04463

# ICSNP包 — F值不同但P值相同
with(data14_3, ICSNP::HotellingsT2(cbind(weight, height) ~ group))
## T.2 = 4.3119, df1 = 2, df2 = 10, p-value = 0.04463
## P < 0.05，孕期保健教育对婴儿生长发育有促进作用
```

### 示例5：多组MANOVA（课本例14-4/14-5）

三组慢性胃炎儿童外周血T3、T4、T8细胞百分比（group=1对照, 2治疗I组, 3治疗II组）。

```r
data14_4 <- haven::read_sav("datasets/例14-04.sav")

fit <- with(data14_4, manova(cbind(t3, t4, t8) ~ as.factor(group)))

summary(fit, test = "Wilks")
##                  Df    Wilks approx F num Df den Df   Pr(>F)
## as.factor(group)  2 0.088735   5.4997      6     14 0.004104 **
## Wilks Lambda=0.089, P<0.01，三组T细胞免疫功能有差异

summary.aov(fit)  # 单变量ANOVA作为补充
## Response t3: F=32.55, P=7.58e-05 ***
## Response t4: F=8.72,  P=0.0078 **
## Response t8: F=14.81, P=0.0014 **

# 可选统计量: "Pillai", "Wilks", "Hotelling-Lawley", "Roy"
summary(fit, test = "Pillai")
## Pillai=1.0492, approx F=2.9424, P=0.0394
```

### 示例6：多变量vs单变量对比（课本例14-6）

两组新生儿体重、身长数据——关键演示：单变量不显著而多变量显著。

```r
data14_6 <- haven::read_sav("datasets/例14-06.sav")

t.test(weight ~ group, data = data14_6)
## t = -1.6219, df = 13.822, p-value = 0.1274  ← 体重不显著

t.test(height ~ group, data = data14_6)
## t = -0.036693, df = 11.938, p-value = 0.9713  ← 身长不显著

with(data14_6, Hotelling::hotelling.test(cbind(weight, height) ~ group))
## Test stat: 9.8669, Numerator df: 2, Denominator df: 13, P-value: 0.0312
## ← 多变量显著！单变量与多变量回答不同统计问题，不能相互替代
```

### 示例7：重复测量数据的多变量分析（课本例14-7）

10名肥胖患者服药前（t1）及服药后1-4周（t2-t5）体重。

```r
data14_7 <- haven::read_sav("datasets/例14-07.sav")
data14_7 <- as.data.frame(data14_7)

weight_change <- as.matrix(data14_7[, 3:6] - data14_7[, 2])  # 相对于服药前的改变

ICSNP::HotellingsT2(X = weight_change, mu = c(0, 0, 0, 0))
## T.2 = 41.308, df1 = 4, df2 = 6, p-value = 0.0001676
## P < 0.01，服药后体重比服药前降低
```

### 示例8：轮廓分析（课本例14-8）

50名硕士生和30名博士生对7个健康问卷问题的回答（每题1-4分制）。

```r
library(profileR)

data14_8 <- haven::read_sav("datasets/例14-08.sav")
data14_8 <- haven::zap_formats(data14_8)

mod <- pbg(data = data14_8[, 3:9], group = data14_8$group)
summary(mod)
## H0: Profiles are parallel
##   Wilks=0.962, P=0.821 → 两组轮廓平行
## H0: Profiles have equal levels
##   F=0.276, P=0.601 → 两组轮廓重合
## H0: Profiles are flat
##   F=10.840, P<0.001 → 总轮廓非水平线（各题得分有高有低）
```

## 结果解读指南

| 输出项 | 含义 | 解读要点 |
|--------|------|----------|
| 均值向量 (`colMeans`) | 各变量在各组的平均值 | 判断集中趋势与组间方向 |
| 协方差矩阵 (`cov`) | 对角为方差，非对角为协方差 | 方差大=个体差异大 |
| 相关矩阵 (`cor`) | 变量间的线性相关程度 | r>0.5中等，r>0.8强相关 |
| T²/Hotelling统计量 | 多变量组间差异的综合统计量 | 值越大差异越大 |
| F值（`T.2`/`approx F`） | 对T²转换后的近似F检验统计量 | 结合df得出P值 |
| Wilks Lambda | MANOVA核心统计量 | 越接近0差异越大 |
| P值 | 拒绝H0的假阳性概率 | p<0.05有统计学意义 |

**论文模板：**

- 单样本：「经单样本Hotelling T²检验，怀疑冠心病成年男性的血脂与正常男性差异有统计学意义（F=2389.8, P<0.001）。」
- 两样本：「对两组新生儿的体重和身长行Hotelling T²检验，二维空间分布差异有统计学意义（T²=9.87, F=4.58, P=0.031）。」
- MANOVA：「多变量方差分析显示三组间T细胞免疫功能差异有统计学意义（Wilks λ=0.089, F=5.50, P=0.004）。」
- 轮廓分析：「轮廓分析平行检验表明两组轮廓平行（P>0.05），相合检验表明两组轮廓重合（P>0.05），认定硕士生和博士生对健康问卷7个问题的回答无差异。」

## 常见问题与注意事项

**Q: Hotelling T²检验和多个t检验有何不同？**

A: 做m次单变量t检验会增加I类错误概率（多重比较问题）。更重要的是两者的统计问题不同：t检验反映单一变量在数轴上的组间差别，Hotelling T²反映多个变量在平面/空间上的联合差异。可能出现单变量都不显著而多变量显著的情况（例14-6），不能相互替代。

**Q: ICSNP::HotellingsT2()输出的T.2是什么？**

A: 实际是F值，不是真正的T²。`Hotelling`包输出真正的T²（Test stat），与SPSS一致。`rrcov::T2.test()`同时给出T²和F。不同包的F值可能不同，但P值相同，不影响结论。

**Q: 如何选择MANOVA检验统计量？**

A: Wilks Lambda最常用，四种统计量（Pillai、Wilks、Hotelling-Lawley、Roy）在大多数情况下结论一致。样本量小且组数多时Pillai更稳健。

**Q: 重复测量数据用重复测量ANOVA还是多变量Hotelling T²？**

A: 重复测量ANOVA要求球对称假设（Mauchly检验判断），不满足时多变量方法是推荐替代，不依赖该假设但可能需要更大样本量。

**Q: 多变量数据需满足哪些前提条件？**

A: (1) 多元正态分布 — `mshapiro.test()`或`MVN::mvn()`检验；(2) 组间协方差矩阵齐同（两组或多组比较时）；(3) 反应变量为连续型；(4) 观察对象独立。

**Q: 轮廓分析的检验顺序是什么？**

A: 必须按序：(1) 平行检验（P>0.05表示轮廓平行）→ (2) 相合检验（平行成立才有意义）→ (3) 水平检验（重合成立才有意义，检验轮廓是否为水平直线）。
