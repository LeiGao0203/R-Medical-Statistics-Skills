---
name: medical-stat-three-line-table
description: "R语言医学统计：三线表（基线特征表/Table 1）绘制。使用compareGroups、gtsummary、table1等R包生成符合中文期刊格式要求的三线表，自动按分组描述连续变量和分类变量，可直接导出为Word、Excel、LaTeX、PDF等格式。TRIGGER when user mentions 三线表、Table1、基线特征表、描述性表格、compareGroups、table1包、gtsummary，or asks about creating publication-ready descriptive tables。SKIP for 统计检验本身（如t检验、卡方检验应触发对应skill）、ROC曲线绘制、统计图形绘制。"
---

# 三线表绘制 (Three-line Table / Table 1)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

三线表是医学论文中最常见的表格类型，特别是**基线资料表（Table 1）**，用于展示研究对象的人口学特征和临床特征在各分组之间的分布情况。

**适用条件：**
- 需要制作SCI论文或中文期刊投稿用的基线特征表
- 数据包含连续变量（需报告均值±标准差或中位数/IQR）和分类变量（需报告频数和百分比）
- 有明确的分组变量（如治疗组/对照组、疾病分期）
- 需要在表格中自动附上组间比较的P值
- 需要将表格导出为Word、Excel、LaTeX等格式供进一步编辑

**不适用条件：**
- 只需要做统计检验但不需要制表 → 使用对应的统计检验skill（t检验、卡方检验等）
- 需要复杂的多因素回归结果表格 → 使用Logistic回归或Cox回归对应的skill
- 绘制统计图形（柱状图、箱线图等）→ 使用统计绘图skill
- 需要ROC曲线分析 → 使用ROC曲线skill

**常见R包：** compareGroups、tableone、table1、gtsummary、gt、gtExtras等。本文以 compareGroups 为主，其在多数组况下最便捷。

## 前置条件

### 安装R包

```r
install.packages("compareGroups")
install.packages("table1")       # 备选方案
install.packages("gtsummary")    # 备选方案
install.packages("tableone")     # 备选方案
```

### 数据格式要求

- 数据框格式，行=观测，列=变量
- **分类变量必须预先因子化**（`as.factor()`），否则compareGroups会自动判断（取值数≤5视为分类变量，可通过`min.dis`调整）
- 生存数据（time-to-event）需用`survival::Surv()`包装
- 建议给变量设置`attr(x, "label")`属性，表格中会自动显示label
- 数据框中仅保留需要描述的变量，不需要的变量建议过滤掉

## 方法选择决策树

```
数据情况和需求 →
├── 仅需快速生成基线描述表（有明确分组变量）
│   ├── 偏好简洁语法，支持导出多种格式 → compareGroups::compareGroups() + createTable()
│   ├── 需要精细化控制表格外观（频率/百分比、小数位数等）→ compareGroups（通过createTable参数控制）
│   └── 需要分层表格（如按性别分层）→ strataTable() 或 cbind() 合并
├── 需要同时做单因素分析并展示OR/HR → compareGroups (show.ratio = TRUE)
├── 遵从FDA/CDISC标准，或偏好tidyverse生态 → gtsummary::tbl_summary()
├── 简单数据，不想引入额外依赖 → base R手动计算 + knitr::kable()
└── 连续变量呈现方式选择
    ├── 正态分布 → 均值±标准差（mean±SD，method = 1）
    └── 非正态分布 → 中位数（Q1, Q3），即 median(IQR)，method = 2 或 NA（自动用shapiro.test判定）
```

## 标准工作流

### 步骤1：数据准备

```r
library(compareGroups)
data("regicor")
dim(regicor)  # 2294行 × 25列

# 分类变量确保为factor
# 生存数据需包装为Surv对象
library(survival)
regicor$tmain <- with(regicor, Surv(tocv, cv == "Yes"))
attr(regicor$tmain, "label") <- "Time to CV event or censoring"
```

### 步骤2：执行分组描述

```r
# 基本用法：formula格式，左侧为分组变量，右侧为待描述变量
res <- compareGroups(year ~ age + sex + smoker + bmi, data = regicor)
res
```

compareGroups自动判断变量类型：
- 数值型变量（>5个不同取值）→ 视为连续变量，默认正态分布，做方差分析/ANOVA
- 因子型变量或少数值变量（≤5个不同取值）→ 视为分类变量，做卡方检验
- `method = NA` 时通过 shapiro.test 自动判断是否正态，`alpha`参数设判定阈值

### 步骤3：构建/定制表格

```r
restab <- createTable(res)

# 常用定制参数：
createTable(res,
  hide.no = "no",              # 隐藏含"no"的分类水平（如只显示"有"的类别）
  hide = c(sex = "Male"),      # 隐藏指定的分类水平
  digits = c(age = 2, bmi = 1),# 控制小数位数
  type = 2,                    # 1=仅百分比, 2=计数+百分比(默认), 3=仅计数
  show.all = TRUE,             # 显示全部人群列
  show.p.overall = TRUE,       # 显示组间比较P值
  show.p.trend = TRUE,         # 分组≥3类时显示趋势检验P值
  show.p.mul = TRUE,           # 显示两两比较P值
  show.ratio = TRUE,           # 二分类因变量时展示OR或HR
  show.n = FALSE               # 是否显示每个变量的有效样本量
)
```

### 步骤4：导出表格

```r
# 支持多种格式导出
export2csv(restab, file = "table1.csv")
export2word(restab, file = "table1.docx")
export2xls(restab, file = "table1.xlsx")
export2html(restab, file = "table1.html")
export2latex(restab, file = "table1.tex")
export2pdf(restab, file = "table1.pdf")
export2md(restab, file = "table1.md")
```

### 步骤5：论文中的统计描述写法

```
标准论文表述示例：

"Table 1 shows the baseline characteristics of the study participants.
Continuous variables are presented as mean ± standard deviation or
median (interquartile range), and categorical variables are presented
as frequency (percentage). Group comparisons were performed using
ANOVA or Kruskal-Wallis test for continuous variables and chi-squared
test for categorical variables."
```

## 代码示例

### 示例1：基本基线资料表

```r
library(compareGroups)
data("regicor")

# 以year为分组变量，描述所有变量（排除id列）
compareGroups(year ~ . - id, data = regicor)
## -------- Summary of results by groups of 'Recruitment year'---------
##    var           N    p.value  method            selection
## 1  Age           2294 0.078*   continuous normal ALL
## 2  Sex           2294 0.506    categorical       ALL
## 3  Smoking status 2233 <0.001** categorical       ALL
## ...

# 用createTable呈现完整表格
res <- compareGroups(year ~ age + sex + smoker + bmi + sbp, data = regicor,
                     selec = list(sbp = txhtn == "No"))
restab <- createTable(res)
print(restab, which.table = "descr")
## --------Summary descriptives table by 'Recruitment year'---------
##                               1995        2000        2005     p.overall
##                               N=431       N=786      N=1077
## Age                        54.1 (11.7) 54.3 (11.2) 55.3 (10.6)   0.078
## Sex:                                                             0.506
##     Male                   206 (47.8%) 390 (49.6%) 505 (46.9%)
##     Female                 225 (52.2%) 396 (50.4%) 572 (53.1%)
## Smoking status:                                                 <0.001
##     Never smoker           234 (56.4%) 414 (54.6%) 553 (52.2%)
##     Current or former < 1y 109 (26.3%) 267 (35.2%) 217 (20.5%)
##     Former >= 1y           72 (17.3%)  77 (10.2%)  290 (27.4%)
## Body mass index            27.0 (4.15) 28.1 (4.62) 27.6 (4.63)  <0.001
## Systolic blood pressure    129 (17.4)  130 (20.1)  124 (16.9)   <0.001
```

### 示例2：连续变量方法选择（正态 vs 非正态）

```r
# method: 1=正态(mean±SD), 2=非正态(median[P25;P75]), NA=自动判断
res <- compareGroups(year ~ age + triglyc, data = regicor,
                     method = c(triglyc = NA),  # NA=shapiro.test自动判断
                     alpha = 0.01)
res
##   var          N    p.value  method                selection
## 1 Age          2294 0.078*   continuous normal     ALL
## 2 Triglycerides 2231 0.762    continuous non-normal ALL

# 手动指定triglyc为非正态，并自定义百分位数
res2 <- compareGroups(year ~ age + triglyc, data = regicor,
                      method = c(triglyc = 2),
                      Q1 = 0.025, Q3 = 0.975)  # 显示95%参考值范围
createTable(res2)
##                    1995            2000            2005       p.overall
## Age             54.1 (11.7)     54.3 (11.2)     55.3 (10.6)     0.078
## Triglycerides 94.0 [47.0;292] 98.0 [47.0;278] 98.0 [42.0;293]   0.762
```

### 示例3：数据子集选择

```r
# sub集：只选女性；selec：bmi只选age>50的观测
compareGroups(year ~ age + smoker + bmi, data = regicor,
              selec = list(bmi = age > 50),
              subset = sex == "Female")
##   var             N    p.value  method            selection
## 1 Age             1193 0.351    continuous normal sex == "Female"
## 2 Smoking status  1162 <0.001** categorical       sex == "Female"
## 3 Body mass index  709 0.308    continuous normal (sex == "Female") & (age > 50)
```

### 示例4：展示OR值和调整参考水平

```r
# 二分类因变量，展示OR和p.ratio
res1 <- compareGroups(cv ~ age + sex + bmi + smoker, data = regicor, ref = 1)
createTable(res1, show.ratio = TRUE)
##                                 No          Yes            OR        p.ratio p.overall
## Age                        54.6 (11.1)  57.5 (11.0) 1.02 [1.00;1.04]  0.017    0.018
## Sex:                                                                           0.801
##     Male                   996 (48.1%)  46 (50.0%)        Ref.        Ref.
##     Female                 1075 (51.9%) 46 (50.0%)  0.93 [0.61;1.41]  0.721
## Smoking status:                                                               <0.001
##     Never smoker           1099 (54.3%) 37 (40.2%)        Ref.        Ref.
##     Current or former < 1y 506 (25.0%)  47 (51.1%)  2.75 [1.77;4.32] <0.001

# 自定义参考水平
res2 <- compareGroups(cv ~ age + sex + bmi + smoker, data = regicor,
                      ref = c(smoker = 1, sex = 2))
# ref.no = "NO": 自动以含有"No"的水平为参考
res3 <- compareGroups(cv ~ age + sex + bmi + histhtn + txhtn, data = regicor,
                      ref.no = "NO")
```

### 示例5：连续性变量OR的单位调整

```r
# age每增加10岁、bmi每增加2个单位时的OR
res <- compareGroups(cv ~ age + bmi, data = regicor,
                     fact.ratio = c(age = 10, bmi = 2))
createTable(res, show.ratio = TRUE)
##                     No          Yes            OR        p.ratio p.overall
## Age             54.6 (11.1) 57.5 (11.0) 1.26 [1.04;1.53]  0.017    0.018
## Body mass index 27.6 (4.56) 28.1 (4.48) 1.05 [0.96;1.14]  0.313    0.307
```

### 示例6：生存数据 (HR值)

```r
regicor$tmain <- with(regicor, Surv(tocv, cv == "Yes"))
attr(regicor$tmain, "label") <- "Time to CV event or censoring"

createTable(compareGroups(tmain ~ year + age + sex, data = regicor),
            show.ratio = TRUE)
##                     No event      Event           HR        p.ratio p.overall
## Recruitment year:                                                     0.157
##     1995          388 (18.7%)  10 (10.9%)        Ref.        Ref.
##     2000          706 (34.1%)  35 (38.0%)  1.95 [0.96;3.93]  0.063
##     2005          977 (47.2%)  47 (51.1%)  1.82 [0.92;3.59]  0.087
## Age               54.6 (11.1)  57.5 (11.0) 1.02 [1.00;1.04]  0.021    0.021
## Sex:                                                                  0.696
##     Male          996 (48.1%)  46 (50.0%)        Ref.        Ref.
##     Female        1075 (51.9%) 46 (50.0%)  0.92 [0.61;1.39]  0.696

# p.overall 由 logrank检验计算
survdiff(Surv(tocv, cv == "Yes") ~ sex, data = regicor)
## Chisq= 0.2  on 1 degrees of freedom, p= 0.7

# p.ratio 由 Cox回归计算
aa <- coxph(Surv(tocv, cv == "Yes") ~ sex, data = regicor)
broom::tidy(aa)
## term      estimate std.error statistic p.value
## sexFemale  -0.0814     0.209    -0.390   0.696
```

### 示例7：提取P值进行多重校正

```r
data(SNPs)
tab <- createTable(compareGroups(casco ~ snp10001 + snp10002 + snp10005 +
                                  snp10008 + snp10009, SNPs))
pvals <- getResults(tab, "p.overall")
p.adjust(pvals, method = "BH")
##  snp10001  snp10002  snp10005  snp10008  snp10009
## 0.7051300 0.7072158 0.7583432 0.7583432 0.7072158

# 或在compareGroups阶段直接校正
cg <- compareGroups(casco ~ snp10001 + snp10002 + snp10005 + snp10008 + snp10009, SNPs)
createTable(padjustCompareGroups(cg, method = "BH"))
```

### 示例8：表格合并（行合并 / 列分层）

```r
# 行合并 — 按主题分组展示
restab1 <- createTable(compareGroups(year ~ age + sex, data = regicor))
restab2 <- createTable(compareGroups(year ~ bmi + smoker, data = regicor))
rbind("Non-modifiable risk factors" = restab1,
      "Modifiable risk factors" = restab2)
## Non-modifiable risk factors:
##     Age                        54.1 (11.7) 54.3 (11.2) 55.3 (10.6)   0.078
## Modifiable risk factors:
##     Body mass index            27.0 (4.15) 28.1 (4.62) 27.6 (4.63)  <0.001
##     Smoking status:                                                 <0.001

# 分层表格 — strataTable
res <- compareGroups(year ~ age + bmi + smoker + histchol + histhtn, regicor)
restab <- createTable(res, hide.no = "no")
strataTable(restab, "sex")
##                                    Male                              Female
##                   1995     2000     2005     p.overall  1995     2000     ...
```

### 示例9：gtsummary备选方案

```r
library(gtsummary)
library(dplyr)

# 基础用法
regicor %>%
  select(age, sex, smoker, bmi, sbp, year) %>%
  tbl_summary(
    by = year,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1
  ) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**")
```

### 示例10：table1备选方案

```r
library(table1)

table1(~ age + bmi + sbp + sex + smoker | year, data = regicor,
       overall = "Total",
       render.continuous = c(. = "Mean (SD)", . = "Median [Q1, Q3]"))
```

### 示例11：base R手动构建

```r
# 连续变量：按分组计算mean±SD
aggregate(age ~ year, data = regicor,
          FUN = function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
##   year        age
## 1 1995 54.1 (11.7)
## 2 2000 54.3 (11.2)
## 3 2005 55.3 (10.6)

# 分类变量：按分组计算频率和百分比
tab <- table(regicor$sex, regicor$year)
prop <- prop.table(tab, 2) * 100
cbind(tab, sprintf("n=%d (%.1f%%)", tab[,1], prop[,1]))
```

### 示例12：一键函数descrTable

```r
# descrTable = compareGroups + createTable 的二合一函数
descrTable(year ~ age + bmi + smoker + histchol + histhtn, data = regicor)
##                               1995        2000        2005     p.overall
##                               N=431       N=786      N=1077
## Age                        54.1 (11.7) 54.3 (11.2) 55.3 (10.6)   0.078
## Body mass index            27.0 (4.15) 28.1 (4.62) 27.6 (4.63)  <0.001
## Smoking status:                                                 <0.001
##     Never smoker           234 (56.4%) 414 (54.6%) 553 (52.2%)
##     Current or former < 1y 109 (26.3%) 267 (35.2%) 217 (20.5%)
##     Former >= 1y           72 (17.3%)  77 (10.2%)  290 (27.4%)
## History of hyperchol.:                                          <0.001
## History of hypertension:                                        <0.001
```

## 结果解读指南

### 表格核心列解读

| 列名 | 含义 |
|------|------|
| `[分组名]` 各列 | 每列的格式取决于变量类型：连续变量显示「均值 (标准差)」或「中位数 [P25; P75]」；分类变量显示「频数 (百分比)」 |
| `p.overall` | 组间总体比较的P值，连续变量用ANOVA/Kruskal-Wallis，分类变量用卡方检验 |
| `p.trend` | 趋势检验P值（仅分组≥3类时），正态用Pearson相关、非正态用Spearman相关 |
| `p.ratio` | 单因素logistic/Cox回归中该变量的P值 |
| `OR / HR [95% CI]` | 比值比/风险比及其95%置信区间 |

### P值含义

- **p.overall < 0.05** → 该变量在各组间的分布差异有统计学意义
- **p.trend < 0.05** → 该变量在不同分组之间存在线性趋势
- **p.ratio < 0.05** → 该变量对结局的影响有统计学意义
- **OR > 1** → 危险因素；**OR < 1** → 保护因素（需结合参考水平理解）

### 验证工具

compareGroups的结果可以被标准R函数验证：

```r
# p.overall 对应：连续变量=ANOVA，分类变量=卡方检验
summary(aov(age ~ year, data = regicor))     # p=0.0778 —— 与compareGroups一致
chisq.test(regicor$sex, regicor$year)         # p=0.5056 —— 与compareGroups一致

# p.ratio + OR 对应：二分类=逻辑回归
aa <- glm(cv ~ age, data = regicor, family = binomial())
broom::tidy(aa, exponentiate = TRUE, conf.int = TRUE)
## age: OR=1.02 (1.00-1.04), p=0.017 —— 与compareGroups一致

# 生存数据：p.overall=logrank, p.ratio+HR=Cox回归
```

## 常见问题与注意事项

**Q1: compareGroups vs tableone / table1 / gtsummary 如何选择？**
compareGroups功能最全面（OR/HR、趋势检验、多重校正），语法简洁，导出格式丰富。gtsummary更适合tidyverse生态用户，出图更美观。table1简单轻量。新手推荐compareGroups或gtsummary。

**Q2: 为什么某些变量的method显示为continuous normal但我认为是非正态？**
默认`method = 1`（正态分布），可通过`method = NA`让包自动用Shapiro-Wilk检验判断，或手动指定`method = 2`（非正态）。

**Q3: 分类变量没有被识别为categorical？**
非因子型向量取值≤5个才自动视为分类变量。可通过`min.dis`调整阈值，或直接用`factor()`转换。

**Q4: 如何隐藏二分类变量中无意义的类别（如只显示"Yes"不显示"No"）？**
使用 `hide.no = "no"` 参数（不区分大小写），或 `hide = c(variable = "category_name")`。

**Q5: 表格输出到Word后还需要手动调整吗？**
目前没有任何R包能输出无需修改即可直接发表的表格。`export2word()`输出的表格通常需要在Word中微调对齐、行高、字体等。

**Q6: 变量有缺失值时如何处理？**
compareGroups自动处理缺失值。通过`show.n = TRUE`可查看每个变量的有效样本量。如某分组某类别样本数为0，设置`simplify = FALSE`来抑制warning并跳过检验。

**Q7: SPSS与R的三线表制作有何差异？**
SPSS的`CTABLES`过程可交互式定制表格，对新手更友好但不易复现。R的compareGroups/gtsummary更利于可重复研究，但最终修饰仍需在Word中完成。

**Q8: 如何制作不分组的描述性表格？**
去除formula左侧的分组变量即可，或使用`gtsummary::tbl_summary()`不指定`by`参数。
