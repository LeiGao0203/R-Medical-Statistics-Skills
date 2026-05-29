---
name: medical-stat-cochran-armitage
description: "R语言医学统计：Cochran-Armitage趋势检验。检验R×2列联表中二分类结局随有序分类变量的线性趋势，如药物剂量与有效率的关系。TRIGGER when user mentions 趋势检验、剂量反应关系、有序分类变量的线性趋势、Cochran-Armitage. SKIP for 无序分类变量卡方检验、两组二分类比较、CMH检验、logistic回归纳入多个协变量、p for trend from regression models."
---

# Cochran-Armitage检验 (Cochran-Armitage Test for Trend)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

Cochran-Armitage检验用于检验R×2列联表中，二分类结局变量随有序分类自变量是否存在线性趋势。

**适用情况：**
- 自变量为有序分类变量（如药物剂量等级、年龄分组、疾病严重程度分级）
- 因变量为二分类变量（如有效/无效、阳性/阴性、死亡/存活）
- 研究目的是检验"随着自变量等级增加，结局发生率是否存在线性变化趋势"

**不适用的场景：**
- 自变量为无序分类变量（如血型、民族）→ 使用Pearson卡方检验或Fisher精确检验
- 仅有两组二分类（2×2表）→ 使用四格表卡方检验
- 存在分层/混杂因素 → 使用Cochran-Mantel-Haenszel (CMH) 检验
- 需要同时控制多个协变量 → 使用logistic回归
- 自变量为连续型 → 直接使用logistic回归或线性回归

**医学研究常见应用：**
- 不同药物剂量组的有效/无效率趋势分析
- 疾病严重程度分级与病死率关系
- 年龄分组与高血压患病率关系
- 吸烟严重程度（不吸烟/轻度/中度/重度）与肺癌发病率关系

## 前置条件

**R包：**
```r
install.packages("DescTools")
```

**数据格式：**
- 一个R×2的频数矩阵（matrix/table），行名为有序分类水平，列名为二分类结局
- 或R×2的列联表（table对象）

**统计假设：**
- 观测值相互独立
- 自变量水平有序（ordered categories with natural ordering）
- 检验的是线性趋势，不检测非线性模式

## 方法选择决策树

```
你的数据情况 →
├── 自变量有序分类 + 二分类结局，检验线性趋势
│   └── 使用 Cochran-Armitage 趋势检验
│
├── 自变量有序分类 + 二分类结局，还需控制协变量
│   └── 使用 logistic回归（将自变量编码为有序因子或数值型）
│
├── 自变量无序分类 + 二分类结局
│   └── 使用 Pearson卡方检验 或 Fisher精确检验
│
├── 自变量分类 + 二分类结局，存在分层变量
│   └── 使用 Cochran-Mantel-Haenszel (CMH) 检验
│
└── 自变量连续型 + 二分类结局
    └── 使用 logistic回归，纳入自变量为连续变量
```

**Cochran-Armitage vs CMH vs Logistic回归 的区分：**

| 方法 | 核心问题 | 自变量类型 | 协变量 |
|------|---------|-----------|--------|
| Cochran-Armitage | 是否存在线性趋势？ | 有序分类 | 无 |
| CMH检验 | 控制分层后是否有关联？ | 分类 | 一个分层变量 |
| Logistic回归 | 多因素影响结局？ | 连续/分类 | 任意多个 |

## 标准工作流

### 步骤1：数据准备与探索

将原始数据整理为R×2频数矩阵，行按等级升序排列，并按剂量分组计算粗率以初步判断趋势方向。

```r
df <- matrix(c(13, 136, 17, 125, 16, 104, 32, 149, 9, 45),
             nrow = 5, byrow = TRUE,
             dimnames = list(
               "Dose" = c("50", "100", "200", "300", "500"),
               "effect" = c("Yes", "No")))
```

### 步骤2：描述性分析

计算各等级的阳性率，观察是否存在单调变化趋势：

```r
df[, 1] / rowSums(df)
##         50        100        200        300        500
## 0.08724832 0.11971831 0.13333333 0.17679558 0.16666667
```

### 步骤3：执行Cochran-Armitage检验

```r
library(DescTools)
CochranArmitageTest(df)
```

### 步骤4：结果解读

看Z统计量和p值。p < 0.05表示存在统计学意义上的线性趋势，即随着自变量等级变化，二分类结局的发生率呈现单调递增或递减的变化模式。

### 步骤5：结果报告（论文中的统计描述）

> 采用Cochran-Armitage趋势检验分析不同剂量组间有效率的线性趋势，结果显示Z = 2.2116, p = 0.027，差异有统计学意义，提示药物有效率随剂量增加而升高。

## 代码示例

### 示例：药物剂量与疗效趋势检验

不同剂量水平下药物有效率的Cochran-Armitage趋势检验。

```r
library(DescTools)

# 构建R×2频数矩阵
df <- matrix(c(13, 136, 17, 125, 16, 104, 32, 149, 9, 45),
             nrow = 5, byrow = TRUE,
             dimnames = list(
               "Dose" = c("50", "100", "200", "300", "500"),
               "effect" = c("Yes", "No")))
df
##      effect
## Dose  Yes  No
##   50   13 136
##   100  17 125
##   200  16 104
##   300  32 149
##   500   9  45

# 计算各剂量组有效率
df[, 1] / rowSums(df)
##         50        100        200        300        500
## 0.08724832 0.11971831 0.13333333 0.17679558 0.16666667

# Cochran-Armitage趋势检验（默认双侧检验）
DescTools::CochranArmitageTest(df)
##
##  Cochran-Armitage test for trend
##
## data:  df
## Z = 2.2116, dim = 5, p-value = 0.02699
## alternative hypothesis: two.sided
```

### 与logistic回归的对比

将频数数据展开为个体水平数据后，使用logistic回归验证：

```r
# 频数表转为个案数据
df1 <- rstatix::counts_to_cases(df)

# Dose转为数值型（编码为1,2,3,4,5）
df1$Dose <- as.numeric(factor(df1$Dose))

summary(glm(effect ~ Dose, data = df1, family = binomial()))
##
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)  2.48493    0.29598   8.396   <2e-16 ***
## Dose        -0.21544    0.08985  -2.398   0.0165 *
##
##     Null deviance: 510.57  on 645  degrees of freedom
## Residual deviance: 504.71  on 644  degrees of freedom
## AIC: 508.71
```

logistic回归中Dose的p = 0.0165，与Cochran-Armitage检验的p = 0.027方向一致，均提示剂量与有效性存在线性关联。

## 结果解读指南

`CochranArmitageTest()` 返回内容解读：

| 输出项 | 含义 | 解读 |
|--------|------|------|
| **Z** | 标准正态检验统计量 | 正值表示趋势上升，负值表示趋势下降 |
| **dim** | 有序分类的级别数 | 即自变量的水平数（R） |
| **p-value** | 双侧检验p值 | p < 0.05 表示线性趋势显著 |
| **alternative hypothesis** | 备择假设方向 | 默认"two.sided"；可指定"greater"/"less"用于单侧趋势检验 |

**论文撰写要点：**
- 先报告各水平频数和百分比，再报告趋势检验结果
- 描述趋势方向（递增/递减）
- p值通常报告到小数点后3位
- 若趋势不显著，应避免过度解读，并注明可能原因是样本量不足或趋势确实不存在

## 常见问题与注意事项

**Q1: Cochran-Armitage检验与卡方检验的主要区别是什么？**

卡方检验检验的是"不同组间结局分布是否有差异"，不关心组间顺序。Cochran-Armitage专门检验"结局是否随自变量的有序等级呈线性变化趋势"。卡方检验检出的是任意差异（包括V形、倒V形等非线性模式），CA检验只检线性趋势。通常CA检验的检验效能比卡方检验更高（针对线性趋势）。

**Q2: 行变量的编码方式会影响结果吗？**

会。默认情况下 `CochranArmitageTest()` 使用等间距评分（1, 2, 3, ...）。如果自变量的间距不均匀（如剂量50, 100, 200, 300, 500），可以使用`weights`参数指定自定义评分权重。如果希望使用实际剂量值作为评分，可传入 `scores = list(dose_values)`。

**Q3: 可以指定单侧检验吗？**

可以。使用 `alternative = "greater"` 或 `alternative = "less"` 参数。例如，如果事先假设剂量越高有效率越高，可用单侧检验以增加检验效能。但应在方案阶段确定，不应根据数据选择。

**Q4: CA检验与logistic回归中的趋势检验p值为何不一样？**

CA检验基于列联表数据的趋势卡方分解，LOGISTIC回归基于最大似然估计。两者渐近等价，但小样本时可能有差异。关键区别：logistic回归可纳入协变量，CA检验仅用于单一有序变量。

**Q5: 数据中有0频数细胞如何处理？**

如果某个剂量组的有效/无效为0，CA检验仍然可用（只要不是整行或整列为0）。但样本量过小时检验效能不足。对于总样本量很小的情况，可考虑使用精确趋势检验方法。
