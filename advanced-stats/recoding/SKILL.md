---
name: medical-stat-recoding
description: "R语言医学统计：分类变量（多分类无序/有序变量）在回归模型中的编码方案。涵盖哑变量编码（dummy coding）、效应编码（effect coding）、简单编码（simple coding）、正交多项式编码、Helmert编码、Forward/Backward差分编码，以及如何设置参照组。TRIGGER when user mentions 哑变量、虚拟变量、编码方案、dummy coding、参照组设置、分类变量如何放入回归、contr.treatment、contr.sum、contr.helmert、contr.poly，or asks about encoding categorical predictors in regression models. SKIP for 连续变量转换、数据清洗、变量标准化或归一化。"
---

# 分类变量重编码 (Categorical Variable Recoding)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

- 回归分析（线性回归、logistic回归、Cox回归）中需将多分类自变量（K > 2）放入模型
- 需理解不同类别对应的回归系数含义
- 需自定义参照组与其他类别进行比较
- 需比较某一类别与所有其他类别均值的差异（效应编码）
- 有序分类变量需检验线性趋势（正交多项式编码）
- 需检验相邻类别间的递进差异（差分编码）

不适用场景：
- 二分类变量（仅2个类别），直接设为0/1即可
- 连续变量转换（标准化、归一化、对数变换等）
- 纯粹数据清洗工作（合并类别、重命名标签）

常见医学应用：种族、血型、疾病分期、病理分级、治疗方案分组等作为协变量进入回归模型时。

## 前置条件

**R包**：所有编码方案均使用 R 基础包 stats 的函数，无需安装第三方包。
- `contr.treatment()` — 哑变量编码（普通因子默认）
- `contr.sum()` — 效应编码（Deviation coding）
- `contr.poly()` — 正交多项式编码（有序因子默认）
- `contr.helmert()` — 反 Helmert 编码

**数据格式**：自变量需转换为因子型。有序因子需设 `ordered = TRUE`。

**前提假设**：编码方案本身无统计假设，但所在回归模型需满足对应前提（如线性回归需正态、等方差等）。

## 方法选择决策树

```
分类变量要进入回归模型 →
├── 无序分类（如种族、血型） →
│   ├── 各水平与参照组比较 → 哑变量编码 (contr.treatment)
│   ├── 各水平与总体均值比较 → 效应编码 (contr.sum)
│   ├── 截距为总均值，系数为组间差 → 简单编码 (simple coding)
│   ├── 某水平与后面所有水平比较 → Helmert 编码
│   ├── 某水平与前面所有水平比较 → 反 Helmert (contr.helmert)
│   │   ├── 与相邻下一水平比较 → Forward 差分
│   │   └── 与相邻上一水平比较 → Backward 差分
│   ├── 某水平与相邻水平比较 → Forward/Backward 差分
└── 有序分类（如疾病分期） →
    ├── 检验线性趋势 → 正交多项式编码 (contr.poly，有序因子默认)
    └── 当无序处理 → 选用上述编码之一
```

医学统计 90% 以上的场景使用哑变量编码，R 对数值型和普通因子默认使用此方案。

## 标准工作流

### 步骤1：数据准备与探索

```r
# 转为因子型并查看各水平因变量均值
hsb2$race.f <- factor(hsb2$race, labels = c("Hispanic", "Asian", "African-Am", "Caucasian"))
tapply(hsb2$write, hsb2$race.f, mean)
## Hispanic    Asian African-Am  Caucasian 
## 46.45833  58.00000  48.20000  54.05517
```

### 步骤2：选择并设置编码方案

```r
contrasts(hsb2$race.f)                 # 查看当前编码
contrasts(hsb2$race.f) <- contr.treatment(4)  # 设置为哑变量
contrasts(hsb2$race.f) <- contr.sum(4)        # 设置为效应编码
```

### 步骤3：执行回归分析

```r
summary(lm(write ~ race.f, data = hsb2))
```

### 步骤4-5：解读系数并报告

系数含义取决于编码方案（见下表）。论文中以哑变量为例："以Hispanic为参照，Asian写作成绩显著更高（β=11.54, SE=3.29, p=0.001）。"

## 代码示例

### 演示数据

```r
load(file = "datasets/codingSchemes.rdata")
hsb2$race.f <- factor(hsb2$race, labels = c("Hispanic", "Asian", "African-Am", "Caucasian"))
# 各组均值: Hispanic=46.46, Asian=58.00, African-Am=48.20, Caucasian=54.06
```

### 哑变量编码 (Dummy Coding) — contr.treatment()

K个类别产生K-1个哑变量。截距=参照组均值，系数=该类别均值-参照组均值。

```r
contr.treatment(4)
##   2 3 4
## 1 0 0 0
## 2 1 0 0
## 3 0 1 0
## 4 0 0 1

contrasts(hsb2$race.f) <- contr.treatment(4)
summary(lm(write ~ race.f, data = hsb2))
## (Intercept)   46.458   ...  < 2e-16 ***
## race.f2       11.542   ...  0.000552 ***    # 58 - 46.458
## race.f3        1.742   ...  0.524613        # 48.2 - 46.458
## race.f4        7.597   ...  0.000179 ***    # 54.055 - 46.458
```

### 简单编码 (Simple Coding)

唯一区别：截距 = 总均值（非参照组均值）。系数与哑变量相同。K个类别时参照组编码为 -1/K，比较组编码为 (K-1)/K。

```r
c <- contr.treatment(4)
my.simple <- c - matrix(rep(1/4, 12), ncol = 3)
contrasts(hsb2$race.f) <- my.simple
summary(lm(write ~ race.f, data = hsb2))
## (Intercept)  51.6784  # (46.46+58+48.2+54.06)/4
```

### 效应编码 (Deviation Coding) — contr.sum()

比较每个类别均值与总均值。截距=总均值，系数=该类别均值-总均值。

```r
contr.sum(4)
##   [,1] [,2] [,3]
## 1    1    0    0
## 2    0    1    0
## 3    0    0    1
## 4   -1   -1   -1

contrasts(hsb2$race.f) <- contr.sum(4)
summary(lm(write ~ race.f, data = hsb2))
## (Intercept)  51.6784  ...
## race.f1      -5.2200  0.00160 **    # 46.458 - 51.678
## race.f2       6.3216  0.00384 **    # 58 - 51.678
## race.f3      -3.4784  0.04602 *     # 48.2 - 51.678
```

### 正交多项式编码 (Orthogonal Polynomial) — contr.poly()

有序因子的默认编码。检验线性趋势(.L)、二次(.Q)、三次(.C)。

```r
hsb2$readcat <- cut(hsb2$read, 4, ordered = TRUE)
contrasts(hsb2$readcat)  # 默认即正交多项式编码
##              .L   .Q         .C
## [1,] -0.6708204  0.5 -0.2236068
## [2,] -0.2236068 -0.5  0.6708204
## [3,]  0.2236068 -0.5 -0.6708204
## [4,]  0.6708204  0.5  0.2236068

summary(lm(write ~ readcat, data = hsb2))
## readcat.L    14.2587  <2e-16 ***    # 显著线性趋势
## readcat.Q    -0.9680  0.446         # 无二次趋势
## readcat.C    -0.1554  0.877         # 无三次趋势
```

### Helmert / 反 Helmert 编码

Helmert：当前水平与后面所有水平比较。反 Helmert (`contr.helmert()`)：当前水平与前面所有水平比较。

```r
my.helmert <- matrix(c(3/4,-1/4,-1/4,-1/4, 0,2/3,-1/3,-1/3, 0,0,1/2,-1/2), ncol=3)
contrasts(hsb2$race.f) <- my.helmert
summary(lm(write ~ race.f, hsb2))
## race.f1  -6.96  # 46.46 - (58+48.2+54.06)/3

contr.helmert(4)
##   [,1] [,2] [,3]
## 1   -1   -1   -1
## 2    1   -1   -1
## 3    0    2   -1
## 4    0    0    3
```

### Forward / Backward 差分编码

Forward：当前水平与相邻下一水平比较。Backward：当前水平与相邻上一水平比较（符号相反）。

```r
my.forward.diff <- matrix(c(3/4,-1/4,-1/4,-1/4, 1/2,1/2,-1/2,-1/2, 1/4,1/4,1/4,-3/4), ncol=3)
contrasts(hsb2$race.f) <- my.forward.diff
summary(lm(write ~ race.f, data = hsb2))
## race.f1 -11.542  # 46.458 - 58 (Hispanic vs Asian)
## race.f2   9.800  # 58 - 48.2 (Asian vs African-Am)
## race.f3  -5.855  # 48.2 - 54.055 (African-Am vs Caucasian)
```

## 结果解读指南

| 编码方案 | 截距含义 | 回归系数含义 | R函数 |
|---|---|---|---|
| 哑变量 | 参照组均值 | 该类别 - 参照组 | `contr.treatment()` |
| 简单编码 | 总均值 | 该类别 - 参照组 | 手动构造 |
| 效应编码 | 总均值 | 该类别 - 总均值 | `contr.sum()` |
| 正交多项式 | 总均值 | .L=线性 .Q=二次 .C=三次 | `contr.poly()` |
| Helmert | 总均值 | 该类别 - 后续合并均值 | 手动构造 |
| 反Helmert | 总均值 | 该类别 - 前面合并均值 | `contr.helmert()` |
| Forward差分 | 总均值 | 该类别 - 下一类别 | 手动构造 |
| Backward差分 | 总均值 | 该类别 - 上一类别 | 手动构造 |

## 常见问题与注意事项

**Q1：如何选择编码方案？**
有明确参照组（安慰剂 vs 治疗组）→ 哑变量；无自然参照组、关注偏离总体均值 → 效应编码；有序分类且等距 → 正交多项式。

**Q2：如何改变参照组？**
`hsb2$race.f <- relevel(hsb2$race.f, ref = "Caucasian")` 改变因子水平顺序即可。

**Q3：logistic/Cox回归也适用吗？**
是。编码方案是因子型在回归中的通用机制，`glm()` 和 `coxph()` 中效果相同。

**Q4：SPSS编码对应关系？**
Indicator=哑变量，Simple=简单编码，Deviation=效应编码，Helmert=反Helmert，Difference=Forward差分，Polynomial=正交多项式。原理一致。

**注意事项**：
1. 分类变量务必先转因子型，否则 R 将其当连续值得到错误单系数
2. 编码矩阵列数固定为 K-1（K 为类别数）
3. 同一模型不同变量可使用不同编码方案，互不影响
4. 编码不改变模型整体拟合（R²、F检验），只改变系数含义和检验对象
5. 正交多项式仅适用于有序等距分类变量

参考资料：https://stats.oarc.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
