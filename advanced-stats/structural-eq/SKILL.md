---
name: medical-stat-structural-eq
description: "R语言医学统计：结构方程模型（SEM）。涵盖验证性因子分析（CFA）、路径分析、结构方程建模，使用lavaan包进行模型设定、拟合与评估。TRIGGER when user mentions 结构方程模型、SEM、验证性因子分析、CFA、路径分析、lavaan、潜变量、测量模型、结构模型、路径图、构念效度、模型拟合指数，or asks about testing theoretical models with latent variables. SKIP for 探索性因子分析（EFA）、主成分分析（PCA）、简单回归、t检验、方差分析."
---
# 结构方程模型 (Structural Equation Modeling, SEM)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用结构方程模型：**
- 研究中有**不能直接测量的潜变量**（如学习能力、幸福感、焦虑、工作满意度），需要通过多个观测指标间接测量
- 需要同时处理多个变量之间的复杂关系（多个因变量、中介效应、间接效应）
- 评价量表或问卷的**构念效度**（结构效度），即验证观测指标是否真的测量了预设的潜变量——此时使用**验证性因子分析（CFA）**
- 检验潜变量之间的因果关系或影响路径——此时使用**结构方程模型（SEM）**
- 只有显变量之间的因果关系假设，没有潜变量——使用**路径分析（Path Analysis）**
- 医学研究中的典型场景：量表效度验证、患者报告结局（PRO）量表结构分析、心理健康模型检验

**何时不使用结构方程模型：**
- 不确定潜变量结构，想探索数据中可能隐含的因子 → 使用**探索性因子分析（EFA）**
- 仅仅做数据降维，不考虑潜变量理论 → 使用**主成分分析（PCA）**
- 只有一个因变量、简单的回归关系 → 使用**多元线性回归**或**Logistic回归**
- 两组或多组均数比较 → 使用**t检验**或**方差分析**

## 前置条件

**R包安装：**

```r
install.packages("lavaan")    # 核心包：潜变量分析
install.packages("semPlot")   # 绘制路径图
install.packages("haven")     # 读取SPSS格式数据（.sav文件）
```

**数据格式要求：**
- 每一行是一个观测个体（如一个患者），每一列是一个变量（如量表条目评分）
- 所有观测变量应为**连续型数值变量**（如Likert量表评分 1-5）
- 数据应无缺失值，或使用适当方法处理（如列删、插补）
- 如果使用最大似然估计（ML），数据量建议至少为模型中自由参数数量的 5-10 倍

**统计假设：**
- 使用最大似然估计（ML，默认方法）时，要求观测变量**多元正态分布**。若不满足，可考虑：变量转换、剔除离群值、自助重抽样（bootstrap）、或使用稳健估计方法（如MLR）
- 使用未加权最小二乘法（ULS）时，对分布无特殊要求
- 模型需为**过度识别模型**（free parameters < unique elements in covariance matrix），才能进行模型拟合优度检验

## 方法选择决策树

```
你的研究目的 →
├── 验证预设的因子结构（已知哪些条目归属哪个因子）→ 验证性因子分析（CFA）
│   └── 使用 cfa() 函数，仅写测量模型（=~ 运算符）
├── 检验潜变量之间的因果关系 → 结构方程模型（SEM）
│   └── 使用 sem() 函数，同时写测量模型（=~）和结构模型（~ 运算符）
├── 只有显变量，无潜变量，检验变量间路径关系 → 路径分析（Path Analysis）
│   └── 使用 sem() 函数，仅写回归公式（~ 运算符）
├── 不确定因子结构，想探索数据中的隐藏维度 → 探索性因子分析（EFA）
│   └── 请参考「探索性因子分析」章节
└── 有层次嵌套数据结构（如患者嵌套在医院）→ 多水平结构方程模型
    └── 请参考「多水平模型」章节
```

## 标准工作流

### 步骤1：模型设定

根据理论或研究假设，用 lavaan 公式语法定义变量之间的关系：

| 运算符 | 含义 | 示例 |
|--------|------|------|
| `=~` | 潜变量由显变量测量（测量模型） | `sos =~ PHD1+PHD2+PHD3` |
| `~` | 回归关系（结构模型） | `anx ~ sos` |
| `~~` | 方差或协方差 | `sos ~~ cog` |
| `~1` | 截距项 | `PHD1 ~ 1` |

完整的模型公式必须用单引号 `'...'` 包裹：

```r
model <- '
  # 测量模型（潜变量 =~ 显变量）
  sos =~ PHD1+PHD2+PHD3+PHD4+PHD5+PHD6+PHD7
  cog =~ PHD8+PHD9+PHD10+PHD11

  # 结构模型（潜变量 ~ 潜变量）
  cog ~ sos
'
```

三种模型设定策略：
1. **直接验证模型**：对单一假设模型进行验证（较少见）
2. **选择最优模型**：提出若干备选模型，比较拟合指标选出最优
3. **导出模型**：从初始模型出发，逐步修正直到拟合满意（最常见）

### 步骤2：模型识别

模型识别的本质是**未知参数能否由观测数据得到唯一解**：
- **识别不足模型**（under-identified）：待估参数个数多于方程个数，参数有无穷多解 → 需要增加约束或减少待估参数
- **恰好识别模型**（just-identified）：自由度为0，无法检验拟合优度 → 无实用价值
- **过度识别模型**（over-identified）：自由度>0，可以检验拟合优度 → **这是SEM所追求的目标**

经验法则：每个潜变量至少应有**3个**观测指标；自由度 `df = (p(p+1)/2) - t`，其中p为观测变量数，t为自由参数数，需 `df > 0`。

### 步骤3：模型估计

lavaan 默认使用**最大似然法（ML）**进行参数估计。其他可选方法：

| 估计方法 | 简写 | 特点 |
|----------|------|------|
| 最大似然估计 | ML | 默认方法，要求多元正态分布，大样本下有效 |
| 未加权最小二乘法 | ULS | 对分布无要求 |
| 广义最小二乘法 | GLS | 一致性有效估计 |
| 加权最小二乘法 | WLS | 不要求多元正态性（ADF方法） |
| 对角加权最小二乘法 | DWLS | ML和WLS的折中 |
| 稳健最大似然法 | MLR | 对非正态数据稳健 |

指定估计方法：`sem(model, data = df, estimator = "MLR")`

### 步骤4：模型评价

模型评价分为三个层面：

**（1）参数估计的合理性与显著性检验：**
- 方差应为正值，相关系数的绝对值 ≤ 1
- 各自由参数的P值应 < 0.05（表明该参数设为自由参数是合理的）

**（2）测量模型评价：**
- **因子载荷**（标准化）通常应 > 0.4，越大说明条目与潜变量关系越强
- R² 值表示观测变量被潜变量解释的方差比例

**（3）整体模型拟合评价（拟合指数）：**

| 拟合指数 | 中文名称 | 良好标准 | 类型 |
|----------|----------|----------|------|
| χ²/df | 卡方自由度比 | < 3（宽松）、< 2（严格） | 绝对拟合指数 |
| GFI | 拟合优度指数 | > 0.9 | 绝对拟合指数 |
| AGFI | 调整拟合优度指数 | > 0.9 | 绝对拟合指数 |
| RMSEA | 近似误差均方根 | < 0.08（可接受）、< 0.05（良好） | 绝对拟合指数 |
| SRMR | 标准化残差均方根 | < 0.08（可接受）、< 0.05（良好） | 绝对拟合指数 |
| CFI | 比较拟合指数 | > 0.9 | 相对拟合指数 |
| TLI / NNFI | Tucker-Lewis指数 | > 0.9 | 相对拟合指数 |
| AIC / BIC | 信息准则 | 越小越好（用于模型比较） | 信息标准指数 |

一个理想的拟合指数应满足：①不受样本含量影响；②惩罚复杂模型；③对误设模型敏感。

### 步骤5：模型修正

当拟合不理想时，可通过**修正指数（Modification Index, MI）**指导模型修正：

**修正原则：**
1. 先解决测量模型的设定误差，再处理结构模型
2. 每次只做一个修正，以免影响其他参数估计
3. 先增加有意义的参数，必要时再减少无意义的参数
4. 修正需有实际理论依据，不能仅凭数据驱动

**常见修正操作：**
- 添加误差项之间的协方差（如 `PHD4 ~~ PHD6`）
- 删除不显著的路径
- 增加或释放交叉载荷
- 查看 MI 值：`modificationIndices(fit) %>% arrange(-mi)`

## 代码示例

### 示例1：验证性因子分析（CFA）——量表结构效度验证

孙振球《医学统计学》第5版例23-2。Stroke-PRO量表生理领域的4个维度（躯体症状SOS、认知能力COG、言语交流VEC、自理能力SHS），20个条目。

```r
# 加载R包
library(lavaan)
library(haven)

# 读取数据（295例脑卒中患者，20个条目评分）
df23_2 <- haven::read_sav("datasets/例23-02.sav")
dim(df23_2)
## [1] 295  20

# 模型设定：仅测量模型（潜变量 =~ 显变量）
cfa_models <- ' sos =~ PHD1+PHD2+PHD3+PHD4+PHD5+PHD6+PHD7 
                cog =~ PHD8+PHD9+PHD10+PHD11
                vec =~ PHD12+PHD13+PHD14+PHD15 
                shs =~ PHD16+PHD17+PHD18+PHD19+PHD20 '

# 拟合模型
fit <- cfa(cfa_models, data = df23_2)

# 查看完整结果（含标准化系数和R²）
summary(fit, standardized = TRUE, rsquare = TRUE)
## lavaan 0.6-19 ended normally after 46 iterations
##   Estimator                                         ML
##   Number of model parameters                        46
##   Number of observations                           295
## Model Test User Model:
##   Test statistic                               630.894
##   Degrees of freedom                               164
##   P-value (Chi-square)                           0.000
## Latent Variables: (因子载荷，Std.all为完全标准化系数)
##   sos =~
##     PHD1   1.000   0.675    0.458    (标准化载荷)
##     PHD4   1.415   0.955    0.779
##     PHD5   1.433   0.966    0.831
##     ...
## Covariances: (潜变量之间的相关系数)
##   sos ~~ cog   0.372   0.592    (标准化相关系数)
##   cog ~~ vec   0.737   0.792
##   vec ~~ shs   0.653   0.634
## R-Square: (各条目被潜变量解释的变异比例)
##   PHD5   0.691
##   PHD20  0.858

# 查看拟合指数
fitMeasures(fit, fit.measures = c("chisq","df","aic","gfi","rmsea","cfi"))
##     chisq        df       aic       gfi     rmsea       cfi 
##   630.894   164.000 17132.886     0.812     0.098     0.867

# 获取标准化结果
std_res <- standardizedSolution(fit)
head(std_res)
##   lhs op  rhs est.std    se      z pvalue ci.lower ci.upper
## 1 sos =~ PHD1   0.458 0.050  9.218      0    0.361    0.556
## 2 sos =~ PHD2   0.508 0.047 10.833      0    0.416    0.600
## 3 sos =~ PHD3   0.464 0.049  9.393      0    0.367    0.561

# 绘制路径图
library(semPlot)
semPaths(fit,
         what = 'col',           # 线条用不同颜色表示
         groups = "latents",     # 根据潜变量上色
         pastel = TRUE,          # 柔和色调
         whatLabels = 'std',     # 显示标准化载荷
         style = "lisrel",       # LISREL风格
         rotation = 2,           # 旋转方向
         edge.label.cex = 1,     # 载荷字体大小
         mar = c(1, 6, 1, 6))     # 图形边距
```

### 示例2：简单结构方程模型——单外生潜变量 + 单内生潜变量

孙振球《医学统计学》第5版例23-3。躯体症状（SOS）→ 焦虑（ANX）。

```r
library(lavaan)
library(haven)

df23_3 <- haven::read_sav("datasets/例23-03.sav")
dim(df23_3)
## [1] 295  12

# 模型设定：前两个是测量模型，第三个是结构模型
sem_models <- ' sos =~ PHD1+PHD2+PHD3+PHD4+PHD5+PHD6+PHD7
                anx =~ PSD1+PSD2+PSD3+PSD4+PSD5
                anx ~ sos'         # 结构模型：SOS影响ANX

# 拟合SEM
fit <- sem(sem_models, data = df23_3)

# 查看结果
summary(fit, standardized = TRUE, rsquare = TRUE)
## lavaan 0.6-19 ended normally after 33 iterations
##   Number of model parameters                        25
##   Number of observations                           295
## Model Test User Model:
##   Test statistic                               161.889
##   Degrees of freedom                                53
##   P-value (Chi-square)                           0.000
## Latent Variables: (测量模型—因子载荷)
##   sos =~
##     PHD1   1.000   0.695   0.472
##     PHD5   1.377   0.958   0.824
##   anx =~
##     PSD1   1.000   0.931   0.721
##     PSD5   1.113   1.036   0.859
## Regressions: (结构模型—潜变量间的路径系数)
##   anx ~
##     sos    0.790   0.590   0.590   (p < 0.001)
## R-Square:
##     anx    0.348   (SOS解释了ANX 34.8%的变异)

# 查看拟合指数
fitMeasures(fit, fit.measures = c("chisq","df","aic","gfi","rmsea","cfi"))
##     chisq        df       aic       gfi     rmsea       cfi 
##   161.889    53.000 10261.913     0.908     0.083     0.932

# 绘制路径图
library(semPlot)
semPaths(fit, what = 'est', whatLabels = 'std', style = "lisrel")
```

### 示例3：多潜变量结构方程模型

孙振球《医学统计学》第5版例23-4。躯体症状(SOS) → 焦虑(ANX) → 抑郁(DEP)、回避(AVO)，同时认知功能(COG)也影响抑郁和回避。

```r
library(lavaan)
library(haven)

df23_4 <- haven::read_sav("datasets/例23-04.sav")
dim(df23_4)
## [1] 295  25

# 模型设定：5个潜变量，测量模型 + 结构模型
sem_models <- ' sos =~ PHD1+PHD2+PHD3+PHD4+PHD5+PHD6+PHD7 
                cog =~ PHD8+PHD9+PHD10+PHD11
                anx =~ PSD1+PSD2+PSD3+PSD4+PSD5
                dep =~ PSD6+PSD7+PSD8+PSD9+PSD10
                avo =~ PSD11+PSD12+PSD13+PSD14
                anx ~ sos                   # SOS → ANX
                dep ~ anx + avo + cog       # ANX, AVO, COG → DEP
                avo ~ anx + cog '           # ANX, COG → AVO

# 拟合模型
fit <- sem(sem_models, data = df23_4)

# 查看核心结果
summary(fit, standardized = TRUE, rsquare = TRUE)
## lavaan 0.6-19 ended normally after 41 iterations
##   Number of model parameters                        57
##   Number of observations                           295
## Model Test User Model:
##   Test statistic                               747.369
##   Degrees of freedom                               268
##   P-value (Chi-square)                           0.000
## Regressions: (结构模型—路径系数)
##   anx ~ sos    0.832   0.626   (p < 0.001)
##   dep ~ anx    0.342   0.369   (p < 0.001)
##   dep ~ avo    0.371   0.437   (p < 0.001)
##   dep ~ cog    0.224   0.246   (p < 0.001)
##   avo ~ anx    0.609   0.557   (p < 0.001)
##   avo ~ cog    0.356   0.331   (p < 0.001)
## Covariances:
##   sos ~~ cog   0.401   0.622   (SOS与COG的相关)
## R-Square:
##     anx        0.392
##     dep        0.798
##     avo        0.564

# 绘制路径图
library(semPlot)
semPaths(fit, 'std', 'std', style = "lisrel", exoCov = FALSE)
```

### 模型修正辅助代码

```r
# 查看修正指数，按MI降序排列
library(dplyr)
modificationIndices(fit) %>% arrange(-mi) %>% head(10)

# 当模型拟合欠佳时，根据MI添加参数（需理论支撑）
# 例如：添加误差项间的协方差
model_revised <- ' sos =~ PHD1+PHD2+PHD3+PHD4+PHD5+PHD6+PHD7 
                   cog =~ PHD8+PHD9+PHD10+PHD11
                   vec =~ PHD12+PHD13+PHD14+PHD15 
                   shs =~ PHD16+PHD17+PHD18+PHD19+PHD20
                   PHD4 ~~ PHD6 '    # 根据MI添加误差相关

fit_revised <- cfa(model_revised, data = df23_2)

# 比较修正前后的模型
anova(fit, fit_revised)
```

## 结果解读指南

### lavaan输出结构

`summary(fit, standardized = TRUE, rsquare = TRUE)` 输出四个部分：

**第1部分：模型基本信息（前9行）**
- `ended normally after XX iterations`：模型收敛正常 → 可信任结果
- `Estimator` ML：使用最大似然法估计
- `Number of model parameters`：自由参数个数
- `Number of observations`：实际使用的样本量
- `Test statistic`（卡方值）越大，模型与数据拟合越差（但易受样本量影响，大样本时几乎总是显著）
- `P-value (Chi-square)`：理想情况下希望 p > 0.05（不显著），但在大样本时几乎总是显著，此时需结合其他拟合指数判断

**第2部分：参数估计（Parameter Estimates）**
- `Latent Variables`：测量模型结果，即**因子载荷**
  - `Estimate`：非标准化载荷系数（类似回归系数）
  - `Std.Err`：标准误
  - `z-value`：Z检验统计量
  - `P(>|z|)`：p值，**p < 0.05 表示该载荷具有统计学意义**
  - `Std.lv`：仅潜变量标准化的载荷系数
  - `Std.all`：**完全标准化载荷系数**（潜变量和观测变量都标准化），这是论文中最重要的报告指标
  - 经验标准：标准化载荷 > 0.4 可接受，> 0.7 理想
- `Regressions`：结构模型结果，即**路径系数**
  - 解读同因子载荷，表示潜变量之间影响的大小和方向
- `Covariances`：潜变量之间的协方差/相关系数
- `Variances`：误差方差（`.变量名`表示显变量的残差方差）

**第3部分：R-Square**
- 显变量的R²：潜变量对该观测指标的解释程度（如 PHD5 的 0.691 = 躯体症状解释了该条目69.1%的变异）
- 内生潜变量的R²：结构模型中外生变量对内生变量的解释程度（如 anx 的 0.348 = SOS 解释了 ANX 34.8%的变异）

### 论文中如何报告

国内医学论文中CFA结果的典型报告格式：

> 验证性因子分析结果显示，所有条目在对应因子的标准化载荷均大于0.4（范围：0.46–0.93），差异均有统计学意义（P < 0.001）。模型拟合指标中，χ²/df = 3.85，RMSEA = 0.098，CFI = 0.867，GFI = 0.812，模型拟合效果尚可接受。躯体症状(SOS)与认知能力(COG)呈中等正相关（r = 0.592），与言语交流(VEC)呈强正相关（r = 0.702）。

SEM结果的典型报告格式：

> 结构方程模型分析结果显示，躯体症状对焦虑具有正向预测作用（β = 0.590，P < 0.001），焦虑对回避具有正向预测作用（β = 0.557，P < 0.001）。模型拟合指标为 χ²/df = 3.05、RMSEA = 0.083、CFI = 0.932，表明模型拟合良好。

## 常见问题与注意事项

**Q1: CFA的卡方检验P值总是<0.05，是不是模型一定不好？**

不是。卡方值对样本量极其敏感。样本量较大时（如n > 200），即使模型与数据差异很小，卡方检验也会显著。此时应**重点参考CFI、TLI、RMSEA、SRMR等不受样本量影响的拟合指数**，而非依赖卡方P值。

**Q2: lavaan中 cfa() 和 sem() 有什么区别？**

本质上 `cfa()` 是 `sem()` 的包装函数。`cfa()` 默认所有潜变量之间两两相关（自动添加协方差），适合纯测量模型。`sem()` 需要手动指定所有关系。但实际上两者可以互用。

**Q3: 标准化载荷出现负值或>1怎么办？**

标准化载荷的绝对值>1（如1.2）通常意味着模型存在问题：①出现了负误差方差（Heywood case）；②模型设定错误；③样本量过小导致估计不稳定。应检查模型设定和样本量。

**Q4: 何时使用 lavaan 的 ~ 运算符 vs =~ 运算符？**

- `=~` 用于**测量模型**：潜变量（左侧）由显变量（右侧）测量，如 `anx =~ PSD1+PSD2+PSD3`
- `~` 用于**结构模型**：潜变量（左侧）被其他潜变量或显变量（右侧）影响，如 `anx ~ sos`

**Q5: SEM对样本量有什么要求？**

一般建议：
- 绝对底线：n ≥ 100
- 经验标准：n ≥ 200，或观测变量数的10-15倍
- 自由参数数的5-10倍
- 若使用ML估计但数据非正态，需要更大样本量

**Q6: 修改指数（MI）值多大才需要考虑修正？**

没有绝对标准。通常 MI > 10 或 MI > 3.84（对应于卡方检验临界值）时值得关注。但核心原则是：**修改必须有理论依据，不能纯粹数据驱动**。没有理论支撑的修正可能导致模型过度拟合当前样本，泛化到新样本时拟合很差。

**Q7: 探索性因子分析(EFA)和验证性因子分析(CFA)的根本区别？**

| 维度 | EFA | CFA |
|------|-----|-----|
| 目的 | 探索数据结构，发现潜在因子 | 验证预设的因子结构 |
| 因子数 | 数据驱动 | 理论驱动 |
| 载荷 | 每个条目可在所有因子上有载荷 | 每个条目仅在预设因子上有载荷（其余固定为0） |
| 因子相关 | 通常假设因子正交或斜交 | 允许且预期因子相关 |
| 使用包 | `psych::fa()`, `factanal()` | `lavaan::cfa()` |

**Q8: SPSS中的AMOS和R中的lavaan有何区别？**

- AMOS：图形界面操作，画路径图即为模型设定，适合不熟悉编程的用户，是医学研究常用工具
- lavaan：代码驱动，通过公式语法设定模型，更灵活，完全免费开源，结果可重复性更强
- 两者的统计原理和内置估计方法完全相同，均为协方差结构分析
- 论文中两种工具均可使用，lavaan的结果与AMOS等价

**Q9: 模型不收敛怎么办？**

常见原因与对策：
1. 模型设定有误（变量名写错、关系设定不合理）→ 检查公式语法
2. 数据问题（高度共线性、变量方差过小、缺失值过多）→ 检查数据描述统计
3. 模型过于复杂，样本量不足 → 简化模型或收集更多数据
4. 起始值不合理 → 尝试 `start` 参数提供合理起始值

**Q10: lavaan 公式中注释用什么符号？**

使用 `#` 符号，但必须**单独一行**，不能放在公式行末尾。注释前后不能有空白。

```r
# 正确写法
model <- '
  # 这是测量模型
  sos =~ PHD1+PHD2+PHD3
  # 这是结构模型
  anx ~ sos
'

# 错误写法（注释在公式同一行）
model <- '
  sos =~ PHD1+PHD2+PHD3  # 这样不行
'
```
