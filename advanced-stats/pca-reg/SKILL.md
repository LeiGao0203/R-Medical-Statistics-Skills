---
name: medical-stat-pca-reg
description: "R语言医学统计：主成分回归。将主成分分析提取的主成分作为自变量进行回归分析，解决自变量共线性问题。TRIGGER when user mentions 主成分回归、PCR、共线性处理、主成分替代变量，or asks about handling multicollinearity with PCA. SKIP for 普通多元回归、岭回归、LASSO回归。"
---

# 主成分回归 (Principal Component Regression, PCR)
> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

- 当多元线性回归中自变量之间存在高度共线性（VIF > 10 或条件数过大），普通最小二乘估计不稳定时
- 需要降维后保留原始变量的大部分信息再进行回归建模
- 医学研究中自变量维度较高（如多个生化指标、基因表达数据）而样本量有限时
- 不适用于：自变量之间相关性较弱的情况——此时普通多元回归效率更高；需要变量选择和系数解释的场合——优先考虑 LASSO 回归；既要处理共线性又要控制系数收缩——使用岭回归

## 前置条件

```r
install.packages("pls")
install.packages("tidymodels")
```

```r
library(pls)
suppressMessages(library(tidymodels))
```

- 数据格式：数值型数据框，因变量和自变量均为连续变量
- 所有自变量需标准化（scale = TRUE），消除量纲影响
- 主成分个数 k 应通过交叉验证确定，不可随意指定
- 样本量 n 应大于自变量个数 p，否则主成分个数上限为 n-1

## 方法选择决策树

```
你的需求 →
├── 快速实现，代码简洁 → 使用 pls::pcr()
├── 需要统一的机器学习流程（预处理+调参+预测） → 使用 tidymodels 工作流
│   ├── 有分类变量需虚拟编码 → recipe 中添加 step_dummy()
│   ├── 只需数值标准化 → recipe 中添加 step_normalize()
│   └── 需后续替换为其他模型（如 LASSO） → 替换 linear_reg() 引擎即可
└── 需要偏最小二乘回归（同时利用 X 和 Y 信息降维） → 使用 pls::plsr()
```

## 标准工作流

### 步骤1：数据准备与探索

```r
# 查看数据结构和共线性
str(mtcars)
pairs(mtcars[, c("mpg", "disp", "drat", "wt", "qsec")])
cor(mtcars[, c("mpg", "disp", "drat", "wt", "qsec")])
```

确认自变量间存在较强相关性（如 mpg 与 wt 相关系数 > 0.8），PCR 才比普通回归有优势。

### 步骤2：前提条件检验

- 检查缺失值：`anyNA(mtcars)`
- 确认所有变量为数值型：`sapply(mtcars, is.numeric)`
- 决定是否对因变量做变换（如对数变换以改善残差正态性）

### 步骤3：执行主成分回归

两种方法选其一。方法一（pls 包）适合快速分析，方法二（tidymodels）适合生产级流程。

若用 pls 包，`pcr()` 中 `validation = "CV"` 启用交叉验证，`scale = TRUE` 标准化变量。通过 `summary()` 查看各主成分数对应的 RMSE 和方差解释比例，选择 RMSE 最小的主成分数。

若用 tidymodels，通过 `tune_grid()` 调优 `num_comp` 超参数，`show_best()` 确定最佳主成分个数，`finalize_workflow()` 锁定参数后用全量数据重新拟合。

### 步骤4：结果解读

核心指标：
- **CV RMSE**：交叉验证均方根误差，越小越好，截距项（0 个主成分）的 RMSE 作为基线
- **方差解释比例（X）**：主成分对自变量总方差的累计解释率，通常 80% 以上即可
- **方差解释比例（Y）**：主成分对因变量方差的累计解释率，反映模型的预测能力
- **validationplot()**：可视化 RMSE / R² 随主成分数的变化，选择拐点处的主成分数

### 步骤5：结果报告

论文中可描述为：「采用主成分回归处理自变量共线性问题，通过 10 折交叉验证确定最佳主成分个数为 2，模型交叉验证 RMSE 为 29.1，前 2 个主成分解释了自变量总变异的 89.4%，以及因变量 hp 变异的 81.3%。」

## 代码示例

### 方法一：pls 包

```r
rm(list = ls())
library(pls)

# 查看数据
head(mtcars)
##                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
## Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
## Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
## Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1

set.seed(1)
model <- pcr(hp ~ mpg + disp + drat + wt + qsec,
             data = mtcars, scale = TRUE, validation = "CV")
summary(model)
## Data:    X dimension: 32 5
##  Y dimension: 32 1
## Fit method: svdpc
## Number of components considered: 5
##
## VALIDATION: RMSEP
## Cross-validated using 10 random segments.
##        (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
## CV           69.66    43.74    34.58    34.93    36.34    37.40
## adjCV        69.66    43.65    34.30    34.61    35.95    36.95
##
## TRAINING: % variance explained
##     1 comps  2 comps  3 comps  4 comps  5 comps
## X     69.83    89.35    95.88    98.96   100.00
## hp    62.38    81.31    81.96    81.98    82.03

# 可视化
validationplot(model)               # RMSE
validationplot(model, val.type = "MSEP")  # MSE
validationplot(model, val.type = "R2")    # R²

# 预测新数据（选择2个主成分）
test <- head(mtcars)
predict(model, test, ncomp = 2)
## , , 2 comps
##                         hp
## Mazda RX4         155.2385
## Mazda RX4 Wag     146.6904
## Datsun 710        100.4458
```

### 方法二：tidymodels

```r
suppressMessages(library(tidymodels))
tidymodels_prefer()

# 模型设定
set.seed(994)
lm_spec <- linear_reg() %>% set_engine("lm")

# 数据划分
mtcars_resamples <- vfold_cv(mtcars, v = 10)

# 配方（预处理步骤）
mtcars_pca_recipe <- recipe(hp ~ mpg + disp + drat + wt + qsec,
                            data = mtcars) %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), num_comp = tune())

# 工作流
mtcars_pca_workflow <- workflow() %>%
  add_model(lm_spec) %>%
  add_recipe(mtcars_pca_recipe)

# 超参数调优
num_comp_grid <- grid_regular(num_comp(range = c(0, 5)), levels = 6)
mtcars_pca_tune <- tune_grid(mtcars_pca_workflow,
                             resamples = mtcars_resamples,
                             grid = num_comp_grid)

# 查看调优结果
autoplot(mtcars_pca_tune)
show_best(mtcars_pca_tune)
## # A tibble: 5 x 7
##   num_comp .metric .estimator  mean     n std_err .config
##      <int> <chr>   <chr>      <dbl> <int>   <dbl> <chr>
## 1        2 rmse    standard    29.1    10    5.63 Preprocessor3_Model1
## 2        3 rmse    standard    31.0    10    5.11 Preprocessor4_Model1

# 使用最佳参数重新建模并预测
mtcars_pca_workflow_final <-
  finalize_workflow(mtcars_pca_workflow,
                    select_best(mtcars_pca_tune, metric = "rmse"))
mtcars_pca_fit_final <- fit(mtcars_pca_workflow_final, data = mtcars)
predict(mtcars_pca_fit_final, new_data = head(mtcars))
## # A tibble: 6 x 1
##   .pred
##   <dbl>
## 1  155.
## 2  147.
## 3  100.
```

## 结果解读指南

| 输出项 | 含义 | 判断标准 |
|--------|------|----------|
| `CV` 行各列 | 不同主成分数对应的交叉验证 RMSE | 取最小值对应的主成分数（本例为 2） |
| `adjCV` | 偏差校正后的 CV RMSE | 与 CV 一致时结果更可靠 |
| `X` 行方差解释 | 主成分对自变量方差的累计解释率 | 通常 > 80% 即可，本例 2 个主成分已达 89.35% |
| `hp` 行方差解释 | 主成分对因变量方差的累计解释率 | 体现模型预测力，本例 2 个主成分为 81.31% |
| `show_best()` 的 `mean` | 10 折交叉验证平均 RMSE | 2 个主成分 RMSE = 29.1，比截距项（0 个主成分）的 32.1 有明显改善 |

- 当增加主成分后 RMSE 不降反升（如本例从 2 个到 3 个：34.58 → 34.93），说明额外主成分引入了噪声
- 主成分是原始变量的线性组合，直接解释系数没有实际意义；如需解释单个变量的效应，应回到原始变量做回归或使用 LASSO
- `validationplot()` 中 RMSE 最低点即最优主成分数

## 常见问题与注意事项

**Q: 主成分回归和逐步回归有什么区别？**
A: 主成分回归通过降维处理共线性，不删除变量；逐步回归通过 AIC/p 值筛选变量，会丢弃变量。两者目的不同：PCR 解决共线性，逐步回归做变量选择。

**Q: pls 包和 tidymodels 结果是否一致？**
A: 两种方法的交叉验证划分策略不同（pls 默认 10 段随机分割；tidymodels 默认 10 折），导致 RMSE 数值略有差异，但最优主成分数的结论通常一致。

**Q: 主成分个数如何选择？**
A: （1）交叉验证 RMSE 最小值；（2）累计方差解释率 > 80%；（3）碎石图拐点。综合判断，不要仅凭单一标准。

**Q: 主成分回归的系数如何还原为原始变量的系数？**
A: 主成分是原始变量的线性组合，可通过载荷矩阵将主成分系数反推为原始变量系数，但此时系数已非最小二乘无偏估计，标准误也需重新计算。建议直接在原始变量尺度上使用岭回归。

**Q: 主成分回归 vs 偏最小二乘回归（PLSR）如何选？**
A: PCR 只利用 X 的方差结构降维，选择使 X 方差最大的方向；PLSR 同时利用 X 和 Y 的信息降维，选择使 X 和 Y 协方差最大的方向。当 Y 和 X 的前几个主成分关联较弱时，PLSR 通常优于 PCR。在 R 中 `pls::plsr()` 即可实现 PLSR。

- 主成分分析前必须标准化，否则方差大的变量会主导主成分方向
- 样本量较小时交叉验证的折数不宜过多（5 折或留一法）
- tidymodels 中 `num_comp = tune()` 时需注意 `step_pca()` 默认保留所有主成分，必须通过 `tune()` 配合调优网格限制
- 主成分回归本质上仍是线性模型，需满足线性回归的基本假设（线性、独立、正态、等方差）
