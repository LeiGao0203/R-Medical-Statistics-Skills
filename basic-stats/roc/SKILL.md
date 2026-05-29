---
name: medical-stat-roc
description: "R语言医学统计：ROC曲线分析与诊断试验评价。涵盖ROC曲线绘制、AUC计算、AUC比较（DeLong检验）、最佳截断值确定（Youden指数）、灵敏度/特异度/阳性预测值/阴性预测值计算。TRIGGER when user mentions ROC曲线、诊断试验、AUC、灵敏度、特异度、最佳截断值、pROC，or asks about evaluating diagnostic test performance. SKIP for 生存分析、校准曲线、决策曲线分析。"
---

# ROC曲线 (ROC Curve Analysis)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用场景：**
- 评价单一诊断试验/生物标志物的区分能力（AUC）
- 比较两种或多种诊断方法的诊断效能（DeLong检验）
- 确定连续型生物标志物的最佳截断值（Youden指数）
- 报告诊断试验的灵敏度、特异度、阳性预测值、阴性预测值
- 临床预测模型中模型区分度的评价

**不适用场景：**
- 生存数据的ROC分析 → 使用 survivalROC / timeROC 包（生存分析）
- 模型校准评估 → 使用校准曲线（calibration plot）
- 临床净获益评估 → 使用决策曲线分析（decision curve analysis）
- 多分类诊断评价 → 使用 pairwise / multiclass ROC 方法

## 前置条件

**R包安装：**

```r
install.packages("pROC")      # ROC分析核心包
install.packages("ggplot2")   # 自定义绘图
install.packages("purrr")     # 函数式编程辅助
```

**数据格式要求：**
- 一列真实分类（金标准结果），必须为二分类 factor，第一水平为阴性/对照
- 一列连续型预测变量（生物标志物值、概率预测值等）
- 样本量建议：病例组和对照组各 ≥ 30 例

**统计前提：**
- 必须有明确的金标准作为分类参照
- 病例组和对照组的定义必须基于金标准
- 区间或比率型预测变量有意义（连续型变量）
- 对于AUC比较，假设两ROC曲线来自相同受试者或独立样本（paired vs unpaired）

## 方法选择决策树

```
你的数据情况 →
├── 单一生物标志物/单一诊断试验 → pROC::roc() 绘制ROC曲线，计算AUC
├── 比较两种诊断方法的AUC（同一样本） → pROC::roc.test(method="delong")
├── 比较两种诊断方法的AUC（独立样本） → pROC::roc.test(method="bootstrap")
├── 需要AUC置信区间 → pROC::ci.auc()
├── 寻找最佳截断值 → 使约登指数(Youden)最大的截断值
│   └── print.thres="best" 或 coords(res, "best")
├── 需要完整诊断试验评价指标 → coords() 可一次性输出所有指标
└── 需要平滑ROC曲线 → roc(..., smooth=TRUE)
```

## 标准工作流

### 步骤1：数据准备
- 确保数据框包含分类结果变量（factor，第一水平为阴性）和连续型预测变量
- 确认金标准分类的组别标签

### 步骤2：探索性分析
- 使用 `str()` 检查数据类型和因子水平
- 使用分组箱线图初步比较两组生物标志物的分布差异

### 步骤3：绘制ROC曲线并计算AUC
- 使用 `pROC::roc()` 构建ROC对象
- 使用 `plot.roc()` 绘制曲线，`print.auc=TRUE` 标注AUC值
- 使用 `ci.auc()` 获取AUC的95%置信区间

### 步骤4：确定最佳截断值
- 使用 `coords(res, "best", ret="threshold")` 获取最佳截断值（约登指数最大）
- 同时输出该截断值下的灵敏度、特异度、阳性预测值、阴性预测值
- `coords(res, "best", ret=c("threshold","sensitivity","specificity","ppv","npv"))`

### 步骤5：对比两种诊断方法（如适用）
- 使用 `roc.test()` 进行DeLong检验比较两个AUC
- 报告 AUC差异及其95%置信区间和p值

### 步骤6：结果报告
论文中建议报告的统计量：
- AUC（95% CI）
- 最佳截断值
- 对应的灵敏度、特异度、阳性预测值、阴性预测值
- 如有多方法比较，报告DeLong检验 p 值和AUC差值

## 代码示例

```r
# ---- 加载包 ----
library(pROC)
library(ggplot2)

# ---- 示例1：基础ROC分析 (aSAH数据集) ----
data(aSAH, package = "pROC")
str(aSAH)
## 'data.frame':    113 obs. of  7 variables:
##  $ outcome: Factor w/ 2 levels "Good","Poor": 1 1 1 1 2 2 1 2 1 1 ...
##  $ s100b  : num  0.13 0.14 0.1 0.04 0.13 0.1 0.47 0.16 0.18 0.1 ...
##  $ ndka   : num  3.01 8.54 8.09 10.42 17.4 ...

# 构建ROC对象
res <- roc(outcome ~ s100b, data = aSAH)
## Setting levels: control = Good, case = Poor
## Setting direction: controls < cases

# 查看AUC
res
## 
## Call:
## roc.formula(formula = outcome ~ s100b, data = aSAH)
## 
## Data: s100b in 72 controls (outcome Good) < 41 cases (outcome Poor).
## Area under the curve: 0.7314

# ---- 示例2：绘制ROC曲线并标注AUC ----
plot(res,
     auc.polygon = TRUE,          # 填充AUC面积
     auc.polygon.col = "steelblue",
     print.auc = TRUE,            # 显示AUC值
     print.auc.x = 0.95,
     print.auc.y = 0.9,
     print.auc.col = "firebrick",
     print.auc.cex = 2)

# ---- 示例3：平滑ROC曲线 ----
res_smooth <- roc(outcome ~ s100b, data = aSAH, smooth = TRUE)
plot(res_smooth,
     auc.polygon = TRUE,
     auc.polygon.col = "steelblue",
     print.auc = TRUE,
     print.auc.x = 0.95,
     print.auc.y = 0.9,
     print.auc.col = "firebrick",
     print.auc.cex = 2)

# ---- 示例4：寻找最佳截断值 ----
plot(res,
     auc.polygon = TRUE,
     auc.polygon.col = "steelblue",
     legacy.axes = TRUE,          # X轴从0开始（1-特异度）
     print.thres = "best")        # 标注最佳截断值（约登指数最大）
## 最佳截断值为 0.205，特异度 0.806，灵敏度 0.634

# ---- 示例5：AUC 置信区间 ----
ci.auc(res)
## 95% CI: 0.6301-0.8326 (DeLong)

# 完整诊断评价指标（在最佳截断值处）
coords(res, "best", ret = c("threshold", "sensitivity", "specificity",
                             "ppv", "npv", "tp", "tn", "fp", "fn"))
## threshold: 0.205
## sensitivity specificiy       ppv       npv
##   0.6341463   0.8055556   0.6500000   0.7945205

# ---- 示例6：比较两个AUC（DeLong检验） ----
roc_s100b <- roc(outcome ~ s100b, data = aSAH)
roc_ndka  <- roc(outcome ~ ndka,  data = aSAH)

# DeLong检验
roc.test(roc_s100b, roc_ndka, method = "delong")
## 
## 	DeLong's test for two correlated ROC curves
## 
## data:  roc_s100b and roc_ndka
## Z = -1.4661, p-value = 0.1426
## alternative hypothesis: true difference in AUC is not equal to 0
## 95 percent confidence interval:
##  -0.20152316  0.02944437
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.7313686   0.8174080

# ---- 示例7：自定义 ggplot 绘制ROC ----
# 先获取坐标数据
roc_list <- list(
  s100b = roc(outcome ~ s100b, data = aSAH),
  ndka  = roc(outcome ~ ndka, data = aSAH)
)

# 提取坐标并合并
roc_data <- do.call(rbind, lapply(names(roc_list), function(name) {
  r <- roc_list[[name]]
  data.frame(
    specificity = r$specificities,
    sensitivity = r$sensitivities,
    marker = name
  )
}))

ggplot(roc_data, aes(1 - specificity, sensitivity, color = marker)) +
  geom_path(size = 1.1) +
  geom_abline(linetype = "dashed", color = "gray50") +
  coord_fixed() +
  theme_bw() +
  labs(x = "1 - 特异度 (1 - Specificity)",
       y = "灵敏度 (Sensitivity)",
       color = "生物标志物")
```

## 结果解读指南

| 指标 | 公式 | 含义 |
|------|------|------|
| **灵敏度（Sensitivity）** | a / (a + c) | 患者被正确诊断为阳性的概率（真阳性率） |
| **特异度（Specificity）** | d / (b + d) | 非患者被正确诊断为阴性的概率（真阴性率） |
| **误诊率（假阳性率）** | b / (b + d) | 非患者被错误诊断为阳性的概率 = 1 − 特异度 |
| **漏诊率（假阴性率）** | c / (a + c) | 患者被错误诊断为阴性的概率 = 1 − 灵敏度 |
| **阳性预测值（PPV）** | a / (a + b) | 诊断阳性者中真正患病的概率 |
| **阴性预测值（NPV）** | d / (c + d) | 诊断阴性者中真正未患病的概率 |
| **阳性似然比（LR+）** | Se / (1 − Sp) | 真阳性率与假阳性率之比，>10 强证据 |
| **阴性似然比（LR−）** | (1 − Se) / Sp | 假阴性率与真阴性率之比，<0.1 强证据 |
| **正确率（Accuracy）** | (a + d) / N | 总符合率，观察结果与金标准的一致程度 |
| **约登指数（Youden Index）** | Se + Sp − 1 | 综合诊断真实性指标，最大值处为最佳截断值 |
| **比数积（OP）** | ad / bc | 比值比形式，越大诊断价值越高 |

**AUC解读：**
- AUC = 0.5：诊断能力无异于随机猜测
- 0.5 < AUC < 0.7：诊断能力较低
- 0.7 ≤ AUC < 0.8：诊断能力可接受
- 0.8 ≤ AUC < 0.9：诊断能力优良
- AUC ≥ 0.9：诊断能力极佳
- AUC = 1.0：完美分类（现实中几乎不存在）

**DeLong检验解读：**
- p < 0.05 表示两个AUC之间的差异具有统计学意义
- 差值 95% CI 不包含0同样说明差异显著
- 样本量较小时 DeLong 检验可能保守，可考虑 bootstrap 方法

## 常见问题与注意事项

**问：pROC 的 roc() 如何判断正负方向？**
pROC 默认将 outcome factor 的第一水平作为阴性（control），第二水平作为阳性（case）。函数会自动检测方向并提示 "Setting levels" 和 "Setting direction"，注意检查这两个提示是否符合预期。若不符合，可以使用 `levels` 参数重新指定或手动调换 factor 的水平顺序。

**问：两个诊断试验的 AUC 比较用什么方法？**
同一样本配对设计（paired）：使用 `roc.test(roc1, roc2, method="delong")`。
两独立样本（unpaired）：使用 `roc.test(roc1, roc2, method="bootstrap")`。
默认 DeLong 检验适用于相关ROC曲线比较，是临床最常用的方法。

**问：约登指数、最短距离、最接近(0,1)等方法哪个更好？**
约登指数（Youden Index）最常用，在灵敏度和特异度之间对等权重权衡。默认使用 `coords(res, "best")` 即为约登指数最大处的截断值。如需其他标准，传入 `best.method="closest.topleft"` 等参数。

**问：阳性预测值和阴性预测值与患病率有关吗？**
是的。灵敏度和特异度是诊断试验的固有属性，不随患病率变化。但PPV和NPV强烈依赖于患病率——患病率越高，PPV越高、NPV越低。因此在低患病率人群中筛查时，即使灵敏度很高，PPV也可能很低（假阳性多）。

**问：AUC ≥ 0.8 就说明诊断效果好？**
AUC衡量的是整体区分能力，但一个 AUC=0.8 的模型在某个截断值下可能灵敏度只有0.5。AUC高不等于临床可用的灵敏度和特异度——必须结合最佳截断值处的具体灵敏度和特异度综合判断。

**问：为什么手动计算ROC和pROC的x轴方向有时不同？**
手动计算常以"1−特异度"为x轴（从0到1），pROC 默认以"特异度"为x轴（从1到0）。设置 `legacy.axes = TRUE` 可使x轴从0开始（即1−特异度方向）。两种画法本质相同，无需纠结。

**问：SPSS 能否做 DeLong 检验？**
SPSS 基础模块不直接支持 DeLong 检验。通常需使用 MedCalc 或 R 的 pROC 包完成。DeLong 检验是学术界比较AUC的标准方法。

**问：ROC曲线的 `direction` 参数什么时候需要手动设置？**
当预测变量值越高代表患病风险越低时（如某些负向生物标志物），需要使用 `direction=">"` 或在 roc 中显式指定 `direction`。pROC 会输出提示信息，注意检查即可。
