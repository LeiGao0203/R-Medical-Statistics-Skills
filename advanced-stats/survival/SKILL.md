---
name: medical-stat-survival
description: "R语言医学统计：生存分析。涵盖Kaplan-Meier估计、log-rank检验、Cox比例风险回归模型（单因素与多因素）、PH假设检验、HR值解读。TRIGGER when user mentions 生存分析、Kaplan-Meier、KM曲线、log-rank、Cox回归、Cox比例风险、HR、hazard ratio，or asks about time-to-event data analysis. SKIP for Fine-Gray竞争风险模型、生存曲线可视化、参数化生存模型。"
---

# 生存分析 (Survival Analysis)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**使用生存分析的典型场景：**

- 临床研究中需要分析患者从入组到发生终点事件（死亡、复发、进展等）的时间，且存在删失数据
- 比较不同治疗组之间生存时间（生存曲线）的差异
- 探索影响生存时间的独立危险因素或保护因素
- 肿瘤学、心血管疾病、慢性病流行病学中任何涉及 time-to-event 终点的研究

**不使用标准生存分析的情况：**

| 你的情况 | 应使用的方法 |
|----------|-------------|
| 存在竞争风险（如非肿瘤死亡的竞争事件） | Fine-Gray 竞争风险模型 |
| 仅需可视化生存曲线而非统计分析 | 生存曲线可视化（survminer 绘图） |
| 假定了特定的生存时间分布（如 Weibull、指数分布） | 参数化生存模型（`survreg()`） |
| 纵向数据中指标多次重复测量 | 联合模型（joint model）或多状态模型 |
| 仅比较一个时间点的率（如 5 年生存率） | 单率比较或卡方检验 |

**医学研究常见应用：**

- 比较两种化疗方案的总体生存期（OS）和无进展生存期（PFS）
- 多因素 Cox 回归筛选影响预后的独立危险因素
- 寻找生物标志物的最佳 cut-point 以区分高/低风险组

## 前置条件

**R 包安装：**

```r
install.packages(c("survival", "survminer", "broom"))
```

**核心数据要求：**

- `time`：生存时间（连续型数值，≥0），单位可为天、月、年，需统一
- `status`（或 `event`）：结局事件指示变量，通常 1 = 发生终点事件，0 = 删失
  - 注意 `survival::lung$status` 原始编码为 1 = 删失，2 = 死亡，使用时需先转换
- 分组变量（如需比较）：因子型（factor），多分类需正确设置水平
- 协变量（Cox 回归）：可为连续型或分类型，分类变量建议转为 factor 以避免被当作连续变量

**科学假设：**

- **随机删失假设**：删失的发生独立于终点事件的发生（即删失不提供预后信息）
- **非信息性删失**：删失的发生不应由未来的终点事件决定
- **比例风险假设（PH assumption）**（仅 Cox 回归）：协变量的风险比（HR）在整个随访期间恒定不变
  - 可通过 Schoenfeld 残差检验（`cox.zph()`）或 K-M 曲线目测交叉来判断是否满足

## 方法选择决策树

```
你的研究目标 →
├── 仅描述单组患者的生存过程（估计各时间点生存率）
│   ├── 数据量大且时间记录精确 → K-M 法（`survfit(Surv(time,status) ~ 1)`）
│   └── 按固定时间区间汇总（如按月、按年） → 寿命表法（`surv_summary()`）
│
├── 比较两组或多组生存曲线
│   ├── 无协变量，仅比较组间生存过程 → log-rank 检验（`survdiff()`）或 Breslow 检验
│   ├── 需要呈现生存曲线 + risk table + p 值 → `ggsurvplot()`（survminer 包）
│   └── 需要寻找连续变量最佳 cut-point 再比较分组 → `surv_cutpoint()` + `surv_categorize()`
│
├── 探索影响生存时间的因素，校正混杂
│   ├── 单个自变量 → 单因素 Cox 回归（`coxph(Surv(time,status) ~ x)`）
│   ├── 多个自变量同时纳入 → 多因素 Cox 回归（`coxph(Surv(time,status) ~ x1 + x2 + x3)`）
│   └── 展示多因素结果 → 森林图（`ggforest()`）
│
├── PH 假设不满足时
│   ├── 系数分段恒定 → 对时间分层（`survSplit()` 按拐点分段 + `strata(tgroup)`）
│   ├── 系数连续变化且可变换为线性 → 时依系数变换（`tt()` 函数，如 `x * log(t+20)`）
│   └── 协变量本身随时间变化（如药物剂量变化） → 时依协变量（需要对数据 reshape）
│
└── 需要提取整洁结果用于制表或报告
    └── `broom::tidy(fit.cox, exponentiate = TRUE, conf.int = TRUE)`
```

## 标准工作流

### 步骤 1：数据准备与探索

```r
library(survival)
library(survminer)

# 使用 lung 数据集（228 例肺癌患者）
df <- lung
# 将 status 重编码：1=死亡, 0=删失（习惯用法）
df$status <- ifelse(df$status == 2, 1, 0)
```

**数据概览**：`time`（天），`status`（1 = 死亡，0 = 删失），`sex`（1 = 男，2 = 女），`age`（岁），`ph.ecog`（ECOG 评分 0-3），`ph.karno`（Karnofsky 评分）。

### 步骤 2：构建生存对象

```r
# 用 Surv() 创建生存对象
surv_obj <- Surv(time = df$time, event = df$status)
head(surv_obj)
## 带 "+" 号者为删失数据
```

### 步骤 3：生存过程的描述

```r
# 单组 K-M 生存曲线
fit_km <- survfit(Surv(time, status) ~ 1, data = df)

# 绘制生存曲线（含 95% CI 和中位生存时间）
ggsurvplot(fit_km, conf.int = TRUE, palette = "blue",
           surv.median.line = "hv", ggtheme = theme_bw())
```

### 步骤 4：生存过程的比较

```r
# log-rank 检验（比较 sex 分组）
fit_logrank <- survdiff(Surv(time, status) ~ sex, data = df)
fit_logrank

# 绘制分组 KM 曲线（含 p 值和 risk table）
fit_km_grp <- survfit(Surv(time, status) ~ sex, data = df)
ggsurvplot(fit_km_grp, data = df, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, surv.median.line = "hv",
           legend.title = "Sex", legend.labs = c("Male", "Female"),
           palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw())
```

### 步骤 5：一般 Cox 回归分析

```r
# 分类变量转为因子型（coxph 会自动哑变量编码）
lung$sex <- factor(lung$sex, labels = c("female", "male"))
lung$ph.ecog <- factor(lung$ph.ecog,
  labels = c("asymptomatic", "symptomatic", "in bed <50%", "in bed >50%"))

# 多因素 Cox 回归
fit_cox <- coxph(Surv(time, status) ~ sex + age + ph.karno, data = lung)
summary(fit_cox)
```

### 步骤 6：比例风险（PH）假设检验

```r
# Schoenfeld 残差检验
fit_ph <- cox.zph(fit_cox)
fit_ph

# 图形化查看
ggcoxzph(fit_ph)
# 或查看特定变量的残差
plot(fit_ph[3]); abline(0, 0, col = "red")
```

### 步骤 7：PH 不满足时的处理（时间分层）

```r
# 以 veteran 数据集为例
vet2 <- survSplit(Surv(time, status) ~ ., data = veteran,
                  cut = c(90, 180), episode = "tgroup", id = "id")

fit2 <- coxph(Surv(tstart, time, status) ~
                trt + prior + karno:strata(tgroup), data = vet2)
cox.zph(fit2)  # 再次检验
```

### 步骤 8：结果提取与报告

```r
# 提取整洁结果（HR 及 95% CI）
library(broom)
result <- tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)
result

# 绘制森林图
fit_cox_full <- coxph(Surv(time, status) ~ ., data = lung)
ggforest(fit_cox_full, data = lung, main = "Hazard ratio",
         fontsize = 0.7, noDigits = 2)
```

## 代码示例

### 示例 1：Kaplan-Meier 估计与 log-rank 检验

```r
library(survival)
library(survminer)

# ---------- 数据准备 ----------
df <- lung
df$status <- ifelse(df$status == 2, 1, 0)

# ---------- 构建生存曲线 ----------
fit <- survfit(Surv(time, status) ~ 1, data = df)

# 寿命表
surv_summary(fit)
##     time n.risk n.event n.censor       surv     std.err     upper      lower
## 1      5    228       1        0 0.99561404 0.004395615 1.0000000 0.98707342
## 2     11    227       3        0 0.98245614 0.008849904 0.9996460 0.96556190

# ---------- 单组 KM 曲线 ----------
ggsurvplot(fit, conf.int = TRUE, palette = "blue",
           surv.median.line = "hv", ggtheme = theme_bw())

# ---------- 分组比较（log-rank）----------
fit_lr <- survdiff(Surv(time, status) ~ sex, data = df)
fit_lr
## Call:
## survdiff(formula = Surv(time, status) ~ sex, data = df)
##
##         N Observed Expected (O-E)^2/E (O-E)^2/V
## sex=1 138      112     91.6      4.55      10.3
## sex=2  90       53     73.4      5.68      10.3
##
##  Chisq= 10.3  on 1 degrees of freedom, p= 0.001

# broom 提取
broom::tidy(fit_lr)
## # A tibble: 2 × 4
##   sex       N   obs   exp
## 1 1       138   112  91.6
## 2 2        90    53  73.4

broom::glance(fit_lr)
## # A tibble: 1 × 3
##   statistic    df p.value
## 1      10.3     1 0.00131

# ---------- 分组 KM 曲线（含 p 值和 risk table）----------
fit_km_sex <- survfit(Surv(time, status) ~ sex, data = df)
ggsurvplot(fit_km_sex, data = df,
           surv.median.line = "hv",
           legend.title = "Sex", legend.labs = c("Male", "Female"),
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw(),
           main = "Survival curve",
           font.main = c(16, "bold", "darkblue"))
```

### 示例 2：寻找最佳切点

```r
data(myeloma)

# 寻找 DEPDC1, WHSC1, CRIM1 的最佳 cut-point
res.cut <- surv_cutpoint(myeloma, time = "time", event = "event",
                         variables = c("DEPDC1", "WHSC1", "CRIM1"))
summary(res.cut)
##        cutpoint statistic
## DEPDC1    279.8  4.275452
## WHSC1    3205.6  3.361330
## CRIM1      82.3  1.968317

# 可视化 cut-point
plot(res.cut, "DEPDC1", palette = "npg")

# 根据最佳切点重新分类
res.cat <- surv_categorize(res.cut)
head(res.cat)
##           time event DEPDC1 WHSC1 CRIM1
## GSM50986 69.24     0   high   low  high
## GSM50988 66.43     0    low   low   low

# 绘制分组 KM 曲线
fit_cat <- survfit(Surv(time, event) ~ DEPDC1, data = res.cat)
ggsurvplot(fit_cat, data = res.cat, risk.table = TRUE, conf.int = TRUE)
```

### 示例 3：多因素 Cox 回归

```r
rm(list = ls())
library(survival)
library(survminer)

# ---------- 数据准备 ----------
lung$sex     <- factor(lung$sex, labels = c("female", "male"))
lung$ph.ecog <- factor(lung$ph.ecog,
  labels = c("asymptomatic", "symptomatic", "in bed <50%", "in bed >50%"))

# ---------- 拟合 Cox 模型 ----------
fit_cox <- coxph(Surv(time, status) ~ sex + age + ph.karno, data = lung)
summary(fit_cox)
## Call:
## coxph(formula = Surv(time, status) ~ sex + age + ph.karno, data = lung)
##
##   n= 227, number of events= 164
##    (1 observation deleted due to missingness)
##
##               coef exp(coef)  se(coef)      z Pr(>|z|)
## sexmale  -0.497170  0.608249  0.167713 -2.964  0.00303 **
## age       0.012375  1.012452  0.009405  1.316  0.18821
## ph.karno -0.013322  0.986767  0.005880 -2.266  0.02348 *
##
##          exp(coef) exp(-coef) lower .95 upper .95
## sexmale     0.6082     1.6441    0.4378    0.8450
## age         1.0125     0.9877    0.9940    1.0313
## ph.karno    0.9868     1.0134    0.9755    0.9982
##
## Concordance= 0.637  (se = 0.025 )
## Likelihood ratio test= 18.81  on 3 df,   p=3e-04
## Wald test            = 18.73  on 3 df,   p=3e-04
## Score (logrank) test = 19.05  on 3 df,   p=3e-04

# ---------- 整洁结果提取 ----------
broom::tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)
## # A tibble: 3 × 7
##   term     estimate std.error statistic p.value conf.low conf.high
## 1 sexmale     0.608   0.168       -2.96 0.00303    0.438     0.845
## 2 age         1.01    0.00940      1.32 0.188      0.994     1.03
## 3 ph.karno    0.987   0.00588     -2.27 0.0235     0.975     0.998

# ---------- 森林图 ----------
fit_cox_full <- coxph(Surv(time, status) ~ ., data = lung)
ggforest(fit_cox_full, data = lung, main = "Hazard ratio",
         cpositions = c(0.01, 0.15, 0.35), fontsize = 0.7,
         refLabel = "reference", noDigits = 2)
```

### 示例 4：PH 假设检验与处理

```r
rm(list = ls())
library(survival)

# ---------- 普通 Cox 回归 ----------
fit <- coxph(Surv(time, status) ~ trt + prior + karno, data = veteran)

# PH 检验
zp <- cox.zph(fit)
zp
##         chisq df       p
## trt     0.288  1 0.59125
## prior   2.168  1 0.14087
## karno  12.138  1 0.00049
## GLOBAL 18.073  3 0.00042

# plot(zp[3]); abline(0, 0, col = "red"); abline(h = fit$coef[3], col = "green", lwd = 2, lty = 2)

# ---------- 方法一：对时间分层 ----------
vet2 <- survSplit(Surv(time, status) ~ ., data = veteran,
                  cut = c(90, 180), episode = "tgroup", id = "id")

fit2 <- coxph(Surv(tstart, time, status) ~
                trt + prior + karno:strata(tgroup), data = vet2)
fit2
##                                   coef exp(coef)  se(coef)      z        p
## trt                          -0.011025  0.989035  0.189062 -0.058    0.953
## prior                        -0.006107  0.993912  0.020355 -0.300    0.764
## karno:strata(tgroup)tgroup=1 -0.048755  0.952414  0.006222 -7.836 4.64e-15
## karno:strata(tgroup)tgroup=2  0.008050  1.008083  0.012823  0.628    0.530
## karno:strata(tgroup)tgroup=3 -0.008349  0.991686  0.014620 -0.571    0.568

cox.zph(fit2)
##                      chisq df     p
## karno:strata(tgroup)  3.04  3 0.385
## GLOBAL                8.03  5 0.154

# ---------- 方法二：时依系数变换（tt 函数）----------
fit3 <- coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
              data = veteran,
              tt = function(x, t, ...) x * log(t + 20))
fit3
##                coef exp(coef)  se(coef)      z        p
## trt        0.016478  1.016614  0.190707  0.086  0.93115
## prior     -0.009317  0.990726  0.020296 -0.459  0.64619
## karno     -0.124662  0.882795  0.028785 -4.331 1.49e-05
## tt(karno)  0.021310  1.021538  0.006607  3.225  0.00126
```

## 结果解读指南

### log-rank 检验输出

- **N**：各组人数
- **Observed**：实测死亡数
- **Expected**：在零假设（两组生存无差异）下的期望死亡数
- **(O-E)^2/E** 和 **(O-E)^2/V**：检验统计量的计算中间值
- **Chisq**：卡方值
- **df**：自由度（通常为组数 - 1）
- **p**：p 值，p < 0.05 表示两组生存曲线差异具有统计学意义

### Cox 回归输出（`summary(fit_cox)`）

- **coef**：回归系数（log HR），正值表示风险增加（危险因素），负值表示风险降低（保护因素）
- **exp(coef)**：即 HR（hazard ratio）。HR > 1 危险因素；HR < 1 保护因素；HR = 1 无影响
- **se(coef)**：回归系数的标准误
- **z**：Wald 检验的 z 统计量（coef / se(coef)）
- **Pr(>|z|)**：Wald 检验的 p 值，p < 0.05 表示该变量对生存有独立影响
- **exp(-coef)**：HR 的倒数（有时用于解释反转效应方向）
- **lower .95 / upper .95**：HR 的 95% 置信区间。区间完全在 1 的同一侧 → 有统计学意义

**论文报告示例**：多因素 Cox 回归结果显示，男性患者的死亡风险显著低于女性（HR = 0.61, 95% CI: 0.44–0.85, p = 0.003）。

### 模型整体检验

- **Concordance**：C-index（Harrell's C），反映模型预测区分度，0.5 = 无区分能力，1.0 = 完美区分，一般 > 0.7 可接受
- **Likelihood ratio test / Wald test / Score (logrank) test**：三种整体模型的显著性检验，p < 0.05 表示模型中至少有一个变量显著

### PH 假设检验（`cox.zph()`）

- 每个变量和变量整体（GLOBAL）的卡方值和 p 值
- **p < 0.05**：该变量不满足比例风险假设，需处理
- **p ≥ 0.05**：没有违反 PH 假设的证据

### broom 提取的结果

- **estimate**：当 `exponentiate = TRUE` 时为 HR 值
- **std.error**：仍为回归系数的标准误（未指数化）
- **statistic**：Wald z 值
- **p.value**：p 值
- **conf.low / conf.high**：HR 的 95% 可信区间

## 常见问题与注意事项

### Q1：log-rank 检验和 Breslow 检验有何区别？

log-rank 检验对所有时间点的死亡一视同仁（等权），对晚期差异更敏感；Breslow 检验（广义 Wilcoxon 检验）以 at-risk 样本量加权，对早期差异更敏感。医学研究中 log-rank 更常用，R 中默认即 log-rank。

### Q2：最佳 cut-point 如何确定？用中位数可以吗？

`surv_cutpoint()` 基于 log-rank 统计量最大化的原则寻找最优切点。但需注意：数据驱动的切点可能导致过拟合和 p 值膨胀，且有多个检验假设时需校正。中位数分组更简单、更可重复，被审稿人质疑的风险更低。建议优先使用中位数或公认的临床阈值。

### Q3：Cox 回归能替代 log-rank 检验吗？

单因素 Cox 回归的变量检验（Wald test）近似等价于 log-rank 检验。多因素 Cox 回归可以在校正其他协变量的同时检验某个变量的效应，是 log-rank 检验不能实现的。但 Cox 回归需要满足 PH 假设，且 log-rank 对生存曲线交叉的情况也可使用（仅报告差异有统计学意义即可）。

### Q4：PH 假设不满足怎么办？

1. **时间分层**（`survSplit()`）：将随访时间分段（如 0-90 天、90-180 天、180 天以上），每段内 HR 近似恒定
2. **时依系数（`tt()` 函数）**：在模型中加入变量与时间的交互项，如 `x * log(t + 20)`
3. **分层 Cox 回归（`strata()`）**：在 `coxph()` 公式中对不满足 PH 的变量使用 `strata()`，但该方法不能估计该变量的效应
4. **改用其他模型**：如参数化加速失效时间模型（AFT）、灵活参数模型等

### Q5：时依协变量和时依系数有何区别？

**时依协变量**：协变量本身的值在随访期间发生了变化（如生化指标、药物剂量），需要用计数过程格式（start-stop 格式）的 Surv 数据。**时依系数**：协变量的系数（效应大小）随时间变化，协变量本身不随时间变化（如性别、基线年龄），需要用 `tt()` 函数或在分层后的数据中估计。

### Q6：原始数据中 status 编码是 1＝删失、2＝死亡怎么办？

`suvival` 包内部处理规则是：0 = 删失，1 = 死亡。`survival::lung$status` 原始编码是 1 = censored，2 = dead，需转换：`status <- ifelse(status == 2, 1, 0)`。直接使用会导致结果错误。

### Q7：森林图中各哑变量共用一个参照水平，如何解读？

`ggforest()` 自动以第一个水平为参照。如 `sex` 有两个水平 "female" 和 "male"，输出中 `sexmale` 的 HR 是以 "female" 为参照。若 HR = 0.61，则说明男性死亡风险为女性的 0.61 倍。多分类变量同理，每个非参照水平都与参照比较。

### Q8：SPSS 与 R 的生存分析有何差异？

核心方法（K-M、log-rank、Cox）原理完全一致。以下为常见注意点：
- SPSS 的 Cox 回归默认将分类变量编码为 indicator（参照为第一个或最后一个类别）；R 的 `coxph()` 对 factor 类型变量自动哑变量化
- SPSS 提供了 `Time` 函数用于时依协变量；R 用 `tt()` 或 `survSplit()` 实现
- SPSS 输出更多默认表格；R 可通过 `broom::tidy()` 灵活提取结果
- C-index：SPSS 需要额外菜单操作；R 在 `summary()` 中直接给出 Concordance

### 关键提醒

- cut-point 确定后用 `surv_categorize()` 生成的分组变量进行 log-rank 检验时，p 值可能存在多重比较问题
- `cox.ph()` 检验是对每个变量的每个可能参数组合分别检验，即使 GLOBAL p > 0.05，个别变量可能仍不满足 PH
- `ggsurvplot()` 的 `pval = TRUE` 传入的是 log-rank p 值，适合比较两组；多组比较需考虑多重校正
- Cox 回归的 Concordance（C-index）仅反映区分度，不等于校准度；预测模型还需做校准曲线评估
