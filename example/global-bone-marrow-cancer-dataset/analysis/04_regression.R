# ============================================================================
# Part 4: 多元线性回归 与 Logistic回归
# Skills: medical-stat-multi-reg, medical-stat-logistic-reg
# Dataset: bone_marrow_cancer_by_country.csv
# ============================================================================

library(GGally)
library(performance)
library(car)
library(gvlma)
library(MASS)
library(glmnet)

data_dir <- "data"
fig_dir <- "analysis/figures"

country <- read.csv(file.path(data_dir, "bone_marrow_cancer_by_country.csv"))
country$Continent <- factor(country$Continent)

# ---- 4.1 数据探索 ----
cat("\n========== 4.1 散点图矩阵 ==========\n")
vars_for_reg <- c("Myeloma_5Y_Survival_Pct", "BMT_Access_Score", "Hematologists_Per_Million",
                  "Myeloma_Incidence_Per_100K", "BMT_Centers", "Population_M")
png(file.path(fig_dir, "pairs_plot.png"), width = 1200, height = 1000)
ggpairs(country[, vars_for_reg]) + theme_bw()
dev.off()
cat("已保存: pairs_plot.png\n")

# ---- 4.2 多元线性回归：预测5年生存率 ----
cat("\n========== 4.2 多元线性回归 ==========\n")

# 完整模型
f_all <- lm(Myeloma_5Y_Survival_Pct ~ BMT_Access_Score + Hematologists_Per_Million +
            Myeloma_Incidence_Per_100K + BMT_Centers + Population_M,
            data = country)
s_all <- summary(f_all)
cat("\n--- 完整模型 ---\n")
print(s_all)

# 模型评价
perf <- model_performance(f_all)
cat(sprintf("\nR² = %.3f, 调整R² = %.3f, AIC = %.1f, BIC = %.1f\n",
            perf$R2, perf$R2_adjusted, perf$AIC, perf$BIC))

# ---- 4.3 回归诊断 ----
cat("\n========== 4.3 回归诊断 ==========\n")

# 综合诊断图
png(file.path(fig_dir, "reg_diagnostics.png"), width = 1200, height = 1000)
par(mfrow = c(2, 2))
plot(f_all)
par(mfrow = c(1, 1))
dev.off()
cat("已保存: reg_diagnostics.png\n")

# 统计检验
gvmodel <- gvlma(f_all)
cat("\n--- gvlma 综合诊断 ---\n")
print(summary(gvmodel))

# 正态性
check_normality(f_all)

# 独立性
dw <- durbinWatsonTest(f_all)
cat(sprintf("Durbin-Watson: statistic = %.3f, p = %.3f\n", dw$dw, dw$p))

# 多重共线性
cat("\n--- VIF 多重共线性诊断 ---\n")
vif_vals <- vif(f_all)
print(vif_vals)
cat("VIF > 5 的变量:", names(vif_vals)[vif_vals > 5], "\n")

# ---- 4.4 变量筛选（逐步回归） ----
cat("\n========== 4.4 逐步回归 ==========\n")

# 后退法
step_back <- stepAIC(f_all, direction = "backward", trace = TRUE)
cat("\n--- 后退法最优模型 ---\n")
s_step <- summary(step_back)
print(s_step)

# 整理结果
library(broom)
cat("\n--- 最终模型回归系数 ---\n")
tidy_result <- tidy(step_back, conf.int = TRUE)
print(tidy_result)

# ---- 4.5 含大陆哑变量的扩展模型 ----
cat("\n========== 4.5 含大陆哑变量的扩展模型 ==========\n")

f_continent <- lm(Myeloma_5Y_Survival_Pct ~ BMT_Access_Score + Hematologists_Per_Million + Continent,
                  data = country)
cat("\n--- 含大陆哑变量模型 ---\n")
s_cont <- summary(f_continent)
print(s_cont)
cat(sprintf("调整R² = %.3f (含大陆: %.3f)\n", s_all$adj.r.squared, s_cont$adj.r.squared))

# ---- 4.6 Logistic回归：预测高生存率 vs 低生存率 ----
cat("\n========== 4.6 Logistic回归 ==========\n")

# 创建二分类结局
country$High_Survival <- ifelse(country$Myeloma_5Y_Survival_Pct > 50, 1, 0)
country$High_Survival <- factor(country$High_Survival, levels = c(0, 1), labels = c("低生存率", "高生存率"))
table(country$High_Survival)

# 拟合Logistic回归
logit_fit <- glm(High_Survival ~ BMT_Access_Score + Hematologists_Per_Million +
                 Myeloma_Incidence_Per_100K + BMT_Centers,
                 data = country, family = binomial())
s_logit <- summary(logit_fit)
cat("\n--- Logistic回归结果 ---\n")
print(s_logit)

# OR值和95%CI
cat("\n--- OR值及95% CI ---\n")
or_result <- tidy(logit_fit, exponentiate = TRUE, conf.int = TRUE)
print(or_result)

# 伪R²
library(DescTools)
pseudo <- PseudoR2(logit_fit, which = c("McFadden", "CoxSnell", "Nagelkerke"))
cat(sprintf("McFadden R² = %.3f, CoxSnell R² = %.3f, Nagelkerke R² = %.3f\n",
            pseudo["McFadden"], pseudo["CoxSnell"], pseudo["Nagelkerke"]))

# 似然比检验
logit_null <- glm(High_Survival ~ 1, data = country, family = binomial())
lrt <- anova(logit_null, logit_fit, test = "Chisq")
cat(sprintf("似然比检验: chi² = %.2f, p = %.4f\n", lrt$Deviance[2], lrt$`Pr(>Chi)`[2]))

# 预测概率
country$pred_prob <- predict(logit_fit, type = "response")
country$pred_class <- ifelse(country$pred_prob > 0.5, "高生存率", "低生存率")

# 混淆矩阵
cat("\n--- 混淆矩阵 ---\n")
cm <- table(Actual = country$High_Survival, Predicted = country$pred_class)
print(cm)
accuracy <- sum(diag(cm)) / sum(cm)
cat(sprintf("准确率 = %.1f%%\n", accuracy * 100))

# ---- 4.7 逐步Logistic回归 ----
cat("\n--- 逐步Logistic回归（后退法）---\n")
step_logit <- step(logit_fit, direction = "backward", trace = TRUE)
cat("\n最优Logistic模型:\n")
print(summary(step_logit))

# ---- 4.8 岭回归和LASSO（处理共线性） ----
cat("\n========== 4.8 岭回归和LASSO ==========\n")

x <- as.matrix(country[, c("BMT_Access_Score", "Hematologists_Per_Million",
                           "Myeloma_Incidence_Per_100K", "BMT_Centers", "Population_M")])
y <- country$Myeloma_5Y_Survival_Pct

# 岭回归
cv_ridge <- cv.glmnet(x, y, alpha = 0)
cat(sprintf("岭回归最佳lambda: %.4f\n", cv_ridge$lambda.min))

# LASSO
cv_lasso <- cv.glmnet(x, y, alpha = 1)
cat(sprintf("LASSO最佳lambda: %.4f\n", cv_lasso$lambda.min))

# LASSO系数
lasso_coef <- coef(cv_lasso, s = "lambda.min")
cat("\nLASSO系数:\n")
print(lasso_coef)

# ---- 4.9 结果汇总 ----
cat("\n========== 4.9 多因素回归结果汇总 ==========\n")
cat("最终多元线性回归模型:\n")
print(tidy_result[, c("term", "estimate", "std.error", "statistic", "p.value")])

cat("\n最终Logistic回归模型:\n")
print(tidy(step_logit, exponentiate = TRUE, conf.int = TRUE)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

cat("\n===== Part 4 完成 =====\n")
