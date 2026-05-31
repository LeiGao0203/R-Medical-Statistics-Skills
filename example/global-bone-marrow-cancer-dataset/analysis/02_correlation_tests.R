# ============================================================================
# Part 2: 相关分析、t检验、非参数检验
# Skills: medical-stat-correlation, medical-stat-ttest, medical-stat-nonparametric
# Dataset: bone_marrow_cancer_by_country.csv
# ============================================================================

library(ggplot2)

data_dir <- "data"
fig_dir <- "analysis/figures"

country <- read.csv(file.path(data_dir, "bone_marrow_cancer_by_country.csv"), stringsAsFactors = FALSE)
country$Continent <- factor(country$Continent)

# ---- 2.1 Pearson 相关分析 ----
cat("\n========== 2.1 Pearson 相关分析 ==========\n")

cor_vars <- c("Myeloma_Incidence_Per_100K", "Leukemia_Incidence_Per_100K",
              "Myeloma_5Y_Survival_Pct", "Leukemia_5Y_Survival_Pct",
              "BMT_Access_Score", "BMT_Centers", "Hematologists_Per_Million")

cat("\n--- 相关矩阵 ---\n")
cor_matrix <- cor(country[, cor_vars], use = "complete.obs")
print(round(cor_matrix, 3))

# 关键相关对分析
cat("\n--- 关键相关对假设检验 ---\n")

# 1) 骨髓瘤发病率 vs 5年生存率
ct1 <- cor.test(~ Myeloma_Incidence_Per_100K + Myeloma_5Y_Survival_Pct, data = country)
cat(sprintf("骨髓瘤发病率 vs 5年生存率: r = %.3f, 95%%CI [%.3f, %.3f], p = %.4f\n",
            ct1$estimate, ct1$conf.int[1], ct1$conf.int[2], ct1$p.value))

# 2) BMT可及性 vs 5年生存率
ct2 <- cor.test(~ BMT_Access_Score + Myeloma_5Y_Survival_Pct, data = country)
cat(sprintf("BMT可及性 vs 骨髓瘤5年生存率: r = %.3f, 95%%CI [%.3f, %.3f], p = %.4f\n",
            ct2$estimate, ct2$conf.int[1], ct2$conf.int[2], ct2$p.value))

# 3) 血液专家密度 vs 5年生存率
ct3 <- cor.test(~ Hematologists_Per_Million + Myeloma_5Y_Survival_Pct, data = country)
cat(sprintf("血液专家密度 vs 骨髓瘤5年生存率: r = %.3f, 95%%CI [%.3f, %.3f], p = %.4f\n",
            ct3$estimate, ct3$conf.int[1], ct3$conf.int[2], ct3$p.value))

# 4) 发病率 vs BMT可及性
ct4 <- cor.test(~ Myeloma_Incidence_Per_100K + BMT_Access_Score, data = country)
cat(sprintf("骨髓瘤发病率 vs BMT可及性: r = %.3f, 95%%CI [%.3f, %.3f], p = %.4f\n",
            ct4$estimate, ct4$conf.int[1], ct4$conf.int[2], ct4$p.value))

# ---- 2.2 Spearman 秩相关分析 ----
cat("\n========== 2.2 Spearman 秩相关分析 ==========\n")

# 五大陆之间的秩相关（用大陆均值）
continent_mean <- aggregate(cbind(Myeloma_Incidence_Per_100K, Myeloma_5Y_Survival_Pct) ~ Continent,
                            data = country, mean)
cat("\n大陆均值:\n")
print(continent_mean)

sr <- cor.test(continent_mean$Myeloma_Incidence_Per_100K, continent_mean$Myeloma_5Y_Survival_Pct,
               method = "spearman")
cat(sprintf("大陆水平 Spearman: rho = %.3f, p = %.4f\n", sr$estimate, sr$p.value))

# ---- 2.3 简单线性回归 ----
cat("\n========== 2.3 简单线性回归 ==========\n")

cat("\n--- 模型1: 5年生存率 ~ BMT可及性评分 ---\n")
fit1 <- lm(Myeloma_5Y_Survival_Pct ~ BMT_Access_Score, data = country)
s1 <- summary(fit1)
print(s1)
cat(sprintf("R² = %.3f, F = %.2f, p = %.4f\n", s1$r.squared, s1$fstatistic[1],
            pf(s1$fstatistic[1], s1$fstatistic[2], s1$fstatistic[3], lower.tail = FALSE)))

cat("\n--- 模型2: 5年生存率 ~ 血液专家密度 ---\n")
fit2 <- lm(Myeloma_5Y_Survival_Pct ~ Hematologists_Per_Million, data = country)
s2 <- summary(fit2)
print(s2)
cat(sprintf("R² = %.3f, F = %.2f, p = %.4f\n", s2$r.squared, s2$fstatistic[1],
            pf(s2$fstatistic[1], s2$fstatistic[2], s2$fstatistic[3], lower.tail = FALSE)))

# ---- 2.4 残差诊断 ----
cat("\n--- 残差诊断 ---\n")
shapiro.test(residuals(fit1))
cat(sprintf("残差正态性 Shapiro-Wilk: W = %.4f, p = %.4f\n",
            shapiro.test(residuals(fit1))$statistic, shapiro.test(residuals(fit1))$p.value))

# 残差图
png(file.path(fig_dir, "residual_diagnostics.png"), width = 1000, height = 800)
par(mfrow = c(2, 2))
plot(fit1)
par(mfrow = c(1, 1))
dev.off()
cat("已保存: residual_diagnostics.png\n")

# ---- 2.5 t检验：高收入 vs 低收入国家比较 ----
cat("\n========== 2.5 独立样本t检验 ==========\n")

# 按生存率中位数分组
country$Survival_Group <- ifelse(country$Myeloma_5Y_Survival_Pct >=
                                  median(country$Myeloma_5Y_Survival_Pct, na.rm = TRUE),
                                 "高生存率国家", "低生存率国家")
country$Survival_Group <- factor(country$Survival_Group)

# 1) 两组BMT可及性比较
cat("\n--- BMT可及性: 高生存率 vs 低生存率 ---\n")
shapiro.test(country$BMT_Access_Score[country$Survival_Group == "高生存率国家"])
shapiro.test(country$BMT_Access_Score[country$Survival_Group == "低生存率国家"])

# 方差齐性检验
var_test <- var.test(BMT_Access_Score ~ Survival_Group, data = country)
cat(sprintf("方差齐性F检验: F = %.3f, p = %.4f\n", var_test$statistic, var_test$p.value))

tt1 <- t.test(BMT_Access_Score ~ Survival_Group, data = country, var.equal = TRUE)
cat(sprintf("Student t检验: t = %.3f, df = %d, p = %.4f\n", tt1$statistic, tt1$parameter, tt1$p.value))
cat(sprintf("均值差 95%%CI: [%.2f, %.2f]\n", tt1$conf.int[1], tt1$conf.int[2]))

# 2) 两组血液专家密度比较
cat("\n--- 血液专家密度: 高生存率 vs 低生存率 ---\n")
tt2 <- t.test(Hematologists_Per_Million ~ Survival_Group, data = country)
cat(sprintf("Welch t检验: t = %.3f, df = %.1f, p = %.4f\n", tt2$statistic, tt2$parameter, tt2$p.value))

# ---- 2.6 非参数检验 (Mann-Whitney U) ----
cat("\n========== 2.6 Mann-Whitney U检验 ==========\n")

# 创建高/低收入区域分组（按BMT可及性中位数）
country$Income_Level <- ifelse(country$BMT_Access_Score >= median(country$BMT_Access_Score, na.rm = TRUE),
                               "High_Access", "Low_Access")
country$Income_Level <- factor(country$Income_Level)

mw1 <- wilcox.test(Myeloma_Incidence_Per_100K ~ Income_Level, data = country)
cat(sprintf("Mann-Whitney: 骨髓瘤发病率 按BMT可及性分组, W = %.1f, p = %.4f\n",
            mw1$statistic, mw1$p.value))

mw2 <- wilcox.test(Myeloma_5Y_Survival_Pct ~ Income_Level, data = country)
cat(sprintf("Mann-Whitney: 骨髓瘤5年生存率 按BMT可及性分组, W = %.1f, p = %.6f\n",
            mw2$statistic, mw2$p.value))

# Kruskal-Wallis: 各大陆骨髓瘤生存率比较
cat("\n--- Kruskal-Wallis: 各大陆骨髓瘤生存率 ---\n")
kw1 <- kruskal.test(Myeloma_5Y_Survival_Pct ~ Continent, data = country)
cat(sprintf("Kruskal-Wallis: chi² = %.3f, df = %d, p = %.6f\n",
            kw1$statistic, kw1$parameter, kw1$p.value))

# 各大陆的中位数生存率
cat("\n各大陆骨髓瘤5年生存率中位数:\n")
med <- tapply(country$Myeloma_5Y_Survival_Pct, country$Continent, median, na.rm = TRUE)
print(sort(med, decreasing = TRUE))

# 多组两两比较
library(PMCMRplus)
country$Continent <- factor(country$Continent)
nemenyi <- kwAllPairsNemenyiTest(Myeloma_5Y_Survival_Pct ~ Continent, data = country)
cat("\nNemenyi 事后两两比较:\n")
summary(nemenyi)

cat("\n===== Part 2 完成 =====\n")
