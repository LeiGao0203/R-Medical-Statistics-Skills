# ============================================================================
# Part 3: 方差分析、非参数检验、卡方检验
# Skills: medical-stat-anova, medical-stat-nonparametric, medical-stat-chisq
# Datasets: bone_marrow_cancer_by_country.csv, bone_marrow_cancer_risk_factors.csv,
#           bone_marrow_cancer_survival.csv, bone_marrow_cancer_treatments.csv
# ============================================================================

library(car)
library(PMCMRplus)
library(ggplot2)

data_dir <- "data"
fig_dir <- "analysis/figures"

country <- read.csv(file.path(data_dir, "bone_marrow_cancer_by_country.csv"))
survival <- read.csv(file.path(data_dir, "bone_marrow_cancer_survival.csv"))
treatments <- read.csv(file.path(data_dir, "bone_marrow_cancer_treatments.csv"))
risk_factors <- read.csv(file.path(data_dir, "bone_marrow_cancer_risk_factors.csv"))

country$Continent <- factor(country$Continent)

# ---- 3.1 单因素ANOVA：各大陆比较 ----
cat("\n========== 3.1 单因素方差分析 ==========\n")

# 方差齐性检验
cat("\n--- Levene 方差齐性检验 ---\n")
lev1 <- leveneTest(Myeloma_5Y_Survival_Pct ~ Continent, data = country)
print(lev1)
cat(sprintf("Levene: F = %.3f, p = %.4f\n", lev1$`F value`[1], lev1$`Pr(>F)`[1]))

lev2 <- leveneTest(Myeloma_Incidence_Per_100K ~ Continent, data = country)
cat(sprintf("发病率Levene: F = %.3f, p = %.4f\n", lev2$`F value`[1], lev2$`Pr(>F)`[1]))

# 单因素ANOVA - 生存率
cat("\n--- ANOVA: 骨髓瘤5年生存率 ~ 大陆 ---\n")
fit_anova1 <- aov(Myeloma_5Y_Survival_Pct ~ Continent, data = country)
s_anova1 <- summary(fit_anova1)
print(s_anova1)

# Tukey HSD 两两比较
cat("\n--- Tukey HSD 事后比较 ---\n")
tukey1 <- TukeyHSD(fit_anova1)
print(tukey1)

# 可视化Tukey HSD
png(file.path(fig_dir, "tukey_hsd_continent.png"), width = 800, height = 600)
par(las = 2, mar = c(5, 10, 4, 2))
plot(tukey1)
dev.off()
cat("已保存: tukey_hsd_continent.png\n")

# 单因素ANOVA - BMT可及性
cat("\n--- ANOVA: BMT可及性 ~ 大陆 ---\n")
fit_anova2 <- aov(BMT_Access_Score ~ Continent, data = country)
print(summary(fit_anova2))
tukey2 <- TukeyHSD(fit_anova2)
print(tukey2)

# LSD检验
cat("\n--- LSD 检验 (BMT可及性) ---\n")
lsd_res <- lsdTest(fit_anova2)
summary(lsd_res)

# ---- 3.2 Welch ANOVA (方差不齐时) ----
cat("\n========== 3.2 Welch 方差分析 ==========\n")

welch1 <- oneway.test(Hematologists_Per_Million ~ Continent, data = country, var.equal = FALSE)
cat(sprintf("Welch ANOVA 血液专家密度: F = %.3f, df1 = %.1f, df2 = %.1f, p = %.4f\n",
            welch1$statistic, welch1$parameter[1], welch1$parameter[2], welch1$p.value))

# ---- 3.3 双因素ANOVA：生存数据(Cancer_Type × Income_Region) ----
cat("\n========== 3.3 双因素方差分析 (癌种×收入水平) ==========\n")

survival$Cancer_Type <- factor(survival$Cancer_Type)
survival$Income_Region <- factor(survival$Income_Region,
  levels = c("High-Income", "Upper-Middle", "Lower-Middle", "Low-Income"))

cat("\n--- 双因素ANOVA: 5年生存率 ~ 癌种 + 收入水平 ---\n")
fit_two <- aov(Five_Year_Survival_Pct ~ Cancer_Type + Income_Region, data = survival)
print(summary(fit_two))

# 含交互项
cat("\n--- 含交互项双因素ANOVA ---\n")
fit_two_int <- aov(Five_Year_Survival_Pct ~ Cancer_Type * Income_Region, data = survival)
print(summary(fit_two_int))

# 各收入水平组均值
cat("\n各收入水平5年生存率均值:\n")
tapply(survival$Five_Year_Survival_Pct, survival$Income_Region, mean)

# 可视化
p_anova <- ggplot(survival, aes(Cancer_Type, Five_Year_Survival_Pct, fill = Income_Region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "癌症类型", y = "5年生存率(%)", fill = "收入水平",
       title = "不同癌症类型和收入水平的5年生存率") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "bar_survival_income_cancer.pdf"), p_anova, width = 12, height = 6)
cat("已保存: bar_survival_income_cancer.pdf\n")

# ---- 3.4 治疗方案ORR/CRR比较 ----
cat("\n========== 3.4 治疗应答率比较 ==========\n")

treatments$Treatment_Category <- factor(treatments$Treatment_Category)

cat("\n--- ANOVA: Overall Response Rate ~ 治疗类别 ---\n")
fit_trt <- aov(Overall_Response_Rate_Pct ~ Treatment_Category, data = treatments)
print(summary(fit_trt))

# 各治疗类别的ORR均值
cat("\n各治疗类别总体应答率均值:\n")
orr_mean <- tapply(treatments$Overall_Response_Rate_Pct, treatments$Treatment_Category, mean)
print(sort(orr_mean, decreasing = TRUE))

# Kruskal-Wallis（样本量小备选）
cat("\n--- Kruskal-Wallis: ORR ~ 治疗类别 ---\n")
kw_trt <- kruskal.test(Overall_Response_Rate_Pct ~ Treatment_Category, data = treatments)
cat(sprintf("Kruskal-Wallis: chi² = %.3f, df = %d, p = %.4f\n",
            kw_trt$statistic, kw_trt$parameter, kw_trt$p.value))

# 可视化
p_trt <- ggplot(treatments, aes(reorder(Treatment_Category, Overall_Response_Rate_Pct), Overall_Response_Rate_Pct)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_text(aes(label = Overall_Response_Rate_Pct), hjust = -0.2, size = 3) +
  coord_flip() + labs(x = "治疗类别", y = "总应答率(%)", title = "各类治疗的总体应答率(ORR)") +
  ylim(0, 100) + theme_classic()
ggsave(file.path(fig_dir, "bar_trt_response.pdf"), p_trt, width = 10, height = 6)
cat("已保存: bar_trt_response.pdf\n")

# ---- 3.5 卡方检验：风险因素分类分析 ----
cat("\n========== 3.5 卡方检验 ==========\n")

# 风险因素类别 vs 证据等级
risk_factors$Category <- factor(risk_factors$Category)
risk_factors$Evidence_Level <- factor(risk_factors$Evidence_Level)

cat("\n--- 风险因素类别分布 ---\n")
cat_tab <- table(risk_factors$Category)
print(cat_tab)

cat("\n--- 证据等级分布 ---\n")
evd_tab <- table(risk_factors$Evidence_Level)
print(evd_tab)

# 构建R×C列联表
if (length(unique(risk_factors$Category)) >= 3) {
  cross_tab <- table(risk_factors$Category, risk_factors$Evidence_Level)
  cat("\n--- 风险类别 vs 证据等级 列联表 ---\n")
  print(cross_tab)

  chisq_risk <- chisq.test(cross_tab, correct = FALSE)
  if (!is.null(chisq_risk$p.value)) {
    cat(sprintf("卡方检验: chi² = %.3f, df = %d, p = %.4f\n",
                chisq_risk$statistic, chisq_risk$parameter, chisq_risk$p.value))
  }
}

# 可改变性 vs 相对风险
cat("\n--- Fisher检验: 可改变 vs 高RR (RR>2.0) ---\n")
risk_factors$High_RR <- risk_factors$Relative_Risk > 2.0
tb_fish <- table(risk_factors$Modifiable, risk_factors$High_RR)
print(tb_fish)
fish_res <- fisher.test(tb_fish)
cat(sprintf("Fisher精确检验: p = %.4f, OR = %.3f\n", fish_res$p.value, fish_res$estimate))

# ---- 3.6 列联系数（关联强度） ----
cat("\n--- 列联相关系数 ---\n")
if (exists("chisq_risk") && !is.null(chisq_risk$statistic)) {
  n_total <- sum(cross_tab)
  phi <- sqrt(chisq_risk$statistic / n_total)
  cat(sprintf("Phi系数 = %.3f\n", phi))
  V <- sqrt(chisq_risk$statistic / (n_total * min(nrow(cross_tab) - 1, ncol(cross_tab) - 1)))
  cat(sprintf("Cramer's V = %.3f\n", V))
}

cat("\n===== Part 3 完成 =====\n")
