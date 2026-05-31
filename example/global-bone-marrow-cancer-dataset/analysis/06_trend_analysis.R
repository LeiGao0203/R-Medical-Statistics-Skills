# ============================================================================
# Part 6: 趋势分析、多项式回归、综合可视化和报告生成
# Skills: medical-stat-polynomial, medical-stat-p-for-trend
# Datasets: bone_marrow_cancer_trends_2000_2026.csv
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(broom)

data_dir <- "data"
fig_dir <- "analysis/figures"

trends <- read.csv(file.path(data_dir, "bone_marrow_cancer_trends_2000_2026.csv"))
trends$Year <- as.numeric(trends$Year)

# ---- 6.1 时间序列趋势图 ----
cat("\n========== 6.1 时间序列趋势可视化 ==========\n")

# 全球新发病例趋势
trends_long <- trends %>%
  select(Year, Global_Myeloma_New_Cases, Global_Leukemia_New_Cases, Global_MDS_New_Cases) %>%
  pivot_longer(-Year, names_to = "Cancer", values_to = "Cases") %>%
  mutate(Cancer = recode(Cancer,
    Global_Myeloma_New_Cases = "骨髓瘤",
    Global_Leukemia_New_Cases = "白血病",
    Global_MDS_New_Cases = "MDS"))

p_cases <- ggplot(trends_long, aes(Year, Cases, color = Cancer, linetype = Cancer)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  labs(x = "年份", y = "全球新发病例数", title = "骨髓相关癌症全球新发病例趋势 (2000-2026)") +
  scale_color_manual(values = c("骨髓瘤" = "#2E9FDF", "白血病" = "#E7B800", "MDS" = "#FC4E07")) +
  theme_classic()
ggsave(file.path(fig_dir, "trend_new_cases.pdf"), p_cases, width = 10, height = 6)
cat("已保存: trend_new_cases.pdf\n")

# 生存率趋势
trends_long2 <- trends %>%
  select(Year, Myeloma_5Y_Survival_Pct, Leukemia_5Y_Survival_Pct, MDS_5Y_Survival_Pct) %>%
  pivot_longer(-Year, names_to = "Cancer", values_to = "Survival") %>%
  mutate(Cancer = recode(Cancer,
    Myeloma_5Y_Survival_Pct = "骨髓瘤",
    Leukemia_5Y_Survival_Pct = "白血病",
    MDS_5Y_Survival_Pct = "MDS"))

p_surv <- ggplot(trends_long2, aes(Year, Survival, color = Cancer)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.5) +
  labs(x = "年份", y = "5年生存率(%)", title = "全球骨髓相关癌症5年生存率趋势 (2000-2026)") +
  scale_color_manual(values = c("骨髓瘤" = "#2E9FDF", "白血病" = "#E7B800", "MDS" = "#FC4E07")) +
  theme_classic()
ggsave(file.path(fig_dir, "trend_survival.pdf"), p_surv, width = 10, height = 6)
cat("已保存: trend_survival.pdf\n")

# BMT手术趋势
p_bmt <- ggplot(trends, aes(Year, Global_BMT_Procedures)) +
  geom_line(color = "#00AFBB", linewidth = 1.2) +
  geom_point(color = "#00AFBB", size = 2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "darkblue") +
  labs(x = "年份", y = "全球BMT手术数量", title = "全球骨髓移植手术数量趋势") +
  scale_y_continuous(labels = scales::comma) +
  theme_classic()
ggsave(file.path(fig_dir, "trend_bmt.pdf"), p_bmt, width = 10, height = 6)
cat("已保存: trend_bmt.pdf\n")

# 综合面板图
p_combined <- p_cases / (p_surv + p_bmt) +
  plot_annotation(title = "全球骨髓相关癌症趋势综合分析",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
ggsave(file.path(fig_dir, "trend_combined.pdf"), p_combined, width = 14, height = 14)
cat("已保存: trend_combined.pdf\n")

# ---- 6.2 线性趋势检验 ----
cat("\n========== 6.2 线性趋势检验 ==========\n")

# 骨髓瘤新发病例的线性趋势
cat("\n--- 骨髓瘤新发病例线性回归 ---\n")
fit_myeloma_cases <- lm(Global_Myeloma_New_Cases ~ Year, data = trends)
s_mc <- summary(fit_myeloma_cases)
annual_increase <- coef(fit_myeloma_cases)["Year"]
cat(sprintf("每年增加: %.0f 例\n", annual_increase))
cat(sprintf("R² = %.3f, p = %.4f\n", s_mc$r.squared,
            coef(s_mc)[2, "Pr(>|t|)"]))
cat(sprintf("2000年预测值: %.0f, 2026年预测值: %.0f\n",
            predict(fit_myeloma_cases, data.frame(Year = 2000)),
            predict(fit_myeloma_cases, data.frame(Year = 2026))))

# 生存率趋势
cat("\n--- 骨髓瘤5年生存率线性趋势 ---\n")
fit_surv <- lm(Myeloma_5Y_Survival_Pct ~ Year, data = trends)
s_surv <- summary(fit_surv)
annual_surv_inc <- coef(fit_surv)["Year"]
cat(sprintf("每年生存率增加: %.2f 百分点\n", annual_surv_inc))
cat(sprintf("R² = %.3f, p = %.4f\n", s_surv$r.squared,
            coef(s_surv)[2, "Pr(>|t|)"]))
cat(sprintf("2000年: %.1f%%, 2026年: %.1f%%\n",
            predict(fit_surv, data.frame(Year = 2000)),
            predict(fit_surv, data.frame(Year = 2026))))

# BMT手术趋势
cat("\n--- BMT手术数量线性趋势 ---\n")
fit_bmt <- lm(Global_BMT_Procedures ~ Year, data = trends)
s_bmt <- summary(fit_bmt)
annual_bmt_inc <- coef(fit_bmt)["Year"]
cat(sprintf("每年BMT增加: %.0f 例\n", annual_bmt_inc))
cat(sprintf("R² = %.3f, p = %.4f\n", s_bmt$r.squared,
            coef(s_bmt)[2, "Pr(>|t|)"]))
cat(sprintf("2000年: %.0f, 2026年: %.0f\n",
            predict(fit_bmt, data.frame(Year = 2000)),
            predict(fit_bmt, data.frame(Year = 2026))))

# ---- 6.3 多项式回归（非线性趋势） ----
cat("\n========== 6.3 多项式回归 ==========\n")

# 二次项回归
cat("\n--- 二次项回归: 骨髓瘤发病趋势 ---\n")
fit_quad <- lm(Global_Myeloma_New_Cases ~ poly(Year, 2), data = trends)
s_quad <- summary(fit_quad)
print(s_quad)
cat(sprintf("二次项R² = %.3f (线性R² = %.3f)\n", s_quad$r.squared, s_mc$r.squared))

# 模型比较
anova_comp <- anova(fit_myeloma_cases, fit_quad)
cat(sprintf("二次项 vs 线性: F = %.2f, p = %.4f\n",
            anova_comp$F[2], anova_comp$`Pr(>F)`[2]))

# 生存率的二次项趋势
cat("\n--- 二次项回归: 生存率趋势 ---\n")
fit_quad_surv <- lm(Myeloma_5Y_Survival_Pct ~ poly(Year, 2), data = trends)
s_quad_surv <- summary(fit_quad_surv)
cat(sprintf("二次项R² = %.3f (线性R² = %.3f)\n", s_quad_surv$r.squared, s_surv$r.squared))

# ---- 6.4 对数转换趋势分析 ----
cat("\n========== 6.4 对数转换趋势 ==========\n")

fit_log <- lm(log(Global_Myeloma_New_Cases) ~ Year, data = trends)
s_log <- summary(fit_log)
cat(sprintf("对数线性模型: R² = %.3f, 年增长率 = %.2f%%\n",
            s_log$r.squared, (exp(coef(fit_log)["Year"]) - 1) * 100))

# ---- 6.5 发病率 vs 死亡率比率趋势 ----
cat("\n========== 6.5 发病率/死亡率比趋势 ==========\n")

trends$incidence_mortality_ratio <- trends$Global_Myeloma_New_Cases / trends$Global_Myeloma_Deaths
fit_ratio <- lm(incidence_mortality_ratio ~ Year, data = trends)
s_ratio <- summary(fit_ratio)
cat(sprintf("发病率/死亡率比: 每年增加 %.4f, R² = %.3f, p = %.6f\n",
            coef(fit_ratio)["Year"], s_ratio$r.squared, coef(s_ratio)[2, "Pr(>|t|)"]))

# ---- 6.6 综合分析报告 ----
cat("\n============================================================\n")
cat("========== 骨髓相关癌症全球趋势综合分析报告 ==========\n")
cat("============================================================\n\n")

cat("1. 全球新发病例趋势:\n")
cat(sprintf("   - 骨髓瘤: 每年增加约 %.0f 例 (p = %.4f)\n",
            coef(fit_myeloma_cases)["Year"], coef(s_mc)[2, "Pr(>|t|)"]))
cat(sprintf("   - 从 %.0f (2000) 增至 %.0f (2026), 增长 %.1f%%\n",
            predict(fit_myeloma_cases, data.frame(Year = 2000)),
            predict(fit_myeloma_cases, data.frame(Year = 2026)),
            (predict(fit_myeloma_cases, data.frame(Year = 2026)) /
             predict(fit_myeloma_cases, data.frame(Year = 2000)) - 1) * 100))

cat("\n2. 生存率改善:\n")
cat(sprintf("   - 骨髓瘤5年生存率从 %.1f%% (2000) 提升至 %.1f%% (2026)\n",
            trends$Myeloma_5Y_Survival_Pct[1],
            trends$Myeloma_5Y_Survival_Pct[nrow(trends)]))
cat(sprintf("   - 白血病5年生存率从 %.1f%% (2000) 提升至 %.1f%% (2026)\n",
            trends$Leukemia_5Y_Survival_Pct[1],
            trends$Leukemia_5Y_Survival_Pct[nrow(trends)]))
cat(sprintf("   - MDS 5年生存率从 %.1f%% (2000) 提升至 %.1f%% (2026)\n",
            trends$MDS_5Y_Survival_Pct[1],
            trends$MDS_5Y_Survival_Pct[nrow(trends)]))

cat("\n3. BMT手术增长:\n")
cat(sprintf("   - 全球BMT手术从 %.0f (2000) 增至 %.0f (2026), 增长 %.1f%%\n",
            trends$Global_BMT_Procedures[1],
            trends$Global_BMT_Procedures[nrow(trends)],
            (trends$Global_BMT_Procedures[nrow(trends)] /
             trends$Global_BMT_Procedures[1] - 1) * 100))

cat("\n4. 骨髓瘤病死率变化:\n")
cat(sprintf("   - 从 %.2f (2000) 降至 %.2f (2026), 改善 %.1f%%\n",
            trends$myeloma_mortality_rate[1],
            trends$myeloma_mortality_rate[nrow(trends)],
            (trends$myeloma_mortality_rate[1] - trends$myeloma_mortality_rate[nrow(trends)]) /
             trends$myeloma_mortality_rate[1] * 100))

cat("\n5. 关键里程碑年份:\n")
milestones <- trends[trends$Key_Milestone != "" & !is.na(trends$Key_Milestone),
                     c("Year", "Key_Milestone")]
for (i in 1:min(6, nrow(milestones))) {
  cat(sprintf("   - %d年: %s\n", milestones$Year[i], milestones$Key_Milestone[i]))
}

cat("\n============================================================\n")
cat("===== Part 6 趋势分析完成 =====\n")
cat("============================================================\n")
