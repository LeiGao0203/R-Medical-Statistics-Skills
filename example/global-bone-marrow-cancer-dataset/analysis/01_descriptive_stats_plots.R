# ============================================================================
# Part 1: 描述性统计、三线表、统计绘图
# Skills: medical-stat-three-line-table, medical-stat-stat-plot
# Dataset: bone_marrow_cancer_by_country.csv, bone_marrow_cancer_types.csv
# ============================================================================

library(compareGroups)
library(ggplot2)
library(dplyr)
library(patchwork)

data_dir <- "data"
fig_dir <- "analysis/figures"

# ---- 1.1 Load data ----
cat("\n========== 1.1 数据加载 ==========\n")
country <- read.csv(file.path(data_dir, "bone_marrow_cancer_by_country.csv"), stringsAsFactors = FALSE)
cancer_types <- read.csv(file.path(data_dir, "bone_marrow_cancer_types.csv"), stringsAsFactors = FALSE)
risk_factors <- read.csv(file.path(data_dir, "bone_marrow_cancer_risk_factors.csv"), stringsAsFactors = FALSE)
survival <- read.csv(file.path(data_dir, "bone_marrow_cancer_survival.csv"), stringsAsFactors = FALSE)
treatments <- read.csv(file.path(data_dir, "bone_marrow_cancer_treatments.csv"), stringsAsFactors = FALSE)
trends <- read.csv(file.path(data_dir, "bone_marrow_cancer_trends_2000_2026.csv"), stringsAsFactors = FALSE)

str(country)
summary(country)
cat("各国数据集：", nrow(country), "行 ×", ncol(country), "列\n")

# ---- 1.2 描述性统计 ----
cat("\n========== 1.2 描述性统计 ==========\n")
# Continuous variables
cont_vars <- c("Population_M", "Myeloma_New_Cases", "Leukemia_New_Cases",
               "Myeloma_Deaths", "Leukemia_Deaths", "Myeloma_Incidence_Per_100K",
               "Leukemia_Incidence_Per_100K", "Myeloma_5Y_Survival_Pct",
               "Leukemia_5Y_Survival_Pct", "BMT_Centers", "BMT_Access_Score",
               "Hematologists_Per_Million")
desc_stats <- data.frame(
  Variable = cont_vars,
  Mean = sapply(country[cont_vars], mean, na.rm = TRUE),
  SD = sapply(country[cont_vars], sd, na.rm = TRUE),
  Median = sapply(country[cont_vars], median, na.rm = TRUE),
  Min = sapply(country[cont_vars], min, na.rm = TRUE),
  Max = sapply(country[cont_vars], max, na.rm = TRUE)
)
print(desc_stats, digits = 2)

# ---- 1.3 三线表 (Table 1 by Continent) ----
cat("\n========== 1.3 三线表 (compareGroups) ==========\n")
country$Continent <- factor(country$Continent)

# 选择需要描述的变量
tab_vars <- c("Population_M", "Myeloma_Incidence_Per_100K", "Leukemia_Incidence_Per_100K",
               "Myeloma_5Y_Survival_Pct", "Leukemia_5Y_Survival_Pct",
               "BMT_Centers", "BMT_Access_Score", "Hematologists_Per_Million")

# 使用compareGroups制作三线表
res <- compareGroups(Continent ~ ., data = country[, c("Continent", tab_vars)], method = NA)
restab <- createTable(res, show.all = TRUE)
print(restab, which.table = "descr")

# 导出
export2csv(restab, file = "analysis/table1_by_continent.csv")
cat("三线表已导出: analysis/table1_by_continent.csv\n")

# ---- 1.4 统计绘图 ----
cat("\n========== 1.4 统计绘图 ==========\n")

# 1.4.1 箱线图：各大陆骨髓瘤发病率
p1 <- ggplot(country, aes(Continent, Myeloma_Incidence_Per_100K)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(x = "大陆", y = "骨髓瘤发病率(每10万)", title = "各大陆骨髓瘤年龄标化发病率") +
  theme_classic()
ggsave(file.path(fig_dir, "boxplot_myeloma_by_continent.pdf"), p1, width = 8, height = 6)
cat("已保存: boxplot_myeloma_by_continent.pdf\n")

# 1.4.2 直方图：5年生存率分布
p2 <- ggplot(country, aes(Myeloma_5Y_Survival_Pct)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "skyblue", color = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(country$Myeloma_5Y_Survival_Pct, na.rm = TRUE),
                                          sd = sd(country$Myeloma_5Y_Survival_Pct, na.rm = TRUE)),
                color = "red", linewidth = 1) +
  labs(x = "骨髓瘤5年生存率(%)", y = "密度", title = "全球骨髓瘤5年生存率分布") +
  theme_classic()
ggsave(file.path(fig_dir, "histogram_survival.pdf"), p2, width = 8, height = 6)
cat("已保存: histogram_survival.pdf\n")

# 1.4.3 散点图：发病率 vs 生存率
p3 <- ggplot(country, aes(Myeloma_Incidence_Per_100K, Myeloma_5Y_Survival_Pct)) +
  geom_point(aes(color = Continent), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  labs(x = "骨髓瘤发病率(每10万)", y = "5年生存率(%)",
       title = "骨髓瘤发病率与5年生存率的关系") +
  theme_classic()
ggsave(file.path(fig_dir, "scatter_incidence_vs_survival.pdf"), p3, width = 8, height = 6)
cat("已保存: scatter_incidence_vs_survival.pdf\n")

# 1.4.4 散点图：BMT可及性 vs 生存率
p4 <- ggplot(country, aes(BMT_Access_Score, Myeloma_5Y_Survival_Pct)) +
  geom_point(aes(color = Continent), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  labs(x = "BMT可及性评分", y = "骨髓瘤5年生存率(%)",
       title = "BMT可及性与生存率的关系") +
  theme_classic()
ggsave(file.path(fig_dir, "scatter_bmt_vs_survival.pdf"), p4, width = 8, height = 6)
cat("已保存: scatter_bmt_vs_survival.pdf\n")

# 1.4.5 分组条形图：癌症类型发病率
p5 <- ggplot(cancer_types, aes(reorder(Cancer_Type, Incidence_Per_100K_US), Incidence_Per_100K_US)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(x = "癌症类型", y = "美国发病率(每10万)", title = "骨髓相关癌症美国发病率") +
  theme_classic()
ggsave(file.path(fig_dir, "barplot_cancer_types.pdf"), p5, width = 10, height = 6)
cat("已保存: barplot_cancer_types.pdf\n")

# 1.4.6 Q-Q图：校验骨髓瘤发病率正态性
png(file.path(fig_dir, "qq_myeloma_incidence.png"), width = 800, height = 600)
qqnorm(country$Myeloma_Incidence_Per_100K, main = "骨髓瘤发病率Q-Q图")
qqline(country$Myeloma_Incidence_Per_100K, col = "red", lwd = 2)
dev.off()
cat("已保存: qq_myeloma_incidence.png\n")

# 1.4.7 分面饼图：各大陆国家数量
continent_counts <- as.data.frame(table(country$Continent))
names(continent_counts) <- c("Continent", "Count")
p7 <- ggplot(continent_counts, aes(x = "", y = Count, fill = Continent)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Continent, "\n", Count)), position = position_stack(vjust = 0.5)) +
  theme_void() + labs(title = "各大陆国家数量分布")
ggsave(file.path(fig_dir, "pie_continent.pdf"), p7, width = 8, height = 7)
cat("已保存: pie_continent.pdf\n")

# 1.4.8 点线图：各大陆关键指标对比
continent_sum <- country %>%
  group_by(Continent) %>%
  summarise(across(c(Myeloma_Incidence_Per_100K, Myeloma_5Y_Survival_Pct,
                      BMT_Access_Score), mean, na.rm = TRUE))
cat("\n各大陆均数汇总:\n")
print(as.data.frame(continent_sum))

cat("\n===== Part 1 完成 =====\n")
