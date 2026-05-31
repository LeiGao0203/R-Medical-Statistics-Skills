# ============================================================================
# Part 5: 主成分分析、聚类分析、生存分析
# Skills: medical-stat-pca, medical-stat-pca-vis, medical-stat-cluster,
#         medical-stat-survival, medical-stat-survival-vis
# Datasets: bone_marrow_cancer_by_country.csv, bone_marrow_cancer_survival.csv,
#           bone_marrow_cancer_trends_2000_2026.csv
# ============================================================================

library(factoextra)
library(FactoMineR)
library(psych)
library(performance)
library(parameters)
library(cluster)
library(NbClust)
library(flexclust)
library(survival)
library(survminer)
library(ggplot2)

data_dir <- "data"
fig_dir <- "analysis/figures"

country <- read.csv(file.path(data_dir, "bone_marrow_cancer_by_country.csv"))
survival <- read.csv(file.path(data_dir, "bone_marrow_cancer_survival.csv"))
trends <- read.csv(file.path(data_dir, "bone_marrow_cancer_trends_2000_2026.csv"))

# ---- 5.1 主成分分析 (PCA) ----
cat("\n========== 5.1 主成分分析 ==========\n")

# 选择连续变量
pca_vars <- c("Myeloma_Incidence_Per_100K", "Leukemia_Incidence_Per_100K",
              "Myeloma_5Y_Survival_Pct", "Leukemia_5Y_Survival_Pct",
              "BMT_Access_Score", "BMT_Centers", "Hematologists_Per_Million")
country_pca <- country[, pca_vars]
rownames(country_pca) <- country$Country

# 前提条件检验
cat("\n--- KMO 和 Bartlett 球形检验 ---\n")
check_factorstructure(country_pca)

# 执行PCA
pca_res <- PCA(country_pca, scale.unit = TRUE, graph = FALSE)

# 特征值与方差贡献率
cat("\n--- 特征值与方差贡献率 ---\n")
eig <- get_eigenvalue(pca_res)
print(eig)

# 碎石图
p_scree <- fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 80))
ggsave(file.path(fig_dir, "scree_plot.pdf"), p_scree, width = 8, height = 6)
cat("已保存: scree_plot.pdf\n")

# 变量载荷图
p_var <- fviz_pca_var(pca_res, col.var = "cos2",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       repel = TRUE)
ggsave(file.path(fig_dir, "pca_variables.pdf"), p_var, width = 8, height = 7)
cat("已保存: pca_variables.pdf\n")

# 变量贡献图
p_contrib <- fviz_contrib(pca_res, choice = "var", axes = 1:2)
ggsave(file.path(fig_dir, "pca_contributions.pdf"), p_contrib, width = 10, height = 6)
cat("已保存: pca_contributions.pdf\n")

# 样本得分图（按大陆着色）
p_ind <- fviz_pca_ind(pca_res, geom.ind = "point",
                       col.ind = country$Continent,
                       palette = "jco",
                       addEllipses = TRUE,
                       legend.title = "Continent",
                       title = "PCA: 各国骨髓瘤/白血病卫生指标空间分布")
ggsave(file.path(fig_dir, "pca_individuals.pdf"), p_ind, width = 9, height = 7)
cat("已保存: pca_individuals.pdf\n")

# 双标图
p_biplot <- fviz_pca_biplot(pca_res, col.ind = country$Continent,
                             palette = "jco", addEllipses = TRUE,
                             label = "var", col.var = "black",
                             repel = TRUE, legend.title = "Continent",
                             title = "PCA Biplot: 国家-变量双标图")
ggsave(file.path(fig_dir, "pca_biplot.pdf"), p_biplot, width = 10, height = 8)
cat("已保存: pca_biplot.pdf\n")

# 维度描述
cat("\n--- 维度描述 (dimdesc) ---\n")
dimdesc_res <- dimdesc(pca_res, axes = c(1, 2), proba = 0.05)
print(dimdesc_res$Dim.1)
print(dimdesc_res$Dim.2)

# ---- 5.2 聚类分析 ----
cat("\n========== 5.2 聚类分析 ==========\n")

# 标准化
country_scaled <- scale(country_pca)

# 确定最佳聚类数
set.seed(123)
nc <- NbClust(country_scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
              method = "average")
cat("\n--- NbClust 最佳聚类数 ---\n")
print(table(nc$Best.nc[1,]))

# 最优K值
best_k <- as.integer(names(which.max(table(nc$Best.nc[1,]))))
cat(sprintf("\n最优聚类数: k = %d\n", best_k))

# 层次聚类
h_clust <- hclust(dist(country_scaled, method = "euclidean"), method = "average")

# 树状图
png(file.path(fig_dir, "dendrogram.png"), width = 1200, height = 700)
plot(h_clust, hang = -1, main = "各国骨髓瘤卫生指标层次聚类",
     xlab = "", sub = "", cex = 0.7)
rect.hclust(h_clust, k = best_k, border = 2:6)
dev.off()
cat("已保存: dendrogram.png\n")

# 切分聚类
clusters <- cutree(h_clust, k = best_k)
country$Cluster <- factor(clusters)

cat("\n--- 各类别国家数量 ---\n")
print(table(country$Cluster))

# K-means聚类
set.seed(123)
km_fit <- kmeans(country_scaled, centers = best_k, nstart = 25)
cat(sprintf("K-means 组间SS/总SS = %.1f%%\n",
            km_fit$betweenss / km_fit$totss * 100))

# 可视化K-means聚类
p_km <- fviz_cluster(km_fit, data = country_scaled,
                      ellipse = TRUE, ellipse.type = "t",
                      geom = "point", palette = "lancet",
                      ggtheme = theme_bw(),
                      main = "K-means聚类: 各国卫生指标分组")
ggsave(file.path(fig_dir, "kmeans_cluster.pdf"), p_km, width = 8, height = 7)
cat("已保存: kmeans_cluster.pdf\n")

# 各类别变量均值
cat("\n--- 各类别变量均值 (原始尺度) ---\n")
for (i in 1:best_k) {
  cat(sprintf("\n类别 %d (n=%d):\n", i, sum(clusters == i)))
  print(colMeans(country_pca[clusters == i, ]))
}

# 各类别大陆分布
cat("\n--- 各类别大陆分布 ---\n")
print(table(country$Cluster, country$Continent))

# ---- 5.3 生存分析 ----
cat("\n========== 5.3 生存分析 ==========\n")

# 使用生存率数据构建模拟生存数据
# 利用trends数据做时间趋势的"生存"分析

# 构建全球骨髓瘤累计生存
cat("\n--- 使用趋势数据模拟生存分析 ---\n")

# 创建模拟生存数据（利用发病率、死亡率反推）
# Global_Myeloma_New_Cases vs Global_Myeloma_Deaths 可计算case fatality ratio

trends$myeloma_mortality_rate <- trends$Global_Myeloma_Deaths / trends$Global_Myeloma_New_Cases
trends$leukemia_mortality_rate <- trends$Global_Leukemia_Deaths / trends$Global_Leukemia_New_Cases

cat("\n历年骨髓瘤病死率变化:\n")
for (yr in c(2000, 2005, 2010, 2015, 2020, 2026)) {
  row <- which(trends$Year == yr)
  if (length(row) > 0) {
    cat(sprintf("  年份 %d: 病死率 = %.2f, 5年生存率 = %.1f%%\n",
                yr, trends$myeloma_mortality_rate[row], trends$Myeloma_5Y_Survival_Pct[row]))
  }
}

# 构建模拟生存对象（5年生存率 > 50% 为"存活"）
# 用country数据中的5年生存率构造
country$surv_time <- 60  # 假设随访60个月
country$surv_status <- rbinom(nrow(country), 1,
                               prob = country$Myeloma_5Y_Survival_Pct / 100)

# KM生存曲线（按大陆分组）
fit_km <- survfit(Surv(surv_time, surv_status) ~ Continent, data = country)

# log-rank检验
survdiff_result <- survdiff(Surv(surv_time, surv_status) ~ Continent, data = country)
cat("\n--- Log-rank 检验（大陆）---\n")
print(survdiff_result)

# Cox回归
cox_fit <- coxph(Surv(surv_time, surv_status) ~ BMT_Access_Score + Hematologists_Per_Million,
                 data = country)
cat("\n--- Cox回归 ---\n")
cox_summary <- summary(cox_fit)
print(cox_summary)

# Cox回归整洁结果
library(broom)
cox_tidy <- tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE)
cat("\nHR及其95% CI:\n")
print(cox_tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value")])

# ---- 5.4 生存曲线可视化 ----
cat("\n========== 5.4 生存曲线可视化 ==========\n")

p_km_plot <- ggsurvplot(fit_km, data = country,
                          pval = TRUE, conf.int = TRUE,
                          surv.median.line = "hv",
                          risk.table = TRUE,
                          palette = "jco",
                          legend.title = "大陆",
                          ggtheme = theme_classic2(),
                          title = "模拟Kaplan-Meier生存曲线 (按大陆分组)")
ggsave(file.path(fig_dir, "km_curves.pdf"), plot = print(p_km_plot),
       width = 10, height = 8)
cat("已保存: km_curves.pdf\n")

# 森林图
png(file.path(fig_dir, "forest_plot.png"), width = 800, height = 600)
ggforest(cox_fit, data = country, main = "骨髓瘤生存风险比(HR)",
         fontsize = 0.8, noDigits = 2)
dev.off()
cat("已保存: forest_plot.png\n")

cat("\n===== Part 5 完成 =====\n")
