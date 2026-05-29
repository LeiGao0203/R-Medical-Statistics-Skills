---
name: medical-stat-stat-plot
description: "R语言医学统计：统计绘图基础。涵盖直方图、箱线图、散点图、条形图、误差线图、茎叶图、P-P图/Q-Q图的R语言实现。TRIGGER when user mentions 统计绘图、箱线图、直方图、散点图、R绘图、ggplot2、base plot，or asks about creating statistical graphs in R for medical data. SKIP for ROC曲线、森林图、生存曲线、三线表。"
---

# 统计绘图 (Statistical Graphics)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用统计绘图：**

- 数据探索阶段，了解数据的分布特征（直方图、箱线图、茎叶图）
- 展示分类变量的频数或构成比（条形图、饼图、百分比条形图）
- 展示两个连续变量之间的关系（散点图）
- 展示随时间变化的趋势（折线图、点线图）
- 比较多组的均数及其变异范围（箱线图、误差条图）
- 检查数据是否服从正态分布（P-P图、Q-Q图）

**何时不使用这些图：**

| 情况 | 替代方法 |
|------|----------|
| 诊断试验评价 | ROC曲线 |
| 多因素回归结果展示 | 森林图 |
| 生存分析 | Kaplan-Meier曲线 |
| 论文表格 | 三线表 |
| 高维数据降维可视化 | PCA图、热图 |

## 前置条件

**R 包依赖：**

```r
# 基础绘图（R自带，无需安装）
# graphics  —— hist(), boxplot(), plot(), barplot(), pie(), stem()
# stats     —— qqnorm(), qqplot(), qqline()

# ggplot2 生态
install.packages("ggplot2")
install.packages("patchwork") # 拼图

# 数据处理
install.packages("dplyr")
install.packages("haven")     # 读取SPSS数据
install.packages("scales")    # 百分比格式化
```

**数据格式要求：**

| 图形类型 | 数据格式 |
|---------|---------|
| 直方图 | 一个连续型数值变量（原始数据，非汇总计数） |
| 箱线图 | 一个分类变量（分组）+ 一个连续型变量 |
| 条形图 | 分类变量 + 数值变量（频数/率/均值） |
| 分组条形图 | 两个分类变量 + 数值变量 |
| 饼图 | 分类变量 + 构成数量 |
| 折线图/点线图 | x轴（时间/有序）+ y轴（数值）+ 可选分组变量 |
| 散点图 | 两个连续型数值变量 |
| 茎叶图 | 一个连续型数值变量（原始数据） |
| 误差条图 | 分组变量 + 均值 + 置信区间上下界 |
| Q-Q图 | 一个连续型数值变量（原始数据） |

**注意事项：**

- 直方图要求原始数据，若只有汇总计数需用条形图伪装（手动调width和position）
- 饼图在ggplot2中通过`coord_polar()`将堆叠条形图转为极坐标实现
- 箱线图自动标记离群值（超出1.5倍IQR的点）

## 方法选择决策树

```
你的数据情况 →
├── 单个连续变量，了解分布形态
│   ├── 样本量较小（< 100） → 茎叶图 stem()
│   ├── 样本量较大 → 直方图 hist() / geom_histogram()
│   ├── 检查正态性 → Q-Q图 qqnorm()
│   └── 查看五数概括及离群值 → 箱线图 boxplot()
│
├── 单个分类变量，展示频数或构成比
│   ├── 类别数 ≤ 6，展示构成比 → 饼图 pie() / coord_polar()
│   └── 类别数较多，展示频数 → 条形图 barplot()
│
├── 一个分类变量 × 一个连续变量（分组比较）
│   ├── 比较分布 → 箱线图 geom_boxplot()
│   └── 展示均值 ± CI → 误差条图 geom_errorbar() + geom_point()
│
├── 两个分类变量 × 一个数值变量
│   ├── 展示绝对数值 → 分组条形图 (position="dodge")
│   └── 展示相对构成 → 百分比堆叠条形图 (position="fill")
│
├── 两个连续变量 → 散点图 plot() / geom_point()
│
└── 时间序列数据
    ├── 单条线 → 折线图 geom_line()
    ├── 多条线 + 数据点 → 点线图 geom_line() + geom_point()
    └── log变换对比 → 双面板 (patchwork拼图)
```

## 标准工作流

### 步骤1：数据准备与探索

```r
library(ggplot2)
library(dplyr)

str(mydata)
head(mydata)
summary(mydata)
```

### 步骤2：选择绘图系统 → 步骤3：调整参数 → 步骤4：导出图形

- **Base R**：快速探索（`hist()`、`boxplot()`、`stem()`、`plot()`）
- **ggplot2**：出版级质量，图层叠加，适合论文发表

```r
# base R
png("output.png", width = 800, height = 600, res = 150)
hist(x)
dev.off()

# ggplot2
ggsave("output.png", width = 8, height = 6, dpi = 300)
```

### 步骤5：论文中的统计图描述

论文中应说明：图形类型、x/y轴含义、各组颜色/填充含义、主要发现。

## 代码示例

### 直方图

连续变量分布形态展示。**base R**：`hist()`；**ggplot2**：`geom_histogram()`。

```r
# Base R —— 模拟100例血红蛋白数据
set.seed(123)
hb <- rnorm(100, mean = 140, sd = 10)

hist(hb, breaks = 20, freq = FALSE, col = "skyblue", border = "white",
     main = "血红蛋白分布", xlab = "血红蛋白 (g/L)")
curve(dnorm(x, mean = mean(hb), sd = sd(hb)),
       col = "red", lwd = 2, add = TRUE)

# ggplot2
df_hb <- data.frame(hb = hb)
ggplot(df_hb, aes(x = hb)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20, fill = "skyblue", color = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(hb), sd = sd(hb)),
                color = "red", linewidth = 1) +
  labs(x = "血红蛋白 (g/L)", y = "密度") +
  theme_classic()
```

若只有汇总计数数据（如课本例10-10），用条形图伪装直方图：

```r
library(haven)
data10_10 <- haven::read_sav("datasets/例10-10.sav", encoding = "GBK")
data10_10 <- as_factor(data10_10)

ggplot(data10_10, aes(age, count)) +
  geom_bar(stat = "identity", fill = "white", color = "black",
           width = 1, position = position_dodge(width = 1)) +
  labs(x = "年龄（岁）", y = "每岁病例数") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.title = element_text(color = "black", size = 15))
```

### 箱线图

比较多个组的连续变量分布，显示中位数、四分位数及离群值。

```r
# Base R
boxplot(Sepal.Length ~ Species, data = iris)

# ggplot2
ggplot(iris, aes(Species, Sepal.Length)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  theme_classic()
```

### 条形图与分组条形图

**基础条形图**（例10-4）与**分组条形图**（例10-5）：

```r
# Base R 条形图
rate <- c(17.9, 20.8, 33.3)
names(rate) <- c("方法1", "方法2", "方法3")
barplot(rate, col = "white", border = "black", ylim = c(0,40), ylab = "再发率（%）")

# ggplot2 条形图
data10_4 <- data.frame(`灌注方法` = c("方法1","方法2","方法3"),
                        rate = c(17.9,20.8,33.3), check.names = FALSE)
ggplot(data10_4, aes(`灌注方法`, rate)) +
  geom_bar(stat = "identity", fill = "white", color = "black", width = 0.4) +
  ylab("再发率（%）") + scale_y_continuous(expand = c(0,0)) + theme_classic()

# 分组条形图（position="dodge"）
data10_5 <- haven::read_sav("datasets/例10-05.sav"); data10_5 <- as_factor(data10_5)
ggplot(data10_5, aes(agent, rate)) +
  geom_bar(stat = "identity", aes(fill = year), position = "dodge", color = "black") +
  labs(x = "性别", y = "患龋率（%）", fill = "年份") +
  scale_y_continuous(expand = c(0,0)) + theme_classic()
```

### 折线图与点线图

展示随时间变化的趋势。

```r
# 折线图（例10-8 布氏菌病发病人数）
data10_8 <- haven::read_sav("datasets/例10-08.sav", encoding = "GBK")
data10_8 <- as_factor(data10_8)

ggplot(data10_8, aes(year, counts)) +
  geom_line(aes(group = agent, linetype = agent)) +
  labs(x = "年份", y = "布氏菌病发病人数", linetype = "性别") +
  theme_classic()

# 点线图 + log变换双面板（例10-9）
data10_9 <- haven::read_sav("datasets/例10-09.sav", encoding = "GBK")
data10_9 <- as_factor(data10_9)

p1 <- ggplot(data10_9, aes(year, `发病率`)) +
  geom_line(aes(group = `病型`, linetype = `病型`)) +
  geom_point(aes(group = `病型`, shape = `病型`), size = 4) + theme_classic()

p2 <- ggplot(data10_9, aes(year, log10(`发病率`))) +
  geom_line(aes(group = `病型`, linetype = `病型`)) +
  geom_point(aes(group = `病型`, shape = `病型`), size = 4) + theme_classic()

library(patchwork)
p1 + p2 + plot_layout(guides = "collect")
```

### 散点图

两个连续变量的关系探索。

```r
plot(mtcars$wt, mtcars$mpg, pch = 19, col = "steelblue",
     xlab = "车重", ylab = "mpg")
abline(lm(mpg ~ wt, data = mtcars), col = "red", lwd = 2)

# ggplot2
ggplot(mtcars, aes(wt, mpg)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_classic()
```

### 饼图

展示分类变量的构成比。Base R直接`pie()`；ggplot2通过`coord_polar()`实现。

```r
# Base R
counts <- c(226, 52, 22, 17, 10)
names(counts) <- c("无菌性松动", "感染", "假体周围骨折", "假体不稳定", "其他")
pie(counts)

# ggplot2（例10-6）
data10_6 <- data.frame(
  `失败原因` = c("无菌性松动", "感染", "假体周围骨折", "假体不稳定", "其他"),
  `数量` = c(226, 52, 22, 17, 10), check.names = FALSE
) %>%
  arrange(desc(`数量`)) %>%
  mutate(prop = scales::percent(round(`数量` / sum(`数量`), 2)))

ggplot(data10_6, aes(x = "", y = `数量`, fill = `失败原因`)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(label = prop), position = position_stack(vjust = 0.5)) +
  coord_polar("y", start = 0) + theme_void()
```

### 百分比堆叠条形图

展示多组分类数据的相对构成（例10-7）。

```r
library(scales)
data10_7 <- haven::read_sav("datasets/例10-07.sav", encoding = "GBK")
data10_7 <- as_factor(data10_7)

ggplot(data10_7, aes(year, percent, fill = reason)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_bw() + theme(legend.position = "bottom")
```

### 茎叶图

保留原始数值每一位信息，适合小样本快速探索。

```r
data10_13 <- haven::read_sav("datasets/例10-13.sav", encoding = "GBK")
data10_13 <- as_factor(data10_13)

stem(data10_13$rbc, scale = 1)
##   The decimal point is 1 digit(s) to the left of the |
##   30 | 7
##   32 | 17
##   35 | 299
##   36 | 0124467789
##   37 | 12266679
##   ...
##   52 | 348
##   54 | 6
```
左侧"茎"为整数部分×0.1，"叶"为下一位。如`30 | 7`表示3.07。

### 误差条图

展示各组的均数和95%置信区间（例10-14）。

```r
library(dplyr)
data10_14 <- haven::read_sav("datasets/例10-14.sav", encoding = "GBK")
data10_14 <- as_factor(data10_14)

data10_14_summary <- data10_14 %>%
  group_by(group) %>%
  summarise(
    mm = mean(dmdz),
    lower = mm - 1.96 * (sd(dmdz) / sqrt(n())),
    upper = mm + 1.96 * (sd(dmdz) / sqrt(n()))
  )

ggplot(data10_14_summary) +
  geom_point(aes(group, mm), size = 4, shape = 0) +
  geom_errorbar(aes(x = group, ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "分组", y = "95%CI") +
  theme_classic() +
  theme(axis.title = element_text(color = "black", size = 15))
```

### Q-Q图（正态性诊断）

```r
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

# Base R
qqnorm(x, main = "正态Q-Q图")
qqline(x, col = "red", lwd = 2)

# ggplot2
ggplot(data.frame(x), aes(sample = x)) +
  stat_qq() + stat_qq_line(color = "red") +
  theme_classic()
```

## 结果解读指南

**直方图**：对称钟形→近似正态；左偏→均值<中位数（如住院天数）；右偏→均值>中位数（如医疗费用）；双峰→可能有两个亚群。

**箱线图**：粗线=中位数，箱体上下边=Q1/Q3，箱体高度=IQR；触须=非离群值范围（±1.5×IQR）；触须外=离群值。箱体位置和中位数差异提示组间差异。

**Q-Q图**：点沿对角线→正态；两端偏离→尾部偏厚/薄；S形偏离→偏态分布。

**误差条图**：误差条=95%CI；两组不重叠→差异可能有统计学意义；多组梯度→剂量-反应关系。

**茎叶图**：保留原始数据每一位信息，适合检查精度，较少用于正式发表。

## 常见问题与注意事项

**Q：直方图组距怎么选？**
A：组距过小图形崎岖，过大丢失信息。`hist()`默认使用Sturges公式自动选择，可尝试不同`breaks`值。

**Q：箱线图 vs 直方图 vs 误差条图如何选？**
A：箱线图展示分位数和离群值，适合多组并行；直方图展示完整分布，适合单变量探索；误差条图只展示均值±CI/SE，信息量最少但常见于论文。

**Q：饼图 vs 条形图？**
A：类别≤5且关注构成比时可用饼图；类别多时用条形图。多数统计学家推荐条形图，因人眼难以精确比较扇形面积。

**Q：ggplot2中让y轴从0开始？**
A：`scale_y_continuous(expand = c(0, 0))`，条形图几乎总是需要。

**Q：stem()的scale参数？**
A：`scale=1`默认；`scale=2`展开；`scale=0.5`压缩。输出过于稀疏或密集时调整。

**Q：如何导出发表级图形？**
A：`ggsave("figure.pdf", width=8, height=6)`输出矢量图（无限分辨率）；`ggsave("figure.tiff", width=8, height=6, dpi=300)`输出位图，满足期刊300 dpi要求。

**Q：patchwork拼图用法？**
A：`p1 + p2`即可拼图；`p1 + p2 + plot_layout(guides = "collect")`合并图例；`(p1 | p2) / p3`控制行列布局。
