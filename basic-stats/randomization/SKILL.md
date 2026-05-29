---
name: medical-stat-randomization
description: "R语言医学统计：临床研究随机分组。涵盖简单随机化（完全随机）、分层随机化、区组随机化（block randomization）的R语言实现。TRIGGER when user mentions 随机分组、随机化、临床试验设计、分层随机、区组随机、简单随机，or asks about generating randomization schemes for clinical trials. SKIP for 抽样方法、随机抽样、样本量计算。"
---

# 随机分组 (Randomization)

> 本章内容改编自 阿越就是我 的《R语言实战医学统计》（https://github.com/ayueme/R_medical_stat），采用 CC BY-SA 4.0 许可证发布。

## 适用场景

**何时使用随机分组：**

- 临床试验设计阶段，需要将受试者随机分配到试验组与对照组
- 动物实验或基础研究中需要将实验对象随机分组
- 多组比较研究（如：试验组1、试验组2、阳性对照组、阴性对照组）
- 需要保证组间基线特征均衡可比的研究

**何时不用随机分组：**

| 情况 | 替代方法 |
|------|----------|
| 从已有总体中抽取代表性样本 | 抽样方法（simple sampling / stratified sampling） |
| 计算研究所需样本量 | 样本量计算 / 功效分析（`pwr` 包） |
| 观察性研究中的组间均衡 | 倾向性评分匹配（PSM） |
| 事后分析组间差异 | t检验 / 卡方检验 / 方差分析等 |

**本章涵盖的随机分组方法：**

- 简单随机（simple randomization / 完全随机）— 最基础的随机分组
- 区组随机（block randomization）— 临床研究最常用的随机分组方法
- 分层随机（stratified randomization）— 根据重要预后因素分层后再随机分组

## 前置条件

**R 包依赖：**

```r
# 基础R（无需安装）
# sample() —— 简单随机化
# runif() / rank() —— 自定义随机化

# 需安装的扩展包
install.packages("randomizr")  # simple_ra(), complete_ra(), block_ra()
install.packages("blockrand")  # blockrand() — 适合临床研究的区组随机 + PDF信封
install.packages("showtext")   # plotblockrand() 中文字体支持
```

**设计前提：**

1. **随机化方案应在研究开始前确定**，并妥善保存（随机信封/中心随机系统）
2. **设置随机种子**（`set.seed()`）确保结果可复现
3. **区组长度应保密**，避免研究者预测分组（必要时使用随机区组长度）
4. **分层因子**应选取真正影响结局的重要预后因素（如：性别、年龄组、疾病严重程度）
5. **分配比例**通常为 1:1，也可按需要设计其他比例（如 2:1）

## 方法选择决策树

```
你的研究设计 →
├── 无重要预后因素需控制，样本量较大
│   └── 使用【简单随机化】
│       - R 内置：sample()
│       - randomizr：simple_ra() 可不等例 / complete_ra() 等比例
│
├── 受试者逐个入组，需保证组间人数均衡
│   └── 使用【区组随机化】
│       - 临床研究首选：blockrand::blockrand()（支持逐个入组 + PDF信封）
│       - 动物/基础实验：randomizr::block_ra()（一次性分组）
│
└── 有已知重要预后因素（性别、年龄等），需分层控制
    └── 使用【分层随机化】
        - blockrand::blockrand() 配合 stratum 参数（各层独立区组随机）
```

## 标准工作流

### 步骤1：确定随机化方案类型

根据研究设计、入组方式、是否有重要协变量，选择简单随机 / 区组随机 / 分层随机。

### 步骤2：设置随机种子并生成分组

```r
set.seed(20260412)  # 确保可复现
```

不同方法的核心函数调用见下方代码示例。

### 步骤3：核查分组结果

使用 `table()` 检查各组人数是否符合预期比例。

### 步骤4：保存随机化方案

- 将分组结果保存为 CSV / Excel 文件，供研究记录使用
- 对于临床研究，可利用 `blockrand::plotblockrand()` 生成可打印的随机信封 PDF，实现分组隐匿（allocation concealment）

### 步骤5：论文中的报告

在论文的方法部分报告随机化方法，例如：

> 本研究采用区组随机化方法，利用R语言blockrand包生成随机分组序列，区组长度为4，按1:1分配至试验组与对照组。

## 代码示例

### 11.1 简单随机化

**方法一：R 内置 sample()**

30名受试者，按完全随机化分入2组，每组15人：

```r
set.seed(123)
id <- 1:30
group <- sample(rep(c("试验组", "对照组"), 15), 30, replace = FALSE)
rand_tbl <- data.frame(ID = id, Group = group)
table(rand_tbl$Group)
## 
## 对照组 试验组 
##     15     15
```

**方法二：randomizr::simple_ra() — 抛硬币式随机**

```r
library(randomizr)

# 100人按抛硬币方式分2组（组间人数可能不等）
sim <- simple_ra(100, num_arms = 2, conditions = c("试验组", "对照组"))
table(sim)
## sim
## 试验组 对照组 
##     52     48
```

**方法三：randomizr::complete_ra() — 等比例随机**

```r
# 确保组间人数相等
com <- complete_ra(100, num_arms = 2, conditions = c("试验组", "对照组"))
table(com)
## com
## 试验组 对照组 
##     50     50
```

**方法四：自定义函数（支持多比例）**

```r
simple_random <- function(size, grp = 2, T_2_C = "1:1") {
  set.seed(20210412)
  id_num <- seq(1, size, 1)
  random_seq <- runif(n = size, min = 0, max = 1)
  int_rank <- rank(random_seq)
  
  ratio_T <- as.numeric(substr(T_2_C, 1, 1))
  ratio_C <- as.numeric(substr(T_2_C, 3, 3))
  
  if (grp == 2) {
    group <- ifelse(int_rank <= size / (ratio_T + ratio_C), "T", "C")
  } else if (grp >= 3) {
    group <- cut(int_rank, breaks = grp, labels = paste("Group", 1:grp))
  }
  
  df <- data.frame(ID = id_num, RandomNum = random_seq,
                   Rank = int_rank, Group = group)
  return(df)
}

simple_random(20)
##    ID  RandomNum Rank Group
## 1   1 0.83237492   19     C
## 2   2 0.95522177   20     C
## 3   3 0.59787880   10     T
## 4   4 0.35076793    4     T
## 5   5 0.43157421    6     T
## ...
```

### 11.2 区组随机化

**randomizr::block_ra()** — 适用于一次性分组的实验设计：

```r
library(randomizr)

# 以毛发颜色为区组，每个区组内分为3组
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                         times = HairEyeColor$Freq), 1:3]

Z <- block_ra(blocks = hec$Hair,
              conditions = c("Control", "Placebo", "Treatment"))
table(Z, hec$Hair)
##             
## Z            Black Brown Red Blond
##   Control       36    95  24    42
##   Placebo       36    96  23    42
##   Treatment     36    95  24    43
```

**blockrand::blockrand()** — 推荐用于临床研究，支持逐个受试者入组：

```r
library(blockrand)

set.seed(111)
res <- blockrand(n = 100, num.levels = 2,
                 levels = c("试验组", "对照组"))

head(res)
##   id block.id block.size treatment
## 1  1        1          4    试验组
## 2  2        1          4    对照组
## 3  3        1          4    试验组
## 4  4        1          4    对照组
## 5  5        2          6    试验组
## 6  6        2          6    试验组

table(res$treatment)
## 
## 对照组 试验组 
##     51     51
```

注意：`blockrand()` 生成的总人数为 `n + floor(n / min(block.sizes))`，因为每个区组末尾可能多填充1人以保证区组内的比例平衡。如上例 100人实际输出102人。若需严格 50:50，可截取前100行或适当调整入参。

**生成随机信封 PDF（分组隐匿）：**

```r
library(showtext)
showtext_auto(enable = TRUE)

plotblockrand(res, file = "res.pdf",
  top = list(
    text = c("xxx临床研究", "受试者编号：%ID%", "入组：%TREAT%"),
    col = c("black", "blue", "red"),
    font = c(2, 2, 4)
  ),
  middle = list(
    text = c("xxx临床研究", "受试者编号：%ID%"),
    col = c("black", "blue"),
    font = c(2, 2)
  ),
  bottom = "联系电话：123456789",
  cut.marks = TRUE  # 裁剪标记
)
```

%ID%、%TREAT%、%STRAT% 分别映射至 id、treatment、stratum 列。信封用于分组隐匿：按入组顺序依次拆开信封，研究者无法预测下一位受试者的分组。

### 11.3 分层随机化

120名受试者分为4组（试验组1、试验组2、阳性对照组、阴性对照组），每组30人。按性别分层，男女各60例。

```r
library(blockrand)

# 男性60例
set.seed(123)
res.M <- blockrand(
  n = 60,
  num.levels = 4,
  levels = c("试验组1", "试验组2", "阳性对照组", "阴性对照组"),
  stratum = "男性",
  id.prefix = "男",
  block.sizes = c(4)
)
table(res.M$treatment)
## 
##    试验组1    试验组2 阳性对照组 阴性对照组 
##         15         15         15         15

# 女性60例
set.seed(456)
res.F <- blockrand(
  n = 60,
  num.levels = 4,
  levels = c("试验组1", "试验组2", "阳性对照组", "阴性对照组"),
  stratum = "女性",
  id.prefix = "女",
  block.sizes = c(4)
)
table(res.F$treatment)
## 
##    试验组1    试验组2 阳性对照组 阴性对照组 
##         15         15         15         15

# 合并分层结果
res <- rbind(res.M, res.F)
table(res$stratum, res$treatment)
##        
##         试验组1 试验组2 阳性对照组 阴性对照组
##     女性      15      15         15         15
##     男性      15      15         15         15
```

**成分层随机信封 PDF：**

```r
library(showtext)
showtext_auto(enable = TRUE)

plotblockrand(res, file = "res1.pdf",
  top = list(
    text = c("不得了临床试验", "受试者编号: %ID%", "组别: %TREAT%"),
    col = c('black', 'blue', 'red'), font = c(2, 2, 4)
  ),
  middle = list(
    text = c("不得了临床试验", "性别: %STRAT%", "受试者编号: %ID%"),
    col = c('black', 'blue', 'red'), font = c(1, 2, 3)
  ),
  bottom = "联系电话：123456789",
  cut.marks = TRUE
)
```

## 结果解读指南

| 输出 | 含义 |
|------|------|
| `table(res$treatment)` | 各组最终分配人数，核验比例是否与预期一致 |
| `id` | 受试者编号，按入组顺序分配 |
| `block.id` | 区组编号，同一区组内已知比例的分组排列 |
| `block.size` | 各区组包含的受试者人数（建议保密） |
| `stratum` | 分层变量取值，分层方案中各层独立完成区组随机 |

**重要概念：**

- **随机种子（seed）**：保证随机序列可复现，需在研究方案中记录。同一 seed + 同一参数 = 相同分组结果。
- **randomizr vs blockrand**：`block_ra()` 需要提前有全部受试者名单（适合动物实验），`blockrand()` 按顺序逐个生成分组（适合临床研究实时入组）。
- **区组长度**：通常取组数的倍数（如 2 组则区组长 4/6/8）。随机区组长度可降低研究者预测分组的可能性。

## 常见问题与注意事项

**Q1: 简单随机和区组随机的根本区别是什么？**

简单随机每个受试者独立随机分配，可能出现组间人数不平衡（尤其是小样本）。区组随机将受试者分成若干区组，每个区组内各组人数按比例分配，保证整体均衡和每批受试者内部均衡。临床研究几乎都使用区组随机。

**Q2: blockrand() 输出人数为什么多于入参 n？**

`blockrand()` 保证每个区组末尾也满足比例平衡，因此总人数可能略微增加（通常多 1-2 人）。实际使用时取前 n 行即可，或适当调整入参 n。

**Q3: 分层随机和事后分层调整有什么区别？**

分层随机是在随机分组之前按分层因子分别生成随机序列（如：男女分别随机），保证各层内组间均衡。事后分层调整是在分析阶段通过回归模型控制混杂因素，但无法完全消除组间基线不均衡的风险。应在设计阶段（分层随机）优先处理。

**Q4: 区组随机为什么需要保管好区组长度？**

如果研究者知道区组长度，当区组内前面的受试者分组已知时，可推断剩余受试者的分组方向，破坏分配隐匿（allocation concealment）。使用 `block.sizes = c(4, 6, 8)` 等随机长度可降低此风险。

**Q5: 什么时候用 block_ra() 而不是 blockrand()？**

`randomizr::block_ra()` 要求一次性输入所有受试者及其区组信息，适合动物实验或提前知道所有样本的基础研究。`blockrand::blockrand()` 适合临床研究，按入组序号依次给出分组结果，无需提前知道全部受试者信息。

**Q6: 多组随机（3组以上）如何操作？**

`blockrand()` 通过 `num.levels` 和 `levels` 参数支持任意组数。`simple_random()` 自定义函数中当 `grp >= 3` 时使用 `cut()` 等距分为指定组数。`randomizr::complete_ra()` 也支持多组设计。

**Q7: 随机分组和随机抽样的区别？**

随机分组是在已有受试者的基础上，将受试者分配到不同治疗组。随机抽样是从目标人群中抽取代表性样本。本章覆盖随机分组，如需抽样方法请参考相关统计教材。
