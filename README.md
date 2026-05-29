# R Medical Statistics Skills

[English version](README.en.md)

一套面向医学统计和 R 语言分析场景的 Codex skills。项目包含基础统计、高级统计、医学文献常用统计方法，以及用于生成普通 R 脚本、统计报告和维护 Jupyter Notebook 的辅助 skill。

## Contents

- `basic-stats/`: t 检验、方差分析、卡方检验、相关分析、ROC、样本量、统计绘图等基础医学统计 skill。
- `advanced-stats/`: 协方差分析、多元回归、Logistic 回归、生存分析、PCA、结构方程、多水平模型等进阶统计 skill。
- `literature-stats/`: 倾向评分、Fine-Gray、限制性立方样条、亚组分析、趋势检验等医学文献常见方法 skill。
- `r-script/`: 原创 R 脚本 skill，用于为 RStudio、命令行 R 或非 notebook 用户生成可复现 `.R` 分析脚本。
- `quarto-report/`: 原创 Quarto/R Markdown 报告 skill，用于生成可导出 HTML、Word 或 PDF 的医学统计分析报告。
- `jupyter-notebook/`: 原创 Jupyter Notebook skill，用于创建、整理和验证可复现 notebook。
- `example/`: 示例目录，后续可放入演示数据、notebook 或使用案例。

## 推荐使用方式

- **普通 R 脚本**：不使用 Jupyter 的用户优先使用 `r-script/`，输出 `analysis.R`，可在 RStudio 中 Source，或用 `Rscript analysis.R` 在命令行运行。
- **Jupyter Notebook**：适合交互式探索、教学演示、逐步解释和需要 `.ipynb` 交付的任务。
- **Quarto / R Markdown 报告**：使用 `quarto-report/`，适合正式报告、论文附录、课题汇报，以及可导出的 HTML、Word、PDF。

## Install

将各个 skill 目录复制到本机 Codex skills 目录：

```bash
mkdir -p ~/.codex/skills
find advanced-stats basic-stats literature-stats -mindepth 1 -maxdepth 1 -type d -exec cp -R {} ~/.codex/skills/ \;
cp -R r-script ~/.codex/skills/
cp -R quarto-report ~/.codex/skills/
cp -R jupyter-notebook ~/.codex/skills/
```

重启 Codex 后，相关 skill 会在对应医学统计、R 脚本、统计报告或 notebook 任务中触发。

## License

本项目采用混合许可证：

- `basic-stats/`、`advanced-stats/`、`literature-stats/` 中与《R语言实战医学统计》相关的内容，改编自阿越就是我的开源项目 [R_medical_stat](https://github.com/ayueme/R_medical_stat)，按 CC BY-SA 4.0 发布。
- `r-script/`、`quarto-report/` 和 `jupyter-notebook/` 为原创内容，按 Apache License 2.0 发布。
- `example/` 中未来新增内容请在对应文件或目录中单独注明许可证。

详见 [LICENSE](LICENSE)。

## Attribution

统计类 skill 的部分内容基于《R语言实战医学统计》整理和改写。再次分发、修改或演绎相关内容时，请保留原作者署名和 CC BY-SA 4.0 授权信息。

## Contributing

欢迎补充新的统计方法、修正代码示例、改进方法选择逻辑或增加可复现实例。提交前请阅读 [CONTRIBUTING.md](CONTRIBUTING.md)。
