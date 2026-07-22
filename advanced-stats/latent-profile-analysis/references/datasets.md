# 对话中提到的数据集：获取与 LPA 适配

| 数据集 | 公开性 | 适合的连续指标 | 访问/分析注意点 |
|---|---|---|---|
| `tidyLPA::pisaUSA15` | R 包内置示例 | 兴趣、享受、自我效能 | 非医学数据，最适合先验证代码 |
| NHANES | CDC 公开下载 | BMI、腰围、血压、HbA1c、脂质、炎症指标 | 多阶段复杂抽样；按周期和文件使用权重、PSU、分层 |
| CHARLS | 需注册/按项目申请 | 心理健康、认知、握力、肺功能、生理测量 | 波次、抽样权重、量表计分和缺失机制需核对 |
| WHO SAGE | 需按官方条件获取 | 功能、活动、健康相关生活质量 | 跨国合并前检查测量可比性和国家设计 |
| HRS | 公共部分需注册，部分数据受限 | 认知、功能、心理社会资源 | 纵向权重、死亡/功能结局和代理回答需处理 |
| CLHLS | 申请/注册后获取 | 认知、ADL、心理和健康长寿指标 | 高龄、代理回答、选择性生存、地板/天花板效应 |
| CFPS | 需注册/申请 | 心理健康、社会支持、健康行为 | 更偏社会医学；家庭/个人/社区层级不能混为独立行 |
| MIMIC-IV | 完成培训和数据使用协议后获取 | 入院生理、实验室、器官功能摘要 | 先固定患者/住院和时间窗；重复观测不能直接当独立样本 |
| ADNI | 申请/数据使用协议 | 认知、影像、淀粉样蛋白、Tau、生物标志物 | 多中心、纵向、指标重复和测量不变性需重点评估 |
| UK Biobank | 获批后访问，部分使用收费 | 代谢、生活方式、体能、影像、生物标志物 | 极大样本可放大微小差异；不要无理论地投入大量指标 |

当前仓库实际下载并分析的是 NHANES 2017–2018 的六个公开 XPT 文件，其中当前示例使用 DEMO、BMX 和 BPX 文件完成三指标 LPA。其来源和变量见 `example/nhanes-lpa/metadata/`。其他数据集不能在没有账户、协议或研究审批的情况下由脚本自动下载；模块应提供获取说明，而不是绕过访问限制。

## 推荐下载入口

- [NHANES 2017–2018 data portal](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017)
- [NHANES 2017–2018 overview](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overview.aspx?BeginYear=2017)
- [NHANES laboratory overview](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overviewlab.aspx?BeginYear=2017)
- [NHANES analytic guidance](https://wwwn.cdc.gov/nchs/nhanes/analyticguidelines.aspx)
- [tidyLPA documentation](https://data-edu.github.io/tidyLPA/)
