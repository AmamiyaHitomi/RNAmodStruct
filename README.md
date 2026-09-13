# RNAmodStruct (ModStruct)

> A reproducible computational framework for studying the relationship between RNA modification levels and local RNA structure.
>
> 用可复现的计算流程研究 RNA 修饰水平与局部 RNA 结构之间的关联、预测价值和证据边界。

![Study overview](results/figures/10_data_overview_figure1.png)

## 项目简介

RNAmodStruct 以公开人类转录组数据为基础，围绕一个核心问题展开：**RNA 修饰水平是否与局部 RNA 结构信号稳定相关，以及结构信息能否提供额外的预测价值？**

主线研究以 HEK293T 为发现集、HeLa 为外部复核集，结合 GLORI m6A 定量、icSHAPE 实验结构信号和 ViennaRNA 预测结构。后续阶段进一步评估条件依赖性、跨技术一致性、多修饰扩展以及扰动数据中的方向性证据。

本仓库强调：

- 分阶段冻结配置，避免分析过程中隐式改变口径；
- 区分发现、外部复核、敏感性分析和事后探索；
- 使用基因与相同序列联合分组，降低数据泄漏风险；
- 同时保留阳性和阴性结果；
- 通过 SHA-256 manifest、运行状态表和自动化测试保证可追溯性；
- 严格区分关联、方向性、因果和机制证据。

## 当前状态

| 项目 | 状态 |
|---|---|
| 主线 HEK293T 分析与 HeLa 复核 | 完成 |
| Phase 1 方法扩展（12–17） | `REPRODUCIBLE` |
| Phase 2 条件与跨技术验证（18–22） | `PASS` |
| Phase 3 多修饰扩展（23–27） | `COMPLETE` |
| Phase 4 扰动方向性分析（28–31） | `COMPLETE_DIRECTIONALITY_ONLY` |
| 自动化测试 | 87 passed |
| 总发布清单 | 308 个文件通过 SHA-256 校验 |

状态文件位于 [`results/status/`](results/status/)，阶段总结位于 [`results/reports/`](results/reports/)。

## 主要结论

1. **m6A 与局部结构存在幅度较小的稳定正相关。** HEK293T 与 HeLa 的效应方向一致；两细胞系固定效应汇总为 β/SD = 0.01007（95% CI 0.00713–0.01301）。
2. **实验结构信号没有表现出稳定的额外外部预测增益。** 扩大结构窗口、加入条件依赖特征或使用非线性模型，都没有把该结论提升为稳定的机制性证据。
3. **多修饰扩展没有检出明确关联。** 在当前公开数据、技术分层和覆盖阈值下，m5C、m7G 与 Nm 的主要 HeLa 检验均未通过 FDR 0.05；这不等于证明真实效应严格为零。
4. **METTL3 抑制数据支持方向性描述，但不支持因果归因。** GSE264642 中 24,898 个配对转录本呈现系统性结构位移；由于缺少位点级 m6A 定量、救援或失活化合物对照及独立复核，目前最高证据等级仍是 `Directionality`。

详细统计结果和限制请以阶段报告为准：

- [Phase 1 方法扩展总结](results/reports/phase1_summary_report.md)
- [Phase 2 完成与退出门槛](results/reports/22_phase2_completion_report.md)
- [Phase 3 多修饰分析总结](results/reports/27_phase3_completion_report.md)
- [Phase 4 方向性分析总结](results/reports/31_phase4_completion_report.md)
- [中文完整报告](docs/reports/final/ModStruct_report_zh.md)
- [会议论文提纲](docs/reports/final/ModStruct_conference_paper_outline.md)

## 分析路线

| 阶段 | 脚本编号 | 目标 |
|---|---:|---|
| 输入与坐标审计 | 02–04 | 下载参考资产、核验版本与坐标映射、审计输入完整性 |
| HEK293T 主分析 | 05–08b | 构建分析集、估计关联、预测消融与结构诊断 |
| HeLa 外部复核 | 09–10 | 独立建集、结构计算、效应复核与冻结模型迁移 |
| 发布完整性 | 11 | 生成或验证总发布 manifest |
| Phase 1 | 12–17 | Meta 分析、条件依赖、结构熵、异构体与窗口敏感性、非线性模型 |
| Phase 2 | 18–22 | 公共数据审计、体内/体外配对、跨技术复核与准入门槛 |
| Phase 3 | 23–27 | m5C、m7G、Nm 等多修饰来源筛选、文件审计与分层关联 |
| Phase 4 | 28–31 | 扰动来源审计、配对方向性分析、证据等级与模型闸门 |

每个可执行阶段均由 `config/<stage>_*.yaml` 固定输入、阈值和输出路径。阶段编号是项目的稳定接口，具体规则见 [文件命名规范](docs/FILE_NAMING_CONVENTION.md)。

## 仓库结构

```text
RNAmodStruct/
├── config/                     # 分阶段冻结的 YAML 配置
├── data/
│   ├── raw_processed/          # 不可变下载数据（Git 忽略，manifest 管理）
│   ├── reference/              # 基因组、转录组和注释参考（Git 忽略）
│   ├── interim/                # 可再生中间数据（Git 忽略）
│   └── final/                  # 分析就绪数据（Git 忽略）
├── docs/
│   ├── reports/final/          # 正式报告、答辩和论文提纲
│   ├── reports/progress/       # 阶段性报告
│   └── research_notes/         # 研究方案与路线记录
├── metadata/
│   ├── audits/                 # 输入、参考和坐标审计
│   ├── decisions/              # 冻结决策与预注册边界
│   ├── manifests/              # 下载、来源和发布清单
│   ├── provenance/             # 字段语义、环境和来源验证
│   └── software/GLORI-tools/   # 固定版本的上游工具子模块
├── results/
│   ├── figures/                # 发布图及 QA 记录
│   ├── reports/                # 分阶段技术报告
│   ├── status/                 # 机器可读的运行状态
│   ├── tables/                 # 汇总统计与模型评估表
│   ├── models/                 # 可再生模型产物（Git 忽略）
│   └── logs/                   # 运行日志（Git 忽略）
├── src/                        # 分析脚本与共享模块
├── tests/                      # 阶段测试与完整性测试
├── environment.yml            # Conda 环境定义
└── README.md
```

## 环境安装

推荐使用 Conda/Mamba 从冻结环境创建：

```bash
conda env create -f environment.yml
conda activate rnamodstruct
```

环境基于 Python 3.13.15，主要依赖包括 NumPy、pandas、SciPy、scikit-learn、statsmodels、PyYAML、Matplotlib、Biopython、PyArrow 和 ViennaRNA 2.7.2。完整环境审计见 [`metadata/provenance/environment.txt`](metadata/provenance/environment.txt)。

本项目不需要 GPU；完整流程建议至少准备约 16 GB 内存。

本机当前可直接使用：

```powershell
E:\ancd\envs\my_pytorch\python.exe -m pytest tests -q
```

## 数据准备

大型原始数据、参考文件和可再生产物不会提交到 Git。恢复数据时以 manifest 为准：

- [`metadata/manifests/source/source_manifest.csv`](metadata/manifests/source/source_manifest.csv)：主线数据来源、实验条件和本地路径；
- [`metadata/manifests/download/download_checksums.csv`](metadata/manifests/download/download_checksums.csv)：原始下载文件校验值；
- [`metadata/manifests/source/02_reference_manifest.csv`](metadata/manifests/source/02_reference_manifest.csv)：参考基因组、转录组和注释版本；
- [`metadata/manifests/download/phase3_download_manifest.csv`](metadata/manifests/download/phase3_download_manifest.csv)：多修饰扩展数据下载记录。

参考资产可由以下脚本准备：

```bash
python src/02_download_reference_assets.py
```

HeLa icSHAPE 的 `GSE145805_RAW.tar` 由分析脚本直接读取内部成员 `GSM4333258_HeLa.out.txt.gz`，无需手动解压。

## 运行与验证

主线阶段应按编号顺序执行。每个脚本默认读取对应 YAML 配置，并将产物写入 `data/`、`metadata/` 或 `results/` 的指定分类目录。

最小验证流程：

```bash
python -m pytest tests -q
python src/11_verify_release.py --verify
```

需要重新冻结总发布清单时：

```bash
python src/11_verify_release.py --write-manifest
python src/11_verify_release.py --verify
```

总清单位于 [`metadata/manifests/release/release_manifest.csv`](metadata/manifests/release/release_manifest.csv)。Phase 1–4 另有独立清单，便于验证各阶段冻结范围而不混淆不同证据层级。

## 关键分析约定

- 响应变量：覆盖度加权的 `combined_ratio = (Acov1 + Acov2) / (AGcov1 + AGcov2)`；
- 主结构指标：修饰位点两侧 −10..−1 与 +1..+10 的有效 reactivity 均值 `R_flank10`，每侧有效覆盖率至少 70%；
- 数据划分：以 `analysis_gene` 与相同 201-nt 序列构成的连通分量为分组单位，冻结为 80/20；
- 主预测比较：HEK293T 使用 ΔMAE = MAE(M2) − MAE(M4)，HeLa 迁移使用 MAE(M1) − MAE(M3)；
- Ridge 超参数网格必须包围最优值，否则阶段直接失败，不冻结边界最优模型；
- 观察性关联、预测性能和扰动方向性均不得自动解释为因果机制。

冻结偏差与模型决策记录：

- [`metadata/decisions/08_prediction_modeling_decision_log.md`](metadata/decisions/08_prediction_modeling_decision_log.md)
- [`metadata/decisions/09_hela_replication_decision_log.md`](metadata/decisions/09_hela_replication_decision_log.md)
- [`metadata/decisions/phase4_preregistration.md`](metadata/decisions/phase4_preregistration.md)

## 图表预览

| HEK293T 关联 | 预测消融 | HeLa 复核与迁移 |
|---|---|---|
| ![HEK293T association](results/figures/07_hek293t_association_overview.png) | ![Prediction ablation](results/figures/08_hek293t_prediction_ablation.png) | ![HeLa replication](results/figures/10_hela_replication_transfer_figure5.png) |

矢量版和投稿版文件位于 [`results/figures/`](results/figures/)。

## 证据边界

本仓库当前支持的是小效应关联、外部复核、跨技术敏感性以及扰动后的总体方向性描述。它不支持以下表述：

- RNA 修饰变化已经被证明会导致特定位点结构变化；
- 预测模型中的特征重要性等同于生物学机制；
- 未达到显著性的结果证明真实效应为零；
- 跨细胞系、跨技术或跨研究的结果可以忽略测量尺度直接合并。

未来若获得同一样本的位点级修饰定量、救援/失活化合物对照及独立验证数据，应建立新的确认性阶段，而不是回写本轮探索性结果。
