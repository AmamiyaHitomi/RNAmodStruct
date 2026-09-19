# ModStruct main_v2：投稿前修订阶段 0 证据与分析冻结记录

记录日期：2026-09-18  
性质：在既有 HEK293T 测试和 HeLa 结果已被查看之后作出的事后分析规划；不构成预注册或新外部验证。  
依据：[投稿前审阅](main_v2_presubmission_review_zh.md)、[修订路线图](main_v2_revision_roadmap_zh.md)、项目决策日志、配置、运行状态、论文源码和数据表。

## 1. 已核实的分析时间线

下表的时间使用记录文件原有时区；带 `Z` 的时间为 UTC。当前 Git 历史仅有少量整批整理提交，关键决策文件在 2026-09-13 的整理提交中同时进入版本库，不能用该提交时间冒充各决策的首次发生时间。

| 事件 | 证据与可证实范围 | 精确发生时间 |
|---|---|---|
| 首次查看 HEK293T 内部测试结果 | [阶段 08 决策日志](../../../metadata/decisions/08_prediction_modeling_decision_log.md) 明确记录：有缺陷的缩放版本结果已被查看，随后才做修正；原始 M1/M3 MAE 分别为 0.2068117/0.2068004。 | `unknown`；日志仅写 2026-09-10 |
| 训练折内缩放修正 | 同一日志记录 3-mer 频率与 `R_flank10` 原先未缩放，修正后 M1/M3 MAE 为 0.2053234/0.2052843；内部测试集因此被复用。 | `unknown`；日志仅写 2026-09-10 |
| 早期可用模型与主比较 | 在预测结构缺席时，M0、M1、M3 可用，日志明确称可用主比较为 M1 − M3；M2/M4 当时未完成。 | 顺序明确；首次冻结时点 `unknown` |
| 预测结构加入 | [08A 状态](../../../results/status/08a_hek293t_predicted_structure_run_status.csv) 记录 2026-09-10T14:14:22Z–14:14:25Z 的一次成功生成；阶段 08 日志随后记录 M2/M4 完成。 | 该次运行时间已知；首次提出加入的时间 `unknown` |
| 完整 M0–M4 矩阵的内部主比较 | [当前 08 配置](../../../config/08_prediction_modeling.yaml) 指定 M2 − M4 为主比较，M1 − M3 为次比较；[补充材料](../../../docs/manuscript/supplement.tex)与[作者说明](../../../docs/manuscript/AUTHOR_NOTES.md)沿用该完整矩阵角色。它与较早的“可用主比较”属于不同阶段，不能追溯合并。 | 首次改为完整矩阵主比较的时间 `unknown` |
| alpha 网格扩展 | 阶段 08 日志记载 2026-09-11 基于 development CV 将上界从 100 扩至 3000；当前 M1–M4 选择 300，且为网格内部点。[最终 08 状态](../../../results/status/08_hek293t_prediction_modeling_run_status.csv)为 2026-09-11T04:21:36Z–04:23:53Z。日志声称未查阅 HEK 测试或 HeLa 标签决定扩展；现有提交历史不能独立证明该负面事实。 | 网格决策精确时间 `unknown`；最终重跑时间已知 |
| 首次 HeLa 关联与冻结预测评估 | [09 决策日志](../../../metadata/decisions/09_hela_replication_decision_log.md)称先完成 09B/09C，再记录 G5，且最终 09C 因 alpha 审计在 2026-09-11 重跑。[当前 09C 状态](../../../results/status/09c_hela_model_transfer_run_status.csv)为 2026-09-11T04:27:31Z–04:28:17Z。09B 状态载有无时区的 `2026-09-10T22:55:29`，不据此推断其与 UTC 记录的精确先后。 | 首次 09C 查看时点 `unknown`；最终重跑时间已知 |
| 401nt 后续分析 | [阶段 16 状态](../../../results/status/16_window_401_sensitivity_run_status.csv)为 2026-09-11T13:54:41Z–13:56:19Z，晚于当前 09C 重跑。 | 有记录的该次运行时间已知；首次探索时间 `unknown` |
| HGB／3-mer 后续分析 | [阶段 17 状态](../../../results/status/17_nonlinear_models_run_status.csv)为 2026-09-11T13:56:21Z–13:58:17Z，晚于阶段 16 和当前 09C 重跑。 | 有记录的该次运行时间已知；首次探索时间 `unknown` |

**时间线结论。** 内部 HEK 测试集在缩放修正前已被查看，此后内部结果只能作描述性证据。原始冻结 HeLa 迁移是项目历史中的主外部检查；阶段 16/17 和本次新增比较都发生在 HeLa 结果已知之后，不能再称为新的“未触碰”验证。上述记录只证明可见文件中的叙述与运行，不证明没有未记录的中间尝试。

## 2. 验证角色冲突及本阶段处理

[阶段 17 配置](../../../config/17_nonlinear_models.yaml)把 HeLa 写为 `untouched_external_evaluation`。这与阶段 17 晚于 09C 的运行记录、[主文](../../../docs/manuscript/main_v2.tex)及补充材料将 HGB 称为后续敏感性分析的表述不一致。阶段 17 还改变了序列表示、使用三折 CV，并在完整 4,096 个 HEK 模型位点上重拟合，不能将其与原阶段 08 的 3,271 位点 development 训练直接解释为仅更换学习器。

本阶段在[固定分析规格](../../../metadata/manifests/presubmission_phase0/analysis_spec.json)中正式将其标为 `subsequent_sensitivity_not_untouched_external_evaluation`。现有阶段 17 配置保留为历史输入；在文稿与配置统一修订时，应同步更正该历史标签，不能通过事后重命名改变真实时间线。主文 Fig. 4 当前仍使用 `reject an experimental-structure gain`，其措辞修订安排在最终文稿阶段。

## 3. HeLa 位点集合及重叠账目

从 09C 冻结逐位点预测表取 24,960 个模型合格位点，以其中的 `site_id`、`joint_group`、`analysis_gene` 和 201nt 序列哈希为准。和 HEK gene 或**完全相同的 201nt 序列**相交即排除；不把此规则称为广义同源排除。构建脚本先核对原始 HeLa 模型资格、序列哈希、HEK assignment 与原 `overlaps_development` 标志，再写出 ID 清单。

| 固定集合 | 位点 | joint groups | 位点列表 SHA-256* | 验证角色 |
|---|---:|---:|---|---|
| 全体合格 | 24,960 | 4,685 | `41ed5d879b88ca206428ed7f9f83afac9749eb77e716023def6acd2d8e6def5f` | 原始冻结 HeLa 主外部评估的人群；新增对比为事后分析 |
| 排除 HEK development 重叠 | 22,256 | 4,046 | `93cd7a430eeb3b385009e6ecf5baba226df1e2d31db8ee06634f79230168b018` | 现有 strict 子集的准确名称；新增对比为事后分析 |
| 排除可审计 HEK 主分析 4,409 位点重叠 | 21,351 | 3,816 | `dce8e7d3a29cbfb705fb8a6821e799ee7f46fcf2054707e8973e3bbc7bf909b3` | 新建的更宽排除敏感性集合，仍非未触碰验证 |

\*哈希对象为 `site_id` 字典序排列、每个 UTF-8 ID 后加 LF（末行也加 LF）的文本，不是 CSV 文件字节。三个 CSV 各自的文件 SHA-256 和输入文件 SHA-256 见[冻结 manifest](../../../metadata/manifests/presubmission_phase0/freeze_manifest.json)。[逐位点标志表](../../../metadata/manifests/presubmission_phase0/hela_site_cohort_flags.csv)保留 gene、序列及集合标志；[生成脚本](../../../src/freeze_presubmission_phase0.py)可审计具体计算。

| HEK 对照范围 | gene 重叠 | 完全相同 201nt 重叠 | 两者交集 | 合并排除 |
|---|---:|---:|---:|---:|
| 3,271 个 development 位点 | 2,704 | 1,542 | 1,542 | 2,704 |
| 4,096 个模型位点（含 825 个已查看 test） | 3,496 | 1,951 | 1,951 | 3,496 |
| 4,409 个主关联分析位点（含额外 313 个非模型位点） | 3,609 | 1,951 | 1,951 | 3,609 |

相对现有 strict 集合，4,409 位点排除规则额外移除 905 个 HeLa 位点。本批数据中相同 201nt 序列的重叠均同时属于 gene 重叠；这一计数不说明不存在较弱的序列同源。4,409 位点是当前可审计的 HEK 主分析人群，不能证明已覆盖历史上所有曾被查看的 HEK 候选位点。因此第三组不得简称为“排除所有已查看 HEK”。

## 4. 固定的新增分析规格

[analysis_spec.json](../../../metadata/manifests/presubmission_phase0/analysis_spec.json)在运行新增模型比较之前列明：既有 M0–M4 冻结诊断、只用 HEK development 标签的常数中位数基线、M0 与 M0＋实验结构、去 abundance 后的 M0/M0＋结构/M1/M3 重训，以及相同 3-mer、相同 3,271 位点训练范围和五折分组下的 Ridge/HGB 比较。所有候选模型保存，预处理及调参只在训练折内；HeLa 不参与选择。新比较全部为已知 HeLa 结果后的事后敏感性分析。

主要误差为位点加权 MAE，`ΔMAE = MAE(基线) − MAE(增加特征后模型)`，正值表示误差降低；比例单位乘 100 才是百分点。对同一 HeLa 集合的所有预定比较，使用相同的联合分组配对 bootstrap 抽样索引：joint group 为单位、2,000 次、种子 `20260912`、百分位 95% 区间。不同集合的组数不同，各自有独立的索引空间。报告字段为 `CI_scope=fixed_fitted_models`；区间不包含重抽训练位点、重新调参或改变特征表示的不确定性。HEK 复用测试结果与新增 HeLa 敏感性结果均按其历史角色解释。

## 5. 阶段 0 验收与后续待处理

- [x] 历史主比较的两个阶段、内部测试复用及当前运行记录已核对；无法证实的首次时间标为 `unknown`。
- [x] 当前阶段 17 validation role 的冲突已定位，并在新分析规格中给出正确角色。
- [x] 三组 HeLa `site_id`、gene/相同 201nt 重叠与交集已固定并哈希；现有 24,960/22,256 数目重现。
- [x] 新增比较、训练范围、ΔMAE 方向、配对 bootstrap 和区间范围已在新增模型运行前固定。
- [ ] 历史 `config/17_nonlinear_models.yaml` 文字标签、主文 Fig. 4 等既有文件的修订尚未执行；需在文稿与配置同步修订时处理，且保留这份历史审计记录。

本记录不认定那些只能从决策日志获知的“未查看某标签”叙述已被独立证明，也不认定第三组是完全无历史接触的独立验证集。
