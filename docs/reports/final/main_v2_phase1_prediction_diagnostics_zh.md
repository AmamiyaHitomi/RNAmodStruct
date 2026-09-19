# ModStruct main_v2：投稿前修订阶段 1 冻结预测诊断

日期：2026-09-18  
分析角色：使用已经保存的逐位点预测进行事后诊断；未恢复、重训或重新选择原 M0–M4 模型。HEK 内部测试已被复用；新 HeLa 排除集合不是新的独立验证。  
输入与代码：[阶段 1 脚本](../../../src/presubmission_stage1_diagnostics.py)、[来源与结果哈希](../../../results/presubmission_stage1/provenance.json)、[阶段 0 位点冻结](../../../metadata/manifests/presubmission_phase0/freeze_manifest.json)。

## 1. 冻结结果复核

已保存的 HEK293T 测试预测覆盖 825 个 `site_id`，与固定 test assignment 完全一致；已保存的 HeLa 预测覆盖 24,960 个 `site_id`，与阶段 0 冻结清单及 `joint_group` 完全一致。重算 M0–M4 的 MAE、RMSE、Spearman、R²，并与既有 HEK 测试、HeLa 全体及原 strict 表逐项比较：**60/60 项通过 `1e-8` 容差，最大绝对差约 `7.73e-10`**。个别 M0 差异来自已保存逐位点预测的小数位精度；没有用重新拟合结果替代冻结结果。逐项差值见[复核表](../../../results/presubmission_stage1/published_metric_reproduction.csv)。

## 2. 开发集常数基线与主要误差

只从 3,271 个 HEK293T development 位点的 `combined_ratio` 计算中位数，固定为 **0.29487179**，再用于 HEK 复用 test 和全部 HeLa 集合。该基线标为 `C_dev_median`。下表 MAE 均为修饰比例单位，数值越低越好；完整 RMSE、Spearman、R²、带符号误差、预测分位数和裁剪率见[诊断指标表](../../../results/presubmission_stage1/diagnostic_metrics.csv)。

| 人群（位点数） | C_dev_median | M0 | M1 | M2 | M3 | M4 |
|---|---:|---:|---:|---:|---:|---:|
| HEK 复用测试（825） | 0.211304 | 0.208203 | 0.205464 | 0.211930 | 0.205490 | 0.211878 |
| HeLa 全体合格（24,960） | 0.223185 | 0.205155 | 0.216084 | 0.218914 | 0.216725 | 0.219470 |
| HeLa 排除 development 重叠（22,256） | 0.223860 | 0.205557 | 0.218589 | 0.222095 | 0.219215 | 0.222645 |
| HeLa 排除可审计 HEK 主分析重叠（21,351） | 0.224648 | 0.205697 | 0.218725 | 0.222153 | 0.219354 | 0.222706 |

按 `ΔMAE = MAE(基线) − MAE(增加特征后模型)`，HeLa 全体中 M0→M1 为 **−1.0929 个百分点**，M1→M3 为 **−0.0641 个百分点**。在新增 21,351 位点集合中，M1→M3 为 **−0.0629 个百分点**。这些是既有冻结模型的描述性点估计；最后一个集合与比较属于事后敏感性分析，不由本阶段新增置信区间。

HeLa 全体中，开发集常数基线 MAE 为 0.223185，高于 M0 的 0.205155。因此 M0 的较好迁移表现不是仅因与未经训练的常数预测比较过弱。与此同时，M1 在 HeLa 比 M0 差，后续结构增量的阴性结果必须继续限定在具体冻结模型及特征设置内，不能归因于结构本身无价值。

## 3. 系统误差、裁剪与校准

本报告定义带符号误差为**预测值减实测值**。HeLa 全体 M0 的平均误差为 `−0.02824`，M1 为 `+0.02741`，M3 为 `+0.03073`；即 M0 平均偏低、M1/M3 平均偏高。M0 平均预测为 0.3909，M1 为 0.4465，M3 为 0.4498。图中的 M0 预测分布也比含序列模型更集中。这些现象提示冻结迁移的预测位置和分布改变，尚不能单凭诊断确定改变的原因。

HeLa 全体中，M0 无边界裁剪，M1 和 M3 各有 4/24,960 个预测恰在 0 或 1，M2 有 12 个、M4 有 11 个；裁剪比例均低于 0.05%。裁剪并非这里大多数位点的误差来源。HEK 复用测试的各模型裁剪数量与其比例亦列于指标表。

[校准表](../../../results/presubmission_stage1/calibration.csv)和[校准图](../../../results/presubmission_stage1/calibration.png)把**预测值**按 HEK development 标签的十分位边界分箱，再报告每箱平均预测与平均实测。分箱边界只由 HEK development 标签确定，未用 HeLa 标签重新分箱。HeLa 全体按该固定分箱计算的描述性加权绝对校准差为 M0 `0.03015`、M1 `0.03546`、M3 `0.03801`；该汇总不是置信区间，也不用于 HeLa 模型选择。[预测分布图](../../../results/presubmission_stage1/prediction_distributions.png)显示四组人群的冻结模型分布和开发集中位数常数位置。

常数预测在全部位点相同，Spearman 无定义；指标表将其留空而不是记为 0。其 R² 可为负值，这是相对于各评估人群均值的正常定义结果。

## 4. HeLa oracle 仅作诊断

[单独 oracle 表](../../../results/presubmission_stage1/oracle_diagnostics.csv)使用各 HeLa **评估集合自身的标签**计算中位数：全体为 0.332494，原 strict 为 0.337254，新增更宽排除集合为 0.340708。相应 MAE 为 0.221350、0.221526、0.221978。这些数字泄露了评估标签，只用于描述标签分布，**不是可部署基线，也不参与冻结模型选择**。

## 5. 输出与阶段验收

- [逐位点预测表](../../../results/presubmission_stage1/site_predictions_with_constant.csv)保留 `site_id`、gene、joint group、实测值、M0–M4 冻结预测及开发集常数预测；HeLa 子集行按所属人群重复列出，不当作独立样本相加。
- [诊断指标表](../../../results/presubmission_stage1/diagnostic_metrics.csv)、[校准表](../../../results/presubmission_stage1/calibration.csv)、两张图和[来源与结果哈希](../../../results/presubmission_stage1/provenance.json)均已生成。
- 原 HEK test、HeLa 全体与原 strict 的 M0–M4 结果按预设容差全部重现；阶段 0 三组 HeLa `site_id` 哈希也再次通过核对。

**阶段 1 状态：完成。** 下一阶段进入预先固定的 M0＋实测结构、去 abundance 重训和相同 3-mer／训练范围下的 Ridge 与 HGB 匹配比较。新增结果继续按事后敏感性分析报告。
