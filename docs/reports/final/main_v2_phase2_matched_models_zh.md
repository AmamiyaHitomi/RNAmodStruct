# ModStruct main_v2：投稿前修订阶段 2 匹配模型比较

日期：2026-09-18  
分析角色：在原 HEK293T 测试与 HeLa 结果已知后，按[阶段 0 固定规格](../../../metadata/manifests/presubmission_phase0/analysis_spec.json)运行的**事后敏感性分析**。不构成预注册或新增未触碰外部验证。  
执行代码：[presubmission_stage2_matched_models.py](../../../src/presubmission_stage2_matched_models.py)；完整输入、代码和输出哈希见[provenance.json](../../../results/presubmission_stage2/provenance.json)。

## Material Passport

- Origin：阶段 0 固定的 C1、C2、C3 比较；输入为项目已处理数据与既有分组。
- Verification Status：ANALYZED。运行成功，并独立核对输出哈希、特征清单、训练位点哈希及一个配对 bootstrap 区间；未声称已由另一套实现完整重跑。
- Population：固定 HEK293T development 3,271 位点训练；HEK 825 位点复用 test 仅作描述；HeLa 三组固定集合分别为 24,960、22,256、21,351 位点。

## 1. 方法和事先固定的边界

三组新比较均沿用 HEK 原 development/test 划分。参数只由 development 的五折 `GroupKFold` 选择；各折的类别发现、缺失填补、缩放及结构特征标准化只在该折训练部分拟合。最终模型只在全部 3,271 个 development 位点拟合一次，随后才对复用 HEK test 和 HeLa 作预测。HeLa 标签不参与模型、参数或特征选择。预测值与原流程一样裁剪到 [0, 1]。

| 组别 | 固定比较 | 特征与训练条件 |
|---|---|---|
| C1 | M0 vs M0＋实测 `reactivity_mean_flank10`（M0R） | 原 201nt 位置 one-hot＋3-mer 编码器；仅 M0R 增加一个实测结构特征 |
| C2 | 去 abundance 后的 M0/M0R 与 M1/M3 | 在训练前从四个模型的特征清单中移除 `log1p_icshape_abundance_rpkm`；所有模型重新训练并以 development CV 重新选参，未在 HeLa 端填零代替删除 |
| C3 | 同一 3-mer 表示下 Ridge vs HGB 的 M0–M4 | 两学习器共用 3,271 位点、相同五折和对应模型的相同特征清单；此表示不含 201nt 逐位 one-hot，预测结构块采用阶段 17 的 8 个标量特征 |

Ridge 搜索固定的九个 alpha（0.01 至 3000）；HGB 使用阶段 17 配置已有的两个候选。所有候选的逐折和均值结果保存在[development_cv_candidates.csv](../../../results/presubmission_stage2/development_cv_candidates.csv)，16 个预定模型的选择、参数和特征数保存在[selected_models.csv](../../../results/presubmission_stage2/selected_models.csv)。11 个 Ridge 选择均不是 alpha 网格端点。HGB 只有两个候选，其选择并不证明更广超参数空间的最优性。

以 `ΔMAE = MAE(基线) − MAE(增加特征后模型)` 报告特征增量，正值为误差降低。外部区间对每组 HeLa 集合按原 `joint_group` 配对抽样 2,000 次，种子 `20260912`；同一集合的所有预定比较使用同一组抽样索引。`CI_scope=fixed_fitted_models`：区间不包括重新抽取训练位点、重新调参或改变特征表示的不确定性。HEK 复用 test 只给描述性点值。所有对比及三组集合均在[paired_contrasts.csv](../../../results/presubmission_stage2/paired_contrasts.csv)。

## 2. HeLa 全体的关键结果

下表 MAE 为修饰比例单位；Δ 与区间为百分点。区间是上述固定已拟合模型的事后敏感性区间，不是完整训练流程的重复研究区间。

| 比较 | 基线 MAE | 增加后 MAE | ΔMAE（百分点） | 95% 配对分组区间 |
|---|---:|---:|---:|---:|
| C1：M0 → M0R | 0.205155 | 0.206711 | −0.1556 | [−0.1884, −0.1231] |
| C2 去 abundance：M0 → M0R | 0.205023 | 0.206471 | −0.1448 | [−0.1738, −0.1160] |
| C2 去 abundance：M1 → M3 | 0.217048 | 0.217469 | −0.0421 | [−0.0470, −0.0374] |
| C3 3-mer Ridge：M1 → M3 | 0.211734 | 0.212630 | −0.0896 | [−0.0994, −0.0794] |
| C3 3-mer HGB：M1 → M3 | 0.188412 | 0.190460 | −0.2048 | [−0.2517, −0.1556] |
| C3 相同 M1 特征：Ridge → HGB | 0.211734 | 0.188412 | +2.3322 | [+2.1509, +2.5162] |

**结果解释。** C1 的实测结构在 HEK development CV 中使 M0 平均 MAE 从 0.203320 降至 0.202840，但在 HeLa 全体使误差增加。去 abundance 后，M1→M3 的 HeLa 增量仍为负，幅度比原冻结 M1→M3 的约 −0.0641 个百分点小。相同 3-mer 表示与训练范围下，HGB 的整体 MAE 比 Ridge 低，但在各自 HGB 基线上增加实测结构仍使 HeLa 误差增加。因此“更换学习器改善总体预测”与“实测结构是否增加预测价值”必须分开表述。

其余预定结果也应同时保留。C3 的 3-mer Ridge M1→M2（加入预测结构）在 HeLa 全体为 **+0.0211** 个百分点，区间 [+0.0028, +0.0411]；原 strict 集合为 +0.0165，区间跨零。相同 HGB 对比在 HeLa 全体为 **−0.0738** 个百分点。C3 HGB M0→M1 在全体为 +0.2434 个百分点，但原 strict 集合为 −0.0539，区间跨零。这些差异说明表示、学习器和排除人群会改变某些特征块的表现；不能只摘取有利子集或把某一结果推广至所有结构特征。所有 M0–M4、Ridge/HGB 及三组集合的完整 MAE、RMSE、Spearman、R² 见[evaluation_metrics.csv](../../../results/presubmission_stage2/evaluation_metrics.csv)。

## 3. 三组 HeLa 集合与内部测试的证据角色

对最直接的 M1→M3 实测结构增量：

| 比较 | HeLa 全体 24,960 | 排除 development 重叠 22,256 | 排除可审计 HEK 主分析重叠 21,351 |
|---|---:|---:|---:|
| C2 去 abundance，Ridge | −0.0421 | −0.0418 | −0.0422 |
| C3 3-mer，Ridge | −0.0896 | −0.0885 | −0.0889 |
| C3 3-mer，HGB | −0.2048 | −0.1749 | −0.1754 |

表中单位为 ΔMAE 百分点，完整配对区间在对比表。第三组是事后构造的更宽排除集合；它只排除当前可审计 HEK 主分析位点的相同 gene 或完全相同 201nt 序列，不能证明排除所有曾被查看的 HEK 位点或广义同源序列。

复用 HEK test 中，C1 M0→M0R 为 +0.0755 个百分点，而 HeLa 全体为 −0.1556；C2 去 abundance 的 M1→M3 在复用 test 约 +0.0006、HeLa 全体为 −0.0421。这种方向差异是迁移诊断，不能将已查看的内部 test 重新当作独立确认。阶段 2 未按该 test 或 HeLa 的表现回头选择模型。

## 4. 运行与核验记录

- 五折 development 搜索保存 **654 条候选记录**（逐折及候选均值）；最终保存 **16 个模型文件**，每个模型都有选择参数、训练位点哈希、特征清单及裁剪规则，见[models 目录](../../../results/presubmission_stage2/models)。
- [HEK 逐位点预测](../../../results/presubmission_stage2/HEK_reused_test_predictions.csv)和[HeLa 逐位点预测](../../../results/presubmission_stage2/HeLa_predictions.csv)保留 `site_id`、gene、joint group。生成 **64 条人群×模型指标**和 **64 条人群×预定对比**。
- C1 重拟合的 M0 与原冻结 M0 的 HEK 测试预测最大绝对差为 `2.87e-8`，HeLa 逐位点预测完全一致；这为新比较的原编码器与训练范围提供了锚点。
- 独立复核确认：所有输出文件 SHA-256 与 provenance 相符；16 个模型的训练位点哈希完全相同；C2 特征名不含 abundance；C3 的 Ridge/HGB 对应模型特征名逐项一致；HeLa 全体 C1 对比的 bootstrap 索引哈希与区间可独立重算。

**阶段 2 状态：完成。** 结果支持将外部阴性结论继续限定为已评估的人群、特征表示、学习器和迁移设置。下一阶段按路线图补齐体内／体外估计量尺度、`NormeRatio` 敏感性和 HeLa 结构覆盖人群的平衡诊断，然后依据全套结果修订主文。
