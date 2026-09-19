# ModStruct main_v2 阶段 0 完成记录

日期：2026-09-18

[阶段 0 证据审计](main_v2_phase0_evidence_audit_zh.md)中的时间线、三组 HeLa 位点 ID 与哈希、重叠账目、[新增分析规格](../../../metadata/manifests/presubmission_phase0/analysis_spec.json)和 bootstrap 规则已固定。审计中无法核实的首次事件继续标为 `unknown`。新增比较均为已知 HeLa 结果之后的事后敏感性分析。

经用户确认，`config/17_nonlinear_models.yaml` 中 HeLa 的 `role` 已从 `untouched_external_evaluation` 改为 `subsequent_sensitivity_evaluation`。YAML 解析通过；Git diff 仅有这一行，模型参数、输入及现有结果均未修改。此处只更正验证角色的文字标签，不改变阶段 17 的历史时间线。

[冻结 manifest](../../../metadata/manifests/presubmission_phase0/freeze_manifest.json)保存的是更正前输入快照，因此其中 `nonlinear_config` 的 SHA-256 为 `2471c6626f3a082e475c7ec4e8a88653fa6cb856f1dbf1e89e2525a83997f535`；更正后该配置文件的 SHA-256 为 `89cc8dbfabbd959c826bb9a5df892d1d1ef3fcaba1ef4b0b60f5b30abe93710b`。位点清单及其哈希不受标签更正影响。

**阶段 0 状态：完成。** 主文 Fig. 4 标题等论文措辞将在最终文稿阶段依据匹配比较结果修订；这不改变阶段 0 的分析冻结。下一阶段是从保存的冻结逐位点预测进行诊断和常数基线核对。