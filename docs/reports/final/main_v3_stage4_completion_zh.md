# main_v3 投稿前修订与复现验收

状态：阶段 4 已完成通用稿的文稿、补充材料、图表、发布清单及 A/B 两级技术验收。目标期刊或会议尚未确定，因而没有声称符合特定投稿模板。作者、单位和最终投稿选择仍由作者填写。

## 修订内容

- 新版位于 `docs/manuscript_v3/`；旧 `docs/manuscript/` 与 `main_v2` 原件保留。主文 10 页，补充材料 5 页。
- 引言与讨论明确承认 [SMART-m6A 原文](https://doi.org/10.1371/journal.pcbi.1014649)既做修饰位点分类，也以 GLORI 连续修饰强度为标签做回归。本文的贡献限定为实测 icSHAPE 的条件增量、明确的背景／序列对照、HEK293T 到 HeLa 的冻结迁移及可审计整合，不对两方法做直接性能排名。
- Fig. 4 标题收窄为实验结构未改善**所评估的冻结迁移模型**。正文补入三个事后匹配比较：M0+实测结构、训练前去 abundance、相同 3-mer 表示和开发集范围下的 Ridge/HGB。完整候选结果仍见 `results/presubmission_stage2/`，未按 HeLa 结果筛选。
- 现有 22,256 位点子集明确称为“排除 HEK293T development 重叠”。新增排除可审计的 4,409 个 HEK 关联位点之基因／相同 201nt 序列重叠后，HeLa 剩 21,351 位点，仅作事后敏感性；两者均不声称完全清除历史查看或广义序列同源。
- 保留 count-pooled Ratio 主标签；补入两个重复分别使用 `NormeRatio` 的敏感性结果、65,050 个双重复检出 HeLa 位点的结构覆盖平衡，以及 3 自由度样条诊断。Stage 19 的 `Δ reactivity` 定义为同一位点体内减体外 flank10 均值，区分该变量的回归斜率与各条件独立标准化斜率之差。
- 通过主文与补充材料指向 GLORI 字段公式、M0 特征清单、Ridge alpha 网格、冻结模型和 `CI_scope=fixed_fitted_models` 的区间边界。基因聚类稳健标准误没有被表述为完全控制基因层面的混杂。

## A 级：冻结结果表到图表与 PDF

在独立 `docs/manuscript_v3/` 目录运行其 `build_assets.py`，从原 `results/tables/` 重绘四幅主图、五张 LaTeX 表并复制 22 个源 CSV；用 `build_v3.ps1` 明确编译 `main_v3.tex` 与 `supplement_v3.tex`。`src/validate_presubmission_v3.py` 验证复制源表哈希、图表存在、引文键、阶段 1 的 60 项原指标复算、阶段 3 HeLa 原模型系数、三组 HeLa 人群及 PDF 的页数、文本边界和 TeX 日志。**44 项检查通过**，PDF 未见未定义引用、overfull box 或页面裁切；主文和补充材料页面均已做联系图检查，补充材料末页在计数修正后重新渲染核对。

`docs/manuscript/validate_evidence.py` 只验证旧版主文及随包数据，不作为 main_v3 的验收。目标投稿模板确定后还需按该模板逐页检查字号、页数、图表与参考文献格式。

## B 级：processed 输入到划分、模型与指标

1. `src/verify_presubmission_stage4_split.py` 从 HEK293T processed 主分析位点重建基因及相同 201nt 序列连通组。4,096 位点的 3,271/825 开发测试划分和 1,205 个分组与冻结表完全一致。
2. `src/replay_presubmission_stage4_primary.py` 在新的 `results/presubmission_stage4/primary_replay/` 目录重建五个历史 Ridge 主模型，直接运行冻结 HeLa 迁移。**12 张输出表逐值一致，最大数值差 0**；运行约 183 秒。
3. `src/replay_presubmission_stage4_models.py` 在独立的 `results/presubmission_stage4/replay/` 目录重训阶段 2 的 **16 个**匹配候选模型。CV、选择、指标、配对对比及逐位点预测共 **六张表逐值一致，最大数值差 0**；运行约 65 秒。

B 级从 processed 输入和已保存的预测结构特征表出发，重建组划分与模型；未重新下载 GEO 原始文件或重新计算 ViennaRNA 结构特征。B 级精确重放不改变原 HEK293T 测试集曾被查看的事实，也不把事后比较变成预注册验证。

## 发布包与环境

权威完整清单为 `metadata/manifests/release/presubmission_v3_release_manifest_full.csv`，附带 `presubmission_v3_release_metadata_full.json`。共 **212** 个带 SHA-256 的工件，覆盖全部 `src/*.py`、`config/*.yaml`、论文源文件、图表、PDF、环境、processed 输入与来源清单、样本 ID、划分与特征字典、历史五个冻结模型、阶段 2 的 16 个模型、结果表和 A/B 验收记录。初版 139 项清单作为生成记录保留；完整清单补入所有分析脚本与配置。Git HEAD 与 GLORI-tools 子模块提交记录在 metadata JSON 中；工作区含未提交修订，不能仅以 HEAD 表示最终文件状态，应使用逐文件哈希。

本次实际环境记录于 `results/presubmission_stage4/runtime_environment.json`：`E:\ancd\envs\my_pytorch\python.exe`，Windows 11，20 个逻辑 CPU、15.63 GiB 内存；NumPy 2.5.2、pandas 3.0.5、SciPy 1.18.1、scikit-learn 1.9.0、statsmodels 0.15.0、patsy 1.0.3、matplotlib 3.11.1、joblib 1.5.3；PDF 由 TeX Live 2026 编译。

## 解释边界

外部配对 bootstrap 的 95% 区间仅描述当前已拟合模型在给定 HeLa 人群上的抽样不确定性；它不覆盖重新抽取训练集、调参、特征表示或更换实验条件。结构覆盖不足与合格 HeLa 位点在测序深度、GC 和转录本位置上并不平衡，关联估计不能外推到全部 GLORI 检出位点。Stage 19 为观察性配对分析，不证明结构对修饰的因果作用。投稿前仍需作者补齐身份信息，确定目标模板并按模板进行最终版式审核。
