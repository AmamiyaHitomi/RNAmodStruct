# ModStruct 投稿前修订实施方案

本方案把投稿前审阅转成可执行任务。新增分析均应标记为后续／事后敏感性分析；它们不能被改写成外部预注册结果，也不能恢复已经查看过的内部测试集。

## A. 先完成文稿和记录修订

1. 修改引言和§4.2，承认SMART-m6A包含定量修饰强度预测，把差异改写为实测结构、条件增量、冻结迁移和可审计整合。
2. 在正文引用补充材料已有的GLORI字段公式、M0特征清单、Ridge alpha网格及bootstrap区间边界。
3. 将strict子集命名为“排除development重叠”，说明其不排除全部已查看HEK数据，也不处理广义同源序列。
4. 将Fig.4标题从“reject an experimental-structure gain”改为更窄的“did not improve the evaluated transfer models”。
5. 整理分析时间线：首次测试查看、缩放修正、预测结构加入、alpha扩展、HeLa评估、401nt和HGB后续分析。无法证明的时间点写为`unknown`。
6. 检查`config/17_nonlinear_models.yaml`、正文和补充材料，统一后续HGB分析的validation role，避免使用`untouched_external_evaluation`造成误读。

## B. 先做无需重训练的预测诊断

读取已保存的逐位点预测；如只有汇总表，恢复原冻结模型并用原结果核对，不以重新拟合替代原结果。对HEK293T测试集、HeLa全体和现有strict集合报告：MAE、RMSE、Spearman、R2、误差均值、clip比例、预测分布和按开发集分箱确定的校准曲线。

增加开发集标签中位数常数预测。HeLa中位数只能作为单独的oracle诊断，不能当作可部署基线。

验收条件：重建后的M0–M4结果与现有结果表相符；常数基线只使用HEK293T development数据；所有预测保留site_id。

## C. 补充最关键的匹配模型比较

固定原开发／测试划分和五折分组，新增：

| 比较 | 目的 |
|---|---|
| M0 vs M0＋实验结构 | 检查结构能否改善简单背景模型 |
| M0、M0＋结构、M1、M3去除abundance后重训 | 检查HeLa缺失abundance的影响 |
| 同一3-mer表示、同一训练范围的Ridge与HGB | 检查阴性结果是否依赖学习器 |

去除abundance必须在训练前删除特征并重新调参；不能只在HeLa端填零。所有预处理、类别发现、缩放和alpha选择均在训练折内完成。候选模型全部保存，不根据HeLa结果挑选。

输出：`baseline_ablation_metrics.csv`、匹配学习器CV表、外部评估表和逐位点预测。

## D. 固定外部评估集合和bootstrap规则

预先保存并hash三组HeLa site ID：全体合格集合、现有development-overlap-excluded集合、可审计范围内的all-reviewed-HEK overlap-excluded集合。记录基因排除、相同201nt序列排除及交集。

每个预定对比使用相同的联合分组配对bootstrap索引，采样单位为joint group，重复2,000次。区间字段必须写明`CI_scope=fixed_fitted_models`。不把新增排除集合的结果称作未经查看的验证。

## E. 补强关联分析的解释

1. 从stage 19现有结果提取原始斜率、各条件SD、标准化斜率和区间；定义并单独报告`Δ reactivity`。
2. 保留count-pooled Ratio为主标签；增加`NormeRatio`敏感性，不能用两重复简单均值代替。
3. 为HeLa补充结构覆盖合格／不合格人群的协变量平衡表。
4. 用已有样条诊断绘制调整曲线及反应性分布，不依据曲线反复寻找有利分组。
5. 基因内／基因间分解作为可选增强，不将基因聚类稳健SE误写为对基因层面混杂的完全控制。

## F. 可复现提交包

1. 让构建入口明确指定`main_v2`或新的修订稿；当前构建脚本默认编译`main`，不能直接视为main_v2构建。
2. 扩展release manifest，覆盖论文源码、图表、配置、环境、结果表和必要子模块；保留旧版本文件。
3. 分别验证A级复现（从冻结表重建图表和PDF）与B级复现（从processed输入重建划分、模型和指标）。
4. 检查`validate_evidence.py`是否识别main_v2和新增结果；旧检查通过不代表新增结果已被验证。
5. 记录真实Python环境、依赖、CPU、内存和运行时间；不要因补充材料中的版本说明而假定当前环境相同。

## 推荐执行顺序

1. 完成文献定位、分析时间线和措辞修订。
2. 完成预测诊断、常数基线和三个外部集合。
3. 完成M0＋结构、去abundance和匹配Ridge/HGB比较。
4. 完成体内／体外估计量、NormeRatio和HeLa筛选平衡的最低补强。
5. 更新正文、补充材料、图表、README和发布清单，再编译和验证。

训练稳定性重采样、新细胞系或更复杂深度学习模型属于增强项，应在最低修订包完成后再决定。完成标准不是得到正的ΔMAE，而是让每个结论都对应清楚的人群、模型、验证角色和不确定性范围。
