# ModStruct 重启范围清单

> 本文列出项目从原始构思（chat idea）到执行过程中**所有被主动放弃、推迟或降级**的部分，供重启决策使用。
> 每条给出：出处原文引用 → 原放弃理由 → 重启所需的依赖与数据缺口。
> 依据文档：同目录《ModStruct_本科生研究方案.md》（正式方案优先）、《chat_idea_整理.md》（原始构思）、《ModStruct_完整执行路线图.md》（执行层）、《ModStruct_通俗解读.md》（口径表述）。

## 文档信息（Material Passport）

- Origin Skill: academic-research-suite
- Origin Mode: plan
- Origin Date: 2026-09-11
- Verification Status: UNVERIFIED（本文为范围清单，非研究结论；所有"数据缺口"均需在重启前重新审计确认）
- Version Label: restart_scope_v1
- 说明：本文不改变任何已冻结的配置或已产出的结果，仅作为重启被放弃项的索引。凡提到"未下载/未实现"的资产，以实际目录审计为准。

---

## 使用说明

- 分类共 5 类、24 项。其中第 24 项（HeLa Q3）不是"放弃"而是"尚未完成"，单独列出以免遗漏。
- 每条"重启依赖与数据缺口"只列**重启时需要补的东西**，不重复原方案已具备的能力（如坐标映射脚本、联合表框架可复用）。
- 出处的行号以本文编制时各文件的实际行号为准，引用为逐字摘录。

---

## 第一类：chat idea 里的"招牌愿景"被砍（最核心的放弃）

### 1. 多修饰比较 / structural niches（structural grammar）

- **出处原文**：
  - 《chat_idea_整理.md》§六：*"不要一开始就分析 70 多种修饰。第一步直接做：m6A → RNA structure。然后如果结果有意思，再扩展到：m6A / m1A / m5C / m7G / Ψ / ac4C / A-to-I / 2'-O-Me ... 最终形成一个非常漂亮的 matrix。这实际上就已经开始产生一个非常有意思的 RNA modification structural landscape。"*
  - 《chat_idea_整理.md》§七：*"Different RNA modifications may occupy distinct structural niches."*、*"这就从一个比较普通的 association analysis，变成了一个 RNA modification structural grammar 的问题。"*
  - 《ModStruct_本科生研究方案.md》§一：*"lncRNA、多修饰比较、因果推断、基础模型训练和数据库网站均不作为结题条件。不同修饰的检测技术和作用对象差异较大，初期合并会放大工作量和解释困难。"*
- **原放弃理由**：不同修饰检测技术差异大、作用对象不同，初期合并会放大工作量和解释困难。
- **重启依赖与数据缺口**：
  - 每种修饰需要**单碱基/高分辨率的定量位点数据**，不能退回 RMBase 的二元注释。目前只有 m6A 有 GLORI 级定量；m1A、m5C、m7G、Ψ、ac4C、2'-O-Me、A-to-I 的可用数据质量需逐一侦察（多数是抗体富集 MeRIP，分辨率低，或仅二元位点）。
  - 每种修饰都要重走 01–06 阶段的字段审计 + 坐标对齐（参考基因组版本、链向、坐标起算各自不同）。
  - 需要与 SMART-m6A / StructRMDB 明确差异：它们都没做"多修饰 structural niches"，此处是真正空白。

### 2. 高维条件互信息 I(M;S | sequence)

- **出处原文**：
  - 《chat_idea_整理.md》§11.6：*"不要只计算 correlation，而计算 I(M;S)。然后进一步 I(M;S∣Sequence)……真正想知道的是：RNA modification 与 structure 的关系，是否超越了 sequence 本身所能解释的程度？"*
  - 《ModStruct_本科生研究方案.md》§十一：*"不建议本科生第一版直接估计高维 I(Y;R|X)：高维条件独立检验需要明确假设，估计结果容易受模型误差及稀疏性影响。"*
  - 《ModStruct_完整执行路线图.md》§8：*"高维条件互信息估计不进入本科主线。"*
- **原放弃理由**：高维条件独立检验方法上困难（Shah & Peters 的 hardness 结论），估计不稳。
- **重启依赖与数据缺口**：
  - 方法依赖：GCM（Generalized Covariance Measure，方案参考文献 [^20] Shah & Peters）或模型-based 条件独立检验；不推荐真的去估高维互信息数值。
  - **建议降级路线**：做低维条件依赖检验——控制序列 k-mer + 区域 + 覆盖后，检验"实验结构 R 是否仍与 m6A 比例条件相关"。可直接复用现有联合表（data/final/06_*），无数据缺口。
  - 与现有框架的衔接：M3 vs M1 的 ΔMAE 已是"预测意义上的信息增益"近似，条件互信息把它变成正式统计量。

### 3. Structural association network（二部网络）

- **出处原文**：《chat_idea_整理.md》§11.4：*"建立一个二部网络：RNA modifications → structural features / structural motifs / sequence motifs / transcript regions……最终你可能得到：RNA modifications occupy distinct structural niches。"*
- **原放弃理由**：仅作为"酷炫升级"候选，未正式立项。
- **重启依赖与数据缺口**：
  - 需要结构 motif 注释、序列 motif、转录本区域等节点定义。
  - 方法依赖：二部图构建 + 模块检测（community detection）。
  - 无外部数据缺口，可复用现有特征表；但"多修饰"网络需等第 1 项的数据落地。

### 4. 多数据集合并 + RMBase v3.0

- **出处原文**：《chat_idea_整理.md》§三：核心流程图中结构数据为 *"Human RNA MaP / GSE145805 / GSE50676"*，修饰来源为 *"RMBase v3.0 → Modification sites"*。
- **原放弃理由**：实际执行中被替换为 GLORI（定量）+ icSHAPE（两套）的窄组合。
- **重启依赖与数据缺口**：
  - Human RNA MaP（U2OS）、GSE145805（HeLa icSHAPE，已下载 tar）、GSE50676（PARS）三个结构数据源需重新纳入。
  - RMBase v3.0 是二元位点注释，无法支撑"定量修饰比例"主线，仅可用于多修饰的"位点存在/不存在"层。
  - 若目标是 chat idea 的"atlas"级，数据源需重新扩充并逐一审计。

### 5. 四层因果框架的第 2–4 层（Directionality → Causal → Mechanism）

- **出处原文**：《chat_idea_整理.md》§九：*"第一层：Association……第二层：Directionality……第三层：Causal evidence……第四层才是：Mechanism。"* 项目实际只做了第一层（Association）。
- **原放弃理由**：第一阶段不需要证明 causal，且纯计算数据无法完成第 3–4 层。
- **重启依赖与数据缺口**：
  - 第 3 层（Causal evidence）需"相同序列的修饰 vs 未修饰 RNA 结构对照"或"writer 敲除 → 修饰↓ → 结构变化"的扰动数据。
  - 数据缺口：需实验协作或借用已发表的扰动数据集（如 Mettl3 KO 的 icSHAPE，见方案 §二 引用的 Spitale 2015）。
  - 第 4 层（Mechanism）超出纯计算，需机制实验证据。

### 6. 正式 meta-analysis（per-dataset effect size 汇总）

- **出处原文**：《chat_idea_整理.md》§五：*"Dataset 1 → m6A vs matched controls → effect size₁……meta-analysis → overall effect。但如果你的主要目标只是一个本科生项目，也可以先采用一个更简单的策略：把所有数据合并，然后在分析中使用 transcript-level / dataset-level bootstrap 或 permutation。"*
- **原放弃理由**：本科生项目先采用"合并 + bootstrap"的简化策略。
- **重启依赖与数据缺口**：
  - 方法依赖：per-dataset 独立效应估计 + 固定/随机效应汇总。
  - 数据缺口：需要多于一个"可独立估计效应"的数据集（当前只有 HEK293T 完成；HeLa 完成后即有第二个数据集）。

---

## 第二类：方案明确排除"不作为结题条件"

### 7. lncRNA

- **出处原文**：《ModStruct_本科生研究方案.md》§一：*"本科生主线限定为：人类、mRNA、m6A……"*、*"lncRNA……不作为结题条件。"*
- **原放弃理由**：主线限定为 mRNA。
- **重启依赖与数据缺口**：现有 GLORI 文件为 mRNA（文件名 `*_mRNA_*`）；lncRNA 需要另找 icSHAPE + m6A 数据源并重走审计。

### 8. 因果推断

- **出处原文**：《ModStruct_本科生研究方案.md》§一：*"因果推断……不作为结题条件。"*
- **原放弃理由**：观察性、跨研究数据无法识别干预效应。
- **重启依赖与数据缺口**：需因果图 + 扰动/工具变量数据，同第 5 项。

### 9. 基础模型 / 大模型训练

- **出处原文**：《ModStruct_本科生研究方案.md》§一：*"基础模型训练……不作为结题条件。"*；《ModStruct_通俗解读.md》§十：*"不做实验、不搞大模型。"*
- **原放弃理由**：超出本科范围；"又一个 RNA modification prediction model"，领域拥挤。
- **重启依赖与数据缺口**：需 GPU、大规模语料；与 SMART-m6A 等正面撞车，差异化价值低。

### 10. 数据库 / 网站

- **出处原文**：《ModStruct_本科生研究方案.md》§一：*"数据库网站均不作为结题条件"*；§十三：*"无需制作在线数据库。"*
- **原放弃理由**：与 StructRMDB 竞争，非优先方向。
- **重启依赖与数据缺口**：需 web 开发；内容需建立在多修饰数据（第 1 项）落地之后才有意义。

### 11. 深度学习 / Transformer

- **出处原文**：《ModStruct_本科生研究方案.md》§九：*"主线不需要 PyTorch、Transformer 或 GPU。"*
- **原放弃理由**：简单模型足够回答"信息从哪来"；避免"又一个预测模型"。
- **重启依赖与数据缺口**：需 GPU；已与 SMART-m6A 完全重叠，重启价值取决于能否绑定第 1/2 项的差异化问题。

---

## 第三类：数据 / 技术层面的"选做 / 推迟"

### 12. in vitro icSHAPE

- **出处原文**：《ModStruct_本科生研究方案.md》§三：*"HEK293T 的 in vivo icSHAPE；in vitro 为选做"*；*"GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt.gz，约 7.0 MB，仅在扩展时使用。"*
- **原放弃理由**：选做。
- **重启依赖与数据缺口**：InVitro 文件**未下载**（raw_processed 当前只有 InVivo）；需下载 + 走 03–06 审计。

### 13. Human RNA MaP（U2OS）

- **出处原文**：《ModStruct_本科生研究方案.md》§三：*"Human RNA MaP 可作后续资源，但不进入第一阶段：它主要对应 U2OS……其 DMS-TRAM-seq 工作在本次核查的 PubMed 记录中仍标为预印本，应按该状态引用。"*
- **原放弃理由**：细胞系与修饰样本不匹配；预印本状态。
- **重启依赖与数据缺口**：数据未下载；需按预印本状态引用；U2OS 需匹配的修饰数据。

### 14. GSE50676（PARS）

- **出处原文**：《ModStruct_本科生研究方案.md》§三：*"GSE50676 对应 Wan 等的人类家系 RNA 结构研究，使用 PARS，不能归类成 DMS-seq；本方案不将其混入 icSHAPE 主分析。"*
- **原放弃理由**：方法（PARS）与 icSHAPE 不同，不混入主分析。
- **重启依赖与数据缺口**：数据未下载；PARS 与 icSHAPE 的技术差异需单独建模。

### 15. m6A-SAC-seq 跨技术复核（GSE162356）

- **出处原文**：《ModStruct_本科生研究方案.md》§三：*"跨技术扩展 | GSE162356……HeLa poly(A) m6A-SAC-seq……选做"*；《路线图》§8：*"SAC-seq 跨技术复核"*。
- **原放弃理由**：选做扩展。
- **重启依赖与数据缺口**：数据未下载；需审计 BED 的定量字段、质量门槛与背景位点可用性。

### 16. 401 nt 窗口敏感性

- **出处原文**：《ModStruct_本科生研究方案.md》§七：*"401 nt 的窗口敏感性分析只在有余力时做，并同步扩展序列基线。"*
- **原放弃理由**：有余力才做。
- **重启依赖与数据缺口**：需重新折叠 401 nt 序列（ViennaRNA），并同步扩展序列基线特征；无外部数据缺口。

### 17. 异构体感知分析（isoform-aware）

- **出处原文**：《ModStruct_本科生研究方案.md》§五：*"异构体感知分析不设为本科生必做项。"*
- **原放弃理由**：短读长异构体混合问题，不做必做项；当前用"确定性代表转录本"规则。
- **重启依赖与数据缺口**：需 isoform-aware 建模方法；现有 `06_hek293t_all_isoform_mappings.csv.gz`（33,962 条映射）可复用为起点。

---

## 第四类：分析 / 模型层面的降级（写进方案但实际没跑或降级跑）

### 18. 树模型 HistGradientBoosting

- **出处原文**：《ModStruct_本科生研究方案.md》§九：*"只比较两类模型：Ridge 回归……HistGradientBoostingRegressor……"*；《路线图》§4：*"约 500—2,000 位点……树模型降为选做"*；§11：*"按 Q4/扩展 → 树模型 → 可选图解释依次删减。"*
- **原放弃理由**：数据量不足时降为选做。
- **重启依赖与数据缺口**：方法已备（scikit-learn），无数据缺口；当前 results/models 只有 m0–m4_ridge.joblib，**树模型未跑**。

### 19. SHAP / 置换重要性可解释性

- **出处原文**：《ModStruct_本科生研究方案.md》§九：*"SHAP 不是必做项，也不是因果解释。"*
- **原放弃理由**：非必做；置换会产生不自然的序列—结构组合。
- **重启依赖与数据缺口**：方法已备；主要证据仍应靠重训练消融，SHAP 只作诊断。

### 20. Q4 修饰 / 低修饰背景分类对照

- **出处原文**：《ModStruct_本科生研究方案.md》§十二：*"这部分不是主线的前提。"*；《路线图》§4：*"缺乏有检测能力的低修饰背景 | 取消 Q4，继续连续比例主线。"*
- **原放弃理由**：缺"测过且足够覆盖"的 A 位点背景。
- **重启依赖与数据缺口**：核心数据缺口——需要可靠低修饰背景位点（当前 GLORI 文件是显著位点列表，无完整可测背景）。

### 21. 信息增益 log-loss→bits（信息论模块第三项）

- **出处原文**：《ModStruct_本科生研究方案.md》§十一：*"第三项作为第 9—10 周的扩展。"*（Δ_bits = {LogLoss(M2) − LogLoss(M4)} / ln 2）
- **原放弃理由**：选做扩展。
- **重启依赖与数据缺口**：依赖 Q4 的分类标签或比例分位点分组；方法已备。

### 22. 大规模超参搜索

- **出处原文**：《ModStruct_本科生研究方案.md》§九：*"不做大规模搜索。"*
- **原放弃理由**：避免"为好看而调参"的纪律考量。
- **重启依赖与数据缺口**：无数据缺口；重启需注意不违反"测试前冻结"纪律。

### 23. 结构熵 H_pair 作为独立信息论分析

- **出处原文**：《ModStruct_本科生研究方案.md》§七：*"H_pair | 计算 ensemble 的配对状态不确定性 | 信息论扩展"*；《chat_idea_整理.md》§11.3（structural entropy / "修饰是否偏好结构不确定性高的区域"）。
- **原放弃理由**：降为 M2/M4 的一个特征列，未发展成独立假设检验。
- **重启依赖与数据缺口**：H_pair 已计算（data/final/08a_* 特征里含 pairing-state entropy）；缺的是把"修饰偏好结构动态区"作为可检验假设的正式分析。

---

## 第五类：不是"放弃"但尚未完成（提醒别漏）

### 24. HeLa 复核 Q3（关联复核 + 模型迁移）

- **出处原文**：《ModStruct_本科生研究方案.md》§一：*"Q3：上述结果是否具有跨数据集稳健性？……完整版必做"*；§十（关联复核 + 迁移评估两层次）。
- **状态**：数据已下载（HeLa GLORI `GSM6432595/2596` + `GSE145805_RAW.tar` 均在 data/raw_processed），但 **09 阶段脚本未写、未跑**。
- **重启依赖与数据缺口**：无新数据缺口；需按冻结的 R_flank10、协变量、QC 规则重走 HeLa 建表（09 阶段），并执行关联复核 + 冻结模型迁移。

---

## 汇总表

| # | 类别 | 被放弃项 | 是否被 SMART-m6A 覆盖 | 重启难度 |
|---|---|---|---|---|
| 1 | 一 | 多修饰 structural niches | **否（真空白）** | 高（数据） |
| 2 | 一 | 高维条件互信息 | **否（真空白）** | 中（方法） |
| 3 | 一 | Structural association network | 否 | 中 |
| 4 | 一 | 多数据集 + RMBase | 部分 | 高（数据） |
| 5 | 一 | 因果框架 2–4 层 | 否 | 高（需实验） |
| 6 | 一 | 正式 meta-analysis | 否 | 中 |
| 7 | 二 | lncRNA | 否 | 高（数据） |
| 8 | 二 | 因果推断 | 部分 | 高 |
| 9 | 二 | 基础模型训练 | **是** | 高（且无差异化） |
| 10 | 二 | 数据库/网站 | 部分（StructRMDB） | 中 |
| 11 | 二 | 深度学习/Transformer | **是** | 高（且无差异化） |
| 12 | 三 | in vitro icSHAPE | 否 | 低（未下载） |
| 13 | 三 | Human RNA MaP | 否 | 中（未下载） |
| 14 | 三 | GSE50676（PARS） | 否 | 中（未下载） |
| 15 | 三 | m6A-SAC-seq 跨技术 | 否 | 中（未下载） |
| 16 | 三 | 401 nt 窗口敏感性 | 否 | 低 |
| 17 | 三 | 异构体感知分析 | 否 | 中 |
| 18 | 四 | 树模型 | 部分 | 低 |
| 19 | 四 | SHAP | 部分 | 低 |
| 20 | 四 | Q4 低修饰背景分类 | 部分 | 高（数据） |
| 21 | 四 | log-loss→bits | 部分 | 低 |
| 22 | 四 | 大规模超参搜索 | 部分 | 低 |
| 23 | 四 | H_pair 独立分析 | 否 | 低 |
| 24 | 五 | HeLa Q3 | 部分 | 中（脚本未写） |

## 结束说明

- 本文仅作范围索引，不改变任何已冻结配置或已产出结果。
- 重启任何一项前，须先按原方案 §四 的审计纪律核实"数据缺口"的实际状态（下载与否、字段是否可用），不凭本文的"未下载"标注直接下结论。
- 优先级建议（针对"重启被砍项"）：**第 1 项（多修饰）和第 2 项（条件互信息）是唯二"未被 SMART-m6A 覆盖且可纯计算推进"的差异化方向**，第 24 项（HeLa Q3）是收尾主线、与重启并行不冲突。
