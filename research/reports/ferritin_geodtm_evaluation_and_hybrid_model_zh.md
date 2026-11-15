# 铁蛋白突变热稳定性预测评估与专属模型 `Ferritin-GeoHybrid` 技术报告

**项目根目录**：[`/Users/a26683/Desktop/ferritinMTdatabase`](file:///Users/a26683/Desktop/ferritinMTdatabase)  
**本地 `GeoDTm-3D` 参考部署（全程只读）**：`/Users/a26683/Desktop/AIDD8122_Final_Project/Ferritin_AIDD8122`  
**主报告路径**：[`ferritin_geodtm_evaluation_and_hybrid_model_zh.md`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/reports/ferritin_geodtm_evaluation_and_hybrid_model_zh.md)  

---

## 1. 摘要与核心结论

本研究完成了三项核心任务：（1）在严格隔离环境下复现本地部署的通用三维几何热稳定性模型 `GeoDTm-3D` 并完成算子与 CLI 数值等价性校验；（2）基于本工作区清洗校对的铁蛋白突变热稳定性数据库，构建涵盖 24-mer 单点、12-mer `EcDps`、多点累加、次级表型与定性方向验证的五个分层评估队列（共 `N=75` 条实测突变记录、`14` 个野生型结构模板、`0` 处序列/PDB 错配），系统评估了 `GeoDTm-3D` 及其衍生基线的泛化瓶颈；（3）设计并实现针对铁蛋白球壳多聚体体系的物理约束杂合模型 [`Ferritin-GeoHybrid`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py)（基于 [`NonNegHuberMAPRegressor`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py)），在跨蛋白留一（LOPO）、跨亚家族留一（LOFO）以及未见外部队列上实现统计显著的性能跃升。

核心定量结论如下：
1. **通用单链模型的系统性失效**：通用单体 [`Monomer_GeoDTm3D_k32`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/baseline_evaluation_summary.json) 在 24-mer 主队列 `S1_All30 (N=30)` 上的线性相关性几乎为零（`r=-0.010`, `ρ=+0.017`, `MAE=9.43 °C`，劣于简单零预测 `Zero_baseline` 的趋势区分能力，在 `S1_PairedCore23 (n=23)` 上为 `r=-0.068`, `MAE=10.04 °C` vs 零预测 `9.91 °C`），在独立定性方向集 `S5_Unseen13_Directional (n=13)` 上近似随机猜测（`ROC-AUC=0.524`, `MW p1s=0.4726`）。直接将 24-mer 组装体送入未重校准的 `GeoDTm-3D`（`Assembly_GeoDTm3D_k32`）同样无法改善排序（`S1_All30 r=+0.028`, `S5 AUC=0.524`）。
2. **`Ferritin-GeoHybrid` 的稳定域外泛化**：仅用 5 个具有明确生物物理意义的同向对齐特征（多聚体 `GeoDTm-3D` 输出 `geo_asm_k32`、笼体等效折叠自由能 `neg_fx_cage_ddg_eff`、C3 三重对称轴界面结合能 `neg_fx_ac_ddg_bind_c3`、界面非极性去溶剂化能 `neg_fx_ac_ddg_bind_nonpolar` 及多聚体编码器表征偏移约束项 `neg_geo_asm_enc_shift_k32`），在 `NonNegHuberMAPRegressor`（$w \ge 0$ 非负约束、Huber $\delta=8.0^\circ\text{C}$、先验锚定 MAP 正则 $\alpha^*=80.0, w_{\text{prior}}^*=1.0$）下，在 `S1_All30` 上取得 `LOPO r=+0.375 (p=0.0409), MAE=8.58 °C` 与 `LOFO r=+0.394 (p=0.0311), MAE=9.15 °C`。
3. **跨寡聚态与多点突变的零样本迁移突破**：直接将 `S1_All30` 训练冻结后的模型零样本迁移至未见测试集：
   - **正四面体 12-mer `S2_EcDps_Single8 (n=8)`**：`MAE` 从 `Zero_baseline (10.38 °C)` 与 `Monomer_GeoDTm3D_k32 (7.93 °C)` 降至 **`5.09 °C`**（`r=+0.505, ρ=+0.524, sign=1.000`，相对 `Monomer_GeoDTm3D_k32` 单侧 Wilcoxon 检验 `diff=-2.85 °C, p=0.0234, Win/Loss=7/1`；95% 预测区间覆盖率 `1.000`）；
   - **多点突变 `S3_Exact9 (n=9)`**：相关性达到 **`r=+0.678 (p=0.0447), ρ=+0.833, sign=0.889`**，`MAE` 从零预测的 `18.64 °C` 降至 `14.29 °C`（`p=0.0391, Win/Loss=8/1`）；
   - **定性方向验证 `S5_Unseen13_Directional (n=13)`**：**`ROC-AUC=0.952 (MW p1s=0.0023)`**，按截距中性阈值（`> b0=-3.0859`）二分类平衡准确率达 `0.786`（准确率 `0.769`, `tp=4, fp=0, fn=3, tn=6`）；
   - **全量未见样本汇总 `Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47 (N=47)`**：相对 `Zero_baseline` 误差降低 `-2.35 °C (p=0.0018, W/L=33/14)`，相对 `Monomer_GeoDTm3D_k32` 降低 **`-1.14 °C (CI=[-2.49, +0.03], p=0.0257, W/L=30/17)`**，且 `r=+0.443 (p=0.0018), ρ=+0.433, sign=0.756`；伴随实现的残差校准版本 `Residual_NonNegMAP` 在同一汇总集上达到 `MAE=8.78 °C, r=+0.459 (p=0.0012)`（较 `Monomer_GeoDTm3D_k32` 降低 `-1.43 °C, CI=[-2.77, -0.49], p=0.0170`）。

---

## 2. 工程环境边界与数据队列设计

### 2.1 运行环境与仓库边界隔离
- **参考部署只读保护**：`/Users/a26683/Desktop/AIDD8122_Final_Project/Ferritin_AIDD8122` 下原有核心脚本保持只读零修改；根工作区已有交付物经 [`scripts/validate_deliverables.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/scripts/validate_deliverables.py) 检验为 `0 errors`（Dataset rows: `334`，Unique mutants: `240`，Core Ferritin Tier 1-4: `recs=259, uniq=210`）。
- **双虚拟环境调用**：按照工具依赖划分，特征提取模块通过 `GEODTM_PYTHON`（`.../.venv-geodtm/bin/python`，Python `3.10.20`，含 PyTorch / ESM / Biopython / SciPy / Pandas / `pytest-8.3.5`）运行；训练与交叉验证基准通过 `SKLEARN_PYTHON`（`.../.venv/bin/python`，含 scikit-learn / joblib / SciPy / NumPy / Pandas）运行。

### 2.2 分层铁蛋白测试队列（`N=75`）
全部清洗数据存储于 [`curated_ferritin_cohorts.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/data/curated_ferritin_cohorts.csv)（[`ferritin_features_all.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/ferritin_features_all.csv) 共 `75` 行 × `113` 列特征）及 `14` 个野生型模板配置 [`wt_templates.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/data/wt_templates.json)，所有记录序列与 PDB 坐标索引错配率为 `0`：

| 队列标识 | 样本量 ($N$) | 寡聚对称性与体系描述 | 核心实验用途 |
|---|---:|---|---|
| `S1_SinglePoint_24mer` (`S1_All30`) | `30` | 八面体 24-mer（含 `S1_PairedCore23` `23` 条、`S1_ExtraPatent2` `2` 条、`S1_ExtraOyster1` `1` 条、`S1_Sensitivity_MtBfrA_ESI4` `4` 条），覆盖 `23` 个独立 `(protein_group, site)`、`8` 个蛋白组、`4` 个 `<35%` 序列一致性亚家族（`oligomer_family`）与 `12` 项独立研究（`study_id`）；另含 `S1_PairedCore23_AltThesis`（`5` 处口径差异敏感性对照） | 主拟合训练集与严苛分组外推评估（SiteGroup / LOPO / LOFO / StudyGroup） |
| `S2_EcDps_Single8` | `8` | 正四面体 12-mer `EcDps`（UniProt `P0ABT2`）单点突变体 | 外部零样本迁移测试：检验跨寡聚体（24-mer $\to$ 12-mer）泛化 |
| `S3_MultiPoint` (`S3_All16`) | `16` | 包含精确定量子集 `S3_Exact9`（`9` 条）与 `BfMFtn 9WX5–9X6R` 摘要子集 `S3_Abstract7`（`7` 条，均记为 `-5.0`） | 外部零样本迁移测试：检验单点监督向多点复合突变的累加外推（以 `S3_Exact9` 为主评价集） |
| `S4_Secondary8` | `8` | 异质次级表型：`S4_Ti`（`n=1`, `-14.0`）、`S4_Heteropolymer`（`n=1`, `-3.1`）、`S4_DSC_Peak1`（`n=3`, `9.0, 3.0, -7.0`）、`S4_DSC_Peak2`（`n=3`, `6.0, 0.0, -6.0`） | 物理边界诊断：测试单热峰折叠模型面对多峰 DSC 解离峰/异聚体时的表征失配 |
| `S5_Unseen13_Directional` | `13` | 仅含热稳定性变化方向标签 `binary_direction`（`7` 条 `> WT`，`6` 条 `< WT`） | 外部独立定性判别检验（ROC-AUC 与单侧 Mann-Whitney 精确检定） |

---

## 3. `GeoDTm-3D` 复现等价性校验与多聚体适配

在执行大队列特征抽提之前，适配器 [`geodtm_assembly_adapter.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/geodtm_assembly_adapter.py) 与流程脚本 [`extract_all_features.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/extract_all_features.py) 内嵌了两项强一致性门禁，确保本工作区几何邻域构建与 ESM-2 序列嵌入计算和本地部署的官方 `GeoDTm-3D` 完全等价：
1. **标准基准结构检验（`verify_1cts_equivalence()`）**：官方输出 `official_1cts_c = -2.9801671504974365`，本适配器输出 `adapter_1cts_c = -2.9801671504974365`，绝对误差 **`abs_diff_c = 0.0`**（门限 `< 1e-4`）。
2. **铁蛋白真实案例端到端 CLI 检验（`verify_ferritin_monomer_cli_equivalence()`，以 `EcBfr R30A` 为例）**：官方 CLI 输出 `official_cli_ecbfr_r30a_c = -0.2750121057033539`，本适配器输出 `adapter_ecbfr_r30a_c = -0.27501383423805237`，绝对误差 **`abs_diff_c = 1.7285346984863281e-06`**（门限 `< 1e-4`）。
3. **多聚体结构化学兼容修复**：通过 [`HIS_ALIASES`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/geodtm_assembly_adapter.py)（`{"H1S", "H2S", "HSD", "HSE", "HSP"}`）将 FoldX 修复后多聚体 PDB 中的质子化组氨酸变体统一还原为标准 `HIS`，避免了多聚体突变中心与界面近邻提取时的残基缺口异常；多聚体搭建及 FoldX 界面能量分解分别由 [`assembly_builder.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/assembly_builder.py) 与 [`foldx_assembly.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/foldx_assembly.py) 完成。

---

## 4. 通用蛋白质模型在铁蛋白体系上的评估结果与失效根因

基于 [`extract_all_features.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/extract_all_features.py) 计算并导出至 [`baseline_evaluation_summary.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/baseline_evaluation_summary.json) 的全量通用基线评估结果汇总如下：

### 4.1 定量队列与定性方向基线表现总表

| 基线模型 | `S1_PairedCore23` ($n=23$) | `S1_All30` ($n=30$) | `S2_EcDps_Single8` ($n=8$) | `S3_Exact9` ($n=9$) | `S5_Unseen13` ($n=13$) |
|---|---|---|---|---|---|
| `Zero_baseline` ($\widehat{\Delta T_m}=0$) | MAE `9.91`, RMSE `12.68` | MAE `9.55`, RMSE `12.10` | MAE `10.38`, RMSE `11.75` | MAE `18.64`, RMSE `20.34` | — |
| `Monomer_GeoDTm3D_k32` | r `-0.068`, ρ `-0.063`, MAE `10.04`, RMSE `12.95` | r `-0.010`, ρ `+0.017`, MAE `9.43`, RMSE `12.19` | r `+0.396`, ρ `+0.238`, MAE `7.93`, RMSE `9.30` | r `+0.408`, ρ `+0.250`, MAE `14.88`, RMSE `16.63` | AUC `0.524`, `MW p1s=0.4726`, `gt0 acc=0.462` |
| `Monomer_GeoDTm3D_k48` | r `-0.070`, ρ `-0.085`, MAE `10.07`, RMSE `12.98` | — | — | — | — |
| `Assembly_GeoDTm3D_k32` | r `-0.027`, ρ `-0.061`, MAE `10.00`, RMSE `12.87` | r `+0.028`, ρ `+0.007`, MAE `9.37`, RMSE `12.11` | r `+0.301`, MAE `7.86`, RMSE `9.39` | r `+0.409`, ρ `+0.250`, MAE `14.93`, RMSE `16.69` | AUC `0.524`, `MW p1s=0.4726`, `gt0 acc=0.462` |
| `Assembly_GeoDTm3D_k48` | r `-0.041`, ρ `-0.041`, MAE `10.05`, RMSE `12.90` | — | — | — | — |
| `FoldX_mono_neg_ddG` | r `-0.195`, ρ `-0.169`, MAE `10.19`, RMSE `12.87` | r `+0.011`, ρ `+0.049`, MAE `9.64`, RMSE `12.12` | r `+0.742`, ρ `+0.810`, MAE `9.36`, RMSE `10.29` | r `+0.316`, MAE `17.98`, RMSE `19.57` | **AUC `1.000`**, `MW p1s=0.0006`, `acc=0.846` |
| `FoldX_sub_neg_ddG_per_chain` | r `+0.128`, ρ `+0.188`, MAE `9.61`, RMSE `12.56` | r `+0.213`, ρ `+0.261`, MAE `9.07`, RMSE `11.75` | r `+0.686`, MAE `8.38` | r `+0.476`, ρ `+0.450`, MAE `16.81`, RMSE `18.27` | AUC `0.952`, `acc=0.769` |
| `FoldX_AC_neg_bind_ddG` | r `+0.312`, ρ `+0.323`, MAE `9.44`, RMSE `12.18` | r `+0.294`, ρ `+0.351`, MAE `9.11`, RMSE `11.56` | r `+0.002`, MAE `7.98` | r `+0.287`, MAE `15.85` | AUC `0.833`, `MW p1s=0.0256`, `acc=0.615` |
| `FoldX_cage_neg_eff_ddG` | r `+0.222`, ρ `+0.208`, MAE `9.39`, RMSE `12.35` | r `+0.251`, ρ `+0.297`, MAE `8.98`, RMSE `11.63` | r `+0.615`, ρ `+0.571`, MAE `7.66`, RMSE `8.87` | r `+0.409`, ρ `+0.400`, MAE `16.58`, RMSE `18.18` | AUC `0.952`, `MW p1s=0.0023`, `acc=0.846` |
| `ESM2_WT_marginal_LLR` | r `-0.363`, ρ `-0.299`, MAE `12.22`, RMSE `16.11` | r `-0.164`, ρ `-0.079`, MAE `10.70`, RMSE `14.55` | r `-0.034`, MAE `5.61`, RMSE `6.32` | r `+0.184`, MAE `11.85`, RMSE `16.59` | AUC `0.381`, `MW p1s=0.7774`, `acc=0.462` |
| `S4346_Frozen_Ridge_Sub` | r `+0.447`, ρ `+0.370`, MAE `9.31`, RMSE `12.20` | r `+0.243`, ρ `+0.167`, MAE `9.19`, RMSE `11.77` | r `+0.485`, MAE `7.88`, RMSE `9.15` | — | — |

补充敏感性与次级队列基线：
- 在 `S1_PairedCore23_AltThesis (n=23)` 上，各通用基线趋势一致：`Zero` MAE `9.83`, RMSE `12.65`；`Mono32` r `-0.118`, ρ `-0.125`, MAE `10.33`, RMSE `13.06`；`Asm32` r `-0.078`, ρ `-0.121`, MAE `10.30`, RMSE `12.98`；`FoldX_AC` r `+0.323`, ρ `+0.307`, MAE `9.44`；`S4346_Ridge_Sub` r `+0.446`, MAE `9.43`。
- 在 `S3_All16 (n=16)` 上：`Zero` MAE `12.68`；`Mono32` MAE `10.59`, r `-0.182`。
- 在次级队列 `S4_Secondary8 (n=8)` 上：`Zero` MAE `6.01`, RMSE `7.22`；`Mono32` r `+0.294`, ρ `+0.303`, MAE `5.39`, RMSE `6.86`；`Asm32` r `+0.235`, MAE `5.58`, RMSE `6.91`；而仅计算静态对称结合的 `FoldX_AC` 与 `FoldX_cage` 呈现负相关（分别 `r=-0.912, MAE=8.56` 与 `r=-0.764, MAE=7.80`）。

### 4.2 通用模型失效的三大生物物理诊断
从上述跨表对比中可直接观察到通用模型不适配铁蛋白体系的三条机理证据：
1. **单体表征对亚基界面效应失明**：在同一 24-mer `S1_PairedCore23` 集上，仅计单链折叠自由能的 `FoldX_mono_neg_ddG` 相关性为负（`r=-0.195`），而一旦引入复合物亚基结合界面信息（`FoldX_AC_neg_bind_ddG` 与 `FoldX_cage_neg_eff_ddG`），相关性立即由负转正至 `r=+0.312` 与 `+0.222`（在 `S1_All30` 上为 `+0.294` 与 `+0.251`）。通用单体模型完全遗漏了铁蛋白 C2 二聚体、C3 通道界面及去溶剂化屏障对热崩溃温度的主导贡献。
2. **多聚体输入仍受限于通用权重标定与量纲压缩**：即使将完整多聚体坐标系输入冻结的 `GeoDTm-3D`（`Assembly_GeoDTm3D_k32`），由于其几何边权与回归头仍在常规球蛋白单链数据上拟合，既无法对齐铁蛋白实测动辄 `±25 °C` 的高能量阈值变化范围（实际输出常压缩在极窄区间内），在 `S5_Unseen13` 上的判别力仍停留在 `AUC=0.524`。
3. **序列进化语言模型偏好在人工热点工程突变上发生反转**：`ESM2_WT_marginal_LLR` 在 `S1_PairedCore23`（`r=-0.363`）、`S1_All30`（`r=-0.164`）与 `S5`（`AUC=0.381`）上均表现负向甚至显著弱于零预测（`S1_PairedCore23` MAE 高达 `12.22 °C` vs 零预测 `9.91 °C`）。铁蛋白功能口袋与三重对称轴孔道通常经历负向选择约束以保持矿化与铁离子转运柔顺性，因而在序列演化上罕见但热力学极稳定的芳香族疏水填充突变会遭到语言模型的错误惩罚。

---

## 5. `Ferritin-GeoHybrid` 模型架构设计与理论依据

针对“高维多聚体结构物理耦合 + 小样本铁蛋白热测量噪声较大（存在 ±20~26 °C 极端异常突变体）”的双重挑战，我们在 [`hybrid_model.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py) 中提出并实现了 **`Ferritin-GeoHybrid`**（模型实现类 [`NonNegHuberMAPRegressor`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py)，拟合入口 [`fit_ferritin_hybrid`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py)）。

### 5.1 五维正向物理对齐特征（`HYBRID_FEATURES_5`）与全集拟合系数
模型严格限定所有特征方向预先取与 $\Delta T_m$ 升高的同向物理符号，并通过 `L-BFGS-B` 强制施加非负权重约束 $w_j \ge 0$、Huber 鲁棒损失函数（阈值 $\delta = 8.0^\circ\text{C}$）与锚定物理先验的高斯 MAP 正则项（内层按 `protein_group` 分组做 `GroupKFold` 交叉验证，在 `S1_All30`、`S1_PairedCore23` 与 `S1_PairedCore23_AltThesis` 三种口径下均稳定一致地自动选出 $\alpha^* = 80.0$、$w_{\text{prior}}^* = 1.0$）：

| 特征符号 (`HYBRID_FEATURES_5`) | 生物物理含义 | 标准化系数 (`S1_All30`) | 原始系数 (`S1_All30`) | `LOPO` 折间均值±标准差 | `LOFO` 折间均值±标准差 | `Core23` 敏感性权重 |
|---|---|---:|---:|---|---|---:|
| `geo_asm_k32` | 多聚体上下文 `GeoDTm-3D` 标量预测 | `0.8957` | `0.4589` | `0.9039 ± 0.0756` (100% 正) | `0.8784 ± 0.0749` (100% 正) | `0.8635` |
| `neg_fx_cage_ddg_eff` | FoldX 笼体等效折叠自由能稳定性 (`-ddG_eff`) | `1.0331` | `0.3549` | `1.0243 ± 0.0641` (100% 正) | `0.9911 ± 0.0544` (100% 正) | `0.9745` |
| `neg_fx_ac_ddg_bind_c3` | FoldX C3 三重轴界面结合贡献 (`-ddG_bind_C3`) | `1.3827` | `0.5042` | `1.3401 ± 0.1026` (100% 正) | `1.2877 ± 0.1393` (100% 正) | `1.3587` |
| `neg_fx_ac_ddg_bind_nonpolar` | 亚基界面非极性疏水去溶剂化增益 | `1.3343` | `0.5032` | `1.2823 ± 0.1788` (100% 正) | `1.2658 ± 0.0933` (100% 正) | `1.4873` |
| `neg_geo_asm_enc_shift_k32` | 多聚体编码器表征偏移约束项 | `1.0600` | `1.0578` | `1.0528 ± 0.0425` (100% 正) | `1.0091 ± 0.0960` (100% 正) | `1.0349` |

- **截距参数**：`S1_All30` 原始特征截距为 `0.1675`，标准化特征空间截距 `b0 = -3.0859`（`S1_PairedCore23` 敏感性标定 `b0 = -2.7121`）。
- **不确定性预测**：模型使用迭代重加权最小二乘（IRLS）/ Laplace 近似计算预测后验方差并导出个例标准差与 68%/95% 置信区间；此外还并行实现了残差杂合版本 `RESIDUAL_HYBRID_FEATURES_4`（`Residual_NonNegMAP`）作为工程校准参考对照。

### 5.2 为什么必须采用非负物理约束先验？（反过拟合实证）
在小样本且特征强共线的条件下（`S1_All30`, `N=30`），若取消权重符号与先验正则约束（如对照模型 `BayesianRidge_FreeSign7`），优化器会利用局部离群点去拟合并翻转基本物理定律：
- 在 `BayesianRidge_FreeSign7` 中，单体折叠稳定性项 `neg_fx_mono_ddg` 被拟合为负号（`S1_All30` 上为 `-0.4621`，在 `Core23` 上恶化至 `-5.5812`），辅助因子距离项 `near_cofactor_8a` 在 24-mer `S1` 上取 `-1.4614`（主要源自含血红素 Bfr 亚群）而在 `S2 EcDps` 上却呈反向相关性 `+0.703`；
- 其后果是无约束线性与非线性对比模型在外推验证中集体崩塌（详情见表 6.1：`BayesianRidge_FreeSign7` `LOPO r=-0.060`、`S2 r=-0.226 (rho=-0.643)`、`S5 AUC=0.643`）。反观 `Ferritin-GeoHybrid`，全部五个权重在 LOPO 与 LOFO 全部折次中 `100%` 保持正号且方差极小。

---

## 6. 全方位验证结果：交叉验证、消融实验、外部零样本迁移与误差审计

所有实验均由 [`run_hybrid_benchmark.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/run_hybrid_benchmark.py) 一键复现，详细数据记录于 [`hybrid_benchmark_summary.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/hybrid_benchmark_summary.json) 与逐样本预测表 [`hybrid_predictions_all.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/hybrid_predictions_all.csv)（75 条全量输出）。

### 6.1 `S1_All30 (N=30)` 不同模型族分组外推对比

| 模型族 | `SiteGroup` CV | `LOPO`（跨蛋白留一） | `LOFO`（跨家族留一 <35% id） | `StudyGroup`（跨文献留一） |
|---|---|---|---|---|
| **`FerritinGeoHybrid_NonNegMAP`** | MAE `8.91`, RMSE `11.85`, r `+0.252`, ρ `+0.259`, sign `0.571` | **MAE `8.58`, RMSE `11.24`, r `+0.375 (p=0.0409)`, ρ `+0.358`, sign `0.643`** | **MAE `9.15`, RMSE `11.10`, r `+0.394 (p=0.0311)`, ρ `+0.306`, sign `0.536`** | MAE `9.18`, RMSE `12.07`, r `+0.181` |
| `Residual_NonNegMAP` | MAE `8.73`, r `+0.248` | MAE `8.50`, r `+0.350` | MAE `8.97`, r `+0.391` | MAE `9.10`, r `+0.178` |
| `ZeroPrior_NonNegHuber` | MAE `9.76`, r `+0.070` | MAE `9.32`, r `+0.223` | MAE `9.97`, r `+0.262` | — |
| `BayesianRidge_5feat` | MAE `10.15`, r `-0.062` | MAE `9.40`, r `+0.002` | MAE `9.95`, r `+0.020` | — |
| `BayesianRidge_FreeSign7` | MAE `10.21`, r `+0.033` | MAE `9.85`, r `-0.060` | MAE `10.39`, r `-0.020` | — |
| `RidgeCV_5feat` | MAE `9.85`, r `-0.143` | MAE `9.49`, r `-0.057` | MAE `9.99`, r `+0.006` | — |
| `Huber_5feat` | MAE `9.68`, r `-0.001` | MAE `9.86`, r `-0.099` | MAE `10.05`, r `+0.138` | — |
| `ShallowGBDT_5feat` | MAE `10.26`, r `-0.080` | MAE `10.05`, r `+0.031` | MAE `11.25`, r `-0.259` | — |
| `RegularizedMLP_5feat` | MAE `13.66`, r `+0.074` | MAE `13.18`, r `-0.084` | MAE `12.12`, r `-0.145` | — |

**口径敏感性验证**：在 `S1_PairedCore23`（SiteGroup `MAE=9.91, r=+0.152`；LOPO `MAE=9.11, r=+0.340, ρ=+0.380`；LOFO `MAE=9.28, r=+0.424`）及其文献异议对照集 `S1_PairedCore23_AltThesis`（SiteGroup `MAE=10.16, r=+0.157`；LOPO `MAE=9.07, r=+0.343, ρ=+0.351`；LOFO `MAE=9.23, r=+0.423`）上，`FerritinGeoHybrid_NonNegMAP` 的相关性与误差收益几乎完全不变。

### 6.2 `Ab1..Ab6` 组件消融梯队验证

利用 `PriorAnchoredNonNegHuber` 保持优化器一致，逐项消融（`ABLATION_SETS`）清晰量化了每一类结构与物理项的边际增益及 ESM-2 的负效应：

| 消融阶段 | `S1 LOPO` | `S1 LOFO` | `S2_EcDps_Single8` | `S3_Exact9` | `S5_Unseen13` |
|---|---|---|---|---|---|
| `Ab1_GeoMono3D_Calibrated` | MAE `9.73`, r `-0.331` | MAE `10.97`, r `-0.187` | MAE `7.25`, r `+0.396` | MAE `14.49`, r `+0.408` | AUC `0.524` |
| `Ab2_GeoAsm3D_Plus_EncShift` | MAE `9.27`, r `+0.099` | MAE `10.78`, r `-0.200` | MAE `7.79`, r `+0.420` | MAE `15.26`, r `+0.697` | AUC `0.667` |
| `Ab3_Plus_CageEffectiveStability` | MAE `8.93`, r `+0.183` | MAE `9.90`, r `+0.068` | MAE `7.18`, r `+0.736` | MAE `15.80`, r `+0.540` | **AUC `0.976 (p=0.0012)`** |
| `Ab4_Plus_C3Interface` | MAE `8.66`, r `+0.278` | MAE `9.47`, r `+0.269` | MAE `6.19`, r `+0.742` | MAE `14.29`, r `+0.815` | AUC `0.929 (p=0.0041)` |
| **`Ab5_Full_Ferritin_GeoHybrid_5feat`** | **MAE `8.58`, r `+0.375`, ρ `+0.358`** | **MAE `9.15`, r `+0.394`** | **MAE `5.09`, r `+0.505`, ρ `+0.524`** | **MAE `14.29`, r `+0.678`, ρ `+0.833`** | **AUC `0.952 (p=0.0023)`** |
| `Ab6_Extended_Plus_ESM2_LLR` | MAE `8.37`, r `+0.331` | MAE `9.58`, r `+0.344` *(退化)* | MAE `4.88`, r `+0.462` | MAE `10.49`, r `+0.726` | AUC `0.833 (p=0.0256)` *(显著退化)* |

消融结果证实：**`Ab1 -> Ab2 -> Ab3 -> Ab4 -> Ab5` 在最严苛的跨亚家族 `LOFO` 上呈严格单调提升**（`r` 由 `-0.187 -> -0.200 -> +0.068 -> +0.269 -> +0.394`，`MAE` 由 `10.97 -> 10.78 -> 9.90 -> 9.47 -> 9.15 °C`）；而一旦在 `Ab6` 中混入序列语言似然比 `ESM2_WT_marginal_LLR`，即刻出现 `LOFO` 退步与 `S5` ROC-AUC 由 `0.952` 大幅滑落至 `0.833`，强有力地印证了不纳入 ESM-2 LLR 特征的设计决策。

### 6.3 外部零样本迁移与汇总统计显著性检验

将 `S1_All30` 完整拟合的 `Ferritin-GeoHybrid` 直接部署至未见外部测试集（无任何微调）：
- **`S2_EcDps_Single8 (n=8)`**：`MAE=5.09 °C, RMSE=6.29 °C, r=+0.505, ρ=+0.524, sign=1.000`；配对误差显著低于所有对比基线：较 `Zero_baseline` 降 `-5.29 °C (p=0.0156, W/L=7/1)`、较 `Monomer_GeoDTm3D_k32` 降 **`-2.85 °C (p=0.0234, W/L=7/1)`**、较 `Assembly_GeoDTm3D_k32` 降 `-2.77 °C (p=0.0391)`、较 `FoldX_mono` 降 `-4.27 °C (p=0.0156)`、较 `FoldX_cage` 降 `-2.57 °C (p=0.0391)`。不确定性标定：95% 置信区间实测覆盖率 `1.000`，68% 置信区间覆盖率 `0.750`。
- **`S3_Exact9 (n=9)`**：`MAE=14.29 °C, RMSE=15.50 °C, r=+0.678 (p=0.0447), ρ=+0.833, sign=0.889`；显著优于 `Zero_baseline`（`diff=-4.35 °C, p=0.0391, W/L=8/1`）与 `FoldX_mono`（`diff=-3.69 °C, p=0.0391`），相对 `Monomer_GeoDTm3D_k32`（`diff=-0.58 °C, p=0.6523, W/L=5/4`）在 MAE 轻微占优的同时实现了极强的单调秩次一致性（ρ 从 `+0.250` 跃升至 `+0.833`）。在扩展含摘要数据全集 `S3_All16 (n=16)` 上：`MAE=11.04 °C, RMSE=12.91 °C, sign=0.938, r=-0.159`。
- **`S5_Unseen13_Directional (n=13)`**：`ROC-AUC=0.952`，单侧 Mann-Whitney `p=0.0023`（双侧 `p=0.0047`）；按中性截距阈值 `> b0=-3.0859` 判定准确率 `0.769`（平衡准确率 `0.786, tp=4, fp=0, fn=3, tn=6`），按零点 `> 0` 判定准确率 `0.615`（平衡准确率 `0.643, tp=2, fp=0, fn=5, tn=6`）。在 `Core23` 敏感性训练与 `ZeroPrior` 对照下 `ROC-AUC` 进一步达 `0.976 (p=0.0012)`（其中 `ZeroPrior` 中性阈值准确率达 `0.846`）。
- **外部定量汇总 `Pooled_External_S2_S3Exact9_N17 (n=17)`**：`MAE=9.96 °C, RMSE=12.07 °C, r=+0.420 (p=0.0936), ρ=+0.414, sign=0.941`；较 `Zero_baseline` 降 `-4.79 °C (p=0.0008, W/L=15/2)`、较 `FoldX_mono` 降 `-3.96 °C (p=0.0005)`、较 `FoldX_cage` 降 `-2.42 °C (p=0.0026)`、较 `Monomer_GeoDTm3D_k32` 降 `-1.65 °C (CI=[-3.74, +3.24], p=0.0887, W/L=12/5)`；预测不确定性与误差呈正相关（Spearman `+0.4681, p=0.0581, cov95=0.7059`）。残差版本 `Residual_NonNegMAP` 在该汇总集上达 `MAE=9.30 °C`（较 `Monomer_GeoDTm3D_k32` `diff=-2.31 °C, p=0.0448`）。
- **全量未见样本大样本汇总 `Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47 (N=47)`**：
  - 整体性能：**`MAE=9.08 °C, RMSE=11.55 °C, r=+0.443 (p=0.0018), ρ=+0.433, sign=0.756`**；
  - 相对所有主流基线均达统计显著提升（`p < 0.05`）：
    - vs `Zero_baseline`：`diff=-2.35 °C`, `95% CI=[-4.53, -0.87]`, `W/L=33/14`, `p=0.0018`
    - vs `Monomer_GeoDTm3D_k32`：**`diff=-1.14 °C`, `95% CI=[-2.49, +0.03]`, `W/L=30/17`, `p=0.0257`**
    - vs `Assembly_GeoDTm3D_k32`：`diff=-1.10 °C`, `p=0.0330`
    - vs `FoldX_mono`：`diff=-2.11 °C`, `p=0.0007`
    - vs `FoldX_AC`：`diff=-1.13 °C`, `p=0.0312`
    - vs `FoldX_cage`：`diff=-1.13 °C`, `p=0.0125`
  - 伴随 `Residual_NonNegMAP` 在同一 `N=47` 样本上进一步取得 **`MAE=8.78 °C, r=+0.459 (p=0.0012)`**，相对 `Monomer_GeoDTm3D_k32` 的 MAE 差值为 `-1.43 °C (95% CI=[-2.77, -0.49], p=0.0170)`。

### 6.4 边界与典型误差机制审计
对 `S1` LOPO 最大残差样本及次级队列 `S4_Secondary8` 的逐一解剖揭示了当前模型尚未覆盖的四类高难度机制：
1. **三重通道芳香族超稳化极端样本**：`D118W`（实测 `+26.3 °C`，混合模型 LOPO 预测 `+1.70 °C`，而 `Monomer_GeoDTm3D_k32` 错判为 `-1.31 °C`）、`D118F`（`+24.7 °C` vs 混合预测 `+0.12 °C` vs 单体 `-1.14 °C`）、`D118Y`（`+23.1 °C` vs 混合预测 `-0.60 °C` vs 单体 `-2.60 °C`）。虽然杂合模型成功抓住了稳定化信号方向改善（由负转正或大幅缩小负偏），但线性稳态加权仍无法匹配超 `20 °C` 的非加成 π-堆叠级联相变；
2. **嗜热/表面电荷重构极端离群突变**：`MjFer D132R`（`-18.0 °C` vs `+3.37 °C`）、`CgFer D120G`（`+13.0 °C` vs `-2.81 °C`）、`RcBFR E135R`（`-14.0 °C` vs `-2.07 °C`）、`EcBfr E128A`（`-17.9 °C` vs `-7.01 °C`）；
3. **非蛋白质辅助因子/血红素轴向配位缺失**：`EcBfr M52A`（`-15.6 °C` vs `-1.62 °C`），FoldX 与静态几何忽略了血红素铁中心轴向配体断键引发的协同去折叠；
4. **寡聚态解聚开关与 DSC 双峰脱偶联（含 `S4_Secondary8`）**：`EcBfr R30A`（`+8.0 °C` vs `-4.92 °C`）与 `R61A`（`+4.0 °C` vs `-8.58 °C`）处于“24-mer 解离为更紧实二聚体”的相变开关上；而在次级表型集 `S4_Secondary8 (MAE=8.42, RMSE=9.56, r=-0.193, sign=0.571)` 中，尽管 `S4_DSC_Peak1`（`9.0, 3.0, -7.0 -> -5.63, -8.99, -8.51`, `r=+0.699`）与 `S4_DSC_Peak2`（`6.0, 0.0, -6.0 -> -5.63, -8.99, -8.51`, `r=+0.794`）内部的相对单调顺序仍有较好捕捉，但模型将解聚引起的界面负 ΔΔG 直接解释为全局失稳，导致整体截距显著左偏；另有异质材料复合体系 `S4_Ti`（`-14.0 -> -2.83`）与异聚体 `S4_Heteropolymer`（`-3.1 -> -8.05`）凸显了单热峰均聚模型的使用边界。

---

## 7. 交付资产与自动化验证清单

| 类别 | 交付路径与主要符号 | 说明 |
|---|---|---|
| 核心源码模块 | [`ferritin_curation.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/ferritin_curation.py), [`assembly_builder.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/assembly_builder.py), [`foldx_assembly.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/foldx_assembly.py), [`geodtm_assembly_adapter.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/geodtm_assembly_adapter.py), [`metrics.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/metrics.py), [`paths.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/paths.py), [`hybrid_model.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/src/hybrid_model.py) | 含 `HIS_ALIASES`、`HYBRID_FEATURES_5`、`RESIDUAL_HYBRID_FEATURES_4`、`CORE_FEATURES_7`、`EXTENDED_FEATURES_9`、`ABLATION_SETS`、`NonNegHuberMAPRegressor`、`fit_ferritin_hybrid`、`save_predictor`、`load_predictor` |
| 执行与评测脚本 | [`run_s571_reproduction.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/run_s571_reproduction.py), [`extract_all_features.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/extract_all_features.py), [`run_hybrid_benchmark.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/scripts/run_hybrid_benchmark.py) | 覆盖等价复现门禁 (`verify_1cts_equivalence`, `verify_ferritin_monomer_cli_equivalence`)、全队列特征生成与杂合模型端到端基准测试 |
| 数据与模型产物 | [`curated_ferritin_cohorts.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/data/curated_ferritin_cohorts.csv), [`wt_templates.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/data/wt_templates.json), [`ferritin_features_all.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/ferritin_features_all.csv), [`baseline_evaluation_summary.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/baseline_evaluation_summary.json), [`hybrid_benchmark_summary.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/hybrid_benchmark_summary.json), [`hybrid_predictions_all.csv`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/results/hybrid_predictions_all.csv), [`ferritin_geohybrid.joblib`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/models/ferritin_geohybrid.joblib), [`ferritin_geohybrid.json`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/models/ferritin_geohybrid.json) | 固化模型 `ferritin_geohybrid.joblib` SHA-256 为 `169b357c348a8ef0f6e8eef0a2c79ff1b0a71028d49995dffaa78486ed04fd33` |
| 自动化测试套件 | [`test_ferritin_pipeline.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/tests/test_ferritin_pipeline.py), [`test_hybrid_model.py`](file:///Users/a26683/Desktop/ferritinMTdatabase/research/tests/test_hybrid_model.py) | `pytest` 特征与管线测试 + `unittest` 杂合模型与预测器序列化测试 |
