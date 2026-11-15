# 铁蛋白（Ferritin）及其超家族突变体温度/热稳定性全网检索汇总报告

> **核心结论（双口径统计 Dual-Count Reconciliation）**：主流蛋白质热稳定性结构化数据库（FireProtDB、ThermoMutDB、GeoStab、ProteinGym、Tsuboyama、MaveDB）中**几近空白（共 0 条天然铁蛋白点突变 Tm 数据，仅 ProThermDB 含 4 条 EBNA1 融合蛋白 DSC Tm 数据且被误标为 wild-type）**；真正可靠的铁蛋白突变热稳定性数据几乎全部散落在**开放获取文献、学位论文（如 NTU 硕士/博士论文 DR-NTU 48015）与付费期刊摘要/全文**中。

---

## 一、 总览与双口径统计核对表（WT 对照严格分离）

| 统计层级 / 口径 | 条件级实验与计算记录数<br>(`ferritin_mutant_thermostability_dataset.csv`) | 去重唯一突变体数<br>(`unique_mutants_summary.csv`) | 说明 |
| :--- | :---: | :---: | :--- |
| **突变体总计（含全部 5 类家族 + 计算预测）** | **292** | **240** | 不含 WT 基线对照 |
| ↳ 其中：**核心铁蛋白家族（Tier 1–4：典型 24-mer、Bfr、Dps、工程融合）** | **259** | **210** | 剔除 Tier 5 扩展超家族（Encapsulin/Artemin/Dodecin/YciF） |
| **湿实验突变体总计（Tier A + B + C，剔除纯计算）** | **174** | **147** | 含定量 Tm/Ti/Tagg/Tdiss、热耐受阈值、定性排序 |
| ↳ 其中：**核心铁蛋白湿实验（Tier 1–4，Tier A–C）** | **141** | **117** | 最严格的核心铁蛋白湿实验数据集 |
| **定量温度测量突变体（Tier A：DSC / CD / nanoDSF / Tycho / DSF / Rayleigh 等）** | **102**（其中精确单值数字 **88**） | **80** | 包括 ΔTm / Tdiss / Tagg_onset / 区间组报告 |
| **纯计算预测条目（Tier D：FoldX ΔΔG + Robetta 虚拟丙氨酸扫描）** | **118** | **93**（纯计算独占；另 **20** 个已在湿实验表重复验证） | 来自 DR-NTU 48015 表 2.3/3.1/4.1 |
| **单独隔离的野生型基线对照（`record_role = WT_baseline`）** | **42** | —（不计入突变体） | pH 扫描、亚基同源寡聚体、FireProtDB P02792 7 条 WT Tm 等 |

---

## 二、 按蛋白家族层级（A1 五大类独立标注）细分统计

| 蛋白家族层级 (`family_tier`) | 突变体全部记录 | 湿实验突变体记录 | 全部唯一突变体 | 湿实验唯一突变体 |
| :--- | :---: | :---: | :---: | :---: |
| `Family_1_Canonical_24mer_Ferritin` | 46 | 46 | 44 | 44 |
| `Family_2_Bacterioferritin_Bfr` | 108 | 44 | 82 | 34 |
| `Family_3_Miniferritin_Dps` | 68 | 14 | 59 | 14 |
| `Family_4_Engineered_Fusion_Display_Ferritin` | 37 | 37 | 25 | 25 |
| `Family_5_Extended_Superfamily_And_Nanocages` | 33 | 33 | 30 | 30 |
| **合计** | **292** | **174** | **240** | **147** |

---

## 三、 按证据质量层级与突变类型拆分

### 3.1 按证据层级 (`evidence_tier`)

| 证据层级 | 条件级记录数 | 唯一突变体数（按最高证据归口） | 数据内涵 |
| :--- | :---: | :---: | :--- |
| `A_WetLab_Quantitative` | 102 | 80 | 定量熔解/解聚/聚集温度（Tm, Ti, Tdiss, Tagg_onset, T_sub_CD, 明确 ΔTm） |
| `B_WetLab_Thermal_Survival` | 40 | 38 | 热澄清 / 耐热存活温度阈值（如 45–100 °C 处理仍保持可溶/笼结构完整） |
| `C_WetLab_Qualitative_Rank` | 32 | 29 | 定性比较或热稳定性排序（>WT / <WT / ~WT / 75 °C 溶解度排名） |
| `D_Computational_Prediction` | 118 | 93 | 计算能量预测（FoldX ΔΔG、Robetta 虚拟 Ala 扫描 ΔΔG_bind） |

### 3.2 按突变类型 (`mutation_category`)

| 突变类型 | 条件级记录数 | 唯一突变体数 |
| :--- | :---: | :---: |
| `Point_Substitution` | 185 | 149 |
| `Multi_Point_Combinatorial` | 43 | 40 |
| `Deletion_Truncation_Frameshift_Swap` | 28 | 21 |
| `Fusion_Display_Heteropolymer` | 36 | 30 |

---

## 四、 各大数据库与文献引擎检索审计汇总（共 17 个库/引擎）

> 完整明细见 [database_and_literature_source_inventory.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/database_and_literature_source_inventory.csv)。

- **结构化稳定性 / DMS / 互作数据库（零值已充分验证）**：
  - `FireProtDB`（2025-09-20 全量 4.83 GB SQL 转储，41,146,172 行）：仅含 `seq_id=8955`（人 `P02792` FTL）**7 条 WT Tm**（与 `4 条 WT ΔG`，来自 `PMID:9681870`），`mutant_id IS NOT NULL` 为 **0**。
  - `ThermoMutDB`（13,337 条 JSON）、`GeoStab`（S4346/S571/S669/S8754 共 14,340 行）、`ProteinGym`（217 substitutions + 66 indels，`A0A1I9GEU1_NEIME` 经核验为 PilE 菌毛素假阳性）、`Tsuboyama Mega-scale`（9,728 条）、`MaveDB`（2,161 条）：全部 **0 条** 铁蛋白突变数据。
  - `SKEMPI 2.0`（7,084 数据行）：含 7 行 `4E6K_AB_G`（铜绿假单胞菌 `BfrB–Bfd`），仅测定 295 K 结合亲和力 `Kd`（SPR），无热变性 `Tm`。
  - `ProThermDB`（缓存检索结果）：含 13+4=17 条相关行（全部为 `WT` 标签），其中 `NO=170, 171, 173, 174` 经检查 `REMARKS` 实为 **人 H-铁蛋白 N/C 端融合 EBV EBNA1（`E1F1`/`F1E1`）的合成构建体 DSC Tm（`70.5/80.1 °C` 与 `63.6/84.6 °C`，`PMID:33825471`）**，全库 ~31k 批量导出列为待补全项。
  - `UniProtKB`（68 个条目 JSON）：含 1 条定性热稳定性下降记录（`P9WNE5` 结核分枝杆菌 `BfrB Δ168–181`，`PMID:21494619`）及多条寡聚体稳定性注释。
  - `RCSB PDB`：检出 **461 个突变铁蛋白家族晶体/冷冻电镜结构**（其中 **38 个** 为 `To Be Published` 或无 PMID/DOI）。
- **文献检索母集与全文深度核验**：
  - `Europe PMC` 10 组检索去重后 **7,083 篇**（原始命中 18,384 条）；`OpenAlex` 10 组去重 **6,365 篇**；`PubMed XML` 11 批去重 **1,856 PMIDs**；`Crossref` 抽样 200 条（含中文食品科学期刊检索）。
  - 完成 **269 篇不重复 PMC 全文**（231 XML + 16 NCBI XML + 22 HTML）、**18 篇 OA PDF/学位论文**（含新加坡南洋理工大学 `DR-NTU 48015` 全文及解码版 `pdf_1007_decoded.txt`）与 **105 篇非 PMC 热稳定性 PubMed 摘要**的逐篇核验；对 260 个 DOI 执行 `Unpaywall` OA 检测（**127 个 OA / 133 个闭源付费**）。

---

## 五、 方法学陷阱与关键数据质量警示（读库必读）

1. **`ΔTm` 符号方向反转陷阱**：Yu Zhang 学位论文（`DR-NTU 48015` 表 2.7/3.3/4.3）及 `PMC6149950`（表 2）采用 **`ΔTm = Tm_WT - Tm_mut`** 定义（正值代表去稳定、负值代表更稳定！例如 `EcBfr R30A` `72.2 °C > WT 64.8 °C` 原文记作 `-7.4`）。本库同时保留原文 `reported_delta_value` 并统一换算 **`normalized_delta_Tm_C = mut - WT`**（正值=更耐热）。
2. **跨源交叉印证与不重复计数**：`PMC2852946`（JBC 2010 `EcBfr`）与 `PMC3267954`（Protein Sci 2011 `EcDps`）在 PMC 中以图片形式嵌表，全文检索不到数值；同实验室 `DR-NTU 48015` 完整保留了表 2.7 与表 4.3 精确数字。已通过 `cross_ref_duplicate_group` 关联并在唯一突变体层面严格去重。
3. **同一突变体跨文献差异**：`EcBfr D56R/E60H` 在 `DR-NTU 48015` 表 3.3 中 `Tm=60.8 °C`（`ΔTm = -7.4 °C`），而同组后续 `J. Biochem. 2015` 论文摘要（`PMID:26115686`）报道为比 WT 低 `16 °C`（体系/缓冲或测定条件不同），两条均已单独入库标注。
4. **测定体系高度异质，严禁直接混池比较绝对 `Tm`**：例如大豆 `SSFH-2` 与 `HuHF-∇EP`（`PMC9069879`）使用**冻干固体粉末 DSC（`HuHF WT = 82 °C`）**；海莴苣 `UpFER`（`PMC3403414`）在 **1.5 M 盐酸胍（GdnHCl）** 中测 DSC；人 `H/L 异源多聚体`（`PMC6800815` vs `PMC9798250`）在水性 NanoDSC 下峰温可达 **92–115 °C**（且载铁 ~700 Fe 会使 L 富集型降低 12–14 °C、H 富集型降低 18–20 °C）。

---

## 六、 存在但暂不可获取 / 待补齐核对来源统计（共 185 项）

> 完整逐条清单见 [inaccessible_and_follow_up_sources.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/inaccessible_and_follow_up_sources.csv)。

| 待补齐类别 / 优先级 (`priority`) | 条目数 | 核心内容示例 |
| :--- | :---: | :--- |
| `P1_High_Closed_ThermalMutant` | 29 | 见 CSV 对应明细及获取建议 |
| `P1_High_DB_Portal` | 1 | 见 CSV 对应明细及获取建议 |
| `P1_High_Quant_Curve_In_ESI` | 1 | 见 CSV 对应明细及获取建议 |
| `P1_High_Quant_FigureOnly` | 1 | 见 CSV 对应明细及获取建议 |
| `P1_High_Unpublished_Thermostab_PDB` | 14 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_DB_Portal` | 2 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_ESI` | 2 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_FigureOnly` | 2 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_PDF_Interleaved` | 1 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_Patents` | 1 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_Preprints` | 1 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_SI_Matrix` | 1 | 见 CSV 对应明细及获取建议 |
| `P2_Medium_ZH_Literature` | 1 | 见 CSV 对应明细及获取建议 |
| `P3_Low_Related_OtherTitle` | 16 | 见 CSV 对应明细及获取建议 |
| `P3_Medium_Closed_FerritinPaper` | 80 | 见 CSV 对应明细及获取建议 |
| `P3_Unpublished_Mutant_PDB` | 24 | 见 CSV 对应明细及获取建议 |
| `P4_OffTopic_Or_NoMeta_ClosedDOI` | 8 | 见 CSV 对应明细及获取建议 |

### 最值得优先补齐的高价值缺口
1. **已有定量摘要/明确突变热稳定性标题的闭源付费论文（29 篇高优先级）**：例如牡蛎铁蛋白 `D120G`（`89.17 °C`, `+13 °C`, `10.1016/j.foodchem.2020.128879`）、虾铁蛋白 `T158R/G159R` 与 `D132R`（`10.1016/j.foodchem.2021.129089`、`10.1016/j.ijbiomac.2021.09.067`）、牛蛙 M 铁蛋白 C4 孔道渐变突变系列（`~ -5 °C`, `10.1021/acs.biochem.5c00744`）、C3 界面重设计 `+20 °C` 笼蛋白（`10.1039/c3cc40886h`）、AF3 设计的无半胱氨酸高热稳 `rHuHF 9LUW`（`10.1021/acs.jafc.5c02232`）、`Dps` 疏水内腔系列 `Tm = 74–90 °C`（`10.1021/ja057069x`）、人 `Ecut-HFn`（`10.1016/j.ijbiomac.2023.126973`）、`EcBfr D56R/E60H` 全文（`10.1093/jb/mvv065`）及 `马脾/HuHF N23/芳环簇/E-helix` 系列（`10.1002/chem.202300488`, `10.1021/bi200207w`, `10.1021/bi301189t`）。
2. **OA 全文的补充材料（ESI）与仅作图未列表数据**：`PMC11776510`（`MtBfrA E121X` RSC ESI Fig S13 精确小数）、`PMC9798250`（Fig 3b 中间 H:L 比例曲线值）、`PMC8816334`（Sci Adv 24-variant 完整 Table S1–S3）、`PMC11452399`（Fig S4g `FTH1-1` DSF 数据）。
3. **RCSB PDB `38` 个未发表（`To Be Published` / 无 PMID）突变结构**：特别是标题直指热稳定的 `9BXP`（`PaFtnA`）、`9BYE`（`EcFtnA`）、`9LUW`（`HuHF`）、`8YKR`（计算辅助二硫键胃稳定性优化）、`7XRG` 等。
4. **中文与专利库深度补全**：知网/万方《食品科学》《生物工程学报》硕博论文全文，以及 `ProThermDB` ~31k 全量包、`SoluProtMutDB`、`BRENDA` 门户检索。

---

## 七、 交付文件清单
- [database_and_literature_source_inventory.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/database_and_literature_source_inventory.csv)
- [ferritin_mutant_thermostability_dataset.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/ferritin_mutant_thermostability_dataset.csv)
- [unique_mutants_summary.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/unique_mutants_summary.csv)
- [inaccessible_and_follow_up_sources.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/inaccessible_and_follow_up_sources.csv)
- [README_database_summary_report.md](file:///Users/a26683/Desktop/ferritinMTdatabase/README_database_summary_report.md)
