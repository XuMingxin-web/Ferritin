# Ferritin-GeoHybrid: Curated Ferritin Thermostability Database & Assembly-Aware ΔTm Predictor

**Author:** Mingxin Xu (Student ID: `P2625571`)  
**Course:** AIDD8122 – Advanced Topics in Artificial Intelligence & Drug Discovery, Faculty of Applied Sciences, Macao Polytechnic University  
**Repository:** https://github.com/XuMingxin-web/Ferritin-AIDD8122

---

## 1. Overview & Motivation

Predicting mutation-induced changes in melting temperature ($\Delta T_m = T_m^{\text{mut}} - T_m^{\text{WT}}$) for oligomeric ferritin nanocages is a critical bottleneck for engineering hyperthermostable iron-storage carriers, food-grade bio-nanomaterials, and subunit vaccine scaffolds. A systematic audit across **17 public biophysical databases and literature corpora** revealed a severe blind spot:

- Major stability repositories—**FireProtDB** (41.1M SQL lines), **ThermoMutDB** (13,337 records), **GeoStab** (`S4346`/`S571`/`S669`/`S8754`, 14,340 rows), **ProteinGym v1.x** (283 DMS assays), **Tsuboyama Mega-scale** (9,728 entries), and **MaveDB** (2,161 score sets)—contain **0 natural ferritin point-mutant $\Delta T_m$ records** (ProThermDB only listed 4 synthetic EBNA1–HuHF fusion rows mislabeled as `wild-type`).
- Existing predictors fail on ferritins for two mechanistic reasons: (i) **monomer-only inference** ignores $C_2$, $C_3$, and $C_4$ oligomeric assembly interfaces that govern nanocage integrity, and (ii) **evolutionary language models (ESM-2 WT-marginal LLR)** invert correlation on functional channel/cavity hotspots ($r = -0.363$ on `S1_PairedCore23`) where stability-enhancing substitutions violate native evolutionary conservation.

To solve both gaps, this repository provides:
1. **A curated, dual-reconciled Ferritin Thermostability Database** (**334 rows**: 292 mutant condition records + 42 WT baselines across **240 unique mutants**, **0 validation errors**).
2. **`Ferritin-GeoHybrid`**: A **5-feature non-negative prior-anchored Huber MAP linear model** (`NonNegHuberMAPRegressor`, $\delta=8.0^\circ\text{C}$, $\alpha=80.0$, $w_{\text{prior}}=1.0$) integrating BIOMT biological assembly GeoDTm-3D inference, FoldX 5.1 cage/interface energy decomposition ($C_3$ and non-polar $A\leftrightarrow C$ binding terms), and evolutionary-context shift (`geo_asm_k32 - geo_asm_noesm_k32`).

---

## 2. Key Benchmark Results

| Evaluation Regime | Metric / Comparison | Baseline Best / Reference | `Ferritin-GeoHybrid` |
|---|---|---|---|
| **S1_All30 Leave-One-Protein-Out (LOPO, $N=30$)** | MAE / Pearson $r$ | Mono GeoDTm-3D $k=32$: MAE $9.43^\circ\text{C}$, $r=-0.010$ | **MAE $8.58^\circ\text{C}$, $r=+0.375$ ($p=0.0409$)** |
| **S1_All30 Leave-One-Family-Out (LOFO $<35\%$ id, $N=30$)** | MAE / Pearson $r$ | Mono GeoDTm-3D $k=32$: $r_{\text{Ab1}}=-0.187$ | **MAE $9.15^\circ\text{C}$, $r=+0.394$ ($p=0.0311$)** |
| **Zero-Shot Transfer: S2 EcDps 12-mer ($N=8$)** | MAE / Sign Acc / $p$ | Mono GeoDTm-3D $k=32$: MAE $7.93^\circ\text{C}$ | **MAE $5.09^\circ\text{C}$, 100% sign ($p=0.0234$, 7/8 wins)** |
| **Zero-Shot Transfer: S3 Exact Multi-Point ($N=9$)** | Spearman $\rho$ / Pearson $r$ | Zero baseline: MAE $18.64^\circ\text{C}$ | **$\rho=+0.833$ ($p=0.0053$), $r=+0.678$ ($p=0.0447$)** |
| **Zero-Shot Transfer: S5 Unseen Directional ($N=13$)** | ROC-AUC (MW $p$) | Mono/Asm GeoDTm-3D: $0.524$; ESM-2: $0.381$ | **ROC-AUC $0.952$ ($p=0.0023$)** |
| **Pooled All-Unseen (S1 LOPO + S2 + S3 Exact, $N=47$)** | MAE / Pearson $r$ | Mono32 $10.22^\circ\text{C}$; FoldX Mono $11.19^\circ\text{C}$ | **MAE $9.08^\circ\text{C}$, $r=+0.443$ ($p=0.0018$)** |

*(Companion `Residual_NonNegMAP` achieves pooled MAE $= 8.78^\circ\text{C}$, $r=+0.459$, $p=0.0012$, Wilcoxon $p=0.0170$ vs. Monomer GeoDTm-3D $k=32$.)*

---

## 3. Repository Structure

```text
ferritinMTdatabase/
├── README.md                                  # Project overview & reproduction guide
├── README_database_summary_report.md          # Detailed curation & database audit report (ZH)
├── requirements.txt                           # Python dependencies
├── ferritin_mutant_thermostability_dataset.csv# Curated 334-row condition-level dataset
├── unique_mutants_summary.csv                 # Deduplicated 240 unique mutants summary
├── database_and_literature_source_inventory.csv # Audit inventory of 17 databases & corpora
├── inaccessible_and_follow_up_sources.csv     # Paywalled / TBP follow-up catalog (165 items)
├── scripts/
│   ├── build_ferritin_mt_db.py                # Reconstructs CSV deliverables from curated JSON/XML
│   └── validate_deliverables.py               # Independent integrity & dual-count validator
└── research/
    ├── src/
    │   ├── paths.py                           # Portable workspace & env-configurable tool paths
    │   ├── ferritin_curation.py               # Curated N=75 cohort & 14 WT template definitions
    │   ├── assembly_builder.py                # BIOMT 24-mer/12-mer symmetry assembly & 7-chain subcomplex builder
    │   ├── foldx_assembly.py                  # Monomer, subcomplex & C2/C3/C4 binding ddG via FoldX 5.1
    │   ├── geodtm_assembly_adapter.py         # Monomer vs assembly GeoDTm-3D & no-ESM encoder shift
    │   ├── hybrid_model.py                    # NonNegHuberMAPRegressor, Laplace std, feature definitions
    │   └── metrics.py                         # Regression, paired bootstrap, Wilcoxon, ROC-AUC metrics
    ├── scripts/
    │   ├── run_s571_reproduction.py           # Benchmark gate: GeoDTm-3D & FoldX 5.1 S571 verification
    │   ├── extract_all_features.py            # End-to-end feature extraction across S1–S5 (N=75)
    │   └── run_hybrid_benchmark.py            # Trains Ferritin-GeoHybrid, baselines, ablations & OOF transfer
    ├── data/
    │   ├── curated_ferritin_cohorts.csv       # 75 benchmark cases across S1–S5 with structural context
    │   ├── wt_templates.json                  # 14 WT ferritin/Bfr/Dps template metadata & symmetry neighbors
    │   ├── fixed_aa_embedding_20x7.json       # Physicochemical 20x7 AA embeddings
    │   ├── ferritin_pdb_raw/                  # Raw WT PDB assemblies
    │   └── s571_pdb/                          # GeoDTm-3D S571 verification structures
    ├── models/
    │   ├── ferritin_geohybrid.joblib          # Serialized deployed Ferritin-GeoHybrid model (SHA256 verified)
    │   └── ferritin_geohybrid.json            # Deployed standardized/raw weights & hyperparameter metadata
    ├── results/
    │   ├── ferritin_features_all.csv          # Complete precomputed feature matrix (75 cases)
    │   ├── baseline_evaluation_summary.json   # Baseline metrics across S1–S5
    │   ├── hybrid_benchmark_summary.json      # Full CV, ablation (Ab1–Ab6), zero-shot & error breakdowns
    │   ├── hybrid_predictions_all.csv         # Per-case OOF & zero-shot predictions across all 75 cases
    │   └── s571_repro/reproduction_summary.json # GeoDTm-3D & FoldX 5.1 parity gate results
    ├── reports/
    │   └── ferritin_geodtm_evaluation_and_hybrid_model_zh.md # Detailed Chinese technical lab note
    └── tests/
        ├── test_ferritin_pipeline.py          # Symmetry builder, GeoDTm 1CTS parity & feature checks
        └── test_hybrid_model.py               # SHA256 determinism, weight sign, zero-leakage & S2/S5 gates
```

---

## 4. Quick Start & Reproduction

### 4.1 Environment Setup
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 4.2 Database Integrity Validation (Zero External Dependencies)
```bash
python scripts/validate_deliverables.py
```
Expected output: `ALL VALIDATION CHECKS PASSED (0 errors)` confirming 334 rows and 240 unique mutants.

### 4.3 Reproducing Model Training, CV, Ablations & Zero-Shot Transfer
Because precomputed structural, FoldX 5.1, and GeoDTm-3D features for all 75 cases are included in `research/results/ferritin_features_all.csv`, the complete modeling benchmark and unit tests run offline in seconds:
```bash
python research/scripts/run_hybrid_benchmark.py
pytest research/tests/test_hybrid_model.py -v
```

### 4.4 Regenerating Raw Structure Features (Optional External Tools)
To rerun FoldX 5.1 and GeoDTm-3D feature extraction from scratch, configure external paths via environment variables before running `research/scripts/extract_all_features.py`:
```bash
export FERRITIN_EXTERNAL_ROOT=/path/to/external_tools
export FOLDX_BIN=/path/to/foldx5.1
export GEOSTAB_ROOT=/path/to/GeoStab
export ESM2_WEIGHTS=/path/to/esm2_t33_650M_UR50D.pt
python research/scripts/extract_all_features.py
```

---

## 5. License & Citation
Released for academic research and coursework evaluation under AIDD8122 (Macao Polytechnic University). When referencing the curated data or model, please cite the primary experimental sources cataloged in `ferritin_mutant_thermostability_dataset.csv` and this repository.
