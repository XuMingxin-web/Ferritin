#!/usr/bin/env python3
"""Generate publication-quality figures for the AIDD8122 Final Report and Slides.

Reads directly from curated CSVs and JSON summaries under ferritinMTdatabase.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = REPO_ROOT / "docs" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_DIR = REPO_ROOT / "research" / "results"
DATA_DIR = REPO_ROOT / "research" / "data"

# Global publication styling
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#455A64",
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

NAVY = "#1B365D"
TEAL = "#00796B"
TEAL_LIGHT = "#4DB6AC"
ORANGE = "#D95F02"
CRIMSON = "#C62828"
SLATE = "#546E7A"
GRAY = "#90A4AE"
GOLD = "#E6A100"
PURPLE = "#6A3D9A"


def _save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.06)
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def compute_roc(y_true: np.ndarray, y_score: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute exact empirical ROC curve and trapezoidal AUC."""
    order = np.argsort(-y_score, kind="mergesort")
    y_sorted = y_true[order]
    n_pos = np.sum(y_true == 1)
    n_neg = np.sum(y_true == 0)
    tpr = [0.0]
    fpr = [0.0]
    tp = 0
    fp = 0
    scores_sorted = y_score[order]
    i = 0
    while i < len(scores_sorted):
        j = i
        while j < len(scores_sorted) and scores_sorted[j] == scores_sorted[i]:
            if y_sorted[j] == 1:
                tp += 1
            else:
                fp += 1
            j += 1
        tpr.append(tp / n_pos)
        fpr.append(fp / n_neg)
        i = j
    auc = float(np.trapezoid(tpr, fpr)) if hasattr(np, "trapezoid") else float(np.trapz(tpr, fpr))
    return np.array(fpr), np.array(tpr), auc


def make_fig1_database_audit() -> None:
    ds = pd.read_csv(REPO_ROOT / "ferritin_mutant_thermostability_dataset.csv")
    muts = ds[ds["record_role"] != "WT_baseline"].copy()
    cohorts = pd.read_csv(DATA_DIR / "curated_ferritin_cohorts.csv")

    fam_order = [
        ("Family_1_Canonical_24mer_Ferritin", "Fam 1:\nCanonical\n24-mer"),
        ("Family_2_Bacterioferritin_Bfr", "Fam 2:\nBacterio-\nferritin"),
        ("Family_3_Miniferritin_Dps", "Fam 3:\nDps\n12-mer"),
        ("Family_4_Engineered_Fusion_Display_Ferritin", "Fam 4:\nFusion /\nDisplay"),
        ("Family_5_Extended_Superfamily_And_Nanocages", "Fam 5:\nExtended\nCages"),
    ]
    labels = [lbl for _, lbl in fam_order]
    keys = [k for k, _ in fam_order]

    wet_rec = [len(muts[(muts["family_tier"] == k) & (muts["evidence_tier"] != "D_Computational_Prediction")]) for k in keys]
    comp_rec = [len(muts[(muts["family_tier"] == k) & (muts["evidence_tier"] == "D_Computational_Prediction")]) for k in keys]

    uniq_df = muts.sort_values("evidence_tier").drop_duplicates("unique_mutant_id")
    wet_uniq = [len(uniq_df[(uniq_df["family_tier"] == k) & (uniq_df["evidence_tier"] != "D_Computational_Prediction")]) for k in keys]
    comp_uniq = [len(uniq_df[(uniq_df["family_tier"] == k) & (uniq_df["evidence_tier"] == "D_Computational_Prediction")]) for k in keys]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.1), gridspec_kw={"width_ratios": [1.18, 1.0]})

    x = np.arange(len(keys))
    w = 0.36
    b1 = ax1.bar(x - w / 2, wet_rec, w, color=NAVY, label="Records: Wet-lab (A–C)")
    b2 = ax1.bar(x - w / 2, comp_rec, w, bottom=wet_rec, color=GRAY, label="Records: Comp. Only (D)")
    b3 = ax1.bar(x + w / 2, wet_uniq, w, color=TEAL, label="Unique Mutants: Wet-lab (A–C)")
    b4 = ax1.bar(x + w / 2, comp_uniq, w, bottom=wet_uniq, color=TEAL_LIGHT, alpha=0.65, label="Unique Mutants: Comp. Only (D)")

    for i in range(len(keys)):
        tot_r = wet_rec[i] + comp_rec[i]
        tot_u = wet_uniq[i] + comp_uniq[i]
        ax1.text(x[i] - w / 2, tot_r + 2, str(tot_r), ha="center", va="bottom", fontsize=8, fontweight="bold", color=NAVY)
        ax1.text(x[i] + w / 2, tot_u + 2, str(tot_u), ha="center", va="bottom", fontsize=8, fontweight="bold", color=TEAL)

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=8.2)
    ax1.set_ylabel("Count (Mutant Rows / Unique Mutants)")
    ax1.set_ylim(0, 125)
    ax1.grid(axis="y", linestyle="--", alpha=0.35)
    ax1.set_title("A. Dual-Count Curated Database (292 Rows / 240 Unique)")
    ax1.legend(frameon=False, loc="upper right", fontsize=7.6)

    # Panel B: Curated benchmark N=75 breakdown by subcohort & oligomer symmetry
    cohort_defs = [
        ("S1_PairedCore23", "S1: Paired Core 24-mer (n=23)"),
        ("S1_SingleOther", "S1: Oyster/Patent/MtBfrA (n=7)"),
        ("S2_EcDps_Single8", "S2: EcDps 12-mer Single (n=8)"),
        ("S3_Exact9", "S3: Exact Multi-Point (n=9)"),
        ("S3_Abstract7", "S3: Abstract Multi-Point (n=7)"),
        ("S4_Secondary", "S4: DSC Dual-Peak / Ti / HL (n=8)"),
        ("S5_TierBC_Directional", "S5: Unseen Directional (n=13)"),
    ]
    c_counts_24 = []
    c_counts_12 = []
    c_labels = []
    for key, lbl in cohort_defs:
        if key == "S1_SingleOther":
            sub = cohorts[(cohorts["cohort"] == "S1_SinglePoint_24mer") & (cohorts["subcohort"] != "S1_PairedCore23")]
        elif key in ("S4_Secondary", "S5_TierBC_Directional"):
            sub = cohorts[cohorts["cohort"] == key]
        else:
            sub = cohorts[cohorts["subcohort"] == key]
        n12 = int((sub["oligomer_family"] == "12mer_dps").sum())
        n24 = int(len(sub) - n12)
        c_counts_24.append(n24)
        c_counts_12.append(n12)
        c_labels.append(lbl)

    y = np.arange(len(c_labels))[::-1]
    ax2.barh(y, c_counts_24, height=0.62, color=NAVY, label="24-mer Octahedral Cage (64)")
    ax2.barh(y, c_counts_12, height=0.62, left=c_counts_24, color=ORANGE, label="12-mer Tetrahedral Dps (9)")
    # Distinguish S4 heteropolymer/24mer count accurately: total S2=8 + S3=1 = 9 are 12mer; 66 are 24mer! Fix label count dynamically:
    ax2.clear()
    n24_tot = sum(c_counts_24)
    n12_tot = sum(c_counts_12)
    ax2.barh(y, c_counts_24, height=0.62, color=NAVY, label=f"24-mer Octahedral Cage (n={n24_tot})")
    ax2.barh(y, c_counts_12, height=0.62, left=c_counts_24, color=ORANGE, label=f"12-mer Tetrahedral Dps (n={n12_tot})")
    for idx, yi in enumerate(y):
        tot = c_counts_24[idx] + c_counts_12[idx]
        role = "Train/CV" if idx <= 1 else ("Zero-Shot Eval" if idx in (2, 3, 6) else "Diagnostic")
        ax2.text(tot + 0.5, yi, f"{tot} ({role})", va="center", ha="left", fontsize=7.8, color="#263238")

    ax2.set_yticks(y)
    ax2.set_yticklabels(c_labels, fontsize=8.2)
    ax2.set_xlabel("Curated Benchmark Cases (Total N = 75 across 14 WT Assemblies)")
    ax2.set_xlim(0, 35)
    ax2.set_ylim(-1.4, 6.65)
    ax2.grid(axis="x", linestyle="--", alpha=0.35)
    ax2.set_title("B. Benchmark Cohorts S1–S5 & Symmetry Split")
    ax2.legend(frameon=False, loc="lower right", ncol=1, fontsize=8.0)

    fig.tight_layout(w_pad=2.2)
    _save(fig, "fig1_database_audit")


def make_fig2_s1_baselines_and_cv() -> None:
    base = json.load(open(RESULTS_DIR / "baseline_evaluation_summary.json"))["baselines"]["S1_All30"]["models"]
    hyb = json.load(open(RESULTS_DIR / "hybrid_benchmark_summary.json"))
    cv = hyb["cv_model_family_comparison"]["S1_All30"]
    pm = hyb["primary_model"]

    # Panel A: Pearson r & MAE comparison on S1_All30
    entries = [
        ("ESM-2 WT LLR", base["ESM2_WT_marginal_LLR"]["pearson_r"], base["ESM2_WT_marginal_LLR"]["mae_c"], CRIMSON),
        ("FoldX Monomer", base["FoldX_mono_neg_ddG"]["pearson_r"], base["FoldX_mono_neg_ddG"]["mae_c"], SLATE),
        ("Mono GeoDTm-3D", base["Monomer_GeoDTm3D_k32"]["pearson_r"], base["Monomer_GeoDTm3D_k32"]["mae_c"], SLATE),
        ("Asm GeoDTm-3D", base["Assembly_GeoDTm3D_k32"]["pearson_r"], base["Assembly_GeoDTm3D_k32"]["mae_c"], SLATE),
        ("FoldX Cage Eff.", base["FoldX_cage_neg_eff_ddG"]["pearson_r"], base["FoldX_cage_neg_eff_ddG"]["mae_c"], TEAL_LIGHT),
        ("FoldX AC Interface", base["FoldX_AC_neg_bind_ddG"]["pearson_r"], base["FoldX_AC_neg_bind_ddG"]["mae_c"], TEAL_LIGHT),
        ("S4346 Frozen Ridge", base["S4346_Frozen_Ridge_Sub"]["pearson_r"], base["S4346_Frozen_Ridge_Sub"]["mae_c"], GOLD),
        ("ShallowGBDT (LOPO)", cv["ShallowGBDT_5feat"]["LOPO"]["pearson_r"], cv["ShallowGBDT_5feat"]["LOPO"]["mae_c"], GRAY),
        ("RidgeCV (LOPO)", cv["RidgeCV_5feat"]["LOPO"]["pearson_r"], cv["RidgeCV_5feat"]["LOPO"]["mae_c"], GRAY),
        ("BayesRidge (LOPO)", cv["BayesianRidge_5feat"]["LOPO"]["pearson_r"], cv["BayesianRidge_5feat"]["LOPO"]["mae_c"], GRAY),
        ("ZeroPrior (LOPO)", cv["ZeroPrior_NonNegHuber"]["LOPO"]["pearson_r"], cv["ZeroPrior_NonNegHuber"]["LOPO"]["mae_c"], TEAL),
        ("Residual MAP (LOPO)", cv["Residual_NonNegMAP"]["LOPO"]["pearson_r"], cv["Residual_NonNegMAP"]["LOPO"]["mae_c"], TEAL),
        ("GeoHybrid (LOPO)", cv["FerritinGeoHybrid_NonNegMAP"]["LOPO"]["pearson_r"], cv["FerritinGeoHybrid_NonNegMAP"]["LOPO"]["mae_c"], NAVY),
        ("GeoHybrid (LOFO)", cv["FerritinGeoHybrid_NonNegMAP"]["LOFO"]["pearson_r"], cv["FerritinGeoHybrid_NonNegMAP"]["LOFO"]["mae_c"], ORANGE),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.5), gridspec_kw={"width_ratios": [1.35, 1.0]})

    y = np.arange(len(entries))[::-1]
    r_vals = [e[1] for e in entries]
    maes = [e[2] for e in entries]
    cols = [e[3] for e in entries]
    names = [e[0] for e in entries]

    ax1.barh(y, r_vals, height=0.66, color=cols, alpha=0.92)
    ax1.axvline(0.0, color="#263238", linewidth=0.9)
    for idx, yi in enumerate(y):
        rv = r_vals[idx]
        mv = maes[idx]
        if rv >= 0:
            ax1.text(rv + 0.012, yi, f"{rv:+.3f} (MAE {mv:.2f}°C)", va="center", ha="left", fontsize=7.6, color="#212121",
                     fontweight="bold" if "GeoHybrid" in names[idx] else "normal")
        else:
            ax1.text(rv - 0.012, yi, f"{rv:+.3f} (MAE {mv:.2f}°C)", va="center", ha="right", fontsize=7.6, color="#212121")

    ax1.set_yticks(y)
    ax1.set_yticklabels(names, fontsize=8.1)
    ax1.set_xlim(-0.46, 0.72)
    ax1.set_xlabel("Pearson Correlation r on S1_All30 (n = 30)")
    ax1.grid(axis="x", linestyle="--", alpha=0.35)
    ax1.set_title("A. S1_All30 Pearson r & MAE: Baselines vs. OOF Hybrid")

    # Panel B: Standardized coefficients stability across Deployed / LOPO / LOFO
    feat_map = [
        ("neg_fx_ac_ddg_bind_c3", "C3 Interface\nBind (-ΔΔG_C3)"),
        ("neg_fx_ac_ddg_bind_nonpolar", "Non-Polar AC\nBind (-ΔΔG_np)"),
        ("neg_geo_asm_enc_shift_k32", "ESM Context\nShift (-ΔGeo_enc)"),
        ("neg_fx_cage_ddg_eff", "Cage Effective\n(-ΔΔG_cage)"),
        ("geo_asm_k32", "Assembly\nGeoDTm-3D"),
    ]
    dep = pm["deployed_weights_s1all30"]["standardized_coefs"]
    lopo_st = pm["lopo_coefficient_stability_s1all30"]
    lofo_st = pm["lofo_coefficient_stability_s1all30"]

    yy = np.arange(len(feat_map))[::-1]
    h = 0.25
    dep_w = [dep[k] for k, _ in feat_map]
    lopo_m = [lopo_st[k]["mean"] for k, _ in feat_map]
    lopo_lo = [lopo_st[k]["mean"] - lopo_st[k]["min"] for k, _ in feat_map]
    lopo_hi = [lopo_st[k]["max"] - lopo_st[k]["mean"] for k, _ in feat_map]
    lofo_m = [lofo_st[k]["mean"] for k, _ in feat_map]
    lofo_lo = [lofo_st[k]["mean"] - lofo_st[k]["min"] for k, _ in feat_map]
    lofo_hi = [lofo_st[k]["max"] - lofo_st[k]["mean"] for k, _ in feat_map]

    ax2.barh(yy + h, dep_w, height=h, color=NAVY, label="Deployed (S1_All30)")
    ax2.barh(yy, lopo_m, height=h, color=TEAL, xerr=[lopo_lo, lopo_hi], capsize=2.5, error_kw={"lw": 1.0}, label="LOPO (mean ± min/max)")
    ax2.barh(yy - h, lofo_m, height=h, color=ORANGE, xerr=[lofo_lo, lofo_hi], capsize=2.5, error_kw={"lw": 1.0}, label="LOFO <35% id (mean ± min/max)")
    ax2.axvline(1.0, color=SLATE, linestyle=":", linewidth=1.1, label="Prior anchor w_prior = 1.0")

    for idx, yi in enumerate(yy):
        ax2.text(max(dep_w[idx], lopo_m[idx] + lopo_hi[idx], lofo_m[idx] + lofo_hi[idx]) + 0.04, yi + 0.05,
                 f"{dep_w[idx]:.2f}", va="center", ha="left", fontsize=7.8, fontweight="bold", color=NAVY)

    ax2.set_yticks(yy)
    ax2.set_yticklabels([lbl for _, lbl in feat_map], fontsize=8.2)
    ax2.set_xlim(0.0, 1.85)
    ax2.set_ylim(-1.75, 4.6)
    ax2.set_xlabel("Standardized Regression Weight (All 100% Positive)")
    ax2.grid(axis="x", linestyle="--", alpha=0.35)
    ax2.set_title("B. Prior-Anchored Weight Stability Across Folds")
    ax2.legend(frameon=False, loc="lower center", ncol=2, fontsize=7.5)

    fig.tight_layout(w_pad=2.0)
    _save(fig, "fig2_s1_baselines_and_cv")


def make_fig3_ablation_and_roc() -> None:
    hyb = json.load(open(RESULTS_DIR / "hybrid_benchmark_summary.json"))["ablation_study"]
    preds = pd.read_csv(RESULTS_DIR / "hybrid_predictions_all.csv")

    ab_keys = [
        ("Ab1_GeoMono3D_Calibrated", "Ab1\nMono3D"),
        ("Ab2_GeoAsm3D_Plus_EncShift", "Ab2\n+Asm/Enc"),
        ("Ab3_Plus_CageEffectiveStability", "Ab3\n+CageEff"),
        ("Ab4_Plus_C3Interface", "Ab4\n+C3 Iface"),
        ("Ab5_Full_Ferritin_GeoHybrid_5feat", "Ab5 (Full)\n+NonPolar"),
        ("Ab6_Extended_Plus_ESM2_LLR", "Ab6\n+ESM2 LLR"),
    ]
    x = np.arange(len(ab_keys))
    lofo_r = [hyb[k]["S1_All30_LOFO"]["pearson_r"] for k, _ in ab_keys]
    lofo_mae = [hyb[k]["S1_All30_LOFO"]["mae_c"] for k, _ in ab_keys]
    s5_auc = [hyb[k]["ZeroShot_S5_Directional13"]["roc_auc"] for k, _ in ab_keys]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax1_r = ax1
    l1 = ax1_r.plot(x, lofo_r, marker="o", lw=2.2, ms=6.5, color=NAVY, label="S1 LOFO Pearson r (left)")
    ax1_r.axhline(0, color=GRAY, lw=0.8, ls="--")
    for xi, rv in zip(x, lofo_r):
        if xi == 2:  # Ab3: place below blue marker so orange 9.90 above does not overlap
            ax1_r.text(xi + 0.08, rv - 0.038, f"{rv:+.3f}", ha="center", va="top", fontsize=7.8, color=NAVY)
        else:
            dy = 0.032 if xi != 4 else 0.036
            ax1_r.text(xi, rv + dy, f"{rv:+.3f}", ha="center", va="bottom", fontsize=7.8, color=NAVY,
                       fontweight="bold" if xi == 4 else "normal")
    ax1_r.set_ylabel("S1 LOFO Pearson r", color=NAVY)
    ax1_r.set_ylim(-0.32, 0.56)
    ax1_r.set_xticks(x)
    ax1_r.set_xticklabels([lbl for _, lbl in ab_keys], fontsize=8.0)
    ax1_r.axvspan(3.65, 4.35, color="#E8F5E9", alpha=0.8, zorder=0)
    ax1_r.grid(axis="x", linestyle=":", alpha=0.4)

    ax1_m = ax1_r.twinx()
    ax1_m.spines["top"].set_visible(False)
    l2 = ax1_m.plot(x, lofo_mae, marker="s", lw=1.9, ls="--", ms=5.5, color=ORANGE, label="S1 LOFO MAE °C (right)")
    for xi, mv in zip(x, lofo_mae):
        ax1_m.text(xi, mv + (0.15 if xi < 4 else -0.28), f"{mv:.2f}", ha="center",
                   va="bottom" if xi < 4 else "top", fontsize=7.6, color=ORANGE)
    ax1_m.set_ylabel("S1 LOFO MAE (°C)", color=ORANGE)
    ax1_m.set_ylim(8.2, 11.7)

    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1_r.legend(lines, labels, frameon=False, loc="upper left", fontsize=7.8)
    ax1_r.set_title("A. Stepwise Mechanistic Ablation (S1 LOFO <35% id)")

    # Panel B: S5 ROC Curves (N=13)
    s5 = preds[preds["cohort"] == "S5_TierBC_Directional"].copy()
    y_true = s5["binary_direction"].to_numpy(int)

    curves = [
        ("FoldX Cage Eff. (-ΔΔG_cage)", s5["neg_fx_cage_ddg_eff"].to_numpy(float), TEAL, "-", 3.6, 0.45),
        ("Ferritin-GeoHybrid", s5["hybrid_deployed_s1all30"].to_numpy(float), NAVY, "-", 2.1, 1.0),
        ("FoldX AC Bind (-ΔΔG_AC)", s5["neg_fx_ac_ddg_bind_total"].to_numpy(float), GOLD, "-.", 1.7, 0.95),
        ("Mono GeoDTm-3D (k=32)", s5["geo_mono_k32"].to_numpy(float), SLATE, ":", 1.8, 0.95),
        ("ESM-2 WT-Marginal LLR", s5["esm2_wt_marginal_llr"].to_numpy(float), CRIMSON, "-", 1.6, 0.9),
    ]
    for cname, scores, col, ls, lw, alp in curves:
        fpr, tpr, auc = compute_roc(y_true, scores)
        ax2.plot(fpr, tpr, color=col, linestyle=ls, linewidth=lw, alpha=alp, label=f"{cname} (AUC={auc:.3f})")

    ax2.plot([0, 1], [0, 1], color=GRAY, linestyle="--", linewidth=0.9, label="Random (AUC=0.500)")
    ax2.set_xlabel("False Positive Rate (1 - Specificity)")
    ax2.set_ylabel("True Positive Rate (Sensitivity)")
    ax2.set_xlim(-0.02, 1.02)
    ax2.set_ylim(-0.02, 1.04)
    ax2.grid(True, linestyle=":", alpha=0.35)
    ax2.set_title("B. Zero-Shot ROC on Unseen S5 Directional (n = 13)")
    ax2.legend(frameon=False, loc="lower right", fontsize=7.6)

    fig.tight_layout(w_pad=2.2)
    _save(fig, "fig3_ablation_and_roc")


def make_fig4_zeroshot_transfer() -> None:
    preds = pd.read_csv(RESULTS_DIR / "hybrid_predictions_all.csv")
    hyb = json.load(open(RESULTS_DIR / "hybrid_benchmark_summary.json"))["external_zero_shot_transfer"]

    s1 = preds[preds["cohort"] == "S1_SinglePoint_24mer"].copy()
    s1["unseen_pred"] = s1["hybrid_oof_lopo"]
    s1["group_lbl"] = "S1 LOPO (24-mer, n=30)"

    s2 = preds[preds["subcohort"] == "S2_EcDps_Single8"].copy()
    s2["unseen_pred"] = s2["hybrid_deployed_s1all30"]
    s2["group_lbl"] = "S2 EcDps (12-mer, n=8)"

    s3 = preds[preds["subcohort"] == "S3_Exact9"].copy()
    s3["unseen_pred"] = s3["hybrid_deployed_s1all30"]
    s3["group_lbl"] = "S3 Multi-Point Exact (n=9)"

    pooled = pd.concat([s1, s2, s3], ignore_index=True)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.1), gridspec_kw={"width_ratios": [1.15, 1.05, 1.05]})
    ax1, ax2, ax3 = axes

    # Panel A: Pooled N=47 Scatter
    styles = [
        ("S1 LOPO (24-mer, n=30)", NAVY, "o"),
        ("S2 EcDps (12-mer, n=8)", ORANGE, "^"),
        ("S3 Multi-Point Exact (n=9)", TEAL, "s"),
    ]
    for glbl, col, mk in styles:
        sub = pooled[pooled["group_lbl"] == glbl]
        ax1.scatter(sub["delta_tm"], sub["unseen_pred"], c=col, marker=mk, s=42, alpha=0.88, edgecolors="white", linewidths=0.5, label=glbl)

    lims = [-34, 30]
    ax1.plot(lims, lims, color="#263238", linestyle="--", linewidth=1.0, label="Identity y = x")
    ax1.axhline(0, color=GRAY, lw=0.7)
    ax1.axvline(0, color=GRAY, lw=0.7)
    # Highlight Mono32 range band
    ax1.axhspan(-12, 5, color="#ECEFF1", alpha=0.45, zorder=0)
    ax1.text(-32, 24.5, "Unseen Pooled N = 47\nr = +0.443 (p = 0.0018)\nρ = +0.433 (p = 0.0023)\nMAE = 9.08 °C | Sign = 75.6%",
             fontsize=7.6, va="top", bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#CFD8DC", alpha=0.92))
    ax1.set_xlim(lims)
    ax1.set_ylim(-20, 28)
    ax1.set_xlabel("Experimental ΔTm (°C)")
    ax1.set_ylabel("Unseen Predicted ΔTm (°C)")
    ax1.set_title("A. Pooled Unseen Parity (N = 47)")
    ax1.legend(frameon=False, loc="lower right", fontsize=7.2)

    # Panel B: S2 EcDps 12-mer Case-by-case Error & Comparison
    s2_sorted = s2.sort_values("delta_tm")
    yy = np.arange(len(s2_sorted))
    h = 0.27
    ax2.barh(yy + h, s2_sorted["delta_tm"], height=h, color="#37474F", label="Experimental ΔTm")
    ax2.barh(yy, s2_sorted["hybrid_deployed_s1all30"], height=h, color=TEAL, label="GeoHybrid (MAE 5.09°C)")
    ax2.barh(yy - h, s2_sorted["geo_mono_k32"], height=h, color=GRAY, label="Mono GeoDTm-3D (MAE 7.93°C)")
    ax2.axvline(0, color="#263238", lw=0.8)
    ax2.set_yticks(yy)
    ax2.set_yticklabels(s2_sorted["mutation_canon"], fontsize=8.0)
    ax2.set_xlim(-22, 2.5)
    ax2.set_ylim(-2.4, 7.6)
    ax2.set_xlabel("ΔTm (°C) — 100% Sign Accuracy, 7/8 Wins")
    ax2.grid(axis="x", linestyle="--", alpha=0.35)
    ax2.set_title("B. Zero-Shot S2: 12-mer EcDps (n = 8)")
    ax2.legend(frameon=False, loc="lower left", fontsize=7.3)

    # Panel C: Paired Wilcoxon MAE Reduction on Pooled N=47
    p47_pri = hyb["Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47"]["deployed_S1All30"]
    p47_res = hyb["Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47"]["residual_companion"]
    refs = [
        ("paired_vs_zero", "Zero (ΔTm=0)"),
        ("paired_vs_fx_mono", "FoldX Mono"),
        ("paired_vs_mono_k32", "Mono GeoDTm-3D"),
        ("paired_vs_fx_cage", "FoldX Cage Eff."),
        ("paired_vs_fx_ac", "FoldX AC Bind"),
        ("paired_vs_asm_k32", "Asm GeoDTm-3D"),
    ]
    y3 = np.arange(len(refs))[::-1]
    h3 = 0.34
    pri_diff = [p47_pri[rk]["mean_abs_err_diff_c"] for rk, _ in refs]
    res_diff = [p47_res[rk]["mean_abs_err_diff_c"] for rk, _ in refs]
    pri_p = [p47_pri[rk]["wilcoxon_p"] for rk, _ in refs]
    pri_w = [p47_pri[rk]["wins"] for rk, _ in refs]

    ax3.barh(y3 + h3 / 2, pri_diff, height=h3, color=NAVY, label="Ferritin-GeoHybrid (MAE 9.08°C)")
    ax3.barh(y3 - h3 / 2, res_diff, height=h3, color=ORANGE, label="Residual NonNegMAP (MAE 8.78°C)")
    ax3.axvline(0, color="#263238", lw=0.8)
    for idx, yi in enumerate(y3):
        p_txt = f"p={pri_p[idx]:.4f}" if pri_p[idx] < 0.01 else f"p={pri_p[idx]:.3f}"
        ax3.text(min(pri_diff[idx], res_diff[idx]) - 0.08, yi, f"{pri_w[idx]}/47 wins, {p_txt}",
                 va="center", ha="right", fontsize=7.3, color="#263238")

    ax3.set_yticks(y3)
    ax3.set_yticklabels([lbl for _, lbl in refs], fontsize=8.0)
    ax3.set_xlim(-5.1, 0.3)
    ax3.set_ylim(-1.75, 5.65)
    ax3.set_xlabel("Paired ΔMAE vs. Baseline (°C, negative = better)")
    ax3.grid(axis="x", linestyle="--", alpha=0.35)
    ax3.set_title("C. Pooled N = 47 Wilcoxon Error Gain")
    ax3.legend(frameon=False, loc="lower left", fontsize=7.3)

    fig.tight_layout(w_pad=1.8)
    _save(fig, "fig4_zeroshot_transfer")


def make_fig5_error_audit() -> None:
    hyb = json.load(open(RESULTS_DIR / "hybrid_benchmark_summary.json"))
    errs = hyb["s1_lopo_largest_errors"][:8]
    short_names = [
        "EcBfr D118W\n(core cavity)",
        "EcBfr D118F\n(core cavity)",
        "EcBfr D118Y\n(core cavity)",
        "MjFer D132R\n(C3 salt bridge)",
        "CgFer D120G\n(3-fold channel)",
        "EcBfr M52A\n(heme pocket)",
        "EcBfr R30A\n(C2+heme)",
        "EcBfr R61A\n(buried/C2)",
    ]
    obs = [e["delta_tm"] for e in errs]
    mono = [e["mono_geo32"] for e in errs]
    oof = [e["oof_lopo"] for e in errs]

    fig, ax = plt.subplots(figsize=(9.6, 3.8))
    x = np.arange(len(errs))
    w = 0.26
    ax.bar(x - w, obs, w, color="#263238", label="Experimental ΔTm (°C)")
    ax.bar(x, mono, w, color=GRAY, label="Monomer GeoDTm-3D k=32")
    ax.bar(x + w, oof, w, color=ORANGE, label="Ferritin-GeoHybrid (LOPO OOF)")
    ax.axhline(0, color="#263238", lw=0.9)

    for xi, (ov, mv, hv) in enumerate(zip(obs, mono, oof)):
        ax.text(xi - w, ov + (0.9 if ov >= 0 else -1.1), f"{ov:+.1f}", ha="center",
                va="bottom" if ov >= 0 else "top", fontsize=7.5, fontweight="bold", color="#263238")
        ax.text(xi + w, hv + (0.8 if hv >= 0 else -1.0), f"{hv:+.1f}", ha="center",
                va="bottom" if hv >= 0 else "top", fontsize=7.4, color=ORANGE)

    ax.set_xticks(x)
    ax.set_xticklabels(short_names, fontsize=8.0)
    ax.set_ylabel("ΔTm (°C)")
    ax.set_ylim(-24, 32)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_title("Diagnostic Audit of Top S1 LOPO Residuals: Amplitude Compression vs. Extreme Non-Linear Transitions")
    ax.legend(frameon=False, loc="upper right", ncol=3, fontsize=8.2)

    fig.tight_layout()
    _save(fig, "fig5_error_audit")


def main() -> None:
    make_fig1_database_audit()
    make_fig2_s1_baselines_and_cv()
    make_fig3_ablation_and_roc()
    make_fig4_zeroshot_transfer()
    make_fig5_error_audit()
    print("Generated 5 figures cleanly in:", FIG_DIR)


if __name__ == "__main__":
    main()
