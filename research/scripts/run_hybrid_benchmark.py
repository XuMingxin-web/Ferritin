#!/usr/bin/env python3
"""Out-of-fold CV, homology transfer, ablation, and zero-shot validation for Ferritin-GeoHybrid."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import LeaveOneGroupOut

warnings.filterwarnings("ignore", category=ConvergenceWarning)

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = RESEARCH_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from hybrid_model import (  # noqa: E402
    ABLATION_SETS,
    CORE_FEATURES_7,
    HYBRID_FEATURES_5,
    RESIDUAL_HYBRID_FEATURES_4,
    compute_derived_feature_frame,
    fit_ferritin_hybrid,
    load_predictor,
    save_predictor,
)
from metrics import paired_comparison, regression_metrics, roc_auc_binary  # noqa: E402
from paths import MODELS_DIR, RESULTS_DIR  # noqa: E402


# (label, model_family, residual_anchor, feature_names)
MODEL_FAMILIES: list[tuple[str, str, bool, list[str]]] = [
    ("FerritinGeoHybrid_NonNegMAP", "PriorAnchoredNonNegHuber", False, HYBRID_FEATURES_5),
    ("Residual_NonNegMAP", "PriorAnchoredNonNegHuber", True, RESIDUAL_HYBRID_FEATURES_4),
    ("ZeroPrior_NonNegHuber", "ZeroPriorNonNegHuber", False, HYBRID_FEATURES_5),
    ("BayesianRidge_5feat", "BayesianRidge", False, HYBRID_FEATURES_5),
    ("BayesianRidge_FreeSign7", "BayesianRidge", False, CORE_FEATURES_7),
    ("RidgeCV_5feat", "RidgeCV", False, HYBRID_FEATURES_5),
    ("Huber_5feat", "Huber", False, HYBRID_FEATURES_5),
    ("ShallowGBDT_5feat", "ShallowGBDT", False, HYBRID_FEATURES_5),
    ("RegularizedMLP_5feat", "RegularizedMLP", False, HYBRID_FEATURES_5),
]


def _make_site_groups(df: pd.DataFrame) -> np.ndarray:
    site_num = df["mutation_canon"].astype(str).str.extract(r"([0-9]+)")[0].fillna("0")
    return (df["protein_group"].astype(str) + "::" + site_num).to_numpy()


def run_oof_cv(
    df: pd.DataFrame,
    *,
    group_kind: str,
    model_family: str,
    residual_anchor: bool = False,
    feature_names: list[str] | None = None,
    target_col: str = "delta_tm",
    train_cohort_label: str = "S1_All30",
) -> tuple[np.ndarray, np.ndarray | None, list[dict[str, Any]]]:
    """Run strict zero-leakage Leave-One-Group-Out CV."""
    work = df.reset_index(drop=True)
    if group_kind == "SiteGroup":
        groups = _make_site_groups(work)
    elif group_kind == "LOPO":
        groups = work["protein_group"].astype(str).to_numpy()
    elif group_kind == "LOFO":
        groups = work["oligomer_family"].astype(str).to_numpy()
    elif group_kind == "StudyGroup":
        groups = work["study_id"].astype(str).to_numpy()
    else:
        raise ValueError(f"Unknown group_kind={group_kind!r}")

    logo = LeaveOneGroupOut()
    oof_pred = np.full(len(work), np.nan, dtype=np.float64)
    oof_std = np.full(len(work), np.nan, dtype=np.float64)
    has_std = False
    fold_coefs: list[dict[str, Any]] = []

    for tr_idx, te_idx in logo.split(work, groups=groups):
        df_tr = work.iloc[tr_idx].reset_index(drop=True)
        df_te = work.iloc[te_idx].reset_index(drop=True)
        model = fit_ferritin_hybrid(
            df_tr,
            model_family=model_family,
            feature_names=feature_names,
            residual_anchor=residual_anchor,
            target_col=target_col,
            group_col="protein_group",
            train_cohort_label=train_cohort_label,
        )
        mean, std = model.predict_with_std(df_te)
        oof_pred[te_idx] = mean
        if std is not None:
            oof_std[te_idx] = std
            has_std = True
        w = model.linear_weights_summary()
        fold_coefs.append(
            {
                "heldout_group": str(groups[te_idx[0]]),
                "n_test": int(len(te_idx)),
                "alpha_selected": model.alpha_selected,
                "w_prior_selected": model.w_prior_selected,
                "standardized_coefs": w["standardized_coefs"] if w else None,
            }
        )

    return oof_pred, (oof_std if has_std else None), fold_coefs


def _summarize_regression(
    df: pd.DataFrame,
    preds: np.ndarray,
    *,
    target_col: str = "delta_tm",
    pred_std: np.ndarray | None = None,
) -> dict[str, Any]:
    y = df[target_col].to_numpy(dtype=np.float64)
    groups = df["protein_group"].astype(str).tolist()
    res = regression_metrics(y, preds)
    res["paired_vs_zero"] = paired_comparison(y, preds, np.zeros_like(y), groups=groups)
    res["paired_vs_mono_k32"] = paired_comparison(y, preds, df["geo_mono_k32"].to_numpy(dtype=np.float64), groups=groups)
    res["paired_vs_asm_k32"] = paired_comparison(y, preds, df["geo_asm_k32"].to_numpy(dtype=np.float64), groups=groups)
    res["paired_vs_fx_mono"] = paired_comparison(y, preds, df["neg_fx_mono_ddg"].to_numpy(dtype=np.float64), groups=groups)
    res["paired_vs_fx_ac"] = paired_comparison(y, preds, df["neg_fx_ac_ddg_bind_total"].to_numpy(dtype=np.float64), groups=groups)
    res["paired_vs_fx_cage"] = paired_comparison(y, preds, df["neg_fx_cage_ddg_eff"].to_numpy(dtype=np.float64), groups=groups)
    if "s4346_ridge_sub" in df.columns and df["s4346_ridge_sub"].notna().all():
        res["paired_vs_s4346_ridge"] = paired_comparison(
            y, preds, df["s4346_ridge_sub"].to_numpy(dtype=np.float64), groups=groups
        )
    if pred_std is not None and np.isfinite(pred_std).all():
        abs_err = np.abs(y - preds)
        rho_u, p_u = spearmanr(pred_std, abs_err)
        cov95 = float(np.mean(abs_err <= 1.96 * pred_std))
        cov68 = float(np.mean(abs_err <= 1.00 * pred_std))
        res["uncertainty_calibration"] = {
            "mean_pred_std": round(float(np.mean(pred_std)), 3),
            "spearman_std_vs_abs_err": round(float(rho_u), 4) if np.isfinite(rho_u) else None,
            "spearman_p": float(p_u) if np.isfinite(p_u) else None,
            "coverage_68pct": round(cov68, 4),
            "coverage_95pct": round(cov95, 4),
        }
    return res


def _evaluate_s5_directional(
    df_s5: pd.DataFrame,
    scores: np.ndarray,
    *,
    neutral_offset: float | None = None,
) -> dict[str, Any]:
    labels = df_s5["binary_direction"].astype(float).astype(int).to_numpy()
    auc = roc_auc_binary(labels, scores)
    pos_scores = scores[labels == 1]
    neg_scores = scores[labels == 0]
    mw_greater = mannwhitneyu(pos_scores, neg_scores, alternative="greater")
    mw_two = mannwhitneyu(pos_scores, neg_scores, alternative="two-sided")

    def _thresh_stats(thr: float) -> dict[str, Any]:
        pred_pos = (scores > thr).astype(int)
        tp = int(np.sum((pred_pos == 1) & (labels == 1)))
        fp = int(np.sum((pred_pos == 1) & (labels == 0)))
        fn = int(np.sum((pred_pos == 0) & (labels == 1)))
        tn = int(np.sum((pred_pos == 0) & (labels == 0)))
        tpr = tp / max(1, tp + fn)
        tnr = tn / max(1, tn + fp)
        return {
            "threshold": round(float(thr), 4),
            "acc": round(float(np.mean(pred_pos == labels)), 4),
            "balanced_acc": round(0.5 * (tpr + tnr), 4),
            "confusion": {"tp": tp, "fp": fp, "fn": fn, "tn": tn},
        }

    out: dict[str, Any] = {
        "n": int(len(labels)),
        "n_pos": int(np.sum(labels == 1)),
        "n_neg": int(np.sum(labels == 0)),
        "roc_auc": auc,
        "mannwhitney_u_stat": float(mw_greater.statistic),
        "mannwhitney_p_one_sided": float(mw_greater.pvalue),
        "mannwhitney_p_two_sided": float(mw_two.pvalue),
        "gt0": _thresh_stats(0.0),
        "pos_mean_score": round(float(np.mean(pos_scores)), 3),
        "neg_mean_score": round(float(np.mean(neg_scores)), 3),
    }
    if neutral_offset is not None:
        out["gt_intercept_neutral"] = _thresh_stats(neutral_offset)
    return out


def _evaluate_external_regression(
    df_sub: pd.DataFrame,
    preds: np.ndarray,
    *,
    pred_std: np.ndarray | None = None,
) -> dict[str, Any]:
    valid = df_sub["delta_tm"].notna().to_numpy()
    sub = df_sub.iloc[valid].reset_index(drop=True)
    p = preds[valid]
    s = pred_std[valid] if pred_std is not None else None
    return _summarize_regression(sub, p, target_col="delta_tm", pred_std=s)


def _coefficient_stability(fold_coefs: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    if not fold_coefs or fold_coefs[0]["standardized_coefs"] is None:
        return {}
    feat_names = list(fold_coefs[0]["standardized_coefs"].keys())
    out: dict[str, dict[str, float]] = {}
    for f in feat_names:
        vals = [float(fc["standardized_coefs"][f]) for fc in fold_coefs]
        signs_nonzero = [np.sign(v) for v in vals if abs(v) > 1e-6]
        ref_sign = np.sign(float(np.mean(vals))) if abs(float(np.mean(vals))) > 1e-6 else 0.0
        sign_cons = float(np.mean([s == ref_sign for s in signs_nonzero])) if signs_nonzero else 1.0
        out[f] = {
            "mean": round(float(np.mean(vals)), 4),
            "std": round(float(np.std(vals)), 4),
            "min": round(float(np.min(vals)), 4),
            "max": round(float(np.max(vals)), 4),
            "sign_consistency": round(sign_cons, 4),
        }
    return out


def main() -> int:
    t0 = time.perf_counter()
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    df_raw = pd.read_csv(RESULTS_DIR / "ferritin_features_all.csv")
    df = compute_derived_feature_frame(df_raw)

    df_s1_all = df[df["cohort"] == "S1_SinglePoint_24mer"].reset_index(drop=True)
    df_s1_core = df[df["subcohort"] == "S1_PairedCore23"].reset_index(drop=True)
    df_s1_alt = df_s1_core.copy()
    has_alt = df_s1_alt["delta_tm_alt_thesis"].notna()
    df_s1_alt.loc[has_alt, "delta_tm"] = df_s1_alt.loc[has_alt, "delta_tm_alt_thesis"].astype(float)

    df_s2 = df[df["cohort"] == "S2_EcDps_12mer"].reset_index(drop=True)
    df_s3_exact = df[df["subcohort"] == "S3_Exact9"].reset_index(drop=True)
    df_s3_all = df[df["cohort"] == "S3_MultiPoint"].reset_index(drop=True)
    df_s4 = df[df["cohort"] == "S4_Secondary"].reset_index(drop=True)
    df_s5 = df[df["cohort"] == "S5_TierBC_Directional"].reset_index(drop=True)

    # 1. Model Family Benchmark on S1_All30, S1_PairedCore23, and AltThesis sensitivity
    cv_schemes = ["SiteGroup", "LOPO", "LOFO", "StudyGroup"]
    model_cv_results: dict[str, Any] = {"S1_All30": {}, "S1_PairedCore23": {}, "S1_PairedCore23_AltThesis": {}}
    oof_store_s1_all: dict[tuple[str, str], np.ndarray] = {}
    oof_std_store_s1_all: dict[tuple[str, str], np.ndarray | None] = {}
    coef_stability_s1_all: dict[tuple[str, str], dict[str, Any]] = {}

    for label, fam_name, res_flag, feats in MODEL_FAMILIES:
        for cohort_label, df_c in [
            ("S1_All30", df_s1_all),
            ("S1_PairedCore23", df_s1_core),
            ("S1_PairedCore23_AltThesis", df_s1_alt),
        ]:
            fam_res: dict[str, Any] = {}
            for scheme in cv_schemes:
                preds, stds, fold_coefs = run_oof_cv(
                    df_c,
                    group_kind=scheme,
                    model_family=fam_name,
                    residual_anchor=res_flag,
                    feature_names=feats,
                    train_cohort_label=cohort_label,
                )
                fam_res[scheme] = _summarize_regression(df_c, preds, pred_std=stds)
                if cohort_label == "S1_All30":
                    oof_store_s1_all[(label, scheme)] = preds
                    oof_std_store_s1_all[(label, scheme)] = stds
                    stab = _coefficient_stability(fold_coefs)
                    if stab:
                        coef_stability_s1_all[(label, scheme)] = stab
            model_cv_results[cohort_label][label] = fam_res

    print("=== Model Family CV Comparison on S1_All30 (N=30) ===")
    for label, _, _, _ in MODEL_FAMILIES:
        sg = model_cv_results["S1_All30"][label]["SiteGroup"]
        lp = model_cv_results["S1_All30"][label]["LOPO"]
        lf = model_cv_results["S1_All30"][label]["LOFO"]
        print(
            f"{label:28s} | SiteGrp MAE={sg['mae_c']:.2f} r={sg['pearson_r']:+.3f} "
            f"| LOPO MAE={lp['mae_c']:.2f} r={lp['pearson_r']:+.3f} "
            f"| LOFO MAE={lf['mae_c']:.2f} r={lf['pearson_r']:+.3f}"
        )

    # 2. Mechanistic Ablation Study (Ab1..Ab6) with PriorAnchoredNonNegHuber
    ablation_results: dict[str, Any] = {}
    for ab_name, ab_feats in ABLATION_SETS.items():
        ab_entry: dict[str, Any] = {"features": ab_feats, "n_features": len(ab_feats)}
        for scheme in ("SiteGroup", "LOPO", "LOFO"):
            preds, stds, _ = run_oof_cv(
                df_s1_all,
                group_kind=scheme,
                model_family="PriorAnchoredNonNegHuber",
                residual_anchor=False,
                feature_names=ab_feats,
                train_cohort_label="S1_All30",
            )
            ab_entry[f"S1_All30_{scheme}"] = _summarize_regression(df_s1_all, preds, pred_std=stds)

        fitted_ab = fit_ferritin_hybrid(
            df_s1_all,
            model_family="PriorAnchoredNonNegHuber",
            feature_names=ab_feats,
            residual_anchor=False,
            train_cohort_label="S1_All30",
        )
        ab_entry["linear_weights"] = fitted_ab.linear_weights_summary()
        ab_entry["alpha_selected"] = fitted_ab.alpha_selected
        ab_entry["w_prior_selected"] = fitted_ab.w_prior_selected
        p_s2, s_s2 = fitted_ab.predict_with_std(df_s2)
        p_s3e, s_s3e = fitted_ab.predict_with_std(df_s3_exact)
        p_s5, _ = fitted_ab.predict_with_std(df_s5)
        b0_ab = fitted_ab.linear_weights_summary()["std_intercept_b0"]
        ab_entry["ZeroShot_S2_EcDps_Single8"] = _evaluate_external_regression(df_s2, p_s2, pred_std=s_s2)
        ab_entry["ZeroShot_S3_Exact9"] = _evaluate_external_regression(df_s3_exact, p_s3e, pred_std=s_s3e)
        ab_entry["ZeroShot_S5_Directional13"] = _evaluate_s5_directional(df_s5, p_s5, neutral_offset=b0_ab)
        ablation_results[ab_name] = ab_entry

    print("\n=== Mechanistic Ablation Study (PriorAnchoredNonNegHuber, S1_All30 CV + External Zero-Shot) ===")
    for ab_name, ab_entry in ablation_results.items():
        lp = ab_entry["S1_All30_LOPO"]
        lf = ab_entry["S1_All30_LOFO"]
        s2m = ab_entry["ZeroShot_S2_EcDps_Single8"]
        s3m = ab_entry["ZeroShot_S3_Exact9"]
        s5m = ab_entry["ZeroShot_S5_Directional13"]
        print(
            f"{ab_name:34s} | LOPO {lp['mae_c']:.2f}/{lp['pearson_r']:+.3f} "
            f"| LOFO {lf['mae_c']:.2f}/{lf['pearson_r']:+.3f} "
            f"| S2 {s2m['mae_c']:.2f}/{s2m['pearson_r']:+.3f} "
            f"| S3 {s3m['mae_c']:.2f}/{s3m['pearson_r']:+.3f} "
            f"| S5 AUC={s5m['roc_auc']:.3f}"
        )

    # 3. Deploy Primary Model (FerritinGeoHybrid_NonNegMAP) + Comparator Models
    deployed_model = fit_ferritin_hybrid(
        df_s1_all,
        model_family="PriorAnchoredNonNegHuber",
        feature_names=HYBRID_FEATURES_5,
        residual_anchor=False,
        train_cohort_label="S1_All30",
    )
    model_path = MODELS_DIR / "ferritin_geohybrid.joblib"
    artifact_sha = save_predictor(deployed_model, model_path)
    reloaded_model = load_predictor(model_path)
    b0_dep = float(reloaded_model.linear_weights_summary()["std_intercept_b0"])

    # Comparator models for full zero-shot audit
    sensitivity_core23_model = fit_ferritin_hybrid(
        df_s1_core,
        model_family="PriorAnchoredNonNegHuber",
        feature_names=HYBRID_FEATURES_5,
        residual_anchor=False,
        train_cohort_label="S1_PairedCore23",
    )
    residual_model = fit_ferritin_hybrid(
        df_s1_all,
        model_family="PriorAnchoredNonNegHuber",
        feature_names=RESIDUAL_HYBRID_FEATURES_4,
        residual_anchor=True,
        train_cohort_label="S1_All30",
    )
    zeroprior_model = fit_ferritin_hybrid(
        df_s1_all,
        model_family="ZeroPriorNonNegHuber",
        feature_names=HYBRID_FEATURES_5,
        residual_anchor=False,
        train_cohort_label="S1_All30",
    )
    bayes7_unconstrained = fit_ferritin_hybrid(
        df_s1_all,
        model_family="BayesianRidge",
        feature_names=CORE_FEATURES_7,
        residual_anchor=False,
        train_cohort_label="S1_All30",
    )

    p_s2_dep, std_s2_dep = reloaded_model.predict_with_std(df_s2)
    p_s3e_dep, std_s3e_dep = reloaded_model.predict_with_std(df_s3_exact)
    p_s3a_dep, std_s3a_dep = reloaded_model.predict_with_std(df_s3_all)
    p_s4_dep, std_s4_dep = reloaded_model.predict_with_std(df_s4)
    p_s5_dep, _ = reloaded_model.predict_with_std(df_s5)

    p_s2_c23, _ = sensitivity_core23_model.predict_with_std(df_s2)
    p_s3e_c23, _ = sensitivity_core23_model.predict_with_std(df_s3_exact)
    p_s5_c23, _ = sensitivity_core23_model.predict_with_std(df_s5)
    b0_c23 = float(sensitivity_core23_model.linear_weights_summary()["std_intercept_b0"])

    p_s2_res, std_s2_res = residual_model.predict_with_std(df_s2)
    p_s3e_res, std_s3e_res = residual_model.predict_with_std(df_s3_exact)
    p_s5_res, _ = residual_model.predict_with_std(df_s5)
    b0_res = float(residual_model.linear_weights_summary()["std_intercept_b0"])

    p_s2_zp, std_s2_zp = zeroprior_model.predict_with_std(df_s2)
    p_s3e_zp, std_s3e_zp = zeroprior_model.predict_with_std(df_s3_exact)
    p_s5_zp, _ = zeroprior_model.predict_with_std(df_s5)
    b0_zp = float(zeroprior_model.linear_weights_summary()["std_intercept_b0"])

    s4_sub_eval: dict[str, Any] = {}
    for subc, df_subc in df_s4.groupby("subcohort", sort=False):
        p_subc, std_subc = reloaded_model.predict_with_std(df_subc)
        s4_sub_eval[str(subc)] = {
            "n": int(len(df_subc)),
            "delta_tm_obs": [float(x) for x in df_subc["delta_tm"]],
            "pred_hybrid": [round(float(x), 3) for x in p_subc],
            "metrics": _evaluate_external_regression(df_subc, p_subc, pred_std=std_subc) if len(df_subc) >= 3 else None,
        }

    # Pooled validation summaries:
    # (A) External quantitative cohorts S2 + S3_Exact9 (n=17)
    # (B) All independent quantitative evaluations S1_All30(LOPO) + S2 + S3_Exact9 (n=47)
    df_ext_quant = pd.concat([df_s2, df_s3_exact], ignore_index=True)
    p_ext_quant = np.concatenate([p_s2_dep, p_s3e_dep])
    std_ext_quant = np.concatenate([std_s2_dep, std_s3e_dep])

    df_pool47 = pd.concat([df_s1_all, df_s2, df_s3_exact], ignore_index=True)
    p_pool47 = np.concatenate([oof_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")], p_s2_dep, p_s3e_dep])
    std_pool47 = np.concatenate([oof_std_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")], std_s2_dep, std_s3e_dep])

    p_ext_res = np.concatenate([p_s2_res, p_s3e_res])
    p_pool47_res = np.concatenate([oof_store_s1_all[("Residual_NonNegMAP", "LOPO")], p_s2_res, p_s3e_res])

    s5_baselines: dict[str, Any] = {
        "Ferritin_GeoHybrid_Deployed_S1All30": _evaluate_s5_directional(df_s5, p_s5_dep, neutral_offset=b0_dep),
        "Ferritin_GeoHybrid_Sens_Core23": _evaluate_s5_directional(df_s5, p_s5_c23, neutral_offset=b0_c23),
        "Residual_NonNegMAP_S1All30": _evaluate_s5_directional(df_s5, p_s5_res, neutral_offset=b0_res),
        "ZeroPrior_NonNegHuber_S1All30": _evaluate_s5_directional(df_s5, p_s5_zp, neutral_offset=b0_zp),
        "BayesianRidge_FreeSign7_Unconstrained": _evaluate_s5_directional(df_s5, bayes7_unconstrained.predict(df_s5)),
        "Monomer_GeoDTm3D_k32": _evaluate_s5_directional(df_s5, df_s5["geo_mono_k32"].to_numpy(dtype=np.float64)),
        "Assembly_GeoDTm3D_k32": _evaluate_s5_directional(df_s5, df_s5["geo_asm_k32"].to_numpy(dtype=np.float64)),
        "FoldX_mono_neg_ddG": _evaluate_s5_directional(df_s5, df_s5["neg_fx_mono_ddg"].to_numpy(dtype=np.float64)),
        "FoldX_AC_neg_bind_ddG": _evaluate_s5_directional(
            df_s5, df_s5["neg_fx_ac_ddg_bind_total"].to_numpy(dtype=np.float64)
        ),
        "FoldX_cage_neg_eff_ddG": _evaluate_s5_directional(df_s5, df_s5["neg_fx_cage_ddg_eff"].to_numpy(dtype=np.float64)),
        "ESM2_WT_marginal_LLR": _evaluate_s5_directional(df_s5, df_s5["esm2_wt_marginal_llr"].to_numpy(dtype=np.float64)),
    }

    # 4. Build per-case prediction CSV across all 75 rows
    p_all_dep, std_all_dep = reloaded_model.predict_with_std(df)
    p_all_res, _ = residual_model.predict_with_std(df)
    p_all_c23, _ = sensitivity_core23_model.predict_with_std(df)
    s1_case_to_idx = {cid: i for i, cid in enumerate(df_s1_all["case_id"].tolist())}
    oof_sg = np.full(len(df), np.nan, dtype=np.float64)
    oof_lopo = np.full(len(df), np.nan, dtype=np.float64)
    oof_lofo = np.full(len(df), np.nan, dtype=np.float64)
    oof_lopo_std = np.full(len(df), np.nan, dtype=np.float64)
    for idx, cid in enumerate(df["case_id"].tolist()):
        if cid in s1_case_to_idx:
            j = s1_case_to_idx[cid]
            oof_sg[idx] = oof_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "SiteGroup")][j]
            oof_lopo[idx] = oof_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")][j]
            oof_lofo[idx] = oof_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOFO")][j]
            if oof_std_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")] is not None:
                oof_lopo_std[idx] = oof_std_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")][j]

    out_pred_df = pd.DataFrame(
        {
            "case_id": df["case_id"],
            "cohort": df["cohort"],
            "subcohort": df["subcohort"],
            "template_key": df["template_key"],
            "protein_group": df["protein_group"],
            "oligomer_family": df["oligomer_family"],
            "study_id": df["study_id"],
            "mutation_canon": df["mutation_canon"],
            "pH": df["pH"],
            "delta_tm": df["delta_tm"],
            "delta_tm_alt_thesis": df["delta_tm_alt_thesis"],
            "binary_direction": df["binary_direction"],
            "geo_mono_k32": df["geo_mono_k32"],
            "geo_asm_k32": df["geo_asm_k32"],
            "neg_fx_cage_ddg_eff": df["neg_fx_cage_ddg_eff"],
            "neg_fx_ac_ddg_bind_c3": df["neg_fx_ac_ddg_bind_c3"],
            "neg_fx_ac_ddg_bind_nonpolar": df["neg_fx_ac_ddg_bind_nonpolar"],
            "neg_geo_asm_enc_shift_k32": df["neg_geo_asm_enc_shift_k32"],
            "neg_fx_mono_ddg": df["neg_fx_mono_ddg"],
            "neg_fx_ac_ddg_bind_total": df["neg_fx_ac_ddg_bind_total"],
            "esm2_wt_marginal_llr": df["esm2_wt_marginal_llr"],
            "s4346_ridge_sub": df["s4346_ridge_sub"],
            "hybrid_oof_sitegroup": oof_sg,
            "hybrid_oof_lopo": oof_lopo,
            "hybrid_oof_lofo": oof_lofo,
            "hybrid_oof_lopo_std": oof_lopo_std,
            "hybrid_deployed_s1all30": p_all_dep,
            "hybrid_deployed_std": std_all_dep,
            "hybrid_residual_s1all30": p_all_res,
            "hybrid_sens_core23": p_all_c23,
        }
    )
    pred_csv_path = RESULTS_DIR / "hybrid_predictions_all.csv"
    out_pred_df.to_csv(pred_csv_path, index=False)

    # Error breakdown by S1 biological site class / worst residuals under LOPO
    s1_diag = df_s1_all.copy()
    s1_diag["oof_lopo"] = oof_store_s1_all[("FerritinGeoHybrid_NonNegMAP", "LOPO")]
    s1_diag["abs_err_lopo"] = np.abs(s1_diag["delta_tm"] - s1_diag["oof_lopo"])
    s1_diag["abs_err_mono32"] = np.abs(s1_diag["delta_tm"] - s1_diag["geo_mono_k32"])
    top_errors = []
    for _, row in s1_diag.sort_values("abs_err_lopo", ascending=False).head(10).iterrows():
        top_errors.append(
            {
                "case_id": str(row["case_id"]),
                "protein_group": str(row["protein_group"]),
                "mutation_canon": str(row["mutation_canon"]),
                "delta_tm": float(row["delta_tm"]),
                "mono_geo32": round(float(row["geo_mono_k32"]), 3),
                "asm_geo32": round(float(row["geo_asm_k32"]), 3),
                "oof_lopo": round(float(row["oof_lopo"]), 3),
                "err_lopo": round(float(row["oof_lopo"] - row["delta_tm"]), 3),
                "abs_err_lopo": round(float(row["abs_err_lopo"]), 3),
                "abs_err_mono32": round(float(row["abs_err_mono32"]), 3),
            }
        )

    summary: dict[str, Any] = {
        "elapsed_sec": round(time.perf_counter() - t0, 2),
        "primary_model": {
            "name": "Ferritin-GeoHybrid (PriorAnchoredNonNegHuber, 5 physical pillars)",
            "artifact_path": str(model_path),
            "artifact_sha256": artifact_sha,
            "train_cohort": "S1_All30",
            "train_n": int(len(df_s1_all)),
            "features": HYBRID_FEATURES_5,
            "alpha_selected": reloaded_model.alpha_selected,
            "w_prior_selected": reloaded_model.w_prior_selected,
            "deployed_weights_s1all30": reloaded_model.linear_weights_summary(),
            "sensitivity_weights_core23": sensitivity_core23_model.linear_weights_summary(),
            "residual_companion_weights_s1all30": residual_model.linear_weights_summary(),
            "zeroprior_weights_s1all30": zeroprior_model.linear_weights_summary(),
            "unconstrained_bayes7_weights_s1all30": bayes7_unconstrained.linear_weights_summary(),
            "lopo_coefficient_stability_s1all30": coef_stability_s1_all.get(("FerritinGeoHybrid_NonNegMAP", "LOPO"), {}),
            "lofo_coefficient_stability_s1all30": coef_stability_s1_all.get(("FerritinGeoHybrid_NonNegMAP", "LOFO"), {}),
        },
        "cv_model_family_comparison": model_cv_results,
        "ablation_study": ablation_results,
        "external_zero_shot_transfer": {
            "S2_EcDps_Single8": {
                "deployed_S1All30": _evaluate_external_regression(df_s2, p_s2_dep, pred_std=std_s2_dep),
                "sensitivity_Core23": _evaluate_external_regression(df_s2, p_s2_c23),
                "residual_companion": _evaluate_external_regression(df_s2, p_s2_res, pred_std=std_s2_res),
                "zeroprior_companion": _evaluate_external_regression(df_s2, p_s2_zp, pred_std=std_s2_zp),
                "unconstrained_bayes7": _evaluate_external_regression(df_s2, bayes7_unconstrained.predict(df_s2)),
            },
            "S3_Exact9": {
                "deployed_S1All30": _evaluate_external_regression(df_s3_exact, p_s3e_dep, pred_std=std_s3e_dep),
                "sensitivity_Core23": _evaluate_external_regression(df_s3_exact, p_s3e_c23),
                "residual_companion": _evaluate_external_regression(df_s3_exact, p_s3e_res, pred_std=std_s3e_res),
                "zeroprior_companion": _evaluate_external_regression(df_s3_exact, p_s3e_zp, pred_std=std_s3e_zp),
                "unconstrained_bayes7": _evaluate_external_regression(df_s3_exact, bayes7_unconstrained.predict(df_s3_exact)),
            },
            "S3_All16": {
                "deployed_S1All30": _evaluate_external_regression(df_s3_all, p_s3a_dep, pred_std=std_s3a_dep),
            },
            "S4_Secondary8": {
                "deployed_S1All30": _evaluate_external_regression(df_s4, p_s4_dep, pred_std=std_s4_dep),
                "by_subcohort": s4_sub_eval,
            },
            "S5_Unseen13_Directional": s5_baselines,
            "Pooled_External_S2_S3Exact9_N17": {
                "deployed_S1All30": _summarize_regression(df_ext_quant, p_ext_quant, pred_std=std_ext_quant),
                "residual_companion": _summarize_regression(df_ext_quant, p_ext_res),
            },
            "Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47": {
                "deployed_S1All30": _summarize_regression(df_pool47, p_pool47, pred_std=std_pool47),
                "residual_companion": _summarize_regression(df_pool47, p_pool47_res),
            },
        },
        "s1_lopo_largest_errors": top_errors,
    }

    summary_path = RESULTS_DIR / "hybrid_benchmark_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("\n=== Deployed Ferritin-GeoHybrid (S1_All30) External Zero-Shot & Pooled Summary ===")
    s2_ext = summary["external_zero_shot_transfer"]["S2_EcDps_Single8"]["deployed_S1All30"]
    s3_ext = summary["external_zero_shot_transfer"]["S3_Exact9"]["deployed_S1All30"]
    s4_ext = summary["external_zero_shot_transfer"]["S4_Secondary8"]["deployed_S1All30"]
    s5_ext = s5_baselines["Ferritin_GeoHybrid_Deployed_S1All30"]
    p17 = summary["external_zero_shot_transfer"]["Pooled_External_S2_S3Exact9_N17"]["deployed_S1All30"]
    p47 = summary["external_zero_shot_transfer"]["Pooled_AllUnseen_S1LOPO_S2_S3Exact9_N47"]["deployed_S1All30"]
    print(f"S2_EcDps_Single8 (n={s2_ext['n']}): MAE={s2_ext['mae_c']:.2f}, RMSE={s2_ext['rmse_c']:.2f}, r={s2_ext['pearson_r']:+.3f}, rho={s2_ext['spearman_rho']:+.3f}, sign={s2_ext['sign_accuracy']:.3f}")
    print(f"S3_Exact9        (n={s3_ext['n']}): MAE={s3_ext['mae_c']:.2f}, RMSE={s3_ext['rmse_c']:.2f}, r={s3_ext['pearson_r']:+.3f}, rho={s3_ext['spearman_rho']:+.3f}, sign={s3_ext['sign_accuracy']:.3f}")
    print(f"S4_Secondary8    (n={s4_ext['n']}): MAE={s4_ext['mae_c']:.2f}, RMSE={s4_ext['rmse_c']:.2f}, r={s4_ext['pearson_r']:+.3f}, rho={s4_ext['spearman_rho']:+.3f}")
    print(
        f"S5_Directional   (n={s5_ext['n']}): ROC-AUC={s5_ext['roc_auc']:.3f}, "
        f"MW p(1s)={s5_ext['mannwhitney_p_one_sided']:.4f}, "
        f"acc(gt0)={s5_ext['gt0']['acc']:.3f}, acc(neutral)={s5_ext['gt_intercept_neutral']['acc']:.3f}"
    )
    print(
        f"Pooled External (N=17): MAE={p17['mae_c']:.2f} r={p17['pearson_r']:+.3f} "
        f"| vs Mono32 dMAE={p17['paired_vs_mono_k32']['mean_abs_err_diff_c']:+.2f} (p={p17['paired_vs_mono_k32']['wilcoxon_p']:.4f})"
    )
    print(
        f"Pooled All-Unseen (N=47): MAE={p47['mae_c']:.2f} r={p47['pearson_r']:+.3f} "
        f"| vs Mono32 dMAE={p47['paired_vs_mono_k32']['mean_abs_err_diff_c']:+.2f} (p={p47['paired_vs_mono_k32']['wilcoxon_p']:.4f}) "
        f"| vs Zero p={p47['paired_vs_zero']['wilcoxon_p']:.4f}"
    )
    print(f"\nSaved {model_path}, {summary_path}, and {pred_csv_path} in {summary['elapsed_sec']}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
