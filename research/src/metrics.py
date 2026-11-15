"""Evaluation metrics, paired comparisons, and bootstrap confidence intervals."""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy.stats import pearsonr, spearmanr, wilcoxon


def regression_metrics(y_true: Iterable[float], y_pred: Iterable[float]) -> dict[str, float | int | None]:
    """Compute descriptive regression and directional classification metrics."""
    y = np.asarray(list(y_true), dtype=float)
    p = np.asarray(list(y_pred), dtype=float)
    if len(y) == 0:
        raise ValueError("Empty evaluation vector")
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(p))):
        raise ValueError("Non-finite values in metrics input")

    err = p - y
    mae = float(np.mean(np.abs(err)))
    rmse = float(math.sqrt(float(np.mean(err ** 2))))
    bias = float(np.mean(err))
    zero_mae = float(np.mean(np.abs(y)))
    zero_rmse = float(math.sqrt(float(np.mean(y ** 2))))

    pr = None
    pr_p = None
    if len(y) >= 2 and np.std(y) > 0 and np.std(p) > 0:
        res = pearsonr(y, p)
        pr = float(getattr(res, "statistic", res[0]))
        pr_p = float(getattr(res, "pvalue", res[1]))

    sr = None
    sr_p = None
    if len(y) >= 2 and np.std(y) > 0 and np.std(p) > 0:
        res_s = spearmanr(y, p)
        sr = float(getattr(res_s, "statistic", res_s[0]))
        sr_p = float(getattr(res_s, "pvalue", res_s[1]))

    nz_mask = y != 0.0
    n_nz = int(np.sum(nz_mask))
    sign_correct = int(np.sum((y[nz_mask] > 0.0) == (p[nz_mask] > 0.0))) if n_nz > 0 else 0
    sign_acc = float(sign_correct / n_nz) if n_nz > 0 else None

    # Stabilizing detection: true > 0 vs pred > 0
    pos_true = y > 0.0
    pos_pred = p > 0.0
    tp = int(np.sum(pos_true & pos_pred))
    fp = int(np.sum((~pos_true) & pos_pred))
    fn = int(np.sum(pos_true & (~pos_pred)))

    return {
        "n": int(len(y)),
        "mae_c": mae,
        "rmse_c": rmse,
        "bias_c": bias,
        "zero_mae_c": zero_mae,
        "zero_rmse_c": zero_rmse,
        "mae_skill_vs_zero": float(1.0 - mae / zero_mae) if zero_mae > 0 else None,
        "rmse_skill_vs_zero": float(1.0 - rmse / zero_rmse) if zero_rmse > 0 else None,
        "pearson_r": pr,
        "pearson_p": pr_p,
        "spearman_rho": sr,
        "spearman_p": sr_p,
        "n_nonzero": n_nz,
        "sign_correct": sign_correct,
        "sign_accuracy": sign_acc,
        "pos_tp": tp,
        "pos_fp": fp,
        "pos_fn": fn,
    }


def paired_comparison(
    y_true: Iterable[float],
    pred_new: Iterable[float],
    pred_ref: Iterable[float],
    *,
    groups: Iterable[str] | None = None,
    n_boot: int = 2000,
    seed: int = 20261115,
) -> dict[str, float | int | None]:
    """Compare absolute error of pred_new minus pred_ref (negative = new improves)."""
    y = np.asarray(list(y_true), dtype=float)
    pn = np.asarray(list(pred_new), dtype=float)
    pr = np.asarray(list(pred_ref), dtype=float)
    ae_new = np.abs(pn - y)
    ae_ref = np.abs(pr - y)
    diff = ae_new - ae_ref  # < 0 means new has lower error
    mean_diff = float(np.mean(diff))
    wins = int(np.sum(diff < 0.0))
    losses = int(np.sum(diff > 0.0))
    ties = int(np.sum(diff == 0.0))

    wilcox_p = None
    nonzero = diff[diff != 0.0]
    if len(nonzero) >= 5:
        try:
            wilcox_p = float(wilcoxon(diff).pvalue)
        except Exception:
            wilcox_p = None

    rng = np.random.default_rng(seed)
    boot_means = []
    if groups is None:
        idx_pool = np.arange(len(y))
        for _ in range(n_boot):
            sample = rng.choice(idx_pool, size=len(y), replace=True)
            boot_means.append(float(np.mean(diff[sample])))
    else:
        grp_arr = np.asarray(list(groups))
        unique_g = np.unique(grp_arr)
        grp_indices = [np.where(grp_arr == g)[0] for g in unique_g]
        for _ in range(n_boot):
            chosen_g = rng.integers(0, len(unique_g), size=len(unique_g))
            sample = np.concatenate([grp_indices[g] for g in chosen_g])
            boot_means.append(float(np.mean(diff[sample])))

    ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
    return {
        "n": int(len(y)),
        "mean_abs_err_diff_c": mean_diff,
        "ci95_low_c": float(ci_lo),
        "ci95_high_c": float(ci_hi),
        "wins": wins,
        "losses": losses,
        "ties": ties,
        "wilcoxon_p": wilcox_p,
        "mae_new_c": float(np.mean(ae_new)),
        "mae_ref_c": float(np.mean(ae_ref)),
        "rmse_new_c": float(math.sqrt(float(np.mean((pn - y) ** 2)))),
        "rmse_ref_c": float(math.sqrt(float(np.mean((pr - y) ** 2)))),
    }


def roc_auc_binary(labels: Iterable[int], scores: Iterable[float]) -> float | None:
    """Compute ROC AUC without sklearn dependency via Mann-Whitney rank statistic."""
    y = np.asarray(list(labels), dtype=int)
    s = np.asarray(list(scores), dtype=float)
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    wins = 0.0
    for vp in pos:
        wins += float(np.sum(vp > neg)) + 0.5 * float(np.sum(vp == neg))
    return float(wins / (len(pos) * len(neg)))
