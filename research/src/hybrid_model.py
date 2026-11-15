"""Domain-specialized Ferritin-GeoHybrid model built on Assembly-GeoDTm-3D."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import BayesianRidge, HuberRegressor, Ridge
from sklearn.model_selection import GroupKFold, KFold, LeaveOneOut
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

CORE_FEATURES_7 = [
    "geo_asm_k32",
    "geo_iface_delta_k32",
    "neg_fx_mono_ddg",
    "fx_sub_vs_mono_diff",
    "neg_fx_ac_ddg_bind_c3",
    "neg_fx_ac_ddg_bind_nonpolar",
    "near_cofactor_8a",
]

HYBRID_FEATURES_5 = [
    "geo_asm_k32",
    "neg_fx_cage_ddg_eff",
    "neg_fx_ac_ddg_bind_c3",
    "neg_fx_ac_ddg_bind_nonpolar",
    "neg_geo_asm_enc_shift_k32",
]

RESIDUAL_HYBRID_FEATURES_4 = [
    "neg_fx_cage_ddg_eff",
    "neg_fx_ac_ddg_bind_c3",
    "neg_fx_ac_ddg_bind_nonpolar",
    "neg_geo_asm_enc_shift_k32",
]

EXTENDED_FEATURES_9 = CORE_FEATURES_7 + [
    "esm2_wt_marginal_llr",
    "delta_abs_charge",
]

ABLATION_SETS: dict[str, list[str]] = {
    "Ab1_GeoMono3D_Calibrated": ["geo_mono_k32"],
    "Ab2_GeoAsm3D_Plus_EncShift": ["geo_asm_k32", "neg_geo_asm_enc_shift_k32"],
    "Ab3_Plus_CageEffectiveStability": [
        "geo_asm_k32",
        "neg_geo_asm_enc_shift_k32",
        "neg_fx_cage_ddg_eff",
    ],
    "Ab4_Plus_C3Interface": [
        "geo_asm_k32",
        "neg_geo_asm_enc_shift_k32",
        "neg_fx_cage_ddg_eff",
        "neg_fx_ac_ddg_bind_c3",
    ],
    "Ab5_Full_Ferritin_GeoHybrid_5feat": HYBRID_FEATURES_5,
    "Ab6_Extended_Plus_ESM2_LLR": HYBRID_FEATURES_5 + ["esm2_wt_marginal_llr"],
}


def compute_derived_feature_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Compute canonical derived biophysical columns cleanly on a copy of df."""
    out = df.copy()
    for col in (
        "geo_mono_k32",
        "geo_asm_k32",
        "geo_asm_k48",
        "geo_iface_delta_k32",
        "geo_ph_shift_asm",
        "geo_mono_enc_shift_k32",
        "geo_asm_enc_shift_k32",
        "fx_mono_ddg",
        "neg_fx_mono_ddg",
        "fx_sub_vs_mono_diff",
        "fx_cage_ddg_eff",
        "neg_fx_cage_ddg_eff",
        "fx_ac_ddg_bind_c3",
        "neg_fx_ac_ddg_bind_c3",
        "fx_ac_ddg_bind_solv_hyd",
        "fx_ac_ddg_bind_vdw",
        "near_cofactor_8a",
        "esm2_wt_marginal_llr",
        "delta_abs_charge",
        "inter_ca10_fraction",
    ):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype(float)

    if "neg_fx_mono_ddg" not in out.columns and "fx_mono_ddg" in out.columns:
        out["neg_fx_mono_ddg"] = -out["fx_mono_ddg"]
    if "neg_fx_cage_ddg_eff" not in out.columns and "fx_cage_ddg_eff" in out.columns:
        out["neg_fx_cage_ddg_eff"] = -out["fx_cage_ddg_eff"]
    if "neg_fx_ac_ddg_bind_c3" not in out.columns and "fx_ac_ddg_bind_c3" in out.columns:
        out["neg_fx_ac_ddg_bind_c3"] = -out["fx_ac_ddg_bind_c3"]
    if "fx_ac_ddg_bind_nonpolar" not in out.columns:
        out["fx_ac_ddg_bind_nonpolar"] = out["fx_ac_ddg_bind_solv_hyd"] + out["fx_ac_ddg_bind_vdw"]
    out["neg_fx_ac_ddg_bind_nonpolar"] = -out["fx_ac_ddg_bind_nonpolar"]
    if "geo_asm_enc_shift_k32" in out.columns:
        out["neg_geo_asm_enc_shift_k32"] = -out["geo_asm_enc_shift_k32"]
    return out


def extract_feature_matrix(df: pd.DataFrame, feature_names: list[str]) -> np.ndarray:
    work = compute_derived_feature_frame(df)
    mat = work[feature_names].to_numpy(dtype=np.float64)
    if not np.isfinite(mat).all():
        bad = [c for c in feature_names if not np.isfinite(work[c]).all()]
        raise ValueError(f"Non-finite feature columns encountered: {bad}")
    return mat


def _inner_cv_for_train(train_groups: np.ndarray | None, n_train: int):
    if train_groups is not None:
        uniq = np.unique(train_groups)
        if len(uniq) >= 3:
            return list(GroupKFold(n_splits=min(5, len(uniq))).split(np.zeros((n_train, 1)), groups=train_groups))
    if n_train >= 10:
        return list(KFold(n_splits=min(5, n_train), shuffle=True, random_state=20261115).split(np.zeros((n_train, 1))))
    return list(LeaveOneOut().split(np.zeros((n_train, 1))))


def _huber_loss_and_grad(r: np.ndarray, delta: float) -> tuple[float, np.ndarray]:
    a = np.abs(r)
    quad = a <= delta
    loss = np.where(quad, 0.5 * r * r, delta * (a - 0.5 * delta))
    grad_r = np.where(quad, r, delta * np.sign(r))
    return float(np.sum(loss)), grad_r


class NonNegHuberMAPRegressor(BaseEstimator, RegressorMixin):
    """L2 prior-anchored non-negative Huber regressor with Laplace predictive std."""

    def __init__(
        self,
        *,
        alpha: float = 20.0,
        w_prior: float = 0.0,
        delta: float = 8.0,
    ) -> None:
        self.alpha = alpha
        self.w_prior = w_prior
        self.delta = delta

    def fit(self, X: np.ndarray, y: np.ndarray) -> "NonNegHuberMAPRegressor":
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        n, p = X.shape
        prior_vec = np.full(p, float(self.w_prior), dtype=np.float64)

        def obj_and_grad(params: np.ndarray) -> tuple[float, np.ndarray]:
            b0 = params[0]
            w = params[1:]
            r = (X @ w + b0) - y
            dl, dr = _huber_loss_and_grad(r, self.delta)
            pen = 0.5 * self.alpha * float(np.sum((w - prior_vec) ** 2))
            g_b0 = float(np.sum(dr))
            g_w = X.T @ dr + self.alpha * (w - prior_vec)
            return dl + pen, np.concatenate(([g_b0], g_w))

        init = np.concatenate(([float(np.mean(y))], prior_vec))
        bounds = [(None, None)] + [(0.0, None)] * p
        res = minimize(obj_and_grad, init, jac=True, method="L-BFGS-B", bounds=bounds)
        self.intercept_ = float(res.x[0])
        self.coef_ = np.asarray(res.x[1:], dtype=np.float64)

        # IRLS/Laplace posterior covariance for predictive standard deviation
        r_fit = y - (self.intercept_ + X @ self.coef_)
        irls_w = np.where(np.abs(r_fit) <= self.delta, 1.0, self.delta / np.maximum(np.abs(r_fit), 1e-8))
        X_aug = np.column_stack([np.ones(n, dtype=np.float64), X])
        hess = (X_aug * irls_w[:, None]).T @ X_aug
        pen_diag = np.diag(np.concatenate(([0.0], np.full(p, max(float(self.alpha), 1e-6)))))
        hess_reg = hess + pen_diag + 1e-6 * np.eye(p + 1)
        cov_inv_hess = np.linalg.inv(hess_reg)
        df_eff = float(np.trace(cov_inv_hess @ hess))
        rss_weighted = float(np.sum(irls_w * (r_fit ** 2)))
        mad = float(np.median(np.abs(r_fit - np.median(r_fit)))) * 1.4826
        sigma_hat = max(float(np.sqrt(rss_weighted / max(1.0, n - df_eff))), mad, 1.0)
        self.sigma_hat_ = sigma_hat
        self.param_cov_ = (sigma_hat ** 2) * cov_inv_hess
        return self

    def predict(self, X: np.ndarray, return_std: bool = False) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        X = np.asarray(X, dtype=np.float64)
        mean = self.intercept_ + X @ self.coef_
        if not return_std:
            return mean
        X_aug = np.column_stack([np.ones(len(X), dtype=np.float64), X])
        epistemic_var = np.sum((X_aug @ self.param_cov_) * X_aug, axis=1)
        std = np.sqrt(np.maximum(self.sigma_hat_ ** 2 + epistemic_var, 1e-8))
        return mean, std


def _select_best_huber_alpha(
    X: np.ndarray,
    y: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
    alphas: tuple[float, ...] = (1.0, 3.0, 10.0, 30.0, 100.0),
) -> float:
    best_alpha = alphas[2]
    best_mse = float("inf")
    for alpha in alphas:
        sq_errs: list[float] = []
        for tr, va in splits:
            if len(tr) < 3:
                continue
            pipe = Pipeline([("scaler", StandardScaler()), ("model", HuberRegressor(alpha=alpha, epsilon=1.35, max_iter=1000))])
            pipe.fit(X[tr], y[tr])
            pred = pipe.predict(X[va])
            sq_errs.extend(((y[va] - pred) ** 2).tolist())
        if sq_errs:
            mse = float(np.mean(sq_errs))
            if mse < best_mse:
                best_mse = mse
                best_alpha = alpha
    return float(best_alpha)


@dataclass
class FerritinHybridPredictor:
    """Self-contained Ferritin-GeoHybrid predictor with optional anchor residual & uncertainty."""

    model_name: str
    feature_names: list[str]
    pipeline: Any
    residual_anchor_col: str | None = None
    alpha_selected: float | None = None
    w_prior_selected: float | None = None
    train_cohort: str = ""
    train_n: int = 0
    feature_stats: dict[str, dict[str, float]] | None = None

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        work = compute_derived_feature_frame(df)
        X = extract_feature_matrix(work, self.feature_names)
        raw = self.pipeline.predict(X).astype(np.float64)
        if self.residual_anchor_col is not None:
            anchor = work[self.residual_anchor_col].to_numpy(dtype=np.float64)
            return anchor + raw
        return raw

    def predict_with_std(self, df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray | None]:
        work = compute_derived_feature_frame(df)
        X = extract_feature_matrix(work, self.feature_names)
        if isinstance(self.pipeline, Pipeline):
            step = self.pipeline[-1]
            if isinstance(step, (BayesianRidge, NonNegHuberMAPRegressor)):
                Xs = self.pipeline[:-1].transform(X)
                mean, std = step.predict(Xs, return_std=True)
                if self.residual_anchor_col is not None:
                    mean = work[self.residual_anchor_col].to_numpy(dtype=np.float64) + mean
                return mean.astype(np.float64), std.astype(np.float64)
        return self.predict(work), None

    def linear_weights_summary(self) -> dict[str, Any] | None:
        if not isinstance(self.pipeline, Pipeline):
            return None
        est = self.pipeline[-1]
        if not hasattr(est, "coef_"):
            return None
        scaler = self.pipeline[0]
        coef_std = np.asarray(est.coef_, dtype=np.float64)
        coef_raw = coef_std / np.where(scaler.scale_ > 1e-12, scaler.scale_, 1.0)
        intercept = float(est.intercept_) - float(np.sum(coef_raw * scaler.mean_))
        return {
            "residual_anchor": self.residual_anchor_col,
            "standardized_coefs": {f: round(float(c), 4) for f, c in zip(self.feature_names, coef_std)},
            "raw_coefs": {f: round(float(c), 4) for f, c in zip(self.feature_names, coef_raw)},
            "raw_intercept": round(intercept, 4),
            "std_intercept_b0": round(float(est.intercept_), 4),
        }


def _select_nested_nonneg_huber(
    X: np.ndarray,
    y_fit: np.ndarray,
    y_true: np.ndarray,
    anchor_vals: np.ndarray | None,
    splits: list[tuple[np.ndarray, np.ndarray]],
    *,
    w_prior_grid: tuple[float, ...],
    alpha_grid: tuple[float, ...] = (5.0, 10.0, 20.0, 40.0, 80.0),
    delta: float = 8.0,
) -> tuple[float, float]:
    best_cfg = (float(alpha_grid[len(alpha_grid) // 2]), float(w_prior_grid[0]))
    best_mae = float("inf")
    for wp in w_prior_grid:
        for alpha in alpha_grid:
            abs_errs: list[float] = []
            for tr, va in splits:
                if len(tr) < 3:
                    continue
                pipe = Pipeline(
                    [
                        ("scaler", StandardScaler()),
                        ("model", NonNegHuberMAPRegressor(alpha=float(alpha), w_prior=float(wp), delta=delta)),
                    ]
                )
                pipe.fit(X[tr], y_fit[tr])
                pred = pipe.predict(X[va])
                if anchor_vals is not None:
                    pred = pred + anchor_vals[va]
                abs_errs.extend(np.abs(y_true[va] - pred).tolist())
            if abs_errs:
                mae = float(np.mean(abs_errs))
                if mae < best_mae:
                    best_mae = mae
                    best_cfg = (float(alpha), float(wp))
    return best_cfg


def fit_ferritin_hybrid(
    df_train: pd.DataFrame,
    *,
    model_family: str = "PriorAnchoredNonNegHuber",
    feature_names: list[str] | None = None,
    residual_anchor: bool = False,
    target_col: str = "delta_tm",
    group_col: str | None = "protein_group",
    train_cohort_label: str = "S1_All30",
    w_prior_grid: tuple[float, ...] = (0.0, 0.5, 1.0),
) -> FerritinHybridPredictor:
    """Fit a regularized Ferritin-GeoHybrid model strictly on df_train."""
    work = compute_derived_feature_frame(df_train)
    y = pd.to_numeric(work[target_col], errors="raise").to_numpy(dtype=np.float64)
    if residual_anchor:
        feats = list(feature_names or RESIDUAL_HYBRID_FEATURES_4)
        anchor_col = "geo_asm_k32"
        anchor_vals = work[anchor_col].to_numpy(dtype=np.float64)
        y_fit = y - anchor_vals
    else:
        feats = list(feature_names or HYBRID_FEATURES_5)
        anchor_col = None
        anchor_vals = None
        y_fit = y

    X = extract_feature_matrix(work, feats)
    groups = work[group_col].to_numpy() if group_col and group_col in work.columns else None
    splits = _inner_cv_for_train(groups, len(work))
    alpha_sel: float | None = None
    wp_sel: float | None = None

    if model_family == "PriorAnchoredNonNegHuber":
        alpha_sel, wp_sel = _select_nested_nonneg_huber(
            X, y_fit, y, anchor_vals, splits, w_prior_grid=w_prior_grid
        )
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", NonNegHuberMAPRegressor(alpha=alpha_sel, w_prior=wp_sel, delta=8.0)),
            ]
        )
        pipe.fit(X, y_fit)
    elif model_family == "ZeroPriorNonNegHuber":
        alpha_sel, wp_sel = _select_nested_nonneg_huber(
            X, y_fit, y, anchor_vals, splits, w_prior_grid=(0.0,)
        )
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", NonNegHuberMAPRegressor(alpha=alpha_sel, w_prior=0.0, delta=8.0)),
            ]
        )
        pipe.fit(X, y_fit)
    elif model_family == "BayesianRidge":
        pipe = Pipeline([("scaler", StandardScaler()), ("model", BayesianRidge())])
        pipe.fit(X, y_fit)
    elif model_family == "RidgeCV":
        alphas = np.array([0.5, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0], dtype=float)
        best_alpha = float(alphas[2])
        best_mse = float("inf")
        for a in alphas:
            sq = []
            for tr, va in splits:
                if len(tr) < 3:
                    continue
                p_in = Pipeline([("scaler", StandardScaler()), ("model", Ridge(alpha=float(a)))])
                p_in.fit(X[tr], y_fit[tr])
                sq.extend(((y_fit[va] - p_in.predict(X[va])) ** 2).tolist())
            if sq and float(np.mean(sq)) < best_mse:
                best_mse = float(np.mean(sq))
                best_alpha = float(a)
        alpha_sel = best_alpha
        pipe = Pipeline([("scaler", StandardScaler()), ("model", Ridge(alpha=best_alpha))])
        pipe.fit(X, y_fit)
    elif model_family == "Huber":
        alpha_sel = _select_best_huber_alpha(X, y_fit, splits)
        pipe = Pipeline(
            [("scaler", StandardScaler()), ("model", HuberRegressor(alpha=alpha_sel, epsilon=1.35, max_iter=1000))]
        )
        pipe.fit(X, y_fit)
    elif model_family == "ShallowGBDT":
        pipe = HistGradientBoostingRegressor(
            max_depth=2,
            max_iter=25,
            min_samples_leaf=4,
            l2_regularization=10.0,
            learning_rate=0.08,
            random_state=20261115,
        )
        pipe.fit(X, y_fit)
    elif model_family == "RegularizedMLP":
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    MLPRegressor(
                        hidden_layer_sizes=(8,),
                        alpha=10.0,
                        max_iter=3000,
                        early_stopping=False,
                        random_state=20261115,
                    ),
                ),
            ]
        )
        pipe.fit(X, y_fit)
    else:
        raise ValueError(f"Unsupported model_family={model_family!r}")

    fstats = {
        f: {
            "mean": float(np.mean(X[:, i])),
            "std": float(np.std(X[:, i])),
            "min": float(np.min(X[:, i])),
            "max": float(np.max(X[:, i])),
        }
        for i, f in enumerate(feats)
    }
    return FerritinHybridPredictor(
        model_name=f"{'Residual_' if residual_anchor else ''}{model_family}",
        feature_names=feats,
        pipeline=pipe,
        residual_anchor_col=anchor_col,
        alpha_selected=alpha_sel,
        w_prior_selected=wp_sel,
        train_cohort=train_cohort_label,
        train_n=len(work),
        feature_stats=fstats,
    )


def save_predictor(predictor: FerritinHybridPredictor, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(predictor, path)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    meta = {
        "model_name": predictor.model_name,
        "feature_names": predictor.feature_names,
        "residual_anchor_col": predictor.residual_anchor_col,
        "alpha_selected": predictor.alpha_selected,
        "w_prior_selected": predictor.w_prior_selected,
        "train_cohort": predictor.train_cohort,
        "train_n": predictor.train_n,
        "linear_weights": predictor.linear_weights_summary(),
        "artifact_sha256": sha,
    }
    path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return sha


def load_predictor(path: Path) -> FerritinHybridPredictor:
    obj = joblib.load(path)
    if not isinstance(obj, FerritinHybridPredictor):
        raise TypeError(f"Expected FerritinHybridPredictor, got {type(obj)}")
    return obj
