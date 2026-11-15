"""Unittest suite for Ferritin-GeoHybrid integrity, invariance, and zero-leakage (SKLEARN_PYTHON)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = RESEARCH_ROOT / "src"
for p in (str(SRC_DIR),):
    if p not in sys.path:
        sys.path.insert(0, p)

from hybrid_model import (  # noqa: E402
    HYBRID_FEATURES_5,
    compute_derived_feature_frame,
    fit_ferritin_hybrid,
    load_predictor,
    save_predictor,
)
from metrics import regression_metrics, roc_auc_binary  # noqa: E402
from paths import MODELS_DIR, RESULTS_DIR  # noqa: E402


class TestFerritinGeoHybridModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.df = compute_derived_feature_frame(pd.read_csv(RESULTS_DIR / "ferritin_features_all.csv"))
        cls.s1 = cls.df[cls.df["cohort"] == "S1_SinglePoint_24mer"].reset_index(drop=True)
        cls.s2 = cls.df[cls.df["cohort"] == "S2_EcDps_12mer"].reset_index(drop=True)
        cls.s5 = cls.df[cls.df["cohort"] == "S5_TierBC_Directional"].reset_index(drop=True)
        cls.model_path = MODELS_DIR / "ferritin_geohybrid.joblib"
        cls.meta_path = MODELS_DIR / "ferritin_geohybrid.json"

    def test_artifact_sha_and_reload_determinism(self) -> None:
        self.assertTrue(self.model_path.exists())
        self.assertTrue(self.meta_path.exists())
        meta = json.loads(self.meta_path.read_text(encoding="utf-8"))
        sha = hashlib.sha256(self.model_path.read_bytes()).hexdigest()
        self.assertEqual(sha, meta["artifact_sha256"])

        model = load_predictor(self.model_path)
        self.assertEqual(model.feature_names, HYBRID_FEATURES_5)
        p1, s1 = model.predict_with_std(self.df)
        p2 = model.predict(self.df)
        self.assertEqual(p1.shape, (75,))
        self.assertIsNotNone(s1)
        self.assertTrue(np.allclose(p1, p2, atol=1e-12))
        self.assertTrue(np.all(np.isfinite(p1)))
        self.assertTrue(np.all(np.isfinite(s1)) and np.all(s1 > 0.0))

    def test_thermodynamic_nonnegativity_and_monotonicity(self) -> None:
        model = load_predictor(self.model_path)
        weights = model.linear_weights_summary()["standardized_coefs"]
        for feat, val in weights.items():
            self.assertGreaterEqual(val, 0.0, msg=f"Negative weight on {feat}")
        self.assertTrue(all(v > 0.5 for v in weights.values()))

        probe = self.s2.iloc[:2].copy()
        base_pred = model.predict(probe)
        bumped = probe.copy()
        bumped["neg_fx_ac_ddg_bind_c3"] = bumped["neg_fx_ac_ddg_bind_c3"] + 2.0
        bumped["neg_fx_cage_ddg_eff"] = bumped["neg_fx_cage_ddg_eff"] + 2.0
        bump_pred = model.predict(bumped)
        self.assertTrue(np.all(bump_pred > base_pred))

    def test_zero_leakage_oof_invariant_to_heldout_labels(self) -> None:
        tr = self.s1[self.s1["protein_group"] != "O29424"].reset_index(drop=True)
        te_clean = self.s1[self.s1["protein_group"] == "O29424"].reset_index(drop=True)
        te_corrupt = te_clean.copy()
        te_corrupt["delta_tm"] = 9999.0

        m = fit_ferritin_hybrid(tr, model_family="PriorAnchoredNonNegHuber", feature_names=HYBRID_FEATURES_5)
        p_clean, s_clean = m.predict_with_std(te_clean)
        p_corrupt, s_corrupt = m.predict_with_std(te_corrupt)
        self.assertTrue(np.allclose(p_clean, p_corrupt, atol=1e-12))
        self.assertTrue(np.allclose(s_clean, s_corrupt, atol=1e-12))

    def test_roundtrip_save_in_tempdir(self) -> None:
        import tempfile

        model = load_predictor(self.model_path)
        with tempfile.TemporaryDirectory() as td:
            tmp_job = Path(td) / "copy.joblib"
            sha = save_predictor(model, tmp_job)
            reloaded = load_predictor(tmp_job)
            self.assertEqual(sha, hashlib.sha256(tmp_job.read_bytes()).hexdigest())
            self.assertTrue(np.allclose(model.predict(self.s2), reloaded.predict(self.s2), atol=1e-12))

    def test_predictions_csv_and_external_transfer_gates(self) -> None:
        pred_df = pd.read_csv(RESULTS_DIR / "hybrid_predictions_all.csv")
        self.assertEqual(len(pred_df), 75)
        self.assertTrue(np.all(np.isfinite(pred_df["hybrid_deployed_s1all30"])))

        model = load_predictor(self.model_path)
        m2 = regression_metrics(self.s2["delta_tm"].to_numpy(float), model.predict(self.s2))
        self.assertLess(m2["mae_c"], 6.0)
        self.assertGreater(m2["pearson_r"], 0.4)

        auc5 = roc_auc_binary(self.s5["binary_direction"].to_numpy(float).astype(int), model.predict(self.s5))
        self.assertIsNotNone(auc5)
        self.assertGreater(auc5, 0.90)


if __name__ == "__main__":
    unittest.main()
