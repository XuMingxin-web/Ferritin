"""Pytest suite for metrics, curation, BIOMT assembly, and GeoDTm-3D parity (GEODTM_PYTHON)."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = RESEARCH_ROOT / "src"
for p in (str(RESEARCH_ROOT), str(SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from assembly_builder import build_biological_assembly_and_subcomplex  # noqa: E402
from geodtm_assembly_adapter import (  # noqa: E402
    _extract_chain_aa_seq,
    _split_pdb_by_chain_clean,
    verify_1cts_equivalence,
)
from metrics import paired_comparison, regression_metrics, roc_auc_binary  # noqa: E402
from paths import DATA_DIR, RESULTS_DIR  # noqa: E402


def test_metrics_contracts_and_bootstrap():
    y = np.array([1.0, -2.0, 3.0, -4.0, 5.0])
    p_good = np.array([1.1, -1.9, 2.8, -3.7, 4.9])
    p_bad = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    m = regression_metrics(y, p_good)
    assert m["n"] == 5
    assert m["mae_c"] < 0.3
    assert m["pearson_r"] > 0.99
    assert np.isclose(m["spearman_rho"], 1.0)
    assert np.isclose(m["sign_accuracy"], 1.0)

    cmp_res = paired_comparison(y, p_good, p_bad, groups=["A", "A", "B", "B", "C"], n_boot=300)
    assert cmp_res["mean_abs_err_diff_c"] < -2.0
    assert cmp_res["ci95_high_c"] < 0.0
    assert cmp_res["wins"] == 5 and cmp_res["losses"] == 0

    auc = roc_auc_binary([1, 1, 0, 0], [0.9, 0.4, 0.2, -0.5])
    assert auc == 1.0


def test_curated_cohorts_and_templates_integrity():
    templates = json.loads((DATA_DIR / "wt_templates.json").read_text(encoding="utf-8"))
    assert len(templates) == 14

    with (DATA_DIR / "curated_ferritin_cohorts.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 75

    counts: dict[str, int] = {}
    for r in rows:
        counts[r["cohort"]] = counts.get(r["cohort"], 0) + 1
        assert r["template_key"] in templates
        assert len(r["wt_sequence"]) == len(r["mut_sequence"])
        assert r["wt_sequence"] != r["mut_sequence"]
        assert r["mutation_pdb"]
    assert counts == {
        "S1_SinglePoint_24mer": 30,
        "S2_EcDps_12mer": 8,
        "S3_MultiPoint": 16,
        "S4_Secondary": 8,
        "S5_TierBC_Directional": 13,
    }


def test_assembly_builder_24mer_and_12mer_symmetry(tmp_path: Path):
    templates = json.loads((DATA_DIR / "wt_templates.json").read_text(encoding="utf-8"))

    # 24-mer EcBfr must yield 24 total subunits with C2, C3, and C4 symmetry partners
    t24 = templates["EcBfr_1BFR_A"]
    info24 = build_biological_assembly_and_subcomplex(
        Path(t24["pdb_path"]),
        target_chain_id=t24["target_chain_id"],
        biomt_source_pdb=Path(t24["biomt_source_pdb"]) if t24.get("biomt_source_pdb") else None,
        out_monomer_pdb=tmp_path / "bfr_mono.pdb",
        out_subcomplex_pdb=tmp_path / "bfr_sub.pdb",
    )
    assert info24["total_subunits"] == 24
    assert len(info24["subcomplex_chains"]) == 7
    sym24 = {n["symmetry_class"] for n in info24["neighbors"]}
    assert sym24 == {"C2", "C3", "C4"}
    assert t24["target_chain_id"] in info24["subcomplex_chains"]

    # 12-mer tetrahedral EcDps must yield 12 total subunits and ZERO C4 partners
    t12 = templates["EcDps_1F33_F"]
    info12 = build_biological_assembly_and_subcomplex(
        Path(t12["pdb_path"]),
        target_chain_id=t12["target_chain_id"],
        biomt_source_pdb=Path(t12["biomt_source_pdb"]) if t12.get("biomt_source_pdb") else None,
        out_monomer_pdb=tmp_path / "dps_mono.pdb",
        out_subcomplex_pdb=tmp_path / "dps_sub.pdb",
    )
    assert info12["total_subunits"] == 12
    sym12 = [n["symmetry_class"] for n in info12["neighbors"]]
    assert "C4" not in sym12
    assert "C3" in sym12


def test_geodtm_adapter_official_1cts_equivalence_and_histidine_alias(tmp_path: Path):
    eq = verify_1cts_equivalence()
    assert eq["abs_diff_c"] < 1e-4
    assert np.isfinite(eq["official_1cts_c"])

    # Verify FoldX HIS protonation aliases (H1S/H2S) & HETATM MSE normalize cleanly via _split_pdb_by_chain_clean
    fake_pdb = tmp_path / "his_alias.pdb"
    fake_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  CA  H1S A   1       0.000   0.000   0.000  1.00 20.00           C",
                "ATOM      2  CA  H2S A   2       3.800   0.000   0.000  1.00 20.00           C",
                "HETATM    3  CA  MSE A   3       7.600   0.000   0.000  1.00 20.00           C",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    split_dir = tmp_path / "split"
    split_dir.mkdir()
    split = _split_pdb_by_chain_clean(fake_pdb, split_dir)
    assert len(split) == 1
    assert _extract_chain_aa_seq(split[0][1]) == "HHM"


def test_feature_table_integrity():
    csv_path = RESULTS_DIR / "ferritin_features_all.csv"
    assert csv_path.exists()
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 75
    required = (
        "geo_mono_k32",
        "geo_asm_k32",
        "neg_fx_cage_ddg_eff",
        "fx_ac_ddg_bind_c3",
        "fx_ac_ddg_bind_solv_hyd",
        "fx_ac_ddg_bind_vdw",
        "geo_asm_enc_shift_k32",
    )
    for r in rows:
        for col in required:
            assert np.isfinite(float(r[col]))
