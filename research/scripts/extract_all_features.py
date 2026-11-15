#!/usr/bin/env python3
"""Extract all FoldX, ESM-2, Monomer/Assembly-GeoDTm-3D & S4346 Ridge features across S1-S5."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import warnings

from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = RESEARCH_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from foldx_assembly import run_foldx_monomer_and_subcomplex  # noqa: E402
from geodtm_assembly_adapter import (  # noqa: E402
    _build_node_pool,
    _parse_canon_muts,
    _slice_geodtm_input,
    compute_geodtm_features_for_case,
    get_esm2_sequence_features,
    verify_1cts_equivalence,
)
from metrics import paired_comparison, regression_metrics, roc_auc_binary  # noqa: E402
from paths import DATA_DIR, EXTERNAL_GEODTM_ROOT, GEOSTAB_ROOT, RESULTS_DIR, SKLEARN_PYTHON  # noqa: E402

import torch  # noqa: E402


def verify_ferritin_monomer_cli_equivalence(wt_templates: dict) -> dict[str, float]:
    """Verify EcBfr R30A monomer tensor matches official GeoStab predict.py CLI."""
    tkey = "EcBfr_1BFR_A"
    fx = run_foldx_monomer_and_subcomplex(tkey, "R30A", template_spec=wt_templates[tkey])
    wt_seq = wt_templates[tkey]["pdb_seq"]
    pos0 = 29  # R30A (1-based 30 -> 0-based 29)
    mut_seq = wt_seq[:pos0] + "A" + wt_seq[pos0 + 1 :]
    ph = 7.2

    wt_esm = get_esm2_sequence_features(wt_seq)["repr33"]
    mt_esm = get_esm2_sequence_features(mut_seq)["repr33"]
    wt_pool = _build_node_pool(
        Path(fx["mono_repair_pdb"]),
        target_chain_id="A",
        full_seq=wt_seq,
        esm_repr=wt_esm,
        mut_seq_pos0=[pos0],
    )
    mt_pool = _build_node_pool(
        Path(fx["mono_mutant_pdb"]),
        target_chain_id="A",
        full_seq=mut_seq,
        esm_repr=mt_esm,
        mut_seq_pos0=[pos0],
        wt_seq=wt_seq,
        wt_esm_repr=wt_esm,
    )
    wt_inp, _ = _slice_geodtm_input(wt_pool, k=32, ph=ph)
    mt_inp, _ = _slice_geodtm_input(mt_pool, k=32, ph=ph)

    with tempfile.TemporaryDirectory(prefix="cli_ecbfr_") as tmp_str:
        tmp = Path(tmp_str)
        ens_path = tmp / "ensemble.pt"
        torch.save({"dTm": {"wt": wt_inp, "mut": mt_inp}}, ens_path)
        out_dir = tmp / "cli_out"
        out_dir.mkdir()
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(GEOSTAB_ROOT / "model_dTm_3D/predict.py"),
                "--ensemble_feature",
                str(ens_path),
                "--saved_folder",
                str(out_dir),
            ],
            cwd=str(EXTERNAL_GEODTM_ROOT),
            check=True,
        )
        cli_pred = float((out_dir / "result_dTm.txt").read_text(encoding="utf-8").strip())

    # Compare with adapter prediction
    feats = compute_geodtm_features_for_case(
        case_key=fx["case_key"],
        wt_seq=wt_seq,
        mut_seq=mut_seq,
        mutation_canon="R30A",
        ph=ph,
        mono_wt_pdb=Path(fx["mono_repair_pdb"]),
        mono_mt_pdb=Path(fx["mono_mutant_pdb"]),
        sub_wt_pdb=Path(fx["sub_repair_pdb"]),
        sub_mt_pdb=Path(fx["sub_mutant_pdb"]),
    )
    return {
        "official_cli_ecbfr_r30a_c": cli_pred,
        "adapter_ecbfr_r30a_c": feats["geo_mono_k32"],
        "abs_diff_c": abs(cli_pred - feats["geo_mono_k32"]),
    }


def score_frozen_s4346_ridge(rows: list[dict]) -> dict[str, dict[str, float | None]]:
    """Call SKLEARN_PYTHON to score single-point variants with archived S4346 Ridge baselines."""
    single_payload = []
    for r in rows:
        if int(float(r["n_mutations"])) != 1:
            continue
        single_payload.append(
            {
                "case_id": r["case_id"],
                "wt_seq": r["wt_sequence"],
                "mut_seq": r["mut_sequence"],
            }
        )
    helper_code = f"""
import json, sys, hashlib, joblib, pandas as pd
from pathlib import Path
ROOT = Path({str(EXTERNAL_GEODTM_ROOT)!r})
sys.path.insert(0, str(ROOT))
from scripts.model_probe_ridge import feature_values, FEATURE_SETS
payload = json.loads(sys.stdin.read())
models = {{}}
for name, sha in [
    ("generic_ridge_substitution", "0419e803040116eda0d27b61b5a9b775226c59b16395ade38a69f07f08fb9803"),
    ("generic_ridge_sequence_context", "7aa8eccef55e5511ed58ebfbfb62d545313d0116d9b60a4655f1fc68f32979d2"),
]:
    mp = ROOT / "models" / f"{{name}}.joblib"
    assert hashlib.sha256(mp.read_bytes()).hexdigest() == sha, f"SHA mismatch {{name}}"
    models[name] = joblib.load(mp)
feats = pd.DataFrame([feature_values(x["wt_seq"], x["mut_seq"]) for x in payload])
out = {{}}
p_sub = models["generic_ridge_substitution"]["pipeline"].predict(feats[FEATURE_SETS["generic_ridge_substitution"]])
p_ctx = models["generic_ridge_sequence_context"]["pipeline"].predict(feats[FEATURE_SETS["generic_ridge_sequence_context"]])
for x, a, b in zip(payload, p_sub, p_ctx):
    out[x["case_id"]] = {{"s4346_ridge_sub": float(a), "s4346_ridge_seqctx": float(b)}}
print(json.dumps(out))
"""
    proc = subprocess.run(
        [str(SKLEARN_PYTHON), "-B", "-c", helper_code],
        input=json.dumps(single_payload),
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(proc.stdout)


def _float_or_none(v: Any) -> float | None:
    if v in ("", None, "NA", "nan", "None"):
        return None
    f = float(v)
    return f if math.isfinite(f) else None


def evaluate_baselines(all_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Compute baseline metrics across S1-S5 cohorts."""
    pred_cols = [
        ("Zero_baseline", None),
        ("Monomer_GeoDTm3D_k32", "geo_mono_k32"),
        ("Monomer_GeoDTm3D_k48", "geo_mono_k48"),
        ("Assembly_GeoDTm3D_k32", "geo_asm_k32"),
        ("Assembly_GeoDTm3D_k48", "geo_asm_k48"),
        ("FoldX_mono_neg_ddG", "neg_fx_mono_ddg"),
        ("FoldX_sub_neg_ddG_per_chain", "neg_fx_sub_ddg_per_chain"),
        ("FoldX_AC_neg_bind_ddG", "neg_fx_ac_ddg_bind_total"),
        ("FoldX_cage_neg_eff_ddG", "neg_fx_cage_ddg_eff"),
        ("ESM2_WT_marginal_LLR", "esm2_wt_marginal_llr"),
        ("S4346_Frozen_Ridge_Sub", "s4346_ridge_sub"),
        ("S4346_Frozen_Ridge_SeqCtx", "s4346_ridge_seqctx"),
    ]

    def eval_subset(sub_rows: list[dict[str, Any]], label_key: str = "delta_tm") -> dict[str, Any]:
        valid_rows = [r for r in sub_rows if _float_or_none(r.get(label_key)) is not None]
        y = [float(r[label_key]) for r in valid_rows]
        res: dict[str, Any] = {"n": len(valid_rows), "models": {}}
        mono_preds = [float(r["geo_mono_k32"]) for r in valid_rows]
        asm32_preds = [float(r["geo_asm_k32"]) for r in valid_rows]
        for name, col in pred_cols:
            if col is None:
                p = [0.0] * len(y)
            else:
                pairs = [(float(r[label_key]), _float_or_none(r.get(col))) for r in valid_rows]
                if any(b is None for _, b in pairs):
                    pairs_valid = [(a, b) for a, b in pairs if b is not None]
                    if len(pairs_valid) < 3:
                        continue
                    res["models"][name] = regression_metrics([a for a, _ in pairs_valid], [b for _, b in pairs_valid])
                    continue
                p = [b for _, b in pairs]
            m = regression_metrics(y, p)
            if name in ("Assembly_GeoDTm3D_k32", "Assembly_GeoDTm3D_k48", "FoldX_AC_neg_bind_ddG", "FoldX_cage_neg_eff_ddG", "ESM2_WT_marginal_LLR"):
                groups = [r["protein_group"] for r in valid_rows]
                m["paired_vs_mono_k32"] = paired_comparison(y, p, mono_preds, groups=groups)
            if name == "Assembly_GeoDTm3D_k48":
                groups = [r["protein_group"] for r in valid_rows]
                m["paired_vs_asm_k32"] = paired_comparison(y, p, asm32_preds, groups=groups)
            res["models"][name] = m
        return res

    s1_core = [r for r in all_rows if r["subcohort"] == "S1_PairedCore23"]
    s1_all = [r for r in all_rows if r["cohort"] == "S1_SinglePoint_24mer"]
    s1_conf = [dict(r) for r in s1_core]
    for r in s1_conf:
        if _float_or_none(r.get("delta_tm_alt_thesis")) is not None:
            r["delta_tm"] = r["delta_tm_alt_thesis"]
    s2 = [r for r in all_rows if r["cohort"] == "S2_EcDps_12mer"]
    s3_exact = [r for r in all_rows if r["subcohort"] == "S3_Exact9"]
    s3_all = [r for r in all_rows if r["cohort"] == "S3_MultiPoint"]
    s4 = [r for r in all_rows if r["cohort"] == "S4_Secondary"]
    s5 = [r for r in all_rows if r["cohort"] == "S5_TierBC_Directional"]

    # Binary directional evaluation on S5 (7 < WT vs 6 > WT)
    s5_labels = [int(float(r["binary_direction"])) for r in s5]
    s5_res: dict[str, Any] = {"n": len(s5), "n_pos": sum(s5_labels), "n_neg": len(s5_labels) - sum(s5_labels), "models": {}}
    for name, col in pred_cols:
        if col is None:
            continue
        vals = [_float_or_none(r.get(col)) for r in s5]
        if any(v is None for v in vals):
            continue
        preds = [float(v) for v in vals]
        auc = roc_auc_binary(s5_labels, preds)
        acc_gt0 = float(np.mean([(1 if p > 0 else 0) == y for p, y in zip(preds, s5_labels)]))
        s5_res["models"][name] = {"roc_auc": auc, "sign_acc": acc_gt0}

    return {
        "S1_PairedCore23": eval_subset(s1_core),
        "S1_All30": eval_subset(s1_all),
        "S1_PairedCore23_AltThesis": eval_subset(s1_conf),
        "S2_EcDps_Single8": eval_subset(s2),
        "S3_Exact9": eval_subset(s3_exact),
        "S3_All16": eval_subset(s3_all),
        "S4_Secondary8": eval_subset(s4),
        "S5_Unseen13_Directional": s5_res,
    }


def main() -> int:
    t0 = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    wt_templates = json.loads((DATA_DIR / "wt_templates.json").read_text(encoding="utf-8"))

    eq_1cts = verify_1cts_equivalence()
    eq_ecbfr = verify_ferritin_monomer_cli_equivalence(wt_templates)
    print("Verification Gate 1CTS:", json.dumps(eq_1cts), flush=True)
    print("Verification Gate EcBfr R30A:", json.dumps(eq_ecbfr), flush=True)
    assert eq_1cts["abs_diff_c"] < 1e-4, eq_1cts
    assert eq_ecbfr["abs_diff_c"] < 1e-4, eq_ecbfr

    with (DATA_DIR / "curated_ferritin_cohorts.csv").open(newline="", encoding="utf-8") as f:
        curated_rows = list(csv.DictReader(f))

    ridge_scores = score_frozen_s4346_ridge(curated_rows)
    enriched_rows: list[dict[str, Any]] = []
    for idx, row in enumerate(curated_rows, 1):
        tkey = row["template_key"]
        fx = run_foldx_monomer_and_subcomplex(
            tkey, row["mutation_pdb"], template_spec=wt_templates[tkey]
        )
        geo = compute_geodtm_features_for_case(
            case_key=fx["case_key"],
            wt_seq=row["wt_sequence"],
            mut_seq=row["mut_sequence"],
            mutation_canon=row["mutation_canon"],
            ph=float(row["pH"]),
            mono_wt_pdb=Path(fx["mono_repair_pdb"]),
            mono_mt_pdb=Path(fx["mono_mutant_pdb"]),
            sub_wt_pdb=Path(fx["sub_repair_pdb"]),
            sub_mt_pdb=Path(fx["sub_mutant_pdb"]),
        )

        rec: dict[str, Any] = dict(row)
        # Merge FoldX features cleanly (skip non-scalar dicts/paths)
        for k, v in fx.items():
            if isinstance(v, (int, float)):
                rec[k] = v
        rec["neg_fx_mono_ddg"] = -float(fx["fx_mono_ddg"])
        rec["neg_fx_sub_ddg_per_chain"] = -float(fx["fx_sub_ddg_per_chain"])
        rec["neg_fx_ac_ddg_bind_total"] = -float(fx["fx_ac_ddg_bind_total"])
        rec["neg_fx_cage_ddg_eff"] = -float(fx["fx_cage_ddg_eff"])

        for k, v in geo.items():
            rec[k] = v

        rs = ridge_scores.get(row["case_id"], {})
        rec["s4346_ridge_sub"] = rs.get("s4346_ridge_sub", "")
        rec["s4346_ridge_seqctx"] = rs.get("s4346_ridge_seqctx", "")
        enriched_rows.append(rec)
        print(
            f"[{idx:2d}/{len(curated_rows)}] {row['case_id']:30s} obs={row['delta_tm'] or ('dir:'+row['binary_direction']):>6s} | "
            f"mono32={geo['geo_mono_k32']:+.2f} asm32={geo['geo_asm_k32']:+.2f} asm48={geo['geo_asm_k48']:+.2f} | "
            f"fx_m={fx['fx_mono_ddg']:+.2f} fx_b={fx['fx_ac_ddg_bind_total']:+.2f} esm={geo['esm2_wt_marginal_llr']:+.2f}",
            flush=True,
        )

    out_csv = RESULTS_DIR / "ferritin_features_all.csv"
    fieldnames = list(enriched_rows[0].keys())
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(enriched_rows)

    summary = {
        "verification_gates": {"1CTS_H274R": eq_1cts, "EcBfr_R30A": eq_ecbfr},
        "total_cases": len(enriched_rows),
        "elapsed_sec": round(time.perf_counter() - t0, 2),
        "baselines": evaluate_baselines(enriched_rows),
    }
    summary_path = RESULTS_DIR / "baseline_evaluation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {out_csv} & {summary_path} in {summary['elapsed_sec']}s", flush=True)
    for cname in ("S1_PairedCore23", "S1_All30", "S2_EcDps_Single8", "S3_Exact9", "S4_Secondary8"):
        c = summary["baselines"][cname]
        print(f"\n=== {cname} (n={c['n']}) ===")
        for mname, m in c["models"].items():
            r_val = m.get("pearson_r")
            rho_val = m.get("spearman_rho")
            r_str = f"{r_val:+.3f}" if r_val is not None else "   NA "
            rho_str = f"{rho_val:+.3f}" if rho_val is not None else "   NA "
            print(f"  {mname:30s}: r={r_str}  rho={rho_str}  MAE={m['mae_c']:.2f}  RMSE={m['rmse_c']:.2f}")
    s5 = summary["baselines"]["S5_Unseen13_Directional"]
    print(f"\n=== S5_Unseen13_Directional (n={s5['n']}, {s5['n_pos']}+ / {s5['n_neg']}-) ===")
    for mname, m in s5["models"].items():
        print(f"  {mname:30s}: AUC={m['roc_auc']:.3f}  sign_acc={m['sign_acc']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
