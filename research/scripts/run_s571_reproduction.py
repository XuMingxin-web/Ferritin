#!/usr/bin/env python3
"""Verify GeoDTm-3D reproduction: full S571 cache audit + 36-case stratified local FoldX+ESM2+GeoDTm run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
import urllib.request
from pathlib import Path

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
if str(RESEARCH_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(RESEARCH_ROOT / "src"))

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import numpy as np

from metrics import paired_comparison, regression_metrics
from paths import (
    EXTERNAL_GEODTM_ROOT,
    FOLDX_BIN,
    GEODTM_PYTHON,
    GEOSTAB_ROOT,
    PINNED_SHA256,
    RESULTS_DIR,
    RUN_GEODTM_SCRIPT,
    S4346_CSV,
    S571_CSV,
    DATA_DIR,
)

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
ANCHOR_NAME = "rcsb_1CTS_A_H274R_7"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_substitution(wt_seq: str, mut_seq: str) -> tuple[int, str, str]:
    """Return 1-based mutation position, WT residue, and mutant residue."""
    if len(wt_seq) != len(mut_seq):
        raise ValueError("WT and mutant sequences differ in length")
    diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_seq, mut_seq)) if a != b]
    if len(diffs) != 1:
        raise ValueError(f"Expected 1 substitution, found {len(diffs)}")
    return diffs[0]


def download_pdb(pdb_id: str, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = dest_dir / f"{pdb_id.upper()}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 500:
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    urllib.request.urlretrieve(url, pdb_path)
    return pdb_path


def pdb_chain_aligns_cleanly(pdb_path: Path, chain_id: str, wt_seq: str, pos1: int, wt_res: str) -> bool:
    """Quick pre-check that PDB parses cleanly, chain exists, and mutated position is resolved."""
    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(pdb_path.stem, str(pdb_path))
        model = struct[0]
        if chain_id not in model:
            return False
        chain = model[chain_id]
        res_list = [r for r in chain if r.id[0] == " "]
        if len(res_list) < 20:
            return False
        if any(not all(atom in r for atom in ("N", "CA", "C", "O")) for r in res_list):
            return False
        pdb_seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in res_list)
        if any(c not in AA_SET for c in pdb_seq):
            return False
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
        aln = aligner.align(wt_seq, pdb_seq)[0]
        # Verify the target 0-based position in wt_seq aligns without gap to identical residue
        q_blocks, t_blocks = aln.aligned
        target0 = pos1 - 1
        for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
            if qs <= target0 < qe:
                pdb_idx = ts + (target0 - qs)
                return pdb_seq[pdb_idx] == wt_res
        return False
    except Exception:
        return False


def compute_author_s571_summary(rows: list[dict[str, str]]) -> dict:
    y = [float(r["dTm"]) for r in rows]
    summary = {
        "n_rows": len(rows),
        "unique_wt_sequences": len({r["wt_seq"] for r in rows}),
        "sha256_s571": sha256_file(S571_CSV),
        "sha256_s4346": sha256_file(S4346_CSV),
        "sha256_verified": (
            sha256_file(S571_CSV) == PINNED_SHA256["S571"]
            and sha256_file(S4346_CSV) == PINNED_SHA256["S4346"]
        ),
        "models": {},
    }
    if not summary["sha256_verified"]:
        raise RuntimeError("Pinned S571/S4346 SHA256 mismatch!")
    summary["models"]["Zero_baseline"] = regression_metrics(y, [0.0] * len(y))
    for col in ("GeoDTm-3D_dir", "GeoDTm-Seq_dir", "AUTO-MUTE_dir", "HoTMuSiC_dir"):
        valid = [(float(r["dTm"]), float(r[col])) for r in rows if r[col] not in ("", "NA", "nan")]
        summary["models"][col] = regression_metrics([a for a, _ in valid], [b for _, b in valid])
    return summary


def select_stratified_candidates(
    rows: list[dict[str, str]],
    pdb_dir: Path,
    cases_dir: Path,
    *,
    per_stratum: int = 12,
    seed: int = 20261115,
) -> list[dict]:
    """Select 3 * per_stratum verified cases balanced across dTm > 2, [-2, 2], < -2."""
    rng = random.Random(seed)
    strata = {"pos": [], "neu": [], "neg": []}
    anchor_row = None
    for r in rows:
        if r["name"] == ANCHOR_NAME:
            anchor_row = r
        val = float(r["dTm"])
        if val > 2.0:
            strata["pos"].append(r)
        elif val < -2.0:
            strata["neg"].append(r)
        else:
            strata["neu"].append(r)

    selected_rows: dict[str, list[dict[str, str]]] = {"pos": [], "neu": [], "neg": []}
    seen_names: set[str] = set()
    prot_counts: dict[str, dict[tuple[str, str], int]] = {"pos": {}, "neu": {}, "neg": {}}

    def try_row(r: dict[str, str], stratum: str) -> dict | None:
        parts = r["name"].split("_")
        if len(parts) < 5:
            return None
        pdb_id, chain_id = parts[1].upper(), parts[2]
        ph = float(parts[-1])
        pos1, wt_res, mut_res = parse_substitution(r["wt_seq"], r["mut_seq"])
        try:
            pdb_path = download_pdb(pdb_id, pdb_dir)
        except Exception:
            return None
        if not pdb_chain_aligns_cleanly(pdb_path, chain_id, r["wt_seq"], pos1, wt_res):
            return None
        return {
            "case_id": r["name"],
            "pdb_path": str(pdb_path.resolve()),
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "model_index0": 0,
            "normalize_mse": True,
            "reference_sequence": r["wt_seq"],
            "mutation_position1": pos1,
            "wt_residue": wt_res,
            "mutant_residue": mut_res,
            "pH": ph,
            "pH_status": "recorded_assay_pH",
            "role": "ferritin_external_evaluation",
            "feature_case_directory": str((cases_dir / r["name"]).resolve()),
            "stratum": stratum,
            "observed_dtm_c": float(r["dTm"]),
            "author_geodtm3d_c": float(r["GeoDTm-3D_dir"]),
            "author_geodtmseq_c": float(r["GeoDTm-Seq_dir"]),
        }

    cases: list[dict] = []
    if anchor_row is not None:
        c = try_row(anchor_row, "neu")
        if c is None:
            raise RuntimeError("Anchor 1CTS case failed structural validation")
        cases.append(c)
        selected_rows["neu"].append(anchor_row)
        seen_names.add(anchor_row["name"])
        prot_counts["neu"][("1CTS", "A")] = 1

    for stratum in ("pos", "neu", "neg"):
        pool = list(strata[stratum])
        rng.shuffle(pool)
        # First pass max 1 per protein; second pass max 2 per protein to keep diversity high
        for max_per_prot in (1, 2, 4):
            for r in pool:
                if len(selected_rows[stratum]) >= per_stratum:
                    break
                if r["name"] in seen_names:
                    continue
                parts = r["name"].split("_")
                key = (parts[1].upper(), parts[2])
                if prot_counts[stratum].get(key, 0) >= max_per_prot:
                    continue
                case = try_row(r, stratum)
                if case is not None:
                    cases.append(case)
                    selected_rows[stratum].append(r)
                    seen_names.add(r["name"])
                    prot_counts[stratum][key] = prot_counts[stratum].get(key, 0) + 1
        if len(selected_rows[stratum]) < per_stratum:
            raise RuntimeError(f"Only found {len(selected_rows[stratum])} valid candidates for {stratum}")

    return cases


def run_geodtm_stages(cases_json: Path, out_dir: Path, threads: int) -> None:
    for stage in ("prepare", "embeddings", "predict"):
        cmd = [
            str(GEODTM_PYTHON),
            "-B",
            str(RUN_GEODTM_SCRIPT),
            "--cases",
            str(cases_json),
            "--stage",
            stage,
            "--foldx",
            str(FOLDX_BIN),
            "--out",
            str(out_dir),
            "--threads",
            str(threads),
        ]
        print(f"Running stage={stage} on {cases_json.name} ...", flush=True)
        res = subprocess.run(cmd, cwd=str(EXTERNAL_GEODTM_ROOT), text=True)
        if res.returncode != 0:
            raise RuntimeError(f"run_geodtm.py failed at stage={stage} (exit {res.returncode})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-stratum", type=int, default=12)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-inference", action="store_true")
    parser.add_argument("--force-reselect", action="store_true")
    args = parser.parse_args()

    out_root = RESULTS_DIR / "s571_repro"
    out_root.mkdir(parents=True, exist_ok=True)
    cases_dir = out_root / "cases"
    pdb_dir = DATA_DIR / "s571_pdb"

    with S571_CSV.open(newline="", encoding="utf-8") as f:
        s571_rows = list(csv.DictReader(f))

    author_summary = compute_author_s571_summary(s571_rows)
    cases_json = out_root / f"s571_stratified_{args.per_stratum * 3}_cases.json"
    if cases_json.exists() and not args.force_reselect:
        cases = json.loads(cases_json.read_text())
    else:
        cases = select_stratified_candidates(s571_rows, pdb_dir, cases_dir, per_stratum=args.per_stratum)
        cases_json.write_text(json.dumps(cases, indent=2) + "\n")

    # Strip metadata fields when passing to run_geodtm if strictly checked, or pass core schema
    core_keys = {
        "case_id", "pdb_path", "pdb_id", "chain_id", "model_index0", "normalize_mse",
        "reference_sequence", "mutation_position1", "wt_residue", "mutant_residue",
        "pH", "pH_status", "role", "feature_case_directory",
    }
    runner_cases = [{k: c[k] for k in core_keys} for c in cases]
    runner_cases_json = out_root / "runner_cases.json"
    runner_cases_json.write_text(json.dumps(runner_cases, indent=2) + "\n")

    if not args.skip_inference:
        run_geodtm_stages(runner_cases_json, cases_dir, args.threads)

    # Run upstream official predict.py CLI on the 1CTS anchor ensemble to confirm tensor equivalence
    anchor_ens = cases_dir / ANCHOR_NAME / "ensemble.pt"
    cli_dir = out_root / "official_cli_1CTS"
    cli_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(GEODTM_PYTHON),
            "-B",
            str(GEOSTAB_ROOT / "model_dTm_3D/predict.py"),
            "--ensemble_feature",
            str(anchor_ens),
            "--saved_folder",
            str(cli_dir),
        ],
        cwd=str(EXTERNAL_GEODTM_ROOT),
        check=True,
    )
    anchor_cli_pred = float((cli_dir / "result_dTm.txt").read_text().strip())

    out_rows = []
    for c in cases:
        cid = c["case_id"]
        local_pred = float((cases_dir / cid / "result_dTm.txt").read_text().strip())
        out_rows.append(
            {
                "case_id": cid,
                "pdb_id": c["pdb_id"],
                "chain_id": c["chain_id"],
                "stratum": c["stratum"],
                "mutation": f"{c['wt_residue']}{c['mutation_position1']}{c['mutant_residue']}",
                "pH": c["pH"],
                "observed_dtm_c": c["observed_dtm_c"],
                "author_geodtm3d_c": c["author_geodtm3d_c"],
                "local_geodtm3d_c": local_pred,
                "local_minus_author_c": local_pred - c["author_geodtm3d_c"],
                "abs_diff_local_vs_author_c": abs(local_pred - c["author_geodtm3d_c"]),
                "abs_err_local_c": abs(local_pred - c["observed_dtm_c"]),
                "abs_err_author_c": abs(c["author_geodtm3d_c"] - c["observed_dtm_c"]),
            }
        )

    csv_path = out_root / f"s571_local_vs_author_{len(out_rows)}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    ys = [r["observed_dtm_c"] for r in out_rows]
    authors = [r["author_geodtm3d_c"] for r in out_rows]
    locals_ = [r["local_geodtm3d_c"] for r in out_rows]
    anchor_row = next(r for r in out_rows if r["case_id"] == ANCHOR_NAME)

    diffs = np.abs(np.asarray(locals_) - np.asarray(authors))
    summary = {
        "full_s571_author_cached": author_summary,
        "stratified_sample_size": len(out_rows),
        "unique_proteins_in_sample": len({(r["pdb_id"], r["chain_id"]) for r in out_rows}),
        "anchor_1CTS_H274R": {
            "observed_dtm_c": anchor_row["observed_dtm_c"],
            "author_cached_c": anchor_row["author_geodtm3d_c"],
            "local_runner_c": anchor_row["local_geodtm3d_c"],
            "official_cli_same_tensor_c": anchor_cli_pred,
            "abs_diff_runner_vs_cli_c": abs(anchor_row["local_geodtm3d_c"] - anchor_cli_pred),
            "abs_diff_local_vs_author_c": anchor_row["abs_diff_local_vs_author_c"],
        },
        "local_vs_author_agreement": {
            "pearson_r": float(np.corrcoef(authors, locals_)[0, 1]),
            "mae_c": float(np.mean(diffs)),
            "median_abs_diff_c": float(np.median(diffs)),
            "max_abs_diff_c": float(np.max(diffs)),
            "same_sign_rate": float(np.mean((np.asarray(authors) > 0) == (np.asarray(locals_) > 0))),
        },
        "metrics_on_sample_vs_experiment": {
            "Zero_baseline": regression_metrics(ys, [0.0] * len(ys)),
            "Author_cached_GeoDTm3D": regression_metrics(ys, authors),
            "Local_FoldX51_GeoDTm3D": regression_metrics(ys, locals_),
            "Paired_Local_vs_Author": paired_comparison(ys, locals_, authors),
        },
    }
    summary_path = out_root / "reproduction_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
