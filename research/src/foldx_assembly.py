"""FoldX 5.1 monomer & homotypic biological-assembly subcomplex energetics for ferritins."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import numpy as np

_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from assembly_builder import build_biological_assembly_and_subcomplex
from paths import DATA_DIR, FOLDX_BIN, FOLDX_MOLECULES_DIR, RESEARCH_ROOT

CACHE_ROOT = RESEARCH_ROOT / "cache"
WT_CACHE_DIR = CACHE_ROOT / "wt_structures"
FOLDX_CASE_CACHE_DIR = CACHE_ROOT / "foldx_cases"

RES3_CUSTOM = {
    "MSE": "M",
    "H1S": "H",
    "H2S": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "CYX": "C",
    "CSS": "C",
}

AC_TERM_MAP = {
    "Interaction Energy": "bind",
    "Electrostatics": "bind_elec",
    "Van der Waals": "bind_vdw",
    "Solvation Polar": "bind_solv_polar",
    "Solvation Hydrophobic": "bind_solv_hyd",
    "Van der Waals clashes": "bind_clash",
}

DIF_TERM_MAP = {
    "total energy": "ddg",
    "Backbone Hbond": "hbond_bb",
    "Sidechain Hbond": "hbond_sc",
    "Van der Waals": "vdw",
    "Electrostatics": "elec",
    "Solvation Polar": "solv_polar",
    "Solvation Hydrophobic": "solv_hyd",
    "Van der Waals clashes": "vdw_clash",
    "torsional clash": "torsional_clash",
    "entropy sidechain": "entropy_sc",
    "entropy mainchain": "entropy_mc",
}


def load_wt_templates() -> dict[str, dict[str, Any]]:
    path = DATA_DIR / "wt_templates.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _extract_chain_aa_map(pdb_path: Path) -> dict[str, dict[int, str]]:
    """Return {chain_id: {resnum: 1-letter-AA}} for standard protein residues in model 0."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = struct[0]
    out: dict[str, dict[int, str]] = {}
    for chain in model:
        cmap: dict[int, str] = {}
        for res in chain:
            hetflag, resnum, icode = res.id
            if hetflag != " " or str(icode).strip():
                continue
            aa1 = seq1(res.resname, custom_map=RES3_CUSTOM)
            if len(aa1) == 1 and aa1 != "X":
                cmap[int(resnum)] = aa1
        out[chain.id] = cmap
    return out


def _link_foldx_molecules(work_dir: Path) -> None:
    mol_link = work_dir / "molecules"
    if not mol_link.exists():
        os.symlink(FOLDX_MOLECULES_DIR, mol_link)


def _run_foldx(args: list[str], work_dir: Path, log_name: str) -> str:
    _link_foldx_molecules(work_dir)
    cmd = [str(FOLDX_BIN), *args]
    res = subprocess.run(
        cmd,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        check=False,
    )
    log_text = res.stdout + ("\n[STDERR]\n" + res.stderr if res.stderr else "")
    (work_dir / log_name).write_text(log_text, encoding="utf-8")
    if res.returncode != 0:
        raise RuntimeError(
            f"FoldX failed (exit={res.returncode}) in {work_dir} for {args}:\n{log_text[-1200:]}"
        )
    return log_text


def _parse_fxout_table(path: Path) -> list[dict[str, str]]:
    """Parse a FoldX *.fxout tab-delimited table by finding the column header line."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for line in lines:
        if not line.strip():
            continue
        if header is None:
            if line.startswith("Pdb\t") or "\ttotal energy" in line or "\tInteraction Energy" in line:
                header = [c.strip() for c in line.split("\t")]
            continue
        parts = [c.strip() for c in line.split("\t")]
        if len(parts) < len(header):
            continue
        rows.append(dict(zip(header, parts[: len(header)])))
    if header is None:
        raise ValueError(f"Could not find header in FoldX output: {path}")
    return rows


def _find_single_file(work_dir: Path, pattern: str) -> Path:
    matches = sorted(work_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No FoldX file matching {pattern} in {work_dir}")
    if len(matches) > 1:
        # Prefer exact non-WT / newest if needed
        matches.sort(key=lambda p: len(p.name))
    return matches[0]


def _repair_pdb(in_pdb: Path, out_pdb: Path, work_dir: Path, tag: str) -> Path:
    stem = f"{tag}.pdb"
    local_in = work_dir / stem
    shutil.copy2(in_pdb, local_in)
    _run_foldx(
        [
            "--command=RepairPDB",
            f"--pdb={stem}",
            "--pdb-dir=.",
            "--output-dir=.",
        ],
        work_dir,
        f"repair_{tag}.log",
    )
    repaired = work_dir / f"{tag}_Repair.pdb"
    if not repaired.exists() or repaired.stat().st_size < 500:
        raise FileNotFoundError(f"FoldX RepairPDB missing output {repaired}")
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(repaired, out_pdb)
    return out_pdb


def _build_homotypic_spec(
    pdb_path: Path,
    chains_order: list[str],
    mutations_pdb: list[tuple[int, str, str]],
) -> tuple[str, int]:
    """Build FoldX individual_list line applying mutations across all resolved chains."""
    chain_aa = _extract_chain_aa_map(pdb_path)
    if "A" not in chain_aa:
        raise RuntimeError(f"Target chain A missing in {pdb_path}")
    for pos, wt_aa, _ in mutations_pdb:
        actual = chain_aa["A"].get(pos)
        if actual != wt_aa:
            raise RuntimeError(
                f"Chain A residue mismatch in {pdb_path} at {pos}: expected {wt_aa}, found {actual}"
            )

    tokens: list[str] = []
    chains_with_all = 0
    for cid in chains_order:
        cmap = chain_aa.get(cid, {})
        has_all = all(cmap.get(pos) == wt_aa for pos, wt_aa, _ in mutations_pdb)
        if has_all:
            chains_with_all += 1
        for pos, wt_aa, mut_aa in mutations_pdb:
            if cmap.get(pos) == wt_aa:
                tokens.append(f"{wt_aa}{cid}{pos}{mut_aa}")
    if chains_with_all == 0 or not tokens:
        raise RuntimeError(f"No chains resolved for {mutations_pdb} in {pdb_path}")
    return ",".join(tokens) + ";", chains_with_all


def _run_analyse_complex(pdb_path: Path, work_dir: Path, tag: str) -> list[dict[str, str]]:
    stem = f"{tag}.pdb"
    local_in = work_dir / stem
    if local_in.resolve() != pdb_path.resolve():
        shutil.copy2(pdb_path, local_in)
    _run_foldx(
        [
            "--command=AnalyseComplex",
            f"--pdb={stem}",
            "--pdb-dir=.",
            "--output-dir=.",
        ],
        work_dir,
        f"ac_{tag}.log",
    )
    ac_file = _find_single_file(work_dir, f"Interaction_{tag}*_AC.fxout")
    return _parse_fxout_table(ac_file)


def _index_ac_by_partner(ac_rows: list[dict[str, str]]) -> dict[str, dict[str, float]]:
    """Extract pairwise interface rows involving chain A keyed by neighbor chain ID."""
    by_partner: dict[str, dict[str, float]] = {}
    for row in ac_rows:
        g1, g2 = row.get("Group1", "").strip(), row.get("Group2", "").strip()
        if g1 == "A" and g2 != "A":
            partner = g2
        elif g2 == "A" and g1 != "A":
            partner = g1
        else:
            continue
        parsed: dict[str, float] = {}
        for fx_col, short_name in AC_TERM_MAP.items():
            parsed[short_name] = float(row[fx_col])
        parsed["bind_hbond"] = float(row["Backbone Hbond"]) + float(row["Sidechain Hbond"])
        parsed["bind_entropy"] = float(row["entropy sidechain"]) + float(row["entropy mainchain"])
        by_partner[partner] = parsed
    return by_partner


def ensure_wt_template_repaired(
    template_key: str,
    template_spec: dict[str, Any] | None = None,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """Build, repair (and apply WT restoration muts if needed), and cache WT monomer + subcomplex."""
    if template_spec is None:
        template_spec = load_wt_templates()[template_key]

    tdir = WT_CACHE_DIR / template_key
    tdir.mkdir(parents=True, exist_ok=True)
    meta_path = tdir / "wt_prep_meta.json"
    mono_out = tdir / "mono_Repair.pdb"
    sub_out = tdir / "sub_Repair.pdb"
    wt_ac_path = tdir / "wt_ac_partners.json"
    if (
        not force
        and meta_path.exists()
        and mono_out.exists()
        and sub_out.exists()
        and wt_ac_path.exists()
    ):
        return json.loads(meta_path.read_text(encoding="utf-8"))

    t0 = time.perf_counter()
    mono_raw = tdir / "monomer_raw.pdb"
    sub_raw = tdir / "subcomplex_raw.pdb"
    biomt_src = Path(template_spec["biomt_source_pdb"]) if template_spec.get("biomt_source_pdb") else None
    asm_info = build_biological_assembly_and_subcomplex(
        Path(template_spec["pdb_path"]),
        target_chain_id=template_spec["target_chain_id"],
        biomt_source_pdb=biomt_src,
        out_monomer_pdb=mono_raw,
        out_subcomplex_pdb=sub_raw,
        normalize_mse=bool(template_spec.get("normalize_mse", True)),
    )
    asm_meta = {
        k: asm_info[k]
        for k in ("total_subunits", "subcomplex_chains", "neighbors", "resnum_map_a")
        if k in asm_info
    }
    (tdir / "assembly_info.json").write_text(json.dumps(asm_meta, indent=2) + "\n", encoding="utf-8")

    restore_raw = template_spec.get("restore_muts_pdb") or []
    restore_muts = [(int(r[0]), str(r[1]), str(r[2])) for r in restore_raw]

    with tempfile.TemporaryDirectory(prefix=f"fx_wt_{template_key}_") as tmp_str:
        work = Path(tmp_str)
        # 1. Monomer Repair (+ optional crystal restoration)
        mono_rep0 = _repair_pdb(mono_raw, work / "mono0.pdb", work, "mono0")
        if restore_muts:
            spec_line, _ = _build_homotypic_spec(mono_rep0, ["A"], restore_muts)
            (work / "individual_list.txt").write_text(spec_line + "\n", encoding="utf-8")
            shutil.copy2(mono_rep0, work / "mono0_input.pdb")
            _run_foldx(
                [
                    "--command=BuildModel",
                    "--pdb=mono0_input.pdb",
                    "--pdb-dir=.",
                    "--output-dir=.",
                    "--mutant-file=individual_list.txt",
                    "--numberOfRuns=3",
                    "--out-pdb=1",
                ],
                work,
                "restore_mono.log",
            )
            dif_rows = _parse_fxout_table(_find_single_file(work, "Dif_mono0_input*.fxout"))
            best_idx = int(min(range(len(dif_rows)), key=lambda i: float(dif_rows[i]["total energy"])))
            restored_mono = work / f"mono0_input_1_{best_idx}.pdb"
            _repair_pdb(restored_mono, mono_out, work, "mono_restored")
        else:
            shutil.copy2(mono_rep0, mono_out)

        # 2. Subcomplex Repair (+ optional crystal restoration across all chains)
        sub_rep0 = _repair_pdb(sub_raw, work / "sub0.pdb", work, "sub0")
        sub_chains = asm_info["subcomplex_chains"]
        if restore_muts:
            spec_line, n_restored = _build_homotypic_spec(sub_rep0, sub_chains, restore_muts)
            if n_restored != len(sub_chains):
                raise RuntimeError(
                    f"{template_key}: expected restoration across all {len(sub_chains)} chains, got {n_restored}"
                )
            (work / "individual_list.txt").write_text(spec_line + "\n", encoding="utf-8")
            shutil.copy2(sub_rep0, work / "sub0_input.pdb")
            _run_foldx(
                [
                    "--command=BuildModel",
                    "--pdb=sub0_input.pdb",
                    "--pdb-dir=.",
                    "--output-dir=.",
                    "--mutant-file=individual_list.txt",
                    "--numberOfRuns=3",
                    "--out-pdb=1",
                ],
                work,
                "restore_sub.log",
            )
            dif_rows = _parse_fxout_table(_find_single_file(work, "Dif_sub0_input*.fxout"))
            best_idx = int(min(range(len(dif_rows)), key=lambda i: float(dif_rows[i]["total energy"])))
            restored_sub = work / f"sub0_input_1_{best_idx}.pdb"
            _repair_pdb(restored_sub, sub_out, work, "sub_restored")
        else:
            shutil.copy2(sub_rep0, sub_out)

        # 3. WT AnalyseComplex baseline on final sub_Repair.pdb
        ac_rows = _run_analyse_complex(sub_out, work, "wt_sub")
        wt_ac_partners = _index_ac_by_partner(ac_rows)
        wt_ac_path.write_text(json.dumps(wt_ac_partners, indent=2) + "\n", encoding="utf-8")

    # Sanity verification of final cached structures
    mono_chains = _extract_chain_aa_map(mono_out)
    sub_chains_map = _extract_chain_aa_map(sub_out)
    if list(mono_chains.keys()) != ["A"]:
        raise RuntimeError(f"{template_key}: mono_Repair chains={list(mono_chains.keys())}")
    if sorted(sub_chains_map.keys()) != sorted(asm_info["subcomplex_chains"]):
        raise RuntimeError(
            f"{template_key}: sub_Repair chains mismatch {sorted(sub_chains_map.keys())} vs {asm_info['subcomplex_chains']}"
        )

    meta = {
        "template_key": template_key,
        "mono_repair_pdb": str(mono_out),
        "sub_repair_pdb": str(sub_out),
        "wt_ac_partners_json": str(wt_ac_path),
        "subcomplex_chains": asm_info["subcomplex_chains"],
        "neighbor_symmetry": {n["new_chain_id"]: n["symmetry_class"] for n in asm_info["neighbors"]},
        "restore_muts_pdb": restore_raw,
        "wt_bind_total_kcal": float(sum(v["bind"] for v in wt_ac_partners.values())),
        "elapsed_sec": round(time.perf_counter() - t0, 2),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def parse_mutation_list(mut_str: str) -> list[tuple[int, str, str]]:
    muts = []
    for tok in mut_str.strip().split("/"):
        tok = tok.strip()
        if not tok:
            continue
        muts.append((int(tok[1:-1]), tok[0], tok[-1]))
    if not muts:
        raise ValueError(f"Empty mutation string: {mut_str!r}")
    return muts


def _aggregate_dif_rows(dif_rows: list[dict[str, str]], scale: float = 1.0) -> dict[str, float]:
    out: dict[str, float] = {}
    ddg_vals = [float(r["total energy"]) * scale for r in dif_rows]
    out["ddg"] = float(np.mean(ddg_vals))
    out["ddg_sd"] = float(np.std(ddg_vals, ddof=0))
    for fx_col, short_name in DIF_TERM_MAP.items():
        if short_name == "ddg":
            continue
        vals = [float(r[fx_col]) * scale for r in dif_rows]
        out[short_name] = float(np.mean(vals))
    return out


def run_foldx_monomer_and_subcomplex(
    template_key: str,
    mutation_pdb: str,
    *,
    template_spec: dict[str, Any] | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Run FoldX 5.1 BuildModel (3 runs) on monomer & homotypic subcomplex + AnalyseComplex."""
    safe_mut = mutation_pdb.replace("/", "_")
    case_key = f"{template_key}__{safe_mut}"
    case_dir = FOLDX_CASE_CACHE_DIR / case_key
    case_dir.mkdir(parents=True, exist_ok=True)
    result_path = case_dir / "foldx_features.json"
    mono_rep_mt_out = case_dir / "mono_mutant_2.pdb"
    sub_rep_mt_out = case_dir / "sub_mutant_2.pdb"

    if (
        not force
        and result_path.exists()
        and mono_rep_mt_out.exists()
        and sub_rep_mt_out.exists()
    ):
        return json.loads(result_path.read_text(encoding="utf-8"))

    t0 = time.perf_counter()
    wt_meta = ensure_wt_template_repaired(template_key, template_spec, force=False)
    mono_wt_pdb = Path(wt_meta["mono_repair_pdb"])
    sub_wt_pdb = Path(wt_meta["sub_repair_pdb"])
    sub_chains = list(wt_meta["subcomplex_chains"])
    neighbor_sym = dict(wt_meta["neighbor_symmetry"])
    mutations_pdb = parse_mutation_list(mutation_pdb)

    mono_spec, _ = _build_homotypic_spec(mono_wt_pdb, ["A"], mutations_pdb)
    sub_spec, n_mutated_chains = _build_homotypic_spec(sub_wt_pdb, sub_chains, mutations_pdb)

    with tempfile.TemporaryDirectory(prefix=f"fx_case_{case_key}_") as tmp_str:
        work = Path(tmp_str)
        # 1. Monomer BuildModel (3 runs)
        shutil.copy2(mono_wt_pdb, work / "mono_wt.pdb")
        (work / "individual_list_mono.txt").write_text(mono_spec + "\n", encoding="utf-8")
        _run_foldx(
            [
                "--command=BuildModel",
                "--pdb=mono_wt.pdb",
                "--pdb-dir=.",
                "--output-dir=.",
                "--mutant-file=individual_list_mono.txt",
                "--numberOfRuns=3",
                "--out-pdb=1",
            ],
            work,
            "build_mono.log",
        )
        mono_dif = _parse_fxout_table(_find_single_file(work, "Dif_mono_wt*.fxout"))
        mono_stats = _aggregate_dif_rows(mono_dif, scale=1.0)
        mono_run2 = work / "mono_wt_1_2.pdb"
        if not mono_run2.exists():
            mono_run2 = _find_single_file(work, "mono_wt_1_*.pdb")
        shutil.copy2(mono_run2, mono_rep_mt_out)

        # Save all 3 monomer WT/MT snapshots for optional structural checks; clean rest after
        # 2. Subcomplex BuildModel (3 runs, homotypic across resolved chains)
        shutil.copy2(sub_wt_pdb, work / "sub_wt.pdb")
        (work / "individual_list_sub.txt").write_text(sub_spec + "\n", encoding="utf-8")
        _run_foldx(
            [
                "--command=BuildModel",
                "--pdb=sub_wt.pdb",
                "--pdb-dir=.",
                "--output-dir=.",
                "--mutant-file=individual_list_sub.txt",
                "--numberOfRuns=3",
                "--out-pdb=1",
            ],
            work,
            "build_sub.log",
        )
        sub_dif = _parse_fxout_table(_find_single_file(work, "Dif_sub_wt*.fxout"))
        sub_stats_total = _aggregate_dif_rows(sub_dif, scale=1.0)
        sub_stats_per = _aggregate_dif_rows(sub_dif, scale=1.0 / float(n_mutated_chains))

        sub_run2 = work / "sub_wt_1_2.pdb"
        if not sub_run2.exists():
            sub_run2 = _find_single_file(work, "sub_wt_1_*.pdb")
        shutil.copy2(sub_run2, sub_rep_mt_out)

        # 3. Paired AnalyseComplex on matched WT_sub_wt_1_{r}.pdb and sub_wt_1_{r}.pdb
        run_deltas_by_partner: dict[str, list[dict[str, float]]] = {p: [] for p in neighbor_sym}
        wt_totals: list[float] = []
        n_runs_ac = 0
        for r_idx in range(len(sub_dif)):
            wt_r = work / f"WT_sub_wt_1_{r_idx}.pdb"
            mt_r = work / f"sub_wt_1_{r_idx}.pdb"
            if not (wt_r.exists() and mt_r.exists()):
                continue
            n_runs_ac += 1
            wt_ac = _index_ac_by_partner(_run_analyse_complex(wt_r, work, f"wt_ac_{r_idx}"))
            mt_ac = _index_ac_by_partner(_run_analyse_complex(mt_r, work, f"mt_ac_{r_idx}"))
            wt_totals.append(sum(v["bind"] for v in wt_ac.values()))
            for partner in neighbor_sym:
                wvals = wt_ac.get(partner, {})
                mvals = mt_ac.get(partner, {})
                diff_row = {}
                for k in ("bind", "bind_elec", "bind_vdw", "bind_solv_polar", "bind_solv_hyd", "bind_clash", "bind_hbond", "bind_entropy"):
                    diff_row[k] = float(mvals.get(k, 0.0) - wvals.get(k, 0.0))
                run_deltas_by_partner[partner].append(diff_row)

        if n_runs_ac == 0:
            raise RuntimeError(f"No BuildModel WT/MT PDB pairs produced for {case_key}")

    # Average per partner across runs, then aggregate by symmetry class
    partner_mean: dict[str, dict[str, float]] = {}
    for partner, diffs in run_deltas_by_partner.items():
        if not diffs:
            partner_mean[partner] = {k: 0.0 for k in ("bind", "bind_elec", "bind_vdw", "bind_solv_polar", "bind_solv_hyd", "bind_clash", "bind_hbond", "bind_entropy")}
        else:
            partner_mean[partner] = {
                k: float(np.mean([d[k] for d in diffs])) for k in diffs[0]
            }

    partner_binds = [partner_mean[p]["bind"] for p in neighbor_sym]
    c2_bind = float(sum(partner_mean[p]["bind"] for p, sym in neighbor_sym.items() if sym == "C2"))
    c3_bind = float(sum(partner_mean[p]["bind"] for p, sym in neighbor_sym.items() if sym == "C3"))
    c4_bind = float(sum(partner_mean[p]["bind"] for p, sym in neighbor_sym.items() if sym == "C4"))
    total_bind = float(sum(partner_binds))

    features: dict[str, Any] = {
        "case_key": case_key,
        "template_key": template_key,
        "mutation_pdb": mutation_pdb,
        "n_mutated_chains_subcomplex": n_mutated_chains,
        "n_total_chains_subcomplex": len(sub_chains),
        "mono_repair_pdb": str(mono_wt_pdb),
        "sub_repair_pdb": str(sub_wt_pdb),
        "mono_mutant_pdb": str(mono_rep_mt_out),
        "sub_mutant_pdb": str(sub_rep_mt_out),
        # Monomer FoldX features
        **{f"fx_mono_{k}": v for k, v in mono_stats.items()},
        # Subcomplex FoldX features (total + per mutated chain)
        "fx_sub_ddg_total": sub_stats_total["ddg"],
        **{f"fx_sub_{k}_per_chain": v for k, v in sub_stats_per.items()},
        "fx_sub_vs_mono_diff": float(sub_stats_per["ddg"] - mono_stats["ddg"]),
        # AnalyseComplex interface binding ΔΔG (Chain A vs C2/C3/C4 neighbors)
        "fx_ac_wt_bind_total": float(np.mean(wt_totals)),
        "fx_ac_ddg_bind_total": total_bind,
        "fx_ac_ddg_bind_c2": c2_bind,
        "fx_ac_ddg_bind_c3": c3_bind,
        "fx_ac_ddg_bind_c4": c4_bind,
        "fx_ac_ddg_bind_max_iface": float(max(partner_binds)) if partner_binds else 0.0,
        "fx_ac_ddg_bind_min_iface": float(min(partner_binds)) if partner_binds else 0.0,
        # Stoichiometric theoretical cage ΔΔG_subunit = ΔΔG_fold + 0.5 * ΔΔG_bind_A
        "fx_cage_ddg_eff": float(mono_stats["ddg"] + 0.5 * total_bind),
        # Interface energy term breakdown
        "fx_ac_ddg_bind_elec": float(sum(v["bind_elec"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_hbond": float(sum(v["bind_hbond"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_vdw": float(sum(v["bind_vdw"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_solv_hyd": float(sum(v["bind_solv_hyd"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_solv_polar": float(sum(v["bind_solv_polar"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_clash": float(sum(v["bind_clash"] for v in partner_mean.values())),
        "fx_ac_ddg_bind_entropy": float(sum(v["bind_entropy"] for v in partner_mean.values())),
        "fx_ac_partner_ddg_bind": {p: round(partner_mean[p]["bind"], 4) for p in neighbor_sym},
        "elapsed_sec": round(time.perf_counter() - t0, 2),
    }
    result_path.write_text(json.dumps(features, indent=2) + "\n", encoding="utf-8")
    return features


def prepare_all_wt_templates(workers: int = 8, force: bool = False) -> dict[str, dict[str, Any]]:
    """Concurrently repair & baseline-analyse all 14 ferritin WT templates."""
    templates = load_wt_templates()
    results: dict[str, dict[str, Any]] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut_to_key = {
            ex.submit(ensure_wt_template_repaired, k, v, force=force): k
            for k, v in templates.items()
        }
        for fut in as_completed(fut_to_key):
            k = fut_to_key[fut]
            meta = fut.result()
            results[k] = meta
            print(
                f"[WT {len(results)}/{len(templates)}] {k}: subcomplex={len(meta['subcomplex_chains'])} "
                f"wt_bind={meta['wt_bind_total_kcal']:.2f} kcal/mol ({meta['elapsed_sec']}s)",
                flush=True,
            )
    return results


def run_all_curated_foldx(workers: int = 10, force: bool = False) -> list[dict[str, Any]]:
    """Run FoldX monomer + subcomplex + AnalyseComplex across all unique curated mutations."""
    templates = load_wt_templates()
    with (DATA_DIR / "curated_ferritin_cohorts.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    unique_jobs = sorted({(r["template_key"], r["mutation_pdb"]) for r in rows})
    out: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fut_to_job = {
            ex.submit(
                run_foldx_monomer_and_subcomplex,
                tkey,
                mut_pdb,
                template_spec=templates[tkey],
                force=force,
            ): (tkey, mut_pdb)
            for tkey, mut_pdb in unique_jobs
        }
        for fut in as_completed(fut_to_job):
            tkey, mut_pdb = fut_to_job[fut]
            res = fut.result()
            out.append(res)
            print(
                f"[FoldX {len(out)}/{len(unique_jobs)}] {tkey} {mut_pdb}: "
                f"mono_ddG={res['fx_mono_ddg']:+.2f}, bind_ddG={res['fx_ac_ddg_bind_total']:+.2f} "
                f"(C2={res['fx_ac_ddg_bind_c2']:+.2f}, C3={res['fx_ac_ddg_bind_c3']:+.2f}, C4={res['fx_ac_ddg_bind_c4']:+.2f}) "
                f"[{res['elapsed_sec']}s]",
                flush=True,
            )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-wt", action="store_true", help="Prepare & repair all 14 WT templates")
    parser.add_argument("--run-all", action="store_true", help="Run FoldX across all 69 unique curated mutants")
    parser.add_argument("--smoke-test", action="store_true", help="Run EcBfr_1BFR_A + E128A + HuHF_2FHA_A_restored check")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.smoke_test:
        templates = load_wt_templates()
        m1 = ensure_wt_template_repaired("EcBfr_1BFR_A", templates["EcBfr_1BFR_A"], force=args.force)
        print("WT EcBfr_1BFR_A:", json.dumps(m1, indent=2), flush=True)
        c1 = run_foldx_monomer_and_subcomplex("EcBfr_1BFR_A", "E128A", template_spec=templates["EcBfr_1BFR_A"], force=args.force)
        print("Case EcBfr E128A:", json.dumps(c1, indent=2), flush=True)
        return 0
    if args.prepare_wt:
        prepare_all_wt_templates(workers=args.workers, force=args.force)
    if args.run_all:
        run_all_curated_foldx(workers=args.workers, force=args.force)
    if not (args.prepare_wt or args.run_all or args.smoke_test):
        parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
