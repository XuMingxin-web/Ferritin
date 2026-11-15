"""Curate stratified, sequence-verified ferritin cohorts across S1-S5."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from src.assembly_builder import build_biological_assembly_and_subcomplex, extract_site_interface_features
from src.paths import DATA_DIR, EXTERNAL_GEODTM_ROOT, FERRO_DATASET_CSV

EXPANSION_CURATED_DIR = EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/curated"
RAW_PDB_DIR = DATA_DIR / "ferritin_pdb_raw"

WT_TEMPLATES: dict[str, dict[str, Any]] = {
    "EcBfr_1BFR_A": {
        "protein_group": "P0ABD3",
        "protein_name": "Bacterioferritin (EcBfr)",
        "oligomer_family": "24mer_bfr",
        "pdb_id": "1BFR",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/raw/structures/1BFR.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "AfFtn_1SQ3_A": {
        "protein_group": "O29424",
        "protein_name": "Archaeoglobus fulgidus ferritin (AfFtn)",
        "oligomer_family": "24mer_archaeal",
        "pdb_id": "1SQ3",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/raw/structures/1SQ3.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "RcBFR_1JGC_A": {
        "protein_group": "Q59738",
        "protein_name": "Rhodobacter capsulatus BFR (RcBFR)",
        "oligomer_family": "24mer_bfr",
        "pdb_id": "1JGC",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/raw/structures/1JGC.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "EcFtnA_1EUM_A": {
        "protein_group": "P0A998",
        "protein_name": "E. coli FtnA (EcFtnA)",
        "oligomer_family": "24mer_nonheme_bac",
        "pdb_id": "1EUM",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/followup13_20261115/sources/1EUM.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "HuHF_2FHA_A_restored": {
        "protein_group": "P02794",
        "protein_name": "Human H-ferritin (HuHF)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "2FHA",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/2FHA.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [[86, "Q", "K"]],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "MjFer_6A4U_A_restored": {
        "protein_group": "T2B7E1",
        "protein_name": "Kuruma prawn ferritin (MjFer)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "6A4U",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/6A4U.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [[89, "R", "Q"]],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "CgFer_6LIJ_A": {
        "protein_group": "Q70MM3",
        "protein_name": "Pacific oyster ferritin (CgFer)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "6LIJ",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/6LIJ.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "AjFER_7VHR_B": {
        "protein_group": "Q2QJV4",
        "protein_name": "Sea cucumber ferritin (AjFER)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "7VHR",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/7VHR.pdb"),
        "target_chain_id": "B",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "HuFTL_2FFX_J": {
        "protein_group": "P02792",
        "protein_name": "L-chain ferritin proxy (2FFX_J)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "2FFX",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/2FFX.pdb"),
        "target_chain_id": "J",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "EcDps_1F33_F": {
        "protein_group": "P0ABT2",
        "protein_name": "E. coli Dps miniferritin (EcDps)",
        "oligomer_family": "12mer_dps",
        "pdb_id": "1F33",
        "pdb_path": str(EXTERNAL_GEODTM_ROOT / "data/expansion_audit/mtdatabase_20261115/structure_sources/1F33.pdb"),
        "target_chain_id": "F",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 12,
    },
    "BfMFtn_1MFR_A": {
        "protein_group": "P07798",
        "protein_name": "Bullfrog M-ferritin (BfMFtn)",
        "oligomer_family": "24mer_canonical_euk",
        "pdb_id": "1MFR",
        "pdb_path": str(RAW_PDB_DIR / "1MFR.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "MtBfrA_3UOI_A": {
        "protein_group": "P9WPQ9",
        "protein_name": "M. tuberculosis BfrA (MtBfrA)",
        "oligomer_family": "24mer_bfr",
        "pdb_id": "3UOI",
        "pdb_path": str(RAW_PDB_DIR / "3UOI.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "MtBfrB_3QD8_A": {
        "protein_group": "P9WNE5",
        "protein_name": "M. tuberculosis BfrB (MtBfrB)",
        "oligomer_family": "24mer_canonical_bac",
        "pdb_id": "3QD8",
        "pdb_path": str(RAW_PDB_DIR / "3QD8.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
    "PfFtn_2JD6_A": {
        "protein_group": "Q8U2T8",
        "protein_name": "Pyrococcus furiosus ferritin (PfFtn)",
        "oligomer_family": "24mer_archaeal",
        "pdb_id": "2JD6",
        "pdb_path": str(RAW_PDB_DIR / "2JD6.pdb"),
        "target_chain_id": "A",
        "biomt_source_pdb": None,
        "restore_muts_pdb": [],
        "normalize_mse": True,
        "expected_subunits": 24,
    },
}

PROTEIN_TO_TEMPLATE = {
    "P0ABD3": "EcBfr_1BFR_A",
    "O29424": "AfFtn_1SQ3_A",
    "Q59738": "RcBFR_1JGC_A",
    "P0A998": "EcFtnA_1EUM_A",
    "P02794": "HuHF_2FHA_A_restored",
    "T2B7E1": "MjFer_6A4U_A_restored",
    "Q70MM3": "CgFer_6LIJ_A",
    "Q2QJV4": "AjFER_7VHR_B",
    "P02792": "HuFTL_2FFX_J",
    "P0ABT2": "EcDps_1F33_F",
    "P07798": "BfMFtn_1MFR_A",
    "P9WPQ9": "MtBfrA_3UOI_A",
    "P9WNE5": "MtBfrB_3QD8_A",
    "Q8U2T8": "PfFtn_2JD6_A",
}


def _make_aligner() -> PairwiseAligner:
    aln = PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 2.0
    aln.mismatch_score = -1.0
    aln.open_gap_score = -10.0
    aln.extend_gap_score = -0.5
    return aln


def _extract_chain_residues(pdb_path: Path, chain_id: str, restore_muts_pdb: list[list[Any]]) -> tuple[list[int], str]:
    parser = PDBParser(QUIET=True)
    ch = parser.get_structure(pdb_path.stem, str(pdb_path))[0][chain_id]
    restore_map = {int(pos): mut_aa for pos, _, mut_aa in restore_muts_pdb}
    resnums: list[int] = []
    aas: list[str] = []
    for res in ch:
        if res.id[0] not in (" ", "H_MSE") or "CA" not in res:
            continue
        rnum = int(res.id[1])
        aa = restore_map.get(rnum, "M" if res.resname.strip() == "MSE" else seq1(res.resname.strip()))
        resnums.append(rnum)
        aas.append(aa)
    return resnums, "".join(aas)


def _align_seq_to_pdb(ref_seq: str, pdb_resnums: list[int], pdb_seq: str) -> dict[int, tuple[int, str]]:
    """Map 1-based ref_seq positions to (pdb_resnum, pdb_aa) pairs."""
    aln = _make_aligner().align(ref_seq, pdb_seq)[0]
    mapping: dict[int, tuple[int, str]] = {}
    for (qa, qb), (ta, tb) in zip(*aln.aligned):
        for k in range(qb - qa):
            mapping[qa + k + 1] = (pdb_resnums[ta + k], pdb_seq[ta + k])
    return mapping


def _parse_mutations(mut_str: str) -> list[tuple[int, str, str]]:
    out = []
    for tok in mut_str.strip().split("/"):
        tok = tok.strip()
        if not tok:
            continue
        wt_aa = tok[0]
        mut_aa = tok[-1]
        pos1 = int(tok[1:-1])
        out.append((pos1, wt_aa, mut_aa))
    if not out:
        raise ValueError(f"Empty mutation string: {mut_str!r}")
    return out


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def _apply_mutations_to_seq(seq: str, muts: list[tuple[int, str, str]]) -> str:
    chars = list(seq)
    for pos1, wt_aa, mut_aa in muts:
        if chars[pos1 - 1] != wt_aa:
            raise ValueError(f"Sequence mismatch at pos {pos1}: expected {wt_aa}, found {chars[pos1 - 1]}")
        chars[pos1 - 1] = mut_aa
    return "".join(chars)


def build_curated_cohorts() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Assemble and sequence-verify S1-S5 ferritin benchmark rows."""
    expanded_rows = _read_csv(EXPANSION_CURATED_DIR / "expanded_author_delta.csv")
    dps_rows = _read_csv(EXPANSION_CURATED_DIR / "dps_external.csv")
    ti_rows = _read_csv(EXPANSION_CURATED_DIR / "fluorescence_Ti.csv")
    het_rows = _read_csv(EXPANSION_CURATED_DIR / "heteropolymer.csv")
    dsc1_rows = _read_csv(EXPANSION_CURATED_DIR / "dsc_peak1.csv")
    dsc2_rows = _read_csv(EXPANSION_CURATED_DIR / "dsc_peak2.csv")

    raw_mtdb = {r["record_id"]: r for r in _read_csv(FERRO_DATASET_CSV)}

    # Assemble PDB biological assemblies once per WT template for sequence + structural annotation
    templates_meta: dict[str, Any] = {}
    template_assemblies: dict[str, dict] = {}
    template_pdb_seq: dict[str, tuple[list[int], str]] = {}
    for tkey, tspec in WT_TEMPLATES.items():
        pdb_path = Path(tspec["pdb_path"])
        info = build_biological_assembly_and_subcomplex(
            pdb_path,
            target_chain_id=tspec["target_chain_id"],
            biomt_source_pdb=Path(tspec["biomt_source_pdb"]) if tspec["biomt_source_pdb"] else None,
            normalize_mse=tspec["normalize_mse"],
        )
        assert info["total_subunits"] == tspec["expected_subunits"], (
            f"{tkey}: expected {tspec['expected_subunits']} subunits, got {info['total_subunits']}"
        )
        resnums, pdb_seq = _extract_chain_residues(pdb_path, tspec["target_chain_id"], tspec["restore_muts_pdb"])
        template_assemblies[tkey] = info
        template_pdb_seq[tkey] = (resnums, pdb_seq)
        templates_meta[tkey] = {
            **tspec,
            "total_subunits": info["total_subunits"],
            "subcomplex_chains": info["subcomplex_chains"],
            "neighbors": info["neighbors"],
            "pdb_seq": pdb_seq,
            "pdb_resnum_min": min(resnums),
            "pdb_resnum_max": max(resnums),
        }

    records: list[dict[str, Any]] = []

    def add_record(
        *,
        case_id: str,
        cohort: str,
        subcohort: str,
        uniprot_id: str,
        wt_seq: str,
        mut_canon: str,
        pdb_mut_override: str | None = None,
        delta_tm: float | None,
        delta_tm_alt: float | None = None,
        ph: float,
        study_id: str,
        endpoint_type: str,
        binary_direction: int | None = None,
        source_record_id: str = "",
        notes: str = "",
    ) -> None:
        tkey = PROTEIN_TO_TEMPLATE[uniprot_id]
        pdb_resnums, pdb_seq = template_pdb_seq[tkey]
        seq_muts = _parse_mutations(mut_canon)
        mut_seq = _apply_mutations_to_seq(wt_seq, seq_muts)
        seq_to_pdb = _align_seq_to_pdb(wt_seq, pdb_resnums, pdb_seq)

        if pdb_mut_override:
            pdb_muts = _parse_mutations(pdb_mut_override)
            pdb_map_by_rnum = dict(zip(pdb_resnums, pdb_seq))
            for rnum, w_aa, _ in pdb_muts:
                if pdb_map_by_rnum.get(rnum) != w_aa:
                    raise ValueError(f"{case_id}: PDB mismatch at {tkey} resnum {rnum}: {pdb_map_by_rnum.get(rnum)} != {w_aa}")
        else:
            pdb_muts = []
            for pos1, w_aa, m_aa in seq_muts:
                rnum, p_aa = seq_to_pdb[pos1]
                if p_aa != w_aa:
                    raise ValueError(f"{case_id}: aligned PDB mismatch at {tkey} seq {pos1}: PDB({rnum})={p_aa} != {w_aa}")
                pdb_muts.append((rnum, w_aa, m_aa))

        iface_feats = extract_site_interface_features(
            template_assemblies[tkey],
            [m[0] for m in pdb_muts],
            pdb_muts,
        )
        # Re-normalize directional class for 0 / +/- thresholds
        if binary_direction is None and delta_tm is not None:
            binary_direction = 1 if delta_tm > 0.0 else 0

        row = {
            "case_id": case_id,
            "cohort": cohort,
            "subcohort": subcohort,
            "template_key": tkey,
            "protein_group": WT_TEMPLATES[tkey]["protein_group"],
            "oligomer_family": WT_TEMPLATES[tkey]["oligomer_family"],
            "study_id": study_id,
            "endpoint_type": endpoint_type,
            "source_record_id": source_record_id,
            "wt_sequence": wt_seq,
            "mut_sequence": mut_seq,
            "mutation_canon": mut_canon,
            "mutation_pdb": "/".join(f"{w}{p}{m}" for p, w, m in pdb_muts),
            "n_mutations": len(pdb_muts),
            "pH": float(ph),
            "delta_tm": "" if delta_tm is None else float(delta_tm),
            "delta_tm_alt_thesis": "" if delta_tm_alt is None else float(delta_tm_alt),
            "binary_direction": "" if binary_direction is None else int(binary_direction),
            "notes": notes,
        }
        for k, v in iface_feats.items():
            row[k] = round(float(v), 6)
        records.append(row)

    # Known conflict alternative values from DR-NTU 48015 Tables 2.7 / 3.3
    conflict_alt = {
        "ZHANG2010_R30A": 7.4,
        "ZHANG2015_R30D": -9.0,
        "ZHANG2015_D56R": -3.8,
        "ZHANG2015_E60H": 4.1,
    }
    # 1. S1 Main Single-Point 24-mer cohort (23 paired_core + 2 patents + 1 author cross-study + 4 MtBfrA ESI sensitivity)
    for r in expanded_rows:
        rid = r["record_id"]
        st = r.get("analysis_stratum") or r.get("recommended_analysis_stratum", "")
        if st in ("prior_journal", "thesis_extension", "paired_journal_addition"):
            subcohort = "S1_PairedCore23"
        elif st == "patent":
            subcohort = "S1_ExtraPatent2"
        elif st == "author_cross_study_sensitivity":
            subcohort = "S1_ExtraOyster1"
        else:
            subcohort = f"S1_{st}"
        add_record(
            case_id=rid,
            cohort="S1_SinglePoint_24mer",
            subcohort=subcohort,
            uniprot_id=r["uniprot_id"],
            wt_seq=r["sequence"],
            mut_canon=r["mutation"],
            delta_tm=float(r["delta_tm"]),
            delta_tm_alt=conflict_alt.get(rid),
            ph=float(r.get("prediction_pH") or r.get("pH") or 7.0),
            study_id=r["source_study_id"],
            endpoint_type=r.get("endpoint", "delta_Tm"),
            source_record_id=r.get("source_record_id") or r.get("input_record_id") or "",
            notes=r.get("data_quality_flag", ""),
        )

    # Add MtBfrA E121A/Q/K/F (ESI Fig S13 approx, P9WPQ9)
    seq_mtbfra = template_pdb_seq["MtBfrA_3UOI_A"][1]
    for fid in ["FMT0216", "FMT0217", "FMT0218", "FMT0219"]:
        r = raw_mtdb[fid]
        add_record(
            case_id=f"MTDB_{fid}",
            cohort="S1_SinglePoint_24mer",
            subcohort="S1_Sensitivity_MtBfrA_ESI4",
            uniprot_id="P9WPQ9",
            wt_seq=seq_mtbfra,
            mut_canon=r["mutation_normalized"].split()[0],
            delta_tm=float(r["normalized_delta_Tm_C"]),
            ph=7.0,
            study_id="2025_ChemSci_PMC11776510",
            endpoint_type="delta_Tm_CD_ESI",
            source_record_id=fid,
            notes="ESI Fig S13 approximate transition",
        )

    # 2. S2 External Family Transfer: 12-mer EcDps (8 single mutants)
    for r in dps_rows:
        add_record(
            case_id=r["record_id"],
            cohort="S2_EcDps_12mer",
            subcohort="S2_EcDps_Single8",
            uniprot_id=r["uniprot_id"],
            wt_seq=r["sequence"],
            mut_canon=r["mutation"],
            delta_tm=float(r["delta_tm"]),
            ph=float(r.get("prediction_pH") or r.get("pH") or 7.2),
            study_id=r["source_study_id"],
            endpoint_type="delta_Tm_CD",
            source_record_id=r.get("source_record_id", ""),
        )

    # 3. S3 Combinatorial Multi-Point Tier A (Exact 9 + Abstract series 7 = 16)
    seq_ecbfr = [r["sequence"] for r in expanded_rows if r["uniprot_id"] == "P0ABD3"][0]
    seq_rcbfr = [r["sequence"] for r in expanded_rows if r["uniprot_id"] == "Q59738"][0]
    seq_ecdps = dps_rows[0]["sequence"]
    seq_ajfer = ti_rows[0]["sequence"]
    seq_huhf = [r["sequence"] for r in expanded_rows if r["uniprot_id"] == "P02794"][0]
    seq_bfmftn = template_pdb_seq["BfMFtn_1MFR_A"][1]

    s3_specs = [
        # (case_id, subcohort, uniprot, wt_seq, mut_canon, pdb_override, dTm, alt_dTm, ph, study, ep, fmt_id)
        ("MTDB_FMT0010_EcBfr_D56A_E60A", "S3_Exact9", "P0ABD3", seq_ecbfr, "D56A/E60A", None, -16.2, None, 7.2, "10.1074/jbc.m109.092445", "delta_Tm_CD", "FMT0010"),
        ("MTDB_FMT0017_EcBfr_D56R_E60H", "S3_Exact9", "P0ABD3", seq_ecbfr, "D56R/E60H", None, -16.0, -7.4, 7.2, "10.1093/jb/mvv065", "delta_Tm_CD", "FMT0016/FMT0017"),
        ("MTDB_FMT0176_EcBfr_D132W_N34W", "S3_Exact9", "P0ABD3", seq_ecbfr, "D132W/N34W", None, -12.1, None, 7.2, "10.3390/molecules22122184", "delta_Tm_CD", "FMT0176"),
        ("MTDB_FMT0040_EcDps_R83A_R133A", "S3_Exact9", "P0ABT2", seq_ecdps, "R83A/R133A", None, -25.3, None, 7.2, "ZHANG2011_PROTEIN_SCI", "delta_Tm_CD", "FMT0040"),
        ("MTDB_FMT0165_RcBFR_E128A_E135A", "S3_Exact9", "Q59738", seq_rcbfr, "E128A/E135A", None, -26.1, None, 7.2, "10.1110/ps.0301903", "delta_Tm_CD", "FMT0165"),
        ("MTDB_FMT0166_RcBFR_E128R_E135R", "S3_Exact9", "Q59738", seq_rcbfr, "E128R/E135R", None, -29.9, None, 7.2, "10.1110/ps.0301903", "delta_Tm_CD", "FMT0166"),
        ("MTDB_FMT0205_AjFER_MF", "S3_Exact9", "Q2QJV4", seq_ajfer, "E25A/E60A/E105A", None, -18.9, None, 8.0, "WU2022_AJFER", "delta_Ti_Tycho", "FMT0205"),
        ("MTDB_FMT0206_AjFER_M3", "S3_Exact9", "Q2QJV4", seq_ajfer, "D129A/E132A", None, -1.3, None, 8.0, "WU2022_AJFER", "delta_Ti_Tycho", "FMT0206"),
        # HFt-W6 uses canonical P02794 numbering WITH Met1 -> 2FHA resnums are pos-1 (57,60,128,137,151,170)
        ("MTDB_FMT0221_HuHF_W6", "S3_Exact9", "P02794", seq_huhf, "H58W/H61W/H129W/Y138W/H152W/F171W", "H57W/H60W/H128W/Y137W/H151W/F170W", -22.0, None, 7.4, "2026_CIMB_HFtW6", "delta_Tagg_onset", "FMT0221"),
        # Abstract-derived Bullfrog M-ferritin 4-fold pore series
        ("MTDB_FMT0252_BfMFtn_9WX5", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/M161E", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0252"),
        ("MTDB_FMT0253_BfMFtn_9WY3", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/L165D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0253"),
        ("MTDB_FMT0254_BfMFtn_9X0V", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/M161E/L165D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0254"),
        ("MTDB_FMT0255_BfMFtn_9X15", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/L165D/K168E/H169D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0255"),
        ("MTDB_FMT0256_BfMFtn_9X16", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/M161Q/L165D/K168E/H169D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0256"),
        ("MTDB_FMT0257_BfMFtn_9X1G", "S3_Abstract7", "P07798", seq_bfmftn, "M161Q/L165D/K168E/H169D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0257"),
        ("MTDB_FMT0258_BfMFtn_9X6R", "S3_Abstract7", "P07798", seq_bfmftn, "E130A/K168E/H169D", None, -5.0, None, 7.0, "2026_Biochem_41700891", "delta_Tm_abstract", "FMT0258"),
    ]
    for cid, subc, uid, wseq, mc, po, dt, alt_dt, ph_v, stdy, ep, srid in s3_specs:
        add_record(
            case_id=cid, cohort="S3_MultiPoint", subcohort=subc, uniprot_id=uid, wt_seq=wseq,
            mut_canon=mc, pdb_mut_override=po, delta_tm=dt, delta_tm_alt=alt_dt, ph=ph_v,
            study_id=stdy, endpoint_type=ep, source_record_id=srid,
        )

    # 4. S4 Secondary endpoints (Fluorescence Ti AjFER E168A, Heteropolymer HuFTL L135P, DSC dimeric peaks EcBfr)
    for r in ti_rows:
        add_record(
            case_id=r["record_id"], cohort="S4_Secondary", subcohort="S4_Ti", uniprot_id=r["uniprot_id"],
            wt_seq=r["sequence"], mut_canon=r["mutation"], delta_tm=float(r["delta_tm"]),
            ph=float(r.get("prediction_pH") or 7.0), study_id=r["source_study_id"],
            endpoint_type="delta_Ti_Tycho", source_record_id=r.get("source_record_id", ""),
        )
    for r in het_rows:
        add_record(
            case_id=r["record_id"], cohort="S4_Secondary", subcohort="S4_Heteropolymer", uniprot_id=r["uniprot_id"],
            wt_seq=r["sequence"], mut_canon=r["mutation"], delta_tm=float(r["delta_tm"]),
            ph=float(r.get("prediction_pH") or 7.0), study_id=r["source_study_id"],
            endpoint_type="delta_Tm_NanoDSC_heteropolymer", source_record_id=r.get("source_record_id", ""),
        )
    for idx_tag, rlist in [("S4_DSC_Peak1", dsc1_rows), ("S4_DSC_Peak2", dsc2_rows)]:
        for r in rlist:
            add_record(
                case_id=f"{r['record_id']}_{idx_tag}", cohort="S4_Secondary", subcohort=idx_tag,
                uniprot_id=r["uniprot_id"], wt_seq=r["sequence"], mut_canon=r["mutation"],
                delta_tm=float(r["delta_tm"]), ph=float(r.get("prediction_pH") or 7.2),
                study_id=r["source_study_id"], endpoint_type=f"delta_{idx_tag}_dimer",
                source_record_id=r.get("source_record_id", ""),
            )

    # 5. S5 Tier B/C Unseen Directional Qualitative Validation (unambiguous >WT/intact=1 vs <WT=0)
    seq_mtbfrb_pdb = template_pdb_seq["MtBfrB_3QD8_A"][1]
    seq_pfftn_pdb = template_pdb_seq["PfFtn_2JD6_A"][1]
    seq_huftl = het_rows[0]["sequence"]

    s5_specs = [
        # Stabilizing or WT-like stable (binary_direction = 1)
        ("MTDB_FMT0283_EcBfr_N23L", "P0ABD3", seq_ecbfr, "N23L", None, 7.2, "ARDEJANI2011_BIOCHEM", "qualitative_CD_chem", 1, "FMT0283", "> WT"),
        ("MTDB_FMT0284_EcBfr_N23F", "P0ABD3", seq_ecbfr, "N23F", None, 7.2, "ARDEJANI2011_BIOCHEM", "qualitative_CD_chem", 1, "FMT0284", "> WT"),
        ("MTDB_FMT0285_EcBfr_N23W", "P0ABD3", seq_ecbfr, "N23W", None, 7.2, "ARDEJANI2011_BIOCHEM", "qualitative_CD_chem", 1, "FMT0285", "> WT"),
        # 3QD8 starts at resnum 10 -> use pdb_seq with pdb_mut_override for MtBfrB canonical Q51/Q126/E135
        ("MTDB_FMT0286_MtBfrB_Q51A", "P9WNE5", seq_mtbfrb_pdb, "Q42A", "Q51A", 7.2, "KHARE2013_BIOCHEM", "qualitative_thermal", 1, "FMT0286", "> WT"),
        ("MTDB_FMT0287_MtBfrB_Q126A", "P9WNE5", seq_mtbfrb_pdb, "Q117A", "Q126A", 7.2, "KHARE2013_BIOCHEM", "qualitative_thermal", 1, "FMT0287", "> WT"),
        ("MTDB_FMT0288_MtBfrB_E135A", "P9WNE5", seq_mtbfrb_pdb, "E126A", "E135A", 7.2, "KHARE2013_BIOCHEM", "qualitative_thermal", 1, "FMT0288", "> WT"),
        ("MTDB_FMT0294_HuHF_9LUW_7mut", "P02794", seq_huhf, "C91A/C103A/H106A/E117A/T123H/C131A/L166H", "C90A/C102A/H105A/E116A/T122H/C130A/L165H", 7.0, "2025_JAFC_40377399", "qualitative_thermal", 1, "FMT0294", "> WT"),
        # Destabilizing (binary_direction = 0)
        ("MTDB_FMT0167_RcBFR_E128A", "Q59738", seq_rcbfr, "E128A", None, 7.2, "10.1110/ps.0301903", "heat_clarification_55C", 0, "FMT0167", "dimer clarifies at 55 C < WT"),
        ("MTDB_FMT0168_RcBFR_E135A", "Q59738", seq_rcbfr, "E135A", None, 7.2, "10.1110/ps.0301903", "heat_clarification_55C", 0, "FMT0168", "dimer clarifies at 55 C < WT"),
        ("MTDB_FMT0223_BfMFtn_L134P", "P07798", seq_bfmftn, "L134P", None, 7.0, "TAKAGI2003_PNAS", "pore_CD_and_A280", 0, "FMT0223", "<43 C pore / <85 C global"),
        ("MTDB_FMT0290_BfMFtn_D80K_I81K", "P07798", seq_bfmftn, "D80K/I81K", None, 7.0, "2014_ACSCHEMBIOL_25148224", "qualitative_thermal", 0, "FMT0290", "< WT & inhibited assembly"),
        ("MTDB_FMT0291_HuFTL_A96T", "P02792", seq_huftl, "A96T", None, 7.0, "2019_BIOCHEM_30986045", "qualitative_thermal", 0, "FMT0291", "< HuFTL WT"),
        ("MTDB_FMT0270_PfFtn_R117A", "Q8U2T8", seq_pfftn_pdb, "R117A", "R117A", 7.0, "2015_SCIREP_PMC4646333", "heat_100C_timecourse", 0, "FMT0270", "< WT-PfFtn"),
    ]
    for cid, uid, wseq, mc, po, ph_v, stdy, ep, bdir, srid, note in s5_specs:
        add_record(
            case_id=cid, cohort="S5_TierBC_Directional", subcohort="S5_Unseen13",
            uniprot_id=uid, wt_seq=wseq, mut_canon=mc, pdb_mut_override=po,
            delta_tm=None, ph=ph_v, study_id=stdy, endpoint_type=ep,
            binary_direction=bdir, source_record_id=srid, notes=note,
        )

    return records, templates_meta


def main() -> None:
    records, templates_meta = build_curated_cohorts()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    tpl_path = DATA_DIR / "wt_templates.json"
    tpl_path.write_text(json.dumps(templates_meta, indent=2), encoding="utf-8")

    csv_path = DATA_DIR / "curated_ferritin_cohorts.csv"
    fieldnames = list(records[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)

    from collections import Counter
    print(f"Wrote {len(records)} curated cases -> {csv_path}")
    print("Cohort counts:", Counter(r["cohort"] for r in records))
    print("Subcohort counts:", Counter(r["subcohort"] for r in records))


if __name__ == "__main__":
    main()
