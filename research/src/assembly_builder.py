"""Biological assembly reconstruction and symmetry-interface feature extraction for ferritins."""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import Atom, Chain, Model, PDBIO, PDBParser, Residue, Structure
from Bio.SeqUtils import seq1

AA_STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

FORMAL_CHARGE = {"ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0, "HIS": 0.1}
AA1_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHAIN_IDS_POOL = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


@dataclass
class SymmetryNeighbor:
    """Descriptor for one symmetry-related neighbor subunit in a biological assembly."""
    new_chain_id: str
    source_chain_id: str
    op_index: int
    min_heavy_dist_a: float
    n_contacts_6a: int
    rotation_deg: float | None
    symmetry_class: str  # C2, C3, C4, OTHER


def parse_biomolecule_transforms(pdb_path: Path, target_chain_id: str) -> list[tuple[list[str], list[np.ndarray]]]:
    """Parse REMARK 350 BIOMT transform groups for the biomolecule containing target_chain_id."""
    lines = pdb_path.read_text(errors="ignore").splitlines()
    biomols: list[dict] = []
    cur: dict | None = None
    cur_group: dict | None = None
    for ln in lines:
        if ln.startswith("REMARK 350 BIOMOLECULE:"):
            cur = {"id": ln.split(":")[-1].strip(), "groups": []}
            biomols.append(cur)
            cur_group = None
        elif cur is not None and ln.startswith("REMARK 350 APPLY THE FOLLOWING TO CHAINS:"):
            chains = [c.strip() for c in ln.split(":", 1)[1].split(",") if c.strip()]
            cur_group = {"chains": chains, "mats": {}}
            cur["groups"].append(cur_group)
        elif cur_group is not None and ln.startswith("REMARK 350                    AND CHAINS:"):
            cur_group["chains"].extend([c.strip() for c in ln.split(":", 1)[1].split(",") if c.strip()])
        elif cur is not None and ln.startswith("REMARK 350   BIOMT"):
            if cur_group is None:
                cur_group = {"chains": [target_chain_id], "mats": {}}
                cur["groups"].append(cur_group)
            row_idx = int(ln[18]) - 1
            op_id = int(ln[19:24].strip())
            vals = [float(x) for x in ln[24:].split()]
            cur_group["mats"].setdefault(op_id, np.zeros((3, 4), dtype=float))[row_idx] = vals
        elif ln.startswith("ATOM  "):
            break

    chosen = None
    for bm in biomols:
        if any(target_chain_id in grp["chains"] and grp["mats"] for grp in bm["groups"]):
            chosen = bm
            break
    if chosen is None and biomols:
        chosen = biomols[0]
    if chosen is None or not any(grp["mats"] for grp in chosen["groups"]):
        return [([target_chain_id], [np.hstack([np.eye(3), np.zeros((3, 1))])])]
    return [
        (grp["chains"], [grp["mats"][k] for k in sorted(grp["mats"])])
        for grp in chosen["groups"]
        if grp["mats"]
    ]


def _clean_residue(res: Residue.Residue, normalize_mse: bool = True) -> Residue.Residue | None:
    """Return a cleaned standard residue copy (MSE optionally converted to MET)."""
    rname = res.resname.strip()
    if normalize_mse and rname == "MSE":
        rname = "MET"
    if rname not in AA_STANDARD:
        return None
    if not all(a in res for a in ("N", "CA", "C", "O")):
        return None
    new_res = Residue.Residue((" ", int(res.id[1]), res.id[2]), rname, res.segid)
    for atom in res:
        elem = atom.element.strip().upper()
        if elem == "H" or atom.name.strip().startswith("H"):
            continue
        aname = atom.name.strip()
        if normalize_mse and res.resname.strip() == "MSE" and aname == "SE":
            aname = "SD"
            elem = "S"
        new_atom = Atom.Atom(
            aname,
            atom.coord.copy(),
            atom.bfactor,
            atom.occupancy,
            " ",
            f" {aname:<3s}" if len(aname) < 4 else aname,
            atom.serial_number,
            elem or aname[0],
        )
        new_res.add(new_atom)
    return new_res


def _build_chain_ca_and_seq(chain: Chain.Chain) -> tuple[list[Residue.Residue], str, np.ndarray]:
    res_list = [r for r in chain if r.id[0] == " " and "CA" in r]
    seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in res_list)
    cas = np.asarray([r["CA"].coord for r in res_list], dtype=float)
    return res_list, seq, cas


def _rotation_angle_aligned(seq_a: str, ca_a: np.ndarray, seq_b: str, ca_b: np.ndarray) -> float | None:
    """Compute Kabsch rotation angle (degrees) between sequence-aligned CAs of two subunits."""
    if len(seq_a) < 20 or len(seq_b) < 20:
        return None
    if seq_a == seq_b and len(ca_a) == len(ca_b):
        pts_a, pts_b = ca_a, ca_b
    else:
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
        aln = aligner.align(seq_a, seq_b)[0]
        idx_a, idx_b = [], []
        for (qa, qb), (ta, tb) in zip(*aln.aligned):
            for k in range(qb - qa):
                idx_a.append(qa + k)
                idx_b.append(ta + k)
        if len(idx_a) < 30:
            return None
        pts_a, pts_b = ca_a[idx_a], ca_b[idx_b]

    A = pts_a - pts_a.mean(axis=0)
    B = pts_b - pts_b.mean(axis=0)
    U, _, Vt = np.linalg.svd(B.T @ A)
    rot = U @ Vt
    if np.linalg.det(rot) < 0:
        U[:, -1] *= -1
        rot = U @ Vt
    cos_t = float(np.clip((np.trace(rot) - 1.0) / 2.0, -1.0, 1.0))
    return math.degrees(math.acos(cos_t))


def _classify_symmetry(angle_deg: float | None) -> str:
    if angle_deg is None:
        return "OTHER"
    if abs(angle_deg - 180.0) <= 20.0:
        return "C2"
    if abs(angle_deg - 120.0) <= 20.0:
        return "C3"
    if abs(angle_deg - 90.0) <= 20.0 or abs(angle_deg - 48.0) <= 15.0:
        return "C4"
    return "OTHER"


def _is_identity(mat: np.ndarray) -> bool:
    return bool(np.allclose(mat[:, :3], np.eye(3), atol=1e-3) and np.allclose(mat[:, 3], 0.0, atol=1e-3))


def build_biological_assembly_and_subcomplex(
    pdb_path: Path,
    *,
    target_chain_id: str = "A",
    biomt_source_pdb: Path | None = None,
    out_monomer_pdb: Path | None = None,
    out_subcomplex_pdb: Path | None = None,
    out_full_assembly_pdb: Path | None = None,
    contact_cutoff_a: float = 10.0,
    normalize_mse: bool = True,
) -> dict:
    """Reconstruct full biological assembly and extract minimal touching-interface subcomplex.

    Chain 'A' in both monomer and subcomplex outputs is always the unmodified target subunit frame.
    """
    parser = PDBParser(QUIET=True)
    raw_struct = parser.get_structure(pdb_path.stem, str(pdb_path))
    raw_model = raw_struct[0]

    bm_source = biomt_source_pdb or pdb_path
    bm_groups = parse_biomolecule_transforms(bm_source, target_chain_id)
    all_bm_chains: list[str] = []
    for grp_chains, _ in bm_groups:
        for cid in grp_chains:
            if cid not in all_bm_chains:
                all_bm_chains.append(cid)
    if target_chain_id not in all_bm_chains:
        all_bm_chains = [target_chain_id] + all_bm_chains

    clean_templates: dict[str, list[Residue.Residue]] = {}
    for cid in all_bm_chains:
        if cid not in raw_model:
            continue
        cleaned = []
        for res in raw_model[cid]:
            if res.id[0] != " ":
                continue
            cr = _clean_residue(res, normalize_mse=normalize_mse)
            if cr is not None:
                cleaned.append(cr)
        if cleaned:
            clean_templates[cid] = cleaned

    if target_chain_id not in clean_templates:
        raise ValueError(f"Target chain {target_chain_id} has no standard complete residues in {pdb_path}")

    # Collect raw heteroatom cofactor coords (metals/heme, excluding HOH/buffers) for geometric annotation
    cofactor_coords: list[np.ndarray] = []
    ignored_het = {"HOH", "DOD", "SO4", "PO4", "ACT", "GOL", "EDO", "PEG", "PGE", "MPD", "FMT", "DMS", "CL", "NA"}
    for chain in raw_model:
        for res in chain:
            het = res.id[0]
            rname = res.resname.strip()
            if het.startswith("H_") or het == "W":
                if rname in ignored_het:
                    continue
                for atom in res:
                    if not atom.element.strip().upper().startswith("H"):
                        cofactor_coords.append(atom.coord.copy())
    cofactor_arr = np.asarray(cofactor_coords, dtype=float) if cofactor_coords else np.zeros((0, 3), dtype=float)

    # Expand (mat, src_cid) pairs per group; place target_chain_id's identity (or first) transform at index 0
    raw_pairs: list[tuple[np.ndarray, str, int]] = []
    op_counter = 0
    for grp_chains, grp_mats in bm_groups:
        id_mats = [m for m in grp_mats if _is_identity(m)]
        other_mats = [m for m in grp_mats if not _is_identity(m)]
        for mat in id_mats + other_mats:
            for cid in grp_chains:
                if cid in clean_templates:
                    raw_pairs.append((mat, cid, op_counter))
            op_counter += 1

    primary_idx = next((i for i, (_, cid, _) in enumerate(raw_pairs) if cid == target_chain_id), 0)
    ordered_pairs = [raw_pairs[primary_idx]] + [p for i, p in enumerate(raw_pairs) if i != primary_idx]

    full_model = Model.Model(0)
    subunits_data: list[dict] = []
    seen_centroids: list[np.ndarray] = []

    for sub_idx, (mat, src_cid, op_idx) in enumerate(ordered_pairs):
        if len(subunits_data) >= len(CHAIN_IDS_POOL):
            break
        R, t = mat[:, :3], mat[:, 3]

        # Build transformed chain and check against exact duplicate placement
        new_cid = CHAIN_IDS_POOL[len(subunits_data)]
        new_chain = Chain.Chain(new_cid)
        heavy_coords = []
        for res_tpl in clean_templates[src_cid]:
            r_copy = copy.deepcopy(res_tpl)
            for atom in r_copy:
                atom.set_coord(R @ atom.coord + t)
                heavy_coords.append(atom.coord)
            new_chain.add(r_copy)

        h_arr = np.asarray(heavy_coords, dtype=float)
        centroid = h_arr.mean(axis=0)
        if any(np.linalg.norm(centroid - prev) < 2.0 for prev in seen_centroids):
            continue
        seen_centroids.append(centroid)
        full_model.add(new_chain)
        _, seq_c, cas_c = _build_chain_ca_and_seq(new_chain)
        subunits_data.append(
            {
                "new_cid": new_cid,
                "src_cid": src_cid,
                "op_idx": op_idx,
                "chain": new_chain,
                "heavy": h_arr,
                "seq": seq_c,
                "cas": cas_c,
            }
        )

    ref = subunits_data[0]
    neighbors: list[SymmetryNeighbor] = []
    touching_cids: list[str] = [ref["new_cid"]]

    for other in subunits_data[1:]:
        ca_min = float(np.linalg.norm(ref["cas"][:, None, :] - other["cas"][None, :, :], axis=-1).min())
        if ca_min > contact_cutoff_a + 12.0:
            continue
        d2 = np.sum((ref["heavy"][:, None, :] - other["heavy"][None, :, :]) ** 2, axis=-1)
        min_d = float(math.sqrt(float(d2.min())))
        if min_d <= contact_cutoff_a:
            n_6a = int(np.sum(np.any(d2 <= 36.0, axis=1)))
            ang = _rotation_angle_aligned(ref["seq"], ref["cas"], other["seq"], other["cas"])
            sym = _classify_symmetry(ang)
            neighbors.append(
                SymmetryNeighbor(
                    new_chain_id=other["new_cid"],
                    source_chain_id=other["src_cid"],
                    op_index=other["op_idx"],
                    min_heavy_dist_a=min_d,
                    n_contacts_6a=n_6a,
                    rotation_deg=ang,
                    symmetry_class=sym,
                )
            )
            touching_cids.append(other["new_cid"])

    io = PDBIO()
    if out_full_assembly_pdb is not None:
        out_full_assembly_pdb.parent.mkdir(parents=True, exist_ok=True)
        fs = Structure.Structure("asm")
        fs.add(full_model)
        io.set_structure(fs)
        io.save(str(out_full_assembly_pdb))

    if out_monomer_pdb is not None:
        out_monomer_pdb.parent.mkdir(parents=True, exist_ok=True)
        ms = Structure.Structure("mono")
        mm = Model.Model(0)
        mm.add(copy.deepcopy(ref["chain"]))
        ms.add(mm)
        io.set_structure(ms)
        io.save(str(out_monomer_pdb))

    if out_subcomplex_pdb is not None:
        out_subcomplex_pdb.parent.mkdir(parents=True, exist_ok=True)
        ss = Structure.Structure("sub")
        sm = Model.Model(0)
        for cid in touching_cids:
            sm.add(copy.deepcopy(full_model[cid]))
        ss.add(sm)
        io.set_structure(ss)
        io.save(str(out_subcomplex_pdb))

    return {
        "total_subunits": len(subunits_data),
        "subcomplex_chains": touching_cids,
        "neighbors": [n.__dict__ for n in neighbors],
        "ref_sequence": ref["seq"],
        "cofactor_coords": cofactor_arr,
        "subcomplex_model": {d["new_cid"]: d["chain"] for d in subunits_data if d["new_cid"] in touching_cids},
        "ref_chain_obj": ref["chain"],
    }


def extract_site_interface_features(
    assembly_info: dict,
    resnum_list: list[int],
    mutations: list[tuple[int, str, str]],
) -> dict[str, float]:
    """Compute 3D monomer vs oligomer interface, symmetry, burial, and electrostatic shifts.

    Args:
        assembly_info: output dict from build_biological_assembly_and_subcomplex.
        resnum_list: PDB residue numbers on Chain A corresponding to the mutation sites (1-based).
        mutations: list of (pos1_seq, wt_aa1, mut_aa1) tuples.
    """
    ref_chain: Chain.Chain = assembly_info["ref_chain_obj"]
    sub_chains: dict[str, Chain.Chain] = assembly_info["subcomplex_model"]
    neighbors_meta: dict[str, dict] = {n["new_chain_id"]: n for n in assembly_info["neighbors"]}
    cofactor_arr: np.ndarray = assembly_info["cofactor_coords"]

    ref_res_map = {int(r.id[1]): r for r in ref_chain if r.id[0] == " "}

    # Pre-collect CA and heavy-atom coordinates of same chain vs neighbor chains
    intra_cas = []
    intra_atoms = []
    for r in ref_chain:
        if r.id[0] != " ":
            continue
        if "CA" in r:
            intra_cas.append((int(r.id[1]), r["CA"].coord))
        for a in r:
            intra_atoms.append((int(r.id[1]), a.coord, FORMAL_CHARGE.get(r.resname.strip(), 0.0)))

    partner_cas = []
    partner_atoms = []
    for cid, ch in sub_chains.items():
        if cid == "A":
            continue
        sym_cls = neighbors_meta[cid]["symmetry_class"]
        for r in ch:
            if r.id[0] != " ":
                continue
            if "CA" in r:
                partner_cas.append((sym_cls, r["CA"].coord))
            q = FORMAL_CHARGE.get(r.resname.strip(), 0.0)
            for a in r:
                # Place charge only on sidechain charged headgroups / terminal atom approximation, or uniform split:
                # Using CZ/NZ/CG/CD/NE2-weighted proxy cleanly via charged residue functional atoms
                an = a.name.strip()
                atom_q = 0.0
                rname = r.resname.strip()
                if rname == "ASP" and an in ("OD1", "OD2"):
                    atom_q = -0.5
                elif rname == "GLU" and an in ("OE1", "OE2"):
                    atom_q = -0.5
                elif rname == "LYS" and an == "NZ":
                    atom_q = 1.0
                elif rname == "ARG" and an in ("NH1", "NH2"):
                    atom_q = 0.5
                elif rname == "HIS" and an in ("ND1", "NE2"):
                    atom_q = 0.05
                partner_atoms.append((cid, sym_cls, a.coord, atom_q))

    partner_ca_coords = np.asarray([x[1] for x in partner_cas], dtype=float) if partner_cas else np.zeros((0, 3))
    partner_atom_coords = np.asarray([x[2] for x in partner_atoms], dtype=float) if partner_atoms else np.zeros((0, 3))
    partner_atom_charges = np.asarray([x[3] for x in partner_atoms], dtype=float) if partner_atoms else np.zeros((0,))

    site_features = []
    for resnum, (_, wt_aa, mut_aa) in zip(resnum_list, mutations):
        res = ref_res_map[resnum]
        res_atoms = np.asarray([a.coord for a in res], dtype=float)
        site_ca = res["CA"].coord

        # Intra-chain neighbor CAs (excluding self) within 10 A
        intra_ca_other = np.asarray([c for rn, c in intra_cas if rn != resnum], dtype=float)
        n_intra_ca10 = int(np.sum(np.linalg.norm(intra_ca_other - site_ca, axis=-1) <= 10.0)) if len(intra_ca_other) else 0
        n_inter_ca10 = int(np.sum(np.linalg.norm(partner_ca_coords - site_ca, axis=-1) <= 10.0)) if len(partner_ca_coords) else 0

        if len(partner_atom_coords) > 0:
            d_heavy = np.sqrt(np.sum((res_atoms[:, None, :] - partner_atom_coords[None, :, :]) ** 2, axis=-1))
            min_per_partner_atom = d_heavy.min(axis=0)
            min_inter_dist = float(min_per_partner_atom.min())
            c5 = int(np.sum(min_per_partner_atom <= 5.0))
            c8 = int(np.sum(min_per_partner_atom <= 8.0))

            # Contacts <= 6 A by symmetry class
            sym_counts = {"C2": 0, "C3": 0, "C4": 0, "OTHER": 0}
            touch_chains = set()
            for idx_pa in np.where(min_per_partner_atom <= 6.0)[0]:
                cid, sym_cls, _, _ = partner_atoms[idx_pa]
                sym_counts[sym_cls] = sym_counts.get(sym_cls, 0) + 1
                touch_chains.add(cid)

            # Neighboring electrostatic field around mutated residue centroid
            site_center = res_atoms.mean(axis=0)
            d_center = np.linalg.norm(partner_atom_coords - site_center, axis=-1)
            mask_8a = (d_center <= 8.0) & (np.abs(partner_atom_charges) > 0)
            env_inter_q = float(np.sum(partner_atom_charges[mask_8a] / np.maximum(d_center[mask_8a], 2.5)))
        else:
            min_inter_dist = 50.0
            c5 = 0
            c8 = 0
            sym_counts = {"C2": 0, "C3": 0, "C4": 0, "OTHER": 0}
            touch_chains = set()
            env_inter_q = 0.0

        # Intra-chain electrostatic field around mutated residue
        intra_env_q = 0.0
        site_center = res_atoms.mean(axis=0)
        for r in ref_chain:
            if r.id[0] != " " or int(r.id[1]) == resnum:
                continue
            q = FORMAL_CHARGE.get(r.resname.strip(), 0.0)
            if abs(q) == 0.0 or "CA" not in r:
                continue
            dist_c = float(np.linalg.norm(r["CA"].coord - site_center))
            if dist_c <= 9.0:
                intra_env_q += q / max(dist_c, 3.0)

        q_wt = AA1_CHARGE.get(wt_aa, 0.0)
        q_mut = AA1_CHARGE.get(mut_aa, 0.0)
        dq = q_mut - q_wt
        delta_elec_inter = dq * env_inter_q
        delta_elec_intra = dq * intra_env_q
        delta_hydro = HYDROPATHY[mut_aa] - HYDROPATHY[wt_aa]

        cof_dist = float(np.linalg.norm(cofactor_arr - site_center, axis=-1).min()) if len(cofactor_arr) else 50.0
        total_ca10 = n_intra_ca10 + n_inter_ca10
        inter_frac = float(n_inter_ca10 / total_ca10) if total_ca10 > 0 else 0.0

        site_features.append(
            {
                "n_intra_ca10": float(n_intra_ca10),
                "n_inter_ca10": float(n_inter_ca10),
                "inter_ca10_fraction": inter_frac,
                "min_inter_heavy_dist_a": min_inter_dist,
                "inter_contacts_5a": float(c5),
                "inter_contacts_8a": float(c8),
                "n_touching_partner_chains_6a": float(len(touch_chains)),
                "c2_contacts_6a": float(sym_counts["C2"]),
                "c3_contacts_6a": float(sym_counts["C3"]),
                "c4_contacts_6a": float(sym_counts["C4"]),
                "is_interface_site": float(min_inter_dist <= 5.5),
                "is_c2_site": float(sym_counts["C2"] > 0 and sym_counts["C2"] >= max(sym_counts["C3"], sym_counts["C4"])),
                "is_c3_site": float(sym_counts["C3"] > 0 and sym_counts["C3"] >= max(sym_counts["C2"], sym_counts["C4"])),
                "is_c4_site": float(sym_counts["C4"] > 0 and sym_counts["C4"] >= max(sym_counts["C2"], sym_counts["C3"])),
                "delta_charge": dq,
                "delta_abs_charge": abs(q_mut) - abs(q_wt),
                "delta_hydropathy": delta_hydro,
                "env_inter_q": env_inter_q,
                "delta_elec_inter_energy_proxy": delta_elec_inter,
                "delta_elec_intra_energy_proxy": delta_elec_intra,
                "min_cofactor_dist_a": cof_dist,
                "near_cofactor_8a": float(cof_dist <= 8.0),
            }
        )

    # Aggregate cleanly over single or multi-point mutation sites
    n_sites = float(len(site_features))
    agg = {"n_mutations": n_sites, "total_subunits": float(assembly_info["total_subunits"])}
    sum_keys = {
        "n_inter_ca10", "inter_contacts_5a", "inter_contacts_8a", "c2_contacts_6a", "c3_contacts_6a",
        "c4_contacts_6a", "delta_charge", "delta_abs_charge", "delta_hydropathy",
        "delta_elec_inter_energy_proxy", "delta_elec_intra_energy_proxy",
    }
    max_keys = {
        "inter_ca10_fraction", "n_touching_partner_chains_6a", "is_interface_site",
        "is_c2_site", "is_c3_site", "is_c4_site", "near_cofactor_8a",
    }
    min_keys = {"min_inter_heavy_dist_a", "min_cofactor_dist_a"}
    for key in site_features[0]:
        vals = [sf[key] for sf in site_features]
        if key in min_keys:
            agg[key] = float(min(vals))
        elif key in max_keys:
            agg[key] = float(max(vals))
        elif key in sum_keys:
            agg[key] = float(sum(vals))
        else:
            agg[key] = float(sum(vals) / len(vals))
    return agg
