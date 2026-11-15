"""Assembly-GeoDTm-3D & ESM-2 feature adapter for monomeric and oligomeric ferritins."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from typing import Any

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import numpy as np
import torch

_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from paths import DATA_DIR, GEODTM_CKPT, GEOSTAB_ROOT, ESM2_WEIGHTS, RESEARCH_ROOT

if str(GEOSTAB_ROOT / "generate_features") not in sys.path:
    sys.path.insert(0, str(GEOSTAB_ROOT / "generate_features"))
if str(GEOSTAB_ROOT / "model_dTm_3D") not in sys.path:
    sys.path.insert(0, str(GEOSTAB_ROOT / "model_dTm_3D"))

import coordinate as coord_mod  # noqa: E402
import esm  # noqa: E402
import pair as pair_mod  # noqa: E402

CACHE_ROOT = RESEARCH_ROOT / "cache"
ESM2_CACHE_DIR = CACHE_ROOT / "esm2_seq"
GEODTM_FEATURE_CACHE_DIR = CACHE_ROOT / "geodtm_cases"

AA20 = "ACDEFGHIKLMNPQRSTVWY"
HIS_ALIASES = {"H1S", "H2S", "HSD", "HSE", "HSP"}

_ESM_BUNDLE: tuple[Any, Any, list[int]] | None = None
_GEODTM_MODEL: Any | None = None
_AA_VEC_TABLE: dict[int, torch.Tensor] | None = None
_COORD_MEM_CACHE: dict[str, list[dict[str, Any]]] = {}


def _get_aa_vec_table() -> dict[int, torch.Tensor]:
    global _AA_VEC_TABLE
    if _AA_VEC_TABLE is None:
        raw = json.loads((DATA_DIR / "fixed_aa_embedding_20x7.json").read_text(encoding="utf-8"))
        _AA_VEC_TABLE = {
            int(k): torch.tensor(v, dtype=torch.float32) for k, v in raw["aa_index_to_vec"].items()
        }
    return _AA_VEC_TABLE


def _get_geodtm_model() -> Any:
    global _GEODTM_MODEL
    if _GEODTM_MODEL is None:
        torch.set_num_threads(8)
        model = torch.load(GEODTM_CKPT, map_location="cpu", weights_only=False)
        model.eval()
        _GEODTM_MODEL = model
    return _GEODTM_MODEL


def _get_esm_bundle() -> tuple[Any, Any, list[int]]:
    global _ESM_BUNDLE
    if _ESM_BUNDLE is None:
        torch.set_num_threads(8)
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(ESM2_WEIGHTS))
        model.eval()
        aa20_ids = [alphabet.get_idx(aa) for aa in AA20]
        _ESM_BUNDLE = (model, alphabet, aa20_ids)
    return _ESM_BUNDLE


def get_esm2_sequence_features(seq: str) -> dict[str, torch.Tensor]:
    """Compute/load cached ESM-2 layer-33 representations and AA20 log-probabilities."""
    ESM2_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    seq_key = hashlib.sha256(seq.encode("utf-8")).hexdigest()[:16]
    cache_path = ESM2_CACHE_DIR / f"{seq_key}.pt"
    if cache_path.exists():
        data = torch.load(cache_path, map_location="cpu", weights_only=False)
        if data.get("seq") == seq and "repr33" in data and "log_probs20" in data:
            return data

    model, alphabet, aa20_ids = _get_esm_bundle()
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("protein", seq)])
    with torch.no_grad():
        res = model(tokens, repr_layers=[33], return_contacts=False)
        repr33 = res["representations"][33][0, 1:-1, :].to(torch.float32).cpu()
        logits_full = res["logits"][0, 1:-1, :]
        log_probs_full = torch.log_softmax(logits_full, dim=-1)
        log_probs20 = log_probs_full[:, aa20_ids].to(torch.float32).cpu()

    data = {"seq": seq, "repr33": repr33, "log_probs20": log_probs20}
    torch.save(data, cache_path)
    return data


def compute_esm2_evolutionary_scores(
    wt_seq: str,
    mutations_canon: list[tuple[int, str, str]],
) -> dict[str, float]:
    """Compute additive WT-marginal LLR, site probabilities, and Shannon entropy."""
    wt_esm = get_esm2_sequence_features(wt_seq)
    lp20 = wt_esm["log_probs20"]
    llrs: list[float] = []
    wt_probs: list[float] = []
    mut_probs: list[float] = []
    entropies: list[float] = []
    for pos1, wt_aa, mut_aa in mutations_canon:
        idx0 = pos1 - 1
        assert wt_seq[idx0] == wt_aa, f"Expected {wt_aa} at {pos1}, got {wt_seq[idx0]}"
        w_i = AA20.index(wt_aa)
        m_i = AA20.index(mut_aa)
        row_lp = lp20[idx0]
        llrs.append(float((row_lp[m_i] - row_lp[w_i]).item()))
        probs = torch.exp(row_lp)
        wt_probs.append(float(probs[w_i].item()))
        mut_probs.append(float(probs[m_i].item()))
        # Normalized entropy over the 20 standard amino acids
        p_norm = probs / probs.sum().clamp_min(1e-12)
        entropies.append(float((-(p_norm * torch.log(p_norm.clamp_min(1e-12))).sum()).item()))
    return {
        "esm2_wt_marginal_llr": float(sum(llrs)),
        "esm2_mean_llr": float(np.mean(llrs)),
        "esm2_min_llr": float(min(llrs)),
        "esm2_site_wt_prob": float(np.mean(wt_probs)),
        "esm2_site_mut_prob": float(np.mean(mut_probs)),
        "esm2_site_entropy": float(np.mean(entropies)),
    }


def _build_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.query_end_gap_score = 0.0
    aligner.target_end_gap_score = 0.0
    return aligner


def _align_chain_to_sequence(chain_seq: str, full_seq: str) -> dict[int, int]:
    """Map 0-based chain residue index -> 0-based sequence index in full_seq."""
    if chain_seq == full_seq:
        return {i: i for i in range(len(chain_seq))}
    idx_direct = full_seq.find(chain_seq)
    if idx_direct >= 0:
        return {i: idx_direct + i for i in range(len(chain_seq))}
    aligner = _build_aligner()
    aln = aligner.align(full_seq, chain_seq)[0]
    q_blocks, t_blocks = aln.aligned
    mapping: dict[int, int] = {}
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        length = min(qe - qs, te - ts)
        for k in range(length):
            mapping[ts + k] = qs + k
    return mapping


def _split_pdb_by_chain_clean(pdb_path: Path, tmp_root: Path) -> list[tuple[str, Path]]:
    """Split PDB into single-chain files normalizing Histidine aliases and MSE."""
    lines_by_chain: dict[str, list[str]] = {}
    order: list[str] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            resname = line[17:20].strip()
            if rec == "HETATM" and resname == "MSE":
                line = "ATOM  " + line[6:17] + "MET" + line[20:]
                resname = "MET"
            elif rec == "HETATM":
                continue
            if resname in HIS_ALIASES:
                line = line[:17] + "HIS" + line[20:]
            cid = line[21]
            if cid not in lines_by_chain:
                lines_by_chain[cid] = []
                order.append(cid)
            lines_by_chain[cid].append(line)

    out: list[tuple[str, Path]] = []
    for cid in order:
        cpath = tmp_root / f"chain_{cid}.pdb"
        cpath.write_text("".join(lines_by_chain[cid]) + "END\n", encoding="utf-8")
        out.append((cid, cpath))
    return out


def _extract_chain_aa_seq(chain_pdb: Path) -> str:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(chain_pdb.stem, str(chain_pdb))
    chain = next(struct[0].get_chains())
    return "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in chain if r.id[0] == " ")


def _extract_pdb_coords_by_chain(pdb_path: Path) -> list[dict[str, Any]]:
    """Extract coordinate.pt dicts per chain (cached in-process by resolved path)."""
    cache_key = str(pdb_path.resolve())
    if cache_key in _COORD_MEM_CACHE:
        return _COORD_MEM_CACHE[cache_key]

    chains_out: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="geodtm_coord_") as tmp_str:
        tmp_root = Path(tmp_str)
        for cid, cpath in _split_pdb_by_chain_clean(pdb_path, tmp_root):
            cdir = tmp_root / f"out_{cid}"
            cdir.mkdir(parents=True, exist_ok=True)
            coord_mod.main.callback(str(cpath), str(cdir))
            coord = torch.load(cdir / "coordinate.pt", map_location="cpu", weights_only=False)
            cseq = _extract_chain_aa_seq(cpath)
            if len(cseq) != len(coord["resseq"]):
                raise RuntimeError(
                    f"{pdb_path} chain {cid}: seq len {len(cseq)} != coord len {len(coord['resseq'])}"
                )
            chains_out.append({"chain_id": cid, "coord": coord, "chain_seq": cseq})
    _COORD_MEM_CACHE[cache_key] = chains_out
    return chains_out


def _build_node_pool(
    pdb_path: Path,
    *,
    target_chain_id: str,
    full_seq: str,
    esm_repr: torch.Tensor,
    mut_seq_pos0: list[int],
    wt_seq: str | None = None,
    wt_esm_repr: torch.Tensor | None = None,
) -> dict[str, Any]:
    """Assemble concatenated multi-chain coordinate + ESM2 pool with Chain A first."""
    chains = _extract_pdb_coords_by_chain(pdb_path)
    # Put target_chain_id first so node indices 0..N_A-1 correspond to central subunit
    ordered = sorted(chains, key=lambda c: (0 if c["chain_id"] == target_chain_id else 1))
    if not ordered or ordered[0]["chain_id"] != target_chain_id:
        raise RuntimeError(f"Chain {target_chain_id} not found in {pdb_path}")

    aa_vec_table = _get_aa_vec_table()
    pos14_list: list[torch.Tensor] = []
    mask_list: list[torch.Tensor] = []
    dyn_list: list[torch.Tensor] = []
    fix7_list: list[torch.Tensor] = []
    resseq_list: list[torch.Tensor] = []
    aa_idx_list: list[torch.Tensor] = []
    is_target_chain_flags: list[bool] = []
    seq_pos0_all: list[int] = []

    for centry in ordered:
        coord = centry["coord"]
        cseq = centry["chain_seq"]
        use_seq = full_seq
        use_repr = esm_repr
        # Fallback to WT sequence alignment only if a crystal neighbor chain wasn't mutated at terminal pos
        c_to_seq = _align_chain_to_sequence(cseq, use_seq)
        mismatches = sum(
            1 for ci, si in c_to_seq.items() if cseq[ci] != use_seq[si]
        )
        if mismatches > 0 and wt_seq is not None and wt_esm_repr is not None:
            alt_map = _align_chain_to_sequence(cseq, wt_seq)
            alt_mism = sum(1 for ci, si in alt_map.items() if cseq[ci] != wt_seq[si])
            if alt_mism < mismatches:
                c_to_seq = alt_map
                use_seq = wt_seq
                use_repr = wt_esm_repr

        n_res = len(coord["resseq"])
        if len(c_to_seq) != n_res:
            raise RuntimeError(
                f"{pdb_path} chain {centry['chain_id']}: aligned {len(c_to_seq)}/{n_res} residues"
            )

        pos14_list.append(coord["pos14"].to(torch.float32))
        mask_list.append(coord["pos14_mask"].to(torch.bool))
        resseq_list.append(coord["resseq"])
        aa_idx_list.append(coord["aa"])

        seq_idxs = [c_to_seq[i] for i in range(n_res)]
        dyn_list.append(use_repr[seq_idxs])
        fix7 = torch.stack([aa_vec_table[int(a)] for a in coord["aa"]], dim=0)
        fix7_list.append(fix7)
        is_target_chain_flags.extend([centry["chain_id"] == target_chain_id] * n_res)
        seq_pos0_all.extend(seq_idxs)

    n_target = len(ordered[0]["coord"]["resseq"])
    target_seq_to_node = {seq_pos0_all[i]: i for i in range(n_target)}
    mut_node_indices = [target_seq_to_node[p0] for p0 in mut_seq_pos0]

    return {
        "pos14": torch.cat(pos14_list, dim=0),
        "pos14_mask": torch.cat(mask_list, dim=0),
        "aa": torch.cat(aa_idx_list, dim=0),
        "resseq": torch.cat(resseq_list, dim=0),
        "dynamic": torch.cat(dyn_list, dim=0),
        "fixed7": torch.cat(fix7_list, dim=0),
        "is_target_chain": torch.tensor(is_target_chain_flags, dtype=torch.bool),
        "mut_node_indices": mut_node_indices,
        "n_target_nodes": n_target,
    }


def _slice_geodtm_input(
    pool: dict[str, Any],
    *,
    k: int,
    ph: float,
) -> tuple[dict[str, torch.Tensor], float]:
    """Select k-nearest CA nodes around mutated site(s) and build GeoDTm-3D tensors."""
    n_total = pool["pos14"].shape[0]
    k_eff = min(int(k), n_total)
    ca = pool["pos14"][:, 1, :]
    mut_idx = torch.tensor(pool["mut_node_indices"], dtype=torch.long)
    dists = torch.linalg.norm(ca[:, None, :] - ca[mut_idx][None, :, :], dim=-1).min(dim=-1).values
    top_idx = torch.topk(dists, k=k_eff, largest=False).indices.sort().values

    pos14_k = pool["pos14"][top_idx]
    mask_k = pool["pos14_mask"][top_idx]
    sub_coord = {
        "pos14": pos14_k,
        "pos14_mask": mask_k,
        "aa": pool["aa"][top_idx],
        "resseq": pool["resseq"][top_idx],
        "icode": " " * k_eff,
        "seq": pool["aa"][top_idx],
    }
    with tempfile.TemporaryDirectory(prefix="geodtm_pair_") as tmp_str:
        cf = Path(tmp_str) / "coordinate.pt"
        torch.save(sub_coord, cf)
        pair_mod.main.callback(str(cf), tmp_str)
        pair_k = torch.load(Path(tmp_str) / "pair.pt", map_location="cpu", weights_only=False)

    fixed9 = torch.cat(
        [
            pool["fixed7"][top_idx],
            torch.full((k_eff, 1), float(ph), dtype=torch.float32),
            torch.ones((k_eff, 1), dtype=torch.float32),
        ],
        dim=-1,
    )
    mut_pos = torch.isin(top_idx, mut_idx)
    if not bool(mut_pos.any()):
        raise RuntimeError("Mutated node missing from k-NN neighborhood!")

    partner_frac = float((~pool["is_target_chain"][top_idx]).float().mean().item())
    tensor_dict = {
        "dynamic_embedding": pool["dynamic"][top_idx].unsqueeze(0),
        "fixed_embedding": fixed9.unsqueeze(0),
        "pair": pair_k.to(torch.float32).unsqueeze(0),
        "pos14": pos14_k.to(torch.float32).unsqueeze(0),
        "atom_mask": mask_k[:, :, 0].bool().unsqueeze(0),
        "mut_pos": mut_pos.unsqueeze(0),
    }
    return tensor_dict, partner_frac


def predict_geodtm_pair(
    wt_pdb: Path,
    mt_pdb: Path,
    *,
    wt_seq: str,
    mut_seq: str,
    mut_seq_pos0: list[int],
    ph: float,
    target_chain_id: str = "A",
    k: int = 32,
) -> dict[str, float]:
    """Run GeoDTm-3D PretrainModel on a WT/MT PDB pair (monomer or subcomplex)."""
    wt_esm = get_esm2_sequence_features(wt_seq)["repr33"]
    mt_esm = get_esm2_sequence_features(mut_seq)["repr33"]

    wt_pool = _build_node_pool(
        wt_pdb,
        target_chain_id=target_chain_id,
        full_seq=wt_seq,
        esm_repr=wt_esm,
        mut_seq_pos0=mut_seq_pos0,
    )
    mt_pool = _build_node_pool(
        mt_pdb,
        target_chain_id=target_chain_id,
        full_seq=mut_seq,
        esm_repr=mt_esm,
        mut_seq_pos0=mut_seq_pos0,
        wt_seq=wt_seq,
        wt_esm_repr=wt_esm,
    )

    wt_inp, wt_sym_frac = _slice_geodtm_input(wt_pool, k=k, ph=ph)
    mt_inp, mt_sym_frac = _slice_geodtm_input(mt_pool, k=k, ph=ph)

    model = _get_geodtm_model()
    with torch.no_grad():
        pred_dtm = float(model(wt_inp, mt_inp).item())
        enc_wt = model.pretrain_encoder(
            wt_inp["dynamic_embedding"], wt_inp["pair"], wt_inp["atom_mask"]
        )
        enc_mt = model.pretrain_encoder(
            mt_inp["dynamic_embedding"], mt_inp["pair"], mt_inp["atom_mask"]
        )
        # L2 shift of pretrain_encoder representations at mutated node(s)
        wt_mut_vec = enc_wt[0][wt_inp["mut_pos"][0]].mean(dim=0)
        mt_mut_vec = enc_mt[0][mt_inp["mut_pos"][0]].mean(dim=0)
        shift_norm = float(torch.linalg.norm(mt_mut_vec - wt_mut_vec).item())

    return {
        "pred_dtm": pred_dtm,
        "sym_partners_frac": 0.5 * (wt_sym_frac + mt_sym_frac),
        "enc_shift_norm": shift_norm,
    }


def _parse_canon_muts(mut_canon: str) -> list[tuple[int, str, str]]:
    out = []
    for tok in mut_canon.strip().split("/"):
        tok = tok.strip()
        if tok:
            out.append((int(tok[1:-1]), tok[0], tok[-1]))
    return out


def compute_geodtm_features_for_case(
    *,
    case_key: str,
    wt_seq: str,
    mut_seq: str,
    mutation_canon: str,
    ph: float,
    mono_wt_pdb: Path,
    mono_mt_pdb: Path,
    sub_wt_pdb: Path,
    sub_mt_pdb: Path,
    force: bool = False,
) -> dict[str, float]:
    """Compute & cache Monomer & Assembly GeoDTm-3D (k=32/48, assay pH & pH7) + ESM2."""
    GEODTM_FEATURE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ph_tag = f"ph{ph:.2f}".replace(".", "p")
    cache_file = GEODTM_FEATURE_CACHE_DIR / f"{case_key}__{ph_tag}.json"
    if not force and cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    muts_canon = _parse_canon_muts(mutation_canon)
    mut_seq_pos0 = [p1 - 1 for p1, _, _ in muts_canon]

    esm_scores = compute_esm2_evolutionary_scores(wt_seq, muts_canon)
    mono_32 = predict_geodtm_pair(
        mono_wt_pdb, mono_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=ph, target_chain_id="A", k=32,
    )
    mono_32_ph7 = predict_geodtm_pair(
        mono_wt_pdb, mono_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=7.0, target_chain_id="A", k=32,
    )
    mono_48 = predict_geodtm_pair(
        mono_wt_pdb, mono_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=ph, target_chain_id="A", k=48,
    )
    asm_32 = predict_geodtm_pair(
        sub_wt_pdb, sub_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=ph, target_chain_id="A", k=32,
    )
    asm_32_ph7 = predict_geodtm_pair(
        sub_wt_pdb, sub_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=7.0, target_chain_id="A", k=32,
    )
    asm_48 = predict_geodtm_pair(
        sub_wt_pdb, sub_mt_pdb, wt_seq=wt_seq, mut_seq=mut_seq,
        mut_seq_pos0=mut_seq_pos0, ph=ph, target_chain_id="A", k=48,
    )

    feats = {
        **esm_scores,
        "geo_mono_k32": mono_32["pred_dtm"],
        "geo_mono_k32_ph7": mono_32_ph7["pred_dtm"],
        "geo_mono_k48": mono_48["pred_dtm"],
        "geo_asm_k32": asm_32["pred_dtm"],
        "geo_asm_k32_ph7": asm_32_ph7["pred_dtm"],
        "geo_asm_k48": asm_48["pred_dtm"],
        "geo_iface_delta_k32": asm_32["pred_dtm"] - mono_32["pred_dtm"],
        "geo_iface_delta_k48": asm_48["pred_dtm"] - mono_48["pred_dtm"],
        "geo_ph_shift_mono": mono_32["pred_dtm"] - mono_32_ph7["pred_dtm"],
        "geo_ph_shift_asm": asm_32["pred_dtm"] - asm_32_ph7["pred_dtm"],
        "geo_asm_sym_partners_frac_k32": asm_32["sym_partners_frac"],
        "geo_asm_sym_partners_frac_k48": asm_48["sym_partners_frac"],
        "geo_mono_enc_shift_k32": mono_32["enc_shift_norm"],
        "geo_asm_enc_shift_k32": asm_32["enc_shift_norm"],
        "geo_asm_enc_shift_k48": asm_48["enc_shift_norm"],
    }
    cache_file.write_text(json.dumps(feats, indent=2) + "\n", encoding="utf-8")
    return feats


def verify_1cts_equivalence() -> dict[str, float]:
    """Verify exact equivalence against official local GeoDTm-3D run on rcsb_1CTS_A_H274R_7."""
    case_dir = RESEARCH_ROOT / "results/s571_repro/cases/rcsb_1CTS_A_H274R_7"
    prep = json.loads((case_dir / "preparation.json").read_text(encoding="utf-8"))
    official_pred = float((case_dir / "result_dTm.txt").read_text(encoding="utf-8").strip())
    adapter_res = predict_geodtm_pair(
        case_dir / "wt_data/relaxed_repair.pdb",
        case_dir / "mut_data/relaxed_repair.pdb",
        wt_seq=prep["wt_sequence"],
        mut_seq=prep["mut_sequence"],
        mut_seq_pos0=[int(prep["mut_pos0"])],
        ph=float(prep["inputs"]["pH"]),
        target_chain_id="A",
        k=32,
    )
    diff = abs(adapter_res["pred_dtm"] - official_pred)
    return {
        "official_1cts_c": official_pred,
        "adapter_1cts_c": adapter_res["pred_dtm"],
        "abs_diff_c": diff,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Verify 1CTS & EcBfr equivalence")
    args = parser.parse_args()
    if args.verify:
        res = verify_1cts_equivalence()
        print("1CTS equivalence:", json.dumps(res, indent=2))
        assert res["abs_diff_c"] < 1e-4, f"Equivalence check failed: {res}"
        return 0
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
