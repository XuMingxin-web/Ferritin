#!/usr/bin/env python3
"""
Independent integrity, consistency, and dual-count reconciliation validator
for the Ferritin Mutant Thermostability deliverables in:
  /Users/a26683/Desktop/ferritinMTdatabase/
"""

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

WORKSPACE_DIR = Path(__file__).resolve().parents[1]

REQUIRED_FILES = [
    "database_and_literature_source_inventory.csv",
    "ferritin_mutant_thermostability_dataset.csv",
    "unique_mutants_summary.csv",
    "inaccessible_and_follow_up_sources.csv",
    "README_database_summary_report.md",
]

VALID_FAMILIES = {
    "Family_1_Canonical_24mer_Ferritin",
    "Family_2_Bacterioferritin_Bfr",
    "Family_3_Miniferritin_Dps",
    "Family_4_Engineered_Fusion_Display_Ferritin",
    "Family_5_Extended_Superfamily_And_Nanocages",
}
CORE_FAMILIES = VALID_FAMILIES - {"Family_5_Extended_Superfamily_And_Nanocages"}
VALID_EVIDENCE = {
    "A_WetLab_Quantitative",
    "B_WetLab_Thermal_Survival",
    "C_WetLab_Qualitative_Rank",
    "D_Computational_Prediction",
}
VALID_CATEGORIES = {
    "WT",
    "Point_Substitution",
    "Multi_Point_Combinatorial",
    "Deletion_Truncation_Frameshift_Swap",
    "Fusion_Display_Heteropolymer",
}
RANK_TIER = {
    "A_WetLab_Quantitative": 1,
    "B_WetLab_Thermal_Survival": 2,
    "C_WetLab_Qualitative_Rank": 3,
    "D_Computational_Prediction": 4,
}


def load_csv(path: Path) -> tuple[list[str], list[dict]]:
    text = path.read_text(encoding="utf-8")
    assert text.strip(), f"CSV is empty: {path}"
    reader = csv.DictReader(text.splitlines())
    rows = list(reader)
    assert reader.fieldnames and len(rows) > 0, f"CSV has no data rows: {path}"
    return list(reader.fieldnames), rows


def check_unique_ids(rows: list[dict], field: str, pattern: str, label: str) -> set[str]:
    ids = [r[field] for r in rows]
    assert all(re.match(pattern, x) for x in ids), f"Invalid {label} ID format found"
    dups = [k for k, v in Counter(ids).items() if v > 1]
    assert not dups, f"Duplicate {label} IDs found: {dups}"
    return set(ids)


def validate() -> int:
    errors: list[str] = []

    def fail(msg: str):
        errors.append(msg)

    # 1. File existence, UTF-8 decoding, non-empty checks
    for fname in REQUIRED_FILES:
        fpath = WORKSPACE_DIR / fname
        if not fpath.is_file():
            fail(f"Missing required deliverable: {fpath}")
            return 1
        try:
            content = fpath.read_text(encoding="utf-8")
            if not content.strip():
                fail(f"Empty deliverable file: {fpath}")
        except UnicodeDecodeError as exc:
            fail(f"Invalid UTF-8 in {fpath}: {exc}")

    if errors:
        for e in errors:
            print(f"FAIL: {e}")
        return 1

    _, inv_rows = load_csv(WORKSPACE_DIR / "database_and_literature_source_inventory.csv")
    _, ds_rows = load_csv(WORKSPACE_DIR / "ferritin_mutant_thermostability_dataset.csv")
    _, um_rows = load_csv(WORKSPACE_DIR / "unique_mutants_summary.csv")
    _, gap_rows = load_csv(WORKSPACE_DIR / "inaccessible_and_follow_up_sources.csv")
    report_text = (WORKSPACE_DIR / "README_database_summary_report.md").read_text(encoding="utf-8")
    report_lines = report_text.splitlines()

    # 2. Unique IDs checks
    check_unique_ids(inv_rows, "source_id", r"^SRC_[A-Z]+_\d+$", "source_id")
    fmt_ids = check_unique_ids(ds_rows, "record_id", r"^FMT\d{4}$", "record_id")
    um_ids = check_unique_ids(um_rows, "unique_mutant_id", r"^UM\d{4}$", "unique_mutant_id")
    check_unique_ids(gap_rows, "gap_id", r"^GAP\d{4}$", "gap_id")

    # 3. WT exclusion and bidirectional foreign-key integrity
    mut_rows = [r for r in ds_rows if r["record_role"] == "Mutant"]
    wt_rows = [r for r in ds_rows if r["record_role"] == "WT_baseline"]
    if len(mut_rows) + len(wt_rows) != len(ds_rows):
        fail("Unexpected record_role values in dataset")

    for w in wt_rows:
        if w["unique_mutant_id"] or w["unique_mutant_key"]:
            fail(f"WT_baseline row {w['record_id']} must not have unique_mutant_id/key")
        if w["mutation_category"] != "WT":
            fail(f"WT_baseline row {w['record_id']} has mutation_category={w['mutation_category']}")

    ds_by_id = {r["record_id"]: r for r in ds_rows}
    um_by_id = {u["unique_mutant_id"]: u for u in um_rows}
    grouped_muts: dict[str, list[dict]] = defaultdict(list)

    for m in mut_rows:
        uid = m["unique_mutant_id"]
        if uid not in um_ids:
            fail(f"Mutant record {m['record_id']} references unknown unique_mutant_id={uid!r}")
        if m["unique_mutant_key"] != um_by_id[uid]["unique_mutant_key"]:
            fail(f"Key mismatch for {m['record_id']} vs {uid}")
        if m["mutation_category"] == "WT":
            fail(f"Mutant record {m['record_id']} labeled mutation_category=WT")
        if m["family_tier"] not in VALID_FAMILIES:
            fail(f"Invalid family_tier on {m['record_id']}: {m['family_tier']}")
        if m["evidence_tier"] not in VALID_EVIDENCE:
            fail(f"Invalid evidence_tier on {m['record_id']}: {m['evidence_tier']}")
        if m["mutation_category"] not in VALID_CATEGORIES:
            fail(f"Invalid mutation_category on {m['record_id']}: {m['mutation_category']}")
        grouped_muts[uid].append(m)

    # Check unique_mutants_summary rows match grouped dataset rows 100%
    if set(grouped_muts.keys()) != um_ids:
        fail("Mismatch between unique_mutant_ids referenced in dataset and unique_mutants_summary.csv")

    for u in um_rows:
        uid = u["unique_mutant_id"]
        if "WT" in (u["unique_mutant_key"], u["mutation_category"]):
            fail(f"WT leaked into unique_mutants_summary: {uid}")
        grp = grouped_muts[uid]
        rec_ids_summary = [x for x in u["record_ids"].split(";") if x]
        rec_ids_ds = [x["record_id"] for x in grp]
        if rec_ids_summary != rec_ids_ds:
            fail(f"{uid} record_ids mismatch: {rec_ids_summary} vs {rec_ids_ds}")
        if int(u["n_condition_records"]) != len(grp):
            fail(f"{uid} n_condition_records mismatch")
        wet_cnt = sum(1 for x in grp if x["evidence_tier"] != "D_Computational_Prediction")
        quant_cnt = sum(1 for x in grp if x["evidence_tier"] == "A_WetLab_Quantitative" and x["value_numeric_C"])
        if int(u["n_wetlab_records"]) != wet_cnt:
            fail(f"{uid} n_wetlab_records mismatch: {u['n_wetlab_records']} vs {wet_cnt}")
        if int(u["n_quant_temp_records"]) != quant_cnt:
            fail(f"{uid} n_quant_temp_records mismatch: {u['n_quant_temp_records']} vs {quant_cnt}")
        best_ev = min(grp, key=lambda x: RANK_TIER[x["evidence_tier"]])["evidence_tier"]
        if u["highest_evidence_tier"] != best_ev:
            fail(f"{uid} highest_evidence_tier mismatch: {u['highest_evidence_tier']} vs {best_ev}")
        # Verify cross-ref duplicate groups listed on summary exist in dataset rows
        summary_xrefs = {x for x in u["cross_ref_groups"].split(";") if x}
        ds_xrefs = {x["cross_ref_duplicate_group"] for x in grp if x["cross_ref_duplicate_group"]}
        if summary_xrefs != ds_xrefs:
            fail(f"{uid} cross_ref_groups mismatch: {summary_xrefs} vs {ds_xrefs}")

    # 4. Internal numeric, sign-convention, and cross_ref_duplicate_group checks
    xref_groups: dict[str, list[dict]] = defaultdict(list)
    for r in ds_rows:
        rid = r["record_id"]
        val_s = r["value_numeric_C"].strip()
        if val_s:
            try:
                val_f = float(val_s)
                if not (0.0 <= val_f <= 150.0):
                    fail(f"{rid} value_numeric_C out of [0, 150] °C range: {val_f}")
            except ValueError:
                fail(f"{rid} non-numeric value_numeric_C: {val_s!r}")

        wt_s = r["wt_baseline_C"].strip()
        norm_d_s = r["normalized_delta_Tm_C"].strip()
        conv = r["delta_sign_convention"].strip()
        rep_d = r["reported_delta_value"].strip()

        if conv == "WT_minus_mut_inverted":
            m_rep = re.match(r"^([-+]?\d+(?:\.\d+)?)", rep_d)
            if not m_rep or not norm_d_s:
                fail(f"{rid} WT_minus_mut_inverted missing parseable rep/normalized delta: {rep_d!r} / {norm_d_s!r}")
            else:
                rep_num = float(m_rep.group(1))
                norm_num = float(norm_d_s)
                if abs(norm_num + rep_num) > 0.15:
                    fail(f"{rid} inverted delta sign mismatch: reported={rep_num}, normalized={norm_num}")
        elif conv == "mut_minus_WT" and norm_d_s and norm_d_s not in {"~0", "<0"}:
            # Verify normalized_delta_Tm_C matches value_numeric_C - wt_baseline_C when both single floats
            if val_s and wt_s:
                try:
                    diff = float(val_s) - float(wt_s)
                    nd = float(norm_d_s)
                    if abs(diff - nd) > 0.25:
                        fail(f"{rid} delta mismatch: val({val_s})-wt({wt_s})={diff:.2f} vs normalized_delta={nd}")
                except ValueError:
                    pass

        if r["cross_ref_duplicate_group"]:
            xref_groups[r["cross_ref_duplicate_group"]].append(r)

    gap_blob = "\n".join(json_str for g in gap_rows for json_str in g.values())
    for xg, rows in xref_groups.items():
        if len(rows) < 2:
            # Must link to an external entry tracked in inaccessible_and_follow_up_sources.csv (e.g. PDB:6LJG)
            tokens = [t for t in xg.split("_") if len(t) >= 4]
            if not any(t in gap_blob for t in tokens):
                fail(f"cross_ref_duplicate_group {xg!r} has <2 rows ({len(rows)}) and no gap linkage")

    # 5. Exact numerical reconciliation between CSV aggregates and README_database_summary_report.md
    wet_mut_rows = [r for r in mut_rows if r["evidence_tier"] != "D_Computational_Prediction"]
    comp_mut_rows = [r for r in mut_rows if r["evidence_tier"] == "D_Computational_Prediction"]
    quant_mut_rows = [r for r in mut_rows if r["evidence_tier"] == "A_WetLab_Quantitative"]
    quant_exact_rows = [r for r in quant_mut_rows if r["value_numeric_C"]]

    wet_uniqs = [u for u in um_rows if u["highest_evidence_tier"] != "D_Computational_Prediction"]
    comp_only_uniqs = [u for u in um_rows if u["highest_evidence_tier"] == "D_Computational_Prediction"]
    comp_and_wet_overlap = sum(
        1 for u in um_rows if u["highest_evidence_tier"] != "D_Computational_Prediction" and u["computational_predictions"]
    )
    quant_uniqs = [u for u in um_rows if u["highest_evidence_tier"] == "A_WetLab_Quantitative"]

    core_mut_rows = [r for r in mut_rows if r["family_tier"] in CORE_FAMILIES]
    core_wet_rows = [r for r in wet_mut_rows if r["family_tier"] in CORE_FAMILIES]
    core_uniqs = [u for u in um_rows if u["family_tier"] in CORE_FAMILIES]
    core_wet_uniqs = [u for u in wet_uniqs if u["family_tier"] in CORE_FAMILIES]

    by_fam_rec = Counter(r["family_tier"] for r in mut_rows)
    by_fam_wet_rec = Counter(r["family_tier"] for r in wet_mut_rows)
    by_fam_uniq = Counter(u["family_tier"] for u in um_rows)
    by_fam_wet_uniq = Counter(u["family_tier"] for u in wet_uniqs)
    by_ev_rec = Counter(r["evidence_tier"] for r in mut_rows)
    by_ev_uniq = Counter(u["highest_evidence_tier"] for u in um_rows)
    by_cat_rec = Counter(r["mutation_category"] for r in mut_rows)
    by_cat_uniq = Counter(u["mutation_category"] for u in um_rows)
    by_gap_prio = Counter(g["priority"] for g in gap_rows)

    expected_report_snippets = [
        f"| **突变体总计（含全部 5 类家族 + 计算预测）** | **{len(mut_rows)}** | **{len(um_rows)}** |",
        f"| ↳ 其中：**核心铁蛋白家族（Tier 1–4：典型 24-mer、Bfr、Dps、工程融合）** | **{len(core_mut_rows)}** | **{len(core_uniqs)}** |",
        f"| **湿实验突变体总计（Tier A + B + C，剔除纯计算）** | **{len(wet_mut_rows)}** | **{len(wet_uniqs)}** |",
        f"| ↳ 其中：**核心铁蛋白湿实验（Tier 1–4，Tier A–C）** | **{len(core_wet_rows)}** | **{len(core_wet_uniqs)}** |",
        f"| **定量温度测量突变体（Tier A：DSC / CD / nanoDSF / Tycho / DSF / Rayleigh 等）** | **{len(quant_mut_rows)}**（其中精确单值数字 **{len(quant_exact_rows)}**） | **{len(quant_uniqs)}** |",
        f"| **纯计算预测条目（Tier D：FoldX ΔΔG + Robetta 虚拟丙氨酸扫描）** | **{len(comp_mut_rows)}** | **{len(comp_only_uniqs)}**（纯计算独占；另 **{comp_and_wet_overlap}** 个已在湿实验表重复验证） |",
        f"| **单独隔离的野生型基线对照（`record_role = WT_baseline`）** | **{len(wt_rows)}** | —（不计入突变体） |",
        f"| **合计** | **{len(mut_rows)}** | **{len(wet_mut_rows)}** | **{len(um_rows)}** | **{len(wet_uniqs)}** |",
        f"## 四、 各大数据库与文献引擎检索审计汇总（共 {len(inv_rows)} 个库/引擎）",
        f"## 六、 存在但暂不可获取 / 待补齐核对来源统计（共 {len(gap_rows)} 项）",
        f"闭源付费论文（{by_gap_prio['P1_High_Closed_ThermalMutant']} 篇高优先级）",
    ]
    for fam in VALID_FAMILIES:
        expected_report_snippets.append(
            f"| `{fam}` | {by_fam_rec[fam]} | {by_fam_wet_rec[fam]} | {by_fam_uniq[fam]} | {by_fam_wet_uniq[fam]} |"
        )
    for ev in VALID_EVIDENCE:
        expected_report_snippets.append(f"| `{ev}` | {by_ev_rec[ev]} | {by_ev_uniq[ev]} |")
    for cat in ("Point_Substitution", "Multi_Point_Combinatorial", "Deletion_Truncation_Frameshift_Swap", "Fusion_Display_Heteropolymer"):
        expected_report_snippets.append(f"| `{cat}` | {by_cat_rec[cat]} | {by_cat_uniq[cat]} |")
    for pr, cnt in sorted(by_gap_prio.items()):
        expected_report_snippets.append(f"| `{pr}` | {cnt} |")

    for snippet in expected_report_snippets:
        if snippet not in report_text:
            fail(f"Report mismatch — expected snippet not found in README_database_summary_report.md:\n  {snippet}")

    # Verify sum of family/evidence/category rows equals totals
    if sum(by_fam_rec.values()) != len(mut_rows) or sum(by_fam_uniq.values()) != len(um_rows):
        fail("Family sum mismatch")
    if sum(by_ev_rec.values()) != len(mut_rows) or sum(by_ev_uniq.values()) != len(um_rows):
        fail("Evidence tier sum mismatch")
    if sum(by_cat_rec.values()) != len(mut_rows) or sum(by_cat_uniq.values()) != len(um_rows):
        fail("Mutation category sum mismatch")

    if errors:
        print(f"VALIDATION FAILED WITH {len(errors)} ERROR(S):")
        for e in errors:
            print(f" - {e}")
        return 1

    print("ALL VALIDATION CHECKS PASSED (0 errors)")
    print(f"README_database_summary_report.md lines: {len(report_lines)}")
    print(f"Inventory sources: {len(inv_rows)}")
    print(f"Dataset rows: {len(ds_rows)} (Mutant={len(mut_rows)}, WT_baseline={len(wt_rows)})")
    print(f"Unique mutants: {len(um_rows)} (WetLab={len(wet_uniqs)}, CompOnly={len(comp_only_uniqs)}, Comp&WetOverlap={comp_and_wet_overlap}, QuantTierA={len(quant_uniqs)})")
    print(f"Core Ferritin (Tier 1-4): recs={len(core_mut_rows)} (WetLab={len(core_wet_rows)}), uniq={len(core_uniqs)} (WetLab={len(core_wet_uniqs)})")
    print("\n--- Breakdown by family_tier (all_rec / wet_rec / all_uniq / wet_uniq) ---")
    for fam in sorted(VALID_FAMILIES):
        print(f"  {fam}: {by_fam_rec[fam]} / {by_fam_wet_rec[fam]} / {by_fam_uniq[fam]} / {by_fam_wet_uniq[fam]}")
    print("\n--- Breakdown by evidence_tier (rec / uniq) ---")
    for ev in sorted(VALID_EVIDENCE):
        print(f"  {ev}: {by_ev_rec[ev]} / {by_ev_uniq[ev]}")
    print("\n--- Breakdown by mutation_category (rec / uniq) ---")
    for cat in ("Point_Substitution", "Multi_Point_Combinatorial", "Deletion_Truncation_Frameshift_Swap", "Fusion_Display_Heteropolymer"):
        print(f"  {cat}: {by_cat_rec[cat]} / {by_cat_uniq[cat]}")
    print("\n--- Breakdown by gap priority ---")
    for pr, cnt in sorted(by_gap_prio.items()):
        print(f"  {pr}: {cnt}")
    print("\n--- Cross-reference duplicate groups ---")
    for xg, rows in sorted(xref_groups.items()):
        print(f"  {xg}: {len(rows)} rows (unique_mutant_ids={len({r['unique_mutant_id'] for r in rows if r['unique_mutant_id']})})")
    return 0


if __name__ == "__main__":
    sys.exit(validate())
