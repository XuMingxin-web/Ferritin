#!/usr/bin/env python3
"""
Build the exhaustive Ferritin Mutant Thermostability & Source Inventory Database
(dual-count reconciliation: condition-level records vs. deduplicated unique mutants).

Outputs in /Users/a26683/Desktop/ferritinMTdatabase/:
  1. database_and_literature_source_inventory.csv
  2. ferritin_mutant_thermostability_dataset.csv
  3. unique_mutants_summary.csv
  4. inaccessible_and_follow_up_sources.csv
  5. README_database_summary_report.md
"""

import csv
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

WORKSPACE_DIR = Path("/Users/a26683/Desktop/ferritinMTdatabase")
SCRATCH_DIR = Path("/Users/a26683/.gemini/antigravity/brain/5596ef0e-2fbf-49e1-b621-055b635f1086/scratch")

FAM_FERRITIN = "Family_1_Canonical_24mer_Ferritin"
FAM_BFR = "Family_2_Bacterioferritin_Bfr"
FAM_DPS = "Family_3_Miniferritin_Dps"
FAM_FUSION = "Family_4_Engineered_Fusion_Display_Ferritin"
FAM_EXT = "Family_5_Extended_Superfamily_And_Nanocages"

EV_A = "A_WetLab_Quantitative"
EV_B = "B_WetLab_Thermal_Survival"
EV_C = "C_WetLab_Qualitative_Rank"
EV_D = "D_Computational_Prediction"

CAT_WT = "WT"
CAT_POINT = "Point_Substitution"
CAT_MULTI = "Multi_Point_Combinatorial"
CAT_INDEL = "Deletion_Truncation_Frameshift_Swap"
CAT_FUSION = "Fusion_Display_Heteropolymer"


def find_scratch_file(filename: str) -> Path:
    """Locate a cached file under the brain scratch directory."""
    for root, _, files in os.walk(SCRATCH_DIR):
        if filename in files:
            return Path(root) / filename
    raise FileNotFoundError(f"Required scratch file not found: {filename}")


def build_source_inventory() -> list[dict]:
    """Construct verified inventory of all structured DBs and literature corpora queried."""
    rows = [
        # Structured biophysical/ML stability databases
        {
            "source_id": "SRC_DB_01",
            "category": "Structured_Stability_DB",
            "name": "ProThermDB (web query cache + ferritin/RNR TSVs)",
            "scope_or_query": "protein/keyword/UniProt searches (ferritin, FTH, FTL, P02792, P02794, P02791, P02793, P0A998, P0ABD3, P0ABT2, 1FHA, 2FHA)",
            "total_records_in_source_or_hitcount": "13 rows (pt2_prot_ferritin) + 4 rows (pt2_prot_ferr_low); web portal ~31k total",
            "ferritin_related_hits": "17 ferritin/fusion rows (13 WT + 4 EBNA1-HuHF fusion)",
            "ferritin_mutant_temp_records": "4 (all synthetic N/C-term EBNA1-HuHF fusions from PMID:33825471; ProThermDB mislabeled MUTATION='wild-type')",
            "access_status": "Partial_Cached_FollowUp_Bulk",
            "notes": "0 point-substitution ferritin mutants found in targeted queries; bulk ~31k web download requires interactive export.",
        },
        {
            "source_id": "SRC_DB_02",
            "category": "Structured_Stability_DB",
            "name": "FireProtDB (2025-09-20 full PostgreSQL dump)",
            "scope_or_query": "Full 4.83 GB SQL dump stream (41,146,172 lines)",
            "total_records_in_source_or_hitcount": "64 SQL COPY tables scanned in full",
            "ferritin_related_hits": "1 protein (seq_id=8955 P02792 FTL, 40 measurement rows incl. 7 WT Tm & 4 WT dG from PMID:9681870)",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero_Mutants",
            "notes": "Confirmed 0 rows for seq_id=8955 with non-null mutant_id. (1 non-ferritin RNR E. coli alpha & 1 Dfx WT protein excluded.)",
        },
        {
            "source_id": "SRC_DB_03",
            "category": "Structured_Stability_DB",
            "name": "ThermoMutDB",
            "scope_or_query": "Full JSON dataset scan",
            "total_records_in_source_or_hitcount": "13,337 JSON records",
            "ferritin_related_hits": "0",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero",
            "notes": "0 ferritin / Bfr / Dps / encapsulin records across entire JSON.",
        },
        {
            "source_id": "SRC_DB_04",
            "category": "Structured_Stability_DB",
            "name": "GeoStab (S4346, S571, S669, S8754)",
            "scope_or_query": "Sequence/PDB/keyword scan of all 4 benchmark CSVs",
            "total_records_in_source_or_hitcount": "4,346 + 571 + 669 + 8,754 = 14,340 rows",
            "ferritin_related_hits": "0",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero",
            "notes": "0 ferritin hits across dTm (S4346/S571) and ddG (S669/S8754) datasets.",
        },
        {
            "source_id": "SRC_DB_05",
            "category": "Structured_DMS_ML_DB",
            "name": "ProteinGym v1.x (Substitutions + Indels)",
            "scope_or_query": "Metadata table scan (217 substitution + 66 indel assays)",
            "total_records_in_source_or_hitcount": "283 DMS assays",
            "ferritin_related_hits": "0 (A0A1I9GEU1_NEIME is Neisseria pilin PilE false positive)",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero",
            "notes": "No ferritin deep mutational scanning assay in ProteinGym.",
        },
        {
            "source_id": "SRC_DB_06",
            "category": "Structured_DMS_ML_DB",
            "name": "Tsuboyama Mega-scale Folding Stability (Zenodo)",
            "scope_or_query": "Domain metadata JSON scan",
            "total_records_in_source_or_hitcount": "9,728 JSON entries",
            "ferritin_related_hits": "0",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero",
            "notes": "Small single-domain folds only; no multimeric ferritin cages.",
        },
        {
            "source_id": "SRC_DB_07",
            "category": "Structured_DMS_ML_DB",
            "name": "MaveDB",
            "scope_or_query": "Score-set search JSON scan",
            "total_records_in_source_or_hitcount": "2,161 JSON entries",
            "ferritin_related_hits": "0",
            "ferritin_mutant_temp_records": "0",
            "access_status": "Exhaustively_Verified_Zero",
            "notes": "0 ferritin thermostability DMS score sets.",
        },
        {
            "source_id": "SRC_DB_08",
            "category": "Structured_Binding_DB",
            "name": "SKEMPI 2.0",
            "scope_or_query": "Full CSV scan (skempi_v2.csv)",
            "total_records_in_source_or_hitcount": "7,084 data rows (+1 header = 7,085)",
            "ferritin_related_hits": "7 rows (4E6K_AB_G: P. aeruginosa BfrB + Bfd, PMID:26368531)",
            "ferritin_mutant_temp_records": "0 (binding affinity Kd measured at 295 K only; no Tm/thermal unfolding)",
            "access_status": "Exhaustively_Verified_BindingOnly",
            "notes": "Logged in Gap/Clarification notes; excluded from Tm dataset as non-thermostability SPR Kd.",
        },
        {
            "source_id": "SRC_DB_09",
            "category": "Sequence_Feature_DB",
            "name": "UniProtKB (Mutagenesis & Natural Variants)",
            "scope_or_query": "uniprot_ferritin_mutagen.json (32), uniprot_ipr_mutagen_all.json (26), uniprot_ferritin_variant_sp.json (10)",
            "total_records_in_source_or_hitcount": "68 JSON entry records (overlapping sets)",
            "ferritin_related_hits": "Multiple canonical ferritins, Bfr, Dps, encapsulin entries",
            "ferritin_mutant_temp_records": "1 qualitative thermostability annotation (P9WNE5 MtBfrB del 168-181, PMID:21494619) + oligomer stability hits (Q2RVS1, Q9WZP2, Q9HY80)",
            "access_status": "Exhaustively_Verified_Qualitative",
            "notes": "No numeric Tm stored in UniProt mutagenesis annotations; P9WNE5 integrated as Tier C.",
        },
        {
            "source_id": "SRC_DB_10",
            "category": "Structural_DB",
            "name": "RCSB PDB",
            "scope_or_query": "Family queries (ferritin 2,759 total/2,000 fetched; Bfr 161; Dps 1; encapsulin 107) + mutant structure filter",
            "total_records_in_source_or_hitcount": "461 mutant PDB structures in rcsb_mutants_clean.json (70 therm/stab/del/fusion title hits)",
            "ferritin_related_hits": "461 mutant structures (incl. 38 To Be Published / no-PMID-DOI entries)",
            "ferritin_mutant_temp_records": "0 numerical Tm stored directly in PDB JSON; links to closed/TBP thermal studies (e.g., 8YKR, 9BXP, 9BYE, 9LUW, 6LJG, 6LBC/6LBD, 7DQO, 8H8L-O)",
            "access_status": "Verified_Structures_With_TBP_Gaps",
            "notes": "38 TBP/no-PMID-DOI mutant entries logged in inaccessible_and_follow_up_sources.csv.",
        },
        {
            "source_id": "SRC_DB_11",
            "category": "Unqueried_Portal_Gap",
            "name": "SoluProtMutDB / MPTherm / HotProtein / BRENDA (EC 1.16.3.1/1.16.3.2) / BioStudies",
            "scope_or_query": "Static HTTP checks + BioStudies coarse index (91,242 hits; S-EPMC6044891 luciferase false positive)",
            "total_records_in_source_or_hitcount": "Portal-dependent",
            "ferritin_related_hits": "Unresolved due to SPA JS rendering / anti-bot / portal query requirements",
            "ferritin_mutant_temp_records": "Unknown (Gap)",
            "access_status": "Inaccessible_JS_Or_Portal_FollowUp",
            "notes": "Flagged in inaccessible_and_follow_up_sources.csv for manual browser/export verification.",
        },
        # Literature engines & corpora
        {
            "source_id": "SRC_LIT_01",
            "category": "Literature_Index_API",
            "name": "Europe PMC REST API (10 query caches)",
            "scope_or_query": "epmc_ferritin_tm (3,000/9,803), epmc_ferritin_md_topt (3,000/5,994), extra_relatives (1,277), neuroferr (417), author_experts (326), title_thermal (255), nanocage (147), fulltext_dsc_tm (87), preprints (78), patents (0)",
            "total_records_in_source_or_hitcount": "7,083 deduplicated unique articles across 10 query JSONs (18,384 raw hits pooled)",
            "ferritin_related_hits": "249 master biophysical candidates (86 PMC structured candidates)",
            "ferritin_mutant_temp_records": "Bulk source for PMC full-text corpus below",
            "access_status": "Queried_And_Triaged (except 78 preprint full texts partial)",
            "notes": "78 Europe PMC preprints logged as follow-up.",
        },
        {
            "source_id": "SRC_LIT_02",
            "category": "Literature_Index_API",
            "name": "OpenAlex API (10 query caches)",
            "scope_or_query": "openalex_bfr_dps_enc (2,000), openalex_ferritin_tm (2,000), openalex_ferritin_md (1,604), oa_title_ferritin_mut_therm (528), oa_zh_title_ferr (600/3 files), oa_zh_ferritin_tm (100), oa_title_dps_mut_therm (69), oa_zh_ferritin_dsc (1)",
            "total_records_in_source_or_hitcount": "6,365 deduplicated unique works across 10 JSONs",
            "ferritin_related_hits": "Triage merged into master candidates + 600 deduped ZH scan",
            "ferritin_mutant_temp_records": "Contributed closed DOI & OA PDF leads",
            "access_status": "Queried_And_Merged",
            "notes": "Chinese metadata yielded 63 mut/fusion/thermal title hits (mostly non-OA/materials false positives; genuine ferritin fusions in PubMed chi).",
        },
        {
            "source_id": "SRC_LIT_03",
            "category": "Literature_Index_API",
            "name": "Crossref REST API",
            "scope_or_query": "crossref_spkx_ferritin (100 sampled / 493,683 total), crossref_zh_ferritin (100 / 514,592)",
            "total_records_in_source_or_hitcount": "200 sampled metadata works (broad >1M index)",
            "ferritin_related_hits": "Sampled for Chinese Food Science (SPKX) & closed journal discovery",
            "ferritin_mutant_temp_records": "Merged into ZH/closed DOI audit",
            "access_status": "Sampled_CNKI_FullText_Gap",
            "notes": "Chinese journal full texts (SPKX, 生物工程学报 etc.) require CNKI/Wanfang access.",
        },
        {
            "source_id": "SRC_LIT_04",
            "category": "Literature_Index_XML",
            "name": "PubMed XML + Thermal Non-PMC Abstracts",
            "scope_or_query": "11 PubMed XML batches (title/novel thermal queries)",
            "total_records_in_source_or_hitcount": "1,856 unique PMIDs; 105 prioritized non-PMC thermal abstracts parsed in full",
            "ferritin_related_hits": "29 high-value mutant+thermal non-PMC abstracts inspected verbatim",
            "ferritin_mutant_temp_records": "Contributed closed-paper abstract Tier A/B/C records (PMID:33406454, 42902949, 16704261, 23511498, 41700891, 26115686, etc.)",
            "access_status": "Abstracts_Verified_FullText_Paywalled",
            "notes": "Full PDFs of paywalled PubMed papers logged in inaccessible_and_follow_up_sources.csv.",
        },
        {
            "source_id": "SRC_LIT_05",
            "category": "FullText_OpenAccess",
            "name": "PubMed Central (PMC) Full-Text Corpus (XML + NCBI XML + HTML)",
            "scope_or_query": "231 pmc_xml + 16 pmc_ncbi_xml + 22 pmc_html",
            "total_records_in_source_or_hitcount": "269 unique nonempty PMC full texts audited",
            "ferritin_related_hits": "All 269 verified by programmatic table/keyword sweep + manual chunk reading",
            "ferritin_mutant_temp_records": "Core OA experimental wet-lab dataset (PMC2323953, PMC5729161, PMC6149950, PMC6800815, PMC9069879, PMC3403414, PMC9785301, PMC9798250, PMC11776510, PMC13406861, PMC152977, PMC10093968, PMC7746619, PMC11838010, PMC8357837, PMC6408514, PMC5812362, PMC4646333, PMC9814774, PMC1133592, PMC13319912, PMC3059798, PMC10367587, PMC13107372, PMC12124451, PMC11190991, PMC11284452, PMC11550324, PMC8816334, PMC11866011, PMC2852946, PMC3267954, etc.)",
            "access_status": "FullText_Verified (some ESI/Figs flagged)",
            "notes": "Note: PMC2852946/PMC3267954 numeric tables were image-embedded in PMC; exact values verified via co-author Yu Zhang DR-NTU 48015 thesis.",
        },
        {
            "source_id": "SRC_LIT_06",
            "category": "FullText_OA_PDF_Thesis",
            "name": "OA PDFs, Preprints & Institutional Theses (Unpaywall 127 OA -> 18 local PDFs + DR-NTU 48015)",
            "scope_or_query": "Unpaywall check on 260 candidate DOIs (127 OA, 133 closed); 18 local PDFs decoded incl. NTU thesis DR-NTU 48015 (Yu Zhang) & pdf_1007_decoded.txt",
            "total_records_in_source_or_hitcount": "18 local PDFs (17 pdfkit + 1 decoded CMap PDF)",
            "ferritin_related_hits": "DR-NTU 48015 complete thesis + 10.1021/acs.cgd.3c00880 + PMID:33900431 (MhDps Q14E) + 10.1101/2025.04.26.649701 + 10.1038/s41467-024-52908-z (PMC11452399)",
            "ferritin_mutant_temp_records": "Major quantitative source: DR-NTU 48015 Ch2 Table 2.7 (Bfr Ala-scan), Ch3 Table 3.3 & Fig 3.17 & Table 3.5, Ch4 Table 4.3 (Dps Ala-scan), plus FoldX Table 3.1 & Robetta Tables 2.3/4.1",
            "access_status": "FullText_Verified (Unpaywall 133 closed DOIs logged)",
            "notes": "133 Unpaywall-closed DOIs (102 ferritin-title + 31 other-title) explicitly cataloged in inaccessible_and_follow_up_sources.csv.",
        },
    ]
    return rows


def normalize_mutation(mut: str) -> str:
    """Create canonical deduplication representation for mutation strings."""
    s = mut.strip()
    s = re.sub(r"\s+", " ", s)
    return s


def make_records() -> list[dict]:
    """Assemble all strictly verified condition-level records (WT baselines + mutants)."""
    records = []

    def add(
        record_role: str,
        family_tier: str,
        evidence_tier: str,
        mutation_category: str,
        organism: str,
        protein_name: str,
        uniprot_id: str,
        pdb_id: str,
        mutation_raw: str,
        mutation_norm: str,
        interface_domain: str,
        metric_name: str,
        value_numeric_C: str,
        value_display: str,
        value_precision: str,
        wt_baseline_C: str,
        reported_delta_value: str,
        delta_sign_convention: str,
        normalized_delta_Tm_C: str,
        delta_Tm_provenance: str,
        assay_method: str,
        state_apo_holo: str,
        buffer_pH_salt_conc: str,
        oligomeric_state: str,
        source_type: str,
        pmid: str,
        pmcid: str,
        doi: str,
        citation_short: str,
        cross_ref_group: str,
        notes: str,
    ):
        idx = len(records) + 1
        ukey = "" if record_role == "WT_baseline" else f"{organism}||{protein_name}||{mutation_norm}"
        records.append({
            "record_id": f"FMT{idx:04d}",
            "record_role": record_role,
            "family_tier": family_tier,
            "evidence_tier": evidence_tier,
            "mutation_category": mutation_category,
            "organism": organism,
            "protein_name": protein_name,
            "uniprot_id": uniprot_id,
            "pdb_id": pdb_id,
            "mutation_raw": mutation_raw,
            "mutation_normalized": mutation_norm,
            "unique_mutant_key": ukey,
            "unique_mutant_id": "",  # populated cleanly during rollup
            "symmetry_interface_or_domain": interface_domain,
            "metric_name": metric_name,
            "value_numeric_C": value_numeric_C,
            "value_display": value_display,
            "value_precision": value_precision,
            "wt_baseline_C": wt_baseline_C,
            "reported_delta_value": reported_delta_value,
            "delta_sign_convention": delta_sign_convention,
            "normalized_delta_Tm_C": normalized_delta_Tm_C,
            "delta_Tm_provenance": delta_Tm_provenance,
            "assay_method": assay_method,
            "state_apo_holo": state_apo_holo,
            "buffer_pH_salt_conc": buffer_pH_salt_conc,
            "oligomeric_state": oligomeric_state,
            "source_type": source_type,
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "citation_short": citation_short,
            "cross_ref_duplicate_group": cross_ref_group,
            "notes_and_caveats": notes,
        })

    # =========================================================================
    # 1. DR-NTU 48015 (Yu Zhang NTU Thesis) + PMC2852946 (E. coli Bfr)
    # NOTE: Thesis reports Delta Tm = Tm_WT - Tm_mut (positive = destabilized).
    # =========================================================================
    bfr_buf_27 = "50 mM NaH2PO4, pH 7.2; 100 ug/mL protein; 20-90 C (5 C steps, 8 min equil)"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR/2VXI",
        "WT", "WT", "WT_reference", "Tm", "64.8", "64.8 ± 1.1 °C", "Exact_Numeric", "64.8", "0", "NA", "0.0", "NA",
        "CD_222nm", "Apo/As-isolated", bfr_buf_27, "24-mer (45%) + Dimer (55%); DLS 13.3 ± 0.2 nm",
        "Thesis_OpenAccess|PMC_OpenAccess", "20139406", "PMC2852946", "10.1074/jbc.M109.092445",
        "Zhang Y (DR-NTU 48015 Table 2.7) / Zhang et al. 2010 JBC", "DRNTU48015_Ch2_vs_PMC2852946",
        "Baseline for Table 2.7 EcBfr Ala-shaving panel. PMC2852946 has image-embedded Table 1; exact numbers taken from co-author DR-NTU 48015.")

    bfr_t27 = [
        ("E157A", "C4", "64.9", "64.9 ± 0.7 °C", "-0.1 ± 1.8", "+0.1", "24-mer (41%) + Dimer (59%)", "Virtual Ala ddG*=0.82 kcal/mol"),
        ("R61A", "C3", "68.8", "68.8 ± 0.8 °C", "-4.0 ± 1.9", "+4.0", "Dimer (100%); DLS 6.1 ± 0.3 nm", "Virtual Ala ddG*=4.36 kcal/mol; dimer more thermostable than WT"),
        ("Y114A", "C3", "48.4", "48.4 ± 0.8 °C", "16.4 ± 1.9", "-16.4", "Dimer (100%); DLS 6.8 ± 0.3 nm", "Virtual Ala ddG*=2.43 kcal/mol"),
        ("E128A", "C3", "46.9", "46.9 ± 0.2 °C", "17.9 ± 1.3", "-17.9", "Dimer (100%); DLS 6.7 ± 0.6 nm", "Virtual Ala ddG*=5.71 kcal/mol"),
        ("R30A", "C2", "72.2", "72.2 ± 0.6 °C", "-7.4 ± 1.7", "+7.4", "Dimer (100%); DLS 6.0 ± 0.2 nm", "Virtual Ala ddG*=9.22 (Table 2.7) / 9.14 (Table 2.3) kcal/mol"),
        ("K33A", "C2", "56.9", "56.9 ± 1.5 °C", "7.9 ± 2.6", "-7.9", "24-mer (4%) + Dimer (96%)", "Virtual Ala ddG*=0.8 kcal/mol"),
        ("Y45A", "C2", "58.6", "58.6 ± 0.4 °C", "6.2 ± 1.5", "-6.2", "24-mer (22%) + Dimer (78%)", "Virtual Ala ddG*=1.48 kcal/mol"),
        ("D56A", "C2", "51.8", "51.8 ± 0.3 °C", "13.0 ± 1.4", "-13.0", "24-mer (16%) + Dimer (84%)", "Virtual Ala ddG*=4.54 kcal/mol"),
        ("D56A/E60A", "C2", "48.6", "48.6 ± 0.6 °C", "16.2 ± 1.7", "-16.2", "24-mer (24%) + Dimer (76%)", "Double C2 interface mutant"),
    ]
    for mut, iface, val, disp, repd, normd, oligo, note in bfr_t27:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_BFR, EV_A, cat, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR/2VXI",
            mut, mut, iface, "Tm", val, disp, "Exact_Numeric", "64.8", repd, "WT_minus_mut_inverted", normd,
            "Sign_Inverted_From_Reported", "CD_222nm", "Apo/As-isolated", bfr_buf_27, oligo,
            "Thesis_OpenAccess|PMC_OpenAccess", "20139406", "PMC2852946", "10.1074/jbc.M109.092445",
            "Zhang Y (DR-NTU 48015 Table 2.7) / Zhang et al. 2010 JBC", "DRNTU48015_Ch2_vs_PMC2852946", note)

    # E60A in Table 2.7: no transition detected ("-")
    add("Mutant", FAM_BFR, EV_C, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR/2VXI",
        "E60A", "E60A", "C2", "Tm", "", "No transition detected (-)", "Qualitative", "64.8", "-", "NA", "", "NA",
        "CD_222nm", "Apo/As-isolated", bfr_buf_27, "Unknown (no CD transition)",
        "Thesis_OpenAccess|PMC_OpenAccess", "20139406", "PMC2852946", "10.1074/jbc.M109.092445",
        "Zhang Y (DR-NTU 48015 Table 2.7) / Zhang et al. 2010 JBC", "DRNTU48015_Ch2_vs_PMC2852946",
        "Table 2.7 lists '-' (no melting transition observed in CD 222 nm; Virtual Ala ddG*=2.07 kcal/mol).")

    # DR-NTU Ch3 Table 3.3 (charged C2 mutants, 150-200 ug/mL)
    bfr_buf_33 = "50 mM NaH2PO4 (150-200 ug/mL); CD 222 nm"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
        "WT", "WT", "WT_reference", "Tm", "68.2", "68.2 °C", "Exact_Numeric", "68.2", "0", "NA", "0.0", "NA",
        "CD_222nm", "Apo/As-isolated", bfr_buf_33, "24-mer & Dimer; Reversibility 71%",
        "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Table 3.3)", "DRNTU48015_Ch3_vs_PMID26115686",
        "Baseline for Chapter 3 C2 charged interface mutants.")
    bfr_t33 = [
        ("R30D", "C2", "59.2", "9.0", "-9.0", "Dimer (solution); Reversibility 81%", "FoldX ddG=12.5 kcal/mol"),
        ("D56R", "C2", "64.4", "3.8", "-3.8", "Dimer (solution); Reversibility 66%", "FoldX ddG=6.8 kcal/mol"),
        ("E60H", "C2", "72.3", "-4.1", "+4.1", "24-mer & Dimer; Reversibility 57%", "FoldX ddG=1.3 kcal/mol"),
        ("D56R/E60H", "C2", "60.8", "7.4", "-7.4", "Dimer (solution); Reversibility 53%", "FoldX ddG=8.5 kcal/mol. Compare PMID:26115686 abstract (-16 C in J. Biochem. 2015 assay)."),
    ]
    for mut, iface, val, repd, normd, oligo, note in bfr_t33:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_BFR, EV_A, cat, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            mut, mut, iface, "Tm", val, f"{val} °C", "Exact_Numeric", "68.2", repd, "WT_minus_mut_inverted", normd,
            "Sign_Inverted_From_Reported", "CD_222nm", "Apo/As-isolated", bfr_buf_33, oligo,
            "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Table 3.3)", "DRNTU48015_Ch3_vs_PMID26115686", note)

    # PMID:26115686 J Biochem 2015 abstract counterpart
    add("Mutant", FAM_BFR, EV_A, CAT_MULTI, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
        "D56R-E60H", "D56R/E60H", "C2", "delta_Tm", "", "-16 °C vs WT (16 °C less stable)", "Exact_Numeric",
        "", "-16 °C (16 °C less stable)", "mut_minus_WT", "-16.0", "Directly_Reported",
        "CD/Thermal_Abstract", "As-isolated", "J Biochem 2015 (closed full text)", "Dimer in solution; cage in dried TEM",
        "PubMed_Abstract_ClosedPaper", "26115686", "", "10.1093/jb/mvv065", "Ardejani et al. 2015 J Biochem", "DRNTU48015_Ch3_vs_PMID26115686",
        "Abstract explicitly states D56R-E60H is 16 C less stable than WT (differs from DR-NTU Table 3.3 delta=-7.4 C); closed DOI.")
    for mut in ("R30D", "D56R", "E60H"):
        add("Mutant", FAM_BFR, EV_C, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            mut, mut, "C2", "Qualitative_Stability", "", "< WT thermostability", "Qualitative", "", "", "NA", "", "Range_Or_Qualitative",
            "CD/Thermal_Abstract", "As-isolated", "J Biochem 2015 (closed full text)", "R30D/D56R dimer; E60H partial cage",
            "PubMed_Abstract_ClosedPaper", "26115686", "", "10.1093/jb/mvv065", "Ardejani et al. 2015 J Biochem", "DRNTU48015_Ch3_vs_PMID26115686",
            "Abstract states 'they all show a lowered thermostability'; exact J Biochem Tm values behind OUP paywall (note E60H was +4.1 C in thesis Table 3.3).")

    # DR-NTU Ch3 Fig 3.17 VP-DSC isolated dimers
    dsc_buf_317 = "50 mM NaH2PO4, pH 7.2; 4.0-6.0 uM protein; VP-DSC 90 C/h (isolated dimers)"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
        "WT (isolated dimer)", "WT (isolated dimer)", "WT_dimer_reference", "Tm1", "77.0", "77 °C (ΔHcal 43.1, ΔHvH 65.6 kcal/mol)",
        "Exact_Numeric", "77.0", "0", "NA", "0.0", "NA", "DSC", "Apo/As-isolated", dsc_buf_317, "Isolated Dimer",
        "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Fig 3.17)", "", "First DSC transition Tm1 of isolated WT EcBfr dimer.")
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
        "WT (isolated dimer)", "WT (isolated dimer)", "WT_dimer_reference", "Tm2", "86.0", "86 °C (ΔHcal 62.0, ΔHvH 100.0 kcal/mol)",
        "Exact_Numeric", "86.0", "0", "NA", "0.0", "NA", "DSC", "Apo/As-isolated", dsc_buf_317, "Isolated Dimer",
        "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Fig 3.17)", "", "Second DSC transition Tm2 of isolated WT EcBfr dimer.")

    for mut, iface, tm1, h1, vh1, tm2, h2, vh2 in [
        ("R30A", "C2", 86, 65.7, 62.3, 92, 54.4, 130.0),
        ("Y114A", "C3", 80, 98.9, 18.6, 86, 58.0, 82.6),
        ("E128A", "C3", 70, 75.5, 47.9, 80, 75.1, 75.4),
    ]:
        add("Mutant", FAM_BFR, EV_A, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            f"{mut} (dimer)", mut, iface, "Tm1", str(tm1), f"{tm1} °C (ΔHcal {h1}, ΔHvH {vh1} kcal/mol)",
            "Exact_Numeric", "77.0", f"+{tm1-77}" if tm1 >= 77 else str(tm1-77), "mut_minus_WT", f"{tm1-77:+d}" if tm1 != 77 else "0",
            "Derived_mut_minus_WT", "DSC", "Apo/As-isolated", dsc_buf_317, "Isolated Dimer",
            "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Fig 3.17)", "",
            "VP-DSC isolated dimer fraction; transition 1 (Tm1).")
        add("Mutant", FAM_BFR, EV_A, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            f"{mut} (dimer)", mut, iface, "Tm2", str(tm2), f"{tm2} °C (ΔHcal {h2}, ΔHvH {vh2} kcal/mol)",
            "Exact_Numeric", "86.0", f"+{tm2-86}" if tm2 > 86 else str(tm2-86), "mut_minus_WT", f"{tm2-86:+d}" if tm2 != 86 else "0",
            "Derived_mut_minus_WT", "DSC", "Apo/As-isolated", dsc_buf_317, "Isolated Dimer",
            "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Fig 3.17)", "",
            "VP-DSC isolated dimer fraction; transition 2 (Tm2).")

    # DR-NTU Table 3.5 WT Bfr pH profile
    for ph, tm_val, ol in [
        ("3.0", "73.2", "Monomer"),
        ("3.5", "93.9", "24-mer + Dimer"),
        ("4.5", "88.9", "24-mer + Dimer"),
        ("7.0", "67.4", "24-mer"),
        ("9.0", "55.8", "Dimer"),
        ("10.0", "51.4", "Dimer"),
    ]:
        add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            f"WT (pH {ph})", f"WT_pH_{ph}", "WT_pH_series", "Tm", tm_val, f"{tm_val} °C",
            "Exact_Numeric", "", "", "NA", "", "NA", "CD_222nm", "Apo/As-isolated",
            f"50 mM phosphate, pH {ph}; 150-200 ug/mL; 1 C/min", ol,
            "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Table 3.5)", "",
            "WT pH stability sweep (pH 1, 2, 5, 6, 11, 12 had no cooperative transition '-').")

    # =========================================================================
    # 2. DR-NTU 48015 Ch4 Table 4.3 + PMC3267954 (E. coli Dps)
    # =========================================================================
    dps_buf = "50 mM NaH2PO4, pH 7.2; 100 ug/mL protein; 20-90 C (5 C steps, 8 min equil)"
    add("WT_baseline", FAM_DPS, EV_A, CAT_WT, "Escherichia coli", "DNA-binding protein from starved cells (EcDps)", "P0ABT2", "1DPS",
        "WT", "WT", "WT_reference", "Tm", "66.7", "66.7 ± 1.0 °C", "Exact_Numeric", "66.7", "0", "NA", "0.0", "NA",
        "CD_222nm", "Apo/As-isolated", dps_buf, "12-mer (100%)",
        "Thesis_OpenAccess|PMC_OpenAccess", "21898653", "PMC3267954", "10.1002/pro.712",
        "Zhang Y (DR-NTU 48015 Table 4.3) / Zhang et al. 2011 Protein Sci", "DRNTU48015_Ch4_vs_PMC3267954",
        "Baseline for EcDps Ala-shaving panel. PMC3267954 Table 1/Fig 7 image-embedded; exact values taken from DR-NTU 48015 Table 4.3.")

    dps_t43 = [
        ("Y16A", "N-term / C3_vicinity", "54.8", "54.8 ± 0.3 °C", "11.9 ± 1.3", "-11.9", "12-mer (84%) + Dimer (16%)", "Solvent-exposed N-term Tyr"),
        ("R18A", "C3", "55.2", "55.2 ± 0.4 °C", "11.5 ± 1.4", "-11.5", "12-mer (100%)", "Robetta ddG_bind=3.93 kcal/mol"),
        ("R83A", "C3", "50.4", "50.4 ± 0.7 °C", "16.3 ± 1.7", "-16.3", "12-mer (61%) + Dimer (39%)", "Robetta ddG_bind=9.84 kcal/mol"),
        ("R133A", "C3", "48.1", "48.1 ± 0.4 °C", "18.6 ± 1.4", "-18.6", "12-mer (51%) + Dimer (49%)", "Robetta ddG_bind=6.11 kcal/mol"),
        ("R83A/R133A", "C3", "41.4", "41.4 ± 0.4 °C", "25.3 ± 1.4", "-25.3", "Dimer (100%); reversible; no TEM cages", "Double C3 salt-bridge mutant; full 12-mer disassembly"),
        ("D146A", "C3 (Fe_channel)", "66.6", "66.6 ± 0.9 °C", "0.1 ± 1.9", "-0.1", "12-mer (100%)", "Robetta ddG_bind=1.48 kcal/mol"),
        ("W160A", "C3 (ferritin-like & Dps-like)", "57.3", "57.3 ± 0.7 °C", "9.4 ± 1.7", "-9.4", "12-mer (100%)", "Robetta ddG_bind=2.42 kcal/mol"),
        ("W52A", "C2", "56.3", "56.3 ± 0.4 °C", "10.4 ± 1.4", "-10.4", "12-mer (94%) + Monomer (6%)", "Robetta ddG_bind=3.35 kcal/mol"),
        ("D78A", "C2", "61.9", "61.9 ± 0.7 °C", "4.8 ± 1.7", "-4.8", "12-mer (100%)", "Robetta ddG_bind=1.09 kcal/mol"),
    ]
    for mut, iface, val, disp, repd, normd, oligo, note in dps_t43:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_DPS, EV_A, cat, "Escherichia coli", "DNA-binding protein from starved cells (EcDps)", "P0ABT2", "1DPS",
            mut, mut, iface, "Tm", val, disp, "Exact_Numeric", "66.7", repd, "WT_minus_mut_inverted", normd,
            "Sign_Inverted_From_Reported", "CD_222nm", "Apo/As-isolated", dps_buf, oligo,
            "Thesis_OpenAccess|PMC_OpenAccess", "21898653", "PMC3267954", "10.1002/pro.712",
            "Zhang Y (DR-NTU 48015 Table 4.3) / Zhang et al. 2011 Protein Sci", "DRNTU48015_Ch4_vs_PMC3267954", note)

    # =========================================================================
    # 3. Computational Tables in DR-NTU 48015 (FoldX Table 3.1 + Robetta 2.3 & 4.1)
    # =========================================================================
    # FoldX Table 3.1 (17 C2 mutants on EcBfr)
    foldx_t31 = [
        ("R30D", "12.46 ± 0.43"), ("R30E", "11.54 ± 0.39"), ("R30A", "5.45 ± 0.66"),
        ("D56R", "6.76 ± 0.31"), ("D56H", "6.23 ± 0.12"), ("D56K", "5.6 ± 0.17"), ("D56A", "4 ± 0.04"),
        ("E60H", "1.27 ± 0.11"), ("E60R", "-1.23 ± 0.45"), ("E60K", "-0.45 ± 0.00"), ("E60A", "-0.15 ± 0.03"),
        ("D56R/E60H", "8.5 ± 0.52"), ("D56R/E60R", "6.6 ± 0.13"), ("D56K/E60K", "6.72 ± 0.46"),
        ("D56R/E60K", "7.29 ± 0.31"), ("D56K/E60R", "5.86 ± 0.21"), ("D56A/E60A", "4.38 ± 0.16"),
    ]
    for mut, val in foldx_t31:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_BFR, EV_D, cat, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            mut, mut, "C2", "FoldX_ddG", "", f"{val} kcal/mol", "Computational_Energy", "", val, "kcal_per_mol", "", "NA",
            "FoldX", "In_Silico", "FoldX BuildModel at pH 7.0", "In_Silico",
            "Thesis_OpenAccess", "", "", "10.32657/10356/75835", "Zhang Y (DR-NTU 48015 Table 3.1)", "",
            "Positive ddG indicates predicted destabilization.")

    # Parse DR-NTU Table 2.3 (47 EcBfr Robetta Virtual Ala Scan rows) from step3_table23_pmc881.txt
    t23_path = find_scratch_file("step3_table23_pmc881.txt")
    t23_known_bfr_aa = {30: "R30A", 33: "K33A", 56: "D56A", 60: "E60A", 61: "R61A", 102: "R102A", 114: "Y114A", 128: "E128A", 157: "E157A"}
    t23_rows = []
    for ln in open(t23_path, encoding="utf-8"):
        m = re.match(r"^(\d+)\s+A\s+([01])\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+(C[234])\s*$", ln.strip())
        if m:
            res_num = int(m.group(1))
            int_id = m.group(2)
            ddg_bind = m.group(3)
            dg_partner = m.group(4)
            iface = m.group(5)
            t23_rows.append((res_num, int_id, ddg_bind, dg_partner, iface))
    assert len(t23_rows) == 47, f"Expected 47 rows in Table 2.3, got {len(t23_rows)}"
    for res_num, int_id, ddg_bind, dg_partner, iface in t23_rows:
        m_norm = t23_known_bfr_aa.get(res_num, f"pos{res_num}->Ala")
        add("Mutant", FAM_BFR, EV_D, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
            f"res{res_num}->Ala (Chain A)", m_norm, iface, "Robetta_ddG_bind", "",
            f"{ddg_bind} kcal/mol (dG_partner={dg_partner})", "Computational_Energy", "", ddg_bind, "kcal_per_mol", "", "NA",
            "Robetta_AlaScan", "In_Silico", "http://robetta.bakerlab.org/alaninescan", "In_Silico",
            "Thesis_OpenAccess", "20139406", "PMC2852946", "10.32657/10356/75835",
            "Zhang Y (DR-NTU 48015 Table 2.3)", "",
            f"int_id={int_id}; positive ddG(bind) = predicted destabilization upon Ala substitution.")

    # Parse DR-NTU Table 4.1 unambiguous EcDps Robetta Virtual Ala Scan rows directly from drntu_48015_pdfkit.txt
    dr_path = find_scratch_file("drntu_48015_pdfkit.txt")
    dr_lines = open(dr_path, encoding="utf-8").read().splitlines()
    dps_known_aa = {18: "R18A", 20: "D20A", 48: "K48A", 52: "W52A", 70: "R70A", 78: "D78A", 83: "R83A",
                    123: "D123A", 133: "R133A", 143: "D143A", 146: "D146A", 156: "D156A", 160: "W160A", 163: "E163A"}
    t41_rows = []
    # Scan clean single-entry lines in splitlines() range 4365..4450 plus explicitly verified 18 (pos 4371)
    for idx in range(4365, 4450):
        st = dr_lines[idx].strip()
        m = re.match(r"^(\d+)\s+F\s+([01])\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+(C[23])\s*$", st)
        if m:
            t41_rows.append((int(m.group(1)), m.group(2), m.group(3), m.group(4), m.group(5)))
    # Add 18 (dr_lines[4371]: "   17 18 F 1 3.93 -0.21 C3" where 18=R18A matches Ch4 Table 4.3 ddG*=3.93!)
    t41_rows.insert(0, (18, "1", "3.93", "-0.21", "C3"))
    assert len(t41_rows) == 54, f"Expected 54 unambiguous clean rows in Table 4.1, got {len(t41_rows)}"
    for res_num, int_id, ddg_bind, dg_partner, iface in t41_rows:
        m_norm = dps_known_aa.get(res_num, f"pos{res_num}->Ala")
        add("Mutant", FAM_DPS, EV_D, CAT_POINT, "Escherichia coli", "DNA-binding protein from starved cells (EcDps)", "P0ABT2", "1DPS",
            f"res{res_num}->Ala (Chain F)", m_norm, iface, "Robetta_ddG_bind", "",
            f"{ddg_bind} kcal/mol (dG_partner={dg_partner})", "Computational_Energy", "", ddg_bind, "kcal_per_mol", "", "NA",
            "Robetta_AlaScan", "In_Silico", "http://robetta.bakerlab.org/alaninescan", "In_Silico",
            "Thesis_OpenAccess", "21898653", "PMC3267954", "10.32657/10356/75835",
            "Zhang Y (DR-NTU 48015 Table 4.1)", "",
            f"int_id={int_id}. Note: 5 top-of-table Page 150 C3 rows (pos 17/19/21/79/82) had PDF multi-column interleaving in pdfkit text and are logged in follow-up CSV.")

    # =========================================================================
    # 4. PMC2323953 (Rhodobacter capsulatus BFR, Protein Sci 2003)
    # =========================================================================
    rc_buf = "50 mM KPi, pH 7.2; CD 222 nm"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Rhodobacter capsulatus", "Bacterioferritin (RcBFR)", "", "1JGC",
        "WT", "WT", "WT_reference", "Tm", "73.2", "73.2 ± 0.4 °C", "Exact_Numeric", "73.2", "0", "NA", "0.0", "NA",
        "CD_222nm", "Apo/As-isolated", rc_buf, "24-mer (heat-clarified 65 C)",
        "PMC_OpenAccess", "12876316", "PMC2323953", "10.1110/ps.0301903",
        "Kilic et al. 2003 Protein Sci", "", "WT baseline (heat-clarified at 65 C).")
    for mut, iface, tm_v, norm_d in [
        ("E135R", "C3/24-mer_interface", "59.2", "-14.0"),
        ("E128A/E135A", "C3/24-mer_interface", "47.1", "-26.1"),
        ("E128R/E135R", "C3/24-mer_interface", "43.3", "-29.9"),
    ]:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        ol = "24-mer" if mut == "E135R" else "Dimer (heat-clarified 55 C)"
        add("Mutant", FAM_BFR, EV_A, cat, "Rhodobacter capsulatus", "Bacterioferritin (RcBFR)", "", "",
            mut, mut, iface, "Tm", tm_v, f"{tm_v} ± 0.3 °C", "Exact_Numeric", "73.2", norm_d, "mut_minus_WT", norm_d,
            "Derived_mut_minus_WT", "CD_222nm", "Apo/As-isolated", rc_buf, ol,
            "PMC_OpenAccess", "12876316", "PMC2323953", "10.1110/ps.0301903",
            "Kilic et al. 2003 Protein Sci", "", "Derived delta Tm = Tm_mut - 73.2.")
    for mut in ("E128A", "E135A"):
        add("Mutant", FAM_BFR, EV_B, CAT_POINT, "Rhodobacter capsulatus", "Bacterioferritin (RcBFR)", "", "",
            mut, mut, "C3/24-mer_interface", "T_survival", "55.0", "Heat-clarified at 55 °C (CD Tm n.d.)",
            "Threshold_Survival", "65.0", "", "NA", "", "Range_Or_Qualitative",
            "Heat_Clarification_Solubility", "Apo/As-isolated", rc_buf, "Dimer",
            "PMC_OpenAccess", "12876316", "PMC2323953", "10.1110/ps.0301903",
            "Kilic et al. 2003 Protein Sci", "", "Purified via 55 C heat step (vs WT 65 C); CD Tm reported as n.d.")

    # =========================================================================
    # 5. PMC5729161 (Archaeoglobus fulgidus AfFtn, Biochemistry 2017)
    # =========================================================================
    af_buf = "800 mM NaCl (Table 1; 0.3 mg/mL CD or 0.5 mg/mL DSC)"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Archaeoglobus fulgidus", "Ferritin (AfFtn)", "O29424", "1S3Q/1SQ3",
        "WT", "WT", "WT_reference", "Tm", "84.0", "84 °C", "Exact_Numeric", "84.0", "0", "NA", "0.0", "NA",
        "CD_and_DSC", "Apo/As-isolated", af_buf, "24-mer (open tetrahedral pore)",
        "PMC_OpenAccess", "28682599", "PMC5729161", "10.1021/acs.biochem.7b00296",
        "Pulford et al. 2017 Biochemistry", "", "Table 1 baseline in 800 mM NaCl.")
    for mut, iface, tm_v, norm_d, meth, ol, note in [
        ("A127R", "Tetrahedral_pore", "85.0", "+1.0", "CD_222nm", "Mixed/salt-dependent", "Table 1 in 800 mM NaCl."),
        ("D138K", "Tetrahedral_pore", "84.0", "0.0", "CD_222nm", "Disassembles at low salt", "Table 1 in 800 mM NaCl."),
        ("E65R", "C3_like / intersubunit", "98.0", "+14.0", "DSC", ">90% 24-mer across all salt concentrations", "PDB:5V5K; >96 C exceeded CD range so measured by DSC."),
    ]:
        add("Mutant", FAM_FERRITIN, EV_A, CAT_POINT, "Archaeoglobus fulgidus", "Ferritin (AfFtn)", "O29424",
            "5V5K" if mut == "E65R" else "", mut, mut, iface, "Tm", tm_v, f"{tm_v.rstrip('.0')} °C", "Exact_Numeric",
            "84.0", norm_d, "mut_minus_WT", norm_d, "Derived_mut_minus_WT", meth, "Apo/As-isolated", af_buf, ol,
            "PMC_OpenAccess", "28682599", "PMC5729161", "10.1021/acs.biochem.7b00296",
            "Pulford et al. 2017 Biochemistry", "", note)

    # =========================================================================
    # 6. PMC6149950 (E. coli BFR B-site aromatics, Molecules 2017)
    # Sign convention: reported Delta Tm = WT - mutant!
    # =========================================================================
    m61_buf = "Phosphate buffer, pH ~7; CD 222 nm, 25-95 C"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "1BFR",
        "BFR WT", "WT", "WT_reference", "Tm", "69.9", "69.9 °C", "Exact_Numeric", "69.9", "0", "NA", "0.0", "NA",
        "CD_222nm", "Apo/As-isolated", m61_buf, "24-mer & Dimer",
        "PMC_OpenAccess", "29292762", "PMC6149950", "10.3390/molecules22122184",
        "Kamyshinsky et al. 2017 Molecules", "", "Baseline for Table 2 B-site Trp mutants.")
    for mut, val, repd, normd, ol in [
        ("D132W", "69.7", "0.2", "-0.2", "24-mer & Dimer"),
        ("N34W", "67.1", "2.8", "-2.8", "Dimer"),
        ("D132W/N34W", "57.8", "12.1", "-12.1", "Dimer (no cages in TEM)"),
    ]:
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_BFR, EV_A, cat, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "",
            mut, mut, "B_site_pore", "Tm", val, f"{val} °C", "Exact_Numeric", "69.9", repd, "WT_minus_mut_inverted", normd,
            "Sign_Inverted_From_Reported", "CD_222nm", "Apo/As-isolated", m61_buf, ol,
            "PMC_OpenAccess", "29292762", "PMC6149950", "10.3390/molecules22122184",
            "Kamyshinsky et al. 2017 Molecules", "", "Reported Delta Tm is WT-mut.")

    # =========================================================================
    # 7. PMC6800815 (Human H/L heteropolymers neuroferritinopathy, Metallomics 2019)
    # =========================================================================
    m68_buf = "50 mM phosphate, 50 mM NaCl, pH 7.0; 0.5-1.0 mg/mL; NanoDSC 3 atm, 25-120 C, 1 C/min"
    # WT baselines: wt H/L (~10-15% L) and HuHF homopolymer (keep all deconvolved peaks!)
    for mut_raw, mut_norm, peaks, star_val, dS, note in [
        ("wt H/L (~10-15% L)", "WT H/L (~10-15% L)", [("Tm1", "100.7"), ("Tm2", "101.8"), ("Tm3", "104.6")], "101.8", "18.5", "WT H/L heteropolymer control (~10-15% L); Tm2*=101.8 C main peak."),
        ("HuHF (homopolymer)", "HuHF WT", [("Tm1", "97.9"), ("Tm2", "100.7"), ("Tm3", "102.0"), ("Tm4", "104.7")], "102.0", "23.2", "H-homopolymer baseline; Tm3*=102.0 C main peak."),
    ]:
        for mname, val in peaks:
            is_star = (val == star_val)
            add("WT_baseline", FAM_FERRITIN if "homopolymer" in mut_raw else FAM_FUSION, EV_A, CAT_WT,
                "Homo sapiens", "Ferritin H/L heteropolymer" if "H/L" in mut_raw else "Ferritin heavy chain (HuHF)",
                "P02794/P02792", "2FHA", mut_raw, mut_norm, "WT_reference", mname, val,
                f"{val} °C" + (" (main peak*)" if is_star else "") + f" [ΔS={dS} kJ/mol·K]",
                "Exact_Numeric", star_val, "0", "NA", "0.0", "NA", "NanoDSC", "Apo", m68_buf, "24-mer",
                "PMC_OpenAccess", "31513212", "PMC6800815", "10.1039/c9mt00154a",
                "Pulsipher et al. 2019 Metallomics", "", note)

    neuro_muts = [
        ("L135P (H/L heteropolymer)", "FTL L135P (in H/L)", CAT_POINT, "P02792/P02794", "", "C3_pore",
         [("Tm1", "96.4"), ("Tm2", "98.7"), ("Tm3", "100.9")], "98.7", "-3.1", "14.5",
         "Non-pathogenic 3-fold L-chain substitution (~10-15% L); Tm2*=98.7 C (-3.1 vs wt H/L 101.8)."),
        ("Ln1 (FTL 460InsA / L154fs in H/L)", "FTL 460InsA (L154fs in H/L)", CAT_INDEL, "P02792/P02794", "", "E_helix_C4",
         [("Tm1", "84.8"), ("Tm2", "89.4"), ("Tm3", "94.2")], "89.4", "-12.4", "26.7",
         "Pathogenic neuroferritinopathy variant; Tm2*=89.4 C (-12.4 vs wt H/L 101.8)."),
        ("Ln2 (FTL 498InsTC / L167fs in H/L)", "FTL 498InsTC (L167fs / F167SfsX26 in H/L)", CAT_INDEL, "P02792/P02794", "4V6B", "E_helix_C4",
         [("Tm1", "91.8"), ("Tm2", "93.9"), ("Tm3", "97.5"), ("Tm4", "102.3")], "93.9/97.5", "-7.9", "30.5",
         "Two major peaks Tm2*=93.9 C (-7.9 vs 101.8) and Tm3*=97.5 C (-4.3 vs 101.8)."),
        ("Ln3 (FTL 442InsC / L148fs in H/L)", "FTL 442InsC (L148fs in H/L)", CAT_INDEL, "P02792/P02794", "", "E_helix_C4",
         [("Tm1", "78.5"), ("Tm2", "82.1"), ("Tm3", "85.6")], "85.6", "-16.2", "29.8",
         "Pathogenic frameshift; Tm3*=85.6 C (-16.2 vs wt H/L 101.8). Note holo (~200-300 Fe/shell) Tm within 5% of apo."),
    ]
    for mraw, mnorm, cat, up, pdb, iface, peaks, star_lbl, norm_d_main, dS, note in neuro_muts:
        for mname, val in peaks:
            is_star = val in star_lbl
            d_derived = f"{float(val) - 101.8:+.1f}" if is_star else ""
            add("Mutant", FAM_FUSION, EV_A, cat, "Homo sapiens", "Ferritin H/L heteropolymer (~10-15% L)",
                up, pdb, mraw, mnorm, iface, mname, val,
                f"{val} °C" + (" (main peak*)" if is_star else "") + f" [ΔS={dS} kJ/mol·K]",
                "Exact_Numeric", "101.8", d_derived, "mut_minus_WT" if d_derived else "NA", d_derived,
                "Derived_mut_minus_WT" if d_derived else "NA", "NanoDSC", "Apo", m68_buf, "24-mer (85-90% H, 10-15% mutant L)",
                "PMC_OpenAccess", "31513212", "PMC6800815", "10.1039/c9mt00154a",
                "Pulsipher et al. 2019 Metallomics", "", note)

    # =========================================================================
    # 8. PMC9069879 (Soybean SSFH-2 & HuHF-EP, RSC Adv 2019 - Freeze-dried DSC!)
    # =========================================================================
    m90_buf = "Dialyzed freeze-dried solid protein, PerkinElmer DSC 8000, 10 C/min"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Glycine max", "Soybean seed H-2 ferritin (SSFH-2)", "Q94IC4", "6J4J",
        "SSFH-2 WT", "WT", "WT_reference", "Tm", "106.0", "106 °C", "Exact_Numeric", "106.0", "0", "NA", "0.0", "NA",
        "DSC", "Apo/As-isolated (freeze-dried solid)", m90_buf, "24-mer",
        "PMC_OpenAccess", "35528680", "PMC9069879", "10.1039/c9ra04785a", "Zhang et al. 2019 RSC Adv", "",
        "Freeze-dried solid DSC! Heat-treated 100 C 5 min (TTH-2) = 108 C; >=10 min decreases; 121 C 5 min denatures.")
    add("Mutant", FAM_FERRITIN, EV_A, CAT_INDEL, "Glycine max", "Soybean seed H-2 ferritin (SSFH-2)", "Q94IC4", "",
        "SSFH-2-ΔEP", "SSFH-2 ΔEP (N-term extra peptide deleted)", "N_term_EP", "Tm", "88.0", "88 °C", "Exact_Numeric",
        "106.0", "-18.0", "mut_minus_WT", "-18.0", "Derived_mut_minus_WT",
        "DSC", "Apo/As-isolated (freeze-dried solid)", m90_buf, "24-mer",
        "PMC_OpenAccess", "35528680", "PMC9069879", "10.1039/c9ra04785a", "Zhang et al. 2019 RSC Adv", "",
        "Deleting ~30 aa plant extra peptide lowers freeze-dried DSC Tm by 18 C.")
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "2FHA",
        "HuHF WT", "WT", "WT_reference", "Tm", "82.0", "82 °C", "Exact_Numeric", "82.0", "0", "NA", "0.0", "NA",
        "DSC", "Apo/As-isolated (freeze-dried solid)", m90_buf, "24-mer",
        "PMC_OpenAccess", "35528680", "PMC9069879", "10.1039/c9ra04785a", "Zhang et al. 2019 RSC Adv", "",
        "Freeze-dried solid DSC value (lower than aqueous NanoDSC ~92-102 C).")
    add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Homo sapiens", "Human H-ferritin (HuHF) + plant EP", "P02794", "6J4A",
        "HuHF-∇EP", "HuHF-∇EP (~30 aa SSFH-2 N-term EP fused)", "N_term_EP_graft", "Tm", "103.0", "103 °C", "Exact_Numeric",
        "82.0", "+21.0", "mut_minus_WT", "+21.0", "Derived_mut_minus_WT",
        "DSC", "Apo/As-isolated (freeze-dried solid)", m90_buf, "24-mer",
        "PMC_OpenAccess", "35528680", "PMC9069879", "10.1039/c9ra04785a", "Zhang et al. 2019 RSC Adv", "",
        "Grafting ~30 aa plant EP onto HuHF N-term raises freeze-dried solid DSC Tm by +21 C (PDB:6J4A).")

    # =========================================================================
    # 9. PMC3403414 (Ulva pertusa UpFER, Protein Sci 2012 - 1.5 M GdnHCl DSC)
    # =========================================================================
    up_buf = "20 mM NaPi, pH 7.6 + 1.5 M GdnHCl; MCS-DSC"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Ulva pertusa", "Sea lettuce ferritin (UpFER)", "", "",
        "Native UpFER", "WT", "WT_reference", "Tm", "81.5", "81.5 ± 0.161 °C", "Exact_Numeric", "81.5", "0", "NA", "0.0", "NA",
        "MCS_DSC", "Apo/As-isolated", up_buf, "24-mer",
        "PMC_OpenAccess", "22419613", "PMC3403414", "10.1002/pro.2061", "Masuda et al. 2012 Protein Sci", "",
        "1.5 M GdnHCl included in DSC buffer to bring Tm below upper limit (paper also notes HuHF = 59.8 C in 2.0 M GdnHCl).")
    for mraw, mnorm, val, disp, dval in [
        ("delEP1 (ΔA1-T23)", "ΔA1-T23 (delEP1)", "77.8", "77.8 ± 0.0945 °C", "-3.7"),
        ("delEP2 (ΔA1-M7)", "ΔA1-M7 (delEP2)", "82.1", "82.1 ± 0.381 °C", "+0.6"),
    ]:
        add("Mutant", FAM_FERRITIN, EV_A, CAT_INDEL, "Ulva pertusa", "Sea lettuce ferritin (UpFER)", "", "",
            mraw, mnorm, "N_term_EP", "Tm", val, disp, "Exact_Numeric", "81.5", dval, "mut_minus_WT", dval, "Derived_mut_minus_WT",
            "MCS_DSC", "Apo/As-isolated", up_buf, "24-mer",
            "PMC_OpenAccess", "22419613", "PMC3403414", "10.1002/pro.2061", "Masuda et al. 2012 Protein Sci", "",
            "1.5 M GdnHCl present.")

    # =========================================================================
    # 10. PMC9785301 (Apostichopus japonicus AjFER, Polymers 2022)
    # =========================================================================
    aj_buf = "25 mM Tris-HCl pH 8.0, 150 mM NaCl; ~5 mg/mL; Tycho NT.6, 35-95 C at 30 C/min"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Apostichopus japonicus", "Sea cucumber ferritin (AjFER)", "", "",
        "AjFER WT", "WT", "WT_reference", "Ti", "88.7", "88.7 °C", "Exact_Numeric", "88.7", "0", "NA", "0.0", "NA",
        "Tycho_NT6", "Apo/As-isolated", aj_buf, "24-mer (90 C 10 min DLS DH=13.5 ± 1.0 nm)",
        "PMC_OpenAccess", "36559745", "PMC9785301", "10.3390/polym14245378", "Wang et al. 2022 Polymers", "",
        "Baseline inflection temperature Ti.")
    for mraw, mnorm, pdb, val, norm_d, dls_note in [
        ("MF (E25A/E60A/E105A)", "E25A/E60A/E105A (MF)", "", "69.8", "-18.9", "90 C 10 min DLS DH=141.8 ± 8.7 nm (aggregated)"),
        ("M3 (D129A/E132A)", "D129A/E132A (M3)", "7Y74", "87.4", "-1.3", "90 C 10 min DLS DH=190.1 ± 7.9 nm (aggregated)"),
        ("M4 (E168A)", "E168A (M4)", "", "74.7", "-14.0", "90 C 10 min DLS DH=122.4 ± 5.4 nm (aggregated)"),
    ]:
        cat = CAT_MULTI if "/" in mnorm else CAT_POINT
        add("Mutant", FAM_FERRITIN, EV_A, cat, "Apostichopus japonicus", "Sea cucumber ferritin (AjFER)", "", pdb,
            mraw, mnorm, "Ferroxidase / Pore_cavity", "Ti", val, f"{val} °C", "Exact_Numeric", "88.7", norm_d,
            "mut_minus_WT", norm_d, "Derived_mut_minus_WT", "Tycho_NT6", "Apo/As-isolated", aj_buf, dls_note,
            "PMC_OpenAccess", "36559745", "PMC9785301", "10.3390/polym14245378", "Wang et al. 2022 Polymers", "",
            "Derived delta Ti = Ti_mut - 88.7.")

    # =========================================================================
    # 11. PMC9798250 (Human H/L heteropolymers, Protein Sci 2023)
    # =========================================================================
    m97_buf = "50 mM phosphate, 50 mM NaCl, pH 7.4; 0.5 ± 0.1 mg/mL; NanoDSC 3 atm, 25-125 C"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human L-ferritin homopolymer (HuLF)", "P02792", "2FFX",
        "HuLF WT", "WT", "WT_reference", "Tm", "115.0", "115 ± 1 °C", "Exact_Numeric", "115.0", "0", "NA", "0.0", "NA",
        "NanoDSC", "Apo", m97_buf, "24-mer (L24)", "PMC_OpenAccess", "36519270", "PMC9798250", "10.1002/pro.4543",
        "Srivastava et al. 2023 Protein Sci", "", "Table 1 / abstract: iron-loading (~700 Fe) lowers L/L-rich Tm by -12 to -14 C.")
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human H-ferritin homopolymer (HuHF)", "P02794", "2FHA",
        "HuHF WT", "WT", "WT_reference", "Tm", "92.5", "92 °C (Table 1) / 93 ± 1 °C (abstract)", "Exact_Numeric", "92.5", "0", "NA", "0.0", "NA",
        "NanoDSC", "Apo", m97_buf, "24-mer (H24)", "PMC_OpenAccess", "36519270", "PMC9798250", "10.1002/pro.4543",
        "Srivastava et al. 2023 Protein Sci", "", "Iron-loading (~700 Fe) lowers H/H-rich Tm by -18 to -20 C.")
    for h_l, tm_str, tm_num, norm_vs_h, note in [
        ("H2:L22 heteropolymer", "~115 °C (L-rich regime)", "", "~+22.5 vs HuHF", "Co-expressed plasmid H/L heteropolymer; additional intermediate ratios plotted in Fig 3b (figure-only gap)."),
        ("H6:L18 heteropolymer", "~115 °C (L-rich regime)", "", "~+22.5 vs HuHF", "Co-expressed plasmid H/L heteropolymer; retains L-like hyperthermostability."),
        ("H19:L5 heteropolymer", "90-95 °C (H-rich regime); holo CD ~10 C lower", "", "0 to +2 vs HuHF", "H-rich heteropolymer; exact intermediate points in Fig 3b (figure-only gap)."),
    ]:
        add("Mutant", FAM_FUSION, EV_C if "90-95" in tm_str else EV_A, CAT_FUSION, "Homo sapiens", "Human H/L heteropolymer ferritin",
            "P02794/P02792", "", h_l, h_l, "H/L_subunit_stoichiometry", "Tm", "115.0" if tm_str.startswith("~115") else "",
            tm_str, "Approx_From_Curve_or_Text" if "~115" in tm_str else "Range_Group", "92.5", "", "mut_minus_WT",
            "+22.5" if "~115" in tm_str else "", "Range_Or_Qualitative", "NanoDSC", "Apo", m97_buf, "24-mer heteropolymer",
            "PMC_OpenAccess", "36519270", "PMC9798250", "10.1002/pro.4543", "Srivastava et al. 2023 Protein Sci", "", note)
    # Literature citations inside PMC9798250 Table 1 / Discussion (non-duplicate new mutant lead)
    add("Mutant", FAM_FERRITIN, EV_A, CAT_POINT, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "",
        "E134K", "E134K", "C3_pore", "delta_Tm", "", "+5 °C vs WT (cited from Tan et al. 2021)", "Approx_From_Curve_or_Text",
        "", "+5 °C", "mut_minus_WT", "+5.0", "Directly_Reported", "Lit_Citation_PMC9798250", "Unknown",
        "Cited in PMC9798250 Discussion (Tan et al. 2021)", "24-mer", "PMC_OpenAccess", "36519270", "PMC9798250", "10.1002/pro.4543",
        "Srivastava et al. 2023 citing Tan et al. 2021", "", "Secondary citation inside PMC9798250; original Tan et al. 2021 paper flagged for primary verification.")
    add("Mutant", FAM_FERRITIN, EV_A, CAT_MULTI, "Marsupenaeus japonicus", "Kuruma prawn ferritin (MjFer)", "T2B7E1", "",
        "3-fold salt-bridge mutant (Tm=91 C vs WT 109 C)", "MjFer 3-fold salt-bridge mutant", "C3_pore", "Tm", "91.0",
        "91 °C (vs MjFer WT 109 °C)", "Exact_Numeric", "109.0", "-18.0", "mut_minus_WT", "-18.0", "Derived_mut_minus_WT",
        "Lit_Citation_PMC9798250", "Unknown", "PMC9798250 Table 1 literature benchmark", "24-mer",
        "PMC_OpenAccess", "36519270", "PMC9798250", "10.1002/pro.4543", "Srivastava et al. 2023 Table 1", "",
        "Secondary citation in PMC9798250 Table 1 (exact residue labels need primary paper verification).")

    # =========================================================================
    # 12. PMC11776510 (M. tuberculosis BfrA, Chem Sci 2025 - ESI Fig S13 approx)
    # =========================================================================
    bfrA_buf = "10 mM Tris, pH 7.0; 6 uM subunit; CD 222 nm, 20-95 C, 1 C/min (RSC ESI Fig S13)"
    add("WT_baseline", FAM_BFR, EV_A, CAT_WT, "Mycobacterium tuberculosis", "Bacterioferritin A (MtBfrA)", "", "",
        "BfrA WT", "WT", "WT_reference", "Tm", "90.0", "~90 °C (ESI Fig S13 approx)", "Approx_From_Curve_or_Text",
        "90.0", "0", "NA", "0.0", "NA", "CD_222nm", "Apo/As-isolated", bfrA_buf, "24-mer",
        "PMC_OpenAccess", "39886445", "PMC11776510", "10.1039/d4sc07021f", "2025 Chem Sci (PMC11776510)", "",
        "Approximate midpoint from ESI Fig S13 prior reading; exact decimal Tm requires external RSC ESI download.")
    for mut, val, norm_d, ol in [
        ("E121A", "90.0", "0.0", "24-mer (self-assembled)"),
        ("E121Q", "84.0", "-6.0", "24-mer (self-assembled)"),
        ("E121K", "78.0", "-12.0", "Assembly-defective (low-T overexpr tested)"),
        ("E121F", "80.0", "-10.0", "Assembly-defective (low-T overexpr tested)"),
    ]:
        add("Mutant", FAM_BFR, EV_A, CAT_POINT, "Mycobacterium tuberculosis", "Bacterioferritin A (MtBfrA)", "", "",
            mut, mut, "C3_pore", "Tm", val, f"~{val.rstrip('.0')} °C (ESI Fig S13 approx)", "Approx_From_Curve_or_Text",
            "90.0", norm_d, "mut_minus_WT", norm_d, "Derived_mut_minus_WT", "CD_222nm", "Apo/As-isolated", bfrA_buf, ol,
            "PMC_OpenAccess", "39886445", "PMC11776510", "10.1039/d4sc07021f", "2025 Chem Sci (PMC11776510)", "",
            "Approximate curve values from ESI Fig S13 (exact decimals flagged for RSC ESI follow-up).")

    # =========================================================================
    # 13. PMC13406861 (HuHF Trp lumen variants, Curr Issues Mol Biol 2026)
    # =========================================================================
    w_buf = "20 mM NaPi pH 7.4 + 150 mM NaCl; 5 uM protein; Rayleigh scattering"
    add("Mutant", FAM_FERRITIN, EV_A, CAT_MULTI, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "",
        "HFt-W4 (H58W/H61W/H129W/H152W)", "H58W/H61W/H129W/H152W (HFt-W4)", "Inner_cavity_lumen", "Tagg_onset",
        "70.0", "~70 °C (remains alpha-helical after 5 h at 75 C)", "Approx_From_Curve_or_Text", "", "", "NA", "", "Range_Or_Qualitative",
        "Rayleigh_Scattering", "Apo/As-isolated", w_buf, "24-mer",
        "PMC_OpenAccess", "42510950", "PMC13406861", "10.3390/cimb48070710", "2026 Curr Issues Mol Biol (PMC13406861)", "",
        "Designed Trp lumen nanocage; no WT Tagg numeric value quoted alongside in same sentence.")
    add("Mutant", FAM_FERRITIN, EV_A, CAT_MULTI, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "",
        "HFt-W6 (H58W/H61W/H129W/H152W/Y138W/F171W)", "H58W/H61W/H129W/H152W/Y138W/F171W (HFt-W6)", "Inner_cavity_lumen", "Tagg_onset",
        "48.0", "~48 °C (converts to beta-sheet at 58-75 C)", "Approx_From_Curve_or_Text", "70.0", "-22.0 (vs HFt-W4)", "mut_minus_WT", "-22.0",
        "Derived_mut_minus_WT", "Rayleigh_Scattering", "Apo/As-isolated", w_buf, "24-mer -> beta-sheet aggregate >=58 C",
        "PMC_OpenAccess", "42510950", "PMC13406861", "10.3390/cimb48070710", "2026 Curr Issues Mol Biol (PMC13406861)", "",
        "Delta calculated vs HFt-W4 reference (~48 vs ~70 C = -22 C).")

    # =========================================================================
    # 14. PMC152977 (Bullfrog M ferritin pore localized unfolding, PNAS 2003)
    # =========================================================================
    bf_buf = "0.01 M MOPS, pH 7.0; CD (pore subdomain) & A280 (global)"
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Rana catesbeiana", "Bullfrog M-ferritin (H')", "P07798", "",
        "H-WT", "WT", "WT_reference", "T_sub_CD", "53.0", "53 °C (pore CD; +1 mM urea -> 43 °C; global A280 >85 °C)",
        "Exact_Numeric", "53.0", "0", "NA", "0.0", "NA", "CD_222nm", "Apo", bf_buf, "24-mer",
        "PMC_OpenAccess", "12634425", "PMC152977", "10.1073/pnas.0636928100", "Takagi et al. 2003 PNAS", "",
        "Localized 3-fold pore subdomain unfolds ~30 C below global cage denaturation (>85 C by A280).")
    add("Mutant", FAM_FERRITIN, EV_B, CAT_POINT, "Rana catesbeiana", "Bullfrog M-ferritin (H')", "P07798", "3SHX",
        "H-L134P", "L134P", "C3_pore", "T_sub_CD", "", "<43 °C pore CD (loss 45-65 C; global A280 <85 °C)", "Threshold_Survival",
        "53.0", "lower than WT+1 mM urea (<43 C)", "mut_minus_WT", "", "Range_Or_Qualitative",
        "CD_and_A280", "Apo", bf_buf, "24-mer with unfolded 3-fold pore",
        "PMC_OpenAccess", "12634425", "PMC152977", "10.1073/pnas.0636928100", "Takagi et al. 2003 PNAS", "",
        "Destabilizes both localized pore subdomain (<43 C) and lowers global A280 stability (<85 C vs >85 C).")
    for mut in ("L134V", "R72D/D122R"):
        cat = CAT_MULTI if "/" in mut else CAT_POINT
        add("Mutant", FAM_FERRITIN, EV_B, cat, "Rana catesbeiana", "Bullfrog M-ferritin (H')", "P07798", "",
            f"H-{mut}", mut, "C3_pore", "T_survival", "85.0", ">85 °C global A280 (unchanged vs WT)", "Threshold_Survival",
            "", "~0 (unchanged >85 C)", "NA", "0.0", "Range_Or_Qualitative",
            "A280_Thermal", "Apo", bf_buf, "24-mer",
            "PMC_OpenAccess", "12634425", "PMC152977", "10.1073/pnas.0636928100", "Takagi et al. 2003 PNAS", "",
            "Maintains WT global thermal threshold (>85 C).")

    # =========================================================================
    # 15. PMC10093968 (H. pylori G27 HpFtn + NadA fusions, IJMS 2023)
    # =========================================================================
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Helicobacter pylori (G27)", "Non-heme ferritin (HpFtn 5-167)", "", "",
        "HpFtn (5-167) alone", "HpFtn (5-167) N-trunc baseline", "WT_reference", "Tm", "100.0", "100 °C", "Exact_Numeric",
        "100.0", "0", "NA", "0.0", "NA", "DSC", "Apo/As-isolated", "PMC10093968 buffer", "24-mer",
        "PMC_OpenAccess", "37047152", "PMC10093968", "10.3390/ijms24076183", "2023 IJMS (PMC10093968)", "",
        "Carrier baseline for NadA-HpFtn fusions; uses HpFtn residues 5-167 + SGGAGS linker.")
    add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Helicobacter pylori / Neisseria meningitidis", "HpFtn(5-167)-SGGAGS-NadA_Sht", "", "",
        "NadA_Sht-Fe", "HpFtn(5-167)-SGGAGS-NadA_Sht", "N_term_antigen_fusion", "Tm", "87.7", "87.7 °C", "Exact_Numeric",
        "100.0", "-12.3 (vs HpFtn); +42.7 (vs free NadA 45 C)", "mut_minus_WT", "-12.3", "Derived_mut_minus_WT",
        "DSF_SYPRO", "As-isolated", "DSF in PMC10093968", "24-mer nanoparticle",
        "PMC_OpenAccess", "37047152", "PMC10093968", "10.3390/ijms24076183", "2023 IJMS (PMC10093968)", "",
        "Free NadA(24-342) Tm=45 C; HpFtn carrier Tm=100 C.")
    for mname, val in [("Tm1", "47.2"), ("Tm2", "56.0")]:
        add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Helicobacter pylori / Neisseria meningitidis", "HpFtn(5-167)-SGGAGS-NadA_Ext", "", "",
            "NadA_Ext-Fe", "HpFtn(5-167)-SGGAGS-NadA_Ext", "N_term_antigen_fusion", mname, val, f"{val} °C", "Exact_Numeric",
            "100.0", f"{float(val)-100.0:+.1f} (vs HpFtn)", "mut_minus_WT", f"{float(val)-100.0:+.1f}", "Derived_mut_minus_WT",
            "DSF_SYPRO", "As-isolated", "DSF in PMC10093968", "24-mer nanoparticle",
            "PMC_OpenAccess", "37047152", "PMC10093968", "10.3390/ijms24076183", "2023 IJMS (PMC10093968)", "",
            "Bimodal DSF unfolding of extended NadA antigen-ferritin fusion.")

    # =========================================================================
    # 16. PMC7746619 (Pyrococcus furiosus ferritin Trx fusion, Front Immunol 2020)
    # =========================================================================
    add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Pyrococcus furiosus", "Pf-Fe-Trx L2(8mer)", "", "",
        "Pf-Fe-Trx L2(8mer)", "Pf-Fe-Trx L2(8mer)", "Surface_Trx_fusion", "Tm", "89.0", "89 °C", "Exact_Numeric",
        "", "", "NA", "", "NA", "CD_222nm", "As-isolated", "0.1X PBS; 1 uM; 25-100 C, 10 min hold/step", "24-mer nanoparticle",
        "PMC_OpenAccess", "33343580", "PMC7746619", "10.3389/fimmu.2020.606569", "2020 Front Immunol (PMC7746619)", "",
        "CD Tm of L2(8mer) antigen-displaying PfFtn fusion. PfFe/3mer/8mer/Trx3mer/Trx8mer also survived 65 C 30 min.")
    for tag_v in ("Pf-Fe-3mer", "Pf-Fe-8mer", "Pf-Fe-Trx3mer", "Pf-Fe-Trx8mer"):
        add("Mutant", FAM_FUSION, EV_B, CAT_FUSION, "Pyrococcus furiosus", f"{tag_v} fusion", "", "",
            tag_v, tag_v, "Surface_Trx_fusion", "T_survival", "65.0", ">=65 °C (30 min heat-purification survival)",
            "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative", "Heat_Clarification_Solubility",
            "As-isolated", "0.1X PBS / purification lysate 65 C 30 min", "24-mer nanoparticle",
            "PMC_OpenAccess", "33343580", "PMC7746619", "10.3389/fimmu.2020.606569", "2020 Front Immunol (PMC7746619)", "",
            "Survived 65 C 30 min heat treatment.")

    # =========================================================================
    # 17. ProThermDB (PMID:33825471 & PMID:9681870) + FireProtDB (PMID:9681870)
    # =========================================================================
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "2FHA",
        "WT (ProThermDB NO=169)", "WT", "WT_reference", "Tm1", "80.1", "80.1 °C", "Exact_Numeric", "80.1", "0", "NA", "0.0", "NA",
        "DSC", "As-isolated", "20 mM Phosphate buffer, pH 7.0", "24-mer",
        "StructuredDB_ProThermDB", "33825471", "", "10.1021/acs.jpcb.1c00276", "Qu et al. 2021 J Phys Chem B (ProThermDB #169)",
        "ProThermDB_vs_PMID33825471", "WT Tm1 baseline in ProThermDB (pt2_prot_ferritin.tsv).")
    add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "2FHA",
        "WT (ProThermDB NO=172)", "WT", "WT_reference", "Tm2", "101.3", "101.3 °C", "Exact_Numeric", "101.3", "0", "NA", "0.0", "NA",
        "DSC", "As-isolated", "20 mM Phosphate buffer, pH 7.0", "24-mer",
        "StructuredDB_ProThermDB", "33825471", "", "10.1021/acs.jpcb.1c00276", "Qu et al. 2021 J Phys Chem B (ProThermDB #172)",
        "ProThermDB_vs_PMID33825471", "WT Tm2 baseline in ProThermDB (pt2_prot_ferritin.tsv).")

    for pt_no, mraw, mnorm, mname, val, wt_ref, normd, term in [
        ("170", "E1F1 (EBNA1 N-term HuHF fusion)", "EBNA1-HuHF N-term fusion (E1F1)", "Tm1", "70.5", "80.1", "-9.6", "N-terminus"),
        ("173", "E1F1 (EBNA1 N-term HuHF fusion)", "EBNA1-HuHF N-term fusion (E1F1)", "Tm2", "80.1", "101.3", "-21.2", "N-terminus"),
        ("171", "F1E1 (EBNA1 C-term HuHF fusion)", "HuHF-EBNA1 C-term fusion (F1E1)", "Tm1", "63.6", "80.1", "-16.5", "C-terminus"),
        ("174", "F1E1 (EBNA1 C-term HuHF fusion)", "HuHF-EBNA1 C-term fusion (F1E1)", "Tm2", "84.6", "101.3", "-16.7", "C-terminus"),
    ]:
        add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Homo sapiens / EBV", "Human H-ferritin + EBNA1 epitope", "P02794", "",
            mraw, mnorm, f"{term}_fusion", mname, val, f"{val} °C", "Exact_Numeric", wt_ref, normd, "mut_minus_WT", normd,
            "Derived_mut_minus_WT", "DSC", "As-isolated", "20 mM Phosphate buffer, pH 7.0", "24-mer nanovaccine",
            "StructuredDB_ProThermDB", "33825471", "", "10.1021/acs.jpcb.1c00276",
            f"Qu et al. 2021 J Phys Chem B (ProThermDB #{pt_no})", "ProThermDB_vs_PMID33825471",
            f"ProThermDB pt2_prot_ferr_low.tsv mistakenly tagged MUTATION='wild-type'; REMARKS reveals EBNA1 {term} fusion. Abstract notes F1E1 shows higher overall heat/pH/denaturant tolerance than E1F1.")

    # FireProtDB (seq_id=8955 P02792) / ProThermDB (#9017-9023) 7 WT subunit DSC baselines
    for pt_no, fp_id, ph, val, dH, dHvH in [
        ("9017", "67854", "2.0", "37.2", "70.3", "41.7"),
        ("9018", "67858", "2.2", "53.1", "86.3", "86.9"),
        ("9019", "67862", "2.4", "65.1", "103.2", "105.4"),
        ("9020", "67866", "2.8", "79.8", "105.5", "145.1"),
        ("9021", "67870", "2.0", "50.0", "40.4", "40.4"),
        ("9022", "67874", "2.4", "50.3", "40.6", "44.1"),
        ("9023", "67878", "2.8", "50.4", "41.1", "40.3"),
    ]:
        add("WT_baseline", FAM_FERRITIN, EV_A, CAT_WT, "Homo sapiens", "Human ferritin L/H WT subunit (FireProt seq_id=8955 / P02792)",
            "P02792", "2FHA", f"WT pH {ph} (#{pt_no})", f"WT_pH_{ph}_#{pt_no}", "WT_acid_unfolding", "Tm", val,
            f"{val} °C (ΔH={dH}, ΔHvH={dHvH} kcal/mol)", "Exact_Numeric", "", "", "NA", "", "NA", "DSC", "Apo",
            f"0.05 M Sodium phosphate-HCl, pH {ph}; 10 ug/mL", "Subunit/Dissociated",
            "StructuredDB_FireProtDB|StructuredDB_ProThermDB", "9681870", "", "10.1093/protein/11.5.377",
            f"Martsev et al. 1998 Protein Eng (FireProt {fp_id}/ProTherm {pt_no})", "",
            "WT acid-pH baseline only; FireProtDB mutant_id IS NULL (0 mutant rows in FireProtDB).")

    # =========================================================================
    # 18. Closed-Paper Abstract Quantitative & Qualitative Records
    # =========================================================================
    add("Mutant", FAM_FERRITIN, EV_A, CAT_POINT, "Crassostrea gigas", "Pacific oyster ferritin", "Q70MM3", "6LJG",
        "D120G (Asp120->Gly120; PDB 6LJG D119G)", "D120G", "C3_interface", "Tm", "89.17",
        "89.17 °C (morphology/size unchanged <90 C)", "Exact_Numeric", "", "+13 °C", "mut_minus_WT", "+13.0",
        "Directly_Reported", "DSC/CD_Abstract", "As-isolated", "Food Chem 2021 abstract (closed DOI)", "24-mer",
        "PubMed_Abstract_ClosedPaper", "33406454", "", "10.1016/j.foodchem.2020.128879",
        "Food Chem 2021 (PMID:33406454)", "Oyster_D120G_6LJG",
        "Implies WT Tm ~76.17 C; closed Elsevier PDF flagged for exact WT value and assay method.")

    add("Mutant", FAM_FUSION, EV_A, CAT_FUSION, "Homo sapiens", "EGF-FTH1 / eGFP-(G4S)3-FTH1 heteromeric nanocage", "P02794", "",
        "40-60% EGF-FTH1 + eGFP-(G4S)3-FTH1 heteromer", "EGF-FTH1 / eGFP-(G4S)3-FTH1 heteromer", "C3_interface", "Tdiss",
        "79.2", "62.9 -> 79.2 °C (across increasing EGF-FTH1 ratio)", "Range_Group", "", "+16.3 °C (max vs min eGFP load)",
        "mut_minus_WT", "+16.3", "Derived_mut_minus_WT", "nanoDSF_and_qPCR_DSF", "As-isolated", "Bioconjug Chem 2026 abstract (closed DOI)",
        "24-mer heteromeric cage", "PubMed_Abstract_ClosedPaper", "42902949", "", "10.1021/acs.bioconjchem.6c00378",
        "Bioconjug Chem 2026 (PMID:42902949)", "",
        "Tdiss of intact 24-mer rises from 62.9 to 79.2 C as EGF-FTH1 proportion increases; exact per-ratio table behind ACS paywall.")

    add("Mutant", FAM_DPS, EV_A, CAT_MULTI, "Escherichia coli", "Dps hydrophobic-cavity variants (up to 120 muts/cage)", "P0ABT2", "",
        "Hydrophobic-cavity mutant set (up to 10 muts/subunit)", "EcDps hydrophobic-cavity mutant set", "Interior_cavity", "T_range_all_mutants",
        "", "74-90 °C (all cavity mutants)", "Range_Group", "", "", "NA", "", "Range_Or_Qualitative",
        "CD_and_Trp_Fluorescence", "As-isolated", "JACS 2006 abstract (closed DOI)", "12-mer",
        "PubMed_Abstract_ClosedPaper", "16704261", "", "10.1021/ja057069x",
        "Swift et al. 2006 JACS (PMID:16704261)", "",
        "Closed ACS paper: individual mutant sequences and per-mutant Tm within 74-90 C flagged for follow-up.")

    add("Mutant", FAM_BFR, EV_A, CAT_MULTI, "Escherichia coli (dimorphic ferritin/Bfr)", "Computational C3-interface redesigned nanocage", "P0ABD3", "",
        "C3-interface monomorphic cage redesign", "EcBfr C3-interface monomorphic redesign", "C3_interface", "delta_Tm",
        "", ">= +20 °C vs parent", "Threshold_Survival", "", ">= +20 °C", "mut_minus_WT", "+20.0", "Directly_Reported",
        "Thermal_Melt_Abstract", "As-isolated", "Chem Commun 2013 abstract (closed DOI)", "Monomorphic 24-mer",
        "PubMed_Abstract_ClosedPaper", "23511498", "", "10.1039/c3cc40886h",
        "Ardejani et al. 2013 Chem Commun (PMID:23511498)", "",
        "Closed RSC paper: exact mutation identities and Tm values flagged for follow-up.")

    # Bullfrog M-ferritin C4 pore gradual mutants (PMID:41700891 + RCSB 9WX5..9X6R)
    for mnorm, pdb in [
        ("E130A/M161E", "9WX5"), ("E130A/L165D", "9WY3"), ("E130A/M161E/L165D", "9X0V"),
        ("E130A/L165D/K168E/H169D", "9X15"), ("E130A/M161Q/L165D/K168E/H169D", "9X16"),
        ("M161Q/L165D/K168E/H169D", "9X1G"), ("E130A/K168E/H169D", "9X6R"),
    ]:
        cat = CAT_MULTI if "/" in mnorm else CAT_POINT
        add("Mutant", FAM_FERRITIN, EV_A, cat, "Rana catesbeiana", "Bullfrog M-ferritin", "P07798", pdb,
            mnorm, mnorm, "C4_pore", "delta_Tm", "", "~ -5 °C (4-fold pore gradual series)", "Approx_From_Curve_or_Text",
            "", "~ -5 °C", "mut_minus_WT", "-5.0", "Range_Or_Qualitative",
            "Thermal_Abstract", "As-isolated", "Biochemistry 2026 abstract (closed DOI)", "24-mer (widened C4 pore)",
            "PubMed_Abstract_ClosedPaper", "41700891", "", "10.1021/acs.biochem.5c00744",
            "Biochemistry 2026 (PMID:41700891)", "",
            "Abstract reports modifications across this series decreased thermal stability ~5 C; per-variant Tm behind ACS paywall.")

    # =========================================================================
    # 19. Tier B Wet-Lab Thermal Survival Thresholds (Ferritin & Bfr & Fusion)
    # =========================================================================
    for mut, cat, t_val, disp, note in [
        ("HFn(+)", CAT_MULTI, "72.0", ">=72 °C (20 min heat-clarification survival)", "Supercharged HFn variant survives 72 C 20 min."),
        ("tHFn(+) 1-174", CAT_INDEL, "72.0", ">=72 °C (20 min survival)", "C-term truncated to aa 174; intact thermal tolerance."),
        ("tHFn(+) 1-162", CAT_INDEL, "72.0", ">=72 °C (20 min survival)", "C-term truncated to aa 162; still survives 72 C 20 min."),
        ("tHFn(+) 1-159", CAT_INDEL, "60.0", ">=60 °C (10 min survival; lowered vs 72 C)", "C-term truncated to aa 159; heat-clarification step reduced to 60 C 10 min."),
    ]:
        add("Mutant", FAM_FERRITIN, EV_B, cat, "Homo sapiens", "Human H-ferritin (HFn(+))", "P02794", "",
            mut, mut, "C_term_E_helix", "T_survival", t_val, disp, "Threshold_Survival",
            "72.0", "", "NA", "0.0" if t_val == "72.0" else "-12.0", "Range_Or_Qualitative",
            "Heat_Clarification_Solubility", "As-isolated", "Heat-clarification purification step", "24-mer",
            "PMC_OpenAccess", "39970222", "PMC11838010", "10.1016/j.apsb.2024.12.011", "2025 PMC11838010", "",
            note + " Shorter truncations 1-156, 1-153, 1-150 purified without heating.")

    add("Mutant", FAM_FERRITIN, EV_B, CAT_MULTI, "Thermotoga maritima", "Thermotoga maritima ferritin (TmFtn)", "Q9X0L2", "7DY8/7DY9",
        "TmFtn FLAL (E112L/E113A/K114L/N147F)", "E112L/E113A/K114L/N147F (FLAL)", "Exterior_interface", "T_survival",
        "90.0", ">=90 °C (10 min survival)", "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative",
        "Heat_Clarification_Solubility", "As-isolated", "90 C 10 min lysate heating", "1D/2D nanocage assemblies",
        "PMC_OpenAccess", "34381032", "PMC8357837", "10.1038/s41467-021-25199-x", "2021 Nat Commun (PMC8357837)", "",
        "Thermal survival during purification at 90 C 10 min.")

    for mraw, mnorm, prot in [
        ("HuHF Mutant A (Q14A/D15A/R22A)", "Q14A/D15A/R22A (Mutant A)", "Human H-ferritin (HuHF)"),
        ("HuHF Mutant B (F81A/Q83A)", "F81A/Q83A (Mutant B)", "Human H-ferritin (HuHF)"),
        ("HuHF Mutant C (Q14A/D15A/R22A/F81A/Q83A)", "Q14A/D15A/R22A/F81A/Q83A (Mutant C)", "Human H-ferritin (HuHF)"),
        ("HuLF Mutant D (S5T/T14Q/S22R/L81F/D116E/A119K/A123D)", "S5T/T14Q/S22R/L81F/D116E/A119K/A123D (Mutant D)", "Human L-ferritin (HuLF)"),
    ]:
        up = "P02792" if "HuLF" in prot else "P02794"
        add("Mutant", FAM_FERRITIN, EV_B, CAT_MULTI, "Homo sapiens", prot, up, "",
            mraw, mnorm, "Outer_surface_B_pore", "T_survival", "55.0", ">=55 °C (8 min heat-clarification survival)",
            "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative", "Heat_Clarification_Solubility",
            "As-isolated", "55 C 8 min heating", "24-mer",
            "PMC_OpenAccess", "30850661", "PMC6408514", "10.1038/s41467-019-09073-5", "2019 Nat Commun (PMC6408514)", "",
            "Purified with 55 C 8 min heat step.")

    for mraw, mnorm, disp, meth in [
        ("HuHF Mutant A (L56H/R63H/E67H)", "L56H/R63H/E67H (Mutant A)", ">=60 °C (10 min heat-clarification)", "Heat_Clarification_Solubility"),
        ("HuHF Mutant B (A52H/L56H/E67H/K71H)", "A52H/L56H/E67H/K71H (Mutant B)", ">=60 °C (10 min heat-clarification + stable at 60 C DLS)", "DLS_and_Heat_Clarification"),
    ]:
        add("Mutant", FAM_FERRITIN, EV_B, CAT_MULTI, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "",
            mraw, mnorm, "B_helix_His_patch", "T_survival", "60.0", disp, "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative",
            meth, "As-isolated", "60 C 10 min heating (+ DLS at 60 C for Mutant B)", "24-mer",
            "PMC_OpenAccess", "29271453", "PMC5812362", "10.1186/s12951-018-0334-5", "2018 J Nanobiotechnol (PMC5812362)", "",
            "Mutant C (N139H/K143H/K146H/D150H) also constructed in study.")

    add("Mutant", FAM_FERRITIN, EV_B, CAT_POINT, "Pyrococcus furiosus", "Pyrococcus furiosus ferritin (PfFtn)", "", "",
        "R117A-PfFtn", "R117A", "Channel/cavity", "T_survival", "100.0",
        "Thermostable at 100 °C (0-70 min timecourse, lower stability than WT-PfFtn)", "Threshold_Survival",
        "", "< WT-PfFtn at 100 C", "NA", "", "Range_Or_Qualitative",
        "Heat_Treatment_SDS_SEC", "As-isolated", "100 C for 0-70 min", "24-mer",
        "PMC_OpenAccess", "26370076", "PMC4646333", "10.1038/srep13924", "2015 Sci Rep (PMC4646333)", "",
        "Survives 100 C heat treatment with lower kinetic resistance than WT.")

    for ph, t_ok, t_bad in [("4.0", "70.0", "75.0"), ("7.0", "65.0", "70.0"), ("9.0", "60.0", "65.0")]:
        add("Mutant", FAM_FERRITIN, EV_B, CAT_POINT, "Marsupenaeus japonicus", "Kuruma prawn ferritin (MjFer)", "T2B7E1", "",
            f"T158H MjFer 3D crystal (pH {ph})", "T158H", "C4_His_crystal_lattice", "T_crystal_limit",
            t_ok, f"Intact {t_ok} °C / damaged at {t_bad} °C (15 min)", "Threshold_Survival", "", "", "NA", "",
            "Range_Or_Qualitative", "Crystal_Heating", "As-isolated_3D_Crystal", f"pH {ph}, 15 min heating",
            "3D crystalline nanocage array", "PMC_OpenAccess", "36703383", "PMC9814774", "10.1038/s42004-020-00394-x",
            "Tan et al. 2020 Commun Chem (PMC9814774)", "", "Crystalline thermal tolerance threshold.")

    add("Mutant", FAM_FERRITIN, EV_B, CAT_INDEL, "Homo sapiens", "Human H-ferritin (HFn)", "P02794", "",
        "Ecut-HFn (E-helix cut)", "Ecut-HFn (E-helix deleted)", "E_helix_C4", "T_survival", "45.0",
        "Heat-sensitive: Cu2+-assisted Dox loading at 45 °C", "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative",
        "Thermal_Drug_Loading", "As-isolated (+Cu2+)", "45 C with Cu2+ assistance", "24-mer",
        "PubMed_Abstract_ClosedPaper", "37729988", "", "10.1016/j.ijbiomac.2023.126973",
        "2023 Int J Biol Macromol (PMID:37729988)", "", "Heat-sensitized mutant enabling mild 45 C drug encapsulation.")

    add("Mutant", FAM_FERRITIN, EV_B, CAT_POINT, "Homo sapiens", "Human H-ferritin (HFn)", "P02794", "",
        "HFn D44A (and RGD-D44A)", "D44A", "C2_channel_R43_D45", "T_survival", "60.0",
        "Stable intact cage at 60 °C (30 min thermoresponsive Dox loading ~94/cage)", "Threshold_Survival", "", "", "NA", "",
        "Range_Or_Qualitative", "Heat_Treatment_SDS_SEC", "As-isolated", "60 C, 30 min", "24-mer (transiently widened C2 channel)",
        "PubMed_Abstract_ClosedPaper", "41638275", "", "10.1016/j.ijbiomac.2026.150643",
        "2026 Int J Biol Macromol (PMID:41638275)", "", "Disrupts R43-D45 H-bond network while preserving overall thermal stability at 60 C.")

    add("Mutant", FAM_FUSION, EV_B, CAT_FUSION, "Zea mays", "Maize ferritin 1 (ZmFer1)", "", "",
        "LTVS-L-ZmFer1 (16-aa linker + LTVSPWY anti-HER2)", "LTVS-L-ZmFer1 (N-term HER2 fusion)", "N_term_EP", "T_survival",
        "100.0", "Withstands up to 100 °C with only partial denaturation", "Threshold_Survival", "", "", "NA", "",
        "Range_Or_Qualitative", "Heat_Treatment_SDS_SEC", "As-isolated", "Heat treatment up to 100 C", "24-mer (19.01 nm DLS)",
        "PubMed_Abstract_ClosedPaper", "38904710", "", "10.1007/s11033-024-09636-w",
        "2024 Mol Biol Rep (PMID:38904710)", "", "HER2-targeted peptide fusion retaining hyperthermostability.")

    add("Mutant", FAM_FUSION, EV_B, CAT_FUSION, "Unknown_Ferritin_Scaffold", "Ferritin-lichenase fusion", "", "",
        "Ferritin-lichenase N/C fusion", "Ferritin-lichenase fusion", "Terminal_fusion", "t_half_50C", "50.0",
        "t1/2 at 50 °C improved 11.09-fold vs free lichenase", "Threshold_Survival", "", "11.09x t1/2 at 50 C",
        "fold_change", "", "Range_Or_Qualitative", "Enzyme_Thermal_HalfLife", "As-isolated", "50 C thermal inactivation",
        "~12 nm soluble polymeric nanoparticle", "PubMed_Abstract_ClosedPaper", "35470630", "", "10.13345/j.cjb.210913",
        "2022 Sheng Wu Gong Cheng Xue Bao (PMID:35470630)", "", "Chinese journal paper; CNKI full text needed for exact WT vs fusion t1/2 minutes.")

    add("Mutant", FAM_DPS, EV_B, CAT_MULTI, "Listeria innocua", "Listeria innocua Dps (LiDps)", "P80725", "",
        "LiDps-tm (triple ferroxidase mutant)", "LiDps triple ferroxidase mutant (LiDps-tm)", "Ferroxidase_center", "T_survival",
        "65.0", "Intact mineralization nanocage at 65 °C", "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative",
        "Heat_Treatment_TEM", "Holo_MNP", "65 C, alkaline pH, H2O2", "12-mer",
        "PubMed_Abstract_ClosedPaper", "19859920", "", "10.1002/chem.200901138",
        "2009 Chem Eur J (PMID:19859920)", "", "Survives 65 C nanoparticle mineralization.")

    # =========================================================================
    # 20. Tier C Wet-Lab Qualitative / Comparative Records (Ferritin, Bfr, Dps, Fusion)
    # =========================================================================
    qual_entries = [
        # (fam, cat, org, prot, up, pdb, mraw, mnorm, iface, disp, meth, state, buf, oligo, stype, pmid, pmcid, doi, cite, note)
        (FAM_FERRITIN, CAT_INDEL, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "", "HuHF C-terminal 2 variants (channel/loop screen of 22 mutants)", "HuHF 2 C-term mutants", "C_term", "< WT thermostability (diminished thermostability)", "Thermal_Abstract", "As-isolated", "Biochem J 1989 (PMC1133592)", "Assembled", "PMC_OpenAccess", "2690826", "PMC1133592", "", "Luzzago & Cesareni 1989 Biochem J", "PMC1133592 screen of 22 channel/loop mutants; ΔC22 & ΔC7 assembled, ΔC28 & 2 single loop deletions failed to assemble."),
        (FAM_FERRITIN, CAT_POINT, "Homo sapiens", "Human H-ferritin (HuHF)", "P02794", "", "Variant #1 (Trp34/channel mutant)", "HuHF variant #1", "Channel/Ferroxidase", "CD & thermal melting curve ~ WT (Figs S1-S2)", "CD_222nm", "Apo", "PMC2653062", "24-mer", "PMC_OpenAccess", "", "PMC2653062", "", "PMC2653062 (2009)", "SI Figs S1-S2 melting curves essentially unchanged vs WT."),
        (FAM_FUSION, CAT_FUSION, "Homo sapiens", "Human H-ferritin (HFt)", "P02794", "", "HFt-MSH (alpha-MSH N-term fusion)", "HFt-MSH", "N_term_fusion", "> WT HFt thermal stability (CD 222 nm 30-100 C, Fig S7)", "CD_222nm", "As-isolated", "PMC3356193 Fig S7 (30-100 C)", "24-mer", "PMC_OpenAccess", "22619508", "PMC3356193", "10.1038/cgt.2012.22", "Falvo et al. 2012 Cancer Gene Ther", "Quantitative curve in SI Fig S7 only."),
        (FAM_FUSION, CAT_FUSION, "Homo sapiens", "Human H-ferritin (FTH1-1)", "P02794", "", "FTH1-1 (Q15R/N110E/S114E) + GDIII-Fagy tag", "Q15R/N110E/S114E + GDIII-Fagy", "Hydrophobic_core_I + tag", "CFX96 DSF TSA measured (10-95 C, 0.5 C/s, Fig S4g)", "DSF_SYPRO", "As-isolated", "10-95 C, 0.5 C/s (Fig S4g)", "24-mer", "PMC_OpenAccess", "39366932", "PMC11452399", "10.1038/s41467-024-52908-z", "2024 Nat Commun (PMC11452399)", "Numeric Tm not printed in main text; figure-only in Fig S4g."),
        (FAM_BFR, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "", "N23L", "N23L", "C2_water_pocket", "> WT thermostability & chemical stability", "CD_and_Chemical", "As-isolated", "Biochemistry 2011 (PMC3179174 / PMID:21488690)", "24-mer + Dimer", "PubMed_Abstract_ClosedPaper", "21488690", "PMC3179174", "10.1021/bi200207w", "Ardejani et al. 2011 Biochemistry", "Plugs C2 interfacial water pocket; exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "", "N23F", "N23F", "C2_water_pocket", "> WT thermostability; shifts population toward 24-mer", "CD_and_Chemical", "As-isolated", "Biochemistry 2011", "24-mer enriched", "PubMed_Abstract_ClosedPaper", "21488690", "PMC3179174", "10.1021/bi200207w", "Ardejani et al. 2011 Biochemistry", "Strongest 24-mer oligomer shift among N23X series; exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_POINT, "Escherichia coli", "Bacterioferritin (EcBfr)", "P0ABD3", "", "N23W", "N23W", "C2_water_pocket", "> WT thermostability & chemical stability", "CD_and_Chemical", "As-isolated", "Biochemistry 2011", "24-mer + Dimer", "PubMed_Abstract_ClosedPaper", "21488690", "PMC3179174", "10.1021/bi200207w", "Ardejani et al. 2011 Biochemistry", "Exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_POINT, "Mycobacterium tuberculosis", "Bacterioferritin B (MtBfrB)", "P9WNE5", "", "Q51A", "Q51A", "Intersubunit_interface", "> WT thermal stability (without affecting assembly)", "Thermal_Abstract", "As-isolated", "Biochemistry 2013 (closed DOI)", "24-mer", "PubMed_Abstract_ClosedPaper", "23409758", "", "10.1021/bi301189t", "Khare et al. 2013 Biochemistry", "Single-point thermostabilizing mutant; exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_POINT, "Mycobacterium tuberculosis", "Bacterioferritin B (MtBfrB)", "P9WNE5", "", "Q126A", "Q126A", "Intersubunit_interface", "> WT thermal stability (without affecting assembly)", "Thermal_Abstract", "As-isolated", "Biochemistry 2013 (closed DOI)", "24-mer", "PubMed_Abstract_ClosedPaper", "23409758", "", "10.1021/bi301189t", "Khare et al. 2013 Biochemistry", "Exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_POINT, "Mycobacterium tuberculosis", "Bacterioferritin B (MtBfrB)", "P9WNE5", "", "E135A", "E135A", "Intersubunit_interface", "> WT thermal stability (without affecting assembly)", "Thermal_Abstract", "As-isolated", "Biochemistry 2013 (closed DOI)", "24-mer", "PubMed_Abstract_ClosedPaper", "23409758", "", "10.1021/bi301189t", "Khare et al. 2013 Biochemistry", "Exact Tm behind ACS paywall."),
        (FAM_BFR, CAT_INDEL, "Mycobacterium tuberculosis", "Bacterioferritin B (MtBfrB)", "P9WNE5", "", "del 168-181 (C-term C4 extension deleted)", "Δ168-181", "C_term_extension", "< WT thermostability (decreased thermostability)", "UniProt_Curated", "As-isolated", "UniProt P9WNE5 annotation", "24-mer (no large assembly change)", "StructuredDB_UniProt", "21494619", "", "", "Khare et al. 2011 (UniProt P9WNE5)", "3.5-fold lower Fe2+ oxidation; 20% lower Fe3+ release."),
        (FAM_FERRITIN, CAT_MULTI, "Rana catesbeiana", "Bullfrog M-ferritin", "P07798", "", "D80K/I81K (D80KI81K)", "D80K/I81K", "L_loop_C2", "< WT thermal stability & inhibited assembly", "Thermal_Abstract", "As-isolated", "ACS Chem Biol 2014 (closed DOI)", "Impaired 24-mer", "PubMed_Abstract_ClosedPaper", "25148224", "", "10.1021/cb500431r", "2014 ACS Chem Biol (PMID:25148224)", "Introduces contiguous Lys at L-loop pair (D80K single mutant also slows nanocage assembly)."),
        (FAM_FERRITIN, CAT_POINT, "Homo sapiens", "Human L-ferritin (HuFTL)", "P02792", "", "A96T (neuroferritinopathy variant)", "A96T", "Subunit_interior", "< HuFTL WT stability against heat, acid, denaturant", "Thermal_Abstract", "Apo/As-isolated", "Biochemistry 2019 (closed DOI)", "24-mer (lower refolding, aggregates)", "PubMed_Abstract_ClosedPaper", "30986045", "", "10.1021/acs.biochem.8b01068", "2019 Biochemistry (PMID:30986045)", "Closed ACS paper; exact Tm behind paywall."),
        (FAM_FUSION, CAT_FUSION, "Homo sapiens", "Glucagon-ferritin N-term fusion (GFTN)", "P02792/P02794", "", "GFTN (29-aa human glucagon N-term fusion)", "Glucagon(1-29)-FTN (GFTN)", "N_term_fusion", "> FTN WT thermal stability (>50% soluble after 40 min vs >80% FTN ppt in 5 min; k_agg 2.8e-4 vs 6.04e-3 s-1)", "Heat_Clarification_Solubility", "As-isolated", "High conc 1.3 g/L", "24-mer", "PubMed_Abstract_ClosedPaper", "11708788", "", "10.1006/bbrc.2001.5931", "2001 BBRC (PMID:11708788)", "Closed Elsevier paper."),
        (FAM_FUSION, CAT_INDEL, "Escherichia coli", "BFR / DPS helix-swap & deletion chimeras", "P0ABD3/P0ABT2", "", "BFR/DPS BC-helix & E-helix deletions/swaps (incl. DPS+E-helix)", "BFR/DPS helix deletions & DPS+E-helix swap", "E_helix / BC_helix", "Cooperative CD thermal unfolding observed across all helical mutants", "CD_222nm", "As-isolated", "Biochemistry 2009 (closed DOI)", "DPS+E-helix novel 12-mer; most BFR mutants fail 24-mer", "PubMed_Abstract_ClosedPaper", "19405543", "", "10.1021/bi900387t", "Fan et al. 2009 Biochemistry (PMID:19405543)", "Closed ACS paper; exact CD Tm values behind paywall."),
        (FAM_FERRITIN, CAT_MULTI, "Homo sapiens", "Human H-ferritin Cys-free + His mutant", "P02794", "9LUW", "C90A/C102A/H105A/E116A/T122H/C130A/L165H", "C90A/C102A/H105A/E116A/T122H/C130A/L165H", "C3_and_C4_interfaces", "Significantly improved thermal stability & monodispersity vs WT rHuHF", "Thermal_Abstract", "As-isolated", "J Agric Food Chem 2025 (closed DOI)", "24-mer", "PubMed_Abstract_ClosedPaper", "40377399", "", "10.1021/acs.jafc.5c02232", "2025 JAFC (PMID:40377399)", "AlphaFold3 design validated by PDB:9LUW; exact Tm behind ACS paywall."),
        (FAM_FUSION, CAT_MULTI, "Homo sapiens", "F1L3E1 (HuHF-L3-EBNA1) E-helix variants C1, C2, C5", "P02794", "", "F1L3E1 C1 (E167Q/D171Q), C2 (E167R/D171R), C5 (heptad a/d->Q)", "F1L3E1 E-helix variants C1, C2, C5", "C_term_E_helix", "< F1L3E1 parent thermostability (C3/C4 hydrophobic swaps cause severe aggregation)", "MD_and_Experiment", "As-isolated", "ACS Appl Bio Mater 2022 (closed DOI)", "24-mer (C1/C2/C5 folded); C3/C4 aggregated", "PubMed_Abstract_ClosedPaper", "35770389", "", "10.1021/acsabm.2c00154", "2022 ACS Appl Bio Mater (PMID:35770389)", "Exact values behind ACS paywall."),
        (FAM_DPS, CAT_POINT, "Marinobacter hydrocarbonoclasticus", "M. hydrocarbonoclasticus Dps (MhDps)", "", "", "Dps Q14E", "Q14E", "N_terminal_tail", "~ MhDps WT thermal stability by SRCD (aggregates upon thermal denaturation)", "SRCD", "Apo & Holo", "Eur Biophys J 2021 (OA decoded PDF)", "12-mer", "OA_PDF", "33900431", "", "10.1007/s00249-021-01538-0", "2021 Eur Biophys J (PMID:33900431)", "Eliminates DNA binding without altering secondary structure or thermostability vs WT."),
        (FAM_DPS, CAT_POINT, "Escherichia coli / Mycobacterium smegmatis", "Dps trimer-interface mutants (E146A, F47E, E146A/F47E)", "P0ABT2/P0C558", "", "E146A, F47E, E146A/F47E", "Dps E146A, F47E, E146A/F47E", "C3_ferritin_like / Dps_like", "E146A trimer/dodecamer T-dependent equilibrium; F47E trimer; E146A/F47E monomer", "Review_PMC3179174", "As-isolated", "Reviewed in PMC3179174", "E146A temp-dependent 3/12-mer", "PMC_OpenAccess", "21954367", "PMC3179174", "10.1016/j.bbagen.2011.09.005", "Fan et al. 2012 BBA (review)", "See Zhang et al. / Ceci et al. cited inside PMC3179174 [74,75]."),
    ]
    for row in qual_entries:
        fam, cat, org, prot, up, pdb, mraw, mnorm, iface, disp, meth, state, buf, oligo, stype, pmid, pmcid, doi, cite, note = row
        add("Mutant", fam, EV_C, cat, org, prot, up, pdb, mraw, mnorm, iface, "Qualitative_Stability",
            "", disp, "Qualitative", "", "", "NA", "", "Range_Or_Qualitative",
            meth, state, buf, oligo, stype, pmid, pmcid, doi, cite, "", note)

    # =========================================================================
    # 21. Family 3 Dps Quantitative (PMC13319912 DgrDps) & Family 5 Extended Superfamily
    # =========================================================================
    # PMC13319912 DgrDps (SRCD 192 nm)
    dgr_buf = "SRCD 192 nm, 24-87 C"
    add("WT_baseline", FAM_DPS, EV_A, CAT_WT, "Deinococcus grandis", "D. grandis Dps (DgrDps)", "", "",
        "DgrDps WT", "WT", "WT_reference", "Tm", "67.0", "67 ± 1 °C (ΔHm = 0.9 ± 0.1e2 kcal/mol)", "Exact_Numeric",
        "67.0", "0", "NA", "0.0", "NA", "SRCD_192nm", "Apo", dgr_buf, "12-mer",
        "PMC_OpenAccess", "41790146", "PMC13319912", "10.1002/pro.70492", "2026 Protein Sci (PMC13319912)", "",
        "WT + pUC19 DNA shows two transitions: 54 ± 1 C and 69 ± 1 C.")
    add("Mutant", FAM_DPS, EV_A, CAT_INDEL, "Deinococcus grandis", "D. grandis Dps (DgrDps)", "", "",
        "DgrDps ΔN (Δ1-46)", "Δ1-46 (DgrDps ΔN)", "N_terminal_tail", "Tm", "70.0", "70 ± 1 °C (ΔHm = 1.4 ± 0.2e2 kcal/mol)",
        "Exact_Numeric", "67.0", "+3.0", "mut_minus_WT", "+3.0", "Derived_mut_minus_WT",
        "SRCD_192nm", "Apo", dgr_buf, "12-mer", "PMC_OpenAccess", "41790146", "PMC13319912", "10.1002/pro.70492",
        "2026 Protein Sci (PMC13319912)", "", "Removing flexible N-terminal DNA-binding tail increases Tm by +3 C.")

    # Family 5: Artemia franciscana Artemin (PMC3059798 + PMID:21600915)
    for mnorm, rank_note, t_val in [
        ("C172A", "75 °C 13 min soluble survival rank: WT > C172A > C166A > C22A > C61A", "75.0"),
        ("C166A", "75 °C 13 min soluble survival rank: WT > C172A > C166A > C22A > C61A", "75.0"),
        ("C22A", "75 °C 13 min soluble survival rank: WT > C172A > C166A > C22A > C61A", "75.0"),
        ("C61A", "75 °C 13 min soluble survival rank: WT > C172A > C166A > C22A > C61A (most labile)", "75.0"),
    ]:
        add("Mutant", FAM_EXT, EV_C, CAT_POINT, "Artemia franciscana", "Artemin (ferritin-superfamily chaperone nanocage)", "", "",
            f"Art{mnorm}", mnorm, "Cys_disulfide_stability", "Qualitative_Stability", "", rank_note, "Qualitative",
            "", "< WT", "mut_minus_WT", "", "Range_Or_Qualitative", "Heat_Clarification_Solubility",
            "As-isolated", "75 C 13 min heat treatment", "24-mer-like chaperone",
            "PMC_OpenAccess", "20878295", "PMC3059798", "10.1007/s12192-010-0222-x", "Rasti et al. 2011 Cell Stress Chaperones", "",
            "Relative soluble recovery after 75 C 13 min stress.")
    add("Mutant", FAM_EXT, EV_C, CAT_INDEL, "Artemia franciscana", "Artemin (ferritin-superfamily chaperone nanocage)", "", "",
        "ΔC39 (C-term 39 aa deleted)", "ΔC39", "C_terminal_tail", "Qualitative_Stability", "",
        "< WT structural stability & reduced heat resistance", "Qualitative", "", "< WT", "mut_minus_WT", "", "Range_Or_Qualitative",
        "In_vivo_Thermotolerance_And_Chaperone", "As-isolated", "IJBM 2011 abstract (closed DOI)", "Oligomer",
        "PubMed_Abstract_ClosedPaper", "21600915", "", "10.1016/j.ijbiomac.2011.05.005", "2011 Int J Biol Macromol (PMID:21600915)", "",
        "Closed Elsevier paper.")

    # Family 5: S. Typhimurium YciF (DUF892, PMC10367587)
    add("Mutant", FAM_EXT, EV_A, CAT_POINT, "Salmonella Typhimurium", "YciF (DUF892 ferritin-like protein)", "", "",
        "E143D", "E143D", "Dinuclear_metal_site", "Tm", "53.5", "53.5 °C (+Mg2+ +4 C; +Zn2+ ~-5 C; +Fe >95 C)", "Exact_Numeric",
        "", "", "NA", "", "NA", "DSF_SYPRO", "Apo", "SYPRO Orange TSA", "Dimer/oligomer",
        "PMC_OpenAccess", "37367303", "PMC10367587", "10.1016/j.jbc.2023.104930", "2023 JBC (PMC10367587)", "",
        "YciF belongs to extended ferritin-like DUF892 superfamily; E113Q & Q54A show ~WT Tm.")
    for mut in ("E113Q", "Q54A"):
        add("Mutant", FAM_EXT, EV_C, CAT_POINT, "Salmonella Typhimurium", "YciF (DUF892 ferritin-like protein)", "", "",
            mut, mut, "Dinuclear_metal_site", "Qualitative_Stability", "", "~ YciF WT Tm (TSA)", "Qualitative",
            "", "~0", "mut_minus_WT", "0.0", "Range_Or_Qualitative", "DSF_SYPRO", "Apo", "SYPRO Orange TSA", "Dimer/oligomer",
            "PMC_OpenAccess", "37367303", "PMC10367587", "10.1016/j.jbc.2023.104930", "2023 JBC (PMC10367587)", "",
            "Exact decimal WT Tm not reprinted in main text snippet; similar to WT.")

    # Family 5: MrkA-dodecin (PMC13107372)
    add("WT_baseline", FAM_EXT, EV_A, CAT_WT, "Klebsiella pneumoniae / Thermophilic dodecin", "MrkA monomer (baseline)", "", "",
        "MrkA monomer", "WT_MrkA_monomer", "WT_reference", "Ti", "87.5", "87.5 °C", "Exact_Numeric", "87.5", "0", "NA", "0.0", "NA",
        "Tycho_NT6", "As-isolated", "Tycho NT.6 thermal shift", "Monomer", "PMC_OpenAccess", "41956462", "PMC13107372",
        "10.1038/s41541-026-01417-7", "2026 npj Vaccines (PMC13107372)", "", "Baseline unfused MrkA antigen.")
    add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Klebsiella pneumoniae / Thermophilic dodecin", "MrkA-dodecin nanocage fusion", "", "",
        "MrkA-dodecin", "MrkA-dodecin fusion", "Dodecin_scaffold_fusion", "Ti", "90.1", "90.1 °C", "Exact_Numeric",
        "87.5", "+2.6", "mut_minus_WT", "+2.6", "Derived_mut_minus_WT", "Tycho_NT6", "As-isolated", "Tycho NT.6 thermal shift",
        "12-mer dodecin nanocage", "PMC_OpenAccess", "41956462", "PMC13107372", "10.1038/s41541-026-01417-7",
        "2026 npj Vaccines (PMC13107372)", "", "Extended nanocage scaffold (dodecin, not true 4-helix ferritin).")

    # Family 5: Quasibacillus thermotolerans encapsulin QtEnc (PMC12124451)
    add("WT_baseline", FAM_EXT, EV_A, CAT_WT, "Quasibacillus thermotolerans", "Encapsulin (QtEnc, E. coli expressed)", "", "",
        "WT QtEnc (E. coli)", "WT", "WT_reference", "Tm", "84.4", "84.4 °C", "Exact_Numeric", "84.4", "0", "NA", "0.0", "NA",
        "nanoDSF", "As-isolated", "nanoDSF in PMC12124451", "T=4 240-mer encapsulin shell",
        "PMC_OpenAccess", "40139971", "PMC12124451", "10.1038/s42004-025-01493-3", "2025 Commun Chem (PMC12124451)", "",
        "Baseline E. coli expressed QtEnc shell.")
    add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Quasibacillus thermotolerans", "Encapsulin (QtEnc) fusion", "", "",
        "LanM-QtEnc (uncleaved N-term fusion)", "LanM-QtEnc", "N_term_fusion", "Tm", "55.9", "55.9 °C", "Exact_Numeric",
        "84.4", "-28.5", "mut_minus_WT", "-28.5", "Derived_mut_minus_WT", "nanoDSF", "As-isolated", "nanoDSF", "Unassembled/destabilized fusion",
        "PMC_OpenAccess", "40139971", "PMC12124451", "10.1038/s42004-025-01493-3", "2025 Commun Chem (PMC12124451)", "",
        "N-terminal LanM fusion severely destabilizes QtEnc (-28.5 C).")
    add("Mutant", FAM_EXT, EV_A, CAT_INDEL, "Quasibacillus thermotolerans", "Encapsulin (QtEnc, TEV-cleaved reassembled)", "", "",
        "TEV-cleaved in vitro assembled QtEnc", "QtEnc (TEV-cleaved in vitro assembled)", "N_term_scar_reassembled", "Tm",
        "100.0", "100 °C", "Exact_Numeric", "84.4", "+15.6", "mut_minus_WT", "+15.6", "Derived_mut_minus_WT",
        "nanoDSF", "As-isolated", "nanoDSF", "T=4 240-mer reassembled",
        "PMC_OpenAccess", "40139971", "PMC12124451", "10.1038/s42004-025-01493-3", "2025 Commun Chem (PMC12124451)", "",
        "TEV-cleaved reassembled QtEnc reaches 100 C (+15.6 C vs WT).")

    # Family 5: T. maritima encapsulin cargo/chaperone variants (PMC11190991)
    add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Thermotoga maritima", "Transketolase + Tm_encap cargo peptide (TK-TmCLP / WT TK)", "", "",
        "TK / TK-TmCLP (ITF)", "TK-TmCLP", "Cargo_loading_peptide", "Tm", "60.85", "334 K (60.85 °C)", "Exact_Numeric",
        "", "", "NA", "", "NA", "ITF", "As-isolated", "Intrinsic Trp fluorescence (K -> C converted by -273.15)", "Cargo dimer",
        "PMC_OpenAccess", "38835217", "PMC11190991", "10.1038/s41598-024-64256-5", "2024 Sci Rep (PMC11190991)", "",
        "Converted from 334 K (-273.15).")
    add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Thermotoga maritima", "TK-TmCLP K603pAzF", "", "",
        "K603pAzF", "TK-TmCLP K603pAzF", "Noncanonical_aa_site", "Tm", "60.65", "333.8 K (60.65 °C)", "Exact_Numeric",
        "60.85", "-0.2", "mut_minus_WT", "-0.2", "Derived_mut_minus_WT", "ITF", "As-isolated", "ITF (333.8 K - 273.15)", "Cargo",
        "PMC_OpenAccess", "38835217", "PMC11190991", "10.1038/s41598-024-64256-5", "2024 Sci Rep (PMC11190991)", "",
        "Photocrosslinkable noncanonical amino acid variant.")
    add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Thermotoga maritima", "Tm_encap-encapsulated TK-TmCLP (FRET)", "", "",
        "Encapsulated TK-TmCLP vs free (FRET)", "Tm_encap-encapsulated TK-TmCLP", "Encapsulin_cage_confinement", "Tm",
        "62.05", "335.2 ± 2.7 K (62.05 °C) vs free 320.5 ± 2.3 K (47.35 °C)", "Exact_Numeric",
        "47.35", "+14.7 (+15 C)", "mut_minus_WT", "+14.7", "Derived_mut_minus_WT", "FRET", "Encapsulated",
        "FRET thermal transition (K -> C)", "Tm_encap 60-mer + cargo",
        "PMC_OpenAccess", "38835217", "PMC11190991", "10.1038/s41598-024-64256-5", "2024 Sci Rep (PMC11190991)", "",
        "Paper summarizes encapsulation improves Tm by ~15 C (335.2 - 320.5 = 14.7 K).")

    # Family 5: A. baumannii (PMC11284452) & S. griseus (PMC11550324) ligand/cargo complexes
    for lig, dtm in [("thiosulfate", "+10.5"), ("GSH", "+3.3"), ("sulfate", "+2.2")]:
        add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Acinetobacter baumannii", "Encapsulin cargo Ab-RLD", "", "",
            f"Ab-RLD + {lig}", f"Ab-RLD + {lig} (ligand stabilization)", "Cargo_active_site", "delta_Tm",
            "", f"{dtm} °C vs apo Ab-RLD", "Exact_Numeric", "", dtm, "mut_minus_WT", dtm, "Directly_Reported",
            "DSF_SYPRO", f"+{lig}", "PMC11284452 DSF", "Cargo protein",
            "PMC_OpenAccess", "39073218", "PMC11284452", "10.1016/j.jbc.2024.107625", "2024 JBC (PMC11284452)", "",
            "Ligand-induced thermal shift of encapsulin-associated rhodanese cargo (not amino acid mutation of shell).")
    for mname, val, lig in [("Tm1", "55.0 ± 0.49", "apo"), ("Tm2", "70.1 ± 0.17", "apo"), ("Tm1", "55.7 ± 2.27", "+cAMP"), ("Tm2", "70.5 ± 0.33", "+cAMP")]:
        num_v = val.split(" ")[0]
        add("Mutant", FAM_EXT, EV_A, CAT_FUSION, "Streptomyces griseus", "Sg Enc_2-MIBS encapsulin-terpene cyclase nanocomplex", "", "",
            f"Sg Enc_2-MIBS ({lig})", "Sg Enc_2-MIBS nanocomplex", "Encapsulin_cargo_complex", mname,
            num_v, f"{val} °C", "Exact_Numeric", "", "", "NA", "", "NA", "DSF_SYPRO", lig, "DSF in PMC11550324",
            "Encapsulin + MIBS cargo", "PMC_OpenAccess", "39521781", "PMC11550324", "10.1074/jbc.2024", "2024 PMC11550324", "",
            "Two-transition DSF melting of native encapsulin-cargo complex.")

    # Family 5: TmEnc S1X-S7X explicitly named purified variants (PMC8816334)
    enc_buf = "150 mM NaCl, 25 mM Tris pH 8.0 (DLS 25->75 C, 2 C/step + overnight 80 C / 10 h TEM)"
    for mut, pdb, cat, tier, t_val, disp, note in [
        ("S1K", "7LIJ", CAT_POINT, EV_B, "80.0", "Thermostable: no irreversible PDI rise <=75 C; intact cages after 80 C 10 h", "Cryo-EM PDB:7LIJ; explicitly verified in PMC8816334 text & PDB."),
        ("S1R", "7LIK", CAT_POINT, EV_B, "80.0", "Thermostable: no irreversible PDI rise <=75 C; intact cages after 80 C 10 h", "Cryo-EM PDB:7LIK; explicitly verified in PMC8816334 text & PDB."),
        ("S3E", "7LIL", CAT_MULTI, EV_B, "80.0", "Thermostable: no irreversible PDI rise <=75 C; intact cages after 80 C 10 h", "Cryo-EM PDB:7LIL ('SE3'); explicitly verified in PMC8816334 text."),
        ("S5D", "7LIS", CAT_INDEL, EV_B, "80.0", "Thermostable: partial loop deletion; no irreversible PDI rise <=75 C; intact 80 C", "Cryo-EM PDB:7LIS; explicitly verified in PMC8816334 text."),
        ("S6E", "7LIM", CAT_INDEL, EV_B, "80.0", "Thermostable: loop deletion; no irreversible PDI rise <=75 C; intact 80 C", "Cryo-EM PDB:7LIM; explicitly verified in PMC8816334 text."),
        ("S7D", "7LII", CAT_INDEL, EV_B, "80.0", "Thermostable: loop deletion; no irreversible PDI rise <=75 C; intact 80 C", "Cryo-EM PDB:7LII; explicitly verified in PMC8816334 text."),
        ("S2E", "", CAT_MULTI, EV_B, "75.0", "Thermally destabilized <=75 C (irreversible DLS PDI increase; partial TEM cages survive 80 C)", "Explicitly named in Fig 2B text as showing irreversible PDI increase <=75 C."),
        ("S4D", "", CAT_INDEL, EV_B, "75.0", "Thermally destabilized <=75 C (lesser irreversible DLS PDI increase; partial TEM cages survive 80 C)", "Explicitly named in Fig 2B text."),
        ("S7G", "7LIT", CAT_INDEL, EV_B, "75.0", "Positive control; irreversible DLS PDI increase <=75 C", "Explicitly named in Fig 2B text & PDB:7LIT."),
        ("S3D/S3K/S3R/S5K (insoluble/host-toxic subset)", "", CAT_MULTI, EV_C, "", "Insoluble (S3D) or host-toxic (S3K, S3R, S5K); failed purified assembly", "Note: remaining 24-library designs reside in Fig 1C/2B graphic & SI Table S1-S3 (logged as ESI follow-up)."),
    ]:
        add("Mutant", FAM_EXT, tier, cat, "Thermotoga maritima", "Encapsulin shell (TmEnc)", "Q9WZP2", pdb,
            f"TmEnc {mut}", f"TmEnc {mut}", "C5_fivefold_pore", "T_survival" if t_val else "Qualitative_Stability",
            t_val, disp, "Threshold_Survival" if t_val else "Qualitative", "", "", "NA", "", "Range_Or_Qualitative",
            "DLS_and_TEM", "As-isolated", enc_buf, "60-mer T=1 encapsulin",
            "PMC_OpenAccess", "35119930", "PMC8816334", "10.1126/sciadv.abl7346", "Williams / Adamson et al. 2022 Sci Adv", "", note)

    # Family 5: Myxococcus xanthus encapsulin MxENC (PMC11866011)
    mx_buf = "DLS 1 h incubation at 20, 30, 40, 50 C (WT disassembles irreversibly at 60 C)"
    add("Mutant", FAM_EXT, EV_B, CAT_FUSION, "Myxococcus xanthus", "Encapsulin (MxEnc)", "Q1D6H4", "",
        "EncA-ST (SpyTag insertion)", "EncA-ST", "Surface_loop_SpyTag", "T_survival", "50.0",
        "Stable up to 50 °C (1 h DLS)", "Threshold_Survival", "", "", "NA", "", "Range_Or_Qualitative",
        "DLS", "As-isolated", mx_buf, "T=3 180-mer", "PMC_OpenAccess", "40028083", "PMC11866011",
        "10.1002/mabi.202400559", "2025 Macromol Biosci (PMC11866011)", "", "Surface SpyTag fusion stable up to 50 C.")
    add("Mutant", FAM_EXT, EV_B, CAT_FUSION, "Myxococcus xanthus", "Encapsulin (MxEnc)", "Q1D6H4", "",
        "EncA-157ST (internal pos 157 SpyTag)", "EncA-157ST", "Loop_pos157_SpyTag", "T_survival", "40.0",
        "Minimal aggregation at 40 °C; strong diameter increase at 50 °C (< WT 60 °C)", "Threshold_Survival",
        "60.0", "<60 C", "mut_minus_WT", "", "Range_Or_Qualitative", "DLS", "As-isolated", mx_buf, "T=3 180-mer",
        "PMC_OpenAccess", "40028083", "PMC11866011", "10.1002/mabi.202400559", "2025 Macromol Biosci (PMC11866011)", "",
        "Internal pos 157 SpyTag destabilizes cage relative to EncA-ST and WT.")

    return records


def build_unique_mutant_rollup(records: list[dict]) -> tuple[list[dict], list[dict]]:
    """Assign deterministic UM IDs and roll up condition-level records into deduplicated unique mutants."""
    rank_tier = {EV_A: 1, EV_B: 2, EV_C: 3, EV_D: 4}
    grouped = defaultdict(list)
    for r in records:
        if r["record_role"] == "Mutant":
            grouped[r["unique_mutant_key"]].append(r)

    unique_rows = []
    # Sort keys deterministically by family_tier, highest_evidence_tier, organism, protein, mutation
    sorted_keys = sorted(
        grouped.keys(),
        key=lambda k: (
            grouped[k][0]["family_tier"],
            min(rank_tier[x["evidence_tier"]] for x in grouped[k]),
            grouped[k][0]["organism"],
            grouped[k][0]["protein_name"],
            grouped[k][0]["mutation_normalized"],
        ),
    )

    key_to_umid = {}
    for i, k in enumerate(sorted_keys, 1):
        umid = f"UM{i:04d}"
        key_to_umid[k] = umid
        grp = grouped[k]
        for x in grp:
            x["unique_mutant_id"] = umid

        fam = grp[0]["family_tier"]
        best_ev = min(grp, key=lambda x: rank_tier[x["evidence_tier"]])["evidence_tier"]
        cat = grp[0]["mutation_category"]
        org = grp[0]["organism"]
        prot = grp[0]["protein_name"]
        ups = ";".join(sorted({x["uniprot_id"] for x in grp if x["uniprot_id"]}))
        pdbs = ";".join(sorted({x["pdb_id"] for x in grp if x["pdb_id"]}))
        mut_norm = grp[0]["mutation_normalized"]

        wet_recs = [x for x in grp if x["evidence_tier"] != EV_D]
        quant_tm_recs = [x for x in grp if x["evidence_tier"] == EV_A and x["value_numeric_C"]]
        comp_recs = [x for x in grp if x["evidence_tier"] == EV_D]

        best_tm = "; ".join(f"{x['metric_name']}={x['value_numeric_C']}°C({x['assay_method']})" for x in quant_tm_recs[:4])
        best_dtm = "; ".join(f"{x['normalized_delta_Tm_C']}°C" for x in grp if x["normalized_delta_Tm_C"])
        comp_vals = "; ".join(f"{x['metric_name']}={x['reported_delta_value']}" for x in comp_recs[:3])

        unique_rows.append({
            "unique_mutant_id": umid,
            "unique_mutant_key": k,
            "family_tier": fam,
            "highest_evidence_tier": best_ev,
            "mutation_category": cat,
            "organism": org,
            "protein_name": prot,
            "uniprot_id": ups,
            "pdb_ids": pdbs,
            "mutation_normalized": mut_norm,
            "n_condition_records": len(grp),
            "n_wetlab_records": len(wet_recs),
            "n_quant_temp_records": len(quant_tm_recs),
            "has_quantitative_temp": "Yes" if quant_tm_recs or any(x["evidence_tier"] == EV_A for x in grp) else "No",
            "quantitative_temp_values_C": best_tm,
            "normalized_delta_Tm_values_C": best_dtm,
            "computational_predictions": comp_vals,
            "metrics_available": ";".join(sorted({x["metric_name"] for x in grp})),
            "assay_methods": ";".join(sorted({x["assay_method"] for x in grp})),
            "cross_ref_groups": ";".join(sorted({x["cross_ref_duplicate_group"] for x in grp if x["cross_ref_duplicate_group"]})),
            "record_ids": ";".join(x["record_id"] for x in grp),
        })

    return unique_rows, records


def build_inaccessible_sources() -> list[dict]:
    """Compile comprehensive follow-up ledger from closed DOIs, ESI gaps, TBP PDBs, ZH/DB gaps."""
    gaps = []

    def add_gap(priority, gap_type, identifier, year, title_or_name, reason, expected_data, suggested_action):
        gaps.append({
            "gap_id": f"GAP{len(gaps)+1:04d}",
            "priority": priority,
            "gap_type": gap_type,
            "identifier_doi_pmid_pdb": identifier,
            "year": year,
            "title_or_resource": title_or_name,
            "inaccessibility_reason": reason,
            "expected_mutant_thermal_content": expected_data,
            "recommended_follow_up": suggested_action,
        })

    # Priority 1 & 2: OA Supplementary / Figure-only leads
    oa_esi_gaps = [
        ("P1_High_Quant_Curve_In_ESI", "OA_ESI_Supplementary", "PMCID:PMC11776510 | DOI:10.1039/d4sc07021f", "2025",
         "M. tuberculosis BfrA E121A/E121Q/E121K/E121F CD thermal melt curves (Fig S13)",
         "Exact numerical Tm decimals reside inside RSC external ESI PDF (only approximate ~90/84/78/80 C curve values captured).",
         "4 MtBfrA E121X variants exact CD Tm decimals", "Download Chem Sci 2025 ESI PDF (Fig S13) and verify exact Tm values."),
        ("P1_High_Quant_FigureOnly", "OA_FigureOnly_Plot", "PMCID:PMC9798250 | DOI:10.1002/pro.4543", "2023",
         "Human H/L heteropolymer intermediate stoichiometry NanoDSC Tm curve (Fig 3b)",
         "Homopolymers (92-93 & 115 C) and extremes (H2:L22, H6:L18, H19:L5) given in text; intermediate H:L ratio Tm points plotted in Fig 3b.",
         "Full H:L stoichiometric NanoDSC curve", "Read Fig 3b data points or author SI table from Wiley."),
        ("P2_Medium_FigureOnly", "OA_FigureOnly_Plot", "PMCID:PMC11452399 | DOI:10.1038/s41467-024-52908-z", "2024",
         "FTH1-1 (Q15R/N110E/S114E) + GDIII-Fagy tag CFX96 DSF TSA (Fig S4g)",
         "Main text & decoded PDF describe assay (10-95 C, 0.5 C/s) without printing numeric Tm in text; plotted in Fig S4g.",
         "Exact DSF Tm for FTH1-1 variant vs WT", "Check Nature Communications Fig S4g source data file."),
        ("P2_Medium_FigureOnly", "OA_FigureOnly_Plot", "PMCID:PMC3356193 | DOI:10.1038/cgt.2012.22", "2012",
         "HFt-MSH vs WT HFt CD 222 nm thermal curve (Fig S7)",
         "Qualitative improvement stated in text; curve in SI Fig S7.",
         "Exact CD inflection Tm for HFt-MSH", "Check Cancer Gene Therapy SI Fig S7."),
        ("P2_Medium_SI_Matrix", "OA_ESI_Supplementary", "PMCID:PMC8816334 | DOI:10.1126/sciadv.abl7346", "2022",
         "T. maritima encapsulin 24-variant pore library S1X-S7X (Fig 1C, Fig 2B, Table S1-S3)",
         "9 variants named explicitly in XML text + 7 metastable + 2 additional purified thermostable variants shown only in Fig 1C/2B graphic & Sci Adv SI.",
         "Remaining ~15 variants from S1X-S7X matrix (identities & DLS profiles)", "Download Science Advances SI PDF for Table S1-S3."),
        ("P2_Medium_ESI", "OA_ESI_Supplementary", "PMCID:PMC5729161 | DOI:10.1021/acs.biochem.7b00296", "2017",
         "AfFtn A127R, D138K, E65R low-salt CD/DSC Supporting Information",
         "Table 1 reports 800 mM NaCl Tm (84/85/84/98 C); low-salt condition curves reside in ACS SI.",
         "Low-salt AfFtn mutant thermal curves", "Check ACS Biochemistry Supporting Info."),
        ("P2_Medium_PDF_Interleaved", "OA_Thesis_TableLayout", "DR-NTU:48015 (Page 150 Table 4.1)", "2011",
         "EcDps Robetta Virtual Ala Scan top 5 C3 interfacial residues (pos 17/19/21/79/82)",
         "PDF multi-column extraction interleaved columns on lines 4373-4384 of Table 4.1.",
         "5 unambiguous ddG_bind pairs on EcDps C3 (54 clean rows already extracted)", "Inspect Page 150 PDF visually if full Table 4.1 parity is desired."),
        ("P2_Medium_ESI", "OA_ESI_Supplementary", "PMCID:PMC2653062", "2009",
         "Human H-chain ferritin Variant #1 CD melting curves (Figs S1-S2)",
         "Main text states melting curves very similar to WT; raw curves in JBC SI Figs S1-S2.",
         "SI melting curve decimals", "Check JBC Supporting Info Figs S1-S2."),
    ]
    for item in oa_esi_gaps:
        add_gap(*item)

    # Priority 1 & 3: Closed DOIs from clean_closed_dois.txt + unpaywall_results.json
    closed_path = find_scratch_file("clean_closed_dois.txt")
    unpay_path = find_scratch_file("unpaywall_results.json")
    udata = json.load(open(unpay_path, encoding="utf-8"))
    unpay_map = udata if isinstance(udata, dict) else {x["doi"].lower(): x for x in udata if isinstance(x, dict) and x.get("doi")}

    high_therm_mut_keywords = re.compile(
        r"thermal\s+stab|thermostab|hyperthermostab|thermally\s+robust|gastric\s+stab|heat\s+sensit|"
        r"10\.1016/j\.foodchem\.2020\.128879|10\.1093/jb/mvv065|10\.1021/ja057069x|"
        r"10\.1039/c3cc40886h|10\.1021/acs\.biochem\.5c00744|10\.1021/acs\.bioconjchem\.6c00378|"
        r"10\.1021/acs\.jafc\.5c02232|10\.1016/j\.foodchem\.2021\.129089|10\.1016/j\.ijbiomac\.2021\.09\.067|"
        r"10\.1002/chem\.202300488|10\.1006/abbi\.1996\.0007|10\.1006/bbrc\.2001\.5931|"
        r"10\.1093/jb/mvm187|10\.1021/bi200207w|10\.1021/bi301189t|10\.1021/cb500431r|"
        r"10\.1021/bi900387t|10\.1021/acs\.biochem\.8b01068|10\.1021/acs\.jpcb\.1c00276|"
        r"10\.1021/acsabm\.2c00154|10\.1016/j\.ijbiomac\.2023\.126973|10\.1007/s11033-024-09636-w|"
        r"10\.1016/j\.ijbiomac\.2026\.150643|10\.1016/j\.ijbiomac\.2011\.05\.005|"
        r"10\.1021/jacs\.0c07835|10\.1002/anie\.201609517|10\.1002/anie\.201800516|"
        r"10\.1039/c6cc03108k|10\.1039/d1nr01632f|10\.1016/j\.ijbiomac\.2022\.01\.174|"
        r"10\.1038/nchembio\.1163|10\.1021/acsinfecdis\.3c00506",
        re.I,
    )

    seen_dois = set()
    section = "ferritin"
    for line in open(closed_path, encoding="utf-8"):
        st = line.strip()
        if "CLOSED OTHER-TITLE PAPERS" in st:
            section = "other"
            continue
        m = re.match(r"^(CF|CN)\d+:\s*([^\s|]+)\s*\|\s*([^|]*)\|\s*([^|]*)\|\s*(.*)$", st)
        if not m:
            continue
        doi = m.group(2).strip().lower()
        if doi in seen_dois:
            continue
        seen_dois.add(doi)
        meta = unpay_map.get(doi, {}) if isinstance(unpay_map, dict) else {}
        title = (meta.get("title") or m.group(5) or "").replace("\n", " ").strip()
        yr = str(meta.get("year") or m.group(3) or "").strip()
        jr = (meta.get("journal_name") or m.group(4) or "").strip()

        is_high = bool(high_therm_mut_keywords.search(doi + " " + title))
        if section == "other" and not is_high:
            # Check if still protein nanocage stability/ferritin-relevant or pure off-topic perovskite/ceramic
            if re.search(r"perovskite|chitosanase|barstar|endocytosis|protein cage|nanocage|superlattice|diiron|gold|protein", title, re.I):
                prio = "P3_Low_Related_OtherTitle"
            else:
                prio = "P4_OffTopic_Or_NoMeta_ClosedDOI"
        else:
            prio = "P1_High_Closed_ThermalMutant" if is_high else "P3_Medium_Closed_FerritinPaper"

        add_gap(
            prio,
            f"Closed_Paywalled_DOI_{section}",
            f"DOI:{doi}" + (f" ({jr})" if jr else ""),
            yr,
            title or f"[Metadata missing in Unpaywall for {doi}]",
            "Paywalled / closed access (Unpaywall is_oa=False).",
            "Potential mutant Tm / dTm / DSC / CD / nanoDSF data inside full text or SI" if is_high else "Needs screening for embedded mutant thermal characterization",
            "Fetch PDF via institutional subscription or interlibrary loan and inspect Materials/Tables.",
        )

    # Priority 3: 38 To Be Published / No-PMID-DOI PDB Mutant Entries
    rc_path = find_scratch_file("rcsb_mutants_clean.json")
    rc_data = json.load(open(rc_path, encoding="utf-8"))
    tbp_high_kw = re.compile(r"thermostab|thermal stability|gastric stability|4QY|S73G|D119|T4FY|S111|G40E|M52H", re.I)
    tbp_rows = []
    for r in rc_data:
        if not isinstance(r, dict):
            continue
        jr = str(r.get("journal") or "")
        pmid = r.get("pmid")
        doi = r.get("doi")
        if ("to be published" in jr.lower()) or (not pmid and not doi):
            pdb_id = r.get("pdb_id") or r.get("rcsb_id") or ""
            title = r.get("title") or ""
            muts = json.dumps(r.get("pdbx_mutations", []), ensure_ascii=False)
            tbp_rows.append((pdb_id, title, muts, jr or "No PMID/DOI"))
            prio = "P1_High_Unpublished_Thermostab_PDB" if tbp_high_kw.search(title) else "P3_Unpublished_Mutant_PDB"
            add_gap(
                prio,
                "RCSB_PDB_Unpublished_Or_NoPMID",
                f"PDB:{pdb_id}",
                "",
                title,
                f"Crystal/cryo-EM structure deposited in RCSB PDB without linked publication ({jr or 'no PMID/DOI'}).",
                f"Mutant annotation in PDB: {muts[:180]}",
                "Monitor RCSB PDB for linked paper release or contact depositing authors.",
            )
    assert len(tbp_rows) == 38, f"Expected 38 TBP/no-PMID PDB entries, got {len(tbp_rows)}"

    # Priority 2 & 4: Chinese / Patent / Preprint / Interactive DB portal gaps
    portal_gaps = [
        ("P1_High_DB_Portal", "Structured_DB_Portal", "ProThermDB (bulk ~31,000 entries)", "",
         "ProThermDB full release table", "Only targeted API/TSV queries checked (yielded 17 ferritin rows incl. 4 mislabeled EBNA1 fusions); bulk archive requires portal export.",
         "Check if synthetic/fusion ferritin entries under alternative names exist in full ~31k rows", "Download full ProThermDB dataset via browser and grep UniProt/Title for ferritin/Bfr/Dps/encapsulin."),
        ("P2_Medium_DB_Portal", "Structured_DB_Portal", "SoluProtMutDB / MPTherm / HotProtein", "",
         "SoluProtMutDB (JS SPA) & MPTherm & HotProtein", "Requires JS browser execution or full bulk download outside offline/static sandbox.",
         "Possible prokaryotic ferritin / membrane-associated Bfr / thermophilic ferritin entries", "Query portals or download bulk releases directly."),
        ("P2_Medium_DB_Portal", "Structured_DB_Portal", "BRENDA (EC 1.16.3.1 ferroxidase & EC 1.16.3.2)", "",
         "BRENDA Enzyme Thermal Stability & Engineered Mutants", "Requires authenticated portal/SOAP/JS access.",
         "Site-directed mutant temperature optimum / stability records", "Check BRENDA mutant/stability tables for EC 1.16.3.1."),
        ("P2_Medium_ZH_Literature", "Chinese_Journal_And_Thesis_DB", "CNKI / Wanfang / SPKX (食品科学) / 生物工程学报", "",
         "Chinese ferritin mutant & fusion literature (e.g., PMID:35470630, 40328702, 34472304; Crossref SPKX sample)",
         "CNKI/Wanfang full texts paywalled/captcha-gated; OpenAlex/Crossref Chinese titles returned 63 candidates.",
         "Chinese mutant Tm / 50-90 C heat-treatment survival / thesis tables", "Search CNKI/Wanfang with '铁蛋白 AND (突变 OR 融合) AND (热稳定性 OR Tm OR DSC)'."),
        ("P2_Medium_Patents", "Patent_FullText", "CNIPA / WIPO / USPTO / EPO (Europe PMC Patents returned 0)", "",
         "Ferritin engineering & vaccine nanocage thermal-stability patents",
         "Europe PMC patent search returned 0; commercial vaccine/food nanocarrier thermostability claims live in patent offices.",
         "Claims on thermostable ferritin mutant nanocages (e.g., SpFN, soybean/prawn/HuHF fusions)", "Search Lens.org / Google Patents / CNIPA for ferritin mutant Tm/thermostability claims."),
        ("P2_Medium_Preprints", "Preprint_FullText", "Europe PMC Preprints (hitCount=78, metadata only)", "",
         "78 bioRxiv / medRxiv / SSRN ferritin preprints from Europe PMC query",
         "Metadata counted; full PDFs of unreviewed preprints (except decoded 10.1101/2025.04.26.649701 etc.) not all downloaded.",
         "Recent 2024-2026 engineered ferritin & encapsulin thermostability preprints", "Batch-fetch remaining bioRxiv XMLs from Europe PMC list."),
    ]
    for item in portal_gaps:
        add_gap(*item)

    return gaps


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)
    src_rows = build_source_inventory()
    raw_records = make_records()
    unique_rows, records = build_unique_mutant_rollup(raw_records)
    gap_rows = build_inaccessible_sources()

    # Write CSVs
    inv_fields = list(src_rows[0].keys())
    rec_fields = list(records[0].keys())
    uniq_fields = list(unique_rows[0].keys())
    gap_fields = list(gap_rows[0].keys())

    write_csv(WORKSPACE_DIR / "database_and_literature_source_inventory.csv", inv_fields, src_rows)
    write_csv(WORKSPACE_DIR / "ferritin_mutant_thermostability_dataset.csv", rec_fields, records)
    write_csv(WORKSPACE_DIR / "unique_mutants_summary.csv", uniq_fields, unique_rows)
    write_csv(WORKSPACE_DIR / "inaccessible_and_follow_up_sources.csv", gap_fields, gap_rows)

    # Compute exact Dual-Count Reconciliation statistics from actual generated data
    mut_recs = [r for r in records if r["record_role"] == "Mutant"]
    wt_recs = [r for r in records if r["record_role"] == "WT_baseline"]

    wet_mut_recs = [r for r in mut_recs if r["evidence_tier"] != EV_D]
    comp_mut_recs = [r for r in mut_recs if r["evidence_tier"] == EV_D]
    quant_mut_recs = [r for r in mut_recs if r["evidence_tier"] == EV_A]
    quant_exact_num = [r for r in quant_mut_recs if r["value_numeric_C"]]

    wet_uniqs = [u for u in unique_rows if u["highest_evidence_tier"] != EV_D]
    comp_only_uniqs = [u for u in unique_rows if u["highest_evidence_tier"] == EV_D]
    comp_and_wet_overlap = sum(1 for u in unique_rows if u["highest_evidence_tier"] != EV_D and u["computational_predictions"])
    quant_uniqs = [u for u in unique_rows if u["highest_evidence_tier"] == EV_A]

    core_fams = {FAM_FERRITIN, FAM_BFR, FAM_DPS, FAM_FUSION}
    core_mut_recs = [r for r in mut_recs if r["family_tier"] in core_fams]
    core_wet_recs = [r for r in wet_mut_recs if r["family_tier"] in core_fams]
    core_uniqs = [u for u in unique_rows if u["family_tier"] in core_fams]
    core_wet_uniqs = [u for u in wet_uniqs if u["family_tier"] in core_fams]

    by_fam_rec = Counter(r["family_tier"] for r in mut_recs)
    by_fam_wet_rec = Counter(r["family_tier"] for r in wet_mut_recs)
    by_fam_uniq = Counter(u["family_tier"] for u in unique_rows)
    by_fam_wet_uniq = Counter(u["family_tier"] for u in wet_uniqs)

    by_ev_rec = Counter(r["evidence_tier"] for r in mut_recs)
    by_ev_uniq = Counter(u["highest_evidence_tier"] for u in unique_rows)

    by_cat_rec = Counter(r["mutation_category"] for r in mut_recs)
    by_cat_uniq = Counter(u["mutation_category"] for u in unique_rows)

    by_gap_prio = Counter(g["priority"] for g in gap_rows)

    fams_order = [FAM_FERRITIN, FAM_BFR, FAM_DPS, FAM_FUSION, FAM_EXT]
    evs_order = [EV_A, EV_B, EV_C, EV_D]
    cats_order = [CAT_POINT, CAT_MULTI, CAT_INDEL, CAT_FUSION]

    report_lines = [
        "# 铁蛋白（Ferritin）及其超家族突变体温度/热稳定性全网检索汇总报告",
        "",
        "> **核心结论（双口径统计 Dual-Count Reconciliation）**：主流蛋白质热稳定性结构化数据库（FireProtDB、ThermoMutDB、GeoStab、ProteinGym、Tsuboyama、MaveDB）中**几近空白（共 0 条天然铁蛋白点突变 Tm 数据，仅 ProThermDB 含 4 条 EBNA1 融合蛋白 DSC Tm 数据且被误标为 wild-type）**；真正可靠的铁蛋白突变热稳定性数据几乎全部散落在**开放获取文献、学位论文（如 NTU 硕士/博士论文 DR-NTU 48015）与付费期刊摘要/全文**中。",
        "",
        "---",
        "",
        "## 一、 总览与双口径统计核对表（WT 对照严格分离）",
        "",
        "| 统计层级 / 口径 | 条件级实验与计算记录数<br>(`ferritin_mutant_thermostability_dataset.csv`) | 去重唯一突变体数<br>(`unique_mutants_summary.csv`) | 说明 |",
        "| :--- | :---: | :---: | :--- |",
        f"| **突变体总计（含全部 5 类家族 + 计算预测）** | **{len(mut_recs)}** | **{len(unique_rows)}** | 不含 WT 基线对照 |",
        f"| ↳ 其中：**核心铁蛋白家族（Tier 1–4：典型 24-mer、Bfr、Dps、工程融合）** | **{len(core_mut_recs)}** | **{len(core_uniqs)}** | 剔除 Tier 5 扩展超家族（Encapsulin/Artemin/Dodecin/YciF） |",
        f"| **湿实验突变体总计（Tier A + B + C，剔除纯计算）** | **{len(wet_mut_recs)}** | **{len(wet_uniqs)}** | 含定量 Tm/Ti/Tagg/Tdiss、热耐受阈值、定性排序 |",
        f"| ↳ 其中：**核心铁蛋白湿实验（Tier 1–4，Tier A–C）** | **{len(core_wet_recs)}** | **{len(core_wet_uniqs)}** | 最严格的核心铁蛋白湿实验数据集 |",
        f"| **定量温度测量突变体（Tier A：DSC / CD / nanoDSF / Tycho / DSF / Rayleigh 等）** | **{len(quant_mut_recs)}**（其中精确单值数字 **{len(quant_exact_num)}**） | **{len(quant_uniqs)}** | 包括 ΔTm / Tdiss / Tagg_onset / 区间组报告 |",
        f"| **纯计算预测条目（Tier D：FoldX ΔΔG + Robetta 虚拟丙氨酸扫描）** | **{len(comp_mut_recs)}** | **{len(comp_only_uniqs)}**（纯计算独占；另 **{comp_and_wet_overlap}** 个已在湿实验表重复验证） | 来自 DR-NTU 48015 表 2.3/3.1/4.1 |",
        f"| **单独隔离的野生型基线对照（`record_role = WT_baseline`）** | **{len(wt_recs)}** | —（不计入突变体） | pH 扫描、亚基同源寡聚体、FireProtDB P02792 7 条 WT Tm 等 |",
        "",
        "---",
        "",
        "## 二、 按蛋白家族层级（A1 五大类独立标注）细分统计",
        "",
        "| 蛋白家族层级 (`family_tier`) | 突变体全部记录 | 湿实验突变体记录 | 全部唯一突变体 | 湿实验唯一突变体 |",
        "| :--- | :---: | :---: | :---: | :---: |",
    ]
    for fam in fams_order:
        report_lines.append(
            f"| `{fam}` | {by_fam_rec[fam]} | {by_fam_wet_rec[fam]} | {by_fam_uniq[fam]} | {by_fam_wet_uniq[fam]} |"
        )
    report_lines.extend([
        f"| **合计** | **{len(mut_recs)}** | **{len(wet_mut_recs)}** | **{len(unique_rows)}** | **{len(wet_uniqs)}** |",
        "",
        "---",
        "",
        "## 三、 按证据质量层级与突变类型拆分",
        "",
        "### 3.1 按证据层级 (`evidence_tier`)",
        "",
        "| 证据层级 | 条件级记录数 | 唯一突变体数（按最高证据归口） | 数据内涵 |",
        "| :--- | :---: | :---: | :--- |",
    ])
    ev_desc = {
        EV_A: "定量熔解/解聚/聚集温度（Tm, Ti, Tdiss, Tagg_onset, T_sub_CD, 明确 ΔTm）",
        EV_B: "热澄清 / 耐热存活温度阈值（如 45–100 °C 处理仍保持可溶/笼结构完整）",
        EV_C: "定性比较或热稳定性排序（>WT / <WT / ~WT / 75 °C 溶解度排名）",
        EV_D: "计算能量预测（FoldX ΔΔG、Robetta 虚拟 Ala 扫描 ΔΔG_bind）",
    }
    for ev in evs_order:
        report_lines.append(f"| `{ev}` | {by_ev_rec[ev]} | {by_ev_uniq[ev]} | {ev_desc[ev]} |")

    report_lines.extend([
        "",
        "### 3.2 按突变类型 (`mutation_category`)",
        "",
        "| 突变类型 | 条件级记录数 | 唯一突变体数 |",
        "| :--- | :---: | :---: |",
    ])
    for cat in cats_order:
        report_lines.append(f"| `{cat}` | {by_cat_rec[cat]} | {by_cat_uniq[cat]} |")

    report_lines.extend([
        "",
        "---",
        "",
        "## 四、 各大数据库与文献引擎检索审计汇总（共 17 个库/引擎）",
        "",
        "> 完整明细见 [database_and_literature_source_inventory.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/database_and_literature_source_inventory.csv)。",
        "",
        "- **结构化稳定性 / DMS / 互作数据库（零值已充分验证）**：",
        "  - `FireProtDB`（2025-09-20 全量 4.83 GB SQL 转储，41,146,172 行）：仅含 `seq_id=8955`（人 `P02792` FTL）**7 条 WT Tm**（与 `4 条 WT ΔG`，来自 `PMID:9681870`），`mutant_id IS NOT NULL` 为 **0**。",
        "  - `ThermoMutDB`（13,337 条 JSON）、`GeoStab`（S4346/S571/S669/S8754 共 14,340 行）、`ProteinGym`（217 substitutions + 66 indels，`A0A1I9GEU1_NEIME` 经核验为 PilE 菌毛素假阳性）、`Tsuboyama Mega-scale`（9,728 条）、`MaveDB`（2,161 条）：全部 **0 条** 铁蛋白突变数据。",
        "  - `SKEMPI 2.0`（7,084 数据行）：含 7 行 `4E6K_AB_G`（铜绿假单胞菌 `BfrB–Bfd`），仅测定 295 K 结合亲和力 `Kd`（SPR），无热变性 `Tm`。",
        "  - `ProThermDB`（缓存检索结果）：含 13+4=17 条相关行（全部为 `WT` 标签），其中 `NO=170, 171, 173, 174` 经检查 `REMARKS` 实为 **人 H-铁蛋白 N/C 端融合 EBV EBNA1（`E1F1`/`F1E1`）的合成构建体 DSC Tm（`70.5/80.1 °C` 与 `63.6/84.6 °C`，`PMID:33825471`）**，全库 ~31k 批量导出列为待补全项。",
        "  - `UniProtKB`（68 个条目 JSON）：含 1 条定性热稳定性下降记录（`P9WNE5` 结核分枝杆菌 `BfrB Δ168–181`，`PMID:21494619`）及多条寡聚体稳定性注释。",
        "  - `RCSB PDB`：检出 **461 个突变铁蛋白家族晶体/冷冻电镜结构**（其中 **38 个** 为 `To Be Published` 或无 PMID/DOI）。",
        "- **文献检索母集与全文深度核验**：",
        "  - `Europe PMC` 10 组检索去重后 **7,083 篇**（原始命中 18,384 条）；`OpenAlex` 10 组去重 **6,365 篇**；`PubMed XML` 11 批去重 **1,856 PMIDs**；`Crossref` 抽样 200 条（含中文食品科学期刊检索）。",
        "  - 完成 **269 篇不重复 PMC 全文**（231 XML + 16 NCBI XML + 22 HTML）、**18 篇 OA PDF/学位论文**（含新加坡南洋理工大学 `DR-NTU 48015` 全文及解码版 `pdf_1007_decoded.txt`）与 **105 篇非 PMC 热稳定性 PubMed 摘要**的逐篇核验；对 260 个 DOI 执行 `Unpaywall` OA 检测（**127 个 OA / 133 个闭源付费**）。",
        "",
        "---",
        "",
        "## 五、 方法学陷阱与关键数据质量警示（读库必读）",
        "",
        "1. **`ΔTm` 符号方向反转陷阱**：Yu Zhang 学位论文（`DR-NTU 48015` 表 2.7/3.3/4.3）及 `PMC6149950`（表 2）采用 **`ΔTm = Tm_WT - Tm_mut`** 定义（正值代表去稳定、负值代表更稳定！例如 `EcBfr R30A` `72.2 °C > WT 64.8 °C` 原文记作 `-7.4`）。本库同时保留原文 `reported_delta_value` 并统一换算 **`normalized_delta_Tm_C = mut - WT`**（正值=更耐热）。",
        "2. **跨源交叉印证与不重复计数**：`PMC2852946`（JBC 2010 `EcBfr`）与 `PMC3267954`（Protein Sci 2011 `EcDps`）在 PMC 中以图片形式嵌表，全文检索不到数值；同实验室 `DR-NTU 48015` 完整保留了表 2.7 与表 4.3 精确数字。已通过 `cross_ref_duplicate_group` 关联并在唯一突变体层面严格去重。",
        "3. **同一突变体跨文献差异**：`EcBfr D56R/E60H` 在 `DR-NTU 48015` 表 3.3 中 `Tm=60.8 °C`（`ΔTm = -7.4 °C`），而同组后续 `J. Biochem. 2015` 论文摘要（`PMID:26115686`）报道为比 WT 低 `16 °C`（体系/缓冲或测定条件不同），两条均已单独入库标注。",
        "4. **测定体系高度异质，严禁直接混池比较绝对 `Tm`**：例如大豆 `SSFH-2` 与 `HuHF-∇EP`（`PMC9069879`）使用**冻干固体粉末 DSC（`HuHF WT = 82 °C`）**；海莴苣 `UpFER`（`PMC3403414`）在 **1.5 M 盐酸胍（GdnHCl）** 中测 DSC；人 `H/L 异源多聚体`（`PMC6800815` vs `PMC9798250`）在水性 NanoDSC 下峰温可达 **92–115 °C**（且载铁 ~700 Fe 会使 L 富集型降低 12–14 °C、H 富集型降低 18–20 °C）。",
        "",
        "---",
        "",
        "## 六、 存在但暂不可获取 / 待补齐核对来源统计（共 " + str(len(gap_rows)) + " 项）",
        "",
        "> 完整逐条清单见 [inaccessible_and_follow_up_sources.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/inaccessible_and_follow_up_sources.csv)。",
        "",
        "| 待补齐类别 / 优先级 (`priority`) | 条目数 | 核心内容示例 |",
        "| :--- | :---: | :--- |",
    ])
    for pr, cnt in sorted(by_gap_prio.items()):
        report_lines.append(f"| `{pr}` | {cnt} | 见 CSV 对应明细及获取建议 |")

    report_lines.extend([
        "",
        "### 最值得优先补齐的高价值缺口",
        f"1. **已有定量摘要/明确突变热稳定性标题的闭源付费论文（{by_gap_prio['P1_High_Closed_ThermalMutant']} 篇高优先级）**：例如牡蛎铁蛋白 `D120G`（`89.17 °C`, `+13 °C`, `10.1016/j.foodchem.2020.128879`）、虾铁蛋白 `T158R/G159R` 与 `D132R`（`10.1016/j.foodchem.2021.129089`、`10.1016/j.ijbiomac.2021.09.067`）、牛蛙 M 铁蛋白 C4 孔道渐变突变系列（`~ -5 °C`, `10.1021/acs.biochem.5c00744`）、C3 界面重设计 `+20 °C` 笼蛋白（`10.1039/c3cc40886h`）、AF3 设计的无半胱氨酸高热稳 `rHuHF 9LUW`（`10.1021/acs.jafc.5c02232`）、`Dps` 疏水内腔系列 `Tm = 74–90 °C`（`10.1021/ja057069x`）、人 `Ecut-HFn`（`10.1016/j.ijbiomac.2023.126973`）、`EcBfr D56R/E60H` 全文（`10.1093/jb/mvv065`）及 `马脾/HuHF N23/芳环簇/E-helix` 系列（`10.1002/chem.202300488`, `10.1021/bi200207w`, `10.1021/bi301189t`）。",
        "2. **OA 全文的补充材料（ESI）与仅作图未列表数据**：`PMC11776510`（`MtBfrA E121X` RSC ESI Fig S13 精确小数）、`PMC9798250`（Fig 3b 中间 H:L 比例曲线值）、`PMC8816334`（Sci Adv 24-variant 完整 Table S1–S3）、`PMC11452399`（Fig S4g `FTH1-1` DSF 数据）。",
        "3. **RCSB PDB `38` 个未发表（`To Be Published` / 无 PMID）突变结构**：特别是标题直指热稳定的 `9BXP`（`PaFtnA`）、`9BYE`（`EcFtnA`）、`9LUW`（`HuHF`）、`8YKR`（计算辅助二硫键胃稳定性优化）、`7XRG` 等。",
        "4. **中文与专利库深度补全**：知网/万方《食品科学》《生物工程学报》硕博论文全文，以及 `ProThermDB` ~31k 全量包、`SoluProtMutDB`、`BRENDA` 门户检索。",
        "",
        "---",
        "",
        "## 七、 交付文件清单",
        "- [database_and_literature_source_inventory.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/database_and_literature_source_inventory.csv)",
        "- [ferritin_mutant_thermostability_dataset.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/ferritin_mutant_thermostability_dataset.csv)",
        "- [unique_mutants_summary.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/unique_mutants_summary.csv)",
        "- [inaccessible_and_follow_up_sources.csv](file:///Users/a26683/Desktop/ferritinMTdatabase/inaccessible_and_follow_up_sources.csv)",
        "- [README_database_summary_report.md](file:///Users/a26683/Desktop/ferritinMTdatabase/README_database_summary_report.md)",
    ])

    report_text = "\n".join(report_lines) + "\n"
    (WORKSPACE_DIR / "README_database_summary_report.md").write_text(report_text, encoding="utf-8")

    print("=== BUILD COMPLETE ===")
    print(f"Source inventory rows: {len(src_rows)}")
    print(f"Dataset total rows: {len(records)} (Mutant={len(mut_recs)}, WT_baseline={len(wt_recs)})")
    print(f"Unique mutants: {len(unique_rows)} (WetLab={len(wet_uniqs)}, CompOnly={len(comp_only_uniqs)}, QuantTemp={len(quant_uniqs)})")
    print(f"Core Ferritin (Tier 1-4): recs={len(core_mut_recs)} (WetLab={len(core_wet_recs)}), uniq={len(core_uniqs)} (WetLab={len(core_wet_uniqs)})")
    print(f"Gap / follow-up rows: {len(gap_rows)}")
    print(f"Report size bytes: {len(report_text.encode('utf-8'))}")


if __name__ == "__main__":
    main()
