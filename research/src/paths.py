"""Centralized filesystem paths for ferritin thermostability benchmarking & modeling."""
from __future__ import annotations

import os
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
RESEARCH_ROOT = WORKSPACE_ROOT / "research"
DATA_DIR = RESEARCH_ROOT / "data"
RESULTS_DIR = RESEARCH_ROOT / "results"
MODELS_DIR = RESEARCH_ROOT / "models"
REPORTS_DIR = RESEARCH_ROOT / "reports"

FERRO_DATASET_CSV = WORKSPACE_ROOT / "ferritin_mutant_thermostability_dataset.csv"
UNIQUE_MUTANTS_CSV = WORKSPACE_ROOT / "unique_mutants_summary.csv"

EXTERNAL_GEODTM_ROOT = Path(os.environ.get("FERRITIN_EXTERNAL_ROOT", str(WORKSPACE_ROOT / "external")))
GEODTM_PYTHON = Path(os.environ.get("GEODTM_PYTHON", str(EXTERNAL_GEODTM_ROOT / ".venv-geodtm/bin/python")))
SKLEARN_PYTHON = Path(os.environ.get("SKLEARN_PYTHON", str(EXTERNAL_GEODTM_ROOT / ".venv/bin/python")))
RUN_GEODTM_SCRIPT = Path(os.environ.get("RUN_GEODTM_SCRIPT", str(EXTERNAL_GEODTM_ROOT / "scripts/run_geodtm.py")))
FOLDX_BIN = Path(os.environ.get("FOLDX_BIN", str(EXTERNAL_GEODTM_ROOT / "software/foldx/foldx5.1")))
FOLDX_MOLECULES_DIR = Path(os.environ.get("FOLDX_MOLECULES_DIR", str(WORKSPACE_ROOT / "molecules")))

GEOSTAB_ROOT = Path(os.environ.get("GEOSTAB_ROOT", str(EXTERNAL_GEODTM_ROOT / "models/GeoStab")))
S571_CSV = GEOSTAB_ROOT / "data/dTm/S571.csv"
S4346_CSV = GEOSTAB_ROOT / "data/dTm/S4346.csv"
ESM2_WEIGHTS = Path(os.environ.get("ESM2_WEIGHTS", str(EXTERNAL_GEODTM_ROOT / "models/esm2/esm2_t33_650M_UR50D.pt")))
GEODTM_CKPT = Path(os.environ.get("GEODTM_CKPT", str(GEOSTAB_ROOT / "model_dTm_3D/model.pt")))

PINNED_SHA256 = {
    "S4346": "285ab904d18dc73a386a7dbe9b96e1f4752529f53a2e73f9ded528f482cb4c49",
    "S571": "cad1a416cc0e6873e524398eef40dd398152d38c7db65feb05f351c99aec34a8",
}
