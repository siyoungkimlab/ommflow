#!/usr/bin/env bash
# Install ommflow.
#
# Two paths, because the dependency stack is split:
#   (default)   a conda-forge environment from environment.yml: the full
#               install, including automatic GAFF ligand parameterization.
#               This is what you want. The OpenFF packages and AmberTools are
#               not published to PyPI, so conda is the only way to get them.
#   --pip-only  a plain venv from PyPI, for environments without conda at all.
#               Everything except automatic GAFF ligands, so a DMS or MAE input
#               carrying a small molecule will not run.
#
# Usage:
#   bash install.sh [NAME] [--pip-only]
#
# NAME defaults to "ommflow". It names the conda environment, or the venv
# directory under --pip-only.
#
#   bash install.sh                  # conda env "ommflow"
#   bash install.sh myenv            # conda env "myenv"
#   bash install.sh --pip-only       # venv in ./ommflow
#   bash install.sh myenv --pip-only # venv in ./myenv
set -euo pipefail

ENV_NAME=""
PIP_ONLY=0
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [ $# -gt 0 ]; do
    case "$1" in
        --pip-only) PIP_ONLY=1; shift ;;
        -h|--help)  awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' \
                        "${BASH_SOURCE[0]}"; exit 0 ;;
        -*) echo "Unknown option: $1 (try --help)" >&2; exit 2 ;;
        *)
            [ -z "$ENV_NAME" ] || {
                echo "Only one environment name is accepted (got '$ENV_NAME' and '$1')." >&2
                exit 2
            }
            ENV_NAME="$1"; shift ;;
    esac
done
ENV_NAME="${ENV_NAME:-ommflow}"
VENV_DIR="$ENV_NAME"

# ----------------------------------------------------------------------
# 1. Create the environment and install ommflow into it.
# ----------------------------------------------------------------------
if [ "$PIP_ONLY" -eq 1 ]; then
    # requires-python is >=3.12; check before creating a venv pip would reject.
    python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)' || {
        echo "Python 3.12 or later is required; found $(python3 -V 2>&1)." >&2
        exit 1
    }
    echo "==> Creating venv at $VENV_DIR"
    python3 -m venv "$VENV_DIR"
    PYTHON="$VENV_DIR/bin/python"
    "$PYTHON" -m pip install --upgrade --quiet pip
    echo "==> Installing ommflow (protein-only; no automatic GAFF ligands)"
    "$PYTHON" -m pip install -e "$REPO_DIR[dev]"
    ACTIVATE="source $VENV_DIR/bin/activate"
else
    command -v conda >/dev/null 2>&1 || {
        echo "conda not found. Install Miniforge, or rerun with --pip-only." >&2
        exit 1
    }
    # conda activate is a shell function, unavailable in a non-interactive
    # script until this profile script is sourced.
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"

    # openmmforcefields pulls in openff-toolkit, which pulls openff-nagl and
    # its PyTorch stack, so this environment is well over 200 packages. Conda's
    # classic solver takes minutes on a graph that size; libmamba takes
    # seconds. libmamba ships with conda 23.10 and later and needs no separate
    # tool, so prefer it and only fall back to the mamba CLI.
    if conda create --help 2>&1 | grep -q libmamba; then
        # Set even though it is the modern default, in case this conda is
        # configured back to the classic solver.
        export CONDA_SOLVER=libmamba
        SOLVER=conda
        echo "==> Solving with conda's libmamba solver"
    elif command -v mamba >/dev/null 2>&1; then
        SOLVER=mamba
        echo "==> This conda has no libmamba solver; solving with mamba instead"
    else
        SOLVER=conda
        echo "Warning: this conda predates the libmamba solver, so the solve"
        echo "         below may take several minutes. To fix it permanently:"
        echo "             conda update -n base -c conda-forge conda"
        echo "         or, without changing conda's version:"
        echo "             conda install -n base -c conda-forge conda-libmamba-solver"
        echo "             conda config --set solver libmamba"
        echo
    fi

    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "==> Updating existing environment '$ENV_NAME' (this can take several minutes)"
        "$SOLVER" env update -n "$ENV_NAME" -f "$REPO_DIR/environment.yml" --prune
    else
        echo "==> Creating environment '$ENV_NAME' from environment.yml"
        echo "    Over 200 packages to solve and download."
        "$SOLVER" env create -n "$ENV_NAME" -f "$REPO_DIR/environment.yml"
    fi
    conda activate "$ENV_NAME"
    PYTHON="$(command -v python)"
    echo "==> Installing ommflow"
    "$PYTHON" -m pip install -e "$REPO_DIR"
    ACTIVATE="conda activate $ENV_NAME"
fi

# ----------------------------------------------------------------------
# 2. Verify what was actually installed, and say plainly what is missing.
#    A missing ligand stack is reported, not treated as a failed install:
#    --pip-only never has one, and it is optional either way.
# ----------------------------------------------------------------------
echo
echo "==> Verifying"
"$PYTHON" - <<'PY'
import shutil
import sys

import openmm

from ommflow.lib.forcefields import (
    FORCE_FIELD_FAMILIES,
    installed_force_field_families,
)

print(f"  python              {sys.version.split()[0]}")
print(f"  openmm              {openmm.version.version}")

families = installed_force_field_families()
missing = sorted(set(FORCE_FIELD_FAMILIES) - set(families))
print(f"  force fields        {', '.join(families) or 'none'}")
if missing:
    print(f"  !! incomplete       {', '.join(missing)} (needs OpenMM 8.5 or later)")

platforms = sorted(
    (openmm.Platform.getPlatform(i) for i in range(openmm.Platform.getNumPlatforms())),
    key=lambda platform: platform.getSpeed(),
    reverse=True,
)
print(f"  platforms           {', '.join(p.getName() for p in platforms)}")

try:
    from openff.toolkit.topology import Molecule  # noqa: F401
    from openmmforcefields.generators import SystemGenerator  # noqa: F401
    from rdkit import Chem  # noqa: F401
except ImportError as error:
    print(f"  automatic ligands   NOT available ({error.name} is missing)")
    print("                      protein and peptide runs are unaffected;")
    print("                      rerun install.sh without --pip-only for GAFF")
else:
    if shutil.which("antechamber") and shutil.which("sqm"):
        print("  automatic ligands   available (GAFF 2.11 with AM1-BCC charges)")
    else:
        print("  automatic ligands   INCOMPLETE: AmberTools binaries not on PATH")
        print("                      antechamber and sqm are needed for AM1-BCC")

if missing:
    sys.exit(1)
PY

# ----------------------------------------------------------------------
# 3. Run the test suite when pytest is present.
# ----------------------------------------------------------------------
if "$PYTHON" -c "import pytest" >/dev/null 2>&1; then
    echo
    echo "==> Running tests"
    "$PYTHON" -m pytest "$REPO_DIR/tests" -q
fi

echo
echo "Done. Activate with:  $ACTIVATE"
echo "Then try:             ommflow --help"
