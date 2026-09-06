#!/usr/bin/env bash
# Install ommflow.
#
# Two paths, because the dependency stack is split:
#   (default)   the full install, including automatic GAFF ligand
#               parameterization. This is what you want. conda supplies the
#               binary pieces that are not on PyPI (AmberTools, OpenMM) and a
#               CPU-safe numpy; pip --no-deps then adds the OpenFF stack,
#               which conda-forge would otherwise burden with PyTorch.
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
    # tool, so prefer it, then the mamba command, then whatever conda defaults
    # to.
    #
    # Advertising libmamba is not the same as being able to load it: a conda
    # whose libarchive has moved on prints "Error while loading conda entry
    # point" and then refuses CONDA_SOLVER=libmamba outright. Detect that,
    # because forcing the solver there turns a slow install into a failed one.
    SOLVER=conda
    conda_help=$(conda create --help 2>&1 || true)
    if [ -n "${OMMFLOW_SOLVER:-}" ]; then
        export CONDA_SOLVER="$OMMFLOW_SOLVER"
        echo "==> Solving with conda's $OMMFLOW_SOLVER solver (OMMFLOW_SOLVER)"
    elif printf '%s' "$conda_help" | grep -q "Error while loading conda entry point: conda-libmamba-solver"; then
        if command -v mamba >/dev/null 2>&1; then
            SOLVER=mamba
            echo "==> conda cannot load its libmamba solver; using mamba instead"
        else
            echo "Warning: conda advertises the libmamba solver but cannot load it,"
            echo "         so this solve falls back to the classic solver and may"
            echo "         take several minutes. Usually a broken libarchive:"
            echo "             conda install -n base -c conda-forge libarchive"
            echo "         or install mamba, or set OMMFLOW_SOLVER to override."
            echo
        fi
    elif printf '%s' "$conda_help" | grep -q libmamba; then
        # Set even though it is the modern default, in case this conda is
        # configured back to the classic solver.
        export CONDA_SOLVER=libmamba
        echo "==> Solving with conda's libmamba solver"
    elif command -v mamba >/dev/null 2>&1; then
        SOLVER=mamba
        echo "==> This conda has no libmamba solver; solving with mamba instead"
    else
        echo "Warning: this conda predates the libmamba solver, so the solve"
        echo "         below may take several minutes. To fix it permanently:"
        echo "             conda update -n base -c conda-forge conda"
        echo
    fi

    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "==> Updating existing environment '$ENV_NAME' (this can take several minutes)"
        "$SOLVER" env update -n "$ENV_NAME" -f "$REPO_DIR/environment.yml" --prune
    else
        echo "==> Creating environment '$ENV_NAME' from environment.yml"
        echo "    165 conda packages, then the OpenFF stack from pip."
        "$SOLVER" env create -n "$ENV_NAME" -f "$REPO_DIR/environment.yml"
    fi
    conda activate "$ENV_NAME"
    PYTHON="$(command -v python)"

    # The OpenFF half, installed with pip rather than conda. conda-forge's
    # openff-toolkit depends on openff-nagl, a neural-network charge model
    # ommflow never calls, which drags in PyTorch and roughly 80 other
    # packages. --no-deps takes the OpenFF packages without that subtree, so
    # their pure-python requirements are listed here explicitly.
    echo "==> Installing the OpenFF ligand stack"
    "$PYTHON" -m pip install --quiet --upgrade pip wheel setuptools versioningit
    "$PYTHON" -m pip install --quiet \
        networkx cachetools "xmltodict<=1.0.2" python-constraint \
        "pydantic>=2,<2.12" pint lxml pyyaml tinydb validators rdkit

    # Pinned to the versions conda-forge's own solver selects together, so the
    # set is known to be mutually consistent, and an upstream commit cannot
    # change this environment underneath you.
    GH=https://github.com/openforcefield
    for package in \
        "$GH/openff-utilities.git@v0.1.18" \
        "$GH/openff-units.git@0.4.0" \
        "$GH/openff-toolkit.git@0.19.0" \
        "$GH/openff-forcefields.git@2026.01.0" \
        "$GH/openff-interchange.git@v0.5.4" \
        "https://github.com/openmm/openmmforcefields.git@0.16.0"
    do
        "$PYTHON" -m pip install --quiet --no-deps "git+${package}"
    done

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

missing = None
try:
    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import SystemGenerator  # noqa: F401
    from rdkit import Chem  # noqa: F401
except ImportError as error:
    missing = error.name

if missing is not None:
    print(f"  automatic ligands   NOT available ({missing} is missing)")
    print("                      protein and peptide runs are unaffected;")
    print("                      rerun install.sh without --pip-only for GAFF")
elif not (shutil.which("antechamber") and shutil.which("sqm")):
    print("  automatic ligands   INCOMPLETE: AmberTools binaries not on PATH")
    print("                      antechamber and sqm are needed for AM1-BCC")
else:
    # Importing is not proof: assign real AM1-BCC charges through sqm, which
    # is what GAFF parameterization does for every ligand.
    try:
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges("am1bcc")
        total = sum(float(charge.m) for charge in molecule.partial_charges)
        print("  automatic ligands   available (GAFF 2.11 with AM1-BCC charges)")
        print(f"                      verified on ethanol, net charge {total:+.3f}")
    except Exception as error:  # noqa: BLE001 - report whatever went wrong
        print(f"  automatic ligands   BROKEN: AM1-BCC failed ({error})")

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
