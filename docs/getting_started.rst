Getting started
===============

Requirements
------------

* Python 3.12, since ``install.sh`` pins numpy below 2, which has no builds
  for 3.13 or later.
* OpenMM 8.5 or later, with the required execution platform support: CPU,
  CUDA, OpenCL, or Metal.

OpenMM 8.5 is the first release carrying every force field and water model
``ommflow`` offers: 8.2 has no ``amber19``, and 8.3 and 8.4 ship
``charmm36_2024.xml`` without its ``charmm36_2024/`` water directory.

Installation
------------

``install.sh`` creates an environment, installs ``ommflow`` into it, reports
what that environment can do, and runs the tests:

.. code-block:: console

   bash install.sh                   # conda env "ommflow", with automatic ligands
   bash install.sh myenv             # same, under another name
   bash install.sh --pip-only        # venv in ./ommflow, protein and peptide only
   bash install.sh --cuda 12         # conda env "ommflow", OpenMM for CUDA 12

The name is optional and defaults to ``ommflow``.

The standard install
~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   bash install.sh

That builds a conda environment, adds the OpenFF stack with pip, installs
``ommflow``, reports what the environment can do, and runs the tests.

The dependency stack is split deliberately. conda supplies the pieces that are
not on PyPI or that must match the machine's CPU: Python, numpy, OpenMM and
AmberTools. pip then adds the OpenFF packages and RDKit with ``--no-deps``,
because conda-forge's ``openff-toolkit`` depends on ``openff-nagl``, which
pulls PyTorch and roughly eighty packages ommflow never uses; GAFF charges come
from AmberTools AM1-BCC, not from a neural-network model. The split takes the
environment from 234 packages to 165.

``environment.yml`` alone is therefore **not** a complete install. Use
``install.sh``, or run its pip steps by hand after creating the environment.

The OpenFF packages are pinned to the versions conda-forge's own solver selects
together, so the set is known to be mutually consistent and an upstream commit
cannot change the environment underneath you.

``numpy`` is held below 2, and Python therefore at 3.12, because numpy 1.x has
no builds for 3.13 or later. numpy 2.4 and later are built against an
x86-64-v2 baseline and abort on older HPC nodes, and AmberTools' bundled tools
pin ``numpy<2`` themselves.

The install finishes by assigning real AM1-BCC charges to ethanol through
``sqm``, so a broken AmberTools is reported at install time rather than on the
first ligand.

AmberTools' bundled Python tools (``ndfes``, ``fetkutils``, ``proprep``)
declare requirements that AmberTools' conda package does not install, and pip
reports each missing one when it adds the OpenFF stack. ommflow uses none of
these tools, but ``environment.yml`` supplies their requirements anyway
(``netcdf4``, ``pdb2pqr``, ``requests``, and ``biopython<1.86``) so the
install finishes without spurious errors.

CUDA
~~~~

By default ``install.sh`` chooses no CUDA release, because many installs are on
machines without an NVIDIA GPU. conda then picks one from the driver on the
machine running the install, or the newest release when there is no driver.
On a cluster that is usually a login node, and GPU nodes with an older driver
then fail at the first run with
``CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)``.

Pass ``--cuda`` with the CUDA release to build for, the same choice OpenMM's
own instructions make with ``cuda-version=12`` and ``openmm[cuda12]``:

.. code-block:: console

   bash install.sh --cuda 12                # conda: pins cuda-version=12
   bash install.sh --pip-only --cuda 12     # pip: installs openmm[cuda12]

The GPU nodes' driver must support the CUDA that gets installed: compare
``conda list cuda-version`` in the environment with ``CUDA Version`` in
``nvidia-smi`` on a GPU node. The conda path pins ``cuda-version`` and sets
``CONDA_OVERRIDE_CUDA``, so it works from a login node with no GPU; a release
conda-forge has no OpenMM build for is reported by conda's solver. Rerunning
with ``--cuda`` also repins an existing environment.

The check at the end of the install lists the CUDA platform only where an
NVIDIA driver is present. Without one it says why CUDA did not load; confirm on
a GPU node with ``python -m openmm.testInstallation``.

Without conda
~~~~~~~~~~~~~

``install.sh --pip-only`` builds a plain venv from PyPI. It runs everything
except automatic GAFF ligands, so a DMS or MAE input carrying a small molecule
will be rejected. Use it only where conda is unavailable.

It applies the same ``numpy<2`` pin as the conda path, so ``python3`` must be
Python 3.12. ``pyproject.toml``
leaves numpy unbounded, because that is what the library needs rather than what
a given machine can run, so the constraint lives in the installer. Installing
ommflow with bare pip instead of ``install.sh`` therefore takes the newest
numpy, which aborts on hosts without x86-64-v2.

Run the installed CLI:

.. code-block:: console

   ommflow --help

Start a run
-----------

Start a new run with the default 100 ns production target:

.. code-block:: console

   ommflow protein.pdb --workdir protein_md

The defaults use Amber ff19SB with OPC water, 298 K, 1 bar, 0.1 ns each of
NVT and NPT equilibration, and a 2 fs integration timestep.

For a CHARMM36 2024 system:

.. code-block:: console

   ommflow protein.pdb \
     --workdir protein_charmm \
     --proteinff charmm36_2024 \
     --waterff tip3p \
     --production-ns 200
