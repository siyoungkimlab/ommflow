Getting started
===============

Requirements
------------

* Python 3.12 or later.
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

The name is optional and defaults to ``ommflow``.

The standard install
~~~~~~~~~~~~~~~~~~~~

``install.sh`` with no arguments builds the conda environment, which is the
complete installation including automatic GAFF ligand parameterization.
Equivalently, by hand:

.. code-block:: console

   conda env create -f environment.yml
   conda activate ommflow
   python -m pip install -e .

Conda is required because the OpenFF packages (``openff-toolkit``,
``openff-units``, ``openff-utilities``, ``openff-interchange``) and AmberTools
are not published to PyPI at all, and ``openmmforcefields`` declares no
dependencies of its own, so pip cannot assemble a working ligand stack.

The environment is over 200 packages, because ``openmmforcefields`` pulls in
``openff-toolkit``, which pulls ``openff-nagl`` and its PyTorch stack. Conda's
classic solver takes minutes on a dependency graph that size, where the
libmamba solver takes seconds: measured on one machine with the same specs and
cached repodata, 277 s against 4 s.

libmamba ships with conda 23.10 and later and is the default there, so usually
there is nothing to do. ``install.sh`` selects it explicitly in case conda has
been configured back to the classic solver, and it does not require the
separate ``mamba`` command. On an older conda:

.. code-block:: console

   conda update -n base -c conda-forge conda

Those timings cover the solve only. Downloading and extracting the packages
takes its own few minutes, which no solver changes.

Without conda
~~~~~~~~~~~~~

``install.sh --pip-only`` builds a plain venv from PyPI. It runs everything
except automatic GAFF ligands, so a DMS or MAE input carrying a small molecule
will be rejected. Use it only where conda is unavailable.

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
