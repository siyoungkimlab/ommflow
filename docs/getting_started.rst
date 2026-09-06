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

``numpy`` is capped below 2.3 because 2.4 and later are built against an
x86-64-v2 baseline and abort on older HPC nodes.

The install finishes by assigning real AM1-BCC charges to ethanol through
``sqm``, so a broken AmberTools is reported at install time rather than on the
first ligand.

pip reports a dependency conflict naming ``proprep``, ``ndfes``, ``fetkutils``
and ``edgembar``, which are AmberTools' own bundled tools pinned to
``numpy<2``. ommflow uses none of them and the install is unaffected.

Without conda
~~~~~~~~~~~~~

``install.sh --pip-only`` builds a plain venv from PyPI. It runs everything
except automatic GAFF ligands, so a DMS or MAE input carrying a small molecule
will be rejected. Use it only where conda is unavailable.

It applies the same ``numpy<2.3`` cap as the conda path. ``pyproject.toml``
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
