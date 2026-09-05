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

Protein and peptide runs
~~~~~~~~~~~~~~~~~~~~~~~~

Protein-only runs, including early stopping on a peptide binder, need nothing
beyond OpenMM and install with plain pip:

.. code-block:: console

   cd ommflow
   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -e .

Automatic GAFF ligands
~~~~~~~~~~~~~~~~~~~~~~

Automatic GAFF ligand parameterization from DMS or MAE input needs more than
pip can provide. The OpenFF packages (``openff-toolkit``, ``openff-units``,
``openff-utilities``, ``openff-interchange``) and AmberTools are not published
to PyPI, and ``openmmforcefields`` declares no dependencies of its own, so
``pip install -e '.[ligands]'`` cannot assemble a working setup. Use the
bundled conda environment:

.. code-block:: console

   conda env create -f environment.yml
   conda activate ommflow
   python -m pip install -e .

If the OpenFF stack and AmberTools are already present, the ``ligands`` extra
adds the remaining pip-installable pieces.

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
