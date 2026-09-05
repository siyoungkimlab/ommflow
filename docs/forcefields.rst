Force fields and water models
=============================

The selected protein force field determines the only water/ion XML directory
that ``ommflow`` will use.

.. list-table::
   :header-rows: 1

   * - ``proteinff``
     - Protein parameters
     - Water/ion directory
   * - ``amber14sb``
     - ``amber14/protein.ff14SB.xml``
     - ``amber14/``
   * - ``amber15ipq``
     - ``amber14/protein.ff15ipq.xml``
     - ``amber14/``
   * - ``amber19sb``
     - ``amber19/protein.ff19SB.xml``
     - ``amber19/``
   * - ``charmm36``
     - ``charmm36.xml``
     - ``charmm36/``
   * - ``charmm36_2024``
     - ``charmm36_2024.xml``
     - ``charmm36_2024/``

Amber supports ``opc``, ``opc3``, ``spce``, ``tip3p``, ``tip3pfb``,
``tip4pew``, and ``tip4pfb``. CHARMM supports ``tip3p``, ``tip3p-pme-b``,
``tip3p-pme-f``, ``spce``, ``tip4p2005``, ``tip4pew``, ``tip5p``, and
``tip5pew``.

For ``--proteinff charmm36_2024 --waterff tip3p``, the workflow loads
``charmm36_2024/water.xml``: the CHARMM-modified TIP3P water and ion model,
not Amber TIP3P. Incompatible protein/water selections are rejected before
system preparation.

The nonbonded cutoff defaults to 0.9 nm for Amber and 1.2 nm for CHARMM.
Override the selected family default with ``--cutoff-nm`` or ``cutoff_nm`` in
the TOML configuration.

Automatic ligand force fields
-----------------------------

``ligand_mode = "auto"`` uses GAFF 2.11 through ``openmmforcefields``
``SystemGenerator`` alongside the selected Amber protein and water XML files.
It is therefore supported only with ``amber14sb``, ``amber15ipq``, or
``amber19sb``. CHARMM ligand parameterization requires CGenFF and is not
implemented by ``ommflow``. Install the optional ligand stack before using
auto mode; it is not available from PyPI alone, so use the bundled conda
environment (see :doc:`getting_started`):

.. code-block:: console

   conda env create -f environment.yml

GAFF charge assignment also needs the AmberTools ``antechamber`` and ``sqm``
binaries, which ``environment.yml`` installs.

The default ``ligand_mode = "auto"`` retains the normal OpenMM
``ForceField`` workflow without importing these optional packages.
