Command-line reference
======================

Every option below has a TOML equivalent in snake_case, so
``--production-report-interval-ns`` is ``production_report_interval_ns``.
Explicit command-line values override a ``--config`` file, which overrides
the built-in defaults. See :doc:`configuration` for what each setting means.

``ommflow INPUT_STRUCTURE [options]``, where ``INPUT_STRUCTURE`` is a
hydrogen-complete structure in PDB, DMS, MAE, or GRO format. It is optional
when resuming an existing ``--workdir``.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Option
     - Default
     - Description
   * - ``--config``
     - --
     - TOML settings file. Explicit CLI values override its settings.
   * - ``--write-default-config``
     - --
     - Write a default TOML settings template to FILE and exit.
   * - ``--list-components``
     - ``false``
     - Print the components of INPUT_STRUCTURE with the selector that names each one, then exit without building or running anything.
   * - ``--workdir``
     - ``openmm_md``
     - Output directory.
   * - ``--padding-nm``
     - ``1.0``
     - Minimum protein-to-box-edge distance in nm (default: 1.0).
   * - ``--cutoff-nm``
     - --
     - Nonbonded cutoff in nm (default: 0.9 Amber, 1.2 CHARMM).
   * - ``--saltM``
     - ``0.15``
     - Added NaCl concentration in molar (default: 0.15).
   * - ``--proteinff``
     - ``amber19sb``
     - Protein force field (default: amber19sb).
   * - ``--waterff``
     - ``opc``
     - Water and ion parameter set (default: opc).
   * - ``--ligand-mode``
     - ``auto``
     - Ligand handling: auto (default) or disabled for strict standard force-field-only preparation.
   * - ``--ligandff``
     - ``gaff-2.11``
     - Automatic ligand force field (default: gaff-2.11).
   * - ``--temperature``
     - ``298.0``
     - Temperature in K (default: 298).
   * - ``--pressure``
     - ``1.0``
     - Pressure in bar (default: 1).
   * - ``--equilibration-ns``
     - ``0.1``
     - NVT and NPT equilibration duration each in ns (default: 0.1).
   * - ``--equilibration-report-interval-ns``
     - ``0.01``
     - Equilibration reporting interval in ns (default: 0.01).
   * - ``--production-ns``
     - ``100.0``
     - Target total NPT production time in ns (default: 100).
   * - ``--production-report-interval-ns``
     - ``1.0``
     - Production reporting interval in ns (default: 1).
   * - ``--checkpoint-interval-ns``
     - ``0.01``
     - Checkpoint writing interval in ns (default: 0.01).
   * - ``--performance-interval-ns``
     - ``1.0``
     - Timing-breakdown interval in ns for performance.csv (default: 1).
   * - ``--integration-fs``
     - ``2.0``
     - Integration timestep in fs (default: 2.0).
   * - ``--hmr``, ``--no-hmr``
     - ``false``
     - Hydrogen mass repartitioning (default: off). Raises bonded hydrogens to 4 amu, taking the mass from their heavy atom, and makes --integration-fs default to 4 instead of 2. Rigid water is left untouched for every water model.
   * - ``--seed``
     - ``0``
     - Random seed.
   * - ``--early-stop``, ``--no-early-stop``
     - ``false``
     - Stop production after confirmed target detachment (default: disabled).
   * - ``--monitor-ligand``
     - --
     - Exact GAFF ligand-N ID from components.json to monitor.
   * - ``--monitor-chain``
     - --
     - Input-structure chain ID whose one disconnected component to monitor.
   * - ``--monitor-component``
     - --
     - Exact component-N ID from components.json to monitor.
   * - ``--monitor-interval-ns``
     - ``0.1``
     - Target-monitoring interval in ns (default: 0.1).
   * - ``--pocket-cutoff-nm``
     - ``0.5``
     - Initial target-to-protein pocket cutoff in nm (default: 0.5).
   * - ``--contact-cutoff-nm``
     - ``0.5``
     - Target-pocket heavy-atom contact cutoff in nm (default: 0.5).
   * - ``--detach-cutoff-nm``
     - ``0.8``
     - Minimum target-pocket distance for detachment in nm (default: 0.8).
   * - ``--confirmation-checks``
     - ``2``
     - Consecutive detached checks required to stop (default: 2).
   * - ``--dihedral-restraint``
     - ``none``
     - Restrain protein backbone phi and psi to the input structure: bb for the whole backbone, ss for residues in helices and sheets only (needs MDTraj), none to disable (default).
   * - ``--dihedral-restraint-kJ``
     - ``20.0``
     - Dihedral restraint strength in kJ/mol (default: 20).
   * - ``--precision``
     - ``mixed``
     - GPU floating-point precision (default: mixed). Applies only to platforms that expose it; a platform that cannot honor the default falls back to its own precision with a warning.
   * - ``--platform``
     - --
     - Optional OpenMM platform; defaults to OpenMM's automatic selection.
