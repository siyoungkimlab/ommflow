Inputs and outputs
==================

Input structures
----------------

The positional ``input_structure`` argument accepts:

.. list-table::
   :header-rows: 1

   * - Extension
     - Reader
     - Required structural data
   * - ``.pdb``
     - OpenMM ``PDBFile``
     - Standard PDB topology and coordinates
   * - ``.dms``
     - ``DMSReader``
     - DMS SQLite ``particle`` and ``bond`` tables
   * - ``.mae``
     - ``MAEReader``
     - MAE ``m_atom`` and ``m_bond`` tables
   * - ``.gro``
     - ``GROReader``
     - GROMACS coordinates; standard protein bonds are inferred

DMS, MAE, and GRO structures are read as topology and coordinates. Embedded
force-field parameters are ignored, and the selected OpenMM XML force fields
are applied. DMS and MAE additionally retain atom formal charges and bond
orders when their input fields are present.

Inputs need hydrogen atoms, standard residue/atom names compatible with the
selected force field, and appropriate protonation states. Make a protein whole
before using coordinates extracted from a periodic trajectory. The workflow
creates a new solvent box from ``padding_nm`` and does not preserve input PBC
box vectors.

Automatic ligands
-----------------

With ``ligand_mode = "auto"``, PDB and GRO inputs remain supported when all
components match the selected standard force field. A nonstandard ligand
candidate requires ``.dms`` or ``.mae`` chemistry data.
They provide the explicit formal charges and bond orders required to build
ligand OpenFF molecules. PDB and GRO inputs are rejected in auto mode rather
than inferring chemistry. A missing DMS/MAE charge or bond order is an error
only when the corresponding component is selected as an automatic ligand.

The workflow classifies covalently connected components, not individual
residues. Disconnected non-polymer components that fail the selected standard
force-field template match are ligand candidates. Components joining standard
protein residues to nonstandard chemistry are written as
``covalent_candidate`` and rejected: covalent ligand support is detect-only.
See ``components.json`` for original atom indices, residue and chain details,
known formal charge, classification, and parameterization provenance. Every
monitorable disconnected component has a stable ``component-N`` ID. Ligand
candidates additionally have stable ``ligand-N`` IDs; all mapped components
include final solvated production-topology atom indices.

Output files
------------

Each work directory contains:

.. list-table::
   :header-rows: 1

   * - File
     - Contents
   * - ``input.<extension>``
     - Original supplied structure
   * - ``solvated.pdb``
     - Solvated and neutralized system
   * - ``solvated.mae``
     - The same system as Maestro atoms, bonds, and box vectors
   * - ``equilibration.dcd`` / ``equilibration.csv``
     - NVT and NPT equilibration trajectory and state data
   * - ``equilibrated.pdb``
     - Final NPT-equilibrated structure
   * - ``equilibrated.mae``
     - The same structure as Maestro atoms, bonds, and box vectors
   * - ``trajectory.dcd`` / ``state.csv``
     - Production trajectory and state data
   * - ``final.pdb``
     - Latest production coordinates
   * - ``final.mae``
     - The same coordinates as Maestro atoms, bonds, and box vectors
   * - ``checkpoint.chk``
     - Current production checkpoint
   * - ``system.xml`` / ``integrator.xml``
     - Production system and integrator definitions
   * - ``final.toml``
     - Resolved settings and production target
   * - ``components.json``
     - Covalent-component classification and selected force-field provenance
   * - ``performance.csv``
     - Wall-time breakdown per production task
   * - ``dihedral_restraints.csv``
     - Restrained torsions, atom indices and reference angles (only with ``dihedral_restraint``)
   * - ``dihedral_restraints.png``
     - The restraint potential (only with ``dihedral_restraint``)
   * - ``pocket.json``
     - Immutable early-stop target/pocket selection (only with ``early_stop``)
   * - ``monitor.csv``
     - Early-stop heavy-atom measurements and confirmation state (only with ``early_stop``)
   * - ``status.json``
     - Production outcome: ``running``, ``target_reached``, or ``detached``

Each MAE file holds the same atoms in the same order as the PDB written beside
it, plus every bond and the periodic box. A viewer can therefore load
``solvated.mae`` as the topology for ``equilibration.dcd`` and
``equilibrated.mae`` for ``trajectory.dcd``, while ``final.mae`` stands on its
own as the latest production frame. A PDB drops the bonds of anything without
a standard residue template and cannot number past 99,999 atoms, both of which
a solvated box reaches. Bond orders from a ``.dms`` or ``.mae`` input are
carried through solvation into the written files; water, ions, and added
hydrogens are single bonds.

Production time and step both start at zero after equilibration, so
``state.csv`` steps line up with the steps recorded in ``monitor.csv``. The
last time in
``state.csv`` therefore gives the cumulative production time, including all
restart segments.
