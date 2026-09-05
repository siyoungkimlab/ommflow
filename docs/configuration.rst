Configuration
=============

Generate a template
-------------------

Generate a documented TOML settings file:

.. code-block:: console

   ommflow --write-default-config protein_md.toml

Run with that file:

.. code-block:: console

   ommflow --config protein_md.toml

Explicit CLI values override TOML values:

.. code-block:: console

   ommflow \
     --config protein_md.toml \
     --production-ns 200

Setting precedence is:

.. code-block:: text

   Built-in defaults < TOML values < explicit CLI values

Example
-------

.. code-block:: toml

   input_structure = "protein.pdb"
   workdir = "protein_md"

   proteinff = "amber19sb"
   waterff = "opc"
   ligand_mode = "auto"  # Automatically parameterize eligible DMS/MAE ligands.
   ligandff = "gaff-2.11"
   padding_nm = 1.0
   # cutoff_nm defaults to 0.9 for Amber and 1.2 for CHARMM.
   saltM = 0.15
   temperature = 298.0
   pressure = 1.0

   equilibration_ns = 0.1
   equilibration_report_interval_ns = 0.01
   production_ns = 100.0
   production_report_interval_ns = 1.0
   checkpoint_interval_ns = 0.01
   performance_interval_ns = 1.0
   integration_fs = 2.0
   # Hydrogen mass repartitioning. When true, integration_fs defaults to 4.0.
   hmr = false
   seed = 0

   # GPU floating-point precision: mixed (default), single, or double.
   # Left commented so the default can fall back on a platform that cannot
   # honor it; setting it here makes an unsupported precision an error.
   # precision = "mixed"

   # Disabled by default. A single GAFF ligand is selected automatically.
   early_stop = false
   # monitor_ligand = "ligand-0"  # Required when more than one ligand exists.
   # monitor_chain = "B"  # Input chain ID; selects its complete component.
   # monitor_component = "component-2"  # Exact durable component ID.
   monitor_interval_ns = 0.1
   pocket_cutoff_nm = 0.5
   contact_cutoff_nm = 0.5
   detach_cutoff_nm = 0.8
   confirmation_checks = 2

   # platform = "CUDA"

All user-facing durations use nanoseconds except ``integration_fs``, which is
in femtoseconds. The CLI uses kebab-case equivalents, such as
``--production-ns`` and ``--integration-fs``.

``saltM``, ``temperature``, and ``pressure`` are the current TOML names.
For migration, ``salt_molarity``, ``temperature_k``, and ``pressure_bar`` are
accepted as aliases, but an alias cannot appear alongside its current name.

Ligands
-------

Automatic ligand handling is the default: ``ligand_mode = "auto"``. It uses
``ligandff = "gaff-2.11"`` for qualifying ligands. Auto mode requires an
Amber protein force field and chemically annotated DMS or MAE input. It
rejects PDB and GRO input rather than guessing missing formal charges or bond
orders. Its optional dependencies are not available from PyPI alone; see
:doc:`getting_started` for the conda environment that provides them.

Hydrogen mass repartitioning
----------------------------

``hmr`` / ``--hmr`` raises bonded hydrogens to 4 amu, taking the mass from
their heavy atom, which makes a 4 fs timestep stable and roughly doubles
throughput. It is off by default. Enabling it changes the default
``integration_fs`` from 2 to 4; an explicit ``integration_fs`` still wins.
Rigid water is never repartitioned, so every 3-, 4-, and 5-site water model is
unaffected, as are the virtual sites of the four- and five-site models. This
alters the dynamics rather than only the speed: the fastest motions are slowed
to permit the longer step, leaving equilibrium properties intact.

Timing breakdown
----------------

Every ``performance_interval_ns`` of production, a row is appended to
``performance.csv`` recording cumulative seconds per task (trajectory, state,
checkpoint, monitor), the time left in integration, the rate over the interval
in ns/day, and the share of wall time spent outside integration. A summary of
the same breakdown prints when production ends. Timings cover the current
submission only and restart at zero on a resume.

Platform and precision
----------------------

``precision`` / ``--precision`` selects the GPU floating-point precision:
``mixed`` (the default), ``single``, or ``double``. It applies only to
platforms exposing a ``Precision`` property, which excludes CPU and Reference.
A platform need not support every value; Apple's OpenCL rejects ``mixed``. The
default therefore falls back to the platform's own precision and reports the
fallback, while a precision requested explicitly is an error when it cannot be
met.

``platform`` / ``--platform`` pins the OpenMM platform. Left unset, the fastest
platform that can create a context is used.

Before integrating, every run prints the system it built (particles, chains,
residues, waters and ions, constraints, forces, box, nonbonded method and
cutoff) and the platform it runs on with all of its properties, including the
device name and the resolved precision. Check that banner first when
throughput is lower than expected: an unnoticed fallback to a slower platform
costs far more than any other setting.

Early-stop target detachment
----------------------------

Early stopping is opt-in with ``early_stop = true`` or ``--early-stop``. It
monitors one production target during production. Automatic selection is
limited to exactly one GAFF ligand candidate. Peptide binders and receptors
are never inferred automatically. Select exactly one target with the legacy
``monitor_ligand`` / ``--monitor-ligand ligand-N``, with an input chain ID
using ``monitor_chain`` / ``--monitor-chain B``, or with a durable component
ID using ``monitor_component`` / ``--monitor-component component-2``.

A selector without ``early_stop`` is rejected rather than silently ignored.

``component-N`` and ``ligand-N`` are assigned from the input structure alone,
before anything is built: components are numbered by their lowest input atom
index, and ``ligand-N`` numbers the ligand candidates in that same order. List
them, with the selector that names each one, without building or running
anything:

.. code-block:: console

   ommflow complex.dms --list-components

Water and ion components are counted rather than listed, and a selector is
resolved immediately after classification, so a wrong one fails within seconds
and reports the valid choices.

``monitor_chain`` applies only to chain IDs in the input structure. It resolves
the complete disconnected component containing that chain, so selecting one
chain of a linked multi-chain peptide monitors the whole peptide. Zero or
multiple matching components is an error. ``components.json`` lists every
monitorable disconnected component as ``component-N``, including standard
force-field protein components, with its chains, residues, classification, and
production topology atom mapping.

At production start, the monitor saves a protein heavy-atom pocket consisting
of atoms within ``pocket_cutoff_nm`` of the target. The pocket excludes the
target, water/ions, and all GAFF ligand candidates. At each
``monitor_interval_ns``, it uses minimum-image, heavy-atom pair distances.
It stops only after ``confirmation_checks`` consecutive checks have both zero
pairs within ``contact_cutoff_nm`` and a minimum pocket distance greater than
``detach_cutoff_nm``. A target with no initial pocket is rejected as unbound.

The monitor interval must be an exact integration-step duration and an integer
multiple of ``checkpoint_interval_ns``. Thus every measurement has a matching
checkpoint.
