# ommflow

`ommflow` prepares and runs explicit-solvent protein molecular dynamics with
OpenMM. It reads hydrogen-complete PDB structures and structural-only DMS or
MAE files, applies a selected protein and water force field, solvates and
neutralizes the system, equilibrates it, and runs restartable production MD.
Protein-only operation is the default. Optional automatic GAFF 2.11 ligand
parameterization is available for chemically annotated DMS and MAE inputs.

## Requirements

- Python 3.12 or later
- OpenMM 8.5 or later, with the desired platform support (CPU, CUDA, OpenCL,
  or Metal)

The protein force fields you can select depend on the parameter files your
OpenMM ships. OpenMM 8.5 is the first release carrying every combination in the
table below: 8.2 has no `amber19`, and 8.3 and 8.4 ship `charmm36_2024.xml`
without its `charmm36_2024/` water directory. Selecting a family your
installation cannot fully provide is reported by name, along with the families
it can.

## Installation

`install.sh` creates an environment, installs `ommflow` into it, reports what
the environment can do, and runs the tests:

```bash
bash install.sh                   # conda env "ommflow", with automatic ligands
bash install.sh myenv             # same, under another name
bash install.sh --pip-only        # venv in ./ommflow, protein and peptide only
```

The name is optional and defaults to `ommflow`. The two modes exist because the
dependency stack is split; the sections below describe what each one installs
and why.

### Protein and peptide runs

Protein-only runs, including peptide-binder early stopping, need nothing beyond
OpenMM and install with plain pip:

```bash
cd ommflow
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

### Automatic GAFF ligands

Automatic GAFF ligands need more than pip can provide. The OpenFF packages
(`openff-toolkit`, `openff-units`, `openff-utilities`, `openff-interchange`)
and AmberTools are not published to PyPI, and `openmmforcefields` declares no
dependencies of its own, so `pip install -e '.[ligands]'` cannot assemble a
working setup. Use the bundled conda environment instead:

```bash
conda env create -f environment.yml
conda activate ommflow
python -m pip install -e .
```

If you already have the OpenFF stack and AmberTools in your environment, the
`ligands` extra adds the remaining pip-installable pieces:

```bash
python -m pip install -e '.[ligands]'
```

On Linux, OpenMM's PyPI wheels require glibc 2.34 or newer. Enterprise
distributions with an older glibc (RHEL 8 and its derivatives, common on HPC
clusters) can only reach OpenMM 8.3 from PyPI, which is below the 8.5 this
package needs; install OpenMM from conda-forge there instead.

| Component | Needed for | Source |
|---|---|---|
| `openmm`, `numpy` | everything | PyPI or conda-forge |
| `openmmforcefields`, `rdkit` | automatic ligands | PyPI or conda-forge |
| `openff-toolkit`, `openff-interchange` | automatic ligands | conda-forge only |
| AmberTools (`antechamber`, `sqm`) | AM1-BCC charges for GAFF | conda-forge only |

To run the test suite:

```bash
python -m pip install -e '.[dev]'
pytest tests/
```

Run the installed CLI:

```bash
ommflow --help
```

Without installing, the same entry point is reachable as
`python -m ommflow.bin.ommflow`.

## Quick start

Start a new 100 ns production run using the defaults:

```bash
ommflow protein.pdb --workdir protein_md
```

This uses Amber ff19SB, OPC water, 298 K, 1 bar, a 0.1 ns NVT phase followed
by a 0.1 ns NPT phase, then runs production toward a 100 ns total target.

Specify a CHARMM36 2024 system:

```bash
ommflow protein.pdb \
  --workdir protein_charmm \
  --proteinff charmm36_2024 \
  --waterff tip3p \
  --production-ns 200
```

For CHARMM36 2024, `--waterff tip3p` resolves to
`charmm36_2024/water.xml`, the CHARMM-modified TIP3P water and ion parameter
set. It does not use Amber TIP3P parameters.

## Configuration files

Generate a ready-to-edit TOML template:

```bash
ommflow --write-default-config protein_md.toml
```

Run with it:

```bash
ommflow --config protein_md.toml
```

Explicit command-line values override TOML values:

```bash
ommflow \
  --config protein_md.toml \
  --production-ns 200
```

The precedence order is:

```text
Built-in defaults < TOML values < explicit CLI values
```

An example configuration is:

```toml
input_structure = "protein.pdb"
workdir = "protein_md"

proteinff = "amber19sb"
waterff = "opc"
ligand_mode = "auto" # Automatically parameterize eligible DMS/MAE ligands.
ligandff = "gaff-2.11"
padding_nm = 1.0
saltM = 0.15
temperature = 298.0
pressure = 1.0

equilibration_ns = 0.1
equilibration_report_interval_ns = 0.01
production_ns = 100.0
production_report_interval_ns = 1.0
checkpoint_interval_ns = 0.01
integration_fs = 2.0
hmr = false
seed = 0
precision = "mixed"
performance_interval_ns = 1.0

# Disabled by default. A single GAFF ligand is selected automatically.
early_stop = false
# monitor_ligand = "ligand-0" # Required when more than one ligand exists.
# monitor_chain = "B" # Input chain ID; selects its complete disconnected component.
# monitor_component = "component-2" # Exact durable component ID.
monitor_interval_ns = 0.1
pocket_cutoff_nm = 0.5
contact_cutoff_nm = 0.5
detach_cutoff_nm = 0.8
confirmation_checks = 2

# platform = "CUDA"
```

`hmr` enables hydrogen mass repartitioning: bonded hydrogens are raised to
4 amu with the mass taken from their heavy atom, which makes a 4 fs timestep
stable and roughly doubles throughput. It is off by default. Turning it on
changes the default `integration_fs` from 2 to 4; setting `integration_fs`
explicitly still wins. Rigid water is never repartitioned, for any water model,
so 3-, 4-, and 5-site models are all unaffected. Repartitioning is a change to
the dynamics, not just a speed setting: it slows the fastest motions to permit
the longer step, leaving equilibrium properties intact.

`precision` selects the GPU floating-point precision: `mixed` (the default),
`single`, or `double`. It applies only to platforms that expose a `Precision`
property, which excludes CPU and Reference. Not every platform can honor every
value — Apple's OpenCL rejects `mixed`, for instance — so the default falls
back to the platform's own precision and says so, while a precision you asked
for explicitly is an error when it cannot be met.

Every run prints the system it built and the platform it runs on before
integrating, including device name and the resolved precision:

```text
System: 50197 particles (50197 topology atoms, 12000 virtual sites)
  chains: 6, residues: 15000 (14800 water, 30 ion), bonds: 12000
  constraints: 24000, forces: HarmonicBondForce, NonbondedForce, ...
  box: 8.00 x 8.00 x 8.00 nm (512.0 nm^3)
  nonbonded: PME, cutoff 0.9 nm
Platform: CUDA (speed 100)
  DeviceIndex: 0
  DeviceName: NVIDIA A100
  Precision: mixed
```

Check that banner first if throughput is lower than expected: an unnoticed
fallback to `CPU` or `OpenCL` costs far more than any setting in this file.

Every `performance_interval_ns` of production, one row is appended to
`performance.csv` recording where wall time went:

```text
production_time_ns  step   wall_s  interval_s  ns_per_day  md_s    trajectory_s  ...  reported_pct
0.008               2000   7.915   7.915       87.329      7.880   0.003              0.444
0.016               4000   15.849  7.934       87.120      15.801  0.006              0.302
```

`ns_per_day` covers the interval since the previous row, the `*_s` columns are
cumulative seconds per task, and `reported_pct` is the share of wall time spent
outside integration. Timings cover the current submission only, so they restart
at zero on a resume. A summary prints when production ends:

```text
Production wall time: 39.6 s (0.22% in reporting and checkpointing)
  trajectory       0.02 s    0.04%
  state            0.03 s    0.07%
  checkpoint       0.01 s    0.03%
  monitor          0.03 s    0.09%
  integration     39.53 s   99.78%
```

If `reported_pct` is small, reporting is not what is limiting the run.

Time-based user settings use nanoseconds except `integration_fs`, which is in
femtoseconds. The default 2 fs timestep requires the script's hydrogen-bond
constraints and rigid-water settings.

## Input structures

The positional `input_structure` argument accepts one of:

| Extension | Reader | Required structural data |
|---|---|---|
| `.pdb` | OpenMM `PDBFile` | Standard PDB topology and coordinates |
| `.dms` | `DMSReader` | DMS SQLite `particle` and `bond` tables |
| `.mae` | `MAEReader` | MAE `m_atom` and `m_bond` tables |
| `.gro` | `GROReader` | GROMACS coordinates; standard protein bonds are inferred |

DMS and MAE inputs are treated as structural inputs only. Their embedded
force-field data, if any, is ignored; `ommflow` applies the selected OpenMM
XML protein and water force fields.

Crystallographic water and ions in the input are kept and are given any virtual
sites the selected water model needs before solvation, so four- and five-site
models such as the default OPC accept them.

Inputs must have hydrogens already added, standard residue and atom names
compatible with the selected force field, and appropriate protonation states.
For an input extracted from periodic simulation coordinates, make the protein
whole before running `ommflow`. Supplying `padding_nm` creates a new solvent
box and replaces any input PBC box vectors.

### Automatic ligands

By default, `ligand_mode = "auto"` identifies disconnected,
untemplated non-polymer components and parameterize qualifying components with
GAFF 2.11. This mode is available only for Amber protein force fields and
requires a `.dms` or `.mae` input: those formats retain the formal charges and
bond orders needed to construct an OpenFF molecule. PDB and GRO inputs are
rejected in auto mode rather than guessing chemistry. Install its optional
dependencies with `pip install -e '.[ligands]'`.

CHARMM ligand parameterization requires CGenFF and is not implemented.
Components in which a standard protein residue is covalently connected to a
nonstandard chemical residue are detected but rejected; covalent ligand
handling is deliberately detect-only. Inspect `components.json` in every new
run for source atom indices, residue and chain membership, charges,
classification, and standard/GAFF force-field provenance.

### Early-stop target detachment

Early stop is disabled by default. Set `early_stop = true` or pass
`--early-stop` to monitor a production target. A selector without `--early-stop`
is rejected rather than silently ignored. Automatic selection remains limited to
exactly one GAFF ligand candidate. A peptide binder is never inferred
automatically: choose exactly one selector, either the compatible legacy
`monitor_ligand` / `--monitor-ligand ligand-N`, an input chain with
`monitor_chain` / `--monitor-chain B`, or a durable component ID with
`monitor_component` / `--monitor-component component-2`. The chain selector
only accepts chain IDs from the input structure, not a solvated output; it
selects that chain's complete disconnected component (including a linked
multi-chain peptide). A chain mapping to zero or multiple components is an
error. Inspect `components.json` for `component-N`, chain/residue membership,
classification, and final production atom mappings.

`component-N` and `ligand-N` are assigned from the input structure alone,
before anything is built: components are numbered by their lowest input atom
index, and `ligand-N` numbers the ligand candidates in that same order. List
them without building or running anything:

```bash
ommflow complex.dms --list-components
```

```text
ID           LIGAND ID  CLASSIFICATION    CHAINS  RESIDUES  COMPOSITION
component-0  -          standard          A       12        ACE, ALA, NME
    --early-stop --monitor-chain A
component-1  -          standard          B       12        ACE, ALA, NME
    --early-stop --monitor-chain B
component-2  ligand-0   ligand_candidate  L       1         BNZ
    --early-stop --monitor-ligand ligand-0
```

Each entry shows the selector that names it unambiguously, preferring a chain
ID when only one component uses that chain. Water and ion components are
counted rather than listed. The IDs are stable for a given input file, and a
selector is resolved right after classification, so a wrong one fails within
seconds and lists the valid choices.

A peptide binder needs no ligand parameterization, so it works from any
supported input, including `.pdb`:

```bash
ommflow complex.pdb \
  --workdir binder_md \
  --early-stop \
  --monitor-chain B
```

The target's own component is excluded from the pocket, so the receptor
chain(s) form the pocket for a peptide target exactly as the protein does for a
GAFF ligand. Selecting a component that holds every non-water protein atom is
rejected.

The production-start protein pocket excludes the selected target, water/ions,
and all GAFF ligand candidates, and is saved in `pocket.json`. At each
checkpoint-aligned `monitor_interval_ns`, the workflow uses PBC-aware direct
heavy-atom distances. It stops only after `confirmation_checks` consecutive
measurements have zero contacts below `contact_cutoff_nm` *and* a minimum
ligand-pocket distance above `detach_cutoff_nm`. `monitor.csv` records each
measurement and `status.json` records `running`, `target_reached`, or
`detached`.

The `ligand_id` field in `pocket.json` and `status.json` names the target the
way you selected it: `ligand-N` for a GAFF ligand, and the `component-N` ID for
a peptide or any other standard-force-field target.

`contact_fraction` is relative to the single contact count measured at the
start of production, so it exceeds 1 whenever the target settles into more
contacts than it began with. Only `contact_count` and
`min_ligand_pocket_distance_nm` decide detachment.

## Force fields and water models

The protein force-field choice determines the only water/ion XML directory
that will be used:

| `proteinff` | Protein parameters | Water/ion directory |
|---|---|---|
| `amber14sb` | `amber14/protein.ff14SB.xml` | `amber14/` |
| `amber15ipq` | `amber14/protein.ff15ipq.xml` | `amber14/` |
| `amber19sb` | `amber19/protein.ff19SB.xml` | `amber19/` |
| `charmm36` | `charmm36.xml` | `charmm36/` |
| `charmm36_2024` | `charmm36_2024.xml` | `charmm36_2024/` |

Supported Amber water models are `opc`, `opc3`, `spce`, `tip3p`, `tip3pfb`,
`tip4pew`, and `tip4pfb`. Supported CHARMM water models are `tip3p`,
`tip3p-pme-b`, `tip3p-pme-f`, `spce`, `tip4p2005`, `tip4pew`, `tip5p`, and
`tip5pew`. Invalid protein/water combinations are rejected before simulation
setup.

## Output files

A new run creates the following inside `workdir`:

| File | Contents |
|---|---|
| `input.<extension>` | Original supplied structure |
| `solvated.pdb` | Solvated and neutralized system |
| `equilibration.dcd` | NVT and NPT equilibration trajectory |
| `equilibration.csv` | Equilibration state data |
| `equilibrated.pdb` | Final NPT-equilibrated coordinates |
| `trajectory.dcd` | Production trajectory |
| `state.csv` | Production state data; step and time both start at 0 |
| `final.pdb` | Latest production coordinates |
| `checkpoint.chk` | Current production OpenMM checkpoint |
| `system.xml` | Production OpenMM system definition |
| `integrator.xml` | Production integrator definition |
| `final.toml` | Resolved settings and production target |
| `components.json` | Covalent-component classifications and parameter provenance |
| `pocket.json` | Immutable early-stop target/pocket selection, when enabled |
| `performance.csv` | Wall-time breakdown per production task |
| `monitor.csv` | Early-stop measurements and confirmation state, when enabled |
| `status.json` | Early-stop production outcome, when enabled |

`final.toml` has a stable setting order, making it suitable for comparison
with a generated default configuration.

## Restarts and wall time

A non-existent or empty `workdir` starts a new run, so a directory pre-created
by a batch scheduler is not mistaken for an unresumable one. An existing
`workdir` holding a previous run automatically attempts a restart; no
`--resume` option is needed.

For a wall-time-limited job, submit the same command again:

```bash
ommflow protein.pdb \
  --workdir protein_md \
  --production-ns 200
```

The checkpoint restores the production simulation state and the command
appends to `trajectory.dcd` and `state.csv`. `production_ns` is an absolute
target, not the duration of an individual submission. Therefore repeated
submissions run only the remaining time needed to reach 200 ns.

A restart reloads every setting saved in `workdir/final.toml` — reporting and
checkpoint intervals, the platform, and all early-stop settings — unless the
same setting is supplied again on the command line or in `--config`. This keeps
`trajectory.dcd` and `state.csv` on one cadence across submissions. `final.toml`
is rewritten only after a restart passes every check, so a rejected resume
leaves the run directory usable.

When no explicit `--production-ns` is supplied during a restart, the saved
target in `workdir/final.toml` is used. Increase the target with an absolute
value:

```bash
ommflow \
  --workdir protein_md \
  --production-ns 300
```

This updates the saved target to 300 ns and runs only the remaining production
time. Use the same OpenMM version and, where possible, the same execution
platform when restoring a checkpoint.

For safe output appending, `checkpoint_interval_ns` must be no greater than
and divide `production_report_interval_ns` exactly in integration steps. For
example, these are valid:

```toml
production_report_interval_ns = 1.0
checkpoint_interval_ns = 0.01
```

Equilibration and production reporting use independent intervals.

For early-stop runs, `pocket.json` is restored rather than recomputed;
`monitor.csv` appends on restart. A confirmed `detached` run does not continue
automatically under the same early-stop target. Disable `--early-stop` or set
a larger explicit `--production-ns` target to continue deliberately.

## License

MIT. See [LICENSE](LICENSE).
