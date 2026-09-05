"""Build a solvated OpenMM system for a new workflow run."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import openmm as mm
from openmm import app, unit

from ommflow.lib.config import HYDROGEN_MASS_AMU
from ommflow.lib.forcefields import create_forcefield, resolve_forcefields
from ommflow.lib.ligands import (
    classify_components,
    component_summary,
    create_gaff_forcefield,
    create_ligand_modeller,
    ligand_candidates,
    production_component_atom_indices,
    reject_covalent_candidates,
    write_components,
)
from ommflow.lib.readers import read_structure


def _report_repartitioned_masses(topology: app.Topology, system: mm.System) -> None:
    """Confirm repartitioning left every atom with a usable mass.

    OpenMM moves mass onto a hydrogen from its heavy atom, skipping rigid water
    and any hydrogen bonded to a virtual site, so every water model is left
    alone. A heavy atom carrying many hydrogens still loses the most mass, so
    check none was drained to the point of instability.
    """
    atoms = list(topology.atoms())
    repartitioned = 0
    drained = []
    for index in range(system.getNumParticles()):
        if system.isVirtualSite(index):
            continue
        mass = system.getParticleMass(index).value_in_unit(unit.amu)
        atom = atoms[index] if index < len(atoms) else None
        if atom is not None and atom.element is not None:
            if atom.element.symbol == "H" and mass > 1.5:
                repartitioned += 1
            elif mass < 1.0:
                drained.append(
                    f"{atom.residue.name} {atom.residue.id} {atom.name} "
                    f"({mass:.3f} amu)"
                )
    if drained:
        raise ValueError(
            "Hydrogen mass repartitioning left these atoms below 1 amu, which "
            "would be unstable: " + ", ".join(drained[:5])
        )
    print(
        f"Hydrogen mass repartitioning: {repartitioned} hydrogens raised to "
        f"{HYDROGEN_MASS_AMU:g} amu (rigid water unchanged)."
    )


def copy_input_structure(input_structure: Path, output_dir: Path) -> Path:
    """Copy the source structure into the run directory."""
    copied_input = output_dir / f"input{input_structure.suffix.lower()}"
    if input_structure.resolve() != copied_input.resolve():
        shutil.copy2(input_structure, copied_input)
    return copied_input


def describe_input_components(args: argparse.Namespace) -> str:
    """Classify an input structure and report its selectable targets.

    This reads the structure only: nothing is solvated and no work directory is
    created, so the durable IDs can be looked up before committing to a run.
    """
    if args.input_structure is None:
        raise ValueError("input_structure must be supplied to list components.")
    if not args.input_structure.is_file():
        raise FileNotFoundError(
            f"Input structure does not exist: {args.input_structure}"
        )
    structure = read_structure(args.input_structure)
    protein_xml, water_directory, water_xml, _ = resolve_forcefields(
        args.proteinff, args.waterff
    )
    forcefield = create_forcefield(
        protein_xml, f"{water_directory}/{water_xml}", args.proteinff
    )
    components = classify_components(structure, forcefield)
    return (
        f"Components of {args.input_structure} "
        f"({args.proteinff} / {args.waterff}):\n\n"
        + component_summary(components)
    )


def build_solvated_system(
    args: argparse.Namespace,
    output_dir: Path,
    solvated_pdb: Path,
    components_json: Path | None = None,
) -> tuple[app.Modeller, mm.System, mm.LangevinMiddleIntegrator, object, object]:
    """Copy input, solvate it, and create the OpenMM system and integrator."""
    copied_input = copy_input_structure(args.input_structure, output_dir)
    structure = read_structure(copied_input)
    protein_xml, water_directory, water_xml, box_model = resolve_forcefields(
        args.proteinff, args.waterff
    )
    water_forcefield_xml = f"{water_directory}/{water_xml}"
    forcefield = create_forcefield(
        protein_xml, water_forcefield_xml, args.proteinff
    )
    components = classify_components(structure, forcefield)
    write_components(
        components_json or output_dir / "components.json",
        components,
        proteinff=args.proteinff,
        waterff=args.waterff,
        ligand_mode=args.ligand_mode,
        ligandff=args.ligandff,
    )
    modeller = app.Modeller(structure.topology, structure.positions)
    source_order = tuple(range(structure.topology.getNumAtoms()))
    candidates = ()
    if args.ligand_mode == "auto":
        reject_covalent_candidates(components)
        candidates = ligand_candidates(components)
        if candidates:
            if args.input_structure.suffix.lower() not in {".dms", ".mae"}:
                raise ValueError(
                    "Automatic ligand parameterization requires a .dms or .mae "
                    "input because PDB and GRO inputs do not retain formal "
                    "charges and bond orders."
                )
            if not args.proteinff.startswith("amber"):
                raise ValueError(
                    "ligand_mode=auto with GAFF 2.11 is supported only with "
                    "Amber protein force fields. CHARMM ligand parameterization "
                    "requires CGenFF, which is not implemented."
                )
            modeller, molecules, source_order = create_ligand_modeller(
                structure, candidates
            )
            forcefield = create_gaff_forcefield(
                protein_xml, water_forcefield_xml, molecules, args.ligandff
            )
    presolvation_topology = modeller.topology
    # Four- and five-site water models need their virtual sites before
    # addSolvent builds a System from the input, so any crystallographic water
    # already present matches its template.
    modeller.addExtraParticles(forcefield)
    modeller.addSolvent(
        forcefield,
        model=box_model,
        padding=args.padding_nm * unit.nanometer,
        neutralize=True,
        ionicStrength=args.saltM * unit.molar,
    )
    modeller.addExtraParticles(forcefield)
    write_components(
        components_json or output_dir / "components.json",
        components,
        proteinff=args.proteinff,
        waterff=args.waterff,
        ligand_mode=args.ligand_mode,
        ligandff=args.ligandff,
        production_topology_indices=production_component_atom_indices(
            modeller.topology, components, presolvation_topology, source_order
        ),
    )
    with solvated_pdb.open("w") as handle:
        app.PDBFile.writeFile(
            modeller.topology, modeller.positions, handle, keepIds=True
        )

    temperature = args.temperature * unit.kelvin
    pressure = args.pressure * unit.bar
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=args.cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        hydrogenMass=HYDROGEN_MASS_AMU * unit.amu if args.hmr else None,
    )
    if args.hmr:
        _report_repartitioned_masses(modeller.topology, system)
    integrator = mm.LangevinMiddleIntegrator(
        temperature, 1.0 / unit.picosecond, args.integration_fs * unit.femtoseconds
    )
    integrator.setRandomNumberSeed(args.seed)
    return modeller, system, integrator, temperature, pressure
