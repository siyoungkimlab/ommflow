"""Ligand component detection and optional GAFF parameterization helpers."""

from __future__ import annotations

from dataclasses import dataclass
import json
from numbers import Real
from pathlib import Path

from openmm import app, unit

from ommflow.lib.readers import BondMetadata, StructureData


COVALENT_SCOPE = "detect-only"
PROTEIN_RESIDUE_NAMES = frozenset(
    {
        "ACE",
        "ALA",
        "ARG",
        "ASH",
        "ASN",
        "ASP",
        "CYM",
        "CYS",
        "CYX",
        "GLH",
        "GLN",
        "GLU",
        "GLY",
        "HID",
        "HIE",
        "HIP",
        "HIS",
        "ILE",
        "LEU",
        "LYN",
        "LYS",
        "MET",
        "NME",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)
_STANDARD_PROTEIN_RESIDUES = PROTEIN_RESIDUE_NAMES
STANDARD_SOLVENT_AND_ION_RESIDUES = frozenset(
    {
        "BR",
        "CA",
        "CL",
        "CS",
        "F",
        "HOH",
        "I",
        "K",
        "LI",
        "MG",
        "NA",
        "OPC",
        "OPC3",
        "RB",
        "SOL",
        "SPC",
        "SPCE",
        "TIP4",
        "TIP5",
        "TIP3",
        "TP3",
        "TP4",
        "TP5",
        "WAT",
        "ZN",
    }
)
_STANDARD_SOLVENT_AND_IONS = STANDARD_SOLVENT_AND_ION_RESIDUES


@dataclass(frozen=True)
class Component:
    """One covalently connected component of a structure."""

    index: int
    topology_atom_indices: tuple[int, ...]
    original_atom_indices: tuple[int, ...]
    residues: tuple[dict[str, object], ...]
    classification: str
    reason: str
    standard_template_match: bool | None
    net_known_formal_charge: int | float | None
    known_formal_charge_sum: int | float
    formal_charge_complete: bool
    ligand_id: str | None = None
    atom_formal_charges: tuple[tuple[int, int | float | None], ...] = ()

    def as_dict(
        self,
        proteinff: str,
        waterff: str,
        ligand_mode: str,
        ligandff: str,
        production_topology_indices: dict[str, tuple[int, ...]] | None = None,
    ) -> dict[str, object]:
        """Serialize a component with force-field provenance."""
        ligand_forcefield = (
            ligandff
            if ligand_mode == "auto" and self.classification == "ligand_candidate"
            else None
        )
        return {
            "component_index": self.index,
            "component_id": f"component-{self.index}",
            "ligand_id": self.ligand_id,
            "original_atom_indices": list(self.original_atom_indices),
            "topology_atom_indices": list(self.topology_atom_indices),
            "classification": self.classification,
            "residues": list(self.residues),
            "input_residue_labels": list(self.residues),
            "atom_formal_charges": [
                {
                    "original_atom_index": original_index,
                    "formal_charge": formal_charge,
                }
                for original_index, formal_charge in self.atom_formal_charges
            ],
            "chains": sorted({str(residue["chain_id"]) for residue in self.residues}),
            "net_known_formal_charge": self.net_known_formal_charge,
            "known_formal_charge_sum": self.known_formal_charge_sum,
            "formal_charge_complete": self.formal_charge_complete,
            "forcefield": {
                "protein": proteinff,
                "water": waterff,
                "standard_template_match": self.standard_template_match,
                "ligand": ligand_forcefield,
            },
            "monitorable": is_monitorable(self),
            "reason": self.reason,
            "production_topology_atom_indices": (
                list(
                    production_topology_indices.get(
                        f"component-{self.index}",
                        production_topology_indices.get(self.ligand_id, ()),
                    )
                )
                if (
                    production_topology_indices is not None
                    and (
                        f"component-{self.index}" in production_topology_indices
                        or self.ligand_id in production_topology_indices
                    )
                )
                else None
            ),
        }


def classify_components(
    structure: StructureData,
    forcefield=None,
) -> tuple[Component, ...]:
    """Classify covalent components and record standard-template provenance."""
    atoms = {atom.topology_index: atom for atom in structure.atom_metadata}
    adjacency = {atom_index: set() for atom_index in atoms}
    for bond in structure.bond_metadata:
        adjacency[bond.atom1_topology_index].add(bond.atom2_topology_index)
        adjacency[bond.atom2_topology_index].add(bond.atom1_topology_index)

    unmatched_residues: set[int] | None = None
    if forcefield is not None:
        unmatched_residues = {
            id(residue)
            for residue in forcefield.getUnmatchedResidues(structure.topology)
        }
    topology_atoms = list(structure.topology.atoms())
    remaining = set(atoms)
    connected_components = []
    while remaining:
        pending = [min(remaining)]
        component_atoms = set()
        while pending:
            atom_index = pending.pop()
            if atom_index in component_atoms:
                continue
            component_atoms.add(atom_index)
            remaining.discard(atom_index)
            pending.extend(adjacency[atom_index] - component_atoms)
        connected_components.append(tuple(sorted(component_atoms)))

    components = []
    ligand_number = 0
    for component_index, atom_indices in enumerate(connected_components):
        metadata = [atoms[index] for index in atom_indices]
        residue_names = {atom.residue_name.upper() for atom in metadata}
        residue_details = _residue_details(metadata)
        contains_protein = bool(residue_names & _STANDARD_PROTEIN_RESIDUES)
        only_standard_solvent_or_ions = residue_names <= _STANDARD_SOLVENT_AND_IONS
        component_residues = {
            id(topology_atoms[index].residue) for index in atom_indices
        }
        template_match = (
            None
            if unmatched_residues is None
            else not bool(component_residues & unmatched_residues)
        )
        net_charge, known_charge_sum, charge_complete = _formal_charge_summary(metadata)

        if contains_protein and not residue_names <= _STANDARD_PROTEIN_RESIDUES:
            classification = "covalent_candidate"
            reason = (
                "Standard protein residue(s) are covalently connected to "
                "nonstandard chemical residue(s); covalent ligand handling is "
                "detect-only."
            )
        elif contains_protein:
            classification = "standard"
            reason = "Contains only standard protein residue(s)."
        elif only_standard_solvent_or_ions:
            classification = "standard"
            reason = "Contains only standard water and/or ion residue(s)."
        elif template_match:
            classification = "standard"
            reason = "All residues matched the selected standard force field."
        else:
            classification = "ligand_candidate"
            reason = (
                "Disconnected non-polymer component with no selected standard "
                "force-field template match; eligible for automatic GAFF."
            )
        ligand_id = None
        if classification == "ligand_candidate":
            ligand_id = f"ligand-{ligand_number}"
            ligand_number += 1
        components.append(
            Component(
                index=component_index,
                ligand_id=ligand_id,
                topology_atom_indices=atom_indices,
                original_atom_indices=tuple(
                    atom.original_index for atom in metadata
                ),
                residues=residue_details,
                classification=classification,
                reason=reason,
                standard_template_match=template_match,
                net_known_formal_charge=net_charge,
                known_formal_charge_sum=known_charge_sum,
                formal_charge_complete=charge_complete,
                atom_formal_charges=tuple(
                    (atom.original_index, _numeric_value(atom.formal_charge))
                    for atom in metadata
                ),
            )
        )
    return tuple(components)


def write_components(
    path: Path,
    components: tuple[Component, ...],
    *,
    proteinff: str,
    waterff: str,
    ligand_mode: str,
    ligandff: str,
    production_topology_indices: dict[str, tuple[int, ...]] | None = None,
) -> None:
    """Write component classifications and parameterization provenance."""
    report = {
        "schema_version": 3,
        "covalent_ligand_scope": COVALENT_SCOPE,
        "proteinff": proteinff,
        "waterff": waterff,
        "ligand_mode": ligand_mode,
        "ligandff": ligandff,
        "components": [
            component.as_dict(
                proteinff,
                waterff,
                ligand_mode,
                ligandff,
                production_topology_indices,
            )
            for component in components
        ],
    }
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def reject_covalent_candidates(components: tuple[Component, ...]) -> None:
    """Reject explicitly detected protein--ligand covalent connectivity."""
    candidates = [
        component
        for component in components
        if component.classification == "covalent_candidate"
    ]
    if not candidates:
        return
    component_ids = ", ".join(str(component.index) for component in candidates)
    raise ValueError(
        "ligand_mode=auto detected covalent protein/nonstandard components "
        f"({component_ids}). Covalent ligand handling is detect-only and is not "
        "implemented; use a non-covalent input or parameterize the covalent "
        "system externally."
    )


def ligand_candidates(components: tuple[Component, ...]) -> tuple[Component, ...]:
    """Return components eligible for automatic small-molecule parameterization."""
    return tuple(
        component
        for component in components
        if component.classification == "ligand_candidate"
    )


def is_monitorable(component: Component) -> bool:
    """Return whether a component can be an early-stop target."""
    return component.classification != "covalent_candidate" and not (
        _is_solvent_or_ion_component(component)
    )


def component_summary(components: tuple[Component, ...]) -> str:
    """Render the durable IDs and selectors a user needs to pick a target.

    Solvent and ion components are counted rather than listed: an input with
    crystallographic water has hundreds of them and none can be a target.
    """
    chain_counts: dict[str, int] = {}
    for component in components:
        for chain in _component_chains(component):
            chain_counts[chain] = chain_counts.get(chain, 0) + 1

    rows = []
    skipped = 0
    for component in components:
        if not is_monitorable(component) and _is_solvent_or_ion_component(component):
            skipped += 1
            continue
        rows.append(
            (
                f"component-{component.index}",
                component.ligand_id or "-",
                component.classification,
                ", ".join(_component_chains(component)) or "-",
                str(len(component.residues)),
                _residue_composition(component),
                _selector_hint(component, chain_counts),
            )
        )

    headers = ("ID", "LIGAND ID", "CLASSIFICATION", "CHAINS", "RESIDUES", "COMPOSITION")
    widths = [
        max(len(headers[column]), *(len(row[column]) for row in rows))
        if rows
        else len(headers[column])
        for column in range(len(headers))
    ]
    lines = ["  ".join(header.ljust(widths[i]) for i, header in enumerate(headers)).rstrip()]
    for row in rows:
        lines.append(
            "  ".join(value.ljust(widths[i]) for i, value in enumerate(row[:6])).rstrip()
        )
        lines.append(f"    {row[6]}")
    if skipped:
        lines.append(
            f"\n{skipped} water/ion component(s) omitted; they cannot be early-stop targets."
        )
    return "\n".join(lines)


def _residue_composition(component: Component, limit: int = 4) -> str:
    names = sorted({str(residue["residue_name"]) for residue in component.residues})
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + f", +{len(names) - limit} more"


def _component_chains(component: Component) -> list[str]:
    return sorted({str(residue["chain_id"]) for residue in component.residues})


def _selector_hint(component: Component, chain_counts: dict[str, int]) -> str:
    """Suggest the selector that unambiguously names this component."""
    if not is_monitorable(component):
        return "not monitorable: " + component.reason
    if component.ligand_id is not None:
        return f"--early-stop --monitor-ligand {component.ligand_id}"
    unique_chains = [
        chain for chain in _component_chains(component) if chain_counts.get(chain) == 1
    ]
    if unique_chains:
        return f"--early-stop --monitor-chain {unique_chains[0]}"
    return f"--early-stop --monitor-component component-{component.index}"


def _is_solvent_or_ion_component(component: Component) -> bool:
    """Return whether every residue in a component is water or an ion."""
    return bool(component.residues) and all(
        str(residue["residue_name"]).upper() in _STANDARD_SOLVENT_AND_IONS
        for residue in component.residues
    )


def production_component_atom_indices(
    topology: app.Topology,
    components: tuple[Component, ...],
    presolvation_topology: app.Topology,
    source_order: tuple[int, ...],
) -> dict[str, tuple[int, ...]]:
    """Map every input component onto the final solvated topology.

    ``source_order`` lists input topology atom indices in the order the
    pre-solvation modeller holds them, which differs from the input order only
    when GAFF candidates were moved into dedicated ligand residues.
    """
    production_indices = _production_atom_map(topology, presolvation_topology)
    if len(production_indices) != len(source_order):
        raise ValueError(
            "The pre-solvation topology and its atom order disagree; cannot map "
            "components onto the production topology."
        )
    by_input_index = dict(zip(source_order, production_indices, strict=True))
    return {
        f"component-{component.index}": tuple(
            by_input_index[atom_index]
            for atom_index in component.topology_atom_indices
        )
        for component in components
    }


def _production_atom_map(
    topology: app.Topology, presolvation_topology: app.Topology
) -> list[int]:
    """Return the production index of every pre-solvation atom, in order.

    ``addSolvent`` only appends residues and ``addExtraParticles`` only inserts
    virtual sites, so the production topology opens with the pre-solvation atoms
    in their original order.  Walking both at once maps them in one pass and
    fails loudly if that assumption ever stops holding.
    """
    production_atoms = topology.atoms()
    mapped: list[int] = []
    for source in presolvation_topology.atoms():
        for atom in production_atoms:
            if (
                atom.name == source.name
                and atom.element is source.element
                and atom.residue.name == source.residue.name
            ):
                mapped.append(atom.index)
                break
            if atom.element is not None:
                raise ValueError(
                    "Solvation reordered the input atoms; cannot map components "
                    f"onto the production topology at atom {atom.index} "
                    f"({atom.residue.name} {atom.name})."
                )
        else:
            raise ValueError(
                f"Input atom {source.name} of residue {source.residue.name} is "
                "missing from the production topology."
            )
    return mapped


def create_ligand_modeller(
    structure: StructureData, candidates: tuple[Component, ...]
) -> tuple[app.Modeller, list[object], tuple[int, ...]]:
    """Create OpenFF molecules, a one-residue-per-ligand topology, and its order."""
    for component in candidates:
        _validate_ligand_component_chemistry(structure, component)
    chemistry = _require_ligand_dependencies()
    molecules = [
        _component_to_openff_molecule(structure, component, chemistry)
        for component in candidates
    ]
    modeller, source_order = _modeller_with_ligand_residues(structure, candidates)
    return modeller, molecules, source_order


def create_gaff_forcefield(
    protein_xml: str,
    water_xml: str,
    molecules: list[object],
    ligandff: str,
):
    """Create a SystemGenerator-backed force field for Amber plus GAFF."""
    chemistry = _require_ligand_dependencies()
    system_generator = chemistry["SystemGenerator"](
        forcefields=[protein_xml, water_xml],
        small_molecule_forcefield=ligandff,
        molecules=molecules,
    )
    return system_generator.forcefield


def _residue_details(metadata) -> tuple[dict[str, object], ...]:
    grouped: dict[tuple[str, str, str], list[int]] = {}
    for atom in metadata:
        key = (atom.chain_id, atom.residue_name, atom.residue_id)
        grouped.setdefault(key, []).append(atom.original_index)
    return tuple(
        {
            "chain_id": chain_id,
            "residue_name": residue_name,
            "residue_id": residue_id,
            "original_atom_indices": atom_indices,
        }
        for (chain_id, residue_name, residue_id), atom_indices in grouped.items()
    )


def _formal_charge_summary(metadata) -> tuple[int | float | None, int | float, bool]:
    charges = [_numeric_value(atom.formal_charge) for atom in metadata]
    known_charges = [charge for charge in charges if charge is not None]
    known_charge_sum = sum(known_charges)
    complete = len(known_charges) == len(charges)
    return (
        _json_number(known_charge_sum) if complete else None,
        _json_number(known_charge_sum),
        complete,
    )


def _numeric_value(value: object | None) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, Real):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def _json_number(value: float) -> int | float:
    """Narrow a charge sum to int when exact.

    ``sum`` of no known charges returns ``int``, and ``int.is_integer`` only
    exists on Python 3.12+, so normalize before asking.
    """
    number = float(value)
    return int(number) if number.is_integer() else number


def _require_ligand_dependencies() -> dict[str, object]:
    """Import optional ligand packages only when automatic ligands are needed."""
    try:
        from openff.toolkit.topology import Molecule
        from openmmforcefields.generators import SystemGenerator
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ImportError as error:
        raise RuntimeError(
            "Automatic ligand parameterization requires RDKit, the OpenFF "
            "Toolkit, and openmmforcefields. The OpenFF packages are not on "
            "PyPI, so pip alone cannot supply them. Create the documented "
            "environment with: conda env create -f environment.yml"
        ) from error
    return {
        "Chem": Chem,
        "Molecule": Molecule,
        "Point3D": Point3D,
        "SystemGenerator": SystemGenerator,
    }


def _component_to_openff_molecule(
    structure: StructureData,
    component: Component,
    chemistry: dict[str, object],
):
    """Build a sanitized, 3D-stereo OpenFF molecule from stored input chemistry."""
    chem = chemistry["Chem"]
    editable_molecule = chem.RWMol()
    metadata_by_index = {
        atom.topology_index: atom for atom in structure.atom_metadata
    }
    rdkit_indices = {}
    for topology_index in component.topology_atom_indices:
        atom = metadata_by_index[topology_index]
        if atom.atomic_number is None or atom.atomic_number < 1:
            raise ValueError(
                f"Ligand component {component.index} atom {atom.original_index} "
                "has no valid element."
            )
        charge = _formal_charge_for_ligand(
            atom.formal_charge, component, atom.original_index
        )
        rdkit_atom = chem.Atom(atom.atomic_number)
        rdkit_atom.SetFormalCharge(charge)
        rdkit_indices[topology_index] = editable_molecule.AddAtom(rdkit_atom)

    for bond in structure.bond_metadata:
        if (
            bond.atom1_topology_index in rdkit_indices
            and bond.atom2_topology_index in rdkit_indices
        ):
            bond_type, aromatic = _rdkit_bond_type(bond, component, chemistry)
            atom1 = rdkit_indices[bond.atom1_topology_index]
            atom2 = rdkit_indices[bond.atom2_topology_index]
            editable_molecule.AddBond(atom1, atom2, bond_type)
            if aromatic:
                editable_molecule.GetAtomWithIdx(atom1).SetIsAromatic(True)
                editable_molecule.GetAtomWithIdx(atom2).SetIsAromatic(True)
                editable_molecule.GetBondBetweenAtoms(atom1, atom2).SetIsAromatic(True)

    molecule = editable_molecule.GetMol()
    conformer = chem.Conformer(len(component.topology_atom_indices))
    conformer.Set3D(True)
    for topology_index, rdkit_index in rdkit_indices.items():
        coordinates = structure.positions[topology_index].value_in_unit(unit.angstrom)
        conformer.SetAtomPosition(
            rdkit_index,
            chemistry["Point3D"](coordinates.x, coordinates.y, coordinates.z),
        )
    molecule.AddConformer(conformer, assignId=True)
    try:
        chem.SanitizeMol(molecule)
        chem.AssignStereochemistryFrom3D(molecule)
        chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        openff_molecule = chemistry["Molecule"].from_rdkit(
            molecule, allow_undefined_stereo=True
        )
    except Exception as error:
        raise ValueError(
            f"Ligand component {component.index} could not be sanitized and "
            f"converted to an OpenFF molecule: {error}"
        ) from error
    openff_molecule.name = f"LIG{component.index + 1}"
    return openff_molecule


def _formal_charge_for_ligand(
    value: object | None, component: Component, original_index: int
) -> int:
    charge = _numeric_value(value)
    if charge is None:
        raise ValueError(
            f"Ligand component {component.index} atom {original_index} is missing "
            "a formal charge. ligand_mode=auto requires formal_charge data."
        )
    if not charge.is_integer():
        raise ValueError(
            f"Ligand component {component.index} atom {original_index} has "
            f"non-integral formal charge {value!r}."
        )
    return int(charge)


def _rdkit_bond_type(
    bond: BondMetadata, component: Component, chemistry: dict[str, object]
):
    order = _validated_bond_order(bond, component)
    if order == "aromatic":
        return chemistry["Chem"].BondType.AROMATIC, True
    bond_types = {
        1.0: chemistry["Chem"].BondType.SINGLE,
        1.5: chemistry["Chem"].BondType.AROMATIC,
        2.0: chemistry["Chem"].BondType.DOUBLE,
        3.0: chemistry["Chem"].BondType.TRIPLE,
    }
    return bond_types[order], order == 1.5


def _validate_ligand_component_chemistry(
    structure: StructureData, component: Component
) -> None:
    metadata_by_index = {
        atom.topology_index: atom for atom in structure.atom_metadata
    }
    for topology_index in component.topology_atom_indices:
        atom = metadata_by_index[topology_index]
        if atom.atomic_number is None or atom.atomic_number < 1:
            raise ValueError(
                f"Ligand component {component.index} atom {atom.original_index} "
                "has no valid element."
            )
        _formal_charge_for_ligand(
            atom.formal_charge, component, atom.original_index
        )
    component_atom_indices = set(component.topology_atom_indices)
    for bond in structure.bond_metadata:
        if (
            bond.atom1_topology_index in component_atom_indices
            and bond.atom2_topology_index in component_atom_indices
        ):
            _validated_bond_order(bond, component)


def _validated_bond_order(
    bond: BondMetadata, component: Component
) -> float | str:
    if bond.order is None:
        raise ValueError(
            f"Ligand component {component.index} bond "
            f"{bond.atom1_original_index}-{bond.atom2_original_index} is missing "
            "a bond order. ligand_mode=auto requires bond order data."
        )
    order_text = str(bond.order).strip().lower()
    if order_text in {"dative", "coordinate", "zero", "0", "0.0"}:
        raise ValueError(
            f"Ligand component {component.index} bond "
            f"{bond.atom1_original_index}-{bond.atom2_original_index} has "
            f"unsupported zero/dative bond order {bond.order!r}; no bond order "
            "will be guessed."
        )
    if order_text in {"ar", "aromatic"}:
        return "aromatic"
    order = _numeric_value(bond.order)
    if order not in {1.0, 1.5, 2.0, 3.0}:
        raise ValueError(
            f"Ligand component {component.index} bond "
            f"{bond.atom1_original_index}-{bond.atom2_original_index} has "
            f"unsupported bond order {bond.order!r}; no bond order will be guessed."
        )
    return order


def _modeller_with_ligand_residues(
    structure: StructureData, candidates: tuple[Component, ...]
) -> tuple[app.Modeller, tuple[int, ...]]:
    """Copy the input topology while normalizing each ligand to one residue."""
    candidate_by_atom = {
        atom_index: component
        for component in candidates
        for atom_index in component.topology_atom_indices
    }
    source_atoms = list(structure.topology.atoms())
    topology = app.Topology()
    box_vectors = structure.topology.getPeriodicBoxVectors()
    if box_vectors is not None:
        topology.setPeriodicBoxVectors(box_vectors)
    atom_map = {}
    position_order = []

    for source_chain in structure.topology.chains():
        target_chain = None
        for source_residue in source_chain.residues():
            retained_atoms = [
                atom
                for atom in source_residue.atoms()
                if atom.index not in candidate_by_atom
            ]
            if not retained_atoms:
                continue
            if target_chain is None:
                target_chain = topology.addChain(source_chain.id)
            target_residue = topology.addResidue(
                source_residue.name,
                target_chain,
                id=source_residue.id,
                insertionCode=source_residue.insertionCode,
            )
            for source_atom in retained_atoms:
                atom_map[source_atom.index] = topology.addAtom(
                    source_atom.name,
                    source_atom.element,
                    target_residue,
                    id=source_atom.id,
                )
                position_order.append(source_atom.index)

    for component in candidates:
        ligand_chain = topology.addChain(f"L{component.index + 1}")
        ligand_residue = topology.addResidue(
            f"LIG{component.index + 1}", ligand_chain, id=str(component.index + 1)
        )
        for source_atom_index in component.topology_atom_indices:
            source_atom = source_atoms[source_atom_index]
            atom_map[source_atom_index] = topology.addAtom(
                source_atom.name,
                source_atom.element,
                ligand_residue,
                id=source_atom.id,
            )
            position_order.append(source_atom_index)

    for source_atom1, source_atom2 in structure.topology.bonds():
        topology.addBond(atom_map[source_atom1.index], atom_map[source_atom2.index])
    positions = _reordered_positions(structure.positions, position_order)
    return app.Modeller(topology, positions), tuple(position_order)


def _reordered_positions(positions, order: list[int]):
    position_unit = getattr(positions, "unit", None)
    if position_unit is None:
        return [positions[index] for index in order]
    return unit.Quantity(
        [positions[index].value_in_unit(position_unit) for index in order],
        position_unit,
    )
