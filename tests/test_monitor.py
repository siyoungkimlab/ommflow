"""Focused tests for ligand-detachment monitor behavior without running MD."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from openmm import app, unit
import openmm as mm
from openmm.app.element import Element

from ommflow.lib.ligands import (
    _modeller_with_ligand_residues,
    _production_atom_map,
    classify_components,
    production_component_atom_indices,
    write_components,
)
from ommflow.lib.monitor import (
    MONITOR_FIELDS,
    append_monitor_row,
    detachment_update,
    initialize_monitor,
    load_status,
    monitor_measurement,
    restore_consecutive_count,
    select_monitor_ligand,
    select_monitor_target,
    write_status,
)
from ommflow.lib.readers import AtomMetadata, BondMetadata, StructureData


def _components(*ligand_ids: str) -> dict[str, object]:
    return {
        "components": [
            {
                "classification": "ligand_candidate",
                "ligand_id": ligand_id,
                "production_topology_atom_indices": [index * 2, index * 2 + 1],
            }
            for index, ligand_id in enumerate(ligand_ids)
        ]
    }


def _standard_components() -> dict[str, object]:
    return {
        "components": [
            {
                "component_id": "component-0",
                "component_index": 0,
                "classification": "standard",
                "monitorable": True,
                "chains": ["A"],
                "residues": [{"chain_id": "A", "residue_name": "ALA"}],
                "production_topology_atom_indices": [0],
            },
            {
                "component_id": "component-1",
                "component_index": 1,
                "classification": "standard",
                "monitorable": True,
                "chains": ["B"],
                "residues": [{"chain_id": "B", "residue_name": "ALA"}],
                "production_topology_atom_indices": [1],
            },
        ]
    }


def _target_args(**overrides):
    values = {
        "monitor_ligand": None,
        "monitor_chain": None,
        "monitor_component": None,
        "monitor_interval_ns": 0.1,
        "pocket_cutoff_nm": 0.15,
        "contact_cutoff_nm": 0.11,
        "detach_cutoff_nm": 0.8,
        "confirmation_checks": 2,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_monitor_ligand_selection_requires_explicit_choice_when_ambiguous() -> None:
    assert select_monitor_ligand(_components("ligand-0"), None)["ligand_id"] == "ligand-0"
    assert (
        select_monitor_ligand(_components("ligand-0", "ligand-1"), "ligand-1")[
            "ligand_id"
        ]
        == "ligand-1"
    )
    with pytest.raises(ValueError, match="ambiguous.*ligand-0, ligand-1"):
        select_monitor_ligand(_components("ligand-0", "ligand-1"), None)
    with pytest.raises(ValueError, match="Unknown monitor ligand"):
        select_monitor_ligand(_components(), "ligand-0")


def test_chain_and_component_selectors_resolve_standard_peptide_component() -> None:
    components = _standard_components()
    chain_target = select_monitor_target(components, _target_args(monitor_chain="B"))
    assert chain_target["_monitor_component_id"] == "component-1"
    component_target = select_monitor_target(
        components, _target_args(monitor_component="component-1")
    )
    assert component_target["_monitor_component_id"] == "component-1"

    ambiguous = _standard_components()
    ambiguous["components"][1]["chains"] = ["B"]
    ambiguous["components"][0]["chains"] = ["B"]
    with pytest.raises(ValueError, match="multiple monitorable.*component-0.*component-1"):
        select_monitor_target(ambiguous, _target_args(monitor_chain="B"))
    with pytest.raises(ValueError, match="No monitorable.*'Z'.*component-0"):
        select_monitor_target(components, _target_args(monitor_chain="Z"))
    with pytest.raises(ValueError, match="Unknown or non-monitorable component"):
        select_monitor_target(components, _target_args(monitor_component="component-9"))


def test_measurements_use_minimum_image_heavy_atom_distances() -> None:
    minimum_distance, contacts = monitor_measurement(
        (0,),
        (1,),
        [(0.05, 0.0, 0.0), (0.95, 0.0, 0.0)],
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        0.11,
    )
    assert minimum_distance == pytest.approx(0.1)
    assert contacts == 1

    count, detached = detachment_update(0.9, 0, 0.8, 0, 2)
    assert (count, detached) == (1, False)
    count, detached = detachment_update(0.9, 0, 0.8, count, 2)
    assert (count, detached) == (2, True)
    assert detachment_update(0.9, 1, 0.8, count, 2) == (0, False)


def test_initial_pocket_selection_uses_pbc_and_excludes_ligand_atoms() -> None:
    topology = app.Topology()
    chain = topology.addChain("A")
    protein_atom = topology.addAtom(
        "CA", Element.getBySymbol("C"), topology.addResidue("ALA", chain)
    )
    ligand_atom = topology.addAtom(
        "C1", Element.getBySymbol("C"), topology.addResidue("LIG", chain)
    )
    args = SimpleNamespace(
        monitor_ligand=None,
        monitor_interval_ns=0.1,
        pocket_cutoff_nm=0.15,
        contact_cutoff_nm=0.11,
        detach_cutoff_nm=0.8,
        confirmation_checks=2,
    )
    definition = initialize_monitor(
        topology,
        [(0.95, 0.0, 0.0), (0.05, 0.0, 0.0)],
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        {
            "components": [
                {
                    "classification": "ligand_candidate",
                    "ligand_id": "ligand-0",
                    "production_topology_atom_indices": [ligand_atom.index],
                }
            ]
        },
        args,
    )
    assert definition.ligand_atom_indices == (ligand_atom.index,)
    assert definition.pocket_atom_indices == (protein_atom.index,)
    assert definition.initial_contact_count == 1


def test_protein_chain_target_uses_the_other_component_as_its_pocket() -> None:
    topology = app.Topology()
    receptor = topology.addChain("A")
    peptide = topology.addChain("B")
    receptor_atom = topology.addAtom(
        "CA", Element.getBySymbol("C"), topology.addResidue("ALA", receptor)
    )
    peptide_atom = topology.addAtom(
        "CA", Element.getBySymbol("C"), topology.addResidue("GLY", peptide)
    )
    definition = initialize_monitor(
        topology,
        [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0)],
        None,
        _standard_components(),
        _target_args(monitor_chain="B"),
    )
    assert definition.component_id == "component-1"
    assert definition.ligand_atom_indices == (peptide_atom.index,)
    assert definition.pocket_atom_indices == (receptor_atom.index,)


def test_component_ligand_ids_and_production_indices_are_durable() -> None:
    topology = app.Topology()
    chain = topology.addChain("A")
    atoms = [
        topology.addAtom(
            "CA", Element.getBySymbol("C"), topology.addResidue("ALA", chain)
        ),
        topology.addAtom(
            "C1", Element.getBySymbol("C"), topology.addResidue("LIG", chain)
        ),
        topology.addAtom(
            "C2", Element.getBySymbol("C"), topology.addResidue("DRG", chain)
        ),
    ]
    structure = StructureData(
        topology,
        unit.Quantity([mm.Vec3(0, 0, 0)] * 3, unit.nanometer),
        tuple(
            AtomMetadata(
                topology_index=atom.index,
                original_index=100 + atom.index,
                atomic_number=6,
                formal_charge=0,
                name=atom.name,
                residue_name=atom.residue.name,
                residue_id=atom.residue.id,
                chain_id=atom.residue.chain.id,
            )
            for atom in atoms
        ),
        (),
    )
    components = classify_components(structure)
    assert [component.ligand_id for component in components] == [
        None,
        "ligand-0",
        "ligand-1",
    ]
    assert [component.index for component in components] == [0, 1, 2]
    mapped = production_component_atom_indices(
        topology, components, topology, tuple(range(topology.getNumAtoms()))
    )
    assert mapped["component-0"] == (0,)
    assert mapped["component-1"] == (1,)
    assert mapped["component-2"] == (2,)

    artifact_dir = Path(__file__).parent / ".monitor-test-artifacts"
    artifact_dir.mkdir(exist_ok=True)
    output = artifact_dir / "components.json"
    try:
        write_components(
            output,
            components,
            proteinff="amber19sb",
            waterff="opc",
            ligand_mode="auto",
            ligandff="gaff-2.11",
            production_topology_indices={"ligand-0": (8,), "ligand-1": (9,)},
        )
        report = json.loads(output.read_text(encoding="utf-8"))
        assert report["schema_version"] == 3
        assert report["components"][1]["component_id"] == "component-1"
        assert report["components"][1]["monitorable"] is True
        assert report["components"][1]["ligand_id"] == "ligand-0"
        assert report["components"][1]["original_atom_indices"] == [101]
        assert report["components"][1]["atom_formal_charges"] == [
            {"original_atom_index": 101, "formal_charge": 0.0}
        ]
        assert report["components"][1]["production_topology_atom_indices"] == [8]
    finally:
        output.unlink(missing_ok=True)
        artifact_dir.rmdir()


def test_monitor_csv_and_status_append_without_duplicate_header() -> None:
    artifact_dir = Path(__file__).parent / ".monitor-test-artifacts"
    artifact_dir.mkdir(exist_ok=True)
    monitor_csv = artifact_dir / "monitor.csv"
    status_json = artifact_dir / "status.json"
    for path in (monitor_csv, status_json):
        path.unlink(missing_ok=True)
    try:
        row = {
            "production_time_ns": "0.1",
            "step": 50000,
            "min_ligand_pocket_distance_nm": "0.9",
            "contact_count": 0,
            "initial_contact_count": 3,
            "contact_fraction": "0",
            "consecutive_detached_count": 1,
            "detached": "false",
        }
        assert append_monitor_row(monitor_csv, row)
        assert not append_monitor_row(monitor_csv, row)
        row["step"] = 100000
        assert append_monitor_row(monitor_csv, row)
        lines = monitor_csv.read_text(encoding="utf-8").splitlines()
        assert lines[0].split(",") == list(MONITOR_FIELDS)
        assert len(lines) == 3

        status = {
            "outcome": "running",
            "final_production_step": 100000,
            "consecutive_detached_count": 1,
        }
        write_status(status_json, status)
        assert load_status(status_json) == status
        assert restore_consecutive_count(status, monitor_csv, 100000) == 1
    finally:
        for path in (monitor_csv, status_json):
            path.unlink(missing_ok=True)
        artifact_dir.rmdir()


def _atom(topology, chain, residue_name, residue_id, name, symbol="C"):
    element = None if symbol is None else Element.getBySymbol(symbol)
    residue = topology.addResidue(residue_name, chain, id=residue_id)
    return topology.addAtom(name, element, residue)


def test_components_map_positionally_through_duplicate_residue_labels() -> None:
    """Two molecules can share a chain, residue name and residue ID."""
    topology = app.Topology()
    chain = topology.addChain("")
    atoms = [_atom(topology, chain, "LIG", "1", f"C{index}") for index in range(4)]
    topology.addBond(atoms[0], atoms[1])
    topology.addBond(atoms[2], atoms[3])
    structure = StructureData(
        topology,
        unit.Quantity([mm.Vec3(0, 0, 0)] * 4, unit.nanometer),
        tuple(
            AtomMetadata(
                topology_index=atom.index,
                original_index=atom.index,
                atomic_number=6,
                formal_charge=0,
                name=atom.name,
                residue_name=atom.residue.name,
                residue_id=atom.residue.id,
                chain_id=atom.residue.chain.id,
            )
            for atom in atoms
        ),
        tuple(
            BondMetadata(
                atom1_topology_index=first,
                atom2_topology_index=second,
                atom1_original_index=first,
                atom2_original_index=second,
                order=1,
            )
            for first, second in ((0, 1), (2, 3))
        ),
    )
    components = classify_components(structure)
    assert len(components) == 2

    mapped = production_component_atom_indices(
        topology, components, topology, tuple(range(4))
    )
    assert mapped["component-0"] == (0, 1)
    assert mapped["component-1"] == (2, 3)


def test_component_mapping_skips_inserted_virtual_sites() -> None:
    """addExtraParticles inserts massless sites among the input atoms."""
    presolvation = app.Topology()
    chain = presolvation.addChain("A")
    water = presolvation.addResidue("HOH", chain, id="1")
    for name, symbol in (("O", "O"), ("H1", "H"), ("H2", "H")):
        presolvation.addAtom(name, Element.getBySymbol(symbol), water)

    production = app.Topology()
    chain = production.addChain("A")
    water = production.addResidue("HOH", chain, id="1")
    for name, symbol in (("O", "O"), ("H1", "H"), ("H2", "H"), ("M", None)):
        production.addAtom(
            name, None if symbol is None else Element.getBySymbol(symbol), water
        )
    added = production.addResidue("HOH", production.addChain("2"), id="2")
    for name, symbol in (("O", "O"), ("H1", "H"), ("H2", "H"), ("M", None)):
        production.addAtom(
            name, None if symbol is None else Element.getBySymbol(symbol), added
        )

    mapped = _production_atom_map(production, presolvation)
    assert mapped == [0, 1, 2]


def test_gaff_reordering_keeps_every_component_addressable() -> None:
    """Ligand atoms move to their own residue; the rest keep their input order."""
    topology = app.Topology()
    chain = topology.addChain("A")
    layout = (
        ("ALA", "1", ("N", "CA")),
        ("LIG", "2", ("C1", "C2")),
        ("ALA", "3", ("C", "O")),
    )
    atoms = []
    for residue_name, residue_id, names in layout:
        residue = topology.addResidue(residue_name, chain, id=residue_id)
        atoms.extend(
            topology.addAtom(name, Element.getBySymbol(name[0]), residue)
            for name in names
        )
    bonds = ((0, 1), (1, 4), (4, 5), (2, 3))
    for first, second in bonds:
        topology.addBond(atoms[first], atoms[second])
    structure = StructureData(
        topology,
        unit.Quantity([mm.Vec3(index, 0, 0) for index in range(6)], unit.nanometer),
        tuple(
            AtomMetadata(
                topology_index=atom.index,
                original_index=atom.index,
                atomic_number=atom.element.atomic_number,
                formal_charge=0,
                name=atom.name,
                residue_name=atom.residue.name,
                residue_id=atom.residue.id,
                chain_id=atom.residue.chain.id,
            )
            for atom in atoms
        ),
        tuple(
            BondMetadata(
                atom1_topology_index=first,
                atom2_topology_index=second,
                atom1_original_index=first,
                atom2_original_index=second,
                order=1,
            )
            for first, second in bonds
        ),
    )
    components = classify_components(structure)
    candidates = tuple(
        component
        for component in components
        if component.classification == "ligand_candidate"
    )
    assert [component.ligand_id for component in candidates] == ["ligand-0"]

    modeller, source_order = _modeller_with_ligand_residues(structure, candidates)
    assert source_order == (0, 1, 4, 5, 2, 3)

    mapped = production_component_atom_indices(
        modeller.topology, components, modeller.topology, source_order
    )
    assert mapped["component-0"] == (0, 1, 2, 3)
    assert mapped["component-1"] == (4, 5)


def test_components_classify_without_any_formal_charge_data() -> None:
    """PDB and GRO inputs carry no formal charges at all."""
    topology = app.Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue("ALA", chain, id="1")
    atoms = [
        topology.addAtom(name, Element.getBySymbol("C"), residue)
        for name in ("N", "CA")
    ]
    topology.addBond(*atoms)
    structure = StructureData(
        topology,
        unit.Quantity([mm.Vec3(0, 0, 0)] * 2, unit.nanometer),
        tuple(
            AtomMetadata(
                topology_index=atom.index,
                original_index=atom.index,
                atomic_number=6,
                formal_charge=None,
                name=atom.name,
                residue_name=atom.residue.name,
                residue_id=atom.residue.id,
                chain_id=atom.residue.chain.id,
            )
            for atom in atoms
        ),
        (
            BondMetadata(
                atom1_topology_index=0,
                atom2_topology_index=1,
                atom1_original_index=0,
                atom2_original_index=1,
                order=None,
            ),
        ),
    )
    (component,) = classify_components(structure)
    assert component.formal_charge_complete is False
    assert component.net_known_formal_charge is None
    assert component.known_formal_charge_sum == 0


def _mae_block(rows: str, natoms: int) -> str:
    """A minimal MAE m_atom/m_bond pair in Maestro's real layout."""
    return f"""f_m_ct {{
 s_m_title
:::
 test
 m_atom[{natoms}] {{
 i_m_residue_number
 s_m_chain_name
 s_m_pdb_residue_name
 s_m_pdb_atom_name
 i_m_atomic_number
 i_m_formal_charge
 r_m_x_coord
 r_m_y_coord
 r_m_z_coord
:::
{rows}:::
 }}
 m_bond[1] {{
 i_m_from
 i_m_to
 i_m_order
:::
 1 1 2 1
:::
 }}
}}
"""


def test_mae_tables_ignore_the_closing_separator(tmp_path: Path) -> None:
    """A block is "headers ::: rows :::"; the trailing one is not data."""
    from ommflow.lib.readers import MAEReader

    rows = (
        ' 1 1 A "ALA " " N  " 7 0 1.0 2.0 3.0\n'
        ' 2 1 A "ALA " " CA " 6 0 4.0 5.0 6.0\n'
    )
    path = tmp_path / "two.mae"
    path.write_text(_mae_block(rows, 2), encoding="utf-8")

    structure = MAEReader(path)
    atoms = list(structure.topology.atoms())
    assert len(atoms) == 2
    # Padding is a PDB column artifact, not part of the name; leaving it in
    # would stop every force-field template from matching.
    assert [atom.name for atom in atoms] == ["N", "CA"]
    assert [atom.residue.name for atom in atoms] == ["ALA", "ALA"]
    assert atoms[0].residue.chain.id == "A"
    assert structure.topology.getNumBonds() == 1
