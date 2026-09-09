"""Round-trip tests for the MAE structures written beside each stage's PDB."""

from __future__ import annotations

from pathlib import Path

import openmm as mm
import pytest
from openmm import app, unit
from openmm.app.element import Element

from ommflow.lib.readers import MAEReader
from ommflow.lib.writers import write_mae


def _ligand_and_water_topology() -> tuple[app.Topology, object]:
    """A two-residue topology: a ligand with real orders, and one water."""
    topology = app.Topology()
    ligand_chain = topology.addChain("L1")
    ligand = topology.addResidue("LIG", ligand_chain, id="1")
    carbon = topology.addAtom("C1", Element.getBySymbol("C"), ligand)
    oxygen = topology.addAtom("O1", Element.getBySymbol("O"), ligand)
    nitrogen = topology.addAtom("N1", Element.getBySymbol("N"), ligand)
    topology.addBond(carbon, oxygen, order=2)
    topology.addBond(carbon, nitrogen, order=3)

    water_chain = topology.addChain("2")
    water = topology.addResidue("HOH", water_chain, id="2")
    water_o = topology.addAtom("O", Element.getBySymbol("O"), water)
    water_h1 = topology.addAtom("H1", Element.getBySymbol("H"), water)
    water_h2 = topology.addAtom("H2", Element.getBySymbol("H"), water)
    topology.addBond(water_o, water_h1)
    topology.addBond(water_o, water_h2)

    positions = unit.Quantity(
        [
            mm.Vec3(0.1, 0.2, 0.3),
            mm.Vec3(0.22, 0.2, 0.3),
            mm.Vec3(0.34, 0.2, 0.3),
            mm.Vec3(0.5, 0.5, 0.5),
            mm.Vec3(0.59, 0.5, 0.5),
            mm.Vec3(0.47, 0.59, 0.5),
        ],
        unit.nanometer,
    )
    return topology, positions


def test_written_mae_reads_back_as_the_same_structure(tmp_path: Path) -> None:
    """The file has to survive our own reader to be worth writing."""
    topology, positions = _ligand_and_water_topology()
    path = tmp_path / "solvated.mae"

    write_mae(path, topology, positions, title="solvated")
    structure = MAEReader(path)

    atoms = list(structure.topology.atoms())
    assert [atom.name for atom in atoms] == ["C1", "O1", "N1", "O", "H1", "H2"]
    assert [atom.element.symbol for atom in atoms] == ["C", "O", "N", "O", "H", "H"]
    assert [atom.residue.name for atom in atoms] == ["LIG"] * 3 + ["HOH"] * 3
    assert [atom.residue.chain.id for atom in atoms] == ["L1"] * 3 + ["2"] * 3
    assert structure.topology.getNumBonds() == topology.getNumBonds()

    original = positions.value_in_unit(unit.angstrom)
    written = structure.positions.value_in_unit(unit.angstrom)
    for source, target in zip(original, written, strict=True):
        # Coordinates are written to six decimals of an Angstrom, far finer
        # than the single-precision coordinates a DCD frame carries.
        assert list(target) == pytest.approx(list(source), abs=1e-5)


def test_bond_orders_are_kept_and_default_to_single(tmp_path: Path) -> None:
    """Input chemistry survives; what the solvent builder adds is single."""
    topology, positions = _ligand_and_water_topology()
    path = tmp_path / "solvated.mae"

    write_mae(path, topology, positions)

    assert [bond.order for bond in MAEReader(path).topology.bonds()] == [2, 3, 1, 1]


def test_input_mae_orders_survive_a_modeller_copy(tmp_path: Path) -> None:
    """Solvation copies the topology, so orders have to ride on its bonds."""
    source = tmp_path / "input.mae"
    source.write_text(
        """f_m_ct {
 s_m_title
:::
 lig
 m_atom[3] {
 i_m_residue_number
 s_m_chain_name
 s_m_pdb_residue_name
 s_m_pdb_atom_name
 i_m_atomic_number
 r_m_x_coord
 r_m_y_coord
 r_m_z_coord
:::
 1 1 A "LIG " " C1 " 6 1.0 2.0 3.0
 2 1 A "LIG " " O1 " 8 2.2 2.0 3.0
 3 1 A "LIG " " N1 " 7 3.4 2.0 3.0
:::
 }
 m_bond[2] {
 i_m_from
 i_m_to
 i_m_order
:::
 1 1 2 2
 2 1 3 3
:::
 }
}
""",
        encoding="utf-8",
    )
    structure = MAEReader(source)
    assert [bond.order for bond in structure.topology.bonds()] == [2, 3]

    modeller = app.Modeller(structure.topology, structure.positions)
    written = tmp_path / "solvated.mae"
    write_mae(written, modeller.topology, modeller.positions)

    assert [bond.order for bond in MAEReader(written).topology.bonds()] == [2, 3]


def test_periodic_box_is_written_in_angstroms(tmp_path: Path) -> None:
    """The chorus box is what gives a viewer the periodic cell."""
    topology, positions = _ligand_and_water_topology()
    path = tmp_path / "equilibrated.mae"

    write_mae(
        path,
        topology,
        positions,
        box_vectors=unit.Quantity(
            [mm.Vec3(3.0, 0, 0), mm.Vec3(0, 3.0, 0), mm.Vec3(0, 0, 3.0)],
            unit.nanometer,
        ),
    )

    contents = path.read_text(encoding="utf-8")
    assert "r_chorus_box_ax" in contents
    assert "30.000000" in contents


def test_a_position_count_mismatch_is_rejected(tmp_path: Path) -> None:
    """A silently truncated structure would be worse than no file at all."""
    topology, positions = _ligand_and_water_topology()
    path = tmp_path / "solvated.mae"

    try:
        write_mae(path, topology, positions[:4])
    except ValueError as error:
        assert "6 topology atoms but 4 positions" in str(error)
    else:
        raise AssertionError("A short position list must not be written.")
