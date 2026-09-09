"""Readers for structural inputs and their chemistry metadata."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
import re
import sqlite3

import openmm as mm
from openmm import app, unit
from openmm.app.element import Element


@dataclass(frozen=True)
class AtomMetadata:
    """Chemical and source information for one topology atom."""

    topology_index: int
    original_index: int
    atomic_number: int | None
    formal_charge: object | None
    name: str
    residue_name: str
    residue_id: str
    chain_id: str


@dataclass(frozen=True)
class BondMetadata:
    """Chemical and source information for one topology bond."""

    atom1_topology_index: int
    atom2_topology_index: int
    atom1_original_index: int
    atom2_original_index: int
    order: object | None


@dataclass
class StructureData:
    """Topology, coordinates, and chemistry retained from a structural input."""

    topology: app.Topology
    positions: object
    atom_metadata: tuple[AtomMetadata, ...]
    bond_metadata: tuple[BondMetadata, ...]

    @property
    def atom_formal_charges(self) -> tuple[object | None, ...]:
        """Formal charges in topology atom order, or ``None`` when unavailable."""
        return tuple(atom.formal_charge for atom in self.atom_metadata)

    @property
    def bond_orders(self) -> tuple[object | None, ...]:
        """Bond orders in topology bond order, or ``None`` when unavailable."""
        return tuple(bond.order for bond in self.bond_metadata)


def _topology_bond_order(order: object | None) -> int | None:
    """Reduce a source bond order to the integer an OpenMM bond can carry.

    Kekule orders pass through unchanged, so a MAE or DMS input keeps the
    double and triple bonds it was written with. An aromatic order has no
    integer form and is recorded as single: this value is only ever drawn by a
    writer, while ligand parameterization reads the unreduced order from the
    bond metadata.
    """
    if order is None:
        return None
    try:
        number = float(order)
    except (TypeError, ValueError):
        return 1
    return int(number) if number.is_integer() else 1


def _structure_data(
    topology: app.Topology,
    positions,
    original_indices: list[int] | None = None,
    formal_charges: list[object | None] | None = None,
    bond_orders: list[object | None] | None = None,
) -> StructureData:
    """Create metadata for a topology whose atom and bond ordering is known."""
    atoms = list(topology.atoms())
    bonds = list(topology.bonds())
    if original_indices is None:
        original_indices = list(range(len(atoms)))
    if formal_charges is None:
        formal_charges = [None] * len(atoms)
    if bond_orders is None:
        bond_orders = [None] * len(bonds)
    if len(atoms) != len(original_indices) or len(atoms) != len(formal_charges):
        raise ValueError("Atom chemistry metadata does not match the topology.")
    if len(bonds) != len(bond_orders):
        raise ValueError("Bond chemistry metadata does not match the topology.")

    atom_metadata = tuple(
        AtomMetadata(
            topology_index=atom.index,
            original_index=original_index,
            atomic_number=atom.element.atomic_number if atom.element else None,
            formal_charge=formal_charge,
            name=atom.name,
            residue_name=atom.residue.name,
            residue_id=atom.residue.id,
            chain_id=atom.residue.chain.id,
        )
        for atom, original_index, formal_charge in zip(
            atoms, original_indices, formal_charges, strict=True
        )
    )
    metadata_by_topology_index = {
        atom.topology_index: atom for atom in atom_metadata
    }
    bond_metadata = tuple(
        BondMetadata(
            atom1_topology_index=bond[0].index,
            atom2_topology_index=bond[1].index,
            atom1_original_index=metadata_by_topology_index[
                bond[0].index
            ].original_index,
            atom2_original_index=metadata_by_topology_index[
                bond[1].index
            ].original_index,
            order=order,
        )
        for bond, order in zip(bonds, bond_orders, strict=True)
    )
    return StructureData(topology, positions, atom_metadata, bond_metadata)


class DMSReader(StructureData):
    """Read structural data and available chemistry from a Desmond DMS database."""

    def __init__(self, path: Path) -> None:
        topology = app.Topology()
        positions = []
        original_indices = []
        formal_charges = []
        bond_orders = []
        with sqlite3.connect(str(path)) as connection:
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            if not {"particle", "bond"} <= tables:
                raise ValueError(
                    "A structural DMS file must contain 'particle' and 'bond' tables."
                )
            particle_columns = {
                row[1] for row in connection.execute("PRAGMA table_info(particle)")
            }
            required_particle_columns = {
                "id",
                "name",
                "anum",
                "resname",
                "resid",
                "chain",
                "x",
                "y",
                "z",
            }
            missing_columns = required_particle_columns - particle_columns
            if missing_columns:
                raise ValueError(
                    "DMS particle table is missing required column(s): "
                    + ", ".join(sorted(missing_columns))
                )
            bond_columns = {
                row[1] for row in connection.execute("PRAGMA table_info(bond)")
            }
            if not {"p0", "p1"} <= bond_columns:
                raise ValueError(
                    "DMS bond table must contain p0 and p1 columns."
                )

            charge_select = (
                "formal_charge" if "formal_charge" in particle_columns else "NULL"
            )
            order_select = '"order"' if "order" in bond_columns else "NULL"
            atoms = {}
            previous_chain = object()
            previous_residue = None
            chain = None
            residue = None
            query = f"""
                SELECT id, name, anum, resname, resid, chain, x, y, z,
                       {charge_select}
                FROM particle ORDER BY id
            """
            for (
                atom_id,
                name,
                atomic_number,
                residue_name,
                residue_id,
                chain_id,
                x,
                y,
                z,
                formal_charge,
            ) in connection.execute(query):
                if chain_id != previous_chain:
                    chain = topology.addChain(str(chain_id or ""))
                    previous_chain = chain_id
                    previous_residue = None
                residue_key = (residue_id, residue_name)
                if residue_key != previous_residue:
                    residue = topology.addResidue(
                        str(residue_name), chain, id=str(residue_id)
                    )
                    previous_residue = residue_key
                atomic_number = int(atomic_number)
                element = (
                    None
                    if atomic_number == 0
                    else Element.getByAtomicNumber(atomic_number)
                )
                atoms[atom_id] = topology.addAtom(str(name), element, residue)
                positions.append(mm.Vec3(x, y, z))
                original_indices.append(int(atom_id))
                formal_charges.append(formal_charge)

            for atom1, atom2, order in connection.execute(
                f"SELECT p0, p1, {order_select} FROM bond"
            ):
                try:
                    topology.addBond(
                        atoms[atom1],
                        atoms[atom2],
                        order=_topology_bond_order(order),
                    )
                except KeyError as error:
                    raise ValueError(
                        f"DMS bond references unknown particle ID: {error.args[0]}"
                    ) from error
                bond_orders.append(order)
        data = _structure_data(
            topology,
            positions * unit.angstrom,
            original_indices=original_indices,
            formal_charges=formal_charges,
            bond_orders=bond_orders,
        )
        super().__init__(
            data.topology, data.positions, data.atom_metadata, data.bond_metadata
        )


class MAEReader(StructureData):
    """Read topology, coordinates, and available chemistry from a MAE file."""

    _token_pattern = re.compile(r'"(?:\\.|[^"\\])*"|[^\s{}]+')

    def __init__(self, path: Path) -> None:
        contents = path.read_text(encoding="utf-8")
        atom_columns, atom_rows = self._read_table(contents, "m_atom")
        bond_columns, bond_rows = self._read_table(contents, "m_bond")

        required_atom_columns = {
            "i_m_atomic_number",
            "r_m_x_coord",
            "r_m_y_coord",
            "r_m_z_coord",
        }
        missing_columns = required_atom_columns - set(atom_columns)
        if missing_columns:
            raise ValueError(
                "MAE atom table is missing required column(s): "
                + ", ".join(sorted(missing_columns))
            )
        if not {"i_m_from", "i_m_to"} <= set(bond_columns):
            raise ValueError("MAE bond table must contain i_m_from and i_m_to columns.")

        topology = app.Topology()
        positions = []
        original_indices = []
        formal_charges = []
        bond_orders = []
        atoms = {}
        previous_chain = object()
        previous_residue = None
        chain = None
        residue = None
        formal_charge_column = self._first_present(
            atom_columns, "i_m_formal_charge", "r_m_formal_charge"
        )
        bond_order_column = self._first_present(
            bond_columns, "i_m_order", "r_m_order"
        )
        for index, row in enumerate(atom_rows, start=1):
            values = dict(zip(atom_columns, row, strict=True))
            residue_name = str(values.get("s_m_pdb_residue_name", "UNK"))
            residue_id = values.get("i_m_residue_number", str(index))
            chain_id = str(values.get("s_m_chain_name", ""))
            if chain_id != previous_chain:
                chain = topology.addChain(chain_id)
                previous_chain = chain_id
                previous_residue = None
            residue_key = (residue_id, residue_name)
            if residue_key != previous_residue:
                residue = topology.addResidue(residue_name, chain, id=str(residue_id))
                previous_residue = residue_key
            atom_name = str(values.get("s_m_pdb_atom_name", f"X{index}"))
            atomic_number = int(values["i_m_atomic_number"])
            atoms[index] = topology.addAtom(
                atom_name,
                None
                if atomic_number == 0
                else Element.getByAtomicNumber(atomic_number),
                residue,
            )
            positions.append(
                mm.Vec3(
                    float(values["r_m_x_coord"]),
                    float(values["r_m_y_coord"]),
                    float(values["r_m_z_coord"]),
                )
            )
            original_indices.append(index)
            formal_charges.append(
                self._chemical_value(values[formal_charge_column])
                if formal_charge_column
                else None
            )
        for row in bond_rows:
            values = dict(zip(bond_columns, row, strict=True))
            order = (
                self._chemical_value(values[bond_order_column])
                if bond_order_column
                else None
            )
            try:
                topology.addBond(
                    atoms[int(values["i_m_from"])],
                    atoms[int(values["i_m_to"])],
                    order=_topology_bond_order(order),
                )
            except KeyError as error:
                raise ValueError(
                    f"MAE bond references unknown atom index: {error.args[0]}"
                ) from error
            bond_orders.append(order)
        data = _structure_data(
            topology,
            positions * unit.angstrom,
            original_indices=original_indices,
            formal_charges=formal_charges,
            bond_orders=bond_orders,
        )
        super().__init__(
            data.topology, data.positions, data.atom_metadata, data.bond_metadata
        )

    @staticmethod
    def _first_present(columns: list[str], *choices: str) -> str | None:
        return next((choice for choice in choices if choice in columns), None)

    @classmethod
    def _read_table(
        cls, contents: str, table_name: str
    ) -> tuple[list[str], list[list[str]]]:
        match = re.search(
            rf"\b{table_name}\[\d+\]\s*\{{(.*?)\n\s*\}}", contents, flags=re.DOTALL
        )
        if match is None:
            raise ValueError(f"MAE file does not contain a {table_name} table.")
        # A block is "headers ::: rows :::", so the closing separator has to be
        # dropped as well; leaving it in adds one stray token and makes every
        # row-length check fail.
        sections = match.group(1).split(":::")
        if len(sections) < 2:
            raise ValueError(
                f"MAE {table_name} table is missing its ::: separator."
            )
        headers, rows = sections[0], sections[1]
        columns = cls._token_pattern.findall(headers)
        values = cls._token_pattern.findall(rows)
        if not columns:
            raise ValueError(f"MAE {table_name} table has malformed rows.")
        if len(values) % (len(columns) + 1) == 0:
            indexed_rows = [
                values[index : index + len(columns) + 1]
                for index in range(0, len(values), len(columns) + 1)
            ]
            if all(
                row[0] == str(index) for index, row in enumerate(indexed_rows, start=1)
            ):
                return columns, [
                    [cls._parse_value(value) for value in row[1:]]
                    for row in indexed_rows
                ]
        if len(values) % len(columns) != 0:
            raise ValueError(f"MAE {table_name} table has malformed rows.")
        return columns, [
            [
                cls._parse_value(value)
                for value in values[index : index + len(columns)]
            ]
            for index in range(0, len(values), len(columns))
        ]

    @staticmethod
    def _parse_value(value: str) -> str:
        """Unquote a MAE value, dropping the PDB column padding it carries.

        Names arrive padded to their PDB field width, so a residue reads as
        "GLU " and an atom as " N  ". Those spaces are formatting, not part of
        the name, and would stop any force-field template from matching.
        """
        if not value.startswith('"'):
            return value
        return ast.literal_eval(value).strip()

    @staticmethod
    def _chemical_value(value: str) -> int | float | str:
        """Convert numeric MAE chemistry fields while retaining unusual values."""
        try:
            number = float(value)
        except ValueError:
            return value
        return int(number) if number.is_integer() else number


class GROReader(StructureData):
    """Read a GRO structure and infer standard protein bonds."""

    _ion_elements = {
        "CL": "Cl",
        "NA": "Na",
        "K": "K",
        "MG": "Mg",
        "CA": "Ca",
        "ZN": "Zn",
    }

    def __init__(self, path: Path) -> None:
        gro = app.GromacsGroFile(str(path))
        topology = app.Topology()
        topology.setPeriodicBoxVectors(gro.getPeriodicBoxVectors())
        chain = topology.addChain()
        previous_residue = None
        residue = None
        for atom_name, residue_name, residue_id in zip(
            gro.atomNames,
            gro.residueNames,
            gro.residueIds,
            strict=True,
        ):
            residue_key = (residue_id, residue_name)
            if residue_key != previous_residue:
                residue = topology.addResidue(residue_name, chain, id=str(residue_id))
                previous_residue = residue_key
            topology.addAtom(atom_name, self._element(atom_name, residue_name), residue)
        topology.createStandardBonds()
        topology.createDisulfideBonds(gro.getPositions())
        data = _structure_data(topology, gro.getPositions())
        super().__init__(
            data.topology, data.positions, data.atom_metadata, data.bond_metadata
        )

    @classmethod
    def _element(cls, atom_name: str, residue_name: str):
        name = atom_name.strip().lstrip("0123456789")
        residue_element = cls._ion_elements.get(residue_name.strip().upper())
        if residue_element:
            return Element.getBySymbol(residue_element)
        if name[:2].upper() in {"CL", "BR"}:
            return Element.getBySymbol(name[:2].capitalize())
        if name and name[0].upper() in {"H", "C", "N", "O", "P", "S", "F", "I"}:
            return Element.getBySymbol(name[0].upper())
        return None


def read_structure(path: Path) -> StructureData:
    """Read a PDB, structural DMS, structural MAE, or GRO input."""
    if path.suffix.lower() == ".pdb":
        pdb = app.PDBFile(str(path))
        return _structure_data(pdb.topology, pdb.positions)
    if path.suffix.lower() == ".dms":
        return DMSReader(path)
    if path.suffix.lower() == ".mae":
        return MAEReader(path)
    if path.suffix.lower() == ".gro":
        return GROReader(path)
    raise ValueError("Input structure must use a .pdb, .dms, .mae, or .gro extension.")
