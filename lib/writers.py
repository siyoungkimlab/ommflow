"""Writers for structure files that carry bonds alongside coordinates."""

from __future__ import annotations

from pathlib import Path

from openmm import unit

_MAE_VERSION_BLOCK = "{\n s_m_m2io_version\n:::\n 2.0.0\n}\n\n"

_ATOM_COLUMNS = (
    "i_m_residue_number",
    "s_m_chain_name",
    "s_m_pdb_residue_name",
    "s_m_pdb_atom_name",
    "i_m_atomic_number",
    "r_m_x_coord",
    "r_m_y_coord",
    "r_m_z_coord",
)

_BOND_COLUMNS = ("i_m_from", "i_m_to", "i_m_order")

_BOX_COLUMNS = (
    "r_chorus_box_ax",
    "r_chorus_box_ay",
    "r_chorus_box_az",
    "r_chorus_box_bx",
    "r_chorus_box_by",
    "r_chorus_box_bz",
    "r_chorus_box_cx",
    "r_chorus_box_cy",
    "r_chorus_box_cz",
)


def _quote(value: object) -> str:
    """Quote a MAE string value, escaping what the format reserves."""
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def _pdb_atom_name(atom) -> str:
    """Pad an atom name back into its four-character PDB column.

    ``s_m_pdb_atom_name`` is the PDB field, so a one-letter element is indented
    by one space exactly as ``PDBFile`` writes it; an unpadded name puts the
    element in the wrong column for anything that reads the padding back.
    """
    name = atom.name[:4]
    symbol = atom.element.symbol if atom.element is not None else ""
    if len(name) < 4 and name[:1].isalpha() and len(symbol) < 2:
        name = " " + name
    return f"{name:<4}"


def _residue_number(residue, fallback: int) -> int:
    """Use the residue ID when it is numeric, else its position in the file.

    Solvent chains built by ``Modeller`` are numbered, but an input residue ID
    can carry an insertion code, and ``i_m_residue_number`` is an integer.
    """
    try:
        return int(str(residue.id).strip())
    except (TypeError, ValueError):
        return fallback


def _angstrom_rows(positions) -> list[tuple[float, float, float]]:
    """Coordinates as plain Angstrom triples, whatever container they arrive in."""
    values = positions.value_in_unit(unit.angstrom)
    return [(float(row[0]), float(row[1]), float(row[2])) for row in values]


def _box_values(box_vectors) -> list[float] | None:
    if box_vectors is None:
        return None
    vectors = box_vectors.value_in_unit(unit.angstrom)
    return [float(component) for vector in vectors for component in vector]


def _bond_order(bond) -> int:
    """The bond order to record, defaulting to single.

    A MAE or DMS input keeps its orders through solvation, so a ligand still
    shows its double and triple bonds. Water, ions, and the rest of what the
    solvent builder adds carry no order because they are single bonds, and the
    column has to hold an integer either way.
    """
    if bond.order is None:
        return 1
    try:
        return int(bond.order)
    except (TypeError, ValueError):
        return 1


def write_mae(
    path: Path,
    topology,
    positions,
    *,
    title: str = "ommflow",
    box_vectors=None,
) -> None:
    """Write a topology and its coordinates as a Maestro structure file.

    Written alongside the PDB of the same stage because a PDB drops the bonds
    of anything without a standard residue template and cannot number more than
    99,999 atoms, both of which a solvated box routinely exceeds. Loading this
    file as the topology for ``equilibration.dcd`` or ``trajectory.dcd`` gives a
    viewer the real connectivity and, through the chorus box, the periodic cell.
    """
    coordinates = _angstrom_rows(positions)
    atoms = list(topology.atoms())
    if len(coordinates) != len(atoms):
        raise ValueError(
            f"Cannot write {path.name}: {len(atoms)} topology atoms but "
            f"{len(coordinates)} positions."
        )
    if box_vectors is None:
        box_vectors = topology.getPeriodicBoxVectors()
    box = _box_values(box_vectors)
    bonds = list(topology.bonds())
    atom_index = {atom.index: index for index, atom in enumerate(atoms, start=1)}

    with path.open("w", encoding="utf-8") as handle:
        handle.write(_MAE_VERSION_BLOCK)
        handle.write("f_m_ct {\n")
        for column in ("s_m_title",) + (_BOX_COLUMNS if box else ()):
            handle.write(f" {column}\n")
        handle.write(":::\n")
        handle.write(f" {_quote(title)}\n")
        for value in box or ():
            handle.write(f" {value:.6f}\n")

        handle.write(f" m_atom[{len(atoms)}] {{\n")
        for column in _ATOM_COLUMNS:
            handle.write(f"  {column}\n")
        handle.write("  :::\n")
        residue_number = 0
        previous_residue = None
        for index, (atom, (x, y, z)) in enumerate(
            zip(atoms, coordinates, strict=True), start=1
        ):
            residue = atom.residue
            if residue is not previous_residue:
                residue_number += 1
                previous_residue = residue
            # A blank chain is written as a space, the way a PDB writes it, and
            # a virtual site keeps its place in the file with atomic number 0
            # so the atom count still matches the trajectory it belongs to.
            chain_name = residue.chain.id or " "
            atomic_number = atom.element.atomic_number if atom.element else 0
            handle.write(
                f"  {index} {_residue_number(residue, residue_number)} "
                f"{_quote(chain_name)} "
                f"{_quote(f'{residue.name[:4]:<4}')} "
                f"{_quote(_pdb_atom_name(atom))} "
                f"{atomic_number} "
                f"{x:.6f} {y:.6f} {z:.6f}\n"
            )
        handle.write("  :::\n")
        handle.write(" }\n")

        if bonds:
            handle.write(f" m_bond[{len(bonds)}] {{\n")
            for column in _BOND_COLUMNS:
                handle.write(f"  {column}\n")
            handle.write("  :::\n")
            for index, bond in enumerate(bonds, start=1):
                handle.write(
                    f"  {index} {atom_index[bond[0].index]} "
                    f"{atom_index[bond[1].index]} {_bond_order(bond)}\n"
                )
            handle.write("  :::\n")
            handle.write(" }\n")
        handle.write("}\n")
