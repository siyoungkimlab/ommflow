"""Backbone dihedral restraints holding phi and psi near the input structure.

The restraint is a truncated Fourier well built from periodic torsion terms:

    V(theta) = sum_{i=1..N} K (-1)^i / i! * [1 + cos(i (theta - theta0 - pi))]

which is what engines that express torsions only as cosine series need. OpenMM
uses exactly that form, ``E = k (1 + cos(n theta - phase))``, so each term maps
onto one ``PeriodicTorsionForce`` entry with periodicity ``i`` and phase
``i (theta0 + pi)``.

``K`` must be negative for the well to sit on ``theta0``; a positive K puts the
minimum at ``theta0 + 180``, restraining the backbone to the opposite
conformation. Only the magnitude of the configured strength is used, so both
``20`` and ``-20`` give the same, correct well.
"""

from __future__ import annotations

import math

import openmm as mm
from openmm import app, unit

from ommflow.lib.ligands import PROTEIN_RESIDUE_NAMES


DIHEDRAL_RESTRAINTS = ("none", "bb", "ss")
# Six terms reproduce a well that is flat far from theta0 and harmonic near it.
FOURIER_TERMS = 6
# DSSP codes counted as secondary structure, in mdtraj's simplified alphabet.
SECONDARY_STRUCTURE_CODES = frozenset({"H", "E"})
_BACKBONE_ATOMS = ("N", "CA", "C")


def add_dihedral_restraints(
    system: mm.System,
    topology: app.Topology,
    positions,
    mode: str,
    strength_kj: float,
) -> tuple[list[dict[str, object]], str]:
    """Restrain protein phi and psi to their values in ``positions``.

    Returns one record per restrained torsion and a description of what was
    selected. The force is added to ``system``, so it is serialized into
    system.xml and survives restarts without being reapplied.
    """
    if mode not in DIHEDRAL_RESTRAINTS:
        raise ValueError(
            f"Unknown dihedral restraint '{mode}'. Choose one of: "
            + ", ".join(DIHEDRAL_RESTRAINTS)
        )
    if mode == "none":
        return [], "disabled"
    if not math.isfinite(strength_kj) or strength_kj == 0:
        raise ValueError("The dihedral restraint strength must be a nonzero number.")

    coordinates = _coordinates_nm(positions)
    torsions = list(_backbone_torsions(topology))
    if mode == "ss":
        selected = _secondary_structure_residues(topology, coordinates)
        torsions = [t for t in torsions if t[0].index in selected]
        description = (
            f"{len(selected)} residues in helices and sheets"
            if selected
            else "no residues in helices or sheets"
        )
    else:
        description = f"{len({t[0].index for t in torsions})} protein residues"

    if not torsions:
        raise ValueError(
            f"--dihedral-restraint {mode} selected no backbone torsions "
            f"({description}); there is nothing to restrain."
        )

    force = mm.PeriodicTorsionForce()
    force.setName("DihedralRestraint")
    # One pass; a lookup per atom would rescan the whole topology 4n times.
    atom_names = {atom.index: atom.name for atom in topology.atoms()}
    records: list[dict[str, object]] = []
    for residue, kind, atoms in torsions:
        reference = _dihedral(coordinates, atoms)
        for periodicity, k in _fourier_terms(strength_kj):
            phase = (periodicity * (reference + math.pi)) % (2 * math.pi)
            force.addTorsion(*atoms, periodicity, phase, k)
        records.append(
            {
                "angle": kind,
                "chain_id": residue.chain.id,
                "residue_name": residue.name,
                "residue_id": residue.id,
                "atom_indices": list(atoms),
                "atom_names": [atom_names[index] for index in atoms],
                "reference_degrees": round(math.degrees(reference), 3),
            }
        )
    system.addForce(force)
    return records, description



def _fourier_terms(strength_kj: float):
    """Yield (periodicity, force constant) for the well of the given strength."""
    # Negative K places the minimum on theta0 rather than opposite it.
    amplitude = -abs(strength_kj)
    for i in range(1, FOURIER_TERMS + 1):
        yield i, amplitude * (-1) ** i / math.factorial(i)


def _backbone_torsions(topology: app.Topology):
    """Yield (residue, "phi"/"psi", atom indices) for every protein torsion.

    A torsion is only produced where the peptide bond to the neighbouring
    residue actually exists, so chain breaks and termini are skipped rather
    than restrained across a gap.
    """
    bonded = {
        frozenset((bond[0].index, bond[1].index)) for bond in topology.bonds()
    }
    for chain in topology.chains():
        residues = [
            residue
            for residue in chain.residues()
            if residue.name.upper() in PROTEIN_RESIDUE_NAMES
        ]
        backbone = [_named_atoms(residue) for residue in residues]
        for position, residue in enumerate(residues):
            current = backbone[position]
            if current is None:
                continue
            previous = backbone[position - 1] if position > 0 else None
            following = backbone[position + 1] if position + 1 < len(residues) else None
            if previous is not None and frozenset(
                (previous["C"], current["N"])
            ) in bonded:
                yield residue, "phi", (
                    previous["C"],
                    current["N"],
                    current["CA"],
                    current["C"],
                )
            if following is not None and frozenset(
                (current["C"], following["N"])
            ) in bonded:
                yield residue, "psi", (
                    current["N"],
                    current["CA"],
                    current["C"],
                    following["N"],
                )


def _named_atoms(residue) -> dict[str, int] | None:
    """Return N, CA and C indices, or None if the residue lacks any of them."""
    found = {
        atom.name: atom.index
        for atom in residue.atoms()
        if atom.name in _BACKBONE_ATOMS
    }
    return found if len(found) == len(_BACKBONE_ATOMS) else None


def _secondary_structure_residues(topology: app.Topology, coordinates) -> set[int]:
    """Return topology residue indices assigned to a helix or a sheet."""
    mdtraj = _require_mdtraj()
    import numpy

    trajectory = mdtraj["Trajectory"](
        numpy.asarray(coordinates, dtype=numpy.float32)[numpy.newaxis, :, :],
        mdtraj["Topology"].from_openmm(topology),
    )
    codes = mdtraj["compute_dssp"](trajectory, simplified=True)[0]
    return {
        index
        for index, code in enumerate(codes)
        if code in SECONDARY_STRUCTURE_CODES
    }


def _require_mdtraj() -> dict[str, object]:
    """Import mdtraj only when secondary structure is actually requested."""
    try:
        import mdtraj
        from mdtraj import Topology, Trajectory, compute_dssp
    except ImportError as error:
        raise RuntimeError(
            "--dihedral-restraint ss assigns secondary structure with MDTraj, "
            "which is not installed. Install it with 'pip install mdtraj', or "
            "use --dihedral-restraint bb to restrain the whole backbone."
        ) from error
    return {
        "mdtraj": mdtraj,
        "Topology": Topology,
        "Trajectory": Trajectory,
        "compute_dssp": compute_dssp,
    }


def _coordinates_nm(positions):
    import numpy

    vectors = (
        positions.value_in_unit(unit.nanometer)
        if hasattr(positions, "value_in_unit")
        else positions
    )
    return numpy.asarray(vectors, dtype=float)


def _dihedral(coordinates, atoms: tuple[int, int, int, int]) -> float:
    """Return the torsion angle in radians, by the usual cross-product form."""
    import numpy

    p0, p1, p2, p3 = (coordinates[index] for index in atoms)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / numpy.linalg.norm(b1)
    v = b0 - numpy.dot(b0, b1) * b1
    w = b2 - numpy.dot(b2, b1) * b1
    return math.atan2(numpy.dot(numpy.cross(b1, v), w), numpy.dot(v, w))


def write_dihedral_restraints(path, records: list[dict[str, object]]) -> None:
    """Record which torsions were restrained, and to what angle."""
    import csv

    fields = (
        "angle",
        "chain_id",
        "residue_name",
        "residue_id",
        "atom_indices",
        "atom_names",
        "reference_degrees",
    )
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            row = dict(record)
            row["atom_indices"] = " ".join(str(i) for i in record["atom_indices"])
            row["atom_names"] = " ".join(record["atom_names"])
            writer.writerow(row)


def plot_dihedral_restraint(path, strength_kj: float) -> bool:
    """Draw the restraint well against deviation from the reference angle.

    matplotlib is a required dependency, so a failure here means a broken
    installation rather than a missing option. It is still not worth aborting a
    simulation over a diagnostic plot, so this reports and carries on.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy
    except ImportError as error:
        print(
            f"Warning: could not plot the restraint potential ({error}). "
            "matplotlib is a required dependency, so this suggests an "
            "incomplete installation; the restraint itself is unaffected."
        )
        return False

    deviation = numpy.linspace(-180.0, 180.0, 1441)
    radians = numpy.radians(deviation)
    energy = sum(
        k * (1 + numpy.cos(periodicity * (radians - math.pi)))
        for periodicity, k in _fourier_terms(strength_kj)
    )
    minimum = energy.min()
    depth = energy.max() - minimum
    # Curvature at the bottom, as the harmonic this well is equivalent to.
    step = numpy.radians(deviation[1] - deviation[0])
    centre = len(deviation) // 2
    curvature = (energy[centre + 1] - 2 * energy[centre] + energy[centre - 1]) / step**2

    figure, axes = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    axes.plot(deviation, energy - minimum, lw=2, color="#1f77b4",
              label=f"{FOURIER_TERMS}-term Fourier well")
    axes.plot(deviation, 0.5 * curvature * radians**2, ls="--", lw=1.2,
              color="#d62728",
              label=f"harmonic, k = {curvature:.0f} kJ/mol/rad$^2$")
    axes.axvline(0.0, color="0.4", lw=1, ls=":")
    axes.set_ylim(-0.05 * depth, 1.15 * depth)
    axes.set_xlim(-180, 180)
    axes.set_xticks(range(-180, 181, 60))
    axes.set_xlabel(r"deviation from the reference angle, $\theta-\theta_0$ (degrees)")
    axes.set_ylabel("restraint energy above the minimum (kJ/mol)")
    axes.set_title(
        f"Backbone dihedral restraint, {abs(strength_kj):g} kJ/mol "
        f"(K = {-abs(strength_kj):g})"
    )
    for offset in (10, 30, 60):
        index = int(numpy.argmin(numpy.abs(deviation - offset)))
        axes.annotate(
            f"{offset}°: {energy[index] - minimum:.1f}",
            xy=(offset, energy[index] - minimum),
            xytext=(offset + 6, energy[index] - minimum + 0.05 * depth),
            fontsize=8,
            color="0.3",
        )
    axes.legend(loc="upper center", frameon=False, fontsize=9)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return True
