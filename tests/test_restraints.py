"""Backbone dihedral restraints."""

from __future__ import annotations

import math

import openmm as mm
import pytest
from openmm import app, unit
from openmm.app.element import Element

from ommflow.lib.restraints import (
    DIHEDRAL_RESTRAINTS,
    FOURIER_TERMS,
    _backbone_torsions,
    _fourier_terms,
    add_dihedral_restraints,
)


def _torsion_energy(theta0: float, strength: float, theta: float) -> float:
    """Energy of the restraint at one angle, evaluated by OpenMM itself."""
    system = mm.System()
    for _ in range(4):
        system.addParticle(12.0)
    force = mm.PeriodicTorsionForce()
    for periodicity, k in _fourier_terms(strength):
        phase = (periodicity * (theta0 + math.pi)) % (2 * math.pi)
        force.addTorsion(0, 1, 2, 3, periodicity, phase, k)
    system.addForce(force)
    context = mm.Context(
        system, mm.VerletIntegrator(1.0), mm.Platform.getPlatformByName("Reference")
    )
    context.setPositions(
        [
            mm.Vec3(1, 0, 0),
            mm.Vec3(0, 0, 0),
            mm.Vec3(0, 0, 1),
            mm.Vec3(math.cos(theta), math.sin(theta), 1),
        ]
    )
    return (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )


def test_the_well_has_its_minimum_on_the_reference_angle() -> None:
    """A positive strength must not put the minimum 180 degrees away."""
    theta0 = math.radians(-57.0)
    at_reference = _torsion_energy(theta0, 20.0, theta0)
    opposite = _torsion_energy(theta0, 20.0, theta0 + math.pi)
    assert at_reference < opposite

    scanned = [
        (_torsion_energy(theta0, 20.0, theta0 + math.radians(offset)), offset)
        for offset in range(-180, 180, 3)
    ]
    assert min(scanned)[1] == 0

    # Monotonic climb away from the reference, over the useful range.
    for offset in (5, 10, 30, 60):
        assert _torsion_energy(theta0, 20.0, theta0 + math.radians(offset)) > at_reference


def test_fourier_coefficients_follow_the_documented_series() -> None:
    terms = list(_fourier_terms(20.0))
    assert len(terms) == FOURIER_TERMS
    for periodicity, k in terms:
        expected = -20.0 * (-1) ** periodicity / math.factorial(periodicity)
        assert k == pytest.approx(expected)
    # The leading term is positive when the strength is given as a magnitude.
    assert terms[0] == (1, pytest.approx(20.0))


def _peptide(residue_names: tuple[str, ...], bond_chain: bool = True) -> app.Topology:
    topology = app.Topology()
    chain = topology.addChain("A")
    previous = None
    for index, name in enumerate(residue_names):
        residue = topology.addResidue(name, chain, id=str(index + 1))
        atoms = {
            atom_name: topology.addAtom(
                atom_name, Element.getBySymbol(atom_name[0]), residue
            )
            for atom_name in ("N", "CA", "C")
        }
        topology.addBond(atoms["N"], atoms["CA"])
        topology.addBond(atoms["CA"], atoms["C"])
        if previous is not None and bond_chain:
            topology.addBond(previous["C"], atoms["N"])
        previous = atoms
    return topology


def test_termini_and_chain_breaks_are_not_restrained() -> None:
    """phi needs the preceding peptide bond, psi the following one."""
    torsions = list(_backbone_torsions(_peptide(("ALA",) * 4)))
    # Four residues: psi on 1..3 and phi on 2..4, so six torsions, not eight.
    assert len(torsions) == 6

    # With no peptide bonds at all, there is nothing to restrain across.
    assert list(_backbone_torsions(_peptide(("ALA",) * 4, bond_chain=False))) == []

    # Non-protein residues are ignored entirely.
    assert list(_backbone_torsions(_peptide(("LIG", "UNK")))) == []


def test_restraints_are_added_to_the_system_and_cover_every_torsion() -> None:
    topology = _peptide(("ALA",) * 3)
    positions = [mm.Vec3(0.1 * i, 0.05 * (i % 3), 0.02 * i) for i in range(9)]
    system = mm.System()
    for _ in range(9):
        system.addParticle(12.0)

    records, description = add_dihedral_restraints(
        system, topology, positions, "bb", 20.0
    )
    assert len(records) == 4
    assert "3 protein residues" in description
    assert {record["angle"] for record in records} == {"phi", "psi"}
    assert all(len(record["atom_indices"]) == 4 for record in records)
    assert all(record["atom_names"] in (["C", "N", "CA", "C"], ["N", "CA", "C", "N"])
               for record in records)

    force = [f for f in system.getForces() if f.getName() == "DihedralRestraint"][0]
    assert force.getNumTorsions() == len(records) * FOURIER_TERMS


def test_disabled_and_invalid_modes() -> None:
    system = mm.System()
    assert add_dihedral_restraints(system, _peptide(("ALA",)), [], "none", 20.0) == (
        [],
        "disabled",
    )
    assert system.getNumForces() == 0
    assert set(DIHEDRAL_RESTRAINTS) == {"none", "bb", "ss"}

    with pytest.raises(ValueError, match="Unknown dihedral restraint"):
        add_dihedral_restraints(system, _peptide(("ALA",)), [], "helix", 20.0)

    topology = _peptide(("ALA",) * 3)
    positions = [mm.Vec3(0.1 * i, 0.0, 0.0) for i in range(9)]
    for _ in range(9):
        system.addParticle(12.0)
    with pytest.raises(ValueError, match="nonzero"):
        add_dihedral_restraints(system, topology, positions, "bb", 0.0)

    # The correct sign is always negative, so the magnitude is what counts and
    # -20 (as written in a Desmond or GROMACS restraint) is accepted verbatim.
    positive, _ = add_dihedral_restraints(mm.System(), topology, positions, "bb", 20.0)
    negative, _ = add_dihedral_restraints(mm.System(), topology, positions, "bb", -20.0)
    assert [r["reference_degrees"] for r in positive] == [
        r["reference_degrees"] for r in negative
    ]
    assert list(_fourier_terms(20.0)) == list(_fourier_terms(-20.0))


def test_secondary_structure_selection_needs_mdtraj() -> None:
    pytest.importorskip("mdtraj", reason="ss mode assigns structure with MDTraj")
    from ommflow.lib.restraints import _secondary_structure_residues

    # A straight extended chain has no helix or sheet to find.
    topology = _peptide(("ALA",) * 6)
    import numpy

    coordinates = numpy.array(
        [[0.15 * i, 0.0, 0.0] for i in range(topology.getNumAtoms())]
    )
    assert _secondary_structure_residues(topology, coordinates) == set()
