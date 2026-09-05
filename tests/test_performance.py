"""Timing instrumentation for production tasks."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

from openmm import app, unit

from ommflow.lib.performance import (
    PERFORMANCE_FIELDS,
    PerformanceReporter,
    PerformanceTracker,
    TimedReporter,
)


class _Recorder:
    """A minimal object satisfying OpenMM's reporter interface."""

    def __init__(self) -> None:
        self.reports = 0

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM's interface
        return (7, True, False, False, False, None)

    def report(self, simulation, state):
        self.reports += 1


def test_tracker_attributes_time_to_the_task_that_spent_it() -> None:
    tracker = PerformanceTracker()
    tracker.start()
    with tracker.measure("checkpoint"):
        sum(range(20000))
    with tracker.measure("monitor"):
        sum(range(20000))
    assert tracker.totals["checkpoint"] > 0
    assert tracker.totals["monitor"] > 0
    assert tracker.totals["trajectory"] == 0
    assert tracker.elapsed() >= sum(tracker.totals.values())
    assert "integration" in tracker.summary()


def test_timed_reporter_delegates_and_charges_its_task() -> None:
    tracker = PerformanceTracker()
    tracker.start()
    recorder = _Recorder()
    wrapped = TimedReporter(recorder, tracker, "trajectory")

    # The description passes through untouched so OpenMM still builds the
    # State the real reporter asked for.
    assert wrapped.describeNextReport(None) == (7, True, False, False, False, None)
    wrapped.report(None, None)
    assert recorder.reports == 1
    assert tracker.totals["trajectory"] > 0
    assert wrapped.reports == 1  # attribute access falls through


def test_performance_rows_measure_rate_from_where_the_run_resumed(
    tmp_path: Path,
) -> None:
    tracker = PerformanceTracker()
    tracker.start()
    path = tmp_path / "performance.csv"
    reporter = PerformanceReporter(
        path, 1000, tracker, timestep_ns=2e-6, initial_step=10_000
    )
    reporter.report(SimpleNamespace(currentStep=11_000), None)

    rows = list(csv.DictReader(path.open()))
    assert list(rows[0]) == list(PERFORMANCE_FIELDS)
    assert rows[0]["step"] == "11000"
    # 11000 steps x 2e-6 ns; the rate covers only the 1000 steps since resuming,
    # not all 11000, so it stays finite and sane rather than exploding.
    assert float(rows[0]["production_time_ns"]) == 0.022
    assert 0 < float(rows[0]["ns_per_day"]) < 1e7


def test_hydrogen_mass_repartitioning_never_touches_water() -> None:
    """Water is rigid, so HMR must skip it for 3-, 4-, and 5-site models."""
    for protein_xml, water_xml, model in (
        ("amber19/protein.ff19SB.xml", "amber19/opc.xml", "tip4pew"),
        ("charmm36_2024.xml", "charmm36_2024/tip5p.xml", "tip5p"),
    ):
        forcefield = app.ForceField(protein_xml, water_xml)
        modeller = app.Modeller(app.Topology(), [])
        modeller.addSolvent(
            forcefield, model=model, boxSize=(2, 2, 2) * unit.nanometer
        )
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            constraints=app.HBonds,
            rigidWater=True,
            hydrogenMass=4.0 * unit.amu,
        )
        hydrogens = {
            round(system.getParticleMass(atom.index).value_in_unit(unit.amu), 3)
            for atom in modeller.topology.atoms()
            if atom.element is not None and atom.element.symbol == "H"
        }
        assert hydrogens == {1.008}, (model, hydrogens)
        assert any(
            system.isVirtualSite(index)
            for index in range(system.getNumParticles())
        )
