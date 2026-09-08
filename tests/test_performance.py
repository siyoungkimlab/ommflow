"""Timing instrumentation for production tasks."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pytest
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


class _Clock:
    """A Simulation stand-in whose reported time carries realistic drift."""

    def __init__(self, picoseconds: float) -> None:
        self._picoseconds = picoseconds

    class _State:
        def __init__(self, picoseconds: float) -> None:
            self._picoseconds = picoseconds

        def getTime(self):  # noqa: N802 - OpenMM's interface
            return self._picoseconds * unit.picoseconds

    @property
    def context(self):
        clock = self

        class _Context:
            def getState(self, *args, **kwargs):  # noqa: N802 - OpenMM's interface
                return _Clock._State(clock._picoseconds)

        return _Context()


def test_step_count_survives_the_clock_drifting_over_a_long_run() -> None:
    """OpenMM adds the timestep once per step, so its clock loses exactness."""
    from ommflow.lib.simulation import _current_steps

    timestep = 0.002 * unit.picoseconds

    # Measured drift: about -2.5e-6 steps at 400k and -7.3e-6 at 600k. A 10 ns
    # run at 2 fs is 5 million steps, where it is larger still. All of these
    # used to raise "Checkpoint time is not an exact multiple of the timestep".
    for steps, drift in ((400_000, -2.5e-6), (600_000, -7.3e-6), (5_000_000, -6e-5)):
        picoseconds = (steps + drift) * 0.002
        assert _current_steps(_Clock(picoseconds), timestep) == steps

    # A genuinely mismatched timestep is still caught. A clock at 1000 ps read
    # with a 3 fs step is 333333.33 steps, which no drift could explain.
    with pytest.raises(ValueError, match="not a multiple of the timestep"):
        _current_steps(_Clock(1000.0), 0.003 * unit.picoseconds)
