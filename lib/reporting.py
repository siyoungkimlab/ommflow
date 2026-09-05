"""Interval validation and OpenMM reporters."""

from __future__ import annotations

from pathlib import Path

from openmm import app, unit

from ommflow.lib.performance import (
    PerformanceReporter,
    PerformanceTracker,
    TimedReporter,
)
from ommflow.lib.restart import RunPaths


def steps_for(duration, timestep) -> int:
    """Return the integral number of steps for a positive exact duration."""
    raw_steps = duration / timestep
    steps = round(raw_steps)
    if steps < 1:
        raise ValueError("Each simulation duration must be at least one integration step.")
    if abs(raw_steps - steps) > 1e-6:
        raise ValueError("Simulation intervals must be exact multiples of the timestep.")
    return steps


def validate_intervals(args, timestep) -> tuple[int, int, int, int, int]:
    """Validate all durations and return equilibration and production step counts."""
    equilibration_steps = steps_for(
        args.equilibration_ns * unit.nanoseconds, timestep
    )
    equilibration_report_steps = steps_for(
        args.equilibration_report_interval_ns * unit.nanoseconds, timestep
    )
    production_target_steps = steps_for(
        args.production_ns * unit.nanoseconds, timestep
    )
    production_report_steps = steps_for(
        args.production_report_interval_ns * unit.nanoseconds, timestep
    )
    steps_for(args.performance_interval_ns * unit.nanoseconds, timestep)
    checkpoint_steps = steps_for(
        args.checkpoint_interval_ns * unit.nanoseconds, timestep
    )
    if production_report_steps % checkpoint_steps != 0:
        raise ValueError(
            "--checkpoint-interval-ns must divide "
            "--production-report-interval-ns so every reported frame has a "
            "matching checkpoint for safe restart appending."
        )
    validate_monitor_interval(args, timestep, checkpoint_steps)
    return (
        equilibration_steps,
        equilibration_report_steps,
        production_target_steps,
        production_report_steps,
        checkpoint_steps,
    )


def validate_monitor_interval(args, timestep, checkpoint_steps: int) -> int | None:
    """Return monitor steps and require every check to be checkpoint-aligned."""
    if not getattr(args, "early_stop", False):
        return None
    monitor_steps = steps_for(args.monitor_interval_ns * unit.nanoseconds, timestep)
    if monitor_steps % checkpoint_steps != 0:
        raise ValueError(
            "--monitor-interval-ns must equal or be an integer multiple of "
            "--checkpoint-interval-ns so every monitor point has a checkpoint."
        )
    return monitor_steps


def _state_reporter(path: Path, report_steps: int, append: bool = False):
    return app.StateDataReporter(
        str(path),
        report_steps,
        append=append,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        density=True,
        volume=True,
        separator=",",
    )


def add_equilibration_reporters(
    simulation: app.Simulation, paths: RunPaths, report_steps: int
) -> None:
    """Attach fresh reporters for the NVT and NPT equilibration phases."""
    simulation.reporters.append(app.DCDReporter(str(paths.equilibration_dcd), report_steps))
    simulation.reporters.append(_state_reporter(paths.equilibration_csv, report_steps))


def add_production_reporters(
    simulation: app.Simulation,
    paths: RunPaths,
    report_steps: int,
    checkpoint_steps: int,
    append: bool,
    tracker: PerformanceTracker | None = None,
    performance_steps: int | None = None,
    timestep_ns: float | None = None,
) -> None:
    """Attach production reporters, safely appending on restart.

    With a tracker, each reporter is wrapped so the wall time it costs is
    attributed to it, and a performance row is written every
    ``performance_steps``.
    """
    reporters = [
        ("trajectory", app.DCDReporter(str(paths.trajectory_dcd), report_steps, append=append)),
        ("state", _state_reporter(paths.state_csv, report_steps, append=append)),
        ("checkpoint", app.CheckpointReporter(str(paths.checkpoint), checkpoint_steps)),
    ]
    for task, reporter in reporters:
        simulation.reporters.append(
            reporter if tracker is None else TimedReporter(reporter, tracker, task)
        )
    if tracker is not None and performance_steps and timestep_ns:
        simulation.reporters.append(
            PerformanceReporter(
                paths.performance_csv,
                performance_steps,
                tracker,
                timestep_ns,
                append=append,
                initial_step=simulation.currentStep,
            )
        )
