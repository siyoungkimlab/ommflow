"""Simulation phases for the OpenMM protein MD workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import openmm as mm
from openmm import app, unit

from ommflow.lib.monitor import (
    MonitorDefinition,
    append_monitor_row,
    detachment_update,
    initialize_monitor,
    load_components,
    load_pocket,
    load_status,
    monitor_measurement,
    monitor_status,
    restore_consecutive_count,
    select_monitor_target,
    validate_monitor_topology,
    validate_saved_monitor,
    write_pocket,
    write_status,
)
from ommflow.lib.reporting import (
    add_equilibration_reporters,
    add_production_reporters,
    steps_for,
    validate_intervals,
    validate_monitor_interval,
)
from ommflow.lib.config import require_early_stop_for_selectors
from ommflow.lib.performance import PerformanceTracker
from ommflow.lib.platforms import (
    create_simulation,
    describe_platform,
    describe_system,
)
from ommflow.lib.restart import (
    RunPaths,
    commit_restart_settings,
    is_restart,
    prepare_restart,
    save_equilibrated_artifacts,
    save_new_configuration,
)
from ommflow.lib.system_builder import build_solvated_system
from ommflow.lib.writers import write_mae


def _simulation(
    topology: app.Topology,
    system: mm.System,
    integrator: object,
    args: argparse.Namespace,
) -> app.Simulation:
    """Create the simulation and report the platform and system it will use."""
    simulation, notes = create_simulation(
        topology,
        system,
        integrator,
        args.platform,
        args.precision,
        getattr(args, "precision_specified", False),
    )
    print(describe_system(topology, system))
    print(describe_platform(simulation.context))
    for note in notes:
        print(f"Note: {note}")
    return simulation


def _write_positions(topology: app.Topology, positions, path: Path) -> None:
    with path.open("w") as handle:
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)


def run_equilibration(
    simulation: app.Simulation,
    system: mm.System,
    modeller: app.Modeller,
    args: argparse.Namespace,
    temperature,
    pressure,
    paths: RunPaths,
    equilibration_steps: int,
    report_steps: int,
) -> None:
    """Minimize, then run NVT and NPT equilibration for a new system."""
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature, args.seed)

    print(f"Using {args.proteinff} protein parameters with {args.waterff} water.")
    add_equilibration_reporters(simulation, paths, report_steps)
    print(f"Running {args.equilibration_ns:g} ns NVT equilibration...")
    simulation.step(equilibration_steps)
    print(f"Running {args.equilibration_ns:g} ns NPT equilibration...")
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
    simulation.context.reinitialize(preserveState=True)
    simulation.step(equilibration_steps)
    simulation.reporters.clear()
    equilibrated_state = simulation.context.getState(getPositions=True)
    _write_positions(modeller.topology, equilibrated_state.getPositions(), paths.equilibrated_pdb)
    # The production trajectory starts from these coordinates, so this is the
    # structure to load it against; the box comes from the state because NPT
    # has resized the one the topology was solvated with.
    write_mae(
        paths.equilibrated_mae,
        modeller.topology,
        equilibrated_state.getPositions(),
        title="equilibrated",
        box_vectors=equilibrated_state.getPeriodicBoxVectors(),
    )
    simulation.context.setTime(0 * unit.picoseconds)
    # Production reporting counts from the first production step, so state.csv
    # steps line up with monitor.csv instead of carrying an equilibration offset.
    simulation.currentStep = 0
    save_equilibrated_artifacts(system, simulation.integrator, paths)


def run_production(
    simulation: app.Simulation,
    modeller: app.Modeller,
    args: argparse.Namespace,
    paths: RunPaths,
    production_target_steps: int,
    report_steps: int,
    checkpoint_steps: int,
    append: bool,
    monitor: MonitorDefinition | None = None,
    monitor_steps: int | None = None,
    performance_steps: int | None = None,
) -> None:
    """Run only the remaining steps needed to reach the absolute target."""
    timestep = simulation.integrator.getStepSize()
    tracker = PerformanceTracker()
    current_steps = _current_steps(simulation, timestep)
    prior_status = (
        load_status(paths.status_json)
        if monitor is not None or paths.status_json.is_file()
        else None
    )
    if monitor is not None:
        _validate_terminal_status(prior_status, args, current_steps)
        if _is_detached_terminal_status(prior_status, args):
            print(
                "Production remains stopped after confirmed target detachment. "
                "Disable --early-stop or increase --production-ns to continue."
            )
            return
        if _is_target_terminal_status(prior_status, args):
            print(f"Production target of {args.production_ns:g} ns has already been reached.")
            return
    add_production_reporters(
        simulation,
        paths,
        report_steps,
        checkpoint_steps,
        append=append,
        tracker=tracker,
        performance_steps=performance_steps,
        timestep_ns=timestep.value_in_unit(unit.nanoseconds),
    )
    tracker.start()
    production_steps = production_target_steps - current_steps
    detached = False
    consecutive_count = (
        restore_consecutive_count(prior_status, paths.monitor_csv, current_steps)
        if monitor is not None
        else 0
    )
    if production_steps > 0:
        remaining_time = production_steps * timestep
        print(
            f"Running {remaining_time.value_in_unit(unit.nanoseconds):g} ns "
            f"to reach {args.production_ns:g} ns of NPT production..."
        )
        if monitor is None:
            simulation.step(production_steps)
            current_steps += production_steps
        else:
            if monitor_steps is None:
                raise AssertionError("An enabled monitor must have an interval.")
            detached, consecutive_count, current_steps = _run_monitored_production(
                simulation,
                monitor,
                args,
                paths,
                current_steps,
                production_target_steps,
                monitor_steps,
                consecutive_count,
                tracker,
            )
    else:
        print(f"Production target of {args.production_ns:g} ns has already been reached.")
    if tracker.started:
        print(tracker.summary())
    state = simulation.context.getState(getPositions=True)
    _write_positions(modeller.topology, state.getPositions(), paths.final_pdb)
    # Rewritten at the end of every segment, so the box tracks the barostat
    # rather than the one production started from.
    write_mae(
        paths.final_mae,
        modeller.topology,
        state.getPositions(),
        title="final",
        box_vectors=state.getPeriodicBoxVectors(),
    )
    final_steps = current_steps
    final_time_ns = state.getTime().value_in_unit(unit.nanoseconds)
    if monitor is not None:
        simulation.saveCheckpoint(str(paths.checkpoint))
        outcome = "detached" if detached else "target_reached"
        write_status(
            paths.status_json,
            monitor_status(
                outcome,
                monitor,
                args.production_ns,
                final_time_ns,
                final_steps,
                consecutive_count,
            ),
        )
        if detached:
            print(
                f"Confirmed detachment of {monitor.ligand_id} at "
                f"{final_time_ns:g} ns; production stopped early."
            )
    elif prior_status is not None:
        _mark_monitor_disabled_completion(
            paths,
            prior_status,
            args.production_ns,
            final_time_ns,
            final_steps,
        )


# OpenMM advances its clock by adding the timestep once per step, so the time
# it reports drifts from an exact multiple by roughly one part in 1e9 per step:
# about 2.5e-6 steps after 400 thousand, and 7e-6 after 600 thousand. A
# tolerance tight enough to be meaningful for a genuinely mismatched timestep,
# which is off by whole steps or more, still has to be loose enough to survive
# that drift over a long run.
_STEP_COUNT_TOLERANCE = 0.01


def _current_steps(simulation: app.Simulation, timestep) -> int:
    """Convert the context clock into a production integration-step count.

    Only used once per submission, to learn where a restored checkpoint sits.
    After that the count is carried as an integer, because re-deriving it from
    the clock accumulates the drift described above.
    """
    current_time = simulation.context.getState().getTime()
    current_steps_raw = current_time / timestep
    current_steps = round(current_steps_raw)
    if abs(current_steps_raw - current_steps) > _STEP_COUNT_TOLERANCE:
        raise ValueError(
            f"Checkpoint time {current_time} is not a multiple of the timestep "
            f"{timestep}; it is {current_steps_raw:g} steps. The checkpoint was "
            "probably written with a different --integration-fs."
        )
    return current_steps


def _run_monitored_production(
    simulation: app.Simulation,
    monitor: MonitorDefinition,
    args: argparse.Namespace,
    paths: RunPaths,
    current_steps: int,
    production_target_steps: int,
    monitor_steps: int,
    consecutive_count: int,
    tracker: PerformanceTracker,
) -> tuple[bool, int, int]:
    """Step precisely to monitor points, checkpointing before each CSV row."""
    while current_steps < production_target_steps:
        next_monitor_step = (current_steps // monitor_steps + 1) * monitor_steps
        step_count = min(
            next_monitor_step - current_steps, production_target_steps - current_steps
        )
        simulation.step(step_count)
        # Incremented rather than re-derived from the clock: exact, and immune
        # to the drift that made a long run fail at a monitor point.
        current_steps += step_count
        if current_steps % monitor_steps:
            continue
        state = simulation.context.getState(getPositions=True)
        with tracker.measure("checkpoint"):
            simulation.saveCheckpoint(str(paths.checkpoint))
        with tracker.measure("monitor"):
            minimum_distance, contact_count = monitor_measurement(
                monitor.ligand_heavy_atom_indices,
                monitor.pocket_atom_indices,
                state.getPositions(),
                state.getPeriodicBoxVectors(),
                args.contact_cutoff_nm,
            )
        consecutive_count, confirmed = detachment_update(
            minimum_distance,
            contact_count,
            args.detach_cutoff_nm,
            consecutive_count,
            args.confirmation_checks,
        )
        production_time_ns = state.getTime().value_in_unit(unit.nanoseconds)
        append_monitor_row(
            paths.monitor_csv,
            {
                "production_time_ns": f"{production_time_ns:.12g}",
                "step": current_steps,
                "min_ligand_pocket_distance_nm": f"{minimum_distance:.12g}",
                "contact_count": contact_count,
                "initial_contact_count": monitor.initial_contact_count,
                "contact_fraction": (
                    f"{contact_count / monitor.initial_contact_count:.12g}"
                    if monitor.initial_contact_count
                    else ""
                ),
                "consecutive_detached_count": consecutive_count,
                "detached": str(confirmed).lower(),
            },
        )
        if confirmed and current_steps < production_target_steps:
            return True, consecutive_count, current_steps
        write_status(
            paths.status_json,
            monitor_status(
                "running",
                monitor,
                args.production_ns,
                production_time_ns,
                current_steps,
                consecutive_count,
            ),
        )
    return False, consecutive_count, current_steps


def _status_target(status: dict[str, object] | None) -> float | None:
    if status is None:
        return None
    target = status.get("target_production_ns")
    return float(target) if isinstance(target, int | float) else None


def _is_detached_terminal_status(
    status: dict[str, object] | None, args: argparse.Namespace
) -> bool:
    """A terminal detachment needs an explicit larger target to continue."""
    if status is None or status.get("outcome") != "detached":
        return False
    previous_target = _status_target(status)
    return (
        not args.production_ns_specified
        or previous_target is None
        or args.production_ns <= previous_target
    )


def _is_target_terminal_status(
    status: dict[str, object] | None, args: argparse.Namespace
) -> bool:
    """Avoid reopening a completed target when the checkpoint is already complete."""
    return (
        status is not None
        and status.get("outcome") == "target_reached"
        and _status_target(status) is not None
        and args.production_ns <= _status_target(status)
    )


def _validate_terminal_status(
    status: dict[str, object] | None, args: argparse.Namespace, current_steps: int
) -> None:
    """Reject malformed terminal status only when it conflicts with a restart."""
    if status is None:
        return
    outcome = status.get("outcome")
    if outcome not in {"running", "target_reached", "detached"}:
        raise ValueError(f"Unknown monitor status outcome: {outcome!r}")
    status_step = status.get("final_production_step")
    if isinstance(status_step, int) and status_step > current_steps:
        raise ValueError(
            "Monitor status is newer than checkpoint; restore the matching "
            "checkpoint before resuming."
        )


def _mark_monitor_disabled_completion(
    paths: RunPaths,
    status: dict[str, object],
    target_production_ns: float,
    final_time_ns: float,
    final_steps: int,
) -> None:
    """Record that a formerly monitored run reached target with monitoring off."""
    status.update(
        {
            "outcome": "target_reached",
            "target_production_ns": target_production_ns,
            "final_production_time_ns": final_time_ns,
            "final_production_step": final_steps,
            "consecutive_detached_count": 0,
            "early_stop_enabled": False,
        }
    )
    write_status(paths.status_json, status)


def run_workflow(args: argparse.Namespace) -> None:
    """Select new-run or restart flow, then perform production to its target."""
    if (
        args.padding_nm <= 0
        or args.cutoff_nm <= 0
        or args.saltM < 0
        or args.integration_fs <= 0
    ):
        raise ValueError(
            "--padding-nm, --cutoff-nm, and --integration-fs must be positive and "
            "--saltM cannot be negative."
        )
    paths = RunPaths(args.workdir)
    resume = is_restart(paths.workdir)
    if resume and not paths.workdir.is_dir():
        raise NotADirectoryError(
            f"Output path exists but is not a directory: {paths.workdir}"
        )
    if not resume:
        require_early_stop_for_selectors(args)
        if args.input_structure is None:
            raise ValueError("input_structure must be supplied for a new simulation.")
        if not args.input_structure.is_file():
            raise FileNotFoundError(
                f"Input structure does not exist: {args.input_structure}"
            )
        paths.workdir.mkdir(parents=True, exist_ok=True)
        modeller, system, integrator, temperature, pressure = build_solvated_system(
            args, paths.workdir, paths.solvated_pdb, paths.components_json
        )
        if args.early_stop:
            select_monitor_target(load_components(paths.components_json), args)
        save_new_configuration(args, paths)
    else:
        modeller, system, integrator = prepare_restart(args, paths)
        require_early_stop_for_selectors(args)
        temperature = pressure = None

    simulation = _simulation(modeller.topology, system, integrator, args)
    timestep = integrator.getStepSize()
    requested_timestep = args.integration_fs * unit.femtoseconds
    if (
        resume
        and args.integration_fs_specified
        and abs(timestep - requested_timestep).value_in_unit(unit.femtoseconds) > 1e-8
    ):
        raise ValueError(
            "Cannot change --integration-fs when resuming; it must match the "
            "timestep stored in integrator.xml."
        )
    (
        equilibration_steps,
        equilibration_report_steps,
        production_target_steps,
        production_report_steps,
        checkpoint_steps,
    ) = validate_intervals(args, timestep)
    monitor_steps = validate_monitor_interval(args, timestep, checkpoint_steps)
    performance_steps = steps_for(
        args.performance_interval_ns * unit.nanoseconds, timestep
    )

    if resume:
        simulation.loadCheckpoint(str(paths.checkpoint))
        print(f"Resuming from {paths.checkpoint}")
        monitor = None
        if args.early_stop:
            monitor = load_pocket(paths.pocket_json)
            validate_saved_monitor(monitor, args)
            validate_monitor_topology(monitor, modeller.topology)
        commit_restart_settings(args, paths)
    else:
        run_equilibration(
            simulation,
            system,
            modeller,
            args,
            temperature,
            pressure,
            paths,
            equilibration_steps,
            equilibration_report_steps,
        )
        monitor = None
        if args.early_stop:
            components = load_components(paths.components_json)
            production_state = simulation.context.getState(getPositions=True)
            monitor = initialize_monitor(
                modeller.topology,
                production_state.getPositions(),
                production_state.getPeriodicBoxVectors(),
                components,
                args,
            )
            write_pocket(paths.pocket_json, monitor)
            write_status(
                paths.status_json,
                monitor_status(
                    "running",
                    monitor,
                    args.production_ns,
                    0.0,
                    0,
                    0,
                ),
            )
    run_production(
        simulation,
        modeller,
        args,
        paths,
        production_target_steps,
        production_report_steps,
        checkpoint_steps,
        append=resume,
        monitor=monitor,
        monitor_steps=monitor_steps,
        performance_steps=performance_steps,
    )
    print(f"Finished. Results written to {paths.workdir}")
