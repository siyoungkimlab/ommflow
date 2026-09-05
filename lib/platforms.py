"""Platform selection, precision, and the banner describing a run."""

from __future__ import annotations

import openmm as mm
from openmm import app, unit

from ommflow.lib.ligands import STANDARD_SOLVENT_AND_ION_RESIDUES


PRECISIONS = ("mixed", "single", "double")
_PRECISION_PROPERTY = "Precision"


def _candidate_platforms(platform_name: str | None) -> list[mm.Platform]:
    """Return the platforms to try, fastest first, or just the requested one."""
    if platform_name:
        return [mm.Platform.getPlatformByName(platform_name)]
    platforms = [
        mm.Platform.getPlatform(index)
        for index in range(mm.Platform.getNumPlatforms())
    ]
    return sorted(platforms, key=lambda platform: platform.getSpeed(), reverse=True)


def create_simulation(
    topology: app.Topology,
    system: mm.System,
    integrator: object,
    platform_name: str | None,
    precision: str,
    precision_specified: bool = False,
) -> tuple[app.Simulation, list[str]]:
    """Create a Simulation on the best usable platform at the wanted precision.

    Not every platform can honor every precision: Apple's OpenCL, for one,
    rejects ``mixed`` because it lacks the double support mixed accumulation
    needs. An explicitly requested precision is therefore an error when it
    cannot be met, while the default quietly falls back to the platform's own
    precision and says so in the returned notes.
    """
    notes: list[str] = []
    failures: list[str] = []
    for platform in _candidate_platforms(platform_name):
        name = platform.getName()
        if _PRECISION_PROPERTY not in platform.getPropertyNames():
            if precision_specified:
                notes.append(
                    f"The {name} platform has no {_PRECISION_PROPERTY} setting; "
                    f"--precision {precision} does not apply to it."
                )
            simulation = _try_create(topology, system, integrator, platform, {}, failures)
            if simulation is not None:
                return simulation, notes
            continue

        simulation = _try_create(
            topology,
            system,
            integrator,
            platform,
            {_PRECISION_PROPERTY: precision},
            failures,
        )
        if simulation is not None:
            return simulation, notes
        if precision_specified:
            raise RuntimeError(
                f"The {name} platform cannot run at --precision {precision}: "
                + failures[-1]
            )
        simulation = _try_create(topology, system, integrator, platform, {}, failures)
        if simulation is not None:
            notes.append(
                f"The {name} platform rejected {_PRECISION_PROPERTY}={precision}; "
                "using its own default precision instead. Pass --precision "
                "explicitly to make this an error."
            )
            return simulation, notes

    raise RuntimeError(
        "Could not create a simulation on any available platform:\n  "
        + "\n  ".join(failures)
    )


def _try_create(
    topology: app.Topology,
    system: mm.System,
    integrator: object,
    platform: mm.Platform,
    properties: dict[str, str],
    failures: list[str],
) -> app.Simulation | None:
    try:
        return app.Simulation(topology, system, integrator, platform, properties)
    except Exception as error:  # noqa: BLE001 - platforms raise assorted types
        described = ", ".join(f"{k}={v}" for k, v in properties.items()) or "defaults"
        failures.append(f"{platform.getName()} ({described}): {error}")
        return None


def describe_platform(context: mm.Context) -> str:
    """Report the platform a context runs on and every property it exposes."""
    platform = context.getPlatform()
    lines = [f"Platform: {platform.getName()} (speed {platform.getSpeed():g})"]
    for name in platform.getPropertyNames():
        try:
            value = platform.getPropertyValue(context, name)
        except Exception:  # noqa: BLE001 - a property may be unreadable
            value = "<unavailable>"
        lines.append(f"  {name}: {value}")
    return "\n".join(lines)


_WATER_RESIDUES = frozenset(
    {"HOH", "WAT", "SOL", "SPC", "SPCE", "TIP3", "TIP4", "TIP5", "TP3", "TP4",
     "TP5", "OPC", "OPC3"}
)


def describe_system(topology: app.Topology, system: mm.System) -> str:
    """Report the size and shape of the system a run will integrate."""
    residues = list(topology.residues())
    waters = sum(1 for residue in residues if residue.name.upper() in _WATER_RESIDUES)
    ions = sum(
        1
        for residue in residues
        if residue.name.upper() in STANDARD_SOLVENT_AND_ION_RESIDUES
        and residue.name.upper() not in _WATER_RESIDUES
    )
    virtual_sites = sum(
        1 for index in range(system.getNumParticles()) if system.isVirtualSite(index)
    )
    lines = [
        f"System: {system.getNumParticles()} particles "
        f"({topology.getNumAtoms()} topology atoms, {virtual_sites} virtual sites)",
        f"  chains: {topology.getNumChains()}, residues: {len(residues)} "
        f"({waters} water, {ions} ion), bonds: {topology.getNumBonds()}",
        f"  constraints: {system.getNumConstraints()}, "
        f"forces: {', '.join(force.__class__.__name__ for force in system.getForces())}",
    ]
    box = system.getDefaultPeriodicBoxVectors()
    if box is not None:
        lengths = [
            vector[axis].value_in_unit(unit.nanometer)
            for axis, vector in enumerate(box)
        ]
        volume = lengths[0] * lengths[1] * lengths[2]
        lines.append(
            f"  box: {lengths[0]:.2f} x {lengths[1]:.2f} x {lengths[2]:.2f} nm "
            f"({volume:.1f} nm^3)"
        )
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            method = {
                mm.NonbondedForce.NoCutoff: "NoCutoff",
                mm.NonbondedForce.CutoffNonPeriodic: "CutoffNonPeriodic",
                mm.NonbondedForce.CutoffPeriodic: "CutoffPeriodic",
                mm.NonbondedForce.Ewald: "Ewald",
                mm.NonbondedForce.PME: "PME",
                mm.NonbondedForce.LJPME: "LJPME",
            }.get(force.getNonbondedMethod(), str(force.getNonbondedMethod()))
            cutoff = force.getCutoffDistance().value_in_unit(unit.nanometer)
            lines.append(f"  nonbonded: {method}, cutoff {cutoff:g} nm")
            break
    return "\n".join(lines)
