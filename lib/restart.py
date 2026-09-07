"""Run artifacts and restart persistence for the OpenMM workflow."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import openmm as mm
from openmm import app

from ommflow.lib.config import (
    RESTARTABLE_SETTINGS,
    load_configuration,
    update_production_target,
    update_restart_settings,
    write_final_configuration,
)


@dataclass(frozen=True)
class RunPaths:
    """Paths to the durable artifacts in one workflow work directory."""

    workdir: Path

    @property
    def solvated_pdb(self) -> Path:
        return self.workdir / "solvated.pdb"

    @property
    def components_json(self) -> Path:
        return self.workdir / "components.json"

    @property
    def system_xml(self) -> Path:
        return self.workdir / "system.xml"

    @property
    def integrator_xml(self) -> Path:
        return self.workdir / "integrator.xml"

    @property
    def checkpoint(self) -> Path:
        return self.workdir / "checkpoint.chk"

    @property
    def final_configuration(self) -> Path:
        return self.workdir / "final.toml"

    @property
    def equilibrated_pdb(self) -> Path:
        return self.workdir / "equilibrated.pdb"

    @property
    def final_pdb(self) -> Path:
        return self.workdir / "final.pdb"

    @property
    def equilibration_dcd(self) -> Path:
        return self.workdir / "equilibration.dcd"

    @property
    def equilibration_csv(self) -> Path:
        return self.workdir / "equilibration.csv"

    @property
    def trajectory_dcd(self) -> Path:
        return self.workdir / "trajectory.dcd"

    @property
    def state_csv(self) -> Path:
        return self.workdir / "state.csv"

    @property
    def pocket_json(self) -> Path:
        return self.workdir / "pocket.json"

    @property
    def dihedral_restraints_csv(self) -> Path:
        return self.workdir / "dihedral_restraints.csv"

    @property
    def dihedral_restraint_png(self) -> Path:
        return self.workdir / "dihedral_restraint.png"

    @property
    def performance_csv(self) -> Path:
        return self.workdir / "performance.csv"

    @property
    def monitor_csv(self) -> Path:
        return self.workdir / "monitor.csv"

    @property
    def status_json(self) -> Path:
        return self.workdir / "status.json"


def is_restart(workdir: Path) -> bool:
    """An existing work directory with contents denotes an attempted restart.

    A pre-created empty directory, which batch schedulers routinely make, is a
    new run rather than an unresumable one.
    """
    if not workdir.exists():
        return False
    if not workdir.is_dir():
        return True
    return any(workdir.iterdir())


def prepare_restart(
    args: argparse.Namespace, paths: RunPaths
) -> tuple[app.Modeller, mm.System, object]:
    """Restore the saved settings and artifacts a restart needs, writing nothing."""
    required_paths = (
        paths.solvated_pdb,
        paths.system_xml,
        paths.integrator_xml,
        paths.checkpoint,
        paths.final_configuration,
    )
    for path in required_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Cannot resume: missing required file: {path}")

    saved_configuration = load_configuration(paths.final_configuration)
    if not args.production_ns_specified and "production_ns" not in saved_configuration:
        raise ValueError(
            f"Cannot resume: missing production_ns in {paths.final_configuration}"
        )
    for setting, _ in RESTARTABLE_SETTINGS:
        if (
            not getattr(args, f"{setting}_specified", False)
            and setting in saved_configuration
        ):
            setattr(args, setting, saved_configuration[setting])

    structure = app.PDBFile(str(paths.solvated_pdb))
    modeller = app.Modeller(structure.topology, structure.positions)
    system = mm.XmlSerializer.deserialize(paths.system_xml.read_text(encoding="utf-8"))
    integrator = mm.XmlSerializer.deserialize(
        paths.integrator_xml.read_text(encoding="utf-8")
    )
    return modeller, system, integrator


def commit_restart_settings(args: argparse.Namespace, paths: RunPaths) -> None:
    """Record the resolved restart settings once every check has passed.

    Writing only after validation keeps a rejected resume from leaving
    ``final.toml`` in a state that no later resume can load.
    """
    if args.production_ns_specified:
        update_production_target(paths.final_configuration, args.production_ns)
    update_restart_settings(paths.final_configuration, args)


def save_new_configuration(args: argparse.Namespace, paths: RunPaths) -> None:
    """Persist the resolved run settings before equilibration begins."""
    write_final_configuration(args, paths.workdir)


def save_equilibrated_artifacts(
    system: mm.System, integrator: object, paths: RunPaths
) -> None:
    """Persist the post-NPT system and integrator used by production."""
    paths.system_xml.write_text(mm.XmlSerializer.serialize(system), encoding="utf-8")
    paths.integrator_xml.write_text(
        mm.XmlSerializer.serialize(integrator), encoding="utf-8"
    )
