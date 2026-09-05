"""Configuration and interval tests for opt-in early stopping."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from openmm import unit

from ommflow.lib.config import DEFAULTS, load_configuration, parse_arguments
from ommflow.lib.reporting import validate_monitor_interval


def test_early_stop_defaults_are_disabled_and_cli_is_opt_in() -> None:
    defaults = parse_arguments([])
    assert defaults.early_stop is False
    assert defaults.monitor_ligand is None
    assert defaults.monitor_chain is None
    assert defaults.monitor_component is None
    assert defaults.monitor_interval_ns == DEFAULTS["monitor_interval_ns"]

    enabled = parse_arguments(
        [
            "--early-stop",
            "--monitor-ligand",
            "ligand-1",
            "--confirmation-checks",
            "3",
        ]
    )
    assert enabled.early_stop is True
    assert enabled.monitor_ligand == "ligand-1"
    assert enabled.confirmation_checks == 3
    assert enabled.early_stop_specified is True

    chain_target = parse_arguments(["--early-stop", "--monitor-chain", "B"])
    assert chain_target.monitor_chain == "B"
    component_target = parse_arguments(
        ["--early-stop", "--monitor-component", "component-2"]
    )
    assert component_target.monitor_component == "component-2"
    with pytest.raises(SystemExit):
        parse_arguments(
            ["--monitor-ligand", "ligand-0", "--monitor-component", "component-1"]
        )


def test_monitor_config_rejects_invalid_range_and_type() -> None:
    artifact_dir = Path(__file__).parent / ".monitor-test-artifacts"
    artifact_dir.mkdir(exist_ok=True)
    config_path = artifact_dir / "invalid.toml"
    try:
        config_path.write_text("monitor_interval_ns = 0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="must be positive"):
            load_configuration(config_path)
        with pytest.raises(SystemExit):
            parse_arguments(["--confirmation-checks", "0"])
        config_path.write_text(
            'monitor_chain = "A"\nmonitor_component = "component-0"\n',
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="at most one monitor target"):
            load_configuration(config_path)
    finally:
        config_path.unlink(missing_ok=True)
        artifact_dir.rmdir()


def test_monitor_interval_must_align_to_checkpoint_steps() -> None:
    args = SimpleNamespace(early_stop=True, monitor_interval_ns=0.1)
    timestep = 2 * unit.femtoseconds
    assert validate_monitor_interval(args, timestep, 5000) == 50000

    args.monitor_interval_ns = 0.015
    with pytest.raises(ValueError, match="integer multiple"):
        validate_monitor_interval(args, timestep, 5000)
