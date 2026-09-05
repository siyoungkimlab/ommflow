"""Tests for the settings a restart restores, and when it persists them."""

from __future__ import annotations

from pathlib import Path

import pytest

from ommflow.lib.config import (
    RESTARTABLE_SETTINGS,
    load_configuration,
    parse_arguments,
    require_early_stop_for_selectors,
    update_restart_settings,
    write_final_configuration,
)
from ommflow.lib.restart import RunPaths, prepare_restart


def _new_run_arguments(workdir: Path, *extra: str):
    args = parse_arguments(
        [
            "complex.pdb",
            "--workdir",
            str(workdir),
            "--production-report-interval-ns",
            "0.002",
            "--checkpoint-interval-ns",
            "0.002",
            "--platform",
            "CPU",
            *extra,
        ]
    )
    return args


def test_every_restartable_setting_reports_whether_it_was_specified() -> None:
    args = parse_arguments(["--checkpoint-interval-ns", "0.002"])
    for setting, _ in RESTARTABLE_SETTINGS:
        assert hasattr(args, f"{setting}_specified"), setting
    assert args.checkpoint_interval_ns_specified is True
    assert args.production_report_interval_ns_specified is False
    assert args.platform_specified is False


def test_restart_restores_reporting_intervals_and_platform(tmp_path: Path) -> None:
    workdir = tmp_path / "run"
    workdir.mkdir()
    write_final_configuration(
        _new_run_arguments(workdir, "--early-stop", "--monitor-chain", "B"), workdir
    )
    paths = RunPaths(workdir)
    for path in (paths.solvated_pdb, paths.system_xml, paths.integrator_xml, paths.checkpoint):
        path.write_text("not a real artifact", encoding="utf-8")
    original = paths.final_configuration.read_text(encoding="utf-8")

    resumed = parse_arguments(["--workdir", str(workdir), "--production-ns", "0.04"])
    assert resumed.checkpoint_interval_ns == 0.01  # command-line default
    with pytest.raises(Exception):
        prepare_restart(resumed, paths)

    assert resumed.checkpoint_interval_ns == 0.002
    assert resumed.production_report_interval_ns == 0.002
    assert resumed.platform == "CPU"
    assert resumed.early_stop is True
    assert resumed.monitor_chain == "B"
    # A resume that fails validation must leave final.toml loadable as it was.
    assert paths.final_configuration.read_text(encoding="utf-8") == original


def test_committing_restart_settings_replaces_only_the_chosen_selector(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "run"
    workdir.mkdir()
    args = _new_run_arguments(workdir, "--early-stop", "--monitor-chain", "B")
    write_final_configuration(args, workdir)
    final_toml = workdir / "final.toml"

    args.monitor_chain = None
    args.monitor_component = "component-2"
    args.confirmation_checks = 5
    update_restart_settings(final_toml, args)

    saved = load_configuration(final_toml)
    assert saved["monitor_component"] == "component-2"
    assert "monitor_chain" not in saved
    assert saved["confirmation_checks"] == 5
    assert saved["platform"] == "CPU"


def test_an_explicit_selector_requires_early_stop(tmp_path: Path) -> None:
    args = parse_arguments(["--monitor-chain", "B"])
    with pytest.raises(ValueError, match="--monitor-chain selects an early-stop target"):
        require_early_stop_for_selectors(args)

    # A selector restored from final.toml on a --no-early-stop resume is not a
    # user request and must not block the run.
    restored = parse_arguments(["--no-early-stop"])
    restored.monitor_chain = "B"
    require_early_stop_for_selectors(restored)

    require_early_stop_for_selectors(parse_arguments(["--early-stop", "--monitor-chain", "B"]))


def test_an_empty_work_directory_is_a_new_run(tmp_path: Path) -> None:
    from ommflow.lib.restart import is_restart

    empty = tmp_path / "scheduler-made"
    empty.mkdir()
    assert is_restart(empty) is False
    assert is_restart(tmp_path / "absent") is False

    (empty / "checkpoint.chk").write_text("", encoding="utf-8")
    assert is_restart(empty) is True


def test_the_installed_openmm_provides_every_documented_force_field() -> None:
    """The openmm>=8.5 floor exists so no family is half-shipped."""
    from ommflow.lib.forcefields import (
        FORCE_FIELD_FAMILIES,
        installed_force_field_families,
    )

    assert set(installed_force_field_families()) == set(FORCE_FIELD_FAMILIES)


def test_precision_defaults_to_mixed_and_is_restorable() -> None:
    from ommflow.lib.platforms import PRECISIONS

    assert parse_arguments([]).precision == "mixed"
    assert parse_arguments([]).precision_specified is False
    explicit = parse_arguments(["--precision", "single"])
    assert explicit.precision == "single"
    assert explicit.precision_specified is True
    assert ("precision", "--precision") in RESTARTABLE_SETTINGS
    assert set(PRECISIONS) == {"mixed", "single", "double"}


def test_precision_config_rejects_an_unknown_value(tmp_path: Path) -> None:
    config = tmp_path / "c.toml"
    config.write_text('precision = "quad"\n', encoding="utf-8")
    with pytest.raises(ValueError, match="'precision' must be one of"):
        load_configuration(config)
    config.write_text("precision = 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must be a string"):
        load_configuration(config)


def _documentation_sources() -> dict[str, str]:
    root = Path(__file__).resolve().parents[1]
    sources = {"README.md": (root / "README.md").read_text(encoding="utf-8")}
    for path in sorted((root / "docs").glob("*.rst")):
        sources[path.name] = path.read_text(encoding="utf-8")
    return sources


def test_every_setting_is_documented() -> None:
    """A new option is not finished until both doc surfaces mention it."""
    from ommflow.lib.config import CONFIGURATION_KEYS, build_parser

    sources = _documentation_sources()
    readme = sources["README.md"]
    reference = "\n".join(text for name, text in sources.items() if name != "README.md")

    undocumented = sorted(
        key
        for key in CONFIGURATION_KEYS
        if key not in readme or key not in reference
    )
    assert not undocumented, f"settings missing from the docs: {undocumented}"

    # Every command-line option must appear somewhere a reader can find it.
    options = {
        option
        for action in build_parser()._actions
        for option in action.option_strings
        if option.startswith("--") and not option.startswith("--no-")
    }
    everywhere = readme + reference
    missing = sorted(option for option in options if option not in everywhere)
    assert not missing, f"options missing from the docs: {missing}"


def test_the_configuration_examples_match_the_real_defaults() -> None:
    """Documented TOML must parse and must not contradict DEFAULTS."""
    import re
    import tomllib

    from ommflow.lib.config import CONFIGURATION_KEYS, DEFAULTS

    blocks: list[tuple[str, str]] = []
    for name, text in _documentation_sources().items():
        pattern = (
            r"```toml\n(.*?)```"
            if name.endswith(".md")
            else r"code-block:: toml\n\n(.*?)(?=\n\S|\Z)"
        )
        for block in re.findall(pattern, text, re.S):
            blocks.append((name, "\n".join(line.strip() for line in block.splitlines())))
    assert blocks

    for name, block in blocks:
        settings = tomllib.loads(block)
        for key, value in settings.items():
            assert key in CONFIGURATION_KEYS, f"{name}: unknown setting {key}"
            if key in DEFAULTS and key not in {"input_structure", "workdir"}:
                assert DEFAULTS[key] == value, (
                    f"{name}: documents {key} = {value!r} but the default is "
                    f"{DEFAULTS[key]!r}"
                )


def test_documented_version_and_metadata_stay_in_step() -> None:
    """conf.py pins a version by hand, so hold it to pyproject's."""
    import re
    import tomllib

    root = Path(__file__).resolve().parents[1]
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    version = project["project"]["version"]

    conf = (root / "docs" / "conf.py").read_text(encoding="utf-8")
    release = re.search(r'release\s*=\s*"([^"]+)"', conf)
    assert release is not None and release.group(1) == version, (
        f"docs/conf.py release {release and release.group(1)!r} != "
        f"pyproject version {version!r}"
    )

    # The license file and the packaging metadata must name the same license.
    assert project["project"]["license"] == "MIT"
    assert "MIT License" in (root / "LICENSE").read_text(encoding="utf-8")

    # ReadTheDocs must build on a Python the package actually supports.
    readthedocs = (root / ".readthedocs.yaml").read_text(encoding="utf-8")
    built_on = re.search(r'python:\s*"([\d.]+)"', readthedocs)
    minimum = project["project"]["requires-python"].lstrip(">=")
    assert built_on is not None
    assert tuple(map(int, built_on.group(1).split("."))) >= tuple(
        map(int, minimum.split("."))
    ), f"ReadTheDocs builds on {built_on.group(1)}, below requires-python {minimum}"
