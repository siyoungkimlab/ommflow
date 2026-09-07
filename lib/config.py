"""Configuration loading, migration, serialization, and command-line parsing."""

from __future__ import annotations

import argparse
import math
from numbers import Real
from pathlib import Path
import re
import sys
import tomllib

from ommflow.lib.platforms import PRECISIONS
from ommflow.lib.restraints import DIHEDRAL_RESTRAINTS
from ommflow.lib.forcefields import (
    FORCE_FIELD_FAMILIES,
    WATER_MODELS,
    compatible_water_models,
    default_cutoff_nm,
)


PLATFORMS = ("CUDA", "OpenCL", "Metal", "CPU", "Reference")
# Hydrogen mass and timestep used when hydrogen mass repartitioning is on.
HYDROGEN_MASS_AMU = 4.0
HMR_INTEGRATION_FS = 4.0
LIGAND_MODES = ("disabled", "auto")
LIGAND_FORCE_FIELDS = ("gaff-2.11",)
DEFAULTS = {
    "workdir": Path("openmm_md"),
    "padding_nm": 1.0,
    "cutoff_nm": None,
    "saltM": 0.15,
    "proteinff": "amber19sb",
    "waterff": "opc",
    "ligand_mode": "auto",
    "ligandff": "gaff-2.11",
    "temperature": 298.0,
    "pressure": 1.0,
    "equilibration_ns": 0.1,
    "equilibration_report_interval_ns": 0.01,
    "production_ns": 100.0,
    "production_report_interval_ns": 1.0,
    "checkpoint_interval_ns": 0.01,
    "performance_interval_ns": 1.0,
    "integration_fs": 2.0,
    "hmr": False,
    "dihedral_restraint": "none",
    "dihedral_restraint_kJ": 20.0,
    "seed": 0,
    "precision": "mixed",
    "early_stop": False,
    "monitor_ligand": None,
    "monitor_chain": None,
    "monitor_component": None,
    "monitor_interval_ns": 0.1,
    "pocket_cutoff_nm": 0.5,
    "contact_cutoff_nm": 0.5,
    "detach_cutoff_nm": 0.8,
    "confirmation_checks": 2,
}
CONFIGURATION_KEYS = {
    "input_structure",
    *DEFAULTS,
    "platform",
}
_NUMERIC_KEYS = {
    "padding_nm",
    "cutoff_nm",
    "saltM",
    "temperature",
    "pressure",
    "equilibration_ns",
    "production_ns",
    "equilibration_report_interval_ns",
    "production_report_interval_ns",
    "checkpoint_interval_ns",
    "performance_interval_ns",
    "integration_fs",
    "dihedral_restraint_kJ",
    "monitor_interval_ns",
    "pocket_cutoff_nm",
    "contact_cutoff_nm",
    "detach_cutoff_nm",
}

# Settings a restart restores from workdir/final.toml when the user does not
# supply them again. Each entry maps the setting to the option that overrides it.
RESTARTABLE_SETTINGS = (
    ("production_ns", "--production-ns"),
    ("integration_fs", "--integration-fs"),
    ("hmr", "--hmr"),
    ("dihedral_restraint", "--dihedral-restraint"),
    ("dihedral_restraint_kJ", "--dihedral-restraint-kJ"),
    ("equilibration_ns", "--equilibration-ns"),
    ("equilibration_report_interval_ns", "--equilibration-report-interval-ns"),
    ("production_report_interval_ns", "--production-report-interval-ns"),
    ("checkpoint_interval_ns", "--checkpoint-interval-ns"),
    ("performance_interval_ns", "--performance-interval-ns"),
    ("platform", "--platform"),
    ("precision", "--precision"),
    ("early_stop", "--early-stop"),
    ("monitor_ligand", "--monitor-ligand"),
    ("monitor_chain", "--monitor-chain"),
    ("monitor_component", "--monitor-component"),
    ("monitor_interval_ns", "--monitor-interval-ns"),
    ("pocket_cutoff_nm", "--pocket-cutoff-nm"),
    ("contact_cutoff_nm", "--contact-cutoff-nm"),
    ("detach_cutoff_nm", "--detach-cutoff-nm"),
    ("confirmation_checks", "--confirmation-checks"),
)
MONITOR_TARGET_SELECTORS = ("monitor_ligand", "monitor_chain", "monitor_component")

_DEFAULT_CONFIGURATION = """\
# OpenMM protein MD settings. Command-line values override this file.
input_structure = "protein.pdb"
workdir = "openmm_md"

# Protein force field: amber14sb, amber15ipq, amber19sb, charmm36, or charmm36_2024
proteinff = "amber19sb"

# Amber water/ion parameter set: opc, opc3, spce, tip3p, tip3pfb, tip4pew, or tip4pfb
# CHARMM water/ion parameter set: tip3p, tip3p-pme-b, tip3p-pme-f, spce,
# tip4p2005, tip4pew, tip5p, or tip5pew
waterff = "opc"

# Automatically parameterize eligible disconnected ligands with GAFF 2.11.
# Ligand candidates require DMS or MAE chemistry data and Amber protein parameters.
ligand_mode = "auto"
ligandff = "gaff-2.11"

padding_nm = 1.0
# cutoff_nm defaults to 0.9 for Amber and 1.2 for CHARMM.
# cutoff_nm = 0.9
saltM = 0.15
temperature = 298.0
pressure = 1.0
equilibration_ns = 0.1
production_ns = 100.0
equilibration_report_interval_ns = 0.01
production_report_interval_ns = 1.0
checkpoint_interval_ns = 0.01
# Timing breakdown written to performance.csv.
performance_interval_ns = 1.0
integration_fs = 2.0
# Hydrogen mass repartitioning. When true, integration_fs defaults to 4.0.
hmr = false

# Restrain protein backbone phi/psi to the input structure: none, bb, or ss.
# ss restrains only residues in helices and sheets, and needs MDTraj.
dihedral_restraint = "none"
dihedral_restraint_kJ = 20.0
seed = 0

# GPU floating-point precision: mixed (default), single, or double.
# Left commented so the default can fall back on a platform that cannot honor
# it; setting it here makes an unsupported precision an error.
# precision = "mixed"

# Disabled by default.  When enabled, stop production after a selected target
# has no pocket contacts and is farther than detach_cutoff_nm repeatedly.
early_stop = false
# Select at most one target. ligand-N remains supported for GAFF ligands.
# monitor_ligand = "ligand-0"
# monitor_chain = "B"
# monitor_component = "component-2"
monitor_interval_ns = 0.1
pocket_cutoff_nm = 0.5
contact_cutoff_nm = 0.5
detach_cutoff_nm = 0.8
confirmation_checks = 2

# Leave this absent to have OpenMM choose automatically.
# platform = "Metal"
"""


def load_configuration(path: Path) -> dict:
    """Load and validate TOML settings."""
    with path.open("rb") as handle:
        configuration = tomllib.load(handle)

    unknown_keys = configuration.keys() - CONFIGURATION_KEYS
    if unknown_keys:
        raise ValueError(
            "Unknown configuration setting(s): " + ", ".join(sorted(unknown_keys))
        )
    for key in ("input_structure", "workdir"):
        if key in configuration and not isinstance(configuration[key], str):
            raise ValueError(f"Configuration setting '{key}' must be a string path.")
    for key in ("ligand_mode", "ligandff", "precision", "dihedral_restraint"):
        if key in configuration and not isinstance(configuration[key], str):
            raise ValueError(f"Configuration setting '{key}' must be a string.")
    for key in MONITOR_TARGET_SELECTORS:
        if key in configuration and (
            not isinstance(configuration[key], str) or not configuration[key]
        ):
            raise ValueError(
                f"Configuration setting '{key}' must be a non-empty string."
            )
    _validate_target_selectors(configuration)
    for key in ("early_stop", "hmr"):
        if key in configuration and not isinstance(configuration[key], bool):
            raise ValueError(f"Configuration setting '{key}' must be a boolean.")
    for key in _NUMERIC_KEYS:
        if key in configuration and (
            isinstance(configuration[key], bool)
            or not isinstance(configuration[key], Real)
        ):
            raise ValueError(f"Configuration setting '{key}' must be a number.")
    if "seed" in configuration and (
        isinstance(configuration["seed"], bool)
        or not isinstance(configuration["seed"], int)
    ):
        raise ValueError("Configuration setting 'seed' must be an integer.")
    if "confirmation_checks" in configuration and (
        isinstance(configuration["confirmation_checks"], bool)
        or not isinstance(configuration["confirmation_checks"], int)
    ):
        raise ValueError(
            "Configuration setting 'confirmation_checks' must be an integer."
        )
    for key in (
        "dihedral_restraint_kJ",
        "performance_interval_ns",
        "monitor_interval_ns",
        "pocket_cutoff_nm",
        "contact_cutoff_nm",
        "detach_cutoff_nm",
    ):
        if key in configuration and (
            not math.isfinite(configuration[key]) or configuration[key] <= 0
        ):
            raise ValueError(f"Configuration setting '{key}' must be positive.")
    if (
        "confirmation_checks" in configuration
        and configuration["confirmation_checks"] < 1
    ):
        raise ValueError(
            "Configuration setting 'confirmation_checks' must be at least one."
        )
    if (
        "proteinff" in configuration
        and configuration["proteinff"] not in FORCE_FIELD_FAMILIES
    ):
        raise ValueError(
            "'proteinff' must be one of: " + ", ".join(FORCE_FIELD_FAMILIES)
        )
    if "waterff" in configuration and configuration["waterff"] not in WATER_MODELS:
        raise ValueError("'waterff' must be one of: " + ", ".join(WATER_MODELS))
    if (
        "ligand_mode" in configuration
        and configuration["ligand_mode"] not in LIGAND_MODES
    ):
        raise ValueError(
            "'ligand_mode' must be one of: " + ", ".join(LIGAND_MODES)
        )
    if (
        "ligandff" in configuration
        and configuration["ligandff"] not in LIGAND_FORCE_FIELDS
    ):
        raise ValueError(
            "'ligandff' must be one of: " + ", ".join(LIGAND_FORCE_FIELDS)
        )
    if (
        "dihedral_restraint" in configuration
        and configuration["dihedral_restraint"] not in DIHEDRAL_RESTRAINTS
    ):
        raise ValueError(
            "'dihedral_restraint' must be one of: " + ", ".join(DIHEDRAL_RESTRAINTS)
        )
    if "precision" in configuration and configuration["precision"] not in PRECISIONS:
        raise ValueError("'precision' must be one of: " + ", ".join(PRECISIONS))
    if "platform" in configuration and configuration["platform"] not in PLATFORMS:
        raise ValueError("'platform' must be one of: " + ", ".join(PLATFORMS))
    return configuration


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser with the workflow's documented defaults."""
    parser = argparse.ArgumentParser(
        description="Run explicit-solvent protein MD with OpenMM."
    )
    parser.add_argument(
        "input_structure",
        type=Path,
        nargs="?",
        metavar="INPUT_STRUCTURE",
        help="Hydrogen-complete protein structure in PDB, DMS, MAE, or GRO format.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="TOML settings file. Explicit CLI values override its settings.",
    )
    parser.add_argument(
        "--write-default-config",
        type=Path,
        metavar="FILE",
        help="Write a default TOML settings template to FILE and exit.",
    )
    parser.add_argument(
        "--list-components",
        action="store_true",
        help=(
            "Print the components of INPUT_STRUCTURE with the selector that "
            "names each one, then exit without building or running anything."
        ),
    )
    parser.add_argument(
        "--workdir", type=Path, default=DEFAULTS["workdir"], help="Output directory."
    )
    parser.add_argument(
        "--padding-nm",
        type=float,
        default=DEFAULTS["padding_nm"],
        help="Minimum protein-to-box-edge distance in nm (default: 1.0).",
    )
    parser.add_argument(
        "--cutoff-nm",
        type=float,
        default=DEFAULTS["cutoff_nm"],
        help="Nonbonded cutoff in nm (default: 0.9 Amber, 1.2 CHARMM).",
    )
    parser.add_argument(
        "--saltM",
        dest="saltM",
        type=float,
        default=DEFAULTS["saltM"],
        help="Added NaCl concentration in molar (default: 0.15).",
    )
    parser.add_argument(
        "--proteinff",
        choices=FORCE_FIELD_FAMILIES,
        default=DEFAULTS["proteinff"],
        help="Protein force field (default: amber19sb).",
    )
    parser.add_argument(
        "--waterff",
        choices=WATER_MODELS,
        default=DEFAULTS["waterff"],
        help="Water and ion parameter set (default: opc).",
    )
    parser.add_argument(
        "--ligand-mode",
        choices=LIGAND_MODES,
        default=DEFAULTS["ligand_mode"],
        help=(
            "Ligand handling: auto (default) or disabled for strict standard "
            "force-field-only preparation."
        ),
    )
    parser.add_argument(
        "--ligandff",
        choices=LIGAND_FORCE_FIELDS,
        default=DEFAULTS["ligandff"],
        help="Automatic ligand force field (default: gaff-2.11).",
    )
    parser.add_argument(
        "--temperature",
        dest="temperature",
        type=float,
        default=DEFAULTS["temperature"],
        help="Temperature in K (default: 298).",
    )
    parser.add_argument(
        "--pressure",
        dest="pressure",
        type=float,
        default=DEFAULTS["pressure"],
        help="Pressure in bar (default: 1).",
    )
    parser.add_argument(
        "--equilibration-ns",
        type=float,
        default=DEFAULTS["equilibration_ns"],
        help="NVT and NPT equilibration duration each in ns (default: 0.1).",
    )
    parser.add_argument(
        "--equilibration-report-interval-ns",
        type=float,
        default=DEFAULTS["equilibration_report_interval_ns"],
        help="Equilibration reporting interval in ns (default: 0.01).",
    )
    parser.add_argument(
        "--production-ns",
        type=float,
        default=DEFAULTS["production_ns"],
        help="Target total NPT production time in ns (default: 100).",
    )
    parser.add_argument(
        "--production-report-interval-ns",
        type=float,
        default=DEFAULTS["production_report_interval_ns"],
        help="Production reporting interval in ns (default: 1).",
    )
    parser.add_argument(
        "--checkpoint-interval-ns",
        type=float,
        default=DEFAULTS["checkpoint_interval_ns"],
        help="Checkpoint writing interval in ns (default: 0.01).",
    )
    parser.add_argument(
        "--performance-interval-ns",
        type=float,
        default=DEFAULTS["performance_interval_ns"],
        help="Timing-breakdown interval in ns for performance.csv (default: 1).",
    )
    parser.add_argument(
        "--integration-fs",
        type=float,
        default=DEFAULTS["integration_fs"],
        help="Integration timestep in fs (default: 2.0).",
    )
    parser.add_argument(
        "--hmr",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["hmr"],
        help=(
            "Hydrogen mass repartitioning (default: off). Raises bonded "
            "hydrogens to 4 amu, taking the mass from their heavy atom, and "
            "makes --integration-fs default to 4 instead of 2. Rigid water is "
            "left untouched for every water model."
        ),
    )
    parser.add_argument("--seed", type=int, default=DEFAULTS["seed"], help="Random seed.")
    parser.add_argument(
        "--early-stop",
        action=argparse.BooleanOptionalAction,
        default=DEFAULTS["early_stop"],
        help=(
            "Stop production after confirmed target detachment "
            "(default: disabled)."
        ),
    )
    parser.add_argument(
        "--monitor-ligand",
        metavar="LIGAND_ID",
        help="Exact GAFF ligand-N ID from components.json to monitor.",
    )
    parser.add_argument(
        "--monitor-chain",
        metavar="CHAIN",
        help="Input-structure chain ID whose one disconnected component to monitor.",
    )
    parser.add_argument(
        "--monitor-component",
        metavar="COMPONENT_ID",
        help="Exact component-N ID from components.json to monitor.",
    )
    parser.add_argument(
        "--monitor-interval-ns",
        type=float,
        default=DEFAULTS["monitor_interval_ns"],
        help="Target-monitoring interval in ns (default: 0.1).",
    )
    parser.add_argument(
        "--pocket-cutoff-nm",
        type=float,
        default=DEFAULTS["pocket_cutoff_nm"],
        help="Initial target-to-protein pocket cutoff in nm (default: 0.5).",
    )
    parser.add_argument(
        "--contact-cutoff-nm",
        type=float,
        default=DEFAULTS["contact_cutoff_nm"],
        help="Target-pocket heavy-atom contact cutoff in nm (default: 0.5).",
    )
    parser.add_argument(
        "--detach-cutoff-nm",
        type=float,
        default=DEFAULTS["detach_cutoff_nm"],
        help="Minimum target-pocket distance for detachment in nm (default: 0.8).",
    )
    parser.add_argument(
        "--confirmation-checks",
        type=int,
        default=DEFAULTS["confirmation_checks"],
        help="Consecutive detached checks required to stop (default: 2).",
    )
    parser.add_argument(
        "--dihedral-restraint",
        choices=DIHEDRAL_RESTRAINTS,
        default=DEFAULTS["dihedral_restraint"],
        help=(
            "Restrain protein backbone phi and psi to the input structure: "
            "bb for the whole backbone, ss for residues in helices and sheets "
            "only (needs MDTraj), none to disable (default)."
        ),
    )
    parser.add_argument(
        "--dihedral-restraint-kJ",
        dest="dihedral_restraint_kJ",
        type=float,
        default=DEFAULTS["dihedral_restraint_kJ"],
        help="Dihedral restraint strength in kJ/mol (default: 20).",
    )
    parser.add_argument(
        "--precision",
        choices=PRECISIONS,
        default=DEFAULTS["precision"],
        help=(
            "GPU floating-point precision (default: mixed). Applies only to "
            "platforms that expose it; a platform that cannot honor the default "
            "falls back to its own precision with a warning."
        ),
    )
    parser.add_argument(
        "--platform",
        choices=PLATFORMS,
        help="Optional OpenMM platform; defaults to OpenMM's automatic selection.",
    )
    return parser


def parse_arguments(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI options after applying TOML values as defaults."""
    raw_arguments = sys.argv[1:] if argv is None else argv
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=Path)
    config_args, _ = config_parser.parse_known_args(raw_arguments)

    parser = build_parser()
    configuration = {}
    if config_args.config:
        if not config_args.config.is_file():
            parser.error(f"Configuration file does not exist: {config_args.config}")
        try:
            configuration = load_configuration(config_args.config)
        except (OSError, tomllib.TOMLDecodeError, ValueError) as error:
            parser.error(f"Unable to read configuration file: {error}")
        parser.set_defaults(**configuration)

    args = parser.parse_args(raw_arguments)
    for setting, option in RESTARTABLE_SETTINGS:
        setattr(
            args,
            f"{setting}_specified",
            setting in configuration
            or any(
                argument == option
                or argument.startswith(f"{option}=")
                or argument == option.replace("--", "--no-", 1)
                for argument in raw_arguments
            ),
        )
    if args.input_structure is not None:
        args.input_structure = Path(args.input_structure)
    args.workdir = Path(args.workdir)
    compatible_models = compatible_water_models(args.proteinff)
    if args.waterff not in compatible_models:
        parser.error(
            f"Water model '{args.waterff}' is incompatible with {args.proteinff}. "
            f"Choose one of: {', '.join(compatible_models)}"
        )
    if args.cutoff_nm is None:
        args.cutoff_nm = default_cutoff_nm(args.proteinff)
    if args.hmr and not args.integration_fs_specified:
        # Repartitioned hydrogens are the reason a 4 fs step is stable.
        args.integration_fs = HMR_INTEGRATION_FS
    _validate_early_stop_settings(args, parser)
    return args


def _validate_early_stop_settings(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    """Reject invalid monitor settings before an expensive workflow starts."""
    numeric_settings = (
        ("--monitor-interval-ns", args.monitor_interval_ns),
        ("--pocket-cutoff-nm", args.pocket_cutoff_nm),
        ("--contact-cutoff-nm", args.contact_cutoff_nm),
        ("--detach-cutoff-nm", args.detach_cutoff_nm),
    )
    for name, value in numeric_settings:
        if not math.isfinite(value) or value <= 0:
            parser.error(f"{name} must be positive.")
    if args.confirmation_checks < 1:
        parser.error("--confirmation-checks must be at least one.")
    for option, value in (
        ("--monitor-ligand", args.monitor_ligand),
        ("--monitor-chain", args.monitor_chain),
        ("--monitor-component", args.monitor_component),
    ):
        if value == "":
            parser.error(f"{option} must be a non-empty value.")
    selected = [
        option
        for option, value in (
            ("--monitor-ligand", args.monitor_ligand),
            ("--monitor-chain", args.monitor_chain),
            ("--monitor-component", args.monitor_component),
        )
        if value is not None
    ]
    if len(selected) > 1:
        parser.error(
            "Specify at most one monitor target selector: " + ", ".join(selected)
        )


def require_early_stop_for_selectors(args) -> None:
    """Reject an explicitly chosen target that early stopping would ignore."""
    if getattr(args, "early_stop", False):
        return
    selected = [
        option
        for setting, option in RESTARTABLE_SETTINGS
        if setting in MONITOR_TARGET_SELECTORS
        and getattr(args, setting, None) is not None
        and getattr(args, f"{setting}_specified", False)
    ]
    if selected:
        raise ValueError(
            ", ".join(selected)
            + " selects an early-stop target, but early stopping is disabled. "
            "Add --early-stop, or drop the selector."
        )


def toml_string(value: object) -> str:
    """Quote a TOML string value."""
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_final_configuration(args: argparse.Namespace, output_dir: Path) -> None:
    """Write the resolved settings used by a new run."""
    settings = [
        ("input_structure", toml_string(args.input_structure)),
        ("workdir", toml_string(output_dir)),
        ("proteinff", toml_string(args.proteinff)),
        ("waterff", toml_string(args.waterff)),
        ("ligand_mode", toml_string(args.ligand_mode)),
        ("ligandff", toml_string(args.ligandff)),
        ("padding_nm", str(args.padding_nm)),
        ("cutoff_nm", str(args.cutoff_nm)),
        ("saltM", str(args.saltM)),
        ("temperature", str(args.temperature)),
        ("pressure", str(args.pressure)),
        ("equilibration_ns", str(args.equilibration_ns)),
        ("production_ns", str(args.production_ns)),
        (
            "equilibration_report_interval_ns",
            str(args.equilibration_report_interval_ns),
        ),
        ("production_report_interval_ns", str(args.production_report_interval_ns)),
        ("checkpoint_interval_ns", str(args.checkpoint_interval_ns)),
        ("performance_interval_ns", str(args.performance_interval_ns)),
        ("integration_fs", str(args.integration_fs)),
        ("hmr", str(args.hmr).lower()),
        ("dihedral_restraint", toml_string(args.dihedral_restraint)),
        ("dihedral_restraint_kJ", str(args.dihedral_restraint_kJ)),
        ("seed", str(args.seed)),
        ("precision", toml_string(args.precision)),
        ("early_stop", str(args.early_stop).lower()),
        ("monitor_interval_ns", str(args.monitor_interval_ns)),
        ("pocket_cutoff_nm", str(args.pocket_cutoff_nm)),
        ("contact_cutoff_nm", str(args.contact_cutoff_nm)),
        ("detach_cutoff_nm", str(args.detach_cutoff_nm)),
        ("confirmation_checks", str(args.confirmation_checks)),
    ]
    if args.monitor_ligand is not None:
        settings.append(("monitor_ligand", toml_string(args.monitor_ligand)))
    if args.monitor_chain is not None:
        settings.append(("monitor_chain", toml_string(args.monitor_chain)))
    if args.monitor_component is not None:
        settings.append(("monitor_component", toml_string(args.monitor_component)))
    if args.platform:
        settings.append(("platform", toml_string(args.platform)))
    contents = (
        "# Resolved OpenMM MD settings used for this run.\n"
        "# Values from the command line override the input TOML file.\n"
        + "".join(f"{key} = {value}\n" for key, value in settings)
    )
    (output_dir / "final.toml").write_text(contents, encoding="utf-8")


def write_default_configuration(path: Path) -> None:
    """Write the documented starter configuration."""
    path.write_text(_DEFAULT_CONFIGURATION, encoding="utf-8")


def update_production_target(path: Path, production_ns: float) -> None:
    """Replace the one saved total production target."""
    contents = path.read_text(encoding="utf-8")
    updated_contents, replacements = re.subn(
        r"^production_ns = .*$",
        f"production_ns = {production_ns}",
        contents,
        flags=re.MULTILINE,
    )
    if replacements != 1:
        raise ValueError(f"Expected exactly one production_ns setting in {path}")
    path.write_text(updated_contents, encoding="utf-8")


def _toml_value(value: object) -> str:
    """Render one resolved setting as a TOML scalar."""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (str, Path)):
        return toml_string(value)
    return str(value)


def update_restart_settings(path: Path, args: argparse.Namespace) -> None:
    """Persist every resolved restart setting so the next resume matches it.

    ``production_ns`` and ``integration_fs`` are excluded: the first is written
    by :func:`update_production_target` and the second can never change once a
    serialized integrator exists.
    """
    contents = path.read_text(encoding="utf-8")
    for setting, _ in RESTARTABLE_SETTINGS:
        if setting in {"production_ns", "integration_fs"}:
            continue
        value = getattr(args, setting, None)
        if value is None:
            contents = re.sub(
                rf"^{re.escape(setting)} = .*\n?", "", contents, flags=re.MULTILINE
            )
            continue
        replacement = f"{setting} = {_toml_value(value)}"
        contents, replacements = re.subn(
            rf"^{re.escape(setting)} = .*$",
            lambda _, text=replacement: text,
            contents,
            flags=re.MULTILINE,
        )
        if replacements == 0:
            contents += replacement + "\n"
        elif replacements != 1:
            raise ValueError(f"Expected at most one {setting} setting in {path}")
    path.write_text(contents, encoding="utf-8")


def _validate_target_selectors(settings: dict) -> None:
    """Require TOML to identify no more than one early-stop target."""
    selected = [key for key in MONITOR_TARGET_SELECTORS if settings.get(key) is not None]
    if len(selected) > 1:
        raise ValueError(
            "Specify at most one monitor target setting: " + ", ".join(selected)
        )
