"""Protein and water force-field compatibility mappings."""

from __future__ import annotations

from pathlib import Path

import openmm
from openmm import app


AMBER_WATER_MODELS = {
    "opc": ("opc.xml", "tip4pew"),
    "opc3": ("opc3.xml", "tip3p"),
    "spce": ("spce.xml", "spce"),
    "tip3p": ("tip3p.xml", "tip3p"),
    "tip3pfb": ("tip3pfb.xml", "tip3p"),
    "tip4pew": ("tip4pew.xml", "tip4pew"),
    "tip4pfb": ("tip4pfb.xml", "tip4pew"),
}
CHARMM_WATER_MODELS = {
    "tip3p": ("water.xml", "tip3p"),
    "tip3p-pme-b": ("tip3p-pme-b.xml", "tip3p"),
    "tip3p-pme-f": ("tip3p-pme-f.xml", "tip3p"),
    "spce": ("spce.xml", "spce"),
    "tip4p2005": ("tip4p2005.xml", "tip4pew"),
    "tip4pew": ("tip4pew.xml", "tip4pew"),
    "tip5p": ("tip5p.xml", "tip5p"),
    "tip5pew": ("tip5pew.xml", "tip5p"),
}
FORCE_FIELD_FAMILIES = {
    "amber14sb": ("amber14/protein.ff14SB.xml", "amber14", AMBER_WATER_MODELS),
    "amber15ipq": ("amber14/protein.ff15ipq.xml", "amber14", AMBER_WATER_MODELS),
    "amber19sb": ("amber19/protein.ff19SB.xml", "amber19", AMBER_WATER_MODELS),
    "charmm36": ("charmm36.xml", "charmm36", CHARMM_WATER_MODELS),
    "charmm36_2024": ("charmm36_2024.xml", "charmm36_2024", CHARMM_WATER_MODELS),
}
WATER_MODELS = tuple(sorted(AMBER_WATER_MODELS.keys() | CHARMM_WATER_MODELS.keys()))


def default_cutoff_nm(proteinff: str) -> float:
    """Return the conventional nonbonded cutoff for a force-field family."""
    if proteinff.startswith("amber"):
        return 0.9
    return 1.2


def compatible_water_models(proteinff: str) -> dict[str, tuple[str, str]]:
    """Return water models compatible with a protein force-field family."""
    return FORCE_FIELD_FAMILIES[proteinff][2]


def resolve_forcefields(proteinff: str, waterff: str) -> tuple[str, str, str, str]:
    """Resolve compatible OpenMM XML paths and Modeller solvent box model."""
    protein_xml, water_directory, water_models = FORCE_FIELD_FAMILIES[proteinff]
    try:
        water_xml, box_model = water_models[waterff]
    except KeyError as error:
        choices = ", ".join(water_models)
        raise ValueError(
            f"Water model '{waterff}' is incompatible with {proteinff}. "
            f"Choose one of: {choices}"
        ) from error
    return protein_xml, water_directory, water_xml, box_model


def _openmm_data_directory() -> Path:
    """Locate the parameter files shipped with the installed OpenMM."""
    return Path(app.forcefield.__file__).parent / "data"


def installed_force_field_families() -> tuple[str, ...]:
    """Return the families this OpenMM ships complete with their water models.

    A family needs its protein XML *and* every water/ion file in its directory:
    OpenMM 8.3 and 8.4 ship ``charmm36_2024.xml`` without the matching
    ``charmm36_2024/`` water directory, which would fail only at solvation.
    """
    data_directory = _openmm_data_directory()
    return tuple(
        family
        for family, (protein_xml, water_directory, water_models) in (
            FORCE_FIELD_FAMILIES.items()
        )
        if (data_directory / protein_xml).is_file()
        and all(
            (data_directory / water_directory / water_xml).is_file()
            for water_xml, _ in water_models.values()
        )
    )


def create_forcefield(protein_xml: str, water_xml: str, proteinff: str) -> app.ForceField:
    """Build the force field, explaining a parameter set this OpenMM lacks."""
    try:
        return app.ForceField(protein_xml, water_xml)
    except ValueError as error:
        if "Could not locate file" not in str(error):
            raise
        installed = installed_force_field_families()
        raise ValueError(
            f"OpenMM {openmm.version.version} does not ship the parameters for "
            f"proteinff '{proteinff}' ({error}). Older OpenMM releases omit the "
            "amber19 and charmm36_2024 parameter sets. Upgrade OpenMM, or choose "
            "a family this installation provides: "
            + (", ".join(installed) if installed else "none found")
        ) from error
