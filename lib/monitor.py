"""Persistent, PBC-aware ligand-detachment monitoring for production MD."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from openmm import app, unit

from ommflow.lib.ligands import PROTEIN_RESIDUE_NAMES


MONITOR_FIELDS = (
    "production_time_ns",
    "step",
    "min_ligand_pocket_distance_nm",
    "contact_count",
    "initial_contact_count",
    "contact_fraction",
    "consecutive_detached_count",
    "detached",
)


@dataclass(frozen=True)
class MonitorDefinition:
    """The immutable production-start pocket used by a detachment monitor."""

    ligand_id: str
    ligand_atom_indices: tuple[int, ...]
    ligand_heavy_atom_indices: tuple[int, ...]
    pocket_atom_indices: tuple[int, ...]
    pocket_residues: tuple[dict[str, object], ...]
    initial_contact_count: int
    settings: dict[str, object]
    component_id: str | None = None
    selector_kind: str | None = None
    selector_value: str | None = None

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-compatible representation of the saved selection."""
        return {
            "schema_version": 2,
            "ligand_id": self.ligand_id,
            "component_id": self.component_id,
            "ligand_atom_indices": list(self.ligand_atom_indices),
            "ligand_heavy_atom_indices": list(self.ligand_heavy_atom_indices),
            "target_atom_indices": list(self.ligand_atom_indices),
            "target_heavy_atom_indices": list(self.ligand_heavy_atom_indices),
            "target_selection": {
                "component_id": self.component_id,
                "ligand_id": self.ligand_id if self.component_id != self.ligand_id else None,
                "selector": self.selector_kind,
                "selector_value": self.selector_value,
            },
            "pocket_atom_indices": list(self.pocket_atom_indices),
            "pocket_residues": list(self.pocket_residues),
            "initial_contact_count": self.initial_contact_count,
            "settings": self.settings,
        }


def monitor_settings(args) -> dict[str, object]:
    """Return all settings that define a persistent monitor's behavior."""
    return {
        "monitor_interval_ns": args.monitor_interval_ns,
        "pocket_cutoff_nm": args.pocket_cutoff_nm,
        "contact_cutoff_nm": args.contact_cutoff_nm,
        "detach_cutoff_nm": args.detach_cutoff_nm,
        "confirmation_checks": args.confirmation_checks,
    }


def ligand_candidates(components: dict[str, object]) -> list[dict[str, object]]:
    """Return component records that have durable ligand IDs."""
    records = components.get("components", [])
    if not isinstance(records, list):
        raise ValueError("components.json has an invalid 'components' field.")
    return [
        component
        for component in records
        if isinstance(component, dict)
        and component.get("classification") == "ligand_candidate"
        and isinstance(component.get("ligand_id"), str)
    ]


def select_monitor_ligand(
    components: dict[str, object], requested_ligand_id: str | None
) -> dict[str, object]:
    """Select one parameterized ligand, requiring an ID when ambiguous."""
    candidates = ligand_candidates(components)
    available = [
        component
        for component in candidates
        if isinstance(component.get("production_topology_atom_indices"), list)
        and component["production_topology_atom_indices"]
    ]
    candidate_ids = [str(component["ligand_id"]) for component in candidates]
    available_ids = [str(component["ligand_id"]) for component in available]
    if requested_ligand_id is not None:
        matching = [
            component
            for component in available
            if component["ligand_id"] == requested_ligand_id
        ]
        if matching:
            return matching[0]
        if requested_ligand_id in candidate_ids:
            raise ValueError(
                f"Selected ligand '{requested_ligand_id}' was not parameterized "
                "and cannot be monitored."
            )
        raise ValueError(
            f"Unknown monitor ligand '{requested_ligand_id}'. Candidates: "
            + (", ".join(candidate_ids) if candidate_ids else "none")
        )
    if not available:
        raise ValueError(
            "Early-stop monitoring requires one automatically parameterized "
            "ligand candidate, but none are available."
        )
    if len(available) != 1:
        raise ValueError(
            "Early-stop monitoring is ambiguous; specify --monitor-ligand with "
            "one of: " + ", ".join(available_ids)
        )
    return available[0]


def select_monitor_target(components: dict[str, object], args) -> dict[str, object]:
    """Resolve one explicit target, or the legacy single-GAFF automatic target."""
    selectors = _requested_selectors(args)
    if len(selectors) > 1:
        raise ValueError(
            "Specify at most one monitor target selector: "
            + ", ".join(f"--{key.replace('_', '-')}" for key, _ in selectors)
        )
    records = _component_records(components)
    monitorable = [record for record in records if _is_monitorable(record)]
    if not selectors:
        candidates = ligand_candidates(components)
        if len(candidates) == 1:
            return _selected_target(candidates[0], "automatic", None)
        if len(candidates) > 1:
            raise ValueError(
                "Early-stop monitoring is ambiguous; specify --monitor-ligand with "
                "one of: "
                + ", ".join(str(component["ligand_id"]) for component in candidates)
            )
        raise ValueError(
            "Early-stop monitoring has no automatically parameterized GAFF ligand. "
            "Select a peptide or standard component explicitly with --monitor-chain "
            "or --monitor-component. Available monitorable components: "
            + _component_descriptions(monitorable)
        )

    selector, requested = selectors[0]
    if selector == "monitor_ligand":
        selected = select_monitor_ligand(components, requested)
        return _selected_target(selected, selector, requested)
    if selector == "monitor_component":
        if not (
            requested.startswith("component-")
            and requested.removeprefix("component-").isdigit()
        ):
            raise ValueError(
                f"Monitor component '{requested}' must be an exact component-N ID."
            )
        selected = [
            component
            for component in monitorable
            if _component_id(component) == requested
        ]
        if selected:
            return _selected_target(selected[0], selector, requested)
        raise ValueError(
            f"Unknown or non-monitorable component '{requested}'. Candidates: "
            + _component_descriptions(monitorable)
        )

    matching = [
        component
        for component in monitorable
        if requested in _component_chains(component)
    ]
    if len(matching) == 1:
        return _selected_target(matching[0], selector, requested)
    if not matching:
        raise ValueError(
            f"No monitorable disconnected component contains input chain '{requested}'. "
            "Candidates: " + _component_descriptions(monitorable)
        )
    raise ValueError(
        f"Input chain '{requested}' maps to multiple monitorable disconnected "
        "components: " + _component_descriptions(matching)
    )


def _requested_selectors(args) -> list[tuple[str, str]]:
    return [
        (key, value)
        for key in ("monitor_ligand", "monitor_chain", "monitor_component")
        if isinstance((value := getattr(args, key, None)), str) and value
    ]


def _component_records(components: dict[str, object]) -> list[dict[str, object]]:
    records = components.get("components", [])
    if not isinstance(records, list):
        raise ValueError("components.json has an invalid 'components' field.")
    return [record for record in records if isinstance(record, dict)]


def _component_id(component: dict[str, object]) -> str:
    value = component.get("component_id")
    if isinstance(value, str):
        return value
    index = component.get("component_index")
    return f"component-{index}" if isinstance(index, int) else ""


def _component_chains(component: dict[str, object]) -> list[str]:
    chains = component.get("chains")
    if isinstance(chains, list):
        return [str(chain) for chain in chains]
    residues = component.get("residues")
    if isinstance(residues, list):
        return sorted(
            {
                str(residue["chain_id"])
                for residue in residues
                if isinstance(residue, dict) and "chain_id" in residue
            }
        )
    return []


def _is_monitorable(component: dict[str, object]) -> bool:
    if "monitorable" in component:
        return component.get("monitorable") is True
    residues = component.get("residues")
    if not isinstance(residues, list) or not residues:
        return component.get("classification") != "covalent_candidate"
    solvent_and_ions = {
        "BR", "CA", "CL", "CS", "F", "HOH", "I", "K", "LI", "MG", "NA", "OPC",
        "OPC3", "RB", "SOL", "SPC", "SPCE", "TIP3", "TIP4", "TIP5", "TP3",
        "TP4", "TP5", "WAT", "ZN",
    }
    return (
        component.get("classification") != "covalent_candidate"
        and any(
            str(residue.get("residue_name", "")).upper() not in solvent_and_ions
            for residue in residues
            if isinstance(residue, dict)
        )
    )


def _production_indices(component: dict[str, object]) -> tuple[int, ...]:
    indices = component.get("production_topology_atom_indices")
    if not isinstance(indices, list) or not indices:
        return ()
    try:
        return tuple(int(index) for index in indices)
    except (TypeError, ValueError):
        return ()


def _selected_target(
    component: dict[str, object], selector_kind: str, selector_value: str | None
) -> dict[str, object]:
    indices = _production_indices(component)
    if not indices:
        raise ValueError(
            f"Selected component '{_component_id(component)}' has no final production "
            "topology atom mapping and cannot be monitored."
        )
    selected = dict(component)
    selected["_monitor_component_id"] = _component_id(component)
    selected["_monitor_selector_kind"] = selector_kind
    selected["_monitor_selector_value"] = selector_value
    return selected


def _component_descriptions(components: Iterable[dict[str, object]]) -> str:
    descriptions = [
        f"{_component_id(component)} (chains: {', '.join(_component_chains(component)) or 'none'}; "
        f"classification: {component.get('classification', 'unknown')})"
        for component in components
    ]
    return ", ".join(descriptions) if descriptions else "none"


def load_components(path: Path) -> dict[str, object]:
    """Load component provenance needed for durable ligand selection."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Early-stop monitoring requires component provenance: {path}"
        ) from None
    except json.JSONDecodeError as error:
        raise ValueError(f"Could not read component provenance {path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"Component provenance {path} must contain a JSON object.")
    return data


def initialize_monitor(
    topology: app.Topology, positions, box_vectors, components: dict[str, object], args
) -> MonitorDefinition:
    """Select a target and establish its immutable protein heavy-atom pocket."""
    selected = select_monitor_target(components, args)
    target_id = str(selected["_monitor_component_id"])
    ligand_atoms = _production_indices(selected)
    topology_atoms = list(topology.atoms())
    _validate_topology_indices(ligand_atoms, topology_atoms, "selected target")
    ligand_heavy_atoms = tuple(
        index for index in ligand_atoms if _is_heavy_atom(topology_atoms[index])
    )
    if not ligand_heavy_atoms:
        raise ValueError(
            f"Selected target '{target_id}' has no heavy atoms to monitor."
        )

    excluded_atoms = set(ligand_atoms)
    excluded_atoms.update(
        index
        for component in ligand_candidates(components)
        for index in _production_indices(component)
    )
    protein_atoms = {
        atom.index
        for atom in topology_atoms
        if atom.residue.name.upper() in PROTEIN_RESIDUE_NAMES
    }
    if protein_atoms and protein_atoms <= set(ligand_atoms):
        raise ValueError(
            f"Selected target '{target_id}' contains all non-water protein atoms; "
            "choose a proper receptor subset rather than the whole protein."
        )
    protein_heavy_atoms = [
        atom.index
        for atom in topology_atoms
        if atom.index not in excluded_atoms
        and atom.residue.name.upper() in PROTEIN_RESIDUE_NAMES
        and _is_heavy_atom(atom)
    ]
    coordinates = _coordinates_nm(positions)
    box = _box_matrices(_box_vectors_nm(box_vectors))
    pocket_atom_indices = tuple(
        atom_index
        for atom_index, distance in zip(
            protein_heavy_atoms,
            _selection_minimum_distances(
                coordinates, protein_heavy_atoms, ligand_heavy_atoms, box
            ),
            strict=True,
        )
        if distance <= args.pocket_cutoff_nm
    )
    if not pocket_atom_indices:
        raise ValueError(
            f"Selected target '{target_id}' has no protein heavy atoms "
            f"within --pocket-cutoff-nm {args.pocket_cutoff_nm:g}; it is not "
            "bound at the start of production."
        )
    _, initial_contacts = monitor_measurement(
        ligand_heavy_atoms,
        pocket_atom_indices,
        positions,
        box_vectors,
        args.contact_cutoff_nm,
    )
    return MonitorDefinition(
        ligand_id=str(selected.get("ligand_id") or target_id),
        ligand_atom_indices=ligand_atoms,
        ligand_heavy_atom_indices=ligand_heavy_atoms,
        pocket_atom_indices=pocket_atom_indices,
        pocket_residues=_pocket_residues(topology_atoms, pocket_atom_indices),
        initial_contact_count=initial_contacts,
        settings=monitor_settings(args),
        component_id=target_id,
        selector_kind=str(selected["_monitor_selector_kind"]),
        selector_value=(
            str(selected["_monitor_selector_value"])
            if selected["_monitor_selector_value"] is not None
            else None
        ),
    )


def _selection_minimum_distances(coordinates, first_indices, second_indices, box):
    """Return each first-selection atom's shortest distance to the second."""
    distances = []
    for block in _pair_distance_blocks(
        coordinates, tuple(first_indices), tuple(second_indices), box
    ):
        distances.extend(block.min(axis=1).tolist())
    return distances


def write_pocket(path: Path, definition: MonitorDefinition) -> None:
    """Persist the production-start monitor definition."""
    path.write_text(json.dumps(definition.as_dict(), indent=2) + "\n", encoding="utf-8")


def load_pocket(path: Path) -> MonitorDefinition:
    """Restore the immutable pocket without recomputing it on restart."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Early-stop monitoring cannot resume without its pocket definition: {path}"
        ) from None
    except json.JSONDecodeError as error:
        raise ValueError(f"Could not read pocket definition {path}: {error}") from error
    required = {
        "ligand_id",
        "ligand_atom_indices",
        "ligand_heavy_atom_indices",
        "pocket_atom_indices",
        "pocket_residues",
        "initial_contact_count",
        "settings",
    }
    if not isinstance(data, dict) or not required <= data.keys():
        raise ValueError(f"Pocket definition {path} is incomplete.")
    try:
        return MonitorDefinition(
            ligand_id=str(data["ligand_id"]),
            ligand_atom_indices=tuple(int(index) for index in data["ligand_atom_indices"]),
            ligand_heavy_atom_indices=tuple(
                int(index) for index in data["ligand_heavy_atom_indices"]
            ),
            pocket_atom_indices=tuple(int(index) for index in data["pocket_atom_indices"]),
            pocket_residues=tuple(data["pocket_residues"]),
            initial_contact_count=int(data["initial_contact_count"]),
            settings=dict(data["settings"]),
            component_id=(
                str(data["component_id"])
                if isinstance(data.get("component_id"), str)
                else None
            ),
            selector_kind=_saved_selector_value(data, "selector"),
            selector_value=_saved_selector_value(data, "selector_value"),
        )
    except (TypeError, ValueError) as error:
        raise ValueError(f"Pocket definition {path} contains invalid data.") from error


def validate_saved_monitor(definition: MonitorDefinition, args) -> None:
    """Ensure a restart cannot silently change a persisted monitor."""
    requested = _requested_selectors(args)
    legacy_ligand_match = (
        definition.selector_kind is None
        and requested == [("monitor_ligand", definition.ligand_id)]
    )
    if requested and not legacy_ligand_match and (
        len(requested) != 1
        or definition.selector_kind is None
        or requested[0] != (definition.selector_kind, definition.selector_value)
    ):
        raise ValueError(
            "Cannot change the monitored target on restart; the persisted pocket "
            f"uses '{definition.ligand_id}'."
        )
    saved_settings = definition.settings
    for key, value in monitor_settings(args).items():
        saved_value = saved_settings.get(key)
        if isinstance(value, float):
            matches = isinstance(saved_value, (int, float)) and math.isclose(
                value, float(saved_value), rel_tol=0.0, abs_tol=1e-12
            )
        else:
            matches = value == saved_value
        if not matches:
            raise ValueError(
                f"Cannot change {key} on restart; it is fixed by "
                f"{definition.ligand_id}'s persisted pocket definition."
            )


def validate_monitor_topology(
    definition: MonitorDefinition, topology: app.Topology
) -> None:
    """Verify a restored selection still addresses the saved production topology."""
    topology_atoms = list(topology.atoms())
    _validate_topology_indices(
        definition.ligand_atom_indices, topology_atoms, "selected ligand"
    )
    _validate_topology_indices(
        definition.pocket_atom_indices, topology_atoms, "pocket"
    )


def monitor_measurement(
    ligand_heavy_atom_indices: Iterable[int],
    pocket_atom_indices: Iterable[int],
    positions,
    box_vectors,
    contact_cutoff_nm: float,
) -> tuple[float, int]:
    """Return minimum PBC-aware distance and the number of contacting atom pairs."""
    ligand_indices = tuple(ligand_heavy_atom_indices)
    pocket_indices = tuple(pocket_atom_indices)
    coordinates = _coordinates_nm(positions)
    box = _box_matrices(_box_vectors_nm(box_vectors))
    minimum = math.inf
    contacts = 0
    for distances in _pair_distance_blocks(
        coordinates, ligand_indices, pocket_indices, box
    ):
        minimum = min(minimum, float(distances.min()))
        contacts += int(np.count_nonzero(distances <= contact_cutoff_nm))
    return minimum, contacts


def detachment_update(
    minimum_distance_nm: float,
    contact_count: int,
    detach_cutoff_nm: float,
    previous_consecutive_count: int,
    confirmation_checks: int,
) -> tuple[int, bool]:
    """Advance the confirmation counter and report confirmed detachment."""
    detached_now = (
        contact_count == 0 and minimum_distance_nm > detach_cutoff_nm
    )
    consecutive_count = previous_consecutive_count + 1 if detached_now else 0
    return consecutive_count, consecutive_count >= confirmation_checks


def append_monitor_row(path: Path, row: dict[str, object]) -> bool:
    """Append one monitor result, never duplicating a checkpointed step."""
    existing_step = last_monitor_step(path)
    step = int(row["step"])
    if existing_step is not None:
        if existing_step == step:
            return False
        if existing_step > step:
            raise ValueError(
                f"monitor.csv already contains step {existing_step}, later than {step}."
            )
    has_contents = path.is_file() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MONITOR_FIELDS)
        if not has_contents:
            writer.writeheader()
        writer.writerow({field: row[field] for field in MONITOR_FIELDS})
    return True


def last_monitor_step(path: Path) -> int | None:
    """Return the last persisted monitor step, if one exists."""
    row = last_monitor_row(path)
    return None if row is None else int(row["step"])


def last_monitor_row(path: Path) -> dict[str, str] | None:
    """Read only the final monitor CSV row for restart state restoration."""
    if not path.is_file() or path.stat().st_size == 0:
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        last = None
        for row in rows:
            last = row
    return last


def load_status(path: Path) -> dict[str, object] | None:
    """Load the last durable monitor status, if monitoring has begun."""
    if not path.is_file():
        return None
    try:
        status = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Could not read monitor status {path}: {error}") from error
    if not isinstance(status, dict):
        raise ValueError(f"Monitor status {path} must contain a JSON object.")
    return status


def write_status(path: Path, status: dict[str, object]) -> None:
    """Persist monitor outcome and restart state."""
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def monitor_status(
    outcome: str,
    definition: MonitorDefinition,
    target_production_ns: float,
    final_production_time_ns: float,
    final_production_step: int,
    consecutive_detached_count: int,
) -> dict[str, object]:
    """Build the durable outcome record for status.json."""
    return {
        "schema_version": 1,
        "outcome": outcome,
        "target_production_ns": target_production_ns,
        "final_production_time_ns": final_production_time_ns,
        "final_production_step": final_production_step,
        "monitor_selection": {
            "target_id": definition.component_id or definition.ligand_id,
            "component_id": definition.component_id,
            "ligand_id": definition.ligand_id,
            "selector": definition.selector_kind,
            "selector_value": definition.selector_value,
            "target_atom_indices": list(definition.ligand_atom_indices),
            "ligand_atom_indices": list(definition.ligand_atom_indices),
            "pocket_atom_indices": list(definition.pocket_atom_indices),
        },
        "consecutive_detached_count": consecutive_detached_count,
    }


def _saved_selector_value(data: dict[str, object], key: str) -> str | None:
    selection = data.get("target_selection")
    if isinstance(selection, dict) and isinstance(selection.get(key), str):
        return str(selection[key])
    return None


def restore_consecutive_count(
    status: dict[str, object] | None, monitor_csv: Path, current_step: int
) -> int:
    """Restore a counter from status, or reconstruct it from the last CSV row."""
    if status is not None:
        status_step = status.get("final_production_step")
        count = status.get("consecutive_detached_count")
        if (
            isinstance(status_step, int)
            and status_step <= current_step
            and isinstance(count, int)
            and count >= 0
        ):
            return count
    row = last_monitor_row(monitor_csv)
    if row is None:
        return 0
    try:
        if int(row["step"]) <= current_step:
            return int(row["consecutive_detached_count"])
    except (KeyError, ValueError):
        pass
    return 0


def _pocket_residues(
    topology_atoms: list, pocket_atom_indices: tuple[int, ...]
) -> tuple[dict[str, object], ...]:
    residue_atoms: dict[tuple[str, str, str, str], list[int]] = {}
    for atom_index in pocket_atom_indices:
        residue = topology_atoms[atom_index].residue
        key = (
            residue.chain.id,
            residue.name,
            residue.id,
            residue.insertionCode,
        )
        residue_atoms.setdefault(key, []).append(atom_index)
    return tuple(
        {
            "chain_id": chain_id,
            "residue_name": residue_name,
            "residue_id": residue_id,
            "insertion_code": insertion_code,
            "atom_indices": atom_indices,
        }
        for (chain_id, residue_name, residue_id, insertion_code), atom_indices in residue_atoms.items()
    )


def _validate_topology_indices(
    atom_indices: Iterable[int], topology_atoms: list, label: str
) -> None:
    invalid = [index for index in atom_indices if index < 0 or index >= len(topology_atoms)]
    if invalid:
        raise ValueError(
            f"Persisted {label} atom indices do not match the production topology: "
            + ", ".join(str(index) for index in invalid)
        )


def _is_heavy_atom(atom) -> bool:
    return atom.element is not None and atom.element.atomic_number > 1


def _minimum_distance_nm(
    coordinates,
    first_indices: Iterable[int],
    second_indices: Iterable[int],
    box_vectors: tuple[tuple[float, float, float], ...] | None,
) -> float:
    """Return the shortest minimum-image distance between two atom selections."""
    return _minimum_distance_with_box(
        coordinates, first_indices, second_indices, _box_matrices(box_vectors)
    )


def _minimum_distance_with_box(
    coordinates, first_indices: Iterable[int], second_indices: Iterable[int], box
) -> float:
    minimum = math.inf
    for distances in _pair_distance_blocks(
        coordinates, tuple(first_indices), tuple(second_indices), box
    ):
        minimum = min(minimum, float(distances.min()))
    return minimum


def _distance_nm(
    coordinates,
    atom1_index: int,
    atom2_index: int,
    box_vectors: tuple[tuple[float, float, float], ...] | None,
) -> float:
    """Return the minimum-image distance between one pair of atoms."""
    return _minimum_distance_nm(
        coordinates, (atom1_index,), (atom2_index,), box_vectors
    )


# Pair blocks are capped so a large target against a large receptor never
# materializes one huge displacement array.
_PAIR_BLOCK_ATOMS = 512


def _pair_distance_blocks(coordinates, first_indices, second_indices, box):
    """Yield minimum-image distance matrices for chunks of the first selection."""
    coordinates = np.asarray(coordinates, dtype=float)
    if not first_indices or not second_indices:
        return
    second = coordinates[np.asarray(second_indices, dtype=int)]
    first_array = np.asarray(first_indices, dtype=int)
    for start in range(0, len(first_array), _PAIR_BLOCK_ATOMS):
        first = coordinates[first_array[start : start + _PAIR_BLOCK_ATOMS]]
        displacement = second[None, :, :] - first[:, None, :]
        yield np.sqrt(
            np.square(_minimum_image_array(displacement, box)).sum(axis=-1)
        )


def _coordinates_nm(positions):
    vectors = (
        positions.value_in_unit(unit.nanometer)
        if hasattr(positions, "value_in_unit")
        else positions
    )
    return np.asarray(vectors, dtype=float)


def _box_vectors_nm(box_vectors) -> tuple[tuple[float, float, float], ...] | None:
    if box_vectors is None:
        return None
    vectors = (
        box_vectors.value_in_unit(unit.nanometer)
        if hasattr(box_vectors, "value_in_unit")
        else box_vectors
    )
    return tuple(
        (float(vector[0]), float(vector[1]), float(vector[2])) for vector in vectors
    )


def _box_matrices(box_vectors: tuple[tuple[float, float, float], ...] | None):
    """Invert the triclinic box once so pair loops never repeat that work."""
    if box_vectors is None:
        return None
    matrix = np.asarray(box_vectors, dtype=float)
    if abs(float(np.linalg.det(matrix))) < 1e-15:
        return None
    return matrix, np.linalg.inv(matrix)


def _minimum_image_array(displacement, box):
    """Wrap displacements with the full triclinic box matrix."""
    if box is None:
        return displacement
    matrix, inverse = box
    fractional = displacement @ inverse
    fractional -= np.floor(fractional + 0.5)
    return fractional @ matrix
