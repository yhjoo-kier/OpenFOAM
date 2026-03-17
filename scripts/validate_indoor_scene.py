#!/usr/bin/env python3
"""Validate Gemini-generated indoor CFD scene JSON specifications.

This validator targets the indoor_cfd_scene_v1 schema used for early-stage
Gemini -> geometry -> meshing -> OpenFOAM pipeline experiments.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "indoor_cfd_scene_v1"
VALID_WALLS = {"west", "east", "south", "north", "floor", "ceiling"}
VALID_OPENING_TYPES = {"inlet", "outlet"}
EPS = 1e-9


class ValidationReport:
    """Collect validation errors and warnings."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    @property
    def ok(self) -> bool:
        return not self.errors


def is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def require_dict(report: ValidationReport, value: Any, label: str) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        report.error(f"{label} must be an object")
        return None
    return value


def require_list(report: ValidationReport, value: Any, label: str) -> list[Any] | None:
    if not isinstance(value, list):
        report.error(f"{label} must be an array")
        return None
    return value


def require_string(report: ValidationReport, value: Any, label: str) -> str | None:
    if not isinstance(value, str) or not value.strip():
        report.error(f"{label} must be a non-empty string")
        return None
    return value


def require_positive_number(report: ValidationReport, value: Any, label: str) -> float | None:
    if not is_number(value):
        report.error(f"{label} must be a finite number")
        return None
    value = float(value)
    if value <= 0.0:
        report.error(f"{label} must be > 0")
        return None
    return value


def require_nonnegative_number(report: ValidationReport, value: Any, label: str) -> float | None:
    if not is_number(value):
        report.error(f"{label} must be a finite number")
        return None
    value = float(value)
    if value < 0.0:
        report.error(f"{label} must be >= 0")
        return None
    return value


def boxes_overlap(a_min: dict[str, float], a_size: dict[str, float], b_min: dict[str, float], b_size: dict[str, float]) -> bool:
    a_max = {axis: a_min[axis] + a_size[f"d{axis}"] for axis in ("x", "y", "z")}
    b_max = {axis: b_min[axis] + b_size[f"d{axis}"] for axis in ("x", "y", "z")}
    return (
        a_min["x"] < b_max["x"] - EPS and a_max["x"] > b_min["x"] + EPS
        and a_min["y"] < b_max["y"] - EPS and a_max["y"] > b_min["y"] + EPS
        and a_min["z"] < b_max["z"] - EPS and a_max["z"] > b_min["z"] + EPS
    )


def validate_room(report: ValidationReport, scene: dict[str, Any]) -> dict[str, float] | None:
    room = require_dict(report, scene.get("room"), "room")
    if room is None:
        return None

    size = require_dict(report, room.get("size"), "room.size")
    if size is None:
        return None

    room_size: dict[str, float] = {}
    for key in ("Lx", "Ly", "Lz"):
        value = require_positive_number(report, size.get(key), f"room.size.{key}")
        if value is not None:
            room_size[key] = value

    return room_size if len(room_size) == 3 else None


def validate_obstacles(report: ValidationReport, scene: dict[str, Any], room_size: dict[str, float] | None) -> list[dict[str, Any]]:
    obstacles = require_list(report, scene.get("obstacles"), "obstacles")
    if obstacles is None:
        return []

    validated: list[dict[str, Any]] = []
    seen_ids: set[str] = set()

    for idx, obstacle in enumerate(obstacles):
        label = f"obstacles[{idx}]"
        obstacle = require_dict(report, obstacle, label)
        if obstacle is None:
            continue

        obs_id = require_string(report, obstacle.get("id"), f"{label}.id")
        if obs_id:
            if obs_id in seen_ids:
                report.error(f"Duplicate obstacle id: {obs_id}")
            seen_ids.add(obs_id)

        obs_type = require_string(report, obstacle.get("type"), f"{label}.type")
        if obs_type is not None and obs_type != "box":
            report.error(f"{label}.type must be 'box'")

        min_corner = require_dict(report, obstacle.get("min"), f"{label}.min")
        size = require_dict(report, obstacle.get("size"), f"{label}.size")
        if min_corner is None or size is None:
            continue

        min_values: dict[str, float] = {}
        size_values: dict[str, float] = {}
        for axis in ("x", "y", "z"):
            value = require_nonnegative_number(report, min_corner.get(axis), f"{label}.min.{axis}")
            if value is not None:
                min_values[axis] = value

        for axis in ("x", "y", "z"):
            key = f"d{axis}"
            value = require_positive_number(report, size.get(key), f"{label}.size.{key}")
            if value is not None:
                size_values[key] = value

        if len(min_values) != 3 or len(size_values) != 3:
            continue

        if room_size is not None:
            max_corner = {
                "x": min_values["x"] + size_values["dx"],
                "y": min_values["y"] + size_values["dy"],
                "z": min_values["z"] + size_values["dz"],
            }
            if max_corner["x"] > room_size["Lx"] + EPS:
                report.error(f"{label} exceeds room bound in x-direction")
            if max_corner["y"] > room_size["Ly"] + EPS:
                report.error(f"{label} exceeds room bound in y-direction")
            if max_corner["z"] > room_size["Lz"] + EPS:
                report.error(f"{label} exceeds room bound in z-direction")

        validated.append({"id": obs_id, "min": min_values, "size": size_values, "index": idx})

    for i in range(len(validated)):
        for j in range(i + 1, len(validated)):
            a = validated[i]
            b = validated[j]
            if boxes_overlap(a["min"], a["size"], b["min"], b["size"]):
                report.error(
                    f"Obstacle overlap detected between {a['id'] or f'obstacles[{a['index']}]'} and {b['id'] or f'obstacles[{b['index']}]'}"
                )

    return validated


def opening_wall_bounds(room_size: dict[str, float], wall: str) -> tuple[float, float]:
    if wall in {"west", "east"}:
        return room_size["Ly"], room_size["Lz"]
    if wall in {"south", "north"}:
        return room_size["Lx"], room_size["Lz"]
    return room_size["Lx"], room_size["Ly"]


def validate_openings(report: ValidationReport, scene: dict[str, Any], room_size: dict[str, float] | None) -> list[dict[str, Any]]:
    openings = require_list(report, scene.get("openings"), "openings")
    if openings is None:
        return []

    validated: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    inlet_count = 0
    outlet_count = 0

    for idx, opening in enumerate(openings):
        label = f"openings[{idx}]"
        opening = require_dict(report, opening, label)
        if opening is None:
            continue

        opening_id = require_string(report, opening.get("id"), f"{label}.id")
        if opening_id:
            if opening_id in seen_ids:
                report.error(f"Duplicate opening id: {opening_id}")
            seen_ids.add(opening_id)

        opening_type = require_string(report, opening.get("type"), f"{label}.type")
        if opening_type is not None:
            if opening_type not in VALID_OPENING_TYPES:
                report.error(f"{label}.type must be one of {sorted(VALID_OPENING_TYPES)}")
            elif opening_type == "inlet":
                inlet_count += 1
            elif opening_type == "outlet":
                outlet_count += 1

        wall = require_string(report, opening.get("wall"), f"{label}.wall")
        if wall is not None and wall not in VALID_WALLS:
            report.error(f"{label}.wall must be one of {sorted(VALID_WALLS)}")

        center = require_dict(report, opening.get("center"), f"{label}.center")
        size = require_dict(report, opening.get("size"), f"{label}.size")
        if center is None or size is None:
            continue

        u = require_nonnegative_number(report, center.get("u"), f"{label}.center.u")
        v = require_nonnegative_number(report, center.get("v"), f"{label}.center.v")
        du = require_positive_number(report, size.get("du"), f"{label}.size.du")
        dv = require_positive_number(report, size.get("dv"), f"{label}.size.dv")
        if None in (u, v, du, dv) or wall is None or room_size is None or wall not in VALID_WALLS:
            validated.append({"id": opening_id, "type": opening_type, "wall": wall})
            continue

        u_max, v_max = opening_wall_bounds(room_size, wall)
        if u - du / 2.0 < -EPS or u + du / 2.0 > u_max + EPS:
            report.error(f"{label} exceeds wall bounds in local u-direction for wall '{wall}'")
        if v - dv / 2.0 < -EPS or v + dv / 2.0 > v_max + EPS:
            report.error(f"{label} exceeds wall bounds in local v-direction for wall '{wall}'")

        validated.append({"id": opening_id, "type": opening_type, "wall": wall})

    if openings and inlet_count != 1:
        report.error(f"Expected exactly 1 inlet, found {inlet_count}")
    if openings and outlet_count != 1:
        report.error(f"Expected exactly 1 outlet, found {outlet_count}")

    return validated


def validate_meta(report: ValidationReport, scene: dict[str, Any]) -> None:
    meta = require_dict(report, scene.get("meta"), "meta")
    if meta is None:
        return

    require_string(report, meta.get("scene_type"), "meta.scene_type")
    require_string(report, meta.get("description"), "meta.description")


def validate_scene(scene: dict[str, Any]) -> ValidationReport:
    report = ValidationReport()

    schema_version = scene.get("schema_version")
    if schema_version != SCHEMA_VERSION:
        report.error(f"schema_version must be '{SCHEMA_VERSION}'")

    units = scene.get("units")
    if units != "m":
        report.error("units must be 'm'")

    coordinate_system = require_dict(report, scene.get("coordinate_system"), "coordinate_system")
    if coordinate_system is not None:
        expected_coordinate_system = {
            "origin": "room_min_corner",
            "x": "west_to_east",
            "y": "south_to_north",
            "z": "floor_to_ceiling",
        }
        for key, expected in expected_coordinate_system.items():
            if coordinate_system.get(key) != expected:
                report.error(f"coordinate_system.{key} must be '{expected}'")

    room_size = validate_room(report, scene)
    obstacles = validate_obstacles(report, scene, room_size)
    openings = validate_openings(report, scene, room_size)
    validate_meta(report, scene)

    if not obstacles:
        report.warn("No valid obstacles parsed")
    elif len(obstacles) < 3 or len(obstacles) > 5:
        report.warn(f"Prototype prompt expected 3-5 obstacles, found {len(obstacles)}")

    if len(openings) != 2:
        report.warn(f"Prototype prompt expected 2 openings, found {len(openings)}")

    return report


def format_report(report: ValidationReport, path: Path) -> str:
    lines = [f"Validation result for {path}", f"Status: {'OK' if report.ok else 'FAILED'}"]
    if report.errors:
        lines.append("Errors:")
        lines.extend(f"  - {msg}" for msg in report.errors)
    if report.warnings:
        lines.append("Warnings:")
        lines.extend(f"  - {msg}" for msg in report.warnings)
    if report.ok and not report.warnings:
        lines.append("No issues found.")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate indoor_cfd_scene_v1 JSON files for Gemini/OpenFOAM prototyping"
    )
    parser.add_argument("scene_json", type=Path, help="Path to indoor scene JSON file")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable validation output as JSON",
    )
    args = parser.parse_args()

    try:
        with args.scene_json.open("r", encoding="utf-8") as f:
            scene = json.load(f)
    except FileNotFoundError:
        print(f"Error: file not found: {args.scene_json}", file=sys.stderr)
        return 2
    except json.JSONDecodeError as exc:
        print(f"Error: invalid JSON in {args.scene_json}: {exc}", file=sys.stderr)
        return 2

    if not isinstance(scene, dict):
        print("Error: root JSON value must be an object", file=sys.stderr)
        return 2

    report = validate_scene(scene)

    if args.json:
        payload = {
            "ok": report.ok,
            "errors": report.errors,
            "warnings": report.warnings,
            "schema_version": scene.get("schema_version"),
            "path": str(args.scene_json),
        }
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print(format_report(report, args.scene_json))

    return 0 if report.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
