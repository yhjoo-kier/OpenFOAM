#!/usr/bin/env python3
"""Build a compact steady metric library from paint-booth post_plenum_metrics.json files."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

METRIC_COLUMNS = [
    "U_supply_mps",
    "case_dir",
    "supply_inflow_m3s",
    "floor_outflow_m3s",
    "relative_mass_imbalance",
    "work_zone_Uz_mean",
    "near_car_Uz_mean",
    "near_car_reverse_fraction_Uz_positive",
    "filter_layer_panel_Uz_mean",
    "filter_pressure_drop_proxy_m2_per_s2",
]


def u_from_case_name(path: Path, data: dict) -> float:
    meta_u = data.get("metadata", {}).get("supply_velocity_mps")
    if meta_u is not None:
        return float(meta_u)
    m = re.search(r"u(\d+)p(\d+)", path.as_posix())
    if not m:
        raise ValueError(f"Cannot infer U_supply from {path}")
    return float(f"{m.group(1)}.{m.group(2)}")


def region_uz(data: dict, name: str) -> float:
    return float(data["regions"][name]["Uz"]["mean"])


def reverse_fraction(data: dict, name: str = "near_car") -> float:
    du = data["regions"][name]["downdraft_uniformity"]
    # Historical docs/scripts used both spellings at different points.
    if "reverse_flow_fraction_Uz_positive" in du:
        return float(du["reverse_flow_fraction_Uz_positive"])
    return float(du["reverse_fraction_uz_positive"])


def row_from_metrics(path: Path, repo: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    mb = data["mass_balance"]
    case_dir = path.parent
    return {
        "U_supply_mps": u_from_case_name(case_dir, data),
        "case_dir": str(case_dir.relative_to(repo) if case_dir.is_relative_to(repo) else case_dir),
        "supply_inflow_m3s": float(mb["supply_inflow_m3s"]),
        "floor_outflow_m3s": float(mb["floor_outflow_m3s"]),
        "relative_mass_imbalance": float(mb["relative_imbalance_out_minus_in"]),
        "work_zone_Uz_mean": region_uz(data, "work_zone"),
        "near_car_Uz_mean": region_uz(data, "near_car"),
        "near_car_reverse_fraction_Uz_positive": reverse_fraction(data, "near_car"),
        "filter_layer_panel_Uz_mean": region_uz(data, "filter_layer_panel"),
        "filter_pressure_drop_proxy_m2_per_s2": float(data["filter_pressure_drop_proxy_m2_per_s2"]),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep-dir", type=Path, default=Path("cases/paint_booth_l2_u_sweep_v0"))
    ap.add_argument("--output", type=Path, default=Path("cases/paint_booth_l2_u_sweep_v0/post_steady_metric_library.csv"))
    args = ap.parse_args()

    repo = Path.cwd().resolve()
    sweep_dir = args.sweep_dir.resolve() if args.sweep_dir.is_absolute() else (repo / args.sweep_dir).resolve()
    output = args.output.resolve() if args.output.is_absolute() else (repo / args.output).resolve()

    metric_files = sorted(sweep_dir.glob("*/post_plenum_metrics.json"))
    if not metric_files:
        raise FileNotFoundError(f"No post_plenum_metrics.json files found under {sweep_dir}")

    rows = sorted((row_from_metrics(path, repo) for path in metric_files), key=lambda r: r["U_supply_mps"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "output": str(output.relative_to(repo) if output.is_relative_to(repo) else output),
        "n_rows": len(rows),
        "U_supply_mps": [r["U_supply_mps"] for r in rows],
        "columns": METRIC_COLUMNS,
    }
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
