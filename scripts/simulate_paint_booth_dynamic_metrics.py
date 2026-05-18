#!/usr/bin/env python3
"""Simulate a metric-level dynamic quasi-steady surrogate for paint-booth airflow.

The model intentionally remains compact:

1. Interpolate steady metric targets from a steady U-sweep library.
2. Apply metric-specific first-order lags.
3. Add direction-specific decaying overshoot states for mean-downdraft steps.
4. Add optional empirical wake/reverse-fraction impulse kernels.
5. Compare predicted metrics to transient CFD time-series CSV files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from bisect import bisect_left
from pathlib import Path
from statistics import mean

METRICS = [
    "supply_inflow_m3s",
    "floor_outflow_m3s",
    "work_zone_Uz_mean",
    "near_car_Uz_mean",
    "near_car_reverse_fraction_Uz_positive",
    "filter_layer_panel_Uz_mean",
    "filter_pressure_drop_proxy_m2_per_s2",
]

DEFAULT_CONFIG = {
    "model": "metric_level_dynamic_quasi_steady_v1",
    "notes": [
        "steady target from linear interpolation over U_supply_mps",
        "first-order metric-specific lag",
        "direction-specific mean-Uz overshoot correction calibrated from paired 10 s high-resolution runs",
        "reverse fraction uses empirical fast/slow signed impulse kernels; it is not assumed first-order",
    ],
    "tau_s": {
        "supply_inflow_m3s": 0.02,
        "floor_outflow_m3s": 0.16,
        "filter_layer_panel_Uz_mean": 0.16,
        "filter_pressure_drop_proxy_m2_per_s2": 0.16,
        "work_zone_Uz_mean": 0.16,
        "near_car_Uz_mean": 0.14,
        "near_car_reverse_fraction_Uz_positive": 8.0,
    },
    "overshoot": {
        "work_zone_Uz_mean": {"up_fraction": 0.10, "down_fraction": 0.14, "up_tau_decay_s": 1.5, "down_tau_decay_s": 0.45},
        "near_car_Uz_mean": {"up_fraction": 0.34, "down_fraction": 0.15, "up_tau_decay_s": 1.8, "down_tau_decay_s": 2.0},
    },
    "impulse_reference_delta_u_mps": 1.09,
    "impulse_kernels": {
        "near_car_reverse_fraction_Uz_positive": {
            "up": [
                {"amplitude": -0.140, "tau_decay_s": 0.90},
                {"amplitude": 0.040, "tau_decay_s": 3.0}
            ],
            "down": [
                {"amplitude": 0.080, "tau_decay_s": 0.90},
                {"amplitude": -0.020, "tau_decay_s": 4.0}
            ]
        }
    },
}


def read_csv_dicts(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_library(path: Path) -> list[dict[str, float]]:
    rows = []
    for row in read_csv_dicts(path):
        rows.append({k: float(v) for k, v in row.items() if k != "case_dir"})
    rows.sort(key=lambda r: r["U_supply_mps"])
    return rows


def interp_metric(library: list[dict[str, float]], u: float, metric: str) -> float:
    xs = [r["U_supply_mps"] for r in library]
    if u <= xs[0]:
        i0, i1 = 0, 1
    elif u >= xs[-1]:
        i0, i1 = len(xs) - 2, len(xs) - 1
    else:
        i1 = bisect_left(xs, u)
        if xs[i1] == u:
            return library[i1][metric]
        i0 = i1 - 1
    x0, x1 = xs[i0], xs[i1]
    y0, y1 = library[i0][metric], library[i1][metric]
    w = (u - x0) / (x1 - x0)
    return y0 + w * (y1 - y0)


def make_step_command(initial_u: float, target_u: float, end_time: float, dt: float, step_time: float) -> list[tuple[float, float]]:
    n = int(round(end_time / dt))
    out = []
    for i in range(n + 1):
        t = round(i * dt, 10)
        # CFD time-series at t=0 represents the copied steady initial field;
        # the boundary-condition step is applied for t > step_time.
        out.append((t, target_u if t > step_time else initial_u))
    return out


def read_command_csv(path: Path) -> list[tuple[float, float]]:
    rows = read_csv_dicts(path)
    out = []
    for row in rows:
        t = float(row.get("time_s", row.get("t", "nan")))
        u = float(row.get("U_supply_mps", row.get("u", "nan")))
        out.append((t, u))
    return sorted(out)


def simulate(library: list[dict[str, float]], command: list[tuple[float, float]], config: dict) -> list[dict[str, float]]:
    if len(command) < 2:
        raise ValueError("Command needs at least two samples")

    y = {m: interp_metric(library, command[0][1], m) for m in METRICS}
    z = {m: 0.0 for m in METRICS}
    impulse_states: list[dict[str, float | str]] = []
    rows: list[dict[str, float]] = []
    prev_t, prev_u = command[0]

    for idx, (t, u) in enumerate(command):
        if idx == 0:
            dt = 0.0
        else:
            dt = t - prev_t
            if dt <= 0:
                raise ValueError("Command times must be strictly increasing")

            prev_targets = {m: interp_metric(library, prev_u, m) for m in METRICS}
            targets = {m: interp_metric(library, u, m) for m in METRICS}

            # Add impulse-like overshoot states at command changes.
            if abs(u - prev_u) > 1e-12:
                direction = "up" if u > prev_u else "down"
                fraction_key = f"{direction}_fraction"
                for m, spec in config.get("overshoot", {}).items():
                    delta_target = targets[m] - prev_targets[m]
                    z[m] += float(spec.get(fraction_key, 0.0)) * delta_target
                impulse_scale = abs(u - prev_u) / max(float(config.get("impulse_reference_delta_u_mps", abs(u - prev_u))), 1e-12)
                for m, per_direction in config.get("impulse_kernels", {}).items():
                    for kernel in per_direction.get(direction, []):
                        impulse_states.append(
                            {
                                "metric": m,
                                "value": impulse_scale * float(kernel.get("amplitude", 0.0)),
                                "tau_decay_s": max(float(kernel.get("tau_decay_s", 1.0)), 1e-12),
                            }
                        )

            for m in METRICS:
                tau = max(float(config["tau_s"].get(m, 0.0)), 1e-12)
                alpha = 1.0 if tau <= 1e-9 else 1.0 - math.exp(-dt / tau)
                y[m] += alpha * (targets[m] - y[m])

            for m, spec in config.get("overshoot", {}).items():
                tau_decay = max(float(spec.get(f"{direction}_tau_decay_s", spec.get("tau_decay_s", 1.0))), 1e-12) if 'direction' in locals() else max(float(spec.get("tau_decay_s", 1.0)), 1e-12)
                z[m] *= math.exp(-dt / tau_decay)

            for state in impulse_states:
                state["value"] = float(state["value"]) * math.exp(-dt / float(state["tau_decay_s"]))
            impulse_states = [state for state in impulse_states if abs(float(state["value"])) > 1e-12]

        row = {"time_s": t, "U_supply_mps": u}
        impulse_by_metric = {m: 0.0 for m in METRICS}
        for state in impulse_states:
            impulse_by_metric[str(state["metric"])] += float(state["value"])
        for m in METRICS:
            row[m] = y[m] + z.get(m, 0.0) + impulse_by_metric.get(m, 0.0)
        rows.append(row)
        prev_t, prev_u = t, u
    return rows


def write_rows(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["time_s", "U_supply_mps", *METRICS]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def interp_series(rows: list[dict[str, float]], t: float, metric: str) -> float:
    times = [float(r["time_s"]) for r in rows]
    if t <= times[0]:
        return float(rows[0][metric])
    if t >= times[-1]:
        return float(rows[-1][metric])
    i1 = bisect_left(times, t)
    if times[i1] == t:
        return float(rows[i1][metric])
    i0 = i1 - 1
    t0, t1 = times[i0], times[i1]
    y0, y1 = float(rows[i0][metric]), float(rows[i1][metric])
    w = (t - t0) / (t1 - t0)
    return y0 + w * (y1 - y0)


def compare(pred_rows: list[dict[str, float]], cfd_path: Path) -> dict:
    cfd_rows = read_csv_dicts(cfd_path)
    result = {"cfd_csv": str(cfd_path), "n_cfd_times": len(cfd_rows), "metrics": {}}
    for m in METRICS:
        if m not in cfd_rows[0]:
            continue
        errs = []
        abs_errs = []
        for row in cfd_rows:
            t = float(row["time_s"])
            obs = float(row[m])
            pred = interp_series(pred_rows, t, m)
            e = pred - obs
            errs.append(e)
            abs_errs.append(abs(e))
        result["metrics"][m] = {
            "rmse": math.sqrt(mean([e * e for e in errs])),
            "mae": mean(abs_errs),
            "max_abs_error": max(abs_errs),
            "bias": mean(errs),
        }
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steady-library", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path)
    ap.add_argument("--config", type=Path)
    ap.add_argument("--command-csv", type=Path)
    ap.add_argument("--initial-u", type=float, default=4.36)
    ap.add_argument("--target-u", type=float, default=5.45)
    ap.add_argument("--end-time", type=float, default=10.0)
    ap.add_argument("--dt", type=float, default=0.25)
    ap.add_argument("--step-time", type=float, default=0.0)
    ap.add_argument("--compare-cfd-csv", type=Path, action="append", default=[])
    args = ap.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8")) if args.config else DEFAULT_CONFIG
    library = read_library(args.steady_library)
    command = read_command_csv(args.command_csv) if args.command_csv else make_step_command(
        args.initial_u, args.target_u, args.end_time, args.dt, args.step_time
    )
    rows = simulate(library, command, config)
    write_rows(args.output_csv, rows)

    summary = {
        "model_config": config,
        "steady_library": str(args.steady_library),
        "output_csv": str(args.output_csv),
        "n_times": len(rows),
        "time_start_s": rows[0]["time_s"],
        "time_end_s": rows[-1]["time_s"],
        "command": {"initial_u": command[0][1], "final_u": command[-1][1]},
        "comparisons": [compare(rows, p) for p in args.compare_cfd_csv],
    }
    out_json = args.output_json or args.output_csv.with_suffix(".json")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
