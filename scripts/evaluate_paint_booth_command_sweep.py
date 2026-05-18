#!/usr/bin/env python3
"""Evaluate paint-booth command-history scenarios with the dynamic metric surrogate.

For each command CSV in a scenario manifest, this script runs the metric-level
surrogate, computes control-oriented objective metrics, compares dynamic output
to a memoryless steady-library prediction, and classifies the scenario into a
screening category.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

# Allow importing the sibling surrogate module when executed from repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate_paint_booth_dynamic_metrics import (  # noqa: E402
    DEFAULT_CONFIG,
    METRICS,
    interp_metric,
    read_command_csv,
    read_library,
    simulate,
    write_rows,
)

MEAN_UZ_METRICS = ["work_zone_Uz_mean", "near_car_Uz_mean", "filter_layer_panel_Uz_mean"]
REVERSE_METRIC = "near_car_reverse_fraction_Uz_positive"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def integrate_positive(rows: list[dict[str, float]], value_fn) -> tuple[float, float]:
    """Trapezoidal integral of positive value_fn(row) and positive duration."""
    if len(rows) < 2:
        return 0.0, 0.0
    integ = 0.0
    duration = 0.0
    for a, b in zip(rows[:-1], rows[1:]):
        dt = b["time_s"] - a["time_s"]
        va = max(0.0, value_fn(a))
        vb = max(0.0, value_fn(b))
        integ += 0.5 * (va + vb) * dt
        if va > 0.0 or vb > 0.0:
            duration += dt
    return integ, duration


def integrate(rows: list[dict[str, float]], value_fn) -> float:
    if len(rows) < 2:
        return 0.0
    return sum(0.5 * (value_fn(a) + value_fn(b)) * (b["time_s"] - a["time_s"]) for a, b in zip(rows[:-1], rows[1:]))


def steady_memoryless_rows(library: list[dict[str, float]], command: list[tuple[float, float]]) -> list[dict[str, float]]:
    rows = []
    for t, u in command:
        row = {"time_s": t, "U_supply_mps": u}
        for m in METRICS:
            row[m] = interp_metric(library, u, m)
        rows.append(row)
    return rows


def max_abs_diff(dynamic_rows: list[dict[str, float]], steady_rows: list[dict[str, float]], metric: str) -> float:
    return max(abs(a[metric] - b[metric]) for a, b in zip(dynamic_rows, steady_rows))


def scenario_summary(manifest_row: dict, dynamic_rows: list[dict[str, float]], steady_rows: list[dict[str, float]], baseline_reverse: float) -> dict:
    work_mag = [-r["work_zone_Uz_mean"] for r in dynamic_rows]
    near_mag = [-r["near_car_Uz_mean"] for r in dynamic_rows]
    reverse = [r[REVERSE_METRIC] for r in dynamic_rows]
    u_vals = [r["U_supply_mps"] for r in dynamic_rows]

    weak020_int, weak020_dur = integrate_positive(dynamic_rows, lambda r: 0.20 - (-r["work_zone_Uz_mean"]))
    weak024_int, weak024_dur = integrate_positive(dynamic_rows, lambda r: 0.24 - (-r["work_zone_Uz_mean"]))
    weak025_int, weak025_dur = integrate_positive(dynamic_rows, lambda r: 0.25 - (-r["work_zone_Uz_mean"]))
    reverse_excess_int, reverse_excess_dur = integrate_positive(dynamic_rows, lambda r: r[REVERSE_METRIC] - baseline_reverse)
    u_int = integrate(dynamic_rows, lambda r: r["U_supply_mps"])
    u3_int = integrate(dynamic_rows, lambda r: r["U_supply_mps"] ** 3)

    diff_work = max_abs_diff(dynamic_rows, steady_rows, "work_zone_Uz_mean")
    diff_near = max_abs_diff(dynamic_rows, steady_rows, "near_car_Uz_mean")
    diff_reverse = max_abs_diff(dynamic_rows, steady_rows, REVERSE_METRIC)
    diff_filter = max_abs_diff(dynamic_rows, steady_rows, "filter_layer_panel_Uz_mean")
    max_dynamic_steady_delta = max(diff_work, diff_near, diff_filter, diff_reverse)

    high_u = float(manifest_row["high_u_mps"])
    hold_s = float(manifest_row["hold_s"])
    ramp_s = float(manifest_row["ramp_s"])

    # Classification is intentionally conservative. It answers: which prediction
    # tier is justified for the scenario, not whether the scenario is good/bad.
    steady_sufficient = (
        diff_work < 0.004
        and diff_near < 0.004
        and diff_filter < 0.004
        and diff_reverse < 0.010
        and reverse_excess_int < 0.02
    )
    outside_transient_calibration = high_u > 5.45 + 1e-9 or ramp_s > 0.0
    short_cycle_or_risky = hold_s < 15.0 or max(reverse) > baseline_reverse + 0.030 or weak020_dur > 0.0

    if steady_sufficient:
        tier = "steady_library_sufficient"
    elif outside_transient_calibration and short_cycle_or_risky:
        tier = "transient_cfd_priority"
    elif outside_transient_calibration:
        tier = "dynamic_surrogate_with_cfd_spot_check"
    else:
        tier = "dynamic_surrogate_needed"

    return {
        **manifest_row,
        "prediction_tier": tier,
        "work_downdraft_min_mps": min(work_mag),
        "work_downdraft_mean_mps": sum(work_mag) / len(work_mag),
        "near_car_downdraft_min_mps": min(near_mag),
        "near_car_downdraft_mean_mps": sum(near_mag) / len(near_mag),
        "reverse_peak": max(reverse),
        "reverse_mean": sum(reverse) / len(reverse),
        "reverse_excess_integral_above_low_baseline": reverse_excess_int,
        "reverse_excess_duration_s": reverse_excess_dur,
        "weak_work_downdraft_integral_below_0p20": weak020_int,
        "weak_work_downdraft_duration_below_0p20_s": weak020_dur,
        "weak_work_downdraft_integral_below_0p24": weak024_int,
        "weak_work_downdraft_duration_below_0p24_s": weak024_dur,
        "weak_work_downdraft_integral_below_0p25": weak025_int,
        "weak_work_downdraft_duration_below_0p25_s": weak025_dur,
        "u_integral_mps_s": u_int,
        "u3_energy_proxy": u3_int,
        "max_dynamic_minus_steady_work_Uz": diff_work,
        "max_dynamic_minus_steady_near_car_Uz": diff_near,
        "max_dynamic_minus_steady_filter_Uz": diff_filter,
        "max_dynamic_minus_steady_reverse": diff_reverse,
        "max_dynamic_minus_steady_any_selected": max_dynamic_steady_delta,
        "U_min": min(u_vals),
        "U_max": max(u_vals),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--steady-library", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--config", type=Path)
    args = ap.parse_args()

    manifest = read_json(args.manifest)
    library = read_library(args.steady_library)
    config = read_json(args.config) if args.config else DEFAULT_CONFIG
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pred_dir = args.output_dir / "predictions"

    baseline_low_u = float(manifest[0]["low_u_mps"])
    baseline_reverse = interp_metric(library, baseline_low_u, REVERSE_METRIC)

    summaries = []
    for item in manifest:
        sid = item["scenario_id"]
        command = read_command_csv(Path(item["command_csv"]))
        dynamic_rows = simulate(library, command, config)
        steady_rows = steady_memoryless_rows(library, command)
        pred_csv = pred_dir / f"{sid}_pred.csv"
        steady_csv = pred_dir / f"{sid}_steady_memoryless.csv"
        write_rows(pred_csv, dynamic_rows)
        write_rows(steady_csv, steady_rows)
        rec = scenario_summary(item, dynamic_rows, steady_rows, baseline_reverse)
        rec["prediction_csv"] = str(pred_csv)
        rec["steady_memoryless_csv"] = str(steady_csv)
        summaries.append(rec)

    summary_csv = args.output_dir / "scenario_summary.csv"
    summary_json = args.output_dir / "scenario_summary.json"
    write_csv(summary_csv, summaries)
    summary_json.write_text(json.dumps(summaries, indent=2), encoding="utf-8")

    tiers: dict[str, int] = {}
    for row in summaries:
        tiers[row["prediction_tier"]] = tiers.get(row["prediction_tier"], 0) + 1

    def top(key: str, n: int = 5, reverse: bool = True) -> list[dict]:
        return sorted(summaries, key=lambda r: float(r[key]), reverse=reverse)[:n]

    report = {
        "n_scenarios": len(summaries),
        "baseline_low_u_mps": baseline_low_u,
        "baseline_reverse_fraction": baseline_reverse,
        "tier_counts": tiers,
        "summary_csv": str(summary_csv),
        "summary_json": str(summary_json),
        "top_reverse_peak": top("reverse_peak"),
        "top_reverse_excess_integral": top("reverse_excess_integral_above_low_baseline"),
        "top_energy_proxy": top("u3_energy_proxy"),
        "lowest_energy_proxy": top("u3_energy_proxy", reverse=False),
        "largest_dynamic_steady_gap": top("max_dynamic_minus_steady_any_selected"),
    }
    report_json = args.output_dir / "sweep_report.json"
    report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"n_scenarios": len(summaries), "tier_counts": tiers, "summary_csv": str(summary_csv), "report_json": str(report_json)}, indent=2))


if __name__ == "__main__":
    main()
