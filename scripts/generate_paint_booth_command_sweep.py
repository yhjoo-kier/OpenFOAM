#!/usr/bin/env python3
"""Generate blower command-history CSVs for paint-booth dynamic surrogate sweeps.

The generated histories preserve the transient-validation convention used in this
project: a switch time itself remains the pre-change value, and the first sample
strictly after the switch time carries the new command.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_floats(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def scenario_id(low_u: float, high_u: float, hold_s: float, ramp_s: float, cycles: int) -> str:
    def fmt(x: float) -> str:
        return f"{x:.2f}".replace(".", "p")

    return f"u{fmt(low_u)}_u{fmt(high_u)}_hold{hold_s:g}s_ramp{ramp_s:g}s_cyc{cycles}"


def command_value(t: float, low_u: float, high_u: float, hold_s: float, ramp_s: float, cycles: int) -> float:
    """Return low/high/ramped command value at time t.

    Sequence:
      initial low hold, then for each cycle:
        ramp up -> high hold -> ramp down -> low hold

    Switch/ramp boundaries are right-open: at the exact boundary, the previous
    state is retained; the next sample changes. This mirrors the CFD t=0 initial
    field convention used by `simulate_paint_booth_dynamic_metrics.py`.
    """
    eps = 1e-12
    elapsed = hold_s
    if t <= elapsed + eps:
        return low_u

    for _ in range(cycles):
        # Ramp up over (elapsed, elapsed + ramp_s]
        if ramp_s > 0:
            if t <= elapsed + ramp_s + eps:
                frac = max(0.0, min(1.0, (t - elapsed) / ramp_s))
                return low_u + frac * (high_u - low_u)
            elapsed += ramp_s
        else:
            if t <= elapsed + eps:
                return low_u

        # High hold over (elapsed, elapsed + hold_s]
        if t <= elapsed + hold_s + eps:
            return high_u
        elapsed += hold_s

        # Ramp down over (elapsed, elapsed + ramp_s]
        if ramp_s > 0:
            if t <= elapsed + ramp_s + eps:
                frac = max(0.0, min(1.0, (t - elapsed) / ramp_s))
                return high_u + frac * (low_u - high_u)
            elapsed += ramp_s
        else:
            if t <= elapsed + eps:
                return high_u

        # Low hold over (elapsed, elapsed + hold_s]
        if t <= elapsed + hold_s + eps:
            return low_u
        elapsed += hold_s

    return low_u


def scenario_duration(hold_s: float, ramp_s: float, cycles: int) -> float:
    return hold_s + cycles * (2.0 * hold_s + 2.0 * ramp_s)


def write_command_csv(path: Path, low_u: float, high_u: float, hold_s: float, ramp_s: float, cycles: int, dt: float) -> int:
    end_time = scenario_duration(hold_s, ramp_s, cycles)
    n = int(round(end_time / dt))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["time_s", "U_supply_mps"])
        writer.writeheader()
        for i in range(n + 1):
            t = round(i * dt, 10)
            u = command_value(t, low_u, high_u, hold_s, ramp_s, cycles)
            writer.writerow({"time_s": f"{t:.6g}", "U_supply_mps": f"{u:.6g}"})
    return n + 1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--low-u", type=float, default=4.36)
    ap.add_argument("--high-u-list", default="5.45,6.00,6.54")
    ap.add_argument("--hold-times", default="5,10,15,30")
    ap.add_argument("--ramp-times", default="0,2,5,10")
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--dt", type=float, default=0.25)
    args = ap.parse_args()

    scenario_dir = args.output_dir / "scenarios"
    high_values = parse_floats(args.high_u_list)
    hold_values = parse_floats(args.hold_times)
    ramp_values = parse_floats(args.ramp_times)

    manifest = []
    for high_u in high_values:
        for hold_s in hold_values:
            for ramp_s in ramp_values:
                sid = scenario_id(args.low_u, high_u, hold_s, ramp_s, args.cycles)
                csv_path = scenario_dir / f"{sid}.csv"
                n_rows = write_command_csv(csv_path, args.low_u, high_u, hold_s, ramp_s, args.cycles, args.dt)
                manifest.append(
                    {
                        "scenario_id": sid,
                        "command_csv": str(csv_path),
                        "low_u_mps": args.low_u,
                        "high_u_mps": high_u,
                        "hold_s": hold_s,
                        "ramp_s": ramp_s,
                        "cycles": args.cycles,
                        "dt_s": args.dt,
                        "duration_s": scenario_duration(hold_s, ramp_s, args.cycles),
                        "n_rows": n_rows,
                    }
                )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_json = args.output_dir / "scenario_manifest.json"
    manifest_csv = args.output_dir / "scenario_manifest.csv"
    manifest_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    with manifest_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
        writer.writeheader()
        writer.writerows(manifest)
    print(json.dumps({"n_scenarios": len(manifest), "manifest_json": str(manifest_json), "manifest_csv": str(manifest_csv)}, indent=2))


if __name__ == "__main__":
    main()
