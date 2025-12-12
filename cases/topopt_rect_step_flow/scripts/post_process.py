#!/usr/bin/env python3
"""
Post-process key scalar metrics from OpenFOAM functionObjects output.

This script expects the functionObjects in system/controlDict:
  - inletPressure / outletPressure : areaAverage(p)
  - inletFlowRate / outletFlowRate : sum(phi)

It writes:
  - results_summary.txt
"""

from __future__ import annotations

import argparse
from pathlib import Path


# Case geometry/physics (SI units, after scaling mm->m)
H = 0.41e-3   # channel height [m]
T = 10.0e-3   # channel thickness [m]
U_MEAN = 0.3  # inlet mean velocity [m/s]
NU = 1e-3     # kinematic viscosity [m^2/s] (see constant/transportProperties)


def _case_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _latest_time_dir(base: Path) -> Path | None:
    if not base.exists():
        return None
    candidates = []
    for d in base.iterdir():
        if not d.is_dir():
            continue
        try:
            t = float(d.name)
        except ValueError:
            continue
        candidates.append((t, d))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def _read_last_value_from_function_object(func_dir: Path) -> float | None:
    """
    Reads last numeric value from the latest time directory in postProcessing/<func>.
    Supports the common surfaceFieldValue output: surfaceFieldValue.dat
    """
    latest = _latest_time_dir(func_dir)
    if latest is None:
        return None

    dat_files = sorted(latest.glob("*.dat"))
    if not dat_files:
        return None

    dat = None
    for f in dat_files:
        if f.name == "surfaceFieldValue.dat":
            dat = f
            break
    if dat is None:
        dat = dat_files[0]

    lines = dat.read_text(encoding="utf-8", errors="ignore").splitlines()
    # Find last non-comment data line
    for line in reversed(lines):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            return float(parts[1])
        except ValueError:
            continue
    return None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Post-process TopOpt OpenFOAM case outputs.")
    p.add_argument("--output", type=str, default="results_summary.txt", help="Output summary filename")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = _case_dir()
    post_dir = case_dir / "postProcessing"

    if not post_dir.exists():
        print(f"[ERR] postProcessing not found: {post_dir}")
        print("      simpleFoam 실행 후 다시 시도하세요.")
        return 2

    pin = _read_last_value_from_function_object(post_dir / "inletPressure")
    pout = _read_last_value_from_function_object(post_dir / "outletPressure")
    qin = _read_last_value_from_function_object(post_dir / "inletFlowRate")
    qout = _read_last_value_from_function_object(post_dir / "outletFlowRate")

    area = H * T
    re = U_MEAN * H / NU

    lines = []
    lines.append("TopOpt STEP -> OpenFOAM results summary")
    lines.append("")
    lines.append(f"Geometry (SI): H={H:.6e} m, T={T:.6e} m, A=H*T={area:.6e} m^2")
    lines.append(f"Reference: Umean={U_MEAN:.6g} m/s, nu={NU:.6g} m^2/s, Re=Umean*H/nu={re:.6g}")
    lines.append("")

    if pin is not None and pout is not None:
        dp = pin - pout
        lines.append(f"p_inlet_avg   = {pin:.6e} Pa")
        lines.append(f"p_outlet_avg  = {pout:.6e} Pa")
        lines.append(f"delta_p       = {dp:.6e} Pa")
    else:
        lines.append("Pressure: (missing) - check functionObjects output (inletPressure/outletPressure)")

    lines.append("")

    if qin is not None and qout is not None:
        # phi sum sign depends on face normals; report both raw and absolute
        q_avg = 0.5 * (abs(qin) + abs(qout))
        u_est = q_avg / area if area > 0 else float("nan")
        lines.append(f"Q_inlet_sum(phi)  = {qin:.6e} m^3/s (raw)")
        lines.append(f"Q_outlet_sum(phi) = {qout:.6e} m^3/s (raw)")
        lines.append(f"Q_avg_abs         = {q_avg:.6e} m^3/s")
        lines.append(f"Umean_est=Q/A      = {u_est:.6e} m/s")
    else:
        lines.append("Flow rate: (missing) - check functionObjects output (inletFlowRate/outletFlowRate)")

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = case_dir / out_path

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[OK] Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


