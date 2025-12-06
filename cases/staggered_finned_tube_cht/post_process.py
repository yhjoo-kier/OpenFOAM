import argparse
import math
from pathlib import Path
from typing import Dict, Optional

DEFAULT_GEOM = {
    "SL": 0.060,
    "ST": 0.060,
    "P_fin": 0.050,
    "R_tube_out": 0.012,
    "R_fin": 0.022,
    "t_fin": 0.001,
}

DEFAULT_PROPS = {
    "rho": 1.18,  # kg/m3
    "k": 0.026,  # W/m/K
}


def latest_time_dir(base: Path) -> Optional[Path]:
    if not base.exists():
        return None
    time_dirs = []
    for item in base.iterdir():
        try:
            time_dirs.append((float(item.name), item))
        except ValueError:
            continue
    if not time_dirs:
        return None
    return sorted(time_dirs, key=lambda x: x[0])[-1][1]


def read_surface_field_value(case_dir: Path, func_name: str) -> Optional[float]:
    base_dir = case_dir / "postProcessing"
    if not base_dir.exists():
        return None

    candidates = [base_dir / func_name]
    for region_dir in base_dir.iterdir():
        if region_dir.is_dir():
            candidate = region_dir / func_name
            candidates.append(candidate)

    data_file: Optional[Path] = None
    for func_dir in candidates:
        time_dir = latest_time_dir(func_dir)
        if time_dir is None:
            continue
        candidate_file = time_dir / "surfaceFieldValue.dat"
        if candidate_file.exists():
            data_file = candidate_file
            break

    if data_file is None:
        return None

    value = None
    for line in data_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        try:
            value = float(parts[-1])
        except (ValueError, IndexError):
            continue
    return value


def hydraulic_diameter(params: Dict[str, float]) -> float:
    flow_area = params["ST"] * params["P_fin"] - math.pi * params["R_fin"]**2
    wetted_perimeter = 2 * (params["ST"] + params["P_fin"]) + 2 * math.pi * params["R_fin"]
    if wetted_perimeter <= 0 or flow_area <= 0:
        raise ValueError("Invalid geometry for hydraulic diameter calculation.")
    return 4 * flow_area / wetted_perimeter


def interface_area(params: Dict[str, float]) -> float:
    tube_area = 2 * math.pi * params["R_tube_out"] * params["P_fin"]
    fin_cyl_area = 2 * math.pi * params["R_fin"] * params["t_fin"]
    fin_annulus_area = 2 * math.pi * (params["R_fin"]**2 - params["R_tube_out"]**2)
    return tube_area + fin_cyl_area + fin_annulus_area


def compute_lmtd(T_hot_in: float, T_hot_out: float, T_wall: float) -> Optional[float]:
    dT1 = T_wall - T_hot_in
    dT2 = T_wall - T_hot_out
    if dT1 == dT2 or dT1 == 0 or dT2 == 0:
        return None
    try:
        return (dT1 - dT2) / math.log(dT1 / dT2)
    except (ValueError, ZeroDivisionError):
        return None


def compute_metrics(case_dir: Path, params: Dict[str, float], props: Dict[str, float], U_mean: float) -> Dict[str, float]:
    inlet_p = read_surface_field_value(case_dir, "inletPressureAverage")
    outlet_p = read_surface_field_value(case_dir, "outletPressureAverage")
    inlet_T = read_surface_field_value(case_dir, "inletTemperatureAverage")
    outlet_T = read_surface_field_value(case_dir, "outletTemperatureAverage")
    wall_T = read_surface_field_value(case_dir, "interfaceTemperatureAverage")
    q_total = read_surface_field_value(case_dir, "wallHeatFluxIntegral")

    results: Dict[str, float] = {}
    if inlet_p is not None and outlet_p is not None:
        results["delta_p"] = outlet_p - inlet_p

    if inlet_T is not None:
        results["T_inlet"] = inlet_T
    if outlet_T is not None:
        results["T_outlet"] = outlet_T
    if wall_T is not None:
        results["T_wall"] = wall_T
    if q_total is not None:
        results["Q"] = q_total

    dh = hydraulic_diameter(params)
    results["hydraulic_diameter"] = dh

    if all(v is not None for v in (inlet_T, outlet_T, wall_T)):
        lmtd = compute_lmtd(inlet_T, outlet_T, wall_T)
        if lmtd:
            area = interface_area(params)
            h_coeff = q_total / (area * lmtd) if q_total else None
            if h_coeff is not None:
                results["Nu"] = h_coeff * dh / props["k"]
                results["h"] = h_coeff
            results["LMTD"] = lmtd

    if "delta_p" in results:
        dp = results["delta_p"]
        results["friction_factor"] = (dp / params["SL"]) * (dh / (2 * props["rho"] * U_mean**2))

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Post-process staggered finned-tube CHT results.")
    parser.add_argument("case", type=Path, help="Case directory containing postProcessing outputs")
    parser.add_argument("--Umean", type=float, default=1.0, help="Bulk mean velocity used in the simulation (m/s)")
    args = parser.parse_args()

    try:
        metrics = compute_metrics(args.case, DEFAULT_GEOM, DEFAULT_PROPS, args.Umean)
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"Post-processing failed: {exc}")

    if not metrics:
        raise SystemExit("No postProcessing data found. Run the solver first.")

    print("--- Post-processing summary ---")
    for key, value in metrics.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
