#!/usr/bin/env python3
"""Grid independence study for the image-to-CFD benchmark.

For 3 representative cases, re-run the REFERENCE CFD at multiple mesh
resolutions (0.35m, 0.25m, 0.18m, 0.10m) and compare:
  1. Self-convergence: how flow metrics change as mesh is refined
  2. Predicted-vs-reference: how the CFD score changes when the reference
     is computed on a finer mesh

Cases: bench_a1_01, bench_a3_03, bench_a4_03 (all floorplan view)
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"
TOOLS_PYTHON = "/home/yhjoo/.venvs/openclaw-tools/bin/python"

CASES = [
    {"name": "bench_a1_01", "scene": "a1_01.json"},
    {"name": "bench_a3_03", "scene": "a3_03.json"},
    {"name": "bench_a4_03", "scene": "a4_03.json"},
]

MESH_SIZES = [0.35, 0.25, 0.18, 0.10]

# Use the "robust" preset — same as the successful reference runs
PRESET = {
    "name": "robust",
    "mode": "RAS",
    "inlet_velocity": 0.02,
    "internal_velocity": 0.0,
    "k": 1.5e-6,
    "omega": 0.35,
    "p_relax": 0.15,
    "u_relax": 0.2,
    "k_relax": 0.2,
    "omega_relax": 0.2,
    "nNonOrthogonalCorrectors": 4,
}

END_TIME = 1000
SOLVER_TIMEOUT = 1200  # 20 min per solve


def run_cmd(cmd: list[str], cwd: Path | None = None, check: bool = True,
            timeout: float | None = None) -> subprocess.CompletedProcess:
    print(f"  + {' '.join(map(str, cmd[:8]))}{'...' if len(cmd) > 8 else ''}")
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    if check and result.returncode != 0:
        print(f"  STDERR: {result.stderr[:500]}")
        raise RuntimeError(f"Command failed ({result.returncode})")
    return result


def docker_exec(case_dir: Path, command: str, check: bool = True,
                timeout: float | None = None) -> subprocess.CompletedProcess:
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{PROJECT_ROOT}:/app",
        "-w", f"/app/{case_dir.relative_to(PROJECT_ROOT)}",
        DOCKER_IMAGE,
        "bash", "-lc", command,
    ]
    print(f"  + docker: {command[:80]}{'...' if len(command) > 80 else ''}")
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=timeout)
    if check and result.returncode != 0:
        print(f"  STDERR: {result.stderr[:500]}")
        raise RuntimeError(f"Docker command failed ({result.returncode}): {command[:60]}")
    return result


def apply_preset_and_boundary(case_dir: Path) -> None:
    """Apply robust preset and fix boundary types — mirrors run_indoor_stabilized.py logic."""
    import re as _re

    preset = PRESET

    # turbulenceProperties
    (case_dir / "constant" / "turbulenceProperties").write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
""", encoding="utf-8")

    # U
    (case_dir / "0" / "U").write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({preset['internal_velocity']} 0 0);
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({preset['inlet_velocity']} 0 0);
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform ({preset['inlet_velocity']} 0 0);
        value           uniform (0 0 0);
    }}

    roomWalls
    {{
        type            noSlip;
    }}

    obstacleWalls
    {{
        type            noSlip;
    }}

    defaultFaces
    {{
        type            noSlip;
    }}
}}
""", encoding="utf-8")

    # k
    (case_dir / "0" / "k").write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {preset['k']};
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {preset['k']};
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {preset['k']};
        value           uniform {preset['k']};
    }}

    roomWalls
    {{
        type            kqRWallFunction;
        value           uniform {preset['k']};
    }}

    obstacleWalls
    {{
        type            kqRWallFunction;
        value           uniform {preset['k']};
    }}

    defaultFaces
    {{
        type            calculated;
        value           uniform {preset['k']};
    }}
}}
""", encoding="utf-8")

    # omega
    (case_dir / "0" / "omega").write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}}

dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {preset['omega']};
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {preset['omega']};
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {preset['omega']};
        value           uniform {preset['omega']};
    }}

    roomWalls
    {{
        type            omegaWallFunction;
        value           uniform {preset['omega']};
    }}

    obstacleWalls
    {{
        type            omegaWallFunction;
        value           uniform {preset['omega']};
    }}

    defaultFaces
    {{
        type            calculated;
        value           uniform {preset['omega']};
    }}
}}
""", encoding="utf-8")

    # Patch nut defaultFaces
    nut_path = case_dir / "0" / "nut"
    if nut_path.exists():
        txt = nut_path.read_text(encoding="utf-8")
        txt = txt.replace(
            """    defaultFaces
    {
        type            empty;
    }""",
            """    defaultFaces
    {
        type            calculated;
        value           uniform 0;
    }""")
        nut_path.write_text(txt, encoding="utf-8")

    # Patch p defaultFaces
    p_path = case_dir / "0" / "p"
    if p_path.exists():
        txt = p_path.read_text(encoding="utf-8")
        txt = txt.replace("type            empty;", "type            zeroGradient;")
        p_path.write_text(txt, encoding="utf-8")

    # Patch fvSolution relaxation
    fvs_path = case_dir / "system" / "fvSolution"
    if fvs_path.exists():
        fv = fvs_path.read_text(encoding="utf-8")
        fv = _re.sub(r"nNonOrthogonalCorrectors\s+\d+;",
                      f"nNonOrthogonalCorrectors {preset['nNonOrthogonalCorrectors']};", fv)
        fv = _re.sub(r"relTol\s+0\.05;", "relTol          0.01;", fv, count=1)
        fv = _re.sub(r"p\s+0\.[0-9]+;", f"p               {preset['p_relax']};", fv, count=1)
        fv = _re.sub(r"U\s+0\.[0-9]+;", f"U               {preset['u_relax']};", fv, count=1)
        fv = _re.sub(r"k\s+0\.[0-9]+;", f"k               {preset['k_relax']};", fv, count=1)
        fv = _re.sub(r"omega\s+0\.[0-9]+;", f"omega           {preset['omega_relax']};", fv, count=1)
        fvs_path.write_text(fv, encoding="utf-8")

    # Fix boundary types (wall instead of patch for roomWalls/obstacleWalls)
    boundary_path = case_dir / "constant" / "polyMesh" / "boundary"
    if boundary_path.exists():
        txt = boundary_path.read_text(encoding="utf-8")
        for name in ("obstacleWalls", "roomWalls"):
            txt = txt.replace(
                f"""    {name}\n    {{\n        type            patch;\n        physicalType    patch;""",
                f"""    {name}\n    {{\n        type            wall;\n        physicalType    wall;""")
        boundary_path.write_text(txt, encoding="utf-8")


def detect_success(log_text: str) -> bool:
    import re as _re
    lower = log_text.lower()
    has_end = "\nend\n" in lower or lower.rstrip().endswith("end")
    reached_target = _re.search(r"\ntime\s*=\s*1000\s*$", lower, flags=_re.MULTILINE) is not None
    return has_end or reached_target


def extract_flow_stats_from_vtk(case_dir: Path, scene_json: Path) -> dict | None:
    """Extract mean velocity magnitude from VTK using pyvista, sampling on normalized grid."""
    try:
        import numpy as np
        import pyvista as pv
        pv.OFF_SCREEN = True

        scene = json.loads(scene_json.read_text(encoding="utf-8"))
        room = scene.get("room", {})
        if "blocks" in room:
            blocks = room["blocks"]
            bounds = {
                "Lx": max(b["origin"]["x"] + b["size"]["dx"] for b in blocks),
                "Ly": max(b["origin"]["y"] + b["size"]["dy"] for b in blocks),
                "Lz": max(b["origin"]["z"] + b["size"]["dz"] for b in blocks),
            }
        else:
            size = room.get("size", {})
            bounds = {"Lx": size.get("Lx", 0), "Ly": size.get("Ly", 0), "Lz": size.get("Lz", 0)}

        # Find latest VTK
        vtk_dir = case_dir / "VTK"
        candidates = []
        for child in vtk_dir.iterdir():
            if not child.is_dir():
                continue
            try:
                ts = int(child.name.split("_")[-1])
            except ValueError:
                continue
            vtu = child / "internal.vtu"
            if vtu.exists():
                candidates.append((ts, vtu))
        if not candidates:
            return None
        candidates.sort(key=lambda x: x[0])
        vtk_path = candidates[-1][1]

        mesh = pv.read(vtk_path)
        if "U" not in mesh.point_data and "U" in mesh.cell_data:
            mesh = mesh.cell_data_to_point_data()
        if "U" not in mesh.point_data:
            return None

        u = np.asarray(mesh.point_data["U"], dtype=float)
        umag = np.linalg.norm(u, axis=1)

        # Sample on 18x18x10 normalized grid (same as benchmark)
        eps = (0.05, 0.05, 0.08)
        xs = np.linspace(eps[0], 1.0 - eps[0], 18)
        ys = np.linspace(eps[1], 1.0 - eps[1], 18)
        zs = np.linspace(eps[2], 1.0 - eps[2], 10)
        norm_pts = np.array([(x, y, z) for z in zs for y in ys for x in xs], dtype=float)
        phys_pts = np.empty_like(norm_pts)
        phys_pts[:, 0] = norm_pts[:, 0] * bounds["Lx"]
        phys_pts[:, 1] = norm_pts[:, 1] * bounds["Ly"]
        phys_pts[:, 2] = norm_pts[:, 2] * bounds["Lz"]

        sampled = pv.PolyData(phys_pts).sample(mesh)
        valid = np.asarray(sampled["vtkValidPointMask"], dtype=bool)
        if valid.sum() == 0:
            return None

        sampled_u = np.asarray(sampled["U"], dtype=float)[valid]
        sampled_umag = np.linalg.norm(sampled_u, axis=1)

        if "p" in sampled.array_names:
            sampled_p = np.asarray(sampled["p"], dtype=float).reshape(-1)[valid]
            p_mean = float(np.mean(sampled_p))
            p_rms = float(np.sqrt(np.mean(sampled_p ** 2)))
        else:
            p_mean = None
            p_rms = None

        return {
            "mean_Umag": round(float(np.mean(sampled_umag)), 8),
            "rms_Umag": round(float(np.sqrt(np.mean(sampled_umag ** 2))), 8),
            "max_Umag": round(float(np.max(sampled_umag)), 8),
            "mean_p": round(p_mean, 8) if p_mean is not None else None,
            "rms_p": round(p_rms, 8) if p_rms is not None else None,
            "valid_sample_count": int(valid.sum()),
            "total_sample_points": int(len(norm_pts)),
            "cell_count_approx": int(mesh.n_cells),
        }
    except Exception as exc:
        print(f"  Warning: VTK extraction failed: {exc}")
        return None


def run_single_case_mesh(case_name: str, scene_file: str, mesh_size: float,
                         work_dir: Path) -> dict:
    """Mesh + solve one case at one mesh size. Returns status dict."""
    scene_json = PROJECT_ROOT / "benchmark" / "scenes" / scene_file
    case_dir = work_dir / f"{case_name}_ms{str(mesh_size).replace('.', 'p')}"
    geo_path = work_dir / f"{case_name}_ms{str(mesh_size).replace('.', 'p')}.geo_unrolled"
    msh_path = work_dir / f"{case_name}_ms{str(mesh_size).replace('.', 'p')}.msh"

    result_entry = {
        "case": case_name,
        "mesh_size": mesh_size,
        "case_dir": str(case_dir),
        "success": False,
        "stage_failed": None,
        "error": None,
        "flow_stats": None,
        "cell_count": None,
    }

    gmsh_python = TOOLS_PYTHON if Path(TOOLS_PYTHON).exists() else shutil.which("python3") or "python3"

    # 1. Gmsh meshing
    print(f"\n[{case_name} @ {mesh_size}m] Step 1: Generating mesh...")
    try:
        run_cmd([
            gmsh_python, str(SCRIPTS / "scene_to_gmsh.py"), str(scene_json),
            "-o", str(msh_path), "--geo", str(geo_path), "--mesh-size", str(mesh_size)
        ], cwd=PROJECT_ROOT, timeout=300)
    except Exception as exc:
        result_entry["stage_failed"] = "scene_to_gmsh"
        result_entry["error"] = str(exc)[:200]
        return result_entry

    # 2. Create OpenFOAM case
    print(f"[{case_name} @ {mesh_size}m] Step 2: Creating OpenFOAM case...")
    try:
        if case_dir.exists():
            shutil.rmtree(case_dir)
        run_cmd([
            "python3", str(SCRIPTS / "create_indoor_openfoam_case.py"),
            str(msh_path), str(case_dir),
            "--end-time", str(END_TIME)
        ], cwd=PROJECT_ROOT, timeout=120)
    except Exception as exc:
        result_entry["stage_failed"] = "create_case"
        result_entry["error"] = str(exc)[:200]
        return result_entry

    # 3. Import mesh via Docker (gmshToFoam + checkMesh)
    print(f"[{case_name} @ {mesh_size}m] Step 3: Importing mesh into OpenFOAM...")
    try:
        msh_rel = msh_path.relative_to(PROJECT_ROOT)
        docker_exec(
            case_dir,
            f"rm -rf constant/polyMesh && cp /app/{msh_rel} . && gmshToFoam {msh_path.name}",
            check=True, timeout=300)
        docker_exec(case_dir, "checkMesh", check=True, timeout=300)
    except Exception as exc:
        result_entry["stage_failed"] = "foam_import"
        result_entry["error"] = str(exc)[:200]
        return result_entry

    # 4. Apply preset and boundary fixes
    print(f"[{case_name} @ {mesh_size}m] Step 4: Applying solver preset...")
    try:
        apply_preset_and_boundary(case_dir)
    except Exception as exc:
        result_entry["stage_failed"] = "apply_preset"
        result_entry["error"] = str(exc)[:200]
        return result_entry

    # 5. Run simpleFoam
    print(f"[{case_name} @ {mesh_size}m] Step 5: Running simpleFoam (timeout={SOLVER_TIMEOUT}s)...")
    t0 = time.time()
    try:
        docker_exec(case_dir, "rm -rf [1-9]* 0.* core* 2>/dev/null || true", check=False)
        solve_result = docker_exec(
            case_dir, "simpleFoam > log.simpleFoam 2>&1",
            check=False, timeout=SOLVER_TIMEOUT)
    except subprocess.TimeoutExpired:
        result_entry["stage_failed"] = "solver_timeout"
        result_entry["error"] = f"simpleFoam timed out after {SOLVER_TIMEOUT}s"
        return result_entry
    except Exception as exc:
        result_entry["stage_failed"] = "solver_run"
        result_entry["error"] = str(exc)[:200]
        return result_entry

    elapsed = time.time() - t0
    log_path = case_dir / "log.simpleFoam"
    log_text = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else ""

    if not detect_success(log_text):
        result_entry["stage_failed"] = "solver_diverged"
        result_entry["error"] = f"simpleFoam did not converge (elapsed={elapsed:.0f}s, rc={solve_result.returncode})"
        # Still try to extract VTK if partial results exist
        print(f"  Solver did not fully converge, attempting VTK export anyway...")

    # 6. Export VTK
    print(f"[{case_name} @ {mesh_size}m] Step 6: Exporting VTK (elapsed={elapsed:.1f}s)...")
    try:
        docker_exec(case_dir, "foamToVTK -latestTime", check=True, timeout=300)
    except Exception as exc:
        if result_entry["stage_failed"] is None:
            result_entry["stage_failed"] = "vtk_export"
            result_entry["error"] = str(exc)[:200]
        return result_entry

    # 7. Extract flow statistics
    print(f"[{case_name} @ {mesh_size}m] Step 7: Extracting flow statistics...")
    stats = extract_flow_stats_from_vtk(case_dir, scene_json)
    if stats:
        result_entry["flow_stats"] = stats
        result_entry["cell_count"] = stats.get("cell_count_approx")

    if result_entry["stage_failed"] is None:
        result_entry["success"] = True
    elif stats is not None:
        # Partial success: solver didn't fully converge but we got VTK data
        result_entry["success"] = True
        result_entry["note"] = "Partial convergence but VTK data extracted"

    print(f"[{case_name} @ {mesh_size}m] Done. success={result_entry['success']}, "
          f"cells={result_entry.get('cell_count')}, "
          f"mean_Umag={stats['mean_Umag'] if stats else 'N/A'}")
    return result_entry


def compute_cfd_score(reference_scene: Path, reference_case: Path,
                      predicted_scene: Path, predicted_case: Path,
                      output_path: Path) -> dict | None:
    """Run compute_benchmark_cfd_metrics.py and return the result."""
    tools_py = TOOLS_PYTHON if Path(TOOLS_PYTHON).exists() else "python3"
    try:
        result = run_cmd([
            tools_py, str(SCRIPTS / "compute_benchmark_cfd_metrics.py"),
            "--reference-scene", str(reference_scene),
            "--reference-case", str(reference_case),
            "--predicted-scene", str(predicted_scene),
            "--predicted-case", str(predicted_case),
            "-o", str(output_path),
        ], cwd=PROJECT_ROOT, timeout=120, check=False)
        if output_path.exists():
            return json.loads(output_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"  Warning: CFD metrics computation failed: {exc}")
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Grid independence study for image-to-CFD benchmark")
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "benchmark" / "manifests" / "grid_independence_results.json")
    parser.add_argument("--work-dir", type=Path,
                        default=PROJECT_ROOT / "cases")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip mesh sizes where VTK already exists")
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    eval_base = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
    scenes_dir = PROJECT_ROOT / "benchmark" / "scenes"

    all_results: dict[str, dict] = {}

    for case_info in CASES:
        case_name = case_info["name"]
        scene_file = case_info["scene"]
        scene_json = scenes_dir / scene_file
        print(f"\n{'='*60}")
        print(f"CASE: {case_name}")
        print(f"{'='*60}")

        case_results = {
            "case": case_name,
            "scene": str(scene_json),
            "mesh_runs": {},
            "self_convergence": {},
            "predicted_vs_refined_reference": {},
        }

        # Run reference at each mesh size
        for mesh_size in MESH_SIZES:
            ms_key = str(mesh_size)
            case_dir = args.work_dir / f"{case_name}_ms{str(mesh_size).replace('.', 'p')}"

            # Check if already done
            vtk_dir = case_dir / "VTK"
            if args.skip_existing and vtk_dir.exists() and any(vtk_dir.iterdir()):
                print(f"\n[{case_name} @ {mesh_size}m] Skipping (VTK exists)")
                stats = extract_flow_stats_from_vtk(case_dir, scene_json)
                entry = {
                    "case": case_name,
                    "mesh_size": mesh_size,
                    "case_dir": str(case_dir),
                    "success": True,
                    "stage_failed": None,
                    "error": None,
                    "flow_stats": stats,
                    "cell_count": stats.get("cell_count_approx") if stats else None,
                    "note": "skipped (existing VTK)",
                }
                case_results["mesh_runs"][ms_key] = entry
                continue

            entry = run_single_case_mesh(case_name, scene_file, mesh_size, args.work_dir)
            case_results["mesh_runs"][ms_key] = entry

        # Compute self-convergence: compare flow stats across mesh levels
        stats_by_mesh = {}
        for ms_key, run_info in case_results["mesh_runs"].items():
            if run_info.get("flow_stats"):
                stats_by_mesh[float(ms_key)] = run_info["flow_stats"]

        if len(stats_by_mesh) >= 2:
            sorted_sizes = sorted(stats_by_mesh.keys())
            convergence_table = []
            for i in range(1, len(sorted_sizes)):
                prev_ms = sorted_sizes[i - 1]
                curr_ms = sorted_sizes[i]
                prev_umag = stats_by_mesh[prev_ms]["mean_Umag"]
                curr_umag = stats_by_mesh[curr_ms]["mean_Umag"]
                delta = abs(curr_umag - prev_umag)
                rel_change = delta / max(abs(prev_umag), 1e-12)
                convergence_table.append({
                    "from_mesh": prev_ms,
                    "to_mesh": curr_ms,
                    "from_cells": stats_by_mesh[prev_ms].get("cell_count_approx"),
                    "to_cells": stats_by_mesh[curr_ms].get("cell_count_approx"),
                    "mean_Umag_from": prev_umag,
                    "mean_Umag_to": curr_umag,
                    "abs_change": round(delta, 8),
                    "rel_change_pct": round(rel_change * 100, 4),
                })
            case_results["self_convergence"] = {
                "stats_by_mesh_size": {str(k): v for k, v in stats_by_mesh.items()},
                "convergence_steps": convergence_table,
            }

        # Compute predicted-vs-refined-reference CFD scores
        eval_dir = eval_base / case_name / "floorplan"
        predicted_scene = eval_dir / "predicted_scene.json"
        predicted_case_dir = None

        # Find the predicted case directory
        if (eval_dir / "predicted_case").exists():
            predicted_case_dir = eval_dir / "predicted_case"
        else:
            # Fallback: look in cases/
            pred_pattern = f"eval_posthoc_uniform_longest_span_v1_{case_name}_floorplan"
            candidate = PROJECT_ROOT / "cases" / pred_pattern
            if candidate.exists():
                predicted_case_dir = candidate

        if predicted_case_dir and predicted_scene.exists():
            for ms_key, run_info in case_results["mesh_runs"].items():
                if not run_info.get("success"):
                    continue
                ref_case_dir = Path(run_info["case_dir"])
                output_metrics = args.work_dir / f"cfd_metrics_{case_name}_pred_vs_ref_{ms_key.replace('.', 'p')}.json"
                print(f"\n[{case_name}] Computing CFD score: predicted vs reference@{ms_key}m...")
                metrics = compute_cfd_score(
                    scene_json, ref_case_dir,
                    predicted_scene, predicted_case_dir,
                    output_metrics,
                )
                if metrics and "aggregate_score" in metrics:
                    case_results["predicted_vs_refined_reference"][ms_key] = {
                        "cfd_score": metrics["aggregate_score"]["cfd_score"],
                        "components": metrics["aggregate_score"].get("components"),
                        "overlap_count": metrics.get("overlap", {}).get("count"),
                    }
                else:
                    case_results["predicted_vs_refined_reference"][ms_key] = {
                        "cfd_score": None,
                        "error": "metrics computation failed",
                    }
        else:
            print(f"\n[{case_name}] Predicted case not found, skipping CFD score comparison")

        all_results[case_name] = case_results

    # Build final summary
    output = {
        "study": "grid_independence",
        "description": "CFD solution convergence as reference mesh is refined",
        "mesh_sizes": MESH_SIZES,
        "preset": PRESET["name"],
        "cases": all_results,
        "summary_table": [],
    }

    # Build human-readable summary table
    for case_name, cr in all_results.items():
        for ms in MESH_SIZES:
            ms_key = str(ms)
            run_info = cr["mesh_runs"].get(ms_key, {})
            stats = run_info.get("flow_stats", {}) or {}
            cfd_entry = cr.get("predicted_vs_refined_reference", {}).get(ms_key, {})
            output["summary_table"].append({
                "case": case_name,
                "mesh_size": ms,
                "cell_count": stats.get("cell_count_approx"),
                "solver_success": run_info.get("success", False),
                "mean_Umag": stats.get("mean_Umag"),
                "rms_Umag": stats.get("rms_Umag"),
                "max_Umag": stats.get("max_Umag"),
                "cfd_score_vs_predicted": cfd_entry.get("cfd_score"),
            })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(f"\n{'='*60}")
    print(f"Results saved to: {args.output}")
    print(f"{'='*60}")

    # Print summary table
    print(f"\n{'Case':<15} {'Mesh(m)':<10} {'Cells':<10} {'OK?':<5} {'Mean|U|':<12} {'RMS|U|':<12} {'CFD Score':<10}")
    print("-" * 74)
    for row in output["summary_table"]:
        cells = row["cell_count"] or "N/A"
        ok = "Y" if row["solver_success"] else "N"
        mu = f"{row['mean_Umag']:.6f}" if row["mean_Umag"] is not None else "N/A"
        ru = f"{row['rms_Umag']:.6f}" if row["rms_Umag"] is not None else "N/A"
        cs = f"{row['cfd_score_vs_predicted']:.4f}" if row["cfd_score_vs_predicted"] is not None else "N/A"
        print(f"{row['case']:<15} {row['mesh_size']:<10} {str(cells):<10} {ok:<5} {mu:<12} {ru:<12} {cs:<10}")

    # Print convergence
    print(f"\nSelf-convergence (mean |U| change between successive mesh refinements):")
    print(f"{'Case':<15} {'From':<8} {'To':<8} {'Abs Change':<12} {'Rel Change %':<12}")
    print("-" * 55)
    for case_name, cr in all_results.items():
        for step in cr.get("self_convergence", {}).get("convergence_steps", []):
            print(f"{case_name:<15} {step['from_mesh']:<8} {step['to_mesh']:<8} "
                  f"{step['abs_change']:<12.8f} {step['rel_change_pct']:<12.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
