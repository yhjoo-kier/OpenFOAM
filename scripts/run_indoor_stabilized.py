#!/usr/bin/env python3
"""Run indoor CFD pipeline with automatic stabilization retries.

Flow:
1. Generate scene JSON with Gemini
2. Mesh with Gmsh
3. Create OpenFOAM case
4. Try solver presets using mesh-risk-aware ordering
5. If needed, repair/simplify scene and retry
6. Export side-by-side visualization on success
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"
TOOLS_PYTHON = "/home/yhjoo/.venvs/openclaw-tools/bin/python"

PRESETS = [
    {
        "name": "baseline",
        "mode": "RAS",
        "inlet_velocity": 0.5,
        "internal_velocity": 0.5,
        "k": 9.375e-4,
        "omega": 0.55901699,
        "p_relax": 0.3,
        "u_relax": 0.7,
        "k_relax": 0.7,
        "omega_relax": 0.7,
        "nNonOrthogonalCorrectors": 1,
    },
    {
        "name": "conservative",
        "mode": "RAS",
        "inlet_velocity": 0.05,
        "internal_velocity": 0.0,
        "k": 9.375e-6,
        "omega": 0.55901699,
        "p_relax": 0.2,
        "u_relax": 0.3,
        "k_relax": 0.3,
        "omega_relax": 0.3,
        "nNonOrthogonalCorrectors": 3,
    },
    {
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
    },
    {
        "name": "ultra_robust",
        "mode": "RAS",
        "inlet_velocity": 0.005,
        "internal_velocity": 0.0,
        "k": 1.0e-7,
        "omega": 0.1,
        "p_relax": 0.12,
        "u_relax": 0.15,
        "k_relax": 0.15,
        "omega_relax": 0.15,
        "nNonOrthogonalCorrectors": 5,
    },
    {
        "name": "laminar_fallback",
        "mode": "laminar",
        "inlet_velocity": 0.01,
        "internal_velocity": 0.0,
        "k": 1.0e-8,
        "omega": 0.05,
        "p_relax": 0.15,
        "u_relax": 0.2,
        "k_relax": 0.2,
        "omega_relax": 0.2,
        "nNonOrthogonalCorrectors": 5,
    },
]

MESH_SIZE_LADDER = [0.35, 0.25, 0.18]


def extract_default_faces(import_log_text: str) -> int:
    m = re.search(r"Found\s+(\d+)\s+undefined faces in mesh", import_log_text)
    return int(m.group(1)) if m else 0


def parse_checkmesh_metrics(checkmesh_text: str) -> dict:
    metrics = {
        "max_non_ortho": None,
        "avg_non_ortho": None,
        "severe_non_ortho_faces": None,
        "max_aspect_ratio": None,
        "max_skewness": None,
        "mesh_ok": "Mesh OK." in checkmesh_text,
    }
    m = re.search(r"Mesh non-orthogonality Max:\s*([0-9eE+\-.]+)\s+average:\s*([0-9eE+\-.]+)", checkmesh_text)
    if m:
        metrics["max_non_ortho"] = float(m.group(1))
        metrics["avg_non_ortho"] = float(m.group(2))
    m = re.search(r"Number of severely non-orthogonal \(> 70 degrees\) faces:\s*(\d+)", checkmesh_text)
    if m:
        metrics["severe_non_ortho_faces"] = int(m.group(1))
    m = re.search(r"Max aspect ratio =\s*([0-9eE+\-.]+)", checkmesh_text)
    if m:
        metrics["max_aspect_ratio"] = float(m.group(1))
    m = re.search(r"Max skewness =\s*([0-9eE+\-.]+)", checkmesh_text)
    if m:
        metrics["max_skewness"] = float(m.group(1))
    return metrics


def classify_mesh_risk(default_faces: int, checkmesh_metrics: dict) -> dict:
    score = 0
    reasons: list[str] = []
    max_non_ortho = checkmesh_metrics.get("max_non_ortho")
    severe_non_ortho = checkmesh_metrics.get("severe_non_ortho_faces")
    max_aspect_ratio = checkmesh_metrics.get("max_aspect_ratio")
    max_skewness = checkmesh_metrics.get("max_skewness")

    if default_faces > 0:
        score += 3
        reasons.append(f"defaultFaces={default_faces}")
    if max_non_ortho is not None:
        if max_non_ortho > 85:
            score += 3
            reasons.append(f"maxNonOrtho={max_non_ortho:.2f}>85")
        elif max_non_ortho > 80:
            score += 2
            reasons.append(f"maxNonOrtho={max_non_ortho:.2f}>80")
        elif max_non_ortho > 70:
            score += 1
            reasons.append(f"maxNonOrtho={max_non_ortho:.2f}>70")
    if severe_non_ortho is not None:
        if severe_non_ortho > 100:
            score += 2
            reasons.append(f"severeNonOrthoFaces={severe_non_ortho}>100")
        elif severe_non_ortho > 0:
            score += 1
            reasons.append(f"severeNonOrthoFaces={severe_non_ortho}")
    if max_aspect_ratio is not None and max_aspect_ratio > 20:
        score += 1
        reasons.append(f"aspectRatio={max_aspect_ratio:.2f}>20")
    if max_skewness is not None and max_skewness > 0.85:
        score += 1
        reasons.append(f"skewness={max_skewness:.3f}>0.85")

    if score >= 5:
        level = "high"
    elif score >= 2:
        level = "medium"
    else:
        level = "low"
    return {"score": score, "level": level, "reasons": reasons}


def preset_order_for_risk(risk_level: str) -> list[dict]:
    by_name = {p["name"]: p for p in PRESETS}
    if risk_level == "high":
        order = ["robust", "ultra_robust", "laminar_fallback", "conservative", "baseline"]
    elif risk_level == "medium":
        order = ["conservative", "robust", "ultra_robust", "laminar_fallback", "baseline"]
    else:
        order = ["baseline", "conservative", "robust", "ultra_robust", "laminar_fallback"]
    return [by_name[name] for name in order]


def classify_failure(log_text: str, returncode: int | None = None) -> dict:
    lower = log_text.lower()
    failure_type = "unknown"
    hints: list[str] = []

    trapfpe_banner = "trapfpe: floating point exception trapping enabled"
    cleaned_lower = lower.replace("foam_sigfpe", "").replace(trapfpe_banner, "")
    has_real_sigfpe = any(token in cleaned_lower for token in [
        "foam::sigfpe::sighandler",
        "floating point exception (core dumped)",
    ]) or "floating point exception" in cleaned_lower

    if has_real_sigfpe:
        failure_type = "floating_point_exception"
    elif "foam fatal error" in lower:
        failure_type = "foam_fatal_error"
    elif re.search(r"(^|[^a-z])nan([^a-z]|$)", lower):
        failure_type = "nan_detected"
    elif "sigterm" in lower:
        failure_type = "terminated"
    elif returncode not in (None, 0):
        failure_type = f"nonzero_exit_{returncode}"

    if "bounding omega" in lower:
        hints.append("omega_instability")
    if "bounding k" in lower:
        hints.append("k_instability")
    if "bounding u" in lower:
        hints.append("u_instability")
    if re.search(r"time step continuity errors\s*:.*e\+", lower):
        hints.append("continuity_explosion")
    if "solving for p" in lower and "no iterations 1000" in lower:
        hints.append("pressure_solver_stall")
    if returncode == 136:
        hints.append("sigfpe_exit_136")

    if any(h in hints for h in ["omega_instability", "k_instability", "u_instability"]):
        family = "turbulence_or_field_blowup"
    elif any(h in hints for h in ["continuity_explosion", "pressure_solver_stall"]):
        family = "pressure_velocity_coupling"
    elif failure_type in {"floating_point_exception", "foam_fatal_error", "nan_detected"} or failure_type.startswith("nonzero_exit_"):
        family = "generic_solver_failure"
    else:
        family = "completed_or_unknown"

    return {"failure_type": failure_type, "failure_family": family, "hints": hints}


def run(cmd: list[str], cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    print("+", " ".join(map(str, cmd)))
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if result.stdout.strip():
        print(result.stdout.rstrip())
    if result.stderr.strip():
        print(result.stderr.rstrip())
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}")
    return result


def extract_generation_summary(stdout_text: str) -> dict | None:
    marker = "---"
    if marker not in stdout_text:
        return None
    candidate = stdout_text.split(marker)[-1].strip()
    try:
        payload = json.loads(candidate)
    except json.JSONDecodeError:
        return None
    return payload if isinstance(payload, dict) else None


def docker_exec(case_dir: Path, command: str, check: bool = True, timeout: float | None = None) -> subprocess.CompletedProcess:
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{PROJECT_ROOT}:/app",
        "-w", f"/app/{case_dir.relative_to(PROJECT_ROOT)}",
        DOCKER_IMAGE,
        "bash", "-lc", command,
    ]
    print("+", " ".join(map(str, cmd)))
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=timeout)
    if result.stdout.strip():
        print(result.stdout.rstrip())
    if result.stderr.strip():
        print(result.stderr.rstrip())
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}")
    return result


def replace_in_file(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    if old not in text:
        raise RuntimeError(f"Could not find target text in {path}: {old}")
    path.write_text(text.replace(old, new), encoding="utf-8")


def patch_boundary_types(case_dir: Path) -> None:
    boundary_path = case_dir / "constant" / "polyMesh" / "boundary"
    if not boundary_path.exists():
        return
    text = boundary_path.read_text(encoding="utf-8")
    text = text.replace(
        """    obstacleWalls
    {
        type            patch;
        physicalType    patch;""",
        """    obstacleWalls
    {
        type            wall;
        physicalType    wall;""",
    )
    text = text.replace(
        """    roomWalls
    {
        type            patch;
        physicalType    patch;""",
        """    roomWalls
    {
        type            wall;
        physicalType    wall;""",
    )
    boundary_path.write_text(text, encoding="utf-8")


def apply_turbulence_mode(case_dir: Path, mode: str) -> None:
    turb_path = case_dir / "constant" / "turbulenceProperties"
    if mode == "laminar":
        turb_path.write_text(
            """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  laminar;
""",
            encoding="utf-8",
        )
    else:
        turb_path.write_text(
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
""",
            encoding="utf-8",
        )


def apply_preset(case_dir: Path, preset: dict) -> None:
    u_path = case_dir / "0" / "U"
    p_path = case_dir / "0" / "p"
    k_path = case_dir / "0" / "k"
    omega_path = case_dir / "0" / "omega"
    nut_path = case_dir / "0" / "nut"
    fvsolution_path = case_dir / "system" / "fvSolution"

    apply_turbulence_mode(case_dir, preset.get("mode", "RAS"))

    write_u = f"""FoamFile
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
"""
    u_path.write_text(write_u, encoding="utf-8")

    p_text = p_path.read_text(encoding="utf-8").replace("type            empty;", "type            zeroGradient;")
    p_path.write_text(p_text, encoding="utf-8")

    k_text = f"""FoamFile
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
"""
    k_path.write_text(k_text, encoding="utf-8")

    omega_text = f"""FoamFile
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
"""
    omega_path.write_text(omega_text, encoding="utf-8")

    nut_text = nut_path.read_text(encoding="utf-8").replace(
        """    defaultFaces
    {
        type            empty;
    }""",
        """    defaultFaces
    {
        type            calculated;
        value           uniform 0;
    }"""
    )
    nut_path.write_text(nut_text, encoding="utf-8")

    fv = fvsolution_path.read_text(encoding="utf-8")
    import re
    fv = re.sub(r"nNonOrthogonalCorrectors\s+\d+;", f"nNonOrthogonalCorrectors {preset['nNonOrthogonalCorrectors']};", fv)
    fv = re.sub(r"relTol\s+0\.05;", "relTol          0.01;", fv, count=1)
    fv = re.sub(r"p\s+0\.[0-9]+;", f"p               {preset['p_relax']};", fv, count=1)
    fv = re.sub(r"U\s+0\.[0-9]+;", f"U               {preset['u_relax']};", fv, count=1)
    fv = re.sub(r"k\s+0\.[0-9]+;", f"k               {preset['k_relax']};", fv, count=1)
    fv = re.sub(r"omega\s+0\.[0-9]+;", f"omega           {preset['omega_relax']};", fv, count=1)
    fvsolution_path.write_text(fv, encoding="utf-8")


def detect_failure(log_text: str, returncode: int | None = None) -> tuple[bool, str]:
    if detect_success(log_text) and returncode in (None, 0):
        return False, "completed"
    info = classify_failure(log_text, returncode=returncode)
    if info["failure_type"] != "unknown":
        return True, info["failure_type"]
    return returncode not in (None, 0), "completed" if returncode in (None, 0) else f"nonzero_exit_{returncode}"


def detect_success(log_text: str) -> bool:
    lower = log_text.lower()
    has_end = "\nend\n" in lower or lower.rstrip().endswith("end")
    reached_target_time = re.search(r"\ntime\s*=\s*1000\s*$", lower, flags=re.MULTILINE) is not None
    return has_end or reached_target_time


def main() -> int:
    parser = argparse.ArgumentParser(description="Run indoor CFD pipeline with stabilization retries")
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--name", default="indoor_stabilized_run")
    parser.add_argument("--backend", choices=["cli", "api"], default="api")
    parser.add_argument("--model", default="gemini-3.1-pro-preview")
    parser.add_argument("--image", dest="images", action="append", default=[], help="Optional image/photo/rendering path; repeatable")
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--skip-mesh-ladder", action="store_true")
    parser.add_argument("--end-time", type=int, default=1000)
    parser.add_argument("--solver-timeout", type=int, default=900, help="Timeout in seconds for each simpleFoam attempt")
    parser.add_argument("--import-timeout", type=int, default=180, help="Timeout in seconds for mesh import (gmshToFoam)")
    parser.add_argument("--checkmesh-timeout", type=int, default=180, help="Timeout in seconds for checkMesh")
    parser.add_argument("--no-fallback", action="store_true")
    parser.add_argument("--disable-repair", action="store_true", help="Disable repaired-scene retry")
    args = parser.parse_args()

    generated_dir = PROJECT_ROOT / "generated"
    case_dir = PROJECT_ROOT / "cases" / args.name
    results_dir = PROJECT_ROOT / "results" / args.name
    generated_dir.mkdir(exist_ok=True, parents=True)
    results_dir.mkdir(exist_ok=True, parents=True)

    scene_json = generated_dir / f"{args.name}.json"
    repaired_scene_json = generated_dir / f"{args.name}_repaired.json"
    geo_path = generated_dir / f"{args.name}.geo_unrolled"
    msh_path = generated_dir / f"{args.name}.msh"

    scenario_path = Path(args.scenario)
    try:
        scenario_is_json_file = scenario_path.exists() and scenario_path.suffix.lower() == ".json"
    except OSError:
        scenario_is_json_file = False
    input_source_type = "scene_json" if scenario_is_json_file else ("image" if args.images else "text")
    generation_invocation = None
    if scenario_is_json_file:
        shutil.copyfile(scenario_path, scene_json)
        print(f"Using existing scene JSON without Gemini regeneration: {scenario_path} -> {scene_json}")
    else:
        gen_cmd = [
            "python3", str(SCRIPTS / "generate_indoor_scene_with_gemini.py"),
            "--backend", args.backend,
            "--scenario", args.scenario,
            "--model", args.model,
            "-o", str(scene_json),
        ]
        for image_path in args.images:
            gen_cmd.extend(["--image", image_path])
        if args.no_fallback:
            gen_cmd.append("--no-fallback")
        generation_result = run(gen_cmd, cwd=PROJECT_ROOT, check=False)
        generation_invocation = {
            "command": gen_cmd,
            "stdout": generation_result.stdout,
            "stderr": generation_result.stderr,
            "returncode": generation_result.returncode,
        }
        if generation_result.returncode != 0:
            if scene_json.exists():
                print(
                    "Generation step returned nonzero, but scene JSON was still written; "
                    "continuing into repair/salvage path.",
                )
            else:
                raise RuntimeError(f"Generation failed before writing scene JSON: {' '.join(map(str, gen_cmd))}")

    gmsh_python = TOOLS_PYTHON if Path(TOOLS_PYTHON).exists() else shutil.which("python3") or "python3"
    mesh_sizes = [args.mesh_size] if args.skip_mesh_ladder else [args.mesh_size] + [m for m in MESH_SIZE_LADDER if m != args.mesh_size]

    attempts = []
    success_preset = None
    success_mesh_size = None
    used_scene_json = scene_json
    repair_info = None

    scene_candidates = [("original", scene_json)]
    if not args.disable_repair:
        try:
            run([
                "python3", str(SCRIPTS / "repair_indoor_scene.py"),
                str(scene_json), "-o", str(repaired_scene_json)
            ], cwd=PROJECT_ROOT)
            scene_candidates.append(("repaired", repaired_scene_json))
        except Exception as exc:
            repair_info = {
                "repair_attempted": True,
                "repair_available": False,
                "repair_error": str(exc),
            }
    else:
        repair_info = {
            "repair_attempted": False,
            "repair_available": False,
            "repair_error": None,
        }

    for scene_variant, active_scene_json in scene_candidates:
        for mesh_size in mesh_sizes:
            try:
                run([
                    gmsh_python, str(SCRIPTS / "scene_to_gmsh.py"), str(active_scene_json),
                    "-o", str(msh_path), "--geo", str(geo_path), "--mesh-size", str(mesh_size)
                ], cwd=PROJECT_ROOT)
            except Exception as exc:
                attempts.append({
                    "scene_variant": scene_variant,
                    "scene_json": str(active_scene_json),
                    "mesh_size": mesh_size,
                    "stage": "scene_to_gmsh",
                    "reason": "meshing_failed",
                    "failure_family": "geometry_or_meshing",
                    "failure_hints": [str(exc)],
                })
                continue

            try:
                if case_dir.exists():
                    shutil.rmtree(case_dir)
                run([
                    "python3", str(SCRIPTS / "create_indoor_openfoam_case.py"),
                    str(msh_path), str(case_dir),
                    "--end-time", str(args.end_time)
                ], cwd=PROJECT_ROOT)
            except Exception as exc:
                attempts.append({
                    "scene_variant": scene_variant,
                    "scene_json": str(active_scene_json),
                    "mesh_size": mesh_size,
                    "stage": "create_case",
                    "reason": "case_creation_failed",
                    "failure_family": "geometry_or_meshing",
                    "failure_hints": [str(exc)],
                })
                continue

            try:
                import_result = docker_exec(
                    case_dir,
                    f"rm -rf constant/polyMesh && cp /app/generated/{args.name}.msh . && gmshToFoam {args.name}.msh",
                    check=True,
                    timeout=args.import_timeout,
                )
                checkmesh_result = docker_exec(case_dir, "checkMesh", check=True, timeout=args.checkmesh_timeout)
            except Exception as exc:
                attempts.append({
                    "scene_variant": scene_variant,
                    "scene_json": str(active_scene_json),
                    "mesh_size": mesh_size,
                    "stage": "foam_import_or_checkmesh",
                    "reason": "foam_import_failed",
                    "failure_family": "mesh_import",
                    "failure_hints": [str(exc)],
                })
                continue
            patch_boundary_types(case_dir)

            import_log = (import_result.stdout or "") + "\n" + (import_result.stderr or "")
            checkmesh_log = (checkmesh_result.stdout or "") + "\n" + (checkmesh_result.stderr or "")
            default_faces = extract_default_faces(import_log)
            mesh_metrics = parse_checkmesh_metrics(checkmesh_log)
            mesh_risk = classify_mesh_risk(default_faces, mesh_metrics)
            preset_sequence = preset_order_for_risk(mesh_risk["level"])

            for preset in preset_sequence:
                apply_preset(case_dir, preset)
                docker_exec(case_dir, "rm -rf [1-9]* 0.* core* 2>/dev/null || true", check=False)
                result = docker_exec(case_dir, "simpleFoam > log.simpleFoam 2>&1", check=False, timeout=args.solver_timeout)
                log_path = case_dir / "log.simpleFoam"
                log_text = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else result.stdout + result.stderr
                archived_log = results_dir / f"log.simpleFoam_{scene_variant}_mesh{str(mesh_size).replace('.', 'p')}_{preset['name']}.txt"
                archived_log.write_text(log_text, encoding="utf-8")
                failed, reason = detect_failure(log_text, returncode=result.returncode)
                failure_info = classify_failure(log_text, returncode=result.returncode)
                succeeded = detect_success(log_text) and result.returncode == 0 and not failed
                attempts.append({
                    "scene_variant": scene_variant,
                    "scene_json": str(active_scene_json),
                    "mesh_size": mesh_size,
                    "preset": preset["name"],
                    "mode": preset.get("mode", "RAS"),
                    "returncode": result.returncode,
                    "reason": "completed" if succeeded else reason,
                    "failure_family": None if succeeded else failure_info["failure_family"],
                    "failure_hints": [] if succeeded else failure_info["hints"],
                    "mesh_risk_level": mesh_risk["level"],
                    "mesh_risk_score": mesh_risk["score"],
                    "mesh_risk_reasons": mesh_risk["reasons"],
                    "default_faces": default_faces,
                    "checkmesh": mesh_metrics,
                    "log": str(archived_log),
                })
                if succeeded:
                    success_preset = preset
                    success_mesh_size = mesh_size
                    used_scene_json = active_scene_json
                    break

            if success_preset is not None:
                break
        if success_preset is not None:
            break

    generation_summary = None if generation_invocation is None else extract_generation_summary(generation_invocation["stdout"])
    summary = {
        "name": args.name,
        "input_source_type": input_source_type,
        "requested_backend": args.backend,
        "requested_model": args.model,
        "input_images": args.images,
        "scene_json": str(scene_json),
        "repaired_scene_json": None if args.disable_repair else str(repaired_scene_json),
        "used_scene_json": str(used_scene_json),
        "msh": str(msh_path),
        "case_dir": str(case_dir),
        "generation_summary": generation_summary,
        "repair_info": repair_info,
        "attempts": attempts,
        "success": success_preset is not None,
        "successful_preset": None if success_preset is None else success_preset["name"],
        "successful_mode": None if success_preset is None else success_preset.get("mode", "RAS"),
        "successful_mesh_size": success_mesh_size,
        "successful_scene_variant": None if success_preset is None else ("repaired" if used_scene_json == repaired_scene_json else "original"),
        "compass_version": "2026-03-15_openfoam_solver_stabilization_compass",
    }
    summary_path = results_dir / "stabilization_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    if success_preset is not None:
        docker_exec(case_dir, "foamToVTK -latestTime", check=True)
        run([
            gmsh_python,
            str(SCRIPTS / "visualize_indoor_case.py"),
            str(used_scene_json),
            str(case_dir),
            "-o", str(results_dir / "comparison_1x2.png"),
        ], cwd=PROJECT_ROOT)
        run([
            gmsh_python,
            str(SCRIPTS / "render_indoor_pipeline_3d.py"),
            "--scene-json", str(used_scene_json),
            "--case-dir", str(case_dir),
            "--output", str(results_dir / "indoor_pipeline_3d_comparison.png"),
        ], cwd=PROJECT_ROOT)

    print(json.dumps(summary, indent=2))
    return 0 if success_preset is not None else 1


if __name__ == "__main__":
    raise SystemExit(main())
