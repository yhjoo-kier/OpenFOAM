#!/usr/bin/env python3
"""End-to-end indoor CFD prototype pipeline.

Stages:
1. Gemini generates indoor_cfd_scene_v1 JSON
2. Validator checks scene JSON
3. Gmsh generator creates .geo_unrolled / .msh
4. OpenFOAM case template is created
5. Optionally run gmshToFoam and checkMesh if available
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
TOOLS_PYTHON = Path("/home/yhjoo/.venvs/openclaw-tools/bin/python")


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(str(x) for x in cmd))
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if result.stdout.strip():
        print(result.stdout.rstrip())
    if result.stderr.strip():
        print(result.stderr.rstrip())
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with code {result.returncode}: {' '.join(str(x) for x in cmd)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the indoor Gemini->Gmsh->OpenFOAM prototype pipeline")
    parser.add_argument("--scenario", required=True, help="Natural-language indoor scene instruction")
    parser.add_argument("--name", default="indoor_pipeline_run", help="Run name for generated files/case folder")
    parser.add_argument("--backend", choices=["cli", "api"], default="api", help="Gemini backend for scene generation")
    parser.add_argument("--model", default="gemini-3.1-pro-preview", help="Gemini model for scene generation")
    parser.add_argument("--image", dest="images", action="append", default=[], help="Optional image/photo/rendering path; repeatable")
    parser.add_argument("--no-fallback", action="store_true", help="Disable fallback to alternate Gemini models")
    parser.add_argument("--mesh-size", type=float, default=0.35, help="Global gmsh target size in meters")
    parser.add_argument("--inlet-velocity", type=float, default=0.5, help="Inlet velocity for OpenFOAM template")
    parser.add_argument("--end-time", type=int, default=1000, help="simpleFoam controlDict endTime")
    parser.add_argument("--convert", action="store_true", help="Run gmshToFoam if available")
    parser.add_argument("--check-mesh", action="store_true", help="Run checkMesh after gmshToFoam")
    args = parser.parse_args()

    generated_dir = PROJECT_ROOT / "generated"
    case_root = PROJECT_ROOT / "cases"
    generated_dir.mkdir(parents=True, exist_ok=True)
    case_root.mkdir(parents=True, exist_ok=True)

    scene_json = generated_dir / f"{args.name}.json"
    geo_path = generated_dir / f"{args.name}.geo_unrolled"
    msh_path = generated_dir / f"{args.name}.msh"
    case_dir = case_root / args.name

    gen_cmd = [
        "python3",
        str(SCRIPTS_DIR / "generate_indoor_scene_with_gemini.py"),
        "--scenario",
        args.scenario,
        "--model",
        args.model,
        "-o",
        str(scene_json),
    ]
    if args.no_fallback:
        gen_cmd.append("--no-fallback")
    run(gen_cmd, cwd=PROJECT_ROOT)

    gmsh_python = str(TOOLS_PYTHON if TOOLS_PYTHON.exists() else shutil.which("python3") or "python3")
    gmsh_cmd = [
        gmsh_python,
        str(SCRIPTS_DIR / "scene_to_gmsh.py"),
        str(scene_json),
        "-o",
        str(msh_path),
        "--geo",
        str(geo_path),
        "--mesh-size",
        str(args.mesh_size),
    ]
    run(gmsh_cmd, cwd=PROJECT_ROOT)

    case_cmd = [
        "python3",
        str(SCRIPTS_DIR / "create_indoor_openfoam_case.py"),
        str(msh_path),
        str(case_dir),
        "--inlet-velocity",
        str(args.inlet_velocity),
        "--end-time",
        str(args.end_time),
    ]
    if args.convert:
        case_cmd.append("--convert")
    if args.check_mesh:
        case_cmd.append("--check-mesh")
    run(case_cmd, cwd=PROJECT_ROOT)

    print("\nPipeline complete.")
    print(f"  scene_json : {scene_json}")
    print(f"  geo        : {geo_path}")
    print(f"  msh        : {msh_path}")
    print(f"  case_dir   : {case_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
