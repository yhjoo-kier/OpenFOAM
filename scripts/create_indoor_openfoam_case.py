#!/usr/bin/env python3
"""Create a minimal OpenFOAM simpleFoam case from indoor scene mesh.

Pipeline stage:
    indoor scene JSON -> Gmsh .msh -> OpenFOAM case template

This script writes a practical first-pass indoor ventilation case using:
- incompressible steady-state solver: simpleFoam
- RANS turbulence: kOmegaSST
- expected patches: inlet, outlet, roomWalls, obstacleWalls

If OpenFOAM tools are available in PATH, it can also run gmshToFoam and checkMesh.
This script is intended to be used directly or from the end-to-end indoor pipeline runner.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def ensure_dirs(case_dir: Path) -> None:
    for rel in ["0", "constant", "system"]:
        (case_dir / rel).mkdir(parents=True, exist_ok=True)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n", encoding="utf-8")


def foam_header(obj: str, cls: str = "dictionary") -> str:
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}"""


def write_initial_conditions(case_dir: Path, inlet_velocity: float, turbulence_intensity: float, length_scale: float) -> None:
    k = 1.5 * (inlet_velocity * turbulence_intensity) ** 2
    omega = (k ** 0.5) / max(0.09 ** 0.25 * length_scale, 1e-9)

    write_text(
        case_dir / "0" / "U",
        f"""
{foam_header('U', 'volVectorField')}

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({inlet_velocity} 0 0);
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({inlet_velocity} 0 0);
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform ({inlet_velocity} 0 0);
        value           uniform ({inlet_velocity} 0 0);
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
        type            empty;
    }}
}}
""",
    )

    write_text(
        case_dir / "0" / "p",
        f"""
{foam_header('p', 'volScalarField')}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    inlet
    {{
        type            zeroGradient;
    }}

    outlet
    {{
        type            fixedValue;
        value           uniform 0;
    }}

    roomWalls
    {{
        type            zeroGradient;
    }}

    obstacleWalls
    {{
        type            zeroGradient;
    }}

    defaultFaces
    {{
        type            empty;
    }}
}}
""",
    )

    write_text(
        case_dir / "0" / "k",
        f"""
{foam_header('k', 'volScalarField')}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k:.8g};
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {k:.8g};
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {k:.8g};
        value           uniform {k:.8g};
    }}

    roomWalls
    {{
        type            kqRWallFunction;
        value           uniform {k:.8g};
    }}

    obstacleWalls
    {{
        type            kqRWallFunction;
        value           uniform {k:.8g};
    }}

    defaultFaces
    {{
        type            empty;
    }}
}}
""",
    )

    write_text(
        case_dir / "0" / "omega",
        f"""
{foam_header('omega', 'volScalarField')}

dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega:.8g};
boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {omega:.8g};
    }}

    outlet
    {{
        type            inletOutlet;
        inletValue      uniform {omega:.8g};
        value           uniform {omega:.8g};
    }}

    roomWalls
    {{
        type            omegaWallFunction;
        value           uniform {omega:.8g};
    }}

    obstacleWalls
    {{
        type            omegaWallFunction;
        value           uniform {omega:.8g};
    }}

    defaultFaces
    {{
        type            empty;
    }}
}}
""",
    )

    write_text(
        case_dir / "0" / "nut",
        f"""
{foam_header('nut', 'volScalarField')}

dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    inlet
    {{
        type            calculated;
        value           uniform 0;
    }}

    outlet
    {{
        type            calculated;
        value           uniform 0;
    }}

    roomWalls
    {{
        type            nutkWallFunction;
        value           uniform 0;
    }}

    obstacleWalls
    {{
        type            nutkWallFunction;
        value           uniform 0;
    }}

    defaultFaces
    {{
        type            empty;
    }}
}}
""",
    )


def write_constant(case_dir: Path) -> None:
    write_text(
        case_dir / "constant" / "transportProperties",
        f"""
{foam_header('transportProperties')}

transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] 1.5e-05;
""",
    )

    write_text(
        case_dir / "constant" / "turbulenceProperties",
        f"""
{foam_header('turbulenceProperties')}

simulationType  RAS;

RAS
{{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}}
""",
    )


def write_system(case_dir: Path, end_time: int) -> None:
    write_text(
        case_dir / "system" / "controlDict",
        f"""
{foam_header('controlDict')}

application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time};
deltaT          1;
writeControl    timeStep;
writeInterval   200;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""",
    )

    write_text(
        case_dir / "system" / "fvSchemes",
        f"""
{foam_header('fvSchemes')}

ddtSchemes
{{
    default         steadyState;
}}

gradSchemes
{{
    default         Gauss linear;
}}

divSchemes
{{
    default                         none;
    div(phi,U)                      bounded Gauss linearUpwind grad(U);
    div(phi,k)                      bounded Gauss upwind;
    div(phi,omega)                  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U)))))   Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         corrected;
}}

wallDist
{{
    method meshWave;
}}
""",
    )

    write_text(
        case_dir / "system" / "fvSolution",
        f"""
{foam_header('fvSolution')}

solvers
{{
    p
    {{
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.05;
        smoother        GaussSeidel;
    }}

    U
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }}

    "(k|omega)"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors 1;
    consistent      yes;
    residualControl
    {{
        p               1e-4;
        U               1e-5;
        k               1e-5;
        omega           1e-5;
    }}
}}

relaxationFactors
{{
    fields
    {{
        p               0.3;
    }}
    equations
    {{
        U               0.7;
        k               0.7;
        omega           0.7;
    }}
}}
""",
    )


def write_notes(case_dir: Path, msh_path: Path) -> None:
    write_text(
        case_dir / "README_case.md",
        f"""
# Indoor ventilation OpenFOAM case

Generated from:
- Mesh: `{msh_path}`
- Solver: `simpleFoam`
- Turbulence model: `kOmegaSST`

Expected mesh boundary patches:
- `inlet`
- `outlet`
- `roomWalls`
- `obstacleWalls`

Typical run sequence when OpenFOAM is available:
```bash
source /usr/share/openfoam/etc/bashrc  # or project env
cd {case_dir}
gmshToFoam {msh_path.name}
checkMesh
simpleFoam
```
""",
    )


def command_exists(name: str) -> bool:
    return shutil.which(name) is not None


def patch_boundary_types(case_dir: Path) -> None:
    boundary_path = case_dir / "constant" / "polyMesh" / "boundary"
    if not boundary_path.exists():
        return

    text = boundary_path.read_text(encoding="utf-8")
    replacements = {
        """    obstacleWalls
    {
        type            patch;
        physicalType    patch;""": """    obstacleWalls
    {
        type            wall;
        physicalType    wall;""",
        """    roomWalls
    {
        type            patch;
        physicalType    patch;""": """    roomWalls
    {
        type            wall;
        physicalType    wall;""",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    boundary_path.write_text(text, encoding="utf-8")


def run_command(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())


def maybe_convert_mesh(case_dir: Path, msh_path: Path, run_check_mesh: bool) -> None:
    if not command_exists("gmshToFoam"):
        print("gmshToFoam not found in PATH; wrote case template only.")
        return

    local_msh = case_dir / msh_path.name
    if msh_path.resolve() != local_msh.resolve():
        shutil.copy2(msh_path, local_msh)

    run_command(["gmshToFoam", local_msh.name], cwd=case_dir)
    patch_boundary_types(case_dir)

    if run_check_mesh:
        if command_exists("checkMesh"):
            run_command(["checkMesh"], cwd=case_dir)
        else:
            print("checkMesh not found in PATH; skipping mesh check.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a simpleFoam case from an indoor Gmsh mesh")
    parser.add_argument("msh_file", type=Path, help="Path to Gmsh .msh file")
    parser.add_argument("case_dir", type=Path, help="Target OpenFOAM case directory")
    parser.add_argument("--inlet-velocity", type=float, default=0.5, help="Inlet velocity magnitude in m/s")
    parser.add_argument("--turbulence-intensity", type=float, default=0.05, help="Inlet turbulence intensity (0-1)")
    parser.add_argument("--length-scale", type=float, default=0.1, help="Reference turbulence length scale in m")
    parser.add_argument("--end-time", type=int, default=1000, help="Pseudo-time / max iterations for simpleFoam output control")
    parser.add_argument("--convert", action="store_true", help="Run gmshToFoam if available")
    parser.add_argument("--check-mesh", action="store_true", help="Run checkMesh after conversion if available")
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    msh_path = args.msh_file.resolve()
    ensure_dirs(case_dir)
    write_initial_conditions(case_dir, args.inlet_velocity, args.turbulence_intensity, args.length_scale)
    write_constant(case_dir)
    write_system(case_dir, args.end_time)
    write_notes(case_dir, msh_path)

    if args.convert:
        maybe_convert_mesh(case_dir, msh_path, args.check_mesh)
    else:
        local_msh = case_dir / msh_path.name
        if msh_path.resolve() != local_msh.resolve():
            shutil.copy2(msh_path, local_msh)

    print(f"Case prepared at: {case_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
