#!/usr/bin/env python3
"""Setup laminar channel flow (S2) OpenFOAM cases for CFD Visual QA benchmark."""
from __future__ import annotations

import json
import shutil
import subprocess
import textwrap
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]  # OpenFOAM/
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"

# Base geometry
H = 1.0   # channel height (m) — Re length scale
W = 0.1   # z-width (2D, 1 cell)

CASES = {
    # --- Correct cases ---
    "S2_correct_lam": {
        "Re": 100, "U": 1.0, "nu": 0.01,
        "L": 10.0, "nx": 200, "ny": 40,
        "endTime": 5000,
        "swap_bc": False,
        "description": "Laminar channel flow, Re=100, 200x40 mesh, developing parabolic profile",
    },
    "S2_correct_developed": {
        "Re": 100, "U": 1.0, "nu": 0.01,
        "L": 20.0, "nx": 400, "ny": 40,
        "endTime": 5000,
        "swap_bc": False,
        "description": "Laminar channel flow, Re=100, L=20m 400x40 mesh, fully developed at outlet",
    },
    # --- Error cases ---
    "S2_E1_underconverged": {
        "Re": 100, "U": 1.0, "nu": 0.01,
        "L": 10.0, "nx": 200, "ny": 40,
        "endTime": 20,   # early stop — not converged
        "swap_bc": False,
        "description": "Under-converged: stopped at 20 iterations (correct setup, wrong endTime)",
    },
    "S2_E2_bc_swap": {
        "Re": 100, "U": 1.0, "nu": 0.01,
        "L": 10.0, "nx": 200, "ny": 40,
        "endTime": 5000,
        "swap_bc": True,  # fixedValue U on right, zeroGradient p on right
        "description": "BC error: inlet/outlet swapped — fixedValue velocity on outlet face",
    },
    "S2_E5_coarse_mesh": {
        "Re": 100, "U": 1.0, "nu": 0.01,
        "L": 10.0, "nx": 20, "ny": 4,   # very coarse — 20x4
        "endTime": 5000,
        "swap_bc": False,
        "description": "Extremely coarse mesh: 20x4, poor boundary layer resolution",
    },
    "S2_E3_wrong_viscosity": {
        "Re": 100, "U": 1.0, "nu": 0.1,  # 10x too high — effective Re=10
        "L": 10.0, "nx": 200, "ny": 40,
        "endTime": 5000,
        "swap_bc": False,
        "description": "Wrong viscosity: nu=0.1 (10x too high), effective Re=10, over-diffused flow",
    },
}


def docker_exec(case_dir: Path, command: str, timeout: float = 900) -> subprocess.CompletedProcess:
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{PROJECT_ROOT}:/app",
        "-w", f"/app/{case_dir.relative_to(PROJECT_ROOT)}",
        DOCKER_IMAGE,
        "bash", "-lc", command,
    ]
    print("+", " ".join(cmd[:6]), "...", command[:80])
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=timeout)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[:500]}")
    return result


def write_blockmesh(case_dir: Path, cfg: dict) -> None:
    """Single hex block: x=[0, L], y=[0, H], z=[0, W].

    Vertex layout (front face z=0, back face z=W):
        3 --- 2          y
        |     |          ^
        0 --- 1          --> x
    Back: 7 --- 6 / 4 --- 5
    """
    L = cfg["L"]
    nx = cfg["nx"]
    ny = cfg["ny"]

    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }}

    scale 1;

    vertices
    (
        (0    0    0  )   // 0  bottom-left  front
        ({L}  0    0  )   // 1  bottom-right front
        ({L}  {H}  0  )   // 2  top-right    front
        (0    {H}  0  )   // 3  top-left     front
        (0    0    {W})   // 4  bottom-left  back
        ({L}  0    {W})   // 5  bottom-right back
        ({L}  {H}  {W})   // 6  top-right    back
        (0    {H}  {W})   // 7  top-left     back
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)
    );

    edges ( );

    boundary
    (
        inlet
        {{
            type patch;
            faces
            (
                (0 4 7 3)   // left face (x=0)
            );
        }}
        outlet
        {{
            type patch;
            faces
            (
                (1 2 6 5)   // right face (x={L})
            );
        }}
        topWall
        {{
            type wall;
            faces
            (
                (3 7 6 2)   // top face (y={H})
            );
        }}
        bottomWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)   // bottom face (y=0)
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                (0 3 2 1)   // front face (z=0)
                (4 5 6 7)   // back face  (z={W})
            );
        }}
    );

    mergePatchPairs ( );
    """)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "system" / "blockMeshDict").write_text(content)


def write_control_dict(case_dir: Path, cfg: dict) -> None:
    end_time = cfg["endTime"]
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }}

    application     simpleFoam;
    startFrom       startTime;
    startTime       0;
    stopAt          endTime;
    endTime         {end_time};
    deltaT          1;
    writeControl    timeStep;
    writeInterval   {end_time};
    purgeWrite      1;
    writeFormat     ascii;
    writePrecision  8;
    writeCompression off;
    timeFormat      general;
    timePrecision   6;
    runTimeModifiable true;
    """)
    (case_dir / "system" / "controlDict").write_text(content)


def write_fvschemes(case_dir: Path) -> None:
    content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSchemes;
    }

    ddtSchemes
    {
        default         steadyState;
    }

    gradSchemes
    {
        default         Gauss linear;
    }

    divSchemes
    {
        default         none;
        div(phi,U)      bounded Gauss linearUpwind grad(U);
        div((nuEff*dev2(T(grad(U))))) Gauss linear;
    }

    laplacianSchemes
    {
        default         Gauss linear corrected;
    }

    interpolationSchemes
    {
        default         linear;
    }

    snGradSchemes
    {
        default         corrected;
    }

    wallDist
    {
        method meshWave;
    }
    """)
    (case_dir / "system" / "fvSchemes").write_text(content)


def write_fvsolution(case_dir: Path) -> None:
    # Channel flow has inlet/outlet — no pRefCell needed
    content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSolution;
    }

    solvers
    {
        p
        {
            solver          GAMG;
            smoother        GaussSeidel;
            tolerance       1e-7;
            relTol          0.01;
        }
        U
        {
            solver          smoothSolver;
            smoother        GaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }
    }

    SIMPLE
    {
        nNonOrthogonalCorrectors 0;
        consistent      yes;
        residualControl
        {
            p               1e-5;
            U               1e-5;
        }
    }

    relaxationFactors
    {
        equations
        {
            U               0.7;
            p               0.3;
        }
    }
    """)
    (case_dir / "system" / "fvSolution").write_text(content)


def write_initial_conditions(case_dir: Path, cfg: dict) -> None:
    U_in = cfg["U"]
    swap = cfg.get("swap_bc", False)

    # swap_bc: fixedValue U goes on outlet face, zeroGradient p goes on outlet face
    if swap:
        vel_fixed_patch = "outlet"
        vel_zero_patch = "inlet"
        p_fixed_patch = "inlet"
        p_zero_patch = "outlet"
    else:
        vel_fixed_patch = "inlet"
        vel_zero_patch = "outlet"
        p_fixed_patch = "outlet"
        p_zero_patch = "inlet"

    u_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }}

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform (0 0 0);

    boundaryField
    {{
        {vel_fixed_patch}
        {{
            type            fixedValue;
            value           uniform ({U_in} 0 0);
        }}
        {vel_zero_patch}
        {{
            type            zeroGradient;
        }}
        topWall
        {{
            type            noSlip;
        }}
        bottomWall
        {{
            type            noSlip;
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    p_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }}

    dimensions      [0 2 -2 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {{
        {p_zero_patch}
        {{
            type            zeroGradient;
        }}
        {p_fixed_patch}
        {{
            type            fixedValue;
            value           uniform 0;
        }}
        topWall
        {{
            type            zeroGradient;
        }}
        bottomWall
        {{
            type            zeroGradient;
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    zero = case_dir / "0"
    zero.mkdir(parents=True, exist_ok=True)
    (zero / "U").write_text(u_content)
    (zero / "p").write_text(p_content)


def write_transport_properties(case_dir: Path, cfg: dict) -> None:
    nu = cfg["nu"]
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      transportProperties;
    }}

    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] {nu:.6e};
    """)
    (case_dir / "constant").mkdir(parents=True, exist_ok=True)
    (case_dir / "constant" / "transportProperties").write_text(content)


def write_turbulence_properties(case_dir: Path) -> None:
    # All S2 cases are laminar
    content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      turbulenceProperties;
    }

    simulationType  laminar;
    """)
    (case_dir / "constant" / "turbulenceProperties").write_text(content)


def setup_case(name: str, cfg: dict) -> Path:
    case_dir = CASE_BASE / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)

    write_blockmesh(case_dir, cfg)
    write_control_dict(case_dir, cfg)
    write_fvschemes(case_dir)
    write_fvsolution(case_dir)
    write_initial_conditions(case_dir, cfg)
    write_transport_properties(case_dir, cfg)
    write_turbulence_properties(case_dir)

    # Determine error type
    error_type = None
    if "_E1_" in name:
        error_type = "E1_underconverged"
    elif "_E2_" in name:
        error_type = "E2_bc_swap"
    elif "_E3_" in name:
        error_type = "E3_wrong_viscosity"
    elif "_E5_" in name:
        error_type = "E5_coarse_mesh"

    meta = {
        "case_name": name,
        "scenario": "S2",
        "scenario_name": "laminar_channel_flow",
        "Re": cfg["Re"],
        "U_inlet": cfg["U"],
        "nu": cfg["nu"],
        "channel_length": cfg["L"],
        "channel_height": H,
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "endTime": cfg["endTime"],
        "swap_bc": cfg.get("swap_bc", False),
        "turbulence": False,
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name} (Re={cfg['Re']}, nu={cfg['nu']}, L={cfg['L']}, "
          f"mesh={cfg['nx']}x{cfg['ny']}, endTime={cfg['endTime']})")
    return case_dir


def run_case(name: str, case_dir: Path) -> bool:
    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"{'='*60}")

    # blockMesh
    r = docker_exec(case_dir, "blockMesh > log.blockMesh 2>&1", timeout=120)
    if r.returncode != 0:
        print(f"  blockMesh FAILED for {name}")
        return False

    # checkMesh
    docker_exec(case_dir, "checkMesh > log.checkMesh 2>&1", timeout=120)

    # simpleFoam — allow non-zero return for error cases (underconverged, bc_swap)
    r = docker_exec(case_dir, "simpleFoam > log.simpleFoam 2>&1", timeout=900)
    if r.returncode != 0:
        print(f"  simpleFoam returned {r.returncode} for {name} (may be OK for error cases)")

    # foamToVTK — export latest time step
    r = docker_exec(case_dir, "foamToVTK -latestTime > log.foamToVTK 2>&1", timeout=600)
    if r.returncode != 0:
        print(f"  foamToVTK FAILED for {name}")
        return False

    print(f"  Completed: {name}")
    return True


def main():
    CASE_BASE.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("CFD Visual QA — S2: Laminar Channel Flow")
    print("=" * 60)

    results = {}
    for name, cfg in CASES.items():
        case_dir = setup_case(name, cfg)
        success = run_case(name, case_dir)
        results[name] = success

    print("\n" + "=" * 60)
    print("Summary:")
    for name, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  {name}: {status}")
    print("=" * 60)


if __name__ == "__main__":
    main()
