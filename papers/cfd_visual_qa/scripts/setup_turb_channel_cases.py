#!/usr/bin/env python3
"""Setup turbulent channel flow (S3) OpenFOAM cases for CFD Visual QA benchmark.

Re=10000, U=1 m/s, nu=1e-4, k-omega SST, simpleFoam.
Wall-graded blockMesh: symmetric refinement near top and bottom walls.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import textwrap
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]  # OpenFOAM/
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"

# Base geometry
L = 10.0   # channel length (m)
H = 1.0    # channel height (m) — Re length scale
W = 0.1    # z-width (2D, 1 cell)

# Turbulence inlet values for Re=10000, U=1 m/s, H=1m
# Turbulence intensity 5%, length scale 0.07*H
TI = 0.05
C_MU = 0.09
L_MIX = 0.07 * H
K_INLET = 1.5 * (TI * 1.0) ** 2          # = 3.75e-4
OMEGA_INLET = math.sqrt(K_INLET) / (C_MU ** 0.25 * L_MIX)  # ≈ 2.36

CASES = {
    # --- Correct cases ---
    "S3_correct_turb": {
        "Re": 10000, "U": 1.0, "nu": 1e-4,
        "nx": 200, "ny": 60,
        "grading": "symmetric",   # finer near both walls
        "endTime": 10000,
        "turbulence": "kOmegaSST",
        "description": "Turbulent channel Re=10000, k-omega SST, 200x60 wall-graded mesh",
    },
    "S3_correct_fine": {
        "Re": 10000, "U": 1.0, "nu": 1e-4,
        "nx": 300, "ny": 100,
        "grading": "symmetric_strong",   # stronger grading
        "endTime": 10000,
        "turbulence": "kOmegaSST",
        "description": "Turbulent channel Re=10000, k-omega SST, 300x100 finer wall-graded mesh",
    },
    # --- Error cases ---
    "S3_E1_underconverged": {
        "Re": 10000, "U": 1.0, "nu": 1e-4,
        "nx": 200, "ny": 60,
        "grading": "symmetric",
        "endTime": 30,   # early stop — not converged
        "turbulence": "kOmegaSST",
        "description": "Under-converged: stopped at 30 iterations (correct setup, wrong endTime)",
    },
    "S3_E4_wrong_turb_model": {
        "Re": 10000, "U": 1.0, "nu": 1e-4,
        "nx": 200, "ny": 60,
        "grading": "symmetric",
        "endTime": 10000,
        "turbulence": "laminar",   # LAMINAR at turbulent Re — wrong physics
        "description": "Wrong model: turbulent Re=10000 solved with laminar assumption — no turbulent diffusion",
    },
    "S3_E5_coarse_mesh": {
        "Re": 10000, "U": 1.0, "nu": 1e-4,
        "nx": 20, "ny": 6,
        "grading": "none",   # no wall grading — y+ completely unresolved
        "endTime": 10000,
        "turbulence": "kOmegaSST",
        "description": "Extremely coarse mesh: 20x6 no wall grading, y+ unresolved for k-omega SST",
    },
    "S3_E3_wrong_viscosity": {
        "Re": 10000, "U": 1.0, "nu": 0.01,   # 100x too high, effective Re=100
        "nx": 200, "ny": 60,
        "grading": "symmetric",
        "endTime": 10000,
        "turbulence": "kOmegaSST",
        "description": "Wrong viscosity: nu=0.01 (100x too high), effective Re=100 but turbulence model active",
    },
}

# y-grading spec: simpleGrading value for the y-direction
# "symmetric": finer near both walls — use multi-grading list
# "symmetric_strong": even finer near walls
# "none": uniform
GRADING_SPEC = {
    "none": "1",
    "symmetric": "((0.5 0.5 4)(0.5 0.5 0.25))",
    "symmetric_strong": "((0.5 0.5 8)(0.5 0.5 0.125))",
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
    nx = cfg["nx"]
    ny = cfg["ny"]
    y_grading = GRADING_SPEC[cfg["grading"]]

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
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 {y_grading} 1)
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
        grad(U)         cellLimited Gauss linear 1;
        grad(k)         cellLimited Gauss linear 1;
        grad(omega)     cellLimited Gauss linear 1;
    }

    divSchemes
    {
        default         none;
        div(phi,U)      bounded Gauss linearUpwind grad(U);
        div(phi,k)      bounded Gauss linearUpwind grad(k);
        div(phi,omega)  bounded Gauss linearUpwind grad(omega);
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


def write_fvsolution(case_dir: Path, turbulence: str) -> None:
    """fvSolution with turbulence solvers when using k-omega SST."""
    turb_solvers = ""
    turb_relax = ""
    turb_residual = ""
    if turbulence == "kOmegaSST":
        turb_solvers = textwrap.dedent("""\
        k
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }
        omega
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }
        """)
        turb_relax = textwrap.dedent("""\
            k               0.5;
            omega           0.5;
        """)
        turb_residual = textwrap.dedent("""\
            k               1e-5;
            omega           1e-5;
        """)

    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSolution;
    }}

    solvers
    {{
        p
        {{
            solver          GAMG;
            smoother        GaussSeidel;
            tolerance       1e-7;
            relTol          0.01;
        }}
        U
        {{
            solver          smoothSolver;
            smoother        GaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }}
        {turb_solvers}
    }}

    SIMPLE
    {{
        nNonOrthogonalCorrectors 0;
        consistent      yes;
        residualControl
        {{
            p               1e-5;
            U               1e-5;
            {turb_residual}
        }}
    }}

    relaxationFactors
    {{
        equations
        {{
            U               0.7;
            p               0.3;
            {turb_relax}
        }}
    }}
    """)
    (case_dir / "system" / "fvSolution").write_text(content)


def write_initial_conditions(case_dir: Path, cfg: dict) -> None:
    """Write 0/U, 0/p, and (if turbulent) 0/k, 0/omega, 0/nut."""
    U_in = cfg["U"]
    turbulence = cfg["turbulence"]

    zero = case_dir / "0"
    zero.mkdir(parents=True, exist_ok=True)

    # --- U ---
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
        inlet
        {{
            type            fixedValue;
            value           uniform ({U_in} 0 0);
        }}
        outlet
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
    (zero / "U").write_text(u_content)

    # --- p ---
    p_content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }

    dimensions      [0 2 -2 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {
        inlet
        {
            type            zeroGradient;
        }
        outlet
        {
            type            fixedValue;
            value           uniform 0;
        }
        topWall
        {
            type            zeroGradient;
        }
        bottomWall
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)
    (zero / "p").write_text(p_content)

    # --- Turbulent fields (only for kOmegaSST) ---
    if turbulence != "kOmegaSST":
        return

    k_val = f"{K_INLET:.6e}"
    omega_val = f"{OMEGA_INLET:.6f}"

    k_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      k;
    }}

    dimensions      [0 2 -2 0 0 0 0];
    internalField   uniform {k_val};

    boundaryField
    {{
        inlet
        {{
            type            fixedValue;
            value           uniform {k_val};
        }}
        outlet
        {{
            type            zeroGradient;
        }}
        topWall
        {{
            type            kqRWallFunction;
            value           uniform {k_val};
        }}
        bottomWall
        {{
            type            kqRWallFunction;
            value           uniform {k_val};
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)
    (zero / "k").write_text(k_content)

    omega_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      omega;
    }}

    dimensions      [0 0 -1 0 0 0 0];
    internalField   uniform {omega_val};

    boundaryField
    {{
        inlet
        {{
            type            fixedValue;
            value           uniform {omega_val};
        }}
        outlet
        {{
            type            zeroGradient;
        }}
        topWall
        {{
            type            omegaWallFunction;
            value           uniform {omega_val};
        }}
        bottomWall
        {{
            type            omegaWallFunction;
            value           uniform {omega_val};
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)
    (zero / "omega").write_text(omega_content)

    nut_content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      nut;
    }

    dimensions      [0 2 -1 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {
        inlet
        {
            type            calculated;
            value           uniform 0;
        }
        outlet
        {
            type            calculated;
            value           uniform 0;
        }
        topWall
        {
            type            nutkWallFunction;
            value           uniform 0;
        }
        bottomWall
        {
            type            nutkWallFunction;
            value           uniform 0;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)
    (zero / "nut").write_text(nut_content)


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


def write_turbulence_properties(case_dir: Path, turbulence: str) -> None:
    if turbulence == "kOmegaSST":
        content = textwrap.dedent("""\
        FoamFile
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
        """)
    else:
        # laminar — E4 wrong turbulence model case
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
    write_fvsolution(case_dir, cfg["turbulence"])
    write_initial_conditions(case_dir, cfg)
    write_transport_properties(case_dir, cfg)
    write_turbulence_properties(case_dir, cfg["turbulence"])

    # Determine error type
    error_type = None
    if "_E1_" in name:
        error_type = "E1_underconverged"
    elif "_E3_" in name:
        error_type = "E3_wrong_viscosity"
    elif "_E4_" in name:
        error_type = "E4_wrong_turb_model"
    elif "_E5_" in name:
        error_type = "E5_coarse_mesh"

    meta = {
        "case_name": name,
        "scenario": "S3",
        "scenario_name": "turbulent_channel_flow",
        "Re": cfg["Re"],
        "U_inlet": cfg["U"],
        "nu": cfg["nu"],
        "channel_length": L,
        "channel_height": H,
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "mesh_grading": cfg["grading"],
        "endTime": cfg["endTime"],
        "turbulence_model": cfg["turbulence"],
        "k_inlet": K_INLET,
        "omega_inlet": round(OMEGA_INLET, 4),
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name} (Re={cfg['Re']}, nu={cfg['nu']}, "
          f"mesh={cfg['nx']}x{cfg['ny']}, grading={cfg['grading']}, "
          f"turb={cfg['turbulence']}, endTime={cfg['endTime']})")
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

    # simpleFoam — allow non-zero return for error cases (underconverged etc.)
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
    print("CFD Visual QA — S3: Turbulent Channel Flow")
    print(f"  Re=10000, U=1 m/s, k-omega SST")
    print(f"  k_inlet={K_INLET:.3e}, omega_inlet={OMEGA_INLET:.2f}")
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
