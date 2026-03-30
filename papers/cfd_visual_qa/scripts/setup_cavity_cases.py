#!/usr/bin/env python3
"""Setup lid-driven cavity (S5) OpenFOAM cases for CFD Visual QA pilot."""
from __future__ import annotations

import json
import shutil
import subprocess
import textwrap
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]  # OpenFOAM/
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"

# Geometry: square cavity 1m × 1m × 0.1m
L = 1.0   # cavity side length
W = 0.1   # z-width (2D, 1 cell)

CASES = {
    # --- Correct cases ---
    "S5_correct_lam": {
        "Re": 100, "U": 1.0, "nu": 0.01, "turb": False,
        "nx": 60, "ny": 60, "grading_y": 1.0,
        "endTime": 5000,
        "lid": "top",   # movingWall = top (y=1)
        "description": "Laminar lid-driven cavity, Re=100, 60x60 mesh",
    },
    "S5_correct_turb": {
        "Re": 10000, "U": 1.0, "nu": 1e-4, "turb": True,
        "nx": 100, "ny": 100, "grading_y": 4.0,
        "endTime": 10000,
        "lid": "top",
        "description": "Turbulent lid-driven cavity, Re=10000, k-omega SST, 100x100 graded mesh",
    },
    # --- Error cases ---
    "S5_E1_underconverged": {
        "Re": 10000, "U": 1.0, "nu": 1e-4, "turb": True,
        "nx": 100, "ny": 100, "grading_y": 4.0,
        "endTime": 30,   # early stop!
        "lid": "top",
        "description": "Under-converged: stopped at 30 iterations",
    },
    "S5_E2_bc_swap": {
        "Re": 10000, "U": 1.0, "nu": 1e-4, "turb": True,
        "nx": 100, "ny": 100, "grading_y": 4.0,
        "endTime": 10000,
        "lid": "bottom",   # ERROR: moving wall on bottom instead of top
        "description": "BC error: lid velocity applied to BOTTOM wall instead of top",
    },
    "S5_E5_coarse_mesh": {
        "Re": 10000, "U": 1.0, "nu": 1e-4, "turb": True,
        "nx": 10, "ny": 10, "grading_y": 1.0,   # extremely coarse, no grading
        "endTime": 10000,
        "lid": "top",
        "description": "Extremely coarse mesh: 10x10, no wall grading",
    },
    "S5_E8_reversed_lid": {
        "Re": 100, "U": -1.0, "nu": 0.01, "turb": False,  # negative U → mirror image
        "nx": 60, "ny": 60, "grading_y": 1.0,
        "endTime": 5000,
        "lid": "top",
        "description": "Reversed lid direction: U=-1 (mirror image of correct laminar)",
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
    nx = cfg["nx"]
    ny = cfg["ny"]
    gy = cfg["grading_y"]
    # grading: refine near walls in y (bottom→top: inv_gy, top→bottom applied symmetrically
    # For cavity: use edgeGrading to refine near all walls
    # Simple approach: uniform grading in y (refine toward top wall)
    inv_gy = round(1.0 / gy, 6) if gy != 1.0 else 1.0

    # 8 vertices: (x, y, z)
    # Vertex ordering for hex block (OpenFOAM right-hand rule)
    # Front face (z=0): 0-1-2-3 counterclockwise from front view
    # Back face (z=W): 4-5-6-7
    #   3-------2
    #   |       |   y
    #   |       |   ^
    #   0-------1   --> x
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
        (0  0  0)    // 0  bottom-left  front
        ({L} 0  0)    // 1  bottom-right front
        ({L} {L} 0)   // 2  top-right    front
        (0  {L} 0)    // 3  top-left     front
        (0  0  {W})   // 4  bottom-left  back
        ({L} 0  {W})  // 5  bottom-right back
        ({L} {L} {W}) // 6  top-right    back
        (0  {L} {W})  // 7  top-left     back
    );

    blocks
    (
        // hex: vertices in OpenFOAM order (0 1 2 3 4 5 6 7)
        // x-grading=1 (uniform), y-grading refines toward top (grading_y={gy})
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 {gy} 1)
    );

    edges ( );

    boundary
    (
        movingWall
        {{
            type wall;
            faces
            (
                (3 7 6 2)   // top face (y={L})
            );
        }}
        fixedWalls
        {{
            type wall;
            faces
            (
                (0 1 5 4)   // bottom face (y=0)
                (0 4 7 3)   // left face   (x=0)
                (1 2 6 5)   // right face  (x={L})
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
    write_interval = end_time  # write only final state
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
    writeInterval   {write_interval};
    purgeWrite      1;
    writeFormat     ascii;
    writePrecision  8;
    writeCompression off;
    timeFormat      general;
    timePrecision   6;
    runTimeModifiable true;
    """)
    (case_dir / "system" / "controlDict").write_text(content)


def write_fvschemes(case_dir: Path, cfg: dict) -> None:
    turb_schemes = ""
    if cfg["turb"]:
        turb_schemes = textwrap.dedent("""\
            div(phi,k)      bounded Gauss linearUpwind grad(k);
            div(phi,omega)  bounded Gauss linearUpwind grad(omega);
        """)
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSchemes;
    }}

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
        default         none;
        div(phi,U)      bounded Gauss linearUpwind grad(U);
        {turb_schemes}    div(phi,nuTilda) bounded Gauss linearUpwind grad(nuTilda);
        div((nuEff*dev2(T(grad(U))))) Gauss linear;
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
    """)
    (case_dir / "system" / "fvSchemes").write_text(content)


def write_fvsolution(case_dir: Path, cfg: dict) -> None:
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
        "(k|omega|nuTilda)"
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
        pRefCell        0;
        pRefValue       0;
        residualControl
        {
            p               1e-5;
            U               1e-5;
            "(k|omega)"     1e-5;
        }
    }

    relaxationFactors
    {
        equations
        {
            U               0.7;
            p               0.3;
            k               0.7;
            omega           0.7;
        }
    }
    """)
    (case_dir / "system" / "fvSolution").write_text(content)


def write_initial_conditions(case_dir: Path, cfg: dict) -> None:
    U_lid = cfg["U"]   # may be negative for reversed case
    turb = cfg["turb"]
    lid = cfg["lid"]   # "top" or "bottom"

    # Decide which patch gets the moving wall BC
    if lid == "top":
        moving_patch = "movingWall"
        fixed_patch = "fixedWalls"
    else:
        # E2: put the velocity on fixedWalls (bottom), movingWall becomes noSlip
        moving_patch = "fixedWalls"
        fixed_patch = "movingWall"

    # U field
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
        {moving_patch}
        {{
            type            fixedValue;
            value           uniform ({U_lid} 0 0);
        }}
        {fixed_patch}
        {{
            type            noSlip;
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    # p field — for a closed cavity (no inlet/outlet), use a reference pressure
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
        movingWall
        {
            type            zeroGradient;
        }
        fixedWalls
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)

    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)
    (zero_dir / "U").write_text(u_content)
    (zero_dir / "p").write_text(p_content)

    if turb:
        # Turbulent quantities — wall-bounded cavity, use moderate initial values
        U_ref = abs(U_lid)
        k_val = 1.5 * (0.05 * U_ref) ** 2   # 5% turbulence intensity
        # omega based on cavity length scale
        omega_val = k_val ** 0.5 / (0.09 ** 0.25 * 0.1 * L)

        k_content = textwrap.dedent(f"""\
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            object      k;
        }}

        dimensions      [0 2 -2 0 0 0 0];
        internalField   uniform {k_val:.6e};

        boundaryField
        {{
            movingWall
            {{
                type            kqRWallFunction;
                value           uniform {k_val:.6e};
            }}
            fixedWalls
            {{
                type            kqRWallFunction;
                value           uniform {k_val:.6e};
            }}
            frontAndBack
            {{
                type            empty;
            }}
        }}
        """)

        omega_content = textwrap.dedent(f"""\
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            object      omega;
        }}

        dimensions      [0 0 -1 0 0 0 0];
        internalField   uniform {omega_val:.6e};

        boundaryField
        {{
            movingWall
            {{
                type            omegaWallFunction;
                value           uniform {omega_val:.6e};
            }}
            fixedWalls
            {{
                type            omegaWallFunction;
                value           uniform {omega_val:.6e};
            }}
            frontAndBack
            {{
                type            empty;
            }}
        }}
        """)

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
            movingWall
            {
                type            nutkWallFunction;
                value           uniform 0;
            }
            fixedWalls
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

        (zero_dir / "k").write_text(k_content)
        (zero_dir / "omega").write_text(omega_content)
        (zero_dir / "nut").write_text(nut_content)


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


def write_turbulence_properties(case_dir: Path, cfg: dict) -> None:
    if cfg["turb"]:
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
    write_fvschemes(case_dir, cfg)
    write_fvsolution(case_dir, cfg)
    write_initial_conditions(case_dir, cfg)
    write_transport_properties(case_dir, cfg)
    write_turbulence_properties(case_dir, cfg)

    # Determine error type
    error_type = None
    if "_E1_" in name:
        error_type = "E1_underconverged"
    elif "_E2_" in name:
        error_type = "E2_bc_swap"
    elif "_E5_" in name:
        error_type = "E5_coarse_mesh"
    elif "_E8_" in name:
        error_type = "E8_reversed_lid"

    meta = {
        "case_name": name,
        "scenario": "S5",
        "scenario_name": "lid_driven_cavity",
        "Re": cfg["Re"],
        "U_lid": cfg["U"],
        "nu": cfg["nu"],
        "turbulence": cfg["turb"],
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "grading_y": cfg["grading_y"],
        "lid_wall": cfg["lid"],
        "endTime": cfg["endTime"],
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name} (Re={cfg['Re']}, turb={cfg['turb']}, lid={cfg['lid']}, "
          f"mesh={cfg['nx']}x{cfg['ny']})")
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
    print("CFD Visual QA — S5: Lid-Driven Cavity")
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
