#!/usr/bin/env python3
"""Setup heated vertical plate / natural convection boundary layer (S1) OpenFOAM cases
for CFD Visual QA benchmark.

Physics: buoyantBoussinesqSimpleFoam — incompressible with Boussinesq buoyancy.
Geometry: 2D, 2m wide × 4m tall, 0.1m deep (1 cell in z, empty BC).
  - Left wall (x=0): heated plate, T=320K, from y=0 to y=4m
  - Right boundary (x=2): open (pressure outlet — fixedValue p_rgh=0)
  - Top (y=4): open (pressure outlet)
  - Bottom (y=0): wall, adiabatic (noSlip, zeroGradient T)

Rayleigh number (plate height H=4m, deltaT=20K, nu=1e-3):
  Ra = g * beta * dT * H^3 * Pr / nu^2
     = 9.81 * 3.4e-3 * 20 * 64 * 0.71 / 1e-6 ≈ 3.02e6
"""
from __future__ import annotations

import json
import shutil
import subprocess
import textwrap
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]  # OpenFOAM/
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"

# Geometry
W_DOMAIN = 2.0   # domain width  (x-direction, m)
H_DOMAIN = 4.0   # domain height (y-direction, m)
W_DEPTH  = 0.1   # z-width (2D, 1 cell)

# Physical parameters
G     = 9.81
BETA  = 3.4e-3   # thermal expansion coefficient [1/K]
DT    = 20.0     # T_hot - T_ref  [K]
T_HOT = 320.0    # heated plate temperature [K]
T_REF = 300.0    # ambient / reference temperature [K]
PR    = 0.71     # Prandtl number

CASES: dict[str, dict] = {
    # -------------------------------------------------------------------------
    # Correct cases
    # -------------------------------------------------------------------------
    "S1_correct_lam": {
        "nu": 1e-3,
        "nx": 60, "ny": 120,
        "grading_x": 0.2,   # simpleGrading: finer near x=0 (heated wall)
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: laminar natural convection boundary layer, "
            "Ra≈3e6, 60×120 wall-graded mesh, nu=1e-3"
        ),
    },
    "S1_correct_fine": {
        "nu": 1e-3,
        "nx": 100, "ny": 200,
        "grading_x": 0.2,
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: laminar natural convection boundary layer, "
            "Ra≈3e6, 100×200 fine wall-graded mesh, nu=1e-3"
        ),
    },
    # -------------------------------------------------------------------------
    # Error cases
    # -------------------------------------------------------------------------
    "S1_E1_underconverged": {
        "nu": 1e-3,
        "nx": 60, "ny": 120,
        "grading_x": 0.2,
        "endTime": 30,        # ERROR: only 30 iterations — far from steady state
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E1: under-converged — stopped at 30 iterations, "
            "boundary layer not developed"
        ),
    },
    "S1_E8_gravity_flipped": {
        "nu": 1e-3,
        "nx": 60, "ny": 120,
        "grading_x": 0.2,
        "endTime": 10000,
        "gravity": (0.0, +9.81, 0.0),   # ERROR: buoyancy reversed — flow descends
        "description": (
            "Error E8: gravity flipped — g=(0,+9.81,0), "
            "buoyancy-driven flow descends instead of rising"
        ),
    },
    "S1_E5_coarse_mesh": {
        "nu": 1e-3,
        "nx": 10, "ny": 20,   # ERROR: 10×20 far too coarse to resolve BL
        "grading_x": 1.0,     # no grading (coarseness is the error)
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E5: coarse mesh — only 10×20 uniform, "
            "boundary layer entirely unresolved"
        ),
    },
    "S1_E3_wrong_viscosity": {
        "nu": 0.1,             # ERROR: 100× too high — suppresses convective BL
        "nx": 60, "ny": 120,
        "grading_x": 0.2,
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E3: wrong viscosity — nu=0.1 (100× too high), "
            "Ra drops to ~3e1, kills the buoyancy-driven boundary layer"
        ),
    },
}


# ---------------------------------------------------------------------------
# Docker helper
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_blockmesh(case_dir: Path, cfg: dict) -> None:
    nx = cfg["nx"]
    ny = cfg["ny"]
    gx = cfg["grading_x"]  # <1 means finer near x=0 (left/heated wall)

    # simpleGrading (gx, 1, 1):
    #   gx < 1 → cells shrink left→right (finer near x=0, the heated wall)
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
        (0           0           0)          // 0  bottom-left  front
        ({W_DOMAIN}  0           0)          // 1  bottom-right front
        ({W_DOMAIN}  {H_DOMAIN}  0)          // 2  top-right    front
        (0           {H_DOMAIN}  0)          // 3  top-left     front
        (0           0           {W_DEPTH})  // 4  bottom-left  back
        ({W_DOMAIN}  0           {W_DEPTH})  // 5  bottom-right back
        ({W_DOMAIN}  {H_DOMAIN}  {W_DEPTH}) // 6  top-right    back
        (0           {H_DOMAIN}  {W_DEPTH}) // 7  top-left     back
    );

    blocks
    (
        // grading: finer near x=0 (heated wall), uniform in y and z
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading ({gx} 1 1)
    );

    edges ( );

    boundary
    (
        heatedWall
        {{
            type wall;
            faces
            (
                (0 4 7 3)   // left face (x=0) — heated plate
            );
        }}
        openRight
        {{
            type patch;
            faces
            (
                (1 2 6 5)   // right face (x={W_DOMAIN}) — open/outlet
            );
        }}
        openTop
        {{
            type patch;
            faces
            (
                (3 7 6 2)   // top face (y={H_DOMAIN}) — open/outlet
            );
        }}
        bottomWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)   // bottom face (y=0) — adiabatic wall
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                (0 3 2 1)   // front face (z=0)
                (4 5 6 7)   // back face  (z={W_DEPTH})
            );
        }}
    );

    mergePatchPairs ( );
    """)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "system" / "blockMeshDict").write_text(content)


def write_control_dict(case_dir: Path, cfg: dict) -> None:
    end_time = cfg["endTime"]
    write_interval = end_time   # write only the final state
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }}

    application     buoyantBoussinesqSimpleFoam;
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
        div(phi,T)      bounded Gauss linearUpwind grad(T);
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
        p_rgh
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
        T
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
        pRefCell        0;
        pRefValue       0;
        residualControl
        {
            p_rgh           1e-5;
            U               1e-5;
            T               1e-5;
        }
    }

    relaxationFactors
    {
        fields
        {
            p_rgh           0.7;
        }
        equations
        {
            U               0.7;
            T               0.5;
        }
    }
    """)
    (case_dir / "system" / "fvSolution").write_text(content)


def write_initial_conditions(case_dir: Path) -> None:
    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)

    # U — heated wall and bottom wall: noSlip
    #     open boundaries (right, top): pressureInletOutletVelocity (handles both in/out)
    u_content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform (0 0 0);

    boundaryField
    {
        heatedWall
        {
            type            noSlip;
        }
        openRight
        {
            type            pressureInletOutletVelocity;
            value           uniform (0 0 0);
        }
        openTop
        {
            type            pressureInletOutletVelocity;
            value           uniform (0 0 0);
        }
        bottomWall
        {
            type            noSlip;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)

    # p_rgh — walls: fixedFluxPressure (no flux through wall)
    #          open boundaries: fixedValue 0 (acts as pressure outlet)
    p_rgh_content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p_rgh;
    }

    dimensions      [0 2 -2 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {
        heatedWall
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        openRight
        {
            type            fixedValue;
            value           uniform 0;
        }
        openTop
        {
            type            fixedValue;
            value           uniform 0;
        }
        bottomWall
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)

    # T — heated wall: Dirichlet 320K
    #     open boundaries: inletOutlet (ambient 300K if flow enters)
    #     bottom wall: zeroGradient (adiabatic)
    t_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      T;
    }}

    dimensions      [0 0 0 1 0 0 0];
    internalField   uniform {T_REF};

    boundaryField
    {{
        heatedWall
        {{
            type            fixedValue;
            value           uniform {T_HOT};
        }}
        openRight
        {{
            type            inletOutlet;
            inletValue      uniform {T_REF};
            value           uniform {T_REF};
        }}
        openTop
        {{
            type            inletOutlet;
            inletValue      uniform {T_REF};
            value           uniform {T_REF};
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

    # alphat — turbulent thermal diffusivity; required by buoyantBoussinesqSimpleFoam
    # even in laminar mode; set to zero (calculated) everywhere.
    alphat_content = textwrap.dedent("""\
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      alphat;
    }

    dimensions      [0 2 -1 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {
        heatedWall
        {
            type            calculated;
            value           uniform 0;
        }
        openRight
        {
            type            calculated;
            value           uniform 0;
        }
        openTop
        {
            type            calculated;
            value           uniform 0;
        }
        bottomWall
        {
            type            calculated;
            value           uniform 0;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)

    (zero_dir / "U").write_text(u_content)
    (zero_dir / "p_rgh").write_text(p_rgh_content)
    (zero_dir / "T").write_text(t_content)
    (zero_dir / "alphat").write_text(alphat_content)


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

    // Kinematic viscosity [m^2/s]
    nu              [0 2 -1 0 0 0 0] {nu:.6e};

    // Thermal expansion coefficient (Boussinesq) [1/K]
    beta            [0 0 0 -1 0 0 0] {BETA:.4e};

    // Reference temperature [K]
    TRef            [0 0 0 1 0 0 0] {T_REF};

    // Prandtl number [-]
    Pr              [0 0 0 0 0 0 0] {PR};

    // Turbulent Prandtl number [-]
    Prt             [0 0 0 0 0 0 0] 0.85;
    """)
    (case_dir / "constant").mkdir(parents=True, exist_ok=True)
    (case_dir / "constant" / "transportProperties").write_text(content)


def write_gravity(case_dir: Path, cfg: dict) -> None:
    gx, gy, gz = cfg["gravity"]
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       uniformDimensionedVectorField;
        object      g;
    }}

    dimensions      [0 1 -2 0 0 0 0];
    value           ({gx} {gy} {gz});
    """)
    (case_dir / "constant" / "g").write_text(content)


def write_turbulence_properties(case_dir: Path) -> None:
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


# ---------------------------------------------------------------------------
# Case orchestration
# ---------------------------------------------------------------------------

def setup_case(name: str, cfg: dict) -> Path:
    case_dir = CASE_BASE / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)

    write_blockmesh(case_dir, cfg)
    write_control_dict(case_dir, cfg)
    write_fvschemes(case_dir)
    write_fvsolution(case_dir)
    write_initial_conditions(case_dir)
    write_transport_properties(case_dir, cfg)
    write_gravity(case_dir, cfg)
    write_turbulence_properties(case_dir)

    # Determine error type
    error_type = None
    if "_E1_" in name:
        error_type = "E1_underconverged"
    elif "_E3_" in name:
        error_type = "E3_wrong_viscosity"
    elif "_E5_" in name:
        error_type = "E5_coarse_mesh"
    elif "_E8_" in name:
        error_type = "E8_gravity_flipped"

    nu = cfg["nu"]
    Ra = G * BETA * DT * H_DOMAIN**3 * PR / nu**2

    meta = {
        "case_name": name,
        "scenario": "S1",
        "scenario_name": "heated_vertical_plate_natural_convection",
        "solver": "buoyantBoussinesqSimpleFoam",
        "Ra": round(Ra, 1),
        "nu": nu,
        "Pr": PR,
        "beta": BETA,
        "T_ref": T_REF,
        "T_hot": T_HOT,
        "dT": DT,
        "domain_width": W_DOMAIN,
        "domain_height": H_DOMAIN,
        "gravity": list(cfg["gravity"]),
        "turbulence": False,
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "grading_x": cfg["grading_x"],
        "endTime": cfg["endTime"],
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name}  (Ra={Ra:.2e}, nu={nu:.4e}, mesh={cfg['nx']}x{cfg['ny']}, "
          f"g_y={cfg['gravity'][1]}, endTime={cfg['endTime']})")
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

    # buoyantBoussinesqSimpleFoam — allow non-zero for error cases
    r = docker_exec(
        case_dir,
        "buoyantBoussinesqSimpleFoam > log.buoyantBoussinesqSimpleFoam 2>&1",
        timeout=900,
    )
    if r.returncode != 0:
        print(f"  solver returned {r.returncode} for {name} (may be OK for error cases)")

    # foamToVTK — export latest time step
    r = docker_exec(case_dir, "foamToVTK -latestTime > log.foamToVTK 2>&1", timeout=600)
    if r.returncode != 0:
        print(f"  foamToVTK FAILED for {name}")
        return False

    print(f"  Completed: {name}")
    return True


def main() -> None:
    CASE_BASE.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("CFD Visual QA — S1: Heated Vertical Plate, Natural Convection BL")
    print("=" * 60)

    results: dict[str, bool] = {}
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
