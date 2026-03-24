#!/usr/bin/env python3
"""Setup Rayleigh-Bénard convection (S9) OpenFOAM cases
for CFD Visual QA pilot.

Physics: buoyantBoussinesqSimpleFoam — incompressible with Boussinesq buoyancy.
Geometry: 4 m wide × 1 m tall cavity, 0.1 m deep (2-D, 1 cell in z, empty BC).
          Aspect ratio 4:1 allows multiple convection cells to form.
Heating: bottom wall hot (310 K), top wall cold (300 K).  ΔT = 10 K, T_ref = 305 K.

Rayleigh number:
  Ra = g * beta * dT * H^3 * Pr / nu^2
  -> nu = sqrt(g * beta * dT * H^3 * Pr / Ra)
  Ra = 5000 : nu = sqrt(9.81 * 3.4e-3 * 10 * 1^3 * 0.71 / 5000)  ≈ 0.006880
  Ra = 20000: nu = sqrt(9.81 * 3.4e-3 * 10 * 1^3 * 0.71 / 20000) ≈ 0.003440
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

W_CAVITY = 4.0  # cavity width  (x, m)
H_CAVITY = 1.0  # cavity height (y, m)
D_CAVITY = 0.1  # z-depth       (m) — 2-D, 1 cell

G = 9.81
BETA = 3.4e-3   # thermal expansion coefficient [1/K]
DT = 10.0       # T_hot - T_cold  [K]
T_HOT = 310.0   # bottom wall temperature [K]
T_COLD = 300.0  # top wall temperature    [K]
T_REF = 305.0   # reference temperature   [K]
PR = 0.71       # Prandtl number


def _nu_for_Ra(Ra: float) -> float:
    """Return kinematic viscosity that gives the requested Rayleigh number.

    Ra = g * beta * dT * H^3 * Pr / nu^2
    => nu = sqrt(g * beta * dT * H^3 * Pr / Ra)
    """
    return math.sqrt(G * BETA * DT * H_CAVITY**3 * PR / Ra)


CASES: dict[str, dict] = {
    # -------------------------------------------------------------------------
    # Correct cases
    # -------------------------------------------------------------------------
    "S9_correct_Ra5k": {
        "Ra": 5000,
        "nu": _nu_for_Ra(5000),    # ≈ 0.006880
        "nx": 160, "ny": 40,
        "grading_y": 1.0,           # uniform — Ra low enough
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "hot_bottom": True,         # canonical: bottom=hot (310 K), top=cold (300 K)
        "description": "Correct: Ra=5000 laminar Rayleigh-Benard convection, 160x40 mesh",
    },
    "S9_correct_Ra20k": {
        "Ra": 20000,
        "nu": _nu_for_Ra(20000),   # ≈ 0.003440
        "nx": 200, "ny": 50,
        "grading_y": 4.0,           # refine near top/bottom walls (boundary layers)
        "endTime": 15000,
        "gravity": (0.0, -9.81, 0.0),
        "hot_bottom": True,
        "description": "Correct: Ra=20000 laminar Rayleigh-Benard convection, 200x50 wall-graded mesh",
    },
    # -------------------------------------------------------------------------
    # Error cases
    # -------------------------------------------------------------------------
    "S9_E1_underconverged": {
        "Ra": 5000,
        "nu": _nu_for_Ra(5000),
        "nx": 160, "ny": 40,
        "grading_y": 1.0,
        "endTime": 30,              # ERROR: only 30 iterations — far from steady state
        "gravity": (0.0, -9.81, 0.0),
        "hot_bottom": True,
        "description": "Error E1: under-converged — stopped at 30 iterations (Ra=5000)",
    },
    "S9_E8_gravity_flipped": {
        "Ra": 5000,
        "nu": _nu_for_Ra(5000),
        "nx": 160, "ny": 40,
        "grading_y": 1.0,
        "endTime": 10000,
        "gravity": (0.0, +9.81, 0.0),  # ERROR: gravity points up — stable stratification
        "hot_bottom": True,
        "description": "Error E8: gravity flipped — g=(0, +9.81, 0), no convection cells form",
    },
    "S9_E5_coarse_mesh": {
        "Ra": 20000,
        "nu": _nu_for_Ra(20000),
        "nx": 16, "ny": 4,          # ERROR: 16x4 far too coarse for Ra=20000
        "grading_y": 1.0,
        "endTime": 15000,
        "gravity": (0.0, -9.81, 0.0),
        "hot_bottom": True,
        "description": "Error E5: coarse mesh — only 16x4 for Ra=20000",
    },
    "S9_E2_bc_swap": {
        "Ra": 5000,
        "nu": _nu_for_Ra(5000),
        "nx": 160, "ny": 40,
        "grading_y": 1.0,
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "hot_bottom": False,        # ERROR: hot on top, cold on bottom — stable stratification
        "description": "Error E2: BC swap — hot on top (310 K), cold on bottom (300 K)",
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
    gy = cfg["grading_y"]

    # For vertical (y) direction we refine near both top and bottom walls.
    # simpleGrading applies ratio front->back (bottom->top here), so with
    # gy > 1 cells cluster toward the top wall.  For symmetric wall refinement
    # a two-block approach is ideal; we use simpleGrading and accept slight
    # asymmetry — the physics test remains meaningful.
    grading_str = f"simpleGrading (1 {gy} 1)"

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
        (0          0           0)           // 0  bottom-left  front
        ({W_CAVITY} 0           0)           // 1  bottom-right front
        ({W_CAVITY} {H_CAVITY}  0)           // 2  top-right    front
        (0          {H_CAVITY}  0)           // 3  top-left     front
        (0          0           {D_CAVITY})  // 4  bottom-left  back
        ({W_CAVITY} 0           {D_CAVITY})  // 5  bottom-right back
        ({W_CAVITY} {H_CAVITY}  {D_CAVITY})  // 6  top-right    back
        (0          {H_CAVITY}  {D_CAVITY})  // 7  top-left     back
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) {grading_str}
    );

    edges ( );

    boundary
    (
        bottomWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)   // bottom face (y=0)
            );
        }}
        topWall
        {{
            type wall;
            faces
            (
                (3 7 6 2)   // top face (y={H_CAVITY})
            );
        }}
        leftWall
        {{
            type wall;
            faces
            (
                (0 4 7 3)   // left face (x=0)
            );
        }}
        rightWall
        {{
            type wall;
            faces
            (
                (1 2 6 5)   // right face (x={W_CAVITY})
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                (0 3 2 1)   // front face (z=0)
                (4 5 6 7)   // back face  (z={D_CAVITY})
            );
        }}
    );

    mergePatchPairs ( );
    """)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "system" / "blockMeshDict").write_text(content)


def write_control_dict(case_dir: Path, cfg: dict) -> None:
    end_time = cfg["endTime"]
    # Use writeInterval = min(500, endTime) so short error cases (e.g. E1 with
    # endTime=30) still produce at least one output time directory.
    write_interval = min(500, end_time)
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


def write_initial_conditions(case_dir: Path, cfg: dict) -> None:
    hot_bottom = cfg["hot_bottom"]
    T_bottom = T_HOT  if hot_bottom else T_COLD
    T_top    = T_COLD if hot_bottom else T_HOT

    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)

    # U — all walls noSlip (closed cavity, buoyancy-driven only)
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
        bottomWall
        {
            type            noSlip;
        }
        topWall
        {
            type            noSlip;
        }
        leftWall
        {
            type            noSlip;
        }
        rightWall
        {
            type            noSlip;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    """)

    # p_rgh — buoyantBoussinesqSimpleFoam uses p_rgh (reduced pressure)
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
        bottomWall
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        topWall
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        leftWall
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        rightWall
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

    # T — Dirichlet on hot/cold walls (top & bottom), zeroGradient on side walls
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
        bottomWall
        {{
            type            fixedValue;
            value           uniform {T_bottom};
        }}
        topWall
        {{
            type            fixedValue;
            value           uniform {T_top};
        }}
        leftWall
        {{
            type            zeroGradient;
        }}
        rightWall
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
    # even in laminar mode.  Set to zero everywhere (no turbulent transport).
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
        bottomWall
        {
            type            calculated;
            value           uniform 0;
        }
        topWall
        {
            type            calculated;
            value           uniform 0;
        }
        leftWall
        {
            type            calculated;
            value           uniform 0;
        }
        rightWall
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
    write_initial_conditions(case_dir, cfg)
    write_transport_properties(case_dir, cfg)
    write_gravity(case_dir, cfg)
    write_turbulence_properties(case_dir)

    # Determine error type
    error_type = None
    if "_E1_" in name:
        error_type = "E1_underconverged"
    elif "_E2_" in name:
        error_type = "E2_bc_swap"
    elif "_E5_" in name:
        error_type = "E5_coarse_mesh"
    elif "_E8_" in name:
        error_type = "E8_gravity_flipped"

    Ra = cfg["Ra"]
    nu = cfg["nu"]
    alpha = nu / PR
    Ra_check = G * BETA * DT * H_CAVITY**3 / (nu * alpha)

    meta = {
        "case_name": name,
        "scenario": "S9",
        "scenario_name": "rayleigh_benard_convection",
        "solver": "buoyantBoussinesqSimpleFoam",
        "Ra": Ra,
        "Ra_check": round(Ra_check, 1),
        "nu": nu,
        "Pr": PR,
        "beta": BETA,
        "T_ref": T_REF,
        "T_hot": T_HOT,
        "T_cold": T_COLD,
        "dT": DT,
        "cavity_width": W_CAVITY,
        "cavity_height": H_CAVITY,
        "aspect_ratio": W_CAVITY / H_CAVITY,
        "gravity": list(cfg["gravity"]),
        "hot_bottom": cfg["hot_bottom"],
        "turbulence": False,
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "grading_y": cfg["grading_y"],
        "endTime": cfg["endTime"],
        "writeInterval": min(500, cfg["endTime"]),
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name}  (Ra={Ra}, nu={nu:.5f}, mesh={cfg['nx']}x{cfg['ny']}, "
          f"hot_bottom={cfg['hot_bottom']}, g_y={cfg['gravity'][1]})")
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

    # buoyantBoussinesqSimpleFoam — allow non-zero return for error cases
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
    print("CFD Visual QA — S9: Rayleigh-Benard Convection")
    print("Geometry: 4m wide x 1m tall, bottom-heated (310 K), top-cooled (300 K)")
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
