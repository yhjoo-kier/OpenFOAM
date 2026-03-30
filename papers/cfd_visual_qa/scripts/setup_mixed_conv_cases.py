#!/usr/bin/env python3
"""Setup mixed convection channel flow (S7) OpenFOAM cases for CFD Visual QA benchmark.

Physics: buoyantBoussinesqSimpleFoam — incompressible with Boussinesq buoyancy.
Geometry: 10 m × 1 m channel, 0.1 m deep (2-D, 1 cell in z, empty BC).
Heating:  bottom wall hot (320 K), top wall adiabatic; inlet at 300 K.

Richardson number:  Ri = Gr / Re^2   (1 → balanced, <0.1 → forced-dominated)
  Gr = g * beta * dT * H^3 / nu^2 = Re^2 for Ri = 1
  nu = sqrt(g * beta * dT * H^3 / Re^2)

For Ri≈1 (Re=500):
  nu = sqrt(9.81 * 3.4e-3 * 20 * 1^3 / 250000) = sqrt(2.6604e-6) ≈ 0.001633 m^2/s
  U_inlet = Re * nu / H = 500 * 0.001633 / 1 = 0.8165 m/s

For Ri≈0.1 (Re≈1581):
  Ri = Gr/Re^2 = 0.1 → Re = sqrt(Gr/0.1) = sqrt(2500000) ≈ 1581
  U_inlet = Re * nu / H = 1581 * 0.001633 / 1 ≈ 2.582 m/s
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

# Geometry
L = 10.0   # channel length (m)
H = 1.0    # channel height (m)
W = 0.1    # z-width (2-D, 1 cell)

# Thermal / fluid properties
G = 9.81
BETA = 3.4e-3    # thermal expansion coefficient [1/K]
DT = 20.0        # T_bottom - T_inlet [K]
T_INLET = 300.0
T_BOTTOM = 320.0
T_REF = 300.0
PR = 0.71

# Ri=1: Re=500 → nu = sqrt(g*beta*dT*H^3/Re^2)
RE_RI1 = 500.0
NU_BASE = math.sqrt(G * BETA * DT * H**3 / RE_RI1**2)   # ≈ 0.001633
U_RI1 = RE_RI1 * NU_BASE / H                              # ≈ 0.8165

# Ri=0.1: same nu, Re increased by sqrt(10)
RE_RI01 = RE_RI1 * math.sqrt(10.0)    # ≈ 1580.9
U_RI01 = RE_RI01 * NU_BASE / H        # ≈ 2.582

CASES: dict[str, dict] = {
    # -------------------------------------------------------------------------
    # Correct cases
    # -------------------------------------------------------------------------
    "S7_correct_Ri1": {
        "Re": RE_RI1,
        "U_inlet": U_RI1,
        "nu": NU_BASE,
        "nx": 200, "ny": 40,
        "grading_y": 1.0,     # uniform — Ri=1, modest BL
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: Ri≈1 balanced mixed convection, Re=500, 200×40 mesh, "
            "buoyancy and forced flow in equilibrium"
        ),
    },
    "S7_correct_Ri01": {
        "Re": RE_RI01,
        "U_inlet": U_RI01,
        "nu": NU_BASE,
        "nx": 200, "ny": 40,
        "grading_y": 4.0,     # wall-normal grading for stronger BL at higher Re
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: Ri≈0.1 forced-dominated mixed convection, Re≈1581, "
            "200×40 wall-graded mesh, forced convection dominates"
        ),
    },
    # -------------------------------------------------------------------------
    # Error cases
    # -------------------------------------------------------------------------
    "S7_E1_underconverged": {
        "Re": RE_RI1,
        "U_inlet": U_RI1,
        "nu": NU_BASE,
        "nx": 200, "ny": 40,
        "grading_y": 1.0,
        "endTime": 30,        # ERROR: only 30 iterations — far from steady state
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E1: under-converged — stopped at 30 iterations "
            "(Ri≈1 setup, correct BCs, wrong endTime)"
        ),
    },
    "S7_E8_gravity_flipped": {
        "Re": RE_RI1,
        "U_inlet": U_RI1,
        "nu": NU_BASE,
        "nx": 200, "ny": 40,
        "grading_y": 1.0,
        "endTime": 10000,
        "gravity": (0.0, +9.81, 0.0),  # ERROR: g points up — buoyancy pushes hot fluid DOWN
        "description": (
            "Error E8: gravity flipped — g=(0, +9.81, 0); "
            "buoyancy drives hot bottom fluid downward instead of upward"
        ),
    },
    "S7_E5_coarse_mesh": {
        "Re": RE_RI1,
        "U_inlet": U_RI1,
        "nu": NU_BASE,
        "nx": 20, "ny": 4,    # ERROR: 20×4 far too coarse — BL unresolved
        "grading_y": 1.0,
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E5: coarse mesh — only 20×4 cells for Ri≈1; "
            "thermal and velocity BLs completely unresolved"
        ),
    },
    "S7_E3_wrong_viscosity": {
        "Re": RE_RI1,
        "U_inlet": U_RI1,
        "nu": 0.1,            # ERROR: nu=0.1 (≈60× too high); kills forced + buoyancy
        "nx": 200, "ny": 40,
        "grading_y": 1.0,
        "endTime": 10000,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E3: wrong viscosity — nu=0.1 m²/s (≈60× too high); "
            "effective Re≈8 and Gr≈2700, eliminating both forced and buoyancy effects"
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
    """10m × 1m channel with optional wall-normal (y) grading."""
    nx = cfg["nx"]
    ny = cfg["ny"]
    gy = cfg["grading_y"]

    # simpleGrading (x y z): uniform in x, graded in y toward bottom wall
    # gy > 1 means cells get smaller toward the top face; invert to refine at bottom.
    # We want refinement at BOTH walls — use edgeGrading for symmetric refinement.
    # For simplicity (matching S6/S2 approach), use simpleGrading with gy on y.
    # gy > 1 compresses cells toward the top; gy < 1 compresses toward bottom.
    # Bottom heated wall needs refinement → use 1/gy for y so small cells are at y=0.
    inv_gy = round(1.0 / gy, 6) if gy != 1.0 else 1.0
    grading_str = f"simpleGrading (1 {inv_gy} 1)"

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
        (0     0    0  )   // 0  bottom-left  front
        ({L}   0    0  )   // 1  bottom-right front
        ({L}   {H}  0  )   // 2  top-right    front
        (0     {H}  0  )   // 3  top-left     front
        (0     0    {W})   // 4  bottom-left  back
        ({L}   0    {W})   // 5  bottom-right back
        ({L}   {H}  {W})   // 6  top-right    back
        (0     {H}  {W})   // 7  top-left     back
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) {grading_str}
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
                (3 7 6 2)   // top face (y={H}) — adiabatic
            );
        }}
        bottomWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)   // bottom face (y=0) — heated
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
    # Use writeInterval 500 to avoid writing only at endTime issue;
    # purgeWrite 1 keeps only the latest written time step.
    write_interval = 500
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
    # Channel has inlet+outlet — outlet fixes pressure, no pRefCell needed.
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
    U_in = cfg["U_inlet"]

    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)

    # U — fixedValue at inlet, pressureInletOutletVelocity at outlet, noSlip on walls
    u_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }}

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform ({U_in:.6f} 0 0);

    boundaryField
    {{
        inlet
        {{
            type            fixedValue;
            value           uniform ({U_in:.6f} 0 0);
        }}
        outlet
        {{
            type            pressureInletOutletVelocity;
            value           uniform ({U_in:.6f} 0 0);
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

    # p_rgh — fixedFluxPressure at inlet/walls, fixedValue 0 at outlet
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
        inlet
        {
            type            fixedFluxPressure;
            value           uniform 0;
        }
        outlet
        {
            type            fixedValue;
            value           uniform 0;
        }
        topWall
        {
            type            fixedFluxPressure;
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

    # T — fixedValue at inlet (300 K), inletOutlet at outlet, fixedValue at bottom (320 K), zeroGradient at top
    t_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      T;
    }}

    dimensions      [0 0 0 1 0 0 0];
    internalField   uniform {T_INLET};

    boundaryField
    {{
        inlet
        {{
            type            fixedValue;
            value           uniform {T_INLET};
        }}
        outlet
        {{
            type            inletOutlet;
            inletValue      uniform {T_INLET};
            value           uniform {T_INLET};
        }}
        topWall
        {{
            type            zeroGradient;
        }}
        bottomWall
        {{
            type            fixedValue;
            value           uniform {T_BOTTOM};
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    # alphat — turbulent thermal diffusivity; zero for laminar
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
    write_initial_conditions(case_dir, cfg)
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
    Re = cfg["Re"]
    Gr = G * BETA * DT * H**3 / nu**2
    Ri = Gr / Re**2

    meta = {
        "case_name": name,
        "scenario": "S7",
        "scenario_name": "mixed_convection_heated_channel",
        "solver": "buoyantBoussinesqSimpleFoam",
        "Re": round(Re, 2),
        "Gr": round(Gr, 1),
        "Ri": round(Ri, 4),
        "U_inlet": round(cfg["U_inlet"], 6),
        "nu": nu,
        "Pr": PR,
        "beta": BETA,
        "T_ref": T_REF,
        "T_inlet": T_INLET,
        "T_bottom": T_BOTTOM,
        "dT": DT,
        "gravity": list(cfg["gravity"]),
        "turbulence": False,
        "mesh_nx": cfg["nx"],
        "mesh_ny": cfg["ny"],
        "grading_y": cfg["grading_y"],
        "channel_length": L,
        "channel_height": H,
        "endTime": cfg["endTime"],
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name}  (Re={Re:.0f}, Ri={Ri:.3f}, nu={nu:.4e}, "
          f"U={cfg['U_inlet']:.4f}, mesh={cfg['nx']}x{cfg['ny']}, g_y={cfg['gravity'][1]})")
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
    print("CFD Visual QA — S7: Mixed Convection, Heated Channel Flow")
    print(f"  nu_base = {NU_BASE:.6f} m^2/s")
    print(f"  Ri=1:  Re={RE_RI1:.0f}, U={U_RI1:.4f} m/s")
    print(f"  Ri=0.1: Re={RE_RI01:.0f}, U={U_RI01:.4f} m/s")
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
