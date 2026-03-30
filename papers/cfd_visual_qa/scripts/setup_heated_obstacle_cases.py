#!/usr/bin/env python3
"""Setup heated bump/obstacle (S8) OpenFOAM cases for CFD Visual QA benchmark.

Physics: buoyantBoussinesqSimpleFoam — forced convection over a heated bump on
the channel bottom wall, with thermal plume and wake dynamics.

Geometry (2D, z-width=0.1m, 1 cell, empty BC):
  Channel: 10m long × 2m tall
  Bump: 0.5m wide × 0.5m tall, positioned at x=3..3.5m on the bottom wall

5-block blockMesh layout:
  Block 0a: x=[0,3],   y=[0,0.5]    — upstream lower
  Block 0b: x=[0,3],   y=[0.5,2]    — upstream upper
  Block 1:  x=[3,3.5], y=[0.5,2]    — above bump
  Block 2a: x=[3.5,10],y=[0,0.5]    — downstream lower
  Block 2b: x=[3.5,10],y=[0.5,2]    — downstream upper

Boundaries:
  inlet       — left face of 0a + 0b
  outlet      — right face of 2a + 2b
  bottomWall  — bottom of 0a + bottom of 2a  (adiabatic)
  topWall     — top of 0b + top of 1 + top of 2b  (adiabatic)
  bumpSurface — top of Block 1 (y=0.5) + right of 0a (x=3, y=0..0.5)
                + left of 2a (x=3.5, y=0..0.5)  (heated T=320K)
  frontAndBack — empty

Physics parameters:
  U_inlet = 0.5 m/s, T_inlet = 300K, T_bump = 320K
  Re (bump height h=0.5m): Re = U*h/nu
  For Re=500: nu = 0.5*0.5/500 = 5e-4 m^2/s
  beta=3.4e-3 1/K, Pr=0.71
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

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------
L_CH    = 10.0    # channel total length [m]
H_CH    = 2.0     # channel height [m]
X_BL    = 3.0     # bump left x [m]
X_BR    = 3.5     # bump right x [m]
H_BMP   = 0.5     # bump height [m]
W       = 0.1     # z-width (2D, 1 cell) [m]

# Derived lengths
L_UP    = X_BL            # upstream length
L_BUMP  = X_BR - X_BL    # bump footprint length
L_DOWN  = L_CH - X_BR    # downstream length
H_ABOVE = H_CH - H_BMP   # height above bump top

# Physical parameters
G      = 9.81
BETA   = 3.4e-3   # thermal expansion [1/K]
T_REF  = 300.0    # ambient / inlet temperature [K]
T_HOT  = 320.0    # bump surface temperature [K]
DT     = T_HOT - T_REF
PR     = 0.71     # Prandtl number
U_IN   = 0.5      # inlet velocity [m/s]

# ---------------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------------
CASES: dict[str, dict] = {
    # -----------------------------------------------------------------------
    # Correct cases
    # -----------------------------------------------------------------------
    "S8_correct_lam": {
        "nu": U_IN * H_BMP / 500,   # Re=500
        # cell counts per block region
        "nx_up": 40, "nx_bump": 10, "nx_down": 130,
        "ny_low": 20, "ny_high": 60,
        "endTime": 10000, "writeInterval": 500,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: forced convection over heated bump, Re=500, "
            "5-block mesh, laminar, thermal plume and wake"
        ),
    },
    "S8_correct_fine": {
        "nu": U_IN * H_BMP / 500,   # Re=500
        "nx_up": 60, "nx_bump": 15, "nx_down": 200,
        "ny_low": 30, "ny_high": 90,
        "endTime": 10000, "writeInterval": 500,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Correct: forced convection over heated bump, Re=500, "
            "finer 5-block mesh, laminar"
        ),
    },
    # -----------------------------------------------------------------------
    # Error cases
    # -----------------------------------------------------------------------
    "S8_E1_underconverged": {
        "nu": U_IN * H_BMP / 500,
        "nx_up": 40, "nx_bump": 10, "nx_down": 130,
        "ny_low": 20, "ny_high": 60,
        "endTime": 30, "writeInterval": 30,   # ERROR: only 30 iterations
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E1: under-converged — stopped at 30 iterations, "
            "wake and thermal plume not established"
        ),
    },
    "S8_E8_gravity_flipped": {
        "nu": U_IN * H_BMP / 500,
        "nx_up": 40, "nx_bump": 10, "nx_down": 130,
        "ny_low": 20, "ny_high": 60,
        "endTime": 10000, "writeInterval": 500,
        "gravity": (0.0, +9.81, 0.0),   # ERROR: buoyancy reversed
        "description": (
            "Error E8: gravity flipped — g=(0,+9.81,0), "
            "thermal plume descends instead of rising from bump"
        ),
    },
    "S8_E5_coarse_mesh": {
        "nu": U_IN * H_BMP / 500,
        "nx_up": 8, "nx_bump": 2, "nx_down": 26,
        "ny_low": 4, "ny_high": 12,   # ERROR: very coarse — wake unresolved
        "endTime": 10000, "writeInterval": 500,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E5: coarse mesh — 8+2+26 x, 4+12 y, "
            "wake and thermal boundary layer completely unresolved"
        ),
    },
    "S8_E3_wrong_viscosity": {
        "nu": 0.05,    # ERROR: 100× too high → Re_eff ≈ 5, kills wake
        "nx_up": 40, "nx_bump": 10, "nx_down": 130,
        "ny_low": 20, "ny_high": 60,
        "endTime": 10000, "writeInterval": 500,
        "gravity": (0.0, -9.81, 0.0),
        "description": (
            "Error E3: wrong viscosity — nu=0.05 (100× too high), "
            "Re_eff≈5, creeping flow — no wake, no thermal plume structure"
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
    """Write 5-block blockMeshDict for heated-bump geometry.

    Vertex numbering (front face z=0, then back face z=W):
      Block 0a (upstream lower):   v0-v3   + v20-v23
      Block 0b (upstream upper):   v4-v7   + v24-v27
      Block 1  (above bump):       v8-v11  + v28-v31
      Block 2a (downstream lower): v12-v15 + v32-v35
      Block 2b (downstream upper): v16-v19 + v36-v39

    Shared edges between blocks are enforced by using the same vertex indices.
    Block 0b top-right (x=3, y=2): v7/v27 shared with Block 1 top-left.
    Block 0a top-right (x=3, y=0.5): v3/v23 shared with Block 0b bottom-right
      and Block 1 bottom-left.
    """
    nx_up   = cfg["nx_up"]
    nx_bump = cfg["nx_bump"]
    nx_down = cfg["nx_down"]
    ny_low  = cfg["ny_low"]
    ny_high = cfg["ny_high"]

    # Vertex coordinates (front face, z=0):
    # Block 0a: upstream lower (x: 0..L_UP, y: 0..H_BMP)
    #   v0=(0,0,0)  v1=(L_UP,0,0)  v2=(L_UP,H_BMP,0)  v3=(0,H_BMP,0)
    # Block 0b: upstream upper (x: 0..L_UP, y: H_BMP..H_CH)
    #   v4=(0,H_BMP,0)==v3  v5=(L_UP,H_BMP,0)==v2  v6=(L_UP,H_CH,0)  v7=(0,H_CH,0)
    #   → share v2,v3 → we reuse indices
    # Block 1: above bump (x: L_UP..X_BR, y: H_BMP..H_CH)
    #   v8=(L_UP,H_BMP,0)==v2  v9=(X_BR,H_BMP,0)  v10=(X_BR,H_CH,0)  v11=(L_UP,H_CH,0)==v6
    # Block 2a: downstream lower (x: X_BR..L_CH, y: 0..H_BMP)
    #   v12=(X_BR,0,0)  v13=(L_CH,0,0)  v14=(L_CH,H_BMP,0)  v15=(X_BR,H_BMP,0)==v9
    # Block 2b: downstream upper (x: X_BR..L_CH, y: H_BMP..H_CH)
    #   v16=(X_BR,H_BMP,0)==v9==v15  v17=(L_CH,H_BMP,0)==v14  v18=(L_CH,H_CH,0)  v19=(X_BR,H_CH,0)==v10

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
        // Block 0a: upstream lower  x=[0,{L_UP}], y=[0,{H_BMP}]
        ( 0       0        0   )  // 0   BL-front
        ( {L_UP}  0        0   )  // 1   BR-front
        ( {L_UP}  {H_BMP}  0   )  // 2   TR-front  (shared with 0b, Block 1)
        ( 0       {H_BMP}  0   )  // 3   TL-front  (shared with 0b)
        ( 0       0        {W} )  // 4   BL-back
        ( {L_UP}  0        {W} )  // 5   BR-back
        ( {L_UP}  {H_BMP}  {W} )  // 6   TR-back   (shared with 0b, Block 1)
        ( 0       {H_BMP}  {W} )  // 7   TL-back   (shared with 0b)

        // Block 0b: upstream upper  x=[0,{L_UP}], y=[{H_BMP},{H_CH}]
        // bottom vertices reuse 3,2,6,7 above
        ( 0       {H_CH}   0   )  // 8   TL-front
        ( {L_UP}  {H_CH}   0   )  // 9   TR-front  (shared with Block 1)
        ( {L_UP}  {H_CH}   {W} )  // 10  TR-back   (shared with Block 1)
        ( 0       {H_CH}   {W} )  // 11  TL-back

        // Block 1: above bump  x=[{L_UP},{X_BR}], y=[{H_BMP},{H_CH}]
        // left-bottom reuses 2,6; left-top reuses 9,10
        ( {X_BR}  {H_BMP}  0   )  // 12  BR-front  (shared with 2a, 2b)
        ( {X_BR}  {H_CH}   0   )  // 13  TR-front  (shared with 2b)
        ( {X_BR}  {H_CH}   {W} )  // 14  TR-back   (shared with 2b)
        ( {X_BR}  {H_BMP}  {W} )  // 15  BR-back   (shared with 2a, 2b)

        // Block 2a: downstream lower  x=[{X_BR},{L_CH}], y=[0,{H_BMP}]
        // left vertices reuse 12,15
        ( {X_BR}  0        0   )  // 16  BL-front
        ( {L_CH}  0        0   )  // 17  BR-front
        ( {L_CH}  {H_BMP}  0   )  // 18  TR-front  (shared with 2b)
        ( {X_BR}  0        {W} )  // 19  BL-back
        ( {L_CH}  0        {W} )  // 20  BR-back
        ( {L_CH}  {H_BMP}  {W} )  // 21  TR-back   (shared with 2b)

        // Block 2b: downstream upper  x=[{X_BR},{L_CH}], y=[{H_BMP},{H_CH}]
        // left-bottom reuses 12,15; left-top reuses 13,14; bottom-right reuses 18,21
        ( {L_CH}  {H_CH}   0   )  // 22  TR-front
        ( {L_CH}  {H_CH}   {W} )  // 23  TR-back
    );

    blocks
    (
        // Block 0a: upstream lower
        hex (0 1 2 3 4 5 6 7)          ({nx_up}   {ny_low}  1) simpleGrading (1 1 1)
        // Block 0b: upstream upper
        hex (3 2 9 8 7 6 10 11)        ({nx_up}   {ny_high} 1) simpleGrading (1 1 1)
        // Block 1: above bump
        hex (2 12 13 9 6 15 14 10)     ({nx_bump} {ny_high} 1) simpleGrading (1 1 1)
        // Block 2a: downstream lower
        hex (16 17 18 12 19 20 21 15)  ({nx_down} {ny_low}  1) simpleGrading (1 1 1)
        // Block 2b: downstream upper
        // left-bottom=12,15; left-top=13,14 (shared with Block 1 right face)
        hex (12 18 22 13 15 21 23 14)  ({nx_down} {ny_high} 1) simpleGrading (1 1 1)
    );

    edges ( );

    boundary
    (
        inlet
        {{
            type patch;
            faces
            (
                (0 3 7 4)    // left of 0a (x=0, y=0..{H_BMP})
                (3 8 11 7)   // left of 0b (x=0, y={H_BMP}..{H_CH})
            );
        }}
        outlet
        {{
            type patch;
            faces
            (
                (18 17 20 21)   // right of 2a (x={L_CH}, y=0..{H_BMP})
                (22 18 21 23)   // right of 2b (x={L_CH}, y={H_BMP}..{H_CH})
            );
        }}
        bottomWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)     // bottom of 0a (y=0, x=0..{L_UP})
                (16 17 20 19) // bottom of 2a (y=0, x={X_BR}..{L_CH})
            );
        }}
        topWall
        {{
            type wall;
            faces
            (
                (8 9 10 11)   // top of 0b  (y={H_CH}, x=0..{L_UP})
                (9 13 14 10)  // top of Block 1 (y={H_CH}, x={L_UP}..{X_BR})
                (13 22 23 14) // top of 2b  (y={H_CH}, x={X_BR}..{L_CH})
            );
        }}
        bumpSurface
        {{
            type wall;
            faces
            (
                (1 2 6 5)      // left face of bump  (x={L_UP},  y=0..{H_BMP})
                (2 12 15 6)    // top face of bump   (y={H_BMP}, x={L_UP}..{X_BR})
                (12 16 19 15)  // right face of bump (x={X_BR},  y=0..{H_BMP})
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                // front (z=0)
                (0 1 2 3)      // 0a front
                (3 2 9 8)      // 0b front
                (2 12 13 9)    // Block 1 front
                (16 17 18 12)  // 2a front
                (12 18 22 13)  // 2b front
                // back (z={W})
                (4 7 6 5)      // 0a back
                (7 11 10 6)    // 0b back
                (6 10 14 15)   // Block 1 back
                (19 15 21 20)  // 2a back
                (15 14 23 21)  // 2b back
            );
        }}
    );

    mergePatchPairs ( );
    """)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "system" / "blockMeshDict").write_text(content)


def write_control_dict(case_dir: Path, cfg: dict) -> None:
    end_time = cfg["endTime"]
    write_interval = cfg.get("writeInterval", end_time)
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
    zero_dir = case_dir / "0"
    zero_dir.mkdir(parents=True, exist_ok=True)

    # U — inlet: fixedValue; outlet: zeroGradient; walls: noSlip
    u_content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }}

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform ({U_IN} 0 0);

    boundaryField
    {{
        inlet
        {{
            type            fixedValue;
            value           uniform ({U_IN} 0 0);
        }}
        outlet
        {{
            type            zeroGradient;
        }}
        bottomWall
        {{
            type            noSlip;
        }}
        topWall
        {{
            type            noSlip;
        }}
        bumpSurface
        {{
            type            noSlip;
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    # p_rgh — walls: fixedFluxPressure; inlet: fixedFluxPressure; outlet: fixedValue 0
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
        bumpSurface
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

    # T — inlet: fixedValue T_REF; outlet: zeroGradient;
    #     bottomWall/topWall: zeroGradient (adiabatic);
    #     bumpSurface: fixedValue T_HOT
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
        inlet
        {{
            type            fixedValue;
            value           uniform {T_REF};
        }}
        outlet
        {{
            type            zeroGradient;
        }}
        bottomWall
        {{
            type            zeroGradient;
        }}
        topWall
        {{
            type            zeroGradient;
        }}
        bumpSurface
        {{
            type            fixedValue;
            value           uniform {T_HOT};
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    # alphat — turbulent thermal diffusivity (required even in laminar mode)
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
        bumpSurface
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
    Re = U_IN * H_BMP / nu

    meta = {
        "case_name": name,
        "scenario": "S8",
        "scenario_name": "heated_bump_forced_convection",
        "solver": "buoyantBoussinesqSimpleFoam",
        "Re": round(Re, 1),
        "nu": nu,
        "U_inlet": U_IN,
        "Pr": PR,
        "beta": BETA,
        "T_ref": T_REF,
        "T_bump": T_HOT,
        "dT": DT,
        "channel_length": L_CH,
        "channel_height": H_CH,
        "bump_x_start": X_BL,
        "bump_x_end": X_BR,
        "bump_height": H_BMP,
        "gravity": list(cfg["gravity"]),
        "turbulence": False,
        "mesh_nx_up": cfg["nx_up"],
        "mesh_nx_bump": cfg["nx_bump"],
        "mesh_nx_down": cfg["nx_down"],
        "mesh_ny_low": cfg["ny_low"],
        "mesh_ny_high": cfg["ny_high"],
        "endTime": cfg["endTime"],
        "description": cfg["description"],
        "is_correct": error_type is None,
        "error_type": error_type,
    }
    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name}  (Re={Re:.0f}, nu={nu:.4e}, "
          f"mesh={cfg['nx_up']}+{cfg['nx_bump']}+{cfg['nx_down']}x"
          f"{cfg['ny_low']}+{cfg['ny_high']}, "
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
    print("CFD Visual QA — S8: Heated Bump, Forced Convection with Wake")
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
