#!/usr/bin/env python3
"""Setup ventilated room (Nielsen 1990-style) OpenFOAM cases for CFD Visual QA benchmark.

Geometry (2D):
  Room: 9m long × 3m tall
  Slot inlet at top-left:  x=0,   y=2.832..3.0  (width h=0.168m)
  Slot outlet at bottom-right: x=9, y=0..0.168    (width h=0.168m)
  Depth (z): 0.1m, 1 cell, empty BC

4-block layout (clean inlet/outlet patches):
  Block A: x=0..9, y=0.168..2.832  (main body, 180×50×1 or variant)
  Block B: x=0..9, y=2.832..3.0    (top strip, 180×5×1)
  Block C: x=0..9, y=0..0.168      (bottom strip, 180×5×1)
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

# Geometry constants
L = 9.0        # room length
H = 3.0        # room height
h = 0.168      # slot width (inlet & outlet)
W = 0.1        # z-depth (2D)
Y_MID_LO = h          # y=0.168  (bottom strip top / main body bottom)
Y_MID_HI = H - h      # y=2.832  (top strip bottom / main body top)

# Physics
U_INLET = 0.455   # m/s  (Re ≈ 5000 based on slot height h=0.168m, ν=1.5e-5)
NU = 1.5e-5        # kinematic viscosity of air [m²/s]

CASES = {
    # --- Correct cases ---
    "S10_correct_turb": {
        "Re": 5000, "U": U_INLET, "turb": True,
        "nx": 180, "ny_main": 50, "ny_strip": 5,
        "endTime": 10000,
        "description": "Ventilated room, Re=5000, k-ω SST, 180×60 mesh",
    },
    "S10_correct_fine": {
        "Re": 5000, "U": U_INLET, "turb": True,
        "nx": 270, "ny_main": 75, "ny_strip": 7,
        "endTime": 10000,
        "description": "Ventilated room, Re=5000, k-ω SST, fine 270×89 mesh",
    },
    # --- Error cases ---
    "S10_E1_underconverged": {
        "Re": 5000, "U": U_INLET, "turb": True,
        "nx": 180, "ny_main": 50, "ny_strip": 5,
        "endTime": 30,
        "description": "Under-converged: stopped at 30 iterations only",
    },
    "S10_E2_bc_swap": {
        "Re": 5000, "U": U_INLET, "turb": True,
        "nx": 180, "ny_main": 50, "ny_strip": 5,
        "endTime": 10000, "swap_bc": True,
        "description": "Boundary condition error: inlet/outlet swapped",
    },
    "S10_E4_wrong_turb_model": {
        "Re": 5000, "U": U_INLET, "turb": False,
        "nx": 180, "ny_main": 50, "ny_strip": 5,
        "endTime": 10000,
        "description": "Wrong physics: Re=5000 but laminar solver (no turbulence model)",
    },
    "S10_E5_coarse_mesh": {
        "Re": 5000, "U": U_INLET, "turb": True,
        "nx": 30, "ny_main": 8, "ny_strip": 2,
        "endTime": 10000,
        "description": "Extremely coarse mesh: 30×12 only, no boundary layer resolution",
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
    """4-block layout:
       Block 0 (main body):   x=0..9, y=Y_MID_LO..Y_MID_HI
       Block 1 (top strip):   x=0..9, y=Y_MID_HI..H
       Block 2 (bottom strip): x=0..9, y=0..Y_MID_LO

    Vertices (24 total, 8 per block sharing interior faces):
      Block 0:  v0-v7
      Block 1:  v8-v15  (shares y=Y_MID_HI face with Block 0)
      Block 2:  v16-v23 (shares y=Y_MID_LO face with Block 0)

    Full numbering (front z=0 / back z=W):
      y=0:         v16(0,0,0)  v17(L,0,0)  v18(L,0,W)  v19(0,0,W)
      y=Y_MID_LO:  v0(0,lo,0)  v1(L,lo,0)  v2(L,lo,W)  v3(0,lo,W)
      y=Y_MID_HI:  v4(0,hi,0)  v5(L,hi,0)  v6(L,hi,W)  v7(0,hi,W)
      y=H:         v8(0,H,0)   v9(L,H,0)   v10(L,H,W)  v11(0,H,W)

    Blocks (hex uses counter-clockwise bottom then top):
      Block 0 (main): hex(v0 v1 v5 v4  v3 v2 v6 v7)   nx × ny_main × 1
      Block 1 (top):  hex(v4 v5 v9 v8  v7 v6 v10 v11) nx × ny_strip × 1
      Block 2 (bot):  hex(v16 v17 v1 v0  v19 v18 v2 v3) nx × ny_strip × 1

    Patches:
      inlet:       left face of Block 1 (x=0, y=Y_MID_HI..H)   → face (v8 v4 v7 v11)
      outlet:      right face of Block 2 (x=L, y=0..Y_MID_LO)  → face (v17 v1 v2 v18)  -- wait need CW
      upperWall:   top of Block 1 (y=H) + right of Block 1 (x=L)
      lowerWall:   bottom of Block 2 (y=0) + left of Block 2 (x=0)
      leftWall:    left of Block 0 (x=0, y=lo..hi) + left of Block 2 (x=0, y=0..lo)
      rightWall:   right of Block 0 (x=L, y=lo..hi) + right of Block 1 (x=L, y=hi..H)
      frontAndBack: all z faces
    """
    nx = cfg["nx"]
    ny_main = cfg["ny_main"]
    ny_strip = cfg["ny_strip"]
    lo = Y_MID_LO
    hi = Y_MID_HI

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
        // y = {lo} (bottom of main body / top of bottom strip)
        ( 0   {lo}  0 )    // 0
        ( {L}  {lo}  0 )    // 1
        ( {L}  {lo}  {W} )  // 2
        ( 0   {lo}  {W} )  // 3

        // y = {hi} (top of main body / bottom of top strip)
        ( 0   {hi}  0 )    // 4
        ( {L}  {hi}  0 )    // 5
        ( {L}  {hi}  {W} )  // 6
        ( 0   {hi}  {W} )  // 7

        // y = {H} (room ceiling)
        ( 0   {H}   0 )    // 8
        ( {L}  {H}   0 )    // 9
        ( {L}  {H}   {W} )  // 10
        ( 0   {H}   {W} )  // 11

        // y = 0 (room floor)
        ( 0   0    0 )    // 12
        ( {L}  0    0 )    // 13
        ( {L}  0    {W} )  // 14
        ( 0   0    {W} )  // 15
    );

    blocks
    (
        // Block 0: main body  x=0..{L}, y={lo}..{hi}
        hex (0 1 5 4  3 2 6 7)    ({nx} {ny_main} 1) simpleGrading (1 1 1)
        // Block 1: top strip  x=0..{L}, y={hi}..{H}
        hex (4 5 9 8  7 6 10 11)  ({nx} {ny_strip} 1) simpleGrading (1 1 1)
        // Block 2: bot strip  x=0..{L}, y=0..{lo}
        hex (12 13 1 0  15 14 2 3) ({nx} {ny_strip} 1) simpleGrading (1 1 1)
    );

    edges ( );

    boundary
    (
        inlet
        {{
            type patch;
            faces
            (
                (8 4 7 11)
            );
        }}
        outlet
        {{
            type patch;
            faces
            (
                (13 14 2 1)
            );
        }}
        upperWall
        {{
            type wall;
            faces
            (
                (8 9 10 11)
            );
        }}
        lowerWall
        {{
            type wall;
            faces
            (
                (12 13 14 15)
            );
        }}
        leftWall
        {{
            type wall;
            faces
            (
                (0 4 7 3)
                (12 0 3 15)
            );
        }}
        rightWall
        {{
            type wall;
            faces
            (
                (1 5 6 2)
                (5 9 10 6)
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                (0 1 5 4)
                (4 5 9 8)
                (12 13 1 0)
                (3 7 6 2)
                (7 11 10 6)
                (15 3 2 14)
            );
        }}
    );

    mergePatchPairs ( );
    """)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)
    (case_dir / "system" / "blockMeshDict").write_text(content)


def write_control_dict(case_dir: Path, cfg: dict) -> None:
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
    endTime         {cfg["endTime"]};
    deltaT          1;
    writeControl    timeStep;
    writeInterval   {cfg["endTime"]};
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
    U_in = cfg["U"]
    swap = cfg.get("swap_bc", False)

    # For swapped BC: velocity at outlet patch, pressure at inlet patch
    vel_patch = "outlet" if swap else "inlet"
    pres_patch = "inlet" if swap else "outlet"

    # U enters room horizontally from left side of top strip
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
        {vel_patch}
        {{
            type            fixedValue;
            value           uniform ({U_in} 0 0);
        }}
        {pres_patch}
        {{
            type            zeroGradient;
        }}
        "(upperWall|lowerWall|leftWall|rightWall)"
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
        {vel_patch}
        {{
            type            zeroGradient;
        }}
        {pres_patch}
        {{
            type            fixedValue;
            value           uniform 0;
        }}
        "(upperWall|lowerWall|leftWall|rightWall)"
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

    if cfg["turb"]:
        # Turbulent initial conditions: 5% intensity, length scale ~ h/10
        I = 0.05
        k_val = 1.5 * (I * U_in) ** 2
        L_t = 0.07 * h  # mixing-length estimate
        omega_val = (k_val ** 0.5) / (0.09 ** 0.25 * L_t)

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
            {vel_patch}
            {{
                type            fixedValue;
                value           uniform {k_val:.6e};
            }}
            {pres_patch}
            {{
                type            zeroGradient;
            }}
            "(upperWall|lowerWall|leftWall|rightWall)"
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
            {vel_patch}
            {{
                type            fixedValue;
                value           uniform {omega_val:.6e};
            }}
            {pres_patch}
            {{
                type            zeroGradient;
            }}
            "(upperWall|lowerWall|leftWall|rightWall)"
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

        nut_content = textwrap.dedent(f"""\
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            object      nut;
        }}

        dimensions      [0 2 -1 0 0 0 0];
        internalField   uniform 0;

        boundaryField
        {{
            {vel_patch}
            {{
                type            calculated;
                value           uniform 0;
            }}
            {pres_patch}
            {{
                type            calculated;
                value           uniform 0;
            }}
            "(upperWall|lowerWall|leftWall|rightWall)"
            {{
                type            nutkWallFunction;
                value           uniform 0;
            }}
            frontAndBack
            {{
                type            empty;
            }}
        }}
        """)

        (zero / "k").write_text(k_content)
        (zero / "omega").write_text(omega_content)
        (zero / "nut").write_text(nut_content)


def write_transport_properties(case_dir: Path, cfg: dict) -> None:
    # ν fixed at air value; Re and slot width h fix U_inlet
    content = textwrap.dedent(f"""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      transportProperties;
    }}

    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] {NU:.6e};
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

    meta = {
        "case_name": name,
        "scenario": "S10",
        "scenario_name": "ventilated_room",
        "Re": cfg["Re"],
        "U_inlet": cfg["U"],
        "nu": NU,
        "slot_height_m": h,
        "room_length_m": L,
        "room_height_m": H,
        "turbulence": cfg["turb"],
        "mesh_nx": cfg["nx"],
        "mesh_ny_main": cfg["ny_main"],
        "mesh_ny_strip": cfg["ny_strip"],
        "mesh_ny_total": cfg["ny_main"] + 2 * cfg["ny_strip"],
        "endTime": cfg["endTime"],
        "description": cfg["description"],
        "is_correct": not any(k in name for k in ["_E1_", "_E2_", "_E4_", "_E5_"]),
        "error_type": None,
    }
    if "_E1_" in name:
        meta["error_type"] = "E1_underconverged"
    elif "_E2_" in name:
        meta["error_type"] = "E2_bc_swap"
    elif "_E4_" in name:
        meta["error_type"] = "E4_wrong_turb_model"
    elif "_E5_" in name:
        meta["error_type"] = "E5_coarse_mesh"

    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name} (Re={cfg['Re']}, turb={cfg['turb']}, nx={cfg['nx']}, endTime={cfg['endTime']})")
    return case_dir


def run_case(name: str, case_dir: Path) -> bool:
    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"{'='*60}")

    # blockMesh
    r = docker_exec(case_dir, "blockMesh > log.blockMesh 2>&1")
    if r.returncode != 0:
        print(f"  blockMesh FAILED for {name}")
        return False

    # checkMesh
    r = docker_exec(case_dir, "checkMesh > log.checkMesh 2>&1", timeout=120)

    # simpleFoam
    r = docker_exec(case_dir, "simpleFoam > log.simpleFoam 2>&1", timeout=900)
    if r.returncode != 0:
        print(f"  simpleFoam returned {r.returncode} for {name} (may be OK for error cases)")

    # foamToVTK (latest time only)
    r = docker_exec(case_dir, "foamToVTK -latestTime > log.foamToVTK 2>&1", timeout=600)
    if r.returncode != 0:
        print(f"  foamToVTK FAILED for {name}")
        return False

    print(f"  Completed: {name}")
    return True


def main():
    CASE_BASE.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("CFD Visual QA — S10: Ventilated Room (Nielsen 1990 style)")
    print(f"  Room: {L}m × {H}m, slot h={h}m, U_inlet={U_INLET} m/s, Re≈5000")
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
