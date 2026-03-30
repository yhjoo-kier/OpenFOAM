#!/usr/bin/env python3
"""Setup backward-facing step (BFS) OpenFOAM cases for CFD Visual QA pilot."""
from __future__ import annotations

import json
import shutil
import subprocess
import textwrap
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]  # OpenFOAM/
CASE_BASE = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark" / "cases"
DOCKER_IMAGE = "openfoam-pipeline-local:latest"

# Geometry: step height h=1m
H = 1.0          # step height
L_UP = 5.0       # upstream length
L_DOWN = 30.0    # downstream length
W = 0.1          # z-width (2D)

CASES = {
    # --- Correct cases ---
    "S4_correct_lam": {
        "Re": 200, "U": 1.0, "turb": False,
        "nx_up": 30, "nx_down": 180, "ny_upper": 25, "ny_lower": 25,
        "grading_y": 4.0, "endTime": 5000, "description": "Laminar BFS, Re=200, well-converged",
    },
    "S4_correct_turb": {
        "Re": 36000, "U": 1.0, "turb": True,
        "nx_up": 40, "nx_down": 240, "ny_upper": 40, "ny_lower": 40,
        "grading_y": 6.0, "endTime": 5000, "description": "Turbulent BFS, Re=36000, k-ω SST",
    },
    # --- Error cases ---
    "S4_E1_underconverged": {
        "Re": 36000, "U": 1.0, "turb": True,
        "nx_up": 40, "nx_down": 240, "ny_upper": 40, "ny_lower": 40,
        "grading_y": 6.0, "endTime": 50,  # early stop!
        "description": "Under-converged: stopped at 50 iterations",
    },
    "S4_E2_bc_swap": {
        "Re": 36000, "U": 1.0, "turb": True,
        "nx_up": 40, "nx_down": 240, "ny_upper": 40, "ny_lower": 40,
        "grading_y": 6.0, "endTime": 5000, "swap_bc": True,
        "description": "Boundary condition error: inlet/outlet swapped",
    },
    "S4_E5_coarse_mesh": {
        "Re": 36000, "U": 1.0, "turb": True,
        "nx_up": 8, "nx_down": 40, "ny_upper": 6, "ny_lower": 6,
        "grading_y": 1.0, "endTime": 5000,
        "description": "Extremely coarse mesh: no boundary layer resolution",
    },
    "S4_E6_modified": {
        "Re": 36000, "U": 1.0, "turb": True,
        "nx_up": 40, "nx_down": 240, "ny_upper": 40, "ny_lower": 40,
        "grading_y": 6.0, "endTime": 5000, "modify_result": True,
        "description": "Artificial modification: recirculation zone removed post-solve",
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
    nu = cfg["U"] * H / cfg["Re"]
    nx_up = cfg["nx_up"]
    nx_down = cfg["nx_down"]
    ny_u = cfg["ny_upper"]
    ny_l = cfg["ny_lower"]
    gy = cfg["grading_y"]
    inv_gy = round(1.0 / gy, 6)

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
        // Block 0: upstream channel (x: -{L_UP} to 0, y: {H} to {2*H})
        ( -{L_UP}  {H}    0 )    // 0
        (  0       {H}    0 )    // 1
        (  0       {2*H}  0 )    // 2
        ( -{L_UP}  {2*H}  0 )    // 3
        ( -{L_UP}  {H}    {W} )  // 4
        (  0       {H}    {W} )  // 5
        (  0       {2*H}  {W} )  // 6
        ( -{L_UP}  {2*H}  {W} )  // 7

        // Block 1: downstream upper (x: 0 to {L_DOWN}, y: {H} to {2*H})
        ( {L_DOWN}  {H}    0 )   // 8
        ( {L_DOWN}  {2*H}  0 )   // 9
        ( {L_DOWN}  {H}    {W} ) // 10
        ( {L_DOWN}  {2*H}  {W} ) // 11

        // Block 2: downstream lower (x: 0 to {L_DOWN}, y: 0 to {H})
        (  0       0      0 )    // 12
        ( {L_DOWN} 0      0 )    // 13
        ( {L_DOWN} 0      {W} )  // 14
        (  0       0      {W} )  // 15
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7)       ({nx_up}  {ny_u} 1) simpleGrading (1 {gy} 1)
        hex (1 8 9 2 5 10 11 6)     ({nx_down} {ny_u} 1) simpleGrading (1 {gy} 1)
        hex (12 13 8 1 15 14 10 5)  ({nx_down} {ny_l} 1) simpleGrading (1 {inv_gy} 1)
    );

    edges ( );

    boundary
    (
        inlet
        {{
            type patch;
            faces ( (0 3 7 4) );
        }}
        outlet
        {{
            type patch;
            faces
            (
                (8 9 11 10)
                (13 8 10 14)
            );
        }}
        upperWall
        {{
            type wall;
            faces
            (
                (3 2 6 7)
                (2 9 11 6)
            );
        }}
        lowerWall
        {{
            type wall;
            faces ( (12 13 14 15) );
        }}
        stepWall
        {{
            type wall;
            faces
            (
                (0 1 5 4)
                (1 12 15 5)
            );
        }}
        frontAndBack
        {{
            type empty;
            faces
            (
                (0 1 2 3)
                (1 8 9 2)
                (12 13 8 1)
                (4 5 6 7)
                (5 10 11 6)
                (15 14 10 5)
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
    nu = cfg["U"] * H / cfg["Re"]
    U_inlet = cfg["U"]
    swap = cfg.get("swap_bc", False)

    inlet_patch = "outlet" if swap else "inlet"
    outlet_patch = "inlet" if swap else "outlet"

    # U
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
        {inlet_patch}
        {{
            type            fixedValue;
            value           uniform ({U_inlet} 0 0);
        }}
        {outlet_patch}
        {{
            type            zeroGradient;
        }}
        "(upperWall|lowerWall|stepWall)"
        {{
            type            noSlip;
        }}
        frontAndBack
        {{
            type            empty;
        }}
    }}
    """)

    # p
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
        {inlet_patch}
        {{
            type            zeroGradient;
        }}
        {outlet_patch}
        {{
            type            fixedValue;
            value           uniform 0;
        }}
        "(upperWall|lowerWall|stepWall)"
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
        # k
        k_val = 0.5 * (0.05 * U_inlet) ** 2  # 5% turbulence intensity
        omega_val = k_val ** 0.5 / (0.09 ** 0.25 * 0.1 * H)

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
            {inlet_patch}
            {{
                type            fixedValue;
                value           uniform {k_val:.6e};
            }}
            {outlet_patch}
            {{
                type            zeroGradient;
            }}
            "(upperWall|lowerWall|stepWall)"
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
            {inlet_patch}
            {{
                type            fixedValue;
                value           uniform {omega_val:.6e};
            }}
            {outlet_patch}
            {{
                type            zeroGradient;
            }}
            "(upperWall|lowerWall|stepWall)"
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
            {inlet_patch}
            {{
                type            calculated;
                value           uniform 0;
            }}
            {outlet_patch}
            {{
                type            calculated;
                value           uniform 0;
            }}
            "(upperWall|lowerWall|stepWall)"
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
    nu = cfg["U"] * H / cfg["Re"]
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

    # Write case metadata
    meta = {
        "case_name": name,
        "scenario": "S4",
        "scenario_name": "backward_facing_step",
        "Re": cfg["Re"],
        "U_inlet": cfg["U"],
        "turbulence": cfg["turb"],
        "description": cfg["description"],
        "is_correct": not any(k in name for k in ["_E1_", "_E2_", "_E5_", "_E6_"]),
        "error_type": None,
    }
    if "_E1_" in name:
        meta["error_type"] = "E1_underconverged"
    elif "_E2_" in name:
        meta["error_type"] = "E2_bc_swap"
    elif "_E5_" in name:
        meta["error_type"] = "E5_coarse_mesh"
    elif "_E6_" in name:
        meta["error_type"] = "E6_modified"

    (case_dir / "case_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  Setup: {name} (Re={cfg['Re']}, turb={cfg['turb']})")
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

    # Check mesh
    r = docker_exec(case_dir, "checkMesh > log.checkMesh 2>&1", timeout=120)

    # simpleFoam
    r = docker_exec(case_dir, "simpleFoam > log.simpleFoam 2>&1", timeout=900)
    if r.returncode != 0:
        print(f"  simpleFoam returned {r.returncode} for {name} (may be OK for error cases)")

    # foamToVTK
    r = docker_exec(case_dir, "foamToVTK -latestTime > log.foamToVTK 2>&1", timeout=600)
    if r.returncode != 0:
        print(f"  foamToVTK FAILED for {name}")
        return False

    print(f"  Completed: {name}")
    return True


def main():
    CASE_BASE.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("CFD Visual QA — Pilot: Backward-Facing Step")
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
