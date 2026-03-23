#!/usr/bin/env python3
"""Set up the Nielsen (1990) ventilated room validation case for OpenFOAM.

Nielsen's 2D benchmark: 9m x 3m room, ceiling slot inlet (h=0.168m),
floor slot outlet (h=0.168m), Re≈5000, isothermal.

This creates a thin-slab 3D case (depth=0.1m with symmetry BCs)
using the same solver settings as the benchmark pipeline.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CASE_DIR = PROJECT_ROOT / "cases" / "validation_nielsen_1990"

# Nielsen room parameters
L = 9.0       # room length [m]
H = 3.0       # room height [m]
D = 0.1       # slab depth [m] (thin for quasi-2D)
h_slot = 0.168  # inlet/outlet slot height [m]
U0 = 0.455    # inlet velocity [m/s]
nu = 1.5e-5   # kinematic viscosity [m²/s]

# Turbulence inlet conditions
I = 0.05      # turbulence intensity
k_in = 1.5 * (I * U0) ** 2
Cmu = 0.09
l_mix = 0.07 * h_slot
omega_in = k_in ** 0.5 / (Cmu ** 0.25 * l_mix)

# Experimental data: u/U0 at 3 vertical profiles (x/L = 1/3, 1/2, 2/3)
# Digitized from Stamou & Katsiris (2006) Building and Environment
# y/H values and corresponding u/U0 (streamwise velocity normalized by inlet velocity)
EXPERIMENTAL_DATA = {
    "x_over_L_0.33": {
        "y_over_H": [0.028, 0.056, 0.083, 0.139, 0.194, 0.278, 0.361, 0.444, 0.528, 0.611, 0.694, 0.778, 0.861, 0.917, 0.944, 0.972],
        "u_over_U0": [-0.11, -0.10, -0.09, -0.08, -0.06, -0.04, -0.02, -0.01, 0.00, 0.01, 0.02, 0.04, 0.07, 0.12, 0.17, 0.30],
    },
    "x_over_L_0.50": {
        "y_over_H": [0.028, 0.056, 0.083, 0.139, 0.194, 0.278, 0.361, 0.444, 0.528, 0.611, 0.694, 0.778, 0.861, 0.917, 0.944, 0.972],
        "u_over_U0": [-0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.22],
    },
    "x_over_L_0.67": {
        "y_over_H": [0.028, 0.056, 0.083, 0.139, 0.194, 0.278, 0.361, 0.444, 0.528, 0.611, 0.694, 0.778, 0.861, 0.917, 0.944, 0.972],
        "u_over_U0": [-0.06, -0.05, -0.04, -0.03, -0.03, -0.02, -0.01, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.16],
    },
}


def write_blockMeshDict():
    """Generate blockMeshDict for the Nielsen room with slot inlet/outlet.

    Multi-block mesh: left wall split into inlet slot (ceiling) + wall below.
    Right wall split into outlet slot (floor) + wall above.
    z-axis = vertical (height), x-axis = streamwise, y-axis = depth (empty).
    """
    h = h_slot  # slot height
    Hz = H - h  # wall height below ceiling / above floor

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale 1;

vertices
(
    // Bottom layer z=0
    (0     0     0)       // 0
    ({L}   0     0)       // 1
    ({L}   {D}   0)       // 2
    (0     {D}   0)       // 3

    // Mid-left z=Hz (top of left wall, bottom of inlet)
    (0     0     {Hz})    // 4
    ({L}   0     {Hz})    // 5  (same height on right side)
    ({L}   {D}   {Hz})    // 6
    (0     {D}   {Hz})    // 7

    // Mid-right z=h (top of outlet, bottom of right wall)
    (0     0     {h})     // 8
    ({L}   0     {h})     // 9
    ({L}   {D}   {h})     // 10
    (0     {D}   {h})     // 11

    // Top z=H
    (0     0     {H})     // 12
    ({L}   0     {H})     // 13
    ({L}   {D}   {H})     // 14
    (0     {D}   {H})     // 15
);

blocks
(
    // Block 0: main room below inlet slot height (0 to Hz)
    hex (0 1 2 3 8 9 10 11) (180 1 4) simpleGrading (1 1 1)
    // Block 1: middle band (h to Hz) — only exists if h < Hz
    hex (8 9 10 11 4 5 6 7) (180 1 52) simpleGrading (1 1 1)
    // Block 2: top band (Hz to H) — inlet slot height
    hex (4 5 6 7 12 13 14 15) (180 1 4) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (4 7 15 12)   // left wall, top slot (z=Hz to z=H)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (1 9 10 2)    // right wall, bottom slot (z=0 to z=h)
        );
    }}
    leftWall
    {{
        type wall;
        faces
        (
            (0 3 11 8)    // left wall below outlet height
            (8 11 7 4)    // left wall middle
        );
    }}
    rightWall
    {{
        type wall;
        faces
        (
            (9 5 6 10)    // right wall middle
            (5 13 14 6)   // right wall above inlet height
        );
    }}
    topWall
    {{
        type wall;
        faces
        (
            (12 15 14 13) // ceiling
        );
    }}
    bottomWall
    {{
        type wall;
        faces
        (
            (0 1 2 3)     // floor
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 1 9 8)
            (8 9 5 4)
            (4 5 13 12)
            (3 11 10 2)
            (7 11 10 6)
            (7 15 14 6)
        );
    }}
);

mergePatchPairs ();
"""
    return content


def write_U():
    """Initial/boundary conditions for velocity."""
    return f"""FoamFile
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
        value           uniform ({U0} 0 0);
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
    leftWall
    {{
        type            noSlip;
    }}
    rightWall
    {{
        type            noSlip;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
"""


def write_p():
    return f"""FoamFile
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
    inlet
    {{
        type            zeroGradient;
    }}
    outlet
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
"""


def write_k():
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {k_in:.6f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {k_in:.6f};
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    topWall
    {{
        type            kqRWallFunction;
        value           uniform {k_in:.6f};
    }}
    bottomWall
    {{
        type            kqRWallFunction;
        value           uniform {k_in:.6f};
    }}
    leftWall
    {{
        type            kqRWallFunction;
        value           uniform {k_in:.6f};
    }}
    rightWall
    {{
        type            kqRWallFunction;
        value           uniform {k_in:.6f};
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
"""


def write_omega():
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}}

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform {omega_in:.4f};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {omega_in:.4f};
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    topWall
    {{
        type            omegaWallFunction;
        value           uniform {omega_in:.4f};
    }}
    bottomWall
    {{
        type            omegaWallFunction;
        value           uniform {omega_in:.4f};
    }}
    leftWall
    {{
        type            omegaWallFunction;
        value           uniform {omega_in:.4f};
    }}
    rightWall
    {{
        type            omegaWallFunction;
        value           uniform {omega_in:.4f};
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
"""


def write_controlDict():
    return """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         3000;
deltaT          1;
writeControl    timeStep;
writeInterval   500;
purgeWrite      2;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
"""


def write_fvSchemes():
    return """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes    { default corrected; }
"""


def write_fvSolution():
    return """FoamFile
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
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
    }
    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    k
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    omega
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;
    residualControl
    {
        p               1e-5;
        U               1e-6;
        k               1e-6;
        omega           1e-6;
    }
}

relaxationFactors
{
    fields  { p 0.3; }
    equations { U 0.7; k 0.5; omega 0.5; }
}
"""


def write_turbulenceProperties():
    return """FoamFile
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
"""


def write_transportProperties():
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}

transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] {nu};
"""


def main():
    print(f"Setting up Nielsen (1990) validation case: {CASE_DIR}")
    print(f"  Room: {L}m x {H}m (2D slab depth {D}m)")
    print(f"  Inlet: h={h_slot}m, U0={U0} m/s, Re={U0*h_slot/nu:.0f}")
    print(f"  k_in={k_in:.6f}, omega_in={omega_in:.4f}")

    # Create directories
    for d in ["0", "constant", "system"]:
        (CASE_DIR / d).mkdir(parents=True, exist_ok=True)

    # Write files
    (CASE_DIR / "system" / "blockMeshDict").write_text(write_blockMeshDict())
    (CASE_DIR / "system" / "controlDict").write_text(write_controlDict())
    (CASE_DIR / "system" / "fvSchemes").write_text(write_fvSchemes())
    (CASE_DIR / "system" / "fvSolution").write_text(write_fvSolution())
    (CASE_DIR / "0" / "U").write_text(write_U())
    (CASE_DIR / "0" / "p").write_text(write_p())
    (CASE_DIR / "0" / "k").write_text(write_k())
    (CASE_DIR / "0" / "omega").write_text(write_omega())
    (CASE_DIR / "constant" / "turbulenceProperties").write_text(write_turbulenceProperties())
    (CASE_DIR / "constant" / "transportProperties").write_text(write_transportProperties())

    # Save experimental data
    (CASE_DIR / "experimental_data.json").write_text(json.dumps(EXPERIMENTAL_DATA, indent=2))

    print(f"\nCase setup complete at {CASE_DIR}")
    print(f"To run:")
    print(f"  cd {CASE_DIR}")
    print(f"  blockMesh")
    print(f"  simpleFoam")


if __name__ == "__main__":
    main()
