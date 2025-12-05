#!/usr/bin/env python3
"""
Create buoyantSimpleFoam case for heatsink cooling analysis
Simulates air flow over a heated heatsink (base + fins)
"""

import os


def create_heatsink_flow_case(
    case_dir: str,
    # Geometry parameters (m)
    channel_length: float = 0.12,     # 120mm channel length
    channel_height: float = 0.04,     # 40mm channel height
    channel_width: float = 0.06,      # 60mm channel width (depth)
    base_length: float = 0.05,        # 50mm heatsink base length
    base_height: float = 0.005,       # 5mm base thickness
    fin_height: float = 0.020,        # 20mm fin height
    fin_thickness: float = 0.003,     # 3mm fin thickness
    num_fins: int = 5,                # Number of fins
    # Boundary conditions
    inlet_velocity: float = 0.1,      # m/s
    inlet_temp: float = 298.15,       # K (25°C)
    wall_temp: float = 353.15,        # K (80°C) - heatsink surface temp
    # Mesh resolution
    cells_x: int = 60,
    cells_y: int = 40,
    cells_z: int = 30
):
    """Create a buoyantSimpleFoam case for heatsink flow simulation."""

    os.makedirs(case_dir, exist_ok=True)

    # Calculate geometry positions
    hs_start_x = (channel_length - base_length) / 2
    hs_end_x = hs_start_x + base_length
    total_hs_height = base_height + fin_height

    # Create directories
    for d in ['system', 'constant', '0']:
        os.makedirs(os.path.join(case_dir, d), exist_ok=True)

    # ========================================
    # system/blockMeshDict - Simple channel with heatsink block
    # ========================================
    # Simplified: treat heatsink as a solid obstacle with fixed temperature walls
    block_mesh = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale   1;

// Channel: {channel_length*1000:.0f}mm x {channel_height*1000:.0f}mm x {channel_width*1000:.0f}mm
// Heatsink: located at center bottom, {base_length*1000:.0f}mm long, {total_hs_height*1000:.0f}mm tall

vertices
(
    // Bottom plane (y=0)
    (0 0 0)                              // 0
    ({hs_start_x} 0 0)                   // 1
    ({hs_end_x} 0 0)                     // 2
    ({channel_length} 0 0)               // 3
    (0 0 {channel_width})                // 4
    ({hs_start_x} 0 {channel_width})     // 5
    ({hs_end_x} 0 {channel_width})       // 6
    ({channel_length} 0 {channel_width}) // 7

    // Heatsink top plane (y=total_hs_height)
    ({hs_start_x} {total_hs_height} 0)             // 8
    ({hs_end_x} {total_hs_height} 0)               // 9
    ({hs_start_x} {total_hs_height} {channel_width}) // 10
    ({hs_end_x} {total_hs_height} {channel_width})   // 11

    // Channel top plane (y=channel_height)
    (0 {channel_height} 0)                         // 12
    ({hs_start_x} {channel_height} 0)              // 13
    ({hs_end_x} {channel_height} 0)                // 14
    ({channel_length} {channel_height} 0)          // 15
    (0 {channel_height} {channel_width})           // 16
    ({hs_start_x} {channel_height} {channel_width}) // 17
    ({hs_end_x} {channel_height} {channel_width})   // 18
    ({channel_length} {channel_height} {channel_width}) // 19
);

blocks
(
    // Block 0: Inlet section (before heatsink) - full height
    hex (0 1 5 4 12 13 17 16) ({int(cells_x * hs_start_x / channel_length)} {cells_y} {cells_z}) simpleGrading (1 1 1)

    // Block 1: Above heatsink
    hex (8 9 11 10 13 14 18 17) ({int(cells_x * base_length / channel_length)} {int(cells_y * (channel_height - total_hs_height) / channel_height)} {cells_z}) simpleGrading (1 1 1)

    // Block 2: Outlet section (after heatsink) - full height
    hex (2 3 7 6 14 15 19 18) ({int(cells_x * (channel_length - hs_end_x) / channel_length)} {cells_y} {cells_z}) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 4 16 12)
        );
    }}

    outlet
    {{
        type patch;
        faces
        (
            (3 15 19 7)
        );
    }}

    topWall
    {{
        type wall;
        faces
        (
            (12 13 17 16)
            (13 14 18 17)
            (14 15 19 18)
        );
    }}

    bottomWall
    {{
        type wall;
        faces
        (
            (0 1 5 4)
            (2 3 7 6)
        );
    }}

    heatsinkWall
    {{
        type wall;
        faces
        (
            // Front face of heatsink (facing inlet)
            (1 8 10 5)
            // Top of heatsink
            (8 9 11 10)
            // Back face of heatsink (facing outlet)
            (9 2 6 11)
        );
    }}

    frontAndBack
    {{
        type wall;
        faces
        (
            // Front (z=0)
            (0 12 13 1)
            (1 13 8 1)
            (8 13 14 9)
            (9 14 15 2)
            (2 15 3 2)
            // Back (z=channel_width)
            (4 5 17 16)
            (5 10 17 5)
            (10 11 18 17)
            (11 6 18 11)
            (6 7 19 18)
        );
    }}
);

mergePatchPairs ();
"""

    # Simplified mesh - just a channel with the heatsink as a wall BC
    simple_block_mesh = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale   1;

// Simple rectangular channel for flow simulation
// Heatsink represented as a heated wall at the bottom

vertices
(
    (0 0 0)                              // 0
    ({channel_length} 0 0)               // 1
    ({channel_length} {channel_height} 0) // 2
    (0 {channel_height} 0)               // 3
    (0 0 {channel_width})                // 4
    ({channel_length} 0 {channel_width}) // 5
    ({channel_length} {channel_height} {channel_width}) // 6
    (0 {channel_height} {channel_width}) // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({cells_x} {cells_y} {cells_z}) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}

    outlet
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}

    topWall
    {{
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }}

    heatedWall
    {{
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }}

    frontAndBack
    {{
        type wall;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);

mergePatchPairs ();
"""
    with open(os.path.join(case_dir, 'system', 'blockMeshDict'), 'w') as f:
        f.write(simple_block_mesh)

    # ========================================
    # system/controlDict
    # ========================================
    control_dict = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}

application     buoyantSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         2000;
deltaT          1;
writeControl    timeStep;
writeInterval   200;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    fieldAverage1
    {{
        type            fieldAverage;
        libs            ("libfieldFunctionObjects.so");
        enabled         true;
        writeControl    writeTime;
        fields
        (
            U
            {{
                mean        on;
                prime2Mean  off;
                base        time;
            }}
            T
            {{
                mean        on;
                prime2Mean  off;
                base        time;
            }}
        );
    }}

    outletAverage
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        enabled         true;
        writeControl    timeStep;
        writeInterval   10;
        log             true;
        writeFields     false;
        regionType      patch;
        name            outlet;
        operation       areaAverage;
        fields          (T U);
    }}
}}
"""
    with open(os.path.join(case_dir, 'system', 'controlDict'), 'w') as f:
        f.write(control_dict)

    # ========================================
    # system/fvSchemes
    # ========================================
    fv_schemes = """FoamFile
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
    div(phi,h)      bounded Gauss linearUpwind grad(h);
    div(phi,K)      bounded Gauss upwind;
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
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
"""
    with open(os.path.join(case_dir, 'system', 'fvSchemes'), 'w') as f:
        f.write(fv_schemes)

    # ========================================
    # system/fvSolution
    # ========================================
    fv_solution = """FoamFile
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
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|h|k|epsilon|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(U|h|k|epsilon|omega)Final"
    {
        $U;
        relTol          0;
    }
}

SIMPLE
{
    momentumPredictor   yes;
    nNonOrthogonalCorrectors 0;
    pRefCell            0;
    pRefValue           0;

    residualControl
    {
        p_rgh           1e-4;
        U               1e-4;
        h               1e-4;
    }
}

relaxationFactors
{
    fields
    {
        p_rgh           0.7;
        rho             1.0;
    }
    equations
    {
        U               0.3;
        h               0.5;
        k               0.5;
        epsilon         0.5;
        omega           0.5;
    }
}
"""
    with open(os.path.join(case_dir, 'system', 'fvSolution'), 'w') as f:
        f.write(fv_solution)

    # ========================================
    # constant/g
    # ========================================
    g_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}

dimensions      [0 1 -2 0 0 0 0];
value           (0 -9.81 0);
"""
    with open(os.path.join(case_dir, 'constant', 'g'), 'w') as f:
        f.write(g_file)

    # ========================================
    # constant/thermophysicalProperties (Air)
    # ========================================
    thermo_props = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       28.96;
    }
    thermodynamics
    {
        Cp              1005;
        Hf              0;
    }
    transport
    {
        mu              1.85e-5;
        Pr              0.71;
    }
}
"""
    with open(os.path.join(case_dir, 'constant', 'thermophysicalProperties'), 'w') as f:
        f.write(thermo_props)

    # ========================================
    # constant/turbulenceProperties
    # ========================================
    turb_props = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  laminar;
"""
    with open(os.path.join(case_dir, 'constant', 'turbulenceProperties'), 'w') as f:
        f.write(turb_props)

    # ========================================
    # 0/U
    # ========================================
    U_file = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({inlet_velocity} 0 0);

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform ({inlet_velocity} 0 0);
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    topWall
    {{
        type            noSlip;
    }}

    heatedWall
    {{
        type            noSlip;
    }}

    frontAndBack
    {{
        type            noSlip;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'U'), 'w') as f:
        f.write(U_file)

    # ========================================
    # 0/p_rgh
    # ========================================
    p_rgh_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}

dimensions      [1 -1 -2 0 0 0 0];

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

    heatedWall
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }

    frontAndBack
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }
}
"""
    with open(os.path.join(case_dir, '0', 'p_rgh'), 'w') as f:
        f.write(p_rgh_file)

    # ========================================
    # 0/p
    # ========================================
    p_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 101325;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           $internalField;
    }
}
"""
    with open(os.path.join(case_dir, '0', 'p'), 'w') as f:
        f.write(p_file)

    # ========================================
    # 0/T
    # ========================================
    T_file = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {inlet_temp};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {inlet_temp};
    }}

    outlet
    {{
        type            zeroGradient;
    }}

    topWall
    {{
        type            zeroGradient;
    }}

    heatedWall
    {{
        type            fixedValue;
        value           uniform {wall_temp};
    }}

    frontAndBack
    {{
        type            zeroGradient;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'T'), 'w') as f:
        f.write(T_file)

    # ========================================
    # 0/alphat
    # ========================================
    alphat_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alphat;
}

dimensions      [1 -1 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           uniform 0;
    }
}
"""
    with open(os.path.join(case_dir, '0', 'alphat'), 'w') as f:
        f.write(alphat_file)

    # ========================================
    # Allrun script
    # ========================================
    allrun = f"""#!/bin/bash
cd "${{0%/*}}" || exit

# Set OpenFOAM environment
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

echo "=========================================="
echo "Heatsink Flow Simulation"
echo "Inlet velocity: {inlet_velocity} m/s"
echo "Inlet temperature: {inlet_temp - 273.15:.1f} degC"
echo "Wall temperature: {wall_temp - 273.15:.1f} degC"
echo "=========================================="

echo "Creating mesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: blockMesh failed"
    cat log.blockMesh
    exit 1
fi

echo "Checking mesh..."
checkMesh > log.checkMesh 2>&1

echo "Running buoyantSimpleFoam..."
buoyantSimpleFoam > log.buoyantSimpleFoam 2>&1

echo "Converting to VTK..."
foamToVTK > log.foamToVTK 2>&1

echo ""
echo "Simulation complete!"
echo "Check log.buoyantSimpleFoam for solver output"
echo "Results in VTK/ directory"

# Extract outlet temperature
if [ -d postProcessing/outletAverage ]; then
    echo ""
    echo "Outlet average temperature history:"
    tail -5 postProcessing/outletAverage/*/surfaceFieldValue.dat
fi
"""
    with open(os.path.join(case_dir, 'Allrun'), 'w') as f:
        f.write(allrun)
    os.chmod(os.path.join(case_dir, 'Allrun'), 0o755)

    # ========================================
    # Allclean script
    # ========================================
    allclean = """#!/bin/bash
cd "${0%/*}" || exit
rm -rf 0.* [1-9]* constant/polyMesh VTK postProcessing log.*
"""
    with open(os.path.join(case_dir, 'Allclean'), 'w') as f:
        f.write(allclean)
    os.chmod(os.path.join(case_dir, 'Allclean'), 0o755)

    print(f"Heatsink flow case created in {case_dir}")
    print(f"  Channel: {channel_length*1000:.0f} x {channel_height*1000:.0f} x {channel_width*1000:.0f} mm")
    print(f"  Inlet velocity: {inlet_velocity} m/s")
    print(f"  Inlet temperature: {inlet_temp - 273.15:.1f} °C")
    print(f"  Wall temperature: {wall_temp - 273.15:.1f} °C")
    print(f"\nTo run: cd {case_dir} && ./Allrun")

    return case_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create heatsink flow case')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--inlet-velocity', type=float, default=0.1,
                        help='Inlet velocity (m/s)')
    parser.add_argument('--inlet-temp', type=float, default=25.0,
                        help='Inlet temperature (°C)')
    parser.add_argument('--wall-temp', type=float, default=80.0,
                        help='Heated wall temperature (°C)')

    args = parser.parse_args()

    create_heatsink_flow_case(
        args.case_dir,
        inlet_velocity=args.inlet_velocity,
        inlet_temp=args.inlet_temp + 273.15,
        wall_temp=args.wall_temp + 273.15
    )
