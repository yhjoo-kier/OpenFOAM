#!/usr/bin/env python3
"""
Create chtMultiRegionFoam case for heatsink cooling analysis
- Solid region: Aluminum heatsink (base + fin)
- Fluid region: Air flow around the heatsink
"""

import os
import shutil


def create_cht_heatsink_case(
    case_dir: str,
    # Geometry parameters (m)
    channel_length: float = 0.1,      # 100mm channel length
    channel_height: float = 0.03,     # 30mm channel height
    channel_width: float = 0.05,      # 50mm channel width (depth)
    base_length: float = 0.05,        # 50mm heatsink base length
    base_height: float = 0.005,       # 5mm base thickness
    fin_height: float = 0.015,        # 15mm fin height
    fin_thickness: float = 0.003,     # 3mm fin thickness
    # Boundary conditions
    inlet_velocity: float = 0.1,      # m/s
    inlet_temp: float = 298.15,       # K (25°C)
    base_temp: float = 353.15,        # K (80°C)
    # Mesh resolution
    cells_x: int = 50,
    cells_y: int = 25,
    cells_z: int = 25
):
    """Create a complete CHT case for heatsink cooling."""

    os.makedirs(case_dir, exist_ok=True)

    # Calculate geometry
    heatsink_start_x = (channel_length - base_length) / 2
    heatsink_end_x = heatsink_start_x + base_length
    total_solid_height = base_height + fin_height

    # Create directory structure
    for d in ['system', 'system/fluid', 'system/solid',
              'constant', 'constant/fluid', 'constant/solid',
              '0', '0/fluid', '0/solid']:
        os.makedirs(os.path.join(case_dir, d), exist_ok=True)

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

application     chtMultiRegionFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         500;
deltaT          1;
writeControl    timeStep;
writeInterval   50;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    fieldAverage
    {{
        type            fieldAverage;
        libs            ("libfieldFunctionObjects.so");
        enabled         true;
        writeControl    writeTime;
        region          fluid;
        fields
        (
            T
            {{
                mean        on;
                prime2Mean  off;
                base        time;
            }}
            U
            {{
                mean        on;
                prime2Mean  off;
                base        time;
            }}
        );
    }}
}}
"""
    with open(os.path.join(case_dir, 'system', 'controlDict'), 'w') as f:
        f.write(control_dict)

    # ========================================
    # system/fvSchemes (global)
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
}

laplacianSchemes
{
    default         none;
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
    # system/fvSolution (global)
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
}
"""
    with open(os.path.join(case_dir, 'system', 'fvSolution'), 'w') as f:
        f.write(fv_solution)

    # ========================================
    # system/regionProperties
    # ========================================
    region_props = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      regionProperties;
}

regions
(
    fluid       (fluid)
    solid       (solid)
);
"""
    with open(os.path.join(case_dir, 'system', 'regionProperties'), 'w') as f:
        f.write(region_props)

    # ========================================
    # system/fluid/fvSchemes
    # ========================================
    fluid_fv_schemes = """FoamFile
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
    grad(p_rgh)     Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,K)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss linearUpwind grad(h);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
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
    with open(os.path.join(case_dir, 'system', 'fluid', 'fvSchemes'), 'w') as f:
        f.write(fluid_fv_schemes)

    # ========================================
    # system/fluid/fvSolution
    # ========================================
    fluid_fv_solution = """FoamFile
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

    "(U|h|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(U|h|k|epsilon)Final"
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
    }
    equations
    {
        U               0.3;
        h               0.5;
        k               0.5;
        epsilon         0.5;
    }
}
"""
    with open(os.path.join(case_dir, 'system', 'fluid', 'fvSolution'), 'w') as f:
        f.write(fluid_fv_solution)

    # ========================================
    # system/solid/fvSchemes
    # ========================================
    solid_fv_schemes = """FoamFile
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
}

laplacianSchemes
{
    default         none;
    laplacian(alpha,h) Gauss linear corrected;
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
    with open(os.path.join(case_dir, 'system', 'solid', 'fvSchemes'), 'w') as f:
        f.write(solid_fv_schemes)

    # ========================================
    # system/solid/fvSolution
    # ========================================
    solid_fv_solution = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    h
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0.01;
    }

    hFinal
    {
        $h;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
}
"""
    with open(os.path.join(case_dir, 'system', 'solid', 'fvSolution'), 'w') as f:
        f.write(solid_fv_solution)

    # ========================================
    # constant/fluid/thermophysicalProperties (Air)
    # ========================================
    # Air properties at ~50°C
    fluid_thermo = """FoamFile
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
        molWeight       28.96;  // Air molecular weight
    }
    thermodynamics
    {
        Cp              1005;   // J/(kg·K)
        Hf              0;
    }
    transport
    {
        mu              1.85e-5;  // Dynamic viscosity (Pa·s)
        Pr              0.71;     // Prandtl number
    }
}
"""
    with open(os.path.join(case_dir, 'constant', 'fluid', 'thermophysicalProperties'), 'w') as f:
        f.write(fluid_thermo)

    # ========================================
    # constant/fluid/turbulenceProperties
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
    with open(os.path.join(case_dir, 'constant', 'fluid', 'turbulenceProperties'), 'w') as f:
        f.write(turb_props)

    # ========================================
    # constant/fluid/g
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
    with open(os.path.join(case_dir, 'constant', 'fluid', 'g'), 'w') as f:
        f.write(g_file)

    # ========================================
    # constant/solid/thermophysicalProperties (Aluminum)
    # ========================================
    solid_thermo = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}

thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       26.98;  // Aluminum molecular weight
    }
    transport
    {
        kappa           205;    // Thermal conductivity (W/(m·K))
    }
    thermodynamics
    {
        Cp              900;    // Specific heat (J/(kg·K))
        Hf              0;
    }
    equationOfState
    {
        rho             2700;   // Density (kg/m³)
    }
}
"""
    with open(os.path.join(case_dir, 'constant', 'solid', 'thermophysicalProperties'), 'w') as f:
        f.write(solid_thermo)

    # ========================================
    # Create blockMeshDict for multi-region mesh
    # Using a simplified geometry: channel with solid block at bottom
    # ========================================
    # Geometry: Channel with heatsink (base + single fin) at bottom center
    # Fluid flows in x-direction, heatsink in center

    # Simplified: Just base plate (no fin for now) to test CHT
    block_mesh = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale   1;

// Channel: {channel_length*1000:.1f}mm x {channel_height*1000:.1f}mm x {channel_width*1000:.1f}mm
// Heatsink base: {base_length*1000:.1f}mm x {base_height*1000:.1f}mm (at bottom center)

// Calculated positions
// heatsink_start_x = {heatsink_start_x}
// heatsink_end_x = {heatsink_end_x}

vertices
(
    // Bottom of channel (z=0) - y=0 plane
    (0 0 0)                                          // 0
    ({heatsink_start_x} 0 0)                         // 1
    ({heatsink_end_x} 0 0)                           // 2
    ({channel_length} 0 0)                           // 3
    (0 0 {channel_width})                            // 4
    ({heatsink_start_x} 0 {channel_width})           // 5
    ({heatsink_end_x} 0 {channel_width})             // 6
    ({channel_length} 0 {channel_width})             // 7

    // Top of solid / bottom of fluid (y=base_height)
    (0 {base_height} 0)                              // 8
    ({heatsink_start_x} {base_height} 0)             // 9
    ({heatsink_end_x} {base_height} 0)               // 10
    ({channel_length} {base_height} 0)               // 11
    (0 {base_height} {channel_width})                // 12
    ({heatsink_start_x} {base_height} {channel_width}) // 13
    ({heatsink_end_x} {base_height} {channel_width})   // 14
    ({channel_length} {base_height} {channel_width})   // 15

    // Top of channel (y=channel_height)
    (0 {channel_height} 0)                           // 16
    ({heatsink_start_x} {channel_height} 0)          // 17
    ({heatsink_end_x} {channel_height} 0)            // 18
    ({channel_length} {channel_height} 0)            // 19
    (0 {channel_height} {channel_width})             // 20
    ({heatsink_start_x} {channel_height} {channel_width}) // 21
    ({heatsink_end_x} {channel_height} {channel_width})   // 22
    ({channel_length} {channel_height} {channel_width})   // 23
);

blocks
(
    // Solid region (heatsink base) - center bottom
    hex (1 2 6 5 9 10 14 13) solid ({int(cells_x*base_length/channel_length)} {int(cells_y*base_height/channel_height)} {cells_z}) simpleGrading (1 1 1)

    // Fluid region - inlet section (before heatsink)
    hex (0 1 5 4 8 9 13 12) fluid ({int(cells_x*heatsink_start_x/channel_length)} {int(cells_y*base_height/channel_height)} {cells_z}) simpleGrading (1 1 1)
    // Fluid region - outlet section (after heatsink)
    hex (2 3 7 6 10 11 15 14) fluid ({int(cells_x*(channel_length-heatsink_end_x)/channel_length)} {int(cells_y*base_height/channel_height)} {cells_z}) simpleGrading (1 1 1)

    // Fluid region - above inlet
    hex (8 9 13 12 16 17 21 20) fluid ({int(cells_x*heatsink_start_x/channel_length)} {int(cells_y*(channel_height-base_height)/channel_height)} {cells_z}) simpleGrading (1 1 1)
    // Fluid region - above heatsink
    hex (9 10 14 13 17 18 22 21) fluid ({int(cells_x*base_length/channel_length)} {int(cells_y*(channel_height-base_height)/channel_height)} {cells_z}) simpleGrading (1 1 1)
    // Fluid region - above outlet
    hex (10 11 15 14 18 19 23 22) fluid ({int(cells_x*(channel_length-heatsink_end_x)/channel_length)} {int(cells_y*(channel_height-base_height)/channel_height)} {cells_z}) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 4 12 8)
            (8 12 20 16)
        );
    }}

    outlet
    {{
        type patch;
        faces
        (
            (3 11 15 7)
            (11 19 23 15)
        );
    }}

    topWall
    {{
        type wall;
        faces
        (
            (16 17 21 20)
            (17 18 22 21)
            (18 19 23 22)
        );
    }}

    bottomWall_fluid
    {{
        type wall;
        faces
        (
            (0 1 5 4)
            (2 3 7 6)
        );
    }}

    bottomWall_solid
    {{
        type wall;
        faces
        (
            (1 2 6 5)
        );
    }}

    frontAndBack
    {{
        type wall;
        faces
        (
            // Front (z=0)
            (0 8 9 1)
            (1 9 10 2)
            (2 10 11 3)
            (8 16 17 9)
            (9 17 18 10)
            (10 18 19 11)
            // Back (z=channel_width)
            (4 5 13 12)
            (5 6 14 13)
            (6 7 15 14)
            (12 13 21 20)
            (13 14 22 21)
            (14 15 23 22)
        );
    }}

    solid_to_fluid
    {{
        type wall;
        faces
        (
            (9 10 14 13)
        );
    }}

    fluid_to_solid
    {{
        type wall;
        faces
        (
            (9 10 14 13)
        );
    }}
);

mergePatchPairs ();
"""
    with open(os.path.join(case_dir, 'system', 'blockMeshDict'), 'w') as f:
        f.write(block_mesh)

    # ========================================
    # 0/fluid/U
    # ========================================
    fluid_U = f"""FoamFile
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

    bottomWall_fluid
    {{
        type            noSlip;
    }}

    frontAndBack
    {{
        type            noSlip;
    }}

    ".*_to_.*"
    {{
        type            noSlip;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'fluid', 'U'), 'w') as f:
        f.write(fluid_U)

    # ========================================
    # 0/fluid/p_rgh
    # ========================================
    fluid_p_rgh = """FoamFile
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

    bottomWall_fluid
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }

    frontAndBack
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }

    ".*_to_.*"
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }
}
"""
    with open(os.path.join(case_dir, '0', 'fluid', 'p_rgh'), 'w') as f:
        f.write(fluid_p_rgh)

    # ========================================
    # 0/fluid/p
    # ========================================
    fluid_p = """FoamFile
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
    with open(os.path.join(case_dir, '0', 'fluid', 'p'), 'w') as f:
        f.write(fluid_p)

    # ========================================
    # 0/fluid/T
    # ========================================
    fluid_T = f"""FoamFile
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

    bottomWall_fluid
    {{
        type            zeroGradient;
    }}

    frontAndBack
    {{
        type            zeroGradient;
    }}

    "fluid_to_solid"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        value           uniform {inlet_temp};
        Tnbr            T;
        kappaMethod     fluidThermo;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'fluid', 'T'), 'w') as f:
        f.write(fluid_T)

    # ========================================
    # 0/solid/T
    # ========================================
    solid_T = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {base_temp};

boundaryField
{{
    bottomWall_solid
    {{
        type            fixedValue;
        value           uniform {base_temp};
    }}

    frontAndBack
    {{
        type            zeroGradient;
    }}

    "solid_to_fluid"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        value           uniform {base_temp};
        Tnbr            T;
        kappaMethod     solidThermo;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'solid', 'T'), 'w') as f:
        f.write(solid_T)

    # ========================================
    # system/changeDictionaryDict for setting up coupled patches
    # ========================================
    change_dict = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}

dictionaryReplacement
{
}
"""
    with open(os.path.join(case_dir, 'system', 'changeDictionaryDict'), 'w') as f:
        f.write(change_dict)

    # ========================================
    # Allrun script
    # ========================================
    allrun = f"""#!/bin/bash
cd "${{0%/*}}" || exit
. ${{WM_PROJECT_DIR:=/usr/share/openfoam}}/bin/tools/RunFunctions 2>/dev/null || true

# Set OpenFOAM environment
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

echo "Creating mesh..."
blockMesh > log.blockMesh 2>&1

echo "Splitting mesh into regions..."
splitMeshRegions -cellZones -overwrite > log.splitMeshRegions 2>&1

echo "Checking mesh..."
checkMesh -allRegions > log.checkMesh 2>&1

echo "Running chtMultiRegionFoam..."
chtMultiRegionFoam > log.chtMultiRegionFoam 2>&1

echo "Converting to VTK..."
foamToVTK -allRegions > log.foamToVTK 2>&1

echo "Done!"
"""
    with open(os.path.join(case_dir, 'Allrun'), 'w') as f:
        f.write(allrun)
    os.chmod(os.path.join(case_dir, 'Allrun'), 0o755)

    # ========================================
    # Allclean script
    # ========================================
    allclean = """#!/bin/bash
cd "${0%/*}" || exit

rm -rf 0.* [1-9]* constant/polyMesh constant/*/polyMesh
rm -rf VTK postProcessing processor* log.*
rm -rf constant/cellToRegion
"""
    with open(os.path.join(case_dir, 'Allclean'), 'w') as f:
        f.write(allclean)
    os.chmod(os.path.join(case_dir, 'Allclean'), 0o755)

    print(f"CHT case created in {case_dir}")
    print(f"  Inlet velocity: {inlet_velocity} m/s")
    print(f"  Inlet temperature: {inlet_temp - 273.15:.1f} °C")
    print(f"  Base temperature: {base_temp - 273.15:.1f} °C")
    print(f"\nTo run: cd {case_dir} && ./Allrun")

    return case_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create CHT heatsink cooling case')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--inlet-velocity', type=float, default=0.1, help='Inlet velocity (m/s)')
    parser.add_argument('--inlet-temp', type=float, default=25.0, help='Inlet temperature (°C)')
    parser.add_argument('--base-temp', type=float, default=80.0, help='Base temperature (°C)')

    args = parser.parse_args()

    create_cht_heatsink_case(
        args.case_dir,
        inlet_velocity=args.inlet_velocity,
        inlet_temp=args.inlet_temp + 273.15,
        base_temp=args.base_temp + 273.15
    )
