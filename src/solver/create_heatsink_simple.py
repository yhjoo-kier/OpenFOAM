#!/usr/bin/env python3
"""
Create simpleFoam case for heatsink cooling analysis
Uses incompressible flow with temperature as passive scalar
"""

import os


def create_heatsink_simple_case(
    case_dir: str,
    # Geometry parameters (m)
    channel_length: float = 0.12,     # 120mm
    channel_height: float = 0.04,     # 40mm
    channel_width: float = 0.06,      # 60mm
    # Boundary conditions
    inlet_velocity: float = 0.1,      # m/s
    inlet_temp: float = 298.15,       # K (25°C)
    wall_temp: float = 353.15,        # K (80°C)
    # Mesh resolution
    cells_x: int = 60,
    cells_y: int = 20,
    cells_z: int = 30
):
    """Create a simpleFoam case with heat transfer."""

    os.makedirs(case_dir, exist_ok=True)

    # Create directories
    for d in ['system', 'constant', '0']:
        os.makedirs(os.path.join(case_dir, d), exist_ok=True)

    # Air properties at 25-80°C
    nu = 1.6e-5  # kinematic viscosity (m^2/s)
    Pr = 0.71    # Prandtl number
    DT = nu / Pr  # thermal diffusivity (m^2/s)

    # ========================================
    # system/blockMeshDict
    # ========================================
    block_mesh = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale   1;

// Simple rectangular channel
// Bottom wall is heated (heatsink surface at 80°C)

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
    hex (0 1 2 3 4 5 6 7) ({cells_x} {cells_y} {cells_z})
    simpleGrading (
        1  // x: uniform
        (
            (0.3 0.3 4)   // near bottom wall - finer
            (0.4 0.4 1)   // middle
            (0.3 0.3 0.25) // near top wall - finer
        )
        1  // z: uniform
    )
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
        type symmetryPlane;
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
        f.write(block_mesh)

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

application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         3000;
deltaT          1;
writeControl    timeStep;
writeInterval   500;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    scalarTransport
    {{
        type            scalarTransport;
        libs            ("libsolverFunctionObjects.so");

        enabled         true;
        writeControl    writeTime;

        field           T;

        // Diffusivity (thermal diffusivity for air)
        D               {DT};

        // Use phi from simpleFoam
        phi             phi;

        // Boundary conditions are defined in 0/T
        resetOnStartUp  false;

        fvOptions
        {{
        }}

        nCorr           2;
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
    div(phi,T)      bounded Gauss linearUpwind grad(T);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
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
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "(U|k|epsilon|omega|T)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0.1;
    }

    "(U|k|epsilon|omega|T)Final"
    {
        $U;
        relTol          0;
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
        "(k|epsilon|omega)" 1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.7;
        epsilon         0.7;
        omega           0.7;
        T               0.5;
    }
}
"""
    with open(os.path.join(case_dir, 'system', 'fvSolution'), 'w') as f:
        f.write(fv_solution)

    # ========================================
    # constant/transportProperties
    # ========================================
    transport_props = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}

// Air at ~50°C
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {nu};

// Thermal diffusivity for passive scalar transport
DT              [0 2 -1 0 0 0 0] {DT};
"""
    with open(os.path.join(case_dir, 'constant', 'transportProperties'), 'w') as f:
        f.write(transport_props)

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
        type            symmetryPlane;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'U'), 'w') as f:
        f.write(U_file)

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

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }

    topWall
    {
        type            zeroGradient;
    }

    heatedWall
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            symmetryPlane;
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
        type            symmetryPlane;
    }}
}}
"""
    with open(os.path.join(case_dir, '0', 'T'), 'w') as f:
        f.write(T_file)

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
echo "Heatsink Cooling Simulation (simpleFoam)"
echo "Channel: {channel_length*1000:.0f} x {channel_height*1000:.0f} x {channel_width*1000:.0f} mm"
echo "Inlet velocity: {inlet_velocity} m/s"
echo "Inlet temperature: {inlet_temp - 273.15:.1f} degC"
echo "Wall temperature: {wall_temp - 273.15:.1f} degC"
echo "=========================================="

echo ""
echo "Step 1: Creating mesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: blockMesh failed"
    cat log.blockMesh
    exit 1
fi
echo "  Mesh created successfully"

echo ""
echo "Step 2: Checking mesh quality..."
checkMesh > log.checkMesh 2>&1
grep -A 5 "Mesh OK" log.checkMesh || echo "  Check log.checkMesh for details"

echo ""
echo "Step 3: Running simpleFoam with scalarTransport..."
echo "  This may take a few minutes..."
simpleFoam > log.simpleFoam 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: simpleFoam failed. Check log.simpleFoam"
    tail -30 log.simpleFoam
    exit 1
fi

echo ""
echo "Step 4: Converting to VTK format..."
foamToVTK > log.foamToVTK 2>&1

echo ""
echo "=========================================="
echo "Simulation Complete!"
echo "=========================================="

# Extract final residuals
echo ""
echo "Final residuals:"
tail -20 log.simpleFoam | grep -E "(Ux|Uy|Uz|p|T)"

# Calculate outlet temperature
echo ""
echo "Extracting outlet temperature..."
postProcess -func 'patchAverage(name=outlet, T)' -latestTime > log.postProcess 2>&1
if [ -f postProcessing/patchAverage\(name=outlet,T\)/0/surfaceFieldValue.dat ]; then
    echo "Outlet average temperature:"
    tail -1 postProcessing/patchAverage\(name=outlet,T\)/0/surfaceFieldValue.dat
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
rm -rf 0.* [1-9]* constant/polyMesh VTK postProcessing log.* processor*
"""
    with open(os.path.join(case_dir, 'Allclean'), 'w') as f:
        f.write(allclean)
    os.chmod(os.path.join(case_dir, 'Allclean'), 0o755)

    print(f"Heatsink case created in {case_dir}")
    print(f"  Channel: {channel_length*1000:.0f} x {channel_height*1000:.0f} x {channel_width*1000:.0f} mm")
    print(f"  Inlet velocity: {inlet_velocity} m/s")
    print(f"  Inlet temperature: {inlet_temp - 273.15:.1f} °C")
    print(f"  Wall temperature: {wall_temp - 273.15:.1f} °C")
    print(f"  Air thermal diffusivity: {DT:.2e} m²/s")
    print(f"\nTo run: cd {case_dir} && ./Allrun")

    return case_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create heatsink cooling case')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--inlet-velocity', type=float, default=0.1,
                        help='Inlet velocity (m/s)')
    parser.add_argument('--inlet-temp', type=float, default=25.0,
                        help='Inlet temperature (°C)')
    parser.add_argument('--wall-temp', type=float, default=80.0,
                        help='Heated wall temperature (°C)')

    args = parser.parse_args()

    create_heatsink_simple_case(
        args.case_dir,
        inlet_velocity=args.inlet_velocity,
        inlet_temp=args.inlet_temp + 273.15,
        wall_temp=args.wall_temp + 273.15
    )
