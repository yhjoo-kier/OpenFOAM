#!/usr/bin/env python3
"""
Simple Heatsink Geometry using blockMesh
Creates a simple heatsink with one fin for thermal analysis
"""

import os


def create_blockmesh_dict(
    case_dir: str,
    base_length: float = 0.05,    # 50 mm
    base_width: float = 0.05,     # 50 mm
    base_height: float = 0.005,   # 5 mm
    fin_height: float = 0.02,     # 20 mm
    fin_thickness: float = 0.01,  # 10 mm
    cells_x: int = 30,
    cells_y: int = 30,
    cells_z_base: int = 5,
    cells_z_fin: int = 15
):
    """
    Create blockMeshDict for simple heatsink (base + single fin).
    """
    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)

    # Fin position (center)
    fin_y_start = (base_width - fin_thickness) / 2
    fin_y_end = fin_y_start + fin_thickness

    total_height = base_height + fin_height

    blockmesh_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   1;

// Geometry parameters
// Base: {base_length*1000:.0f}mm x {base_width*1000:.0f}mm x {base_height*1000:.0f}mm
// Fin:  {base_length*1000:.0f}mm x {fin_thickness*1000:.0f}mm x {fin_height*1000:.0f}mm (centered)

vertices
(
    // Base vertices (z=0)
    (0 0 0)                              // 0
    ({base_length} 0 0)                  // 1
    ({base_length} {base_width} 0)       // 2
    (0 {base_width} 0)                   // 3

    // Base top / fin bottom (z=base_height)
    (0 0 {base_height})                  // 4
    ({base_length} 0 {base_height})      // 5
    ({base_length} {base_width} {base_height})  // 6
    (0 {base_width} {base_height})       // 7

    // Fin top vertices (z=total_height)
    (0 {fin_y_start} {total_height})                    // 8
    ({base_length} {fin_y_start} {total_height})        // 9
    ({base_length} {fin_y_end} {total_height})          // 10
    (0 {fin_y_end} {total_height})                      // 11

    // Fin bottom corners at base_height (for hex block)
    (0 {fin_y_start} {base_height})                     // 12
    ({base_length} {fin_y_start} {base_height})         // 13
    ({base_length} {fin_y_end} {base_height})           // 14
    (0 {fin_y_end} {base_height})                       // 15
);

blocks
(
    // Base block
    hex (0 1 2 3 4 5 6 7) ({cells_x} {cells_y} {cells_z_base}) simpleGrading (1 1 1)

    // Fin block
    hex (12 13 14 15 8 9 10 11) ({cells_x} {int(cells_y * fin_thickness / base_width)} {cells_z_fin}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    bottom
    {{
        type wall;
        faces
        (
            (0 3 2 1)
        );
    }}

    top
    {{
        type wall;
        faces
        (
            (8 9 10 11)
        );
    }}

    sides
    {{
        type wall;
        faces
        (
            // Base sides
            (0 1 5 4)
            (2 3 7 6)
            (0 4 7 3)
            (1 2 6 5)
            // Base top (excluding fin area) - simplified: include all
            (4 5 6 7)
            // Fin sides
            (12 8 11 15)
            (13 14 10 9)
            (12 13 9 8)
            (15 11 10 14)
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "blockMeshDict"), 'w') as f:
        f.write(blockmesh_dict)

    print(f"Created blockMeshDict in {system_dir}")
    return os.path.join(system_dir, "blockMeshDict")


def create_heat_conduction_case(case_dir: str, heat_flux: float = 5000, ambient_temp: float = 300):
    """
    Create a complete heat conduction case using laplacianFoam.
    """
    os.makedirs(case_dir, exist_ok=True)

    # Create blockMeshDict
    create_blockmesh_dict(case_dir)

    # system/controlDict
    control_dict = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}

application     laplacianFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;
deltaT          1;
writeControl    timeStep;
writeInterval   100;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
"""
    with open(os.path.join(case_dir, "system", "controlDict"), 'w') as f:
        f.write(control_dict)

    # system/fvSchemes
    fv_schemes = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes
{
    default         Euler;
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
    with open(os.path.join(case_dir, "system", "fvSchemes"), 'w') as f:
        f.write(fv_schemes)

    # system/fvSolution
    fv_solution = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    T
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
}
"""
    with open(os.path.join(case_dir, "system", "fvSolution"), 'w') as f:
        f.write(fv_solution)

    # constant/transportProperties (thermal diffusivity for aluminum)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    transport_props = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

// Aluminum thermal diffusivity
// alpha = k / (rho * Cp) = 205 / (2700 * 900) = 8.44e-5 m^2/s
DT              [0 2 -1 0 0 0 0] 8.44e-5;
"""
    with open(os.path.join(case_dir, "constant", "transportProperties"), 'w') as f:
        f.write(transport_props)

    # 0/T (temperature field)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)
    t_field = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {ambient_temp};

boundaryField
{{
    bottom
    {{
        type            fixedGradient;
        gradient        uniform {heat_flux / 205};  // q/k for aluminum
    }}

    top
    {{
        type            fixedValue;
        value           uniform {ambient_temp};
    }}

    sides
    {{
        type            zeroGradient;
    }}
}}
"""
    with open(os.path.join(case_dir, "0", "T"), 'w') as f:
        f.write(t_field)

    print(f"Heat conduction case created in {case_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create simple heatsink case')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--heat-flux', type=float, default=5000, help='Heat flux (W/m²)')

    args = parser.parse_args()

    create_heat_conduction_case(args.case_dir, args.heat_flux)
