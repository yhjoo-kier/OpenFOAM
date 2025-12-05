#!/usr/bin/env python3
"""
Simple Heatsink Geometry using blockMesh
Creates a heatsink with base plate and 3 rectangular fins for thermal analysis
Uses proper multi-block structure with matching nodes at block interfaces
"""

import os


def create_blockmesh_dict(
    case_dir: str,
    base_length: float = 0.060,    # 60 mm (X direction)
    base_width: float = 0.040,     # 40 mm (Y direction)
    base_height: float = 0.005,    # 5 mm base plate thickness
    fin_width: float = 0.006,      # 6 mm fin thickness (Y direction)
    fin_height: float = 0.020,     # 20 mm fin height (Z direction)
    num_fins: int = 3,
    cells_per_mm: float = 1.5      # mesh density
):
    """
    Create blockMeshDict for heatsink with base + 3 fins.

    Uses a structured multi-block approach where all blocks share nodes
    at interfaces. The domain is divided into 7 Y-regions and 2 Z-layers.

    Y-direction: [gap0][fin0][gap1][fin1][gap2][fin2][gap3]
    Z-direction: [base][fin_layer]

    Only fin regions have blocks in the fin_layer.
    """
    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)

    # Calculate fin positions (Y coordinates)
    total_fin_width = num_fins * fin_width
    total_gap = base_width - total_fin_width
    gap = total_gap / (num_fins + 1)

    # Y coordinates for all regions
    y_coords = [0]
    for i in range(num_fins):
        y_coords.append(y_coords[-1] + gap)       # start of fin
        y_coords.append(y_coords[-1] + fin_width) # end of fin
    y_coords.append(base_width)                    # final gap ends at base_width

    total_height = base_height + fin_height

    # Mesh cell counts
    cells_x = max(10, int(base_length * 1000 * cells_per_mm))
    cells_z_base = max(5, int(base_height * 1000 * cells_per_mm))
    cells_z_fin = max(15, int(fin_height * 1000 * cells_per_mm))

    # Y cells for each region
    cells_y = []
    for i in range(len(y_coords) - 1):
        width = y_coords[i+1] - y_coords[i]
        cells_y.append(max(3, int(width * 1000 * cells_per_mm)))

    blockmesh_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
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

// Heatsink geometry:
// Base plate: {base_length*1000:.0f}mm x {base_width*1000:.0f}mm x {base_height*1000:.0f}mm
// Fins: {num_fins} fins, each {fin_width*1000:.0f}mm wide, {fin_height*1000:.0f}mm tall
// Total height: {total_height*1000:.0f}mm

scale   1;  // Already in meters

vertices
(
"""

    # Generate vertices
    # For each Y-coordinate, we need 4 vertices: (0,y,0), (L,y,0), (0,y,base_h), (L,y,base_h)
    # Plus for fin regions at fin_layer: (0,y,total_h), (L,y,total_h)

    vertex_idx = 0
    base_vertices = {}  # (y_idx, z_level) -> vertex indices (front, back at x=0, x=L)
    fin_vertices = {}   # (y_idx, z_level) -> vertex indices for fin layer

    # Z levels: 0 = bottom, 1 = base_top, 2 = fin_top
    z_levels = [0, base_height, total_height]

    # For structured multi-block, all Y boundaries need vertices at z=0 and z=base_height
    # Only fin Y boundaries need vertices at z=total_height

    for y_idx, y in enumerate(y_coords):
        # Bottom (z=0)
        blockmesh_dict += f"    (0 {y} 0)           // {vertex_idx}: y{y_idx}, z=0, x=0\n"
        blockmesh_dict += f"    ({base_length} {y} 0)   // {vertex_idx+1}: y{y_idx}, z=0, x=L\n"
        base_vertices[(y_idx, 0)] = (vertex_idx, vertex_idx+1)
        vertex_idx += 2

        # Base top (z=base_height)
        blockmesh_dict += f"    (0 {y} {base_height})           // {vertex_idx}: y{y_idx}, z=base_h, x=0\n"
        blockmesh_dict += f"    ({base_length} {y} {base_height})   // {vertex_idx+1}: y{y_idx}, z=base_h, x=L\n"
        base_vertices[(y_idx, 1)] = (vertex_idx, vertex_idx+1)
        vertex_idx += 2

        # Fin top (z=total_height) - only for fin boundaries (odd indices: 1,2,3,4,5,6)
        # y_idx 1,2 are fin 0 boundaries, 3,4 are fin 1, 5,6 are fin 2
        if y_idx in [1, 2, 3, 4, 5, 6]:
            blockmesh_dict += f"    (0 {y} {total_height})           // {vertex_idx}: y{y_idx}, z=fin_top, x=0\n"
            blockmesh_dict += f"    ({base_length} {y} {total_height})   // {vertex_idx+1}: y{y_idx}, z=fin_top, x=L\n"
            fin_vertices[(y_idx, 2)] = (vertex_idx, vertex_idx+1)
            vertex_idx += 2

    blockmesh_dict += ");\n\nblocks\n(\n"

    # Base layer blocks (7 blocks for 7 Y-regions)
    for i in range(7):  # 7 Y-regions
        v00 = base_vertices[(i, 0)]     # y=y_i, z=0
        v01 = base_vertices[(i+1, 0)]   # y=y_{i+1}, z=0
        v10 = base_vertices[(i, 1)]     # y=y_i, z=base_h
        v11 = base_vertices[(i+1, 1)]   # y=y_{i+1}, z=base_h

        # hex: (x0y0z0, x1y0z0, x1y1z0, x0y1z0, x0y0z1, x1y0z1, x1y1z1, x0y1z1)
        blockmesh_dict += f"    // Base block {i} (Y region {i})\n"
        blockmesh_dict += f"    hex ({v00[0]} {v00[1]} {v01[1]} {v01[0]} {v10[0]} {v10[1]} {v11[1]} {v11[0]}) "
        blockmesh_dict += f"({cells_x} {cells_y[i]} {cells_z_base}) simpleGrading (1 1 1)\n"

    # Fin layer blocks (3 blocks for 3 fins)
    # Fin 0: Y region 1, Fin 1: Y region 3, Fin 2: Y region 5
    fin_regions = [1, 3, 5]
    for fin_idx, y_region in enumerate(fin_regions):
        y_idx_start = y_region
        y_idx_end = y_region + 1

        v_base_start = base_vertices[(y_idx_start, 1)]  # base top at fin start
        v_base_end = base_vertices[(y_idx_end, 1)]      # base top at fin end
        v_fin_start = fin_vertices[(y_idx_start, 2)]    # fin top at fin start
        v_fin_end = fin_vertices[(y_idx_end, 2)]        # fin top at fin end

        blockmesh_dict += f"\n    // Fin {fin_idx} block (Y region {y_region})\n"
        blockmesh_dict += f"    hex ({v_base_start[0]} {v_base_start[1]} {v_base_end[1]} {v_base_end[0]} "
        blockmesh_dict += f"{v_fin_start[0]} {v_fin_start[1]} {v_fin_end[1]} {v_fin_end[0]}) "
        blockmesh_dict += f"({cells_x} {cells_y[y_region]} {cells_z_fin}) simpleGrading (1 1 1)\n"

    blockmesh_dict += ");\n\nedges\n(\n);\n\nboundary\n(\n"

    # Bottom boundary (heat source - base plate bottom)
    blockmesh_dict += "    bottom\n    {\n        type wall;\n        faces\n        (\n"
    for i in range(7):
        v00 = base_vertices[(i, 0)]
        v01 = base_vertices[(i+1, 0)]
        blockmesh_dict += f"            ({v00[0]} {v01[0]} {v01[1]} {v00[1]})  // Base block {i} bottom\n"
    blockmesh_dict += "        );\n    }\n\n"

    # Top boundary (cooling - fin tops only)
    blockmesh_dict += "    top\n    {\n        type wall;\n        faces\n        (\n"
    for fin_idx, y_region in enumerate(fin_regions):
        y_idx_start = y_region
        y_idx_end = y_region + 1
        v_fin_start = fin_vertices[(y_idx_start, 2)]
        v_fin_end = fin_vertices[(y_idx_end, 2)]
        blockmesh_dict += f"            ({v_fin_start[0]} {v_fin_start[1]} {v_fin_end[1]} {v_fin_end[0]})  // Fin {fin_idx} top\n"
    blockmesh_dict += "        );\n    }\n\n"

    # Sides boundary (all other external faces)
    blockmesh_dict += "    sides\n    {\n        type wall;\n        faces\n        (\n"

    # Base front/back (X boundaries at base layer)
    for i in range(7):
        v00 = base_vertices[(i, 0)]
        v01 = base_vertices[(i+1, 0)]
        v10 = base_vertices[(i, 1)]
        v11 = base_vertices[(i+1, 1)]
        blockmesh_dict += f"            ({v00[0]} {v10[0]} {v11[0]} {v01[0]})  // Base block {i} front (x=0)\n"
        blockmesh_dict += f"            ({v00[1]} {v01[1]} {v11[1]} {v10[1]})  // Base block {i} back (x=L)\n"

    # Fin front/back (X boundaries at fin layer)
    for fin_idx, y_region in enumerate(fin_regions):
        y_idx_start = y_region
        y_idx_end = y_region + 1
        v_base_start = base_vertices[(y_idx_start, 1)]
        v_base_end = base_vertices[(y_idx_end, 1)]
        v_fin_start = fin_vertices[(y_idx_start, 2)]
        v_fin_end = fin_vertices[(y_idx_end, 2)]
        blockmesh_dict += f"            ({v_base_start[0]} {v_fin_start[0]} {v_fin_end[0]} {v_base_end[0]})  // Fin {fin_idx} front (x=0)\n"
        blockmesh_dict += f"            ({v_base_start[1]} {v_base_end[1]} {v_fin_end[1]} {v_fin_start[1]})  // Fin {fin_idx} back (x=L)\n"

    # Y boundaries at base layer (y=0 and y=base_width)
    v_y0_bot = base_vertices[(0, 0)]
    v_y0_top = base_vertices[(0, 1)]
    blockmesh_dict += f"            ({v_y0_bot[0]} {v_y0_bot[1]} {v_y0_top[1]} {v_y0_top[0]})  // Base left (y=0)\n"
    v_ymax_bot = base_vertices[(7, 0)]
    v_ymax_top = base_vertices[(7, 1)]
    blockmesh_dict += f"            ({v_ymax_bot[0]} {v_ymax_top[0]} {v_ymax_top[1]} {v_ymax_bot[1]})  // Base right (y=max)\n"

    # Fin Y-side faces (the 4 faces of each fin perpendicular to Y)
    for fin_idx, y_region in enumerate(fin_regions):
        y_idx_start = y_region
        y_idx_end = y_region + 1
        v_base_start = base_vertices[(y_idx_start, 1)]
        v_base_end = base_vertices[(y_idx_end, 1)]
        v_fin_start = fin_vertices[(y_idx_start, 2)]
        v_fin_end = fin_vertices[(y_idx_end, 2)]
        blockmesh_dict += f"            ({v_base_start[0]} {v_base_start[1]} {v_fin_start[1]} {v_fin_start[0]})  // Fin {fin_idx} front side (y=min)\n"
        blockmesh_dict += f"            ({v_base_end[0]} {v_fin_end[0]} {v_fin_end[1]} {v_base_end[1]})  // Fin {fin_idx} back side (y=max)\n"

    # Base top exposed surfaces (gaps between fins - Y regions 0, 2, 4, 6)
    gap_regions = [0, 2, 4, 6]
    for y_region in gap_regions:
        y_idx_start = y_region
        y_idx_end = y_region + 1
        v_top_start = base_vertices[(y_idx_start, 1)]
        v_top_end = base_vertices[(y_idx_end, 1)]
        blockmesh_dict += f"            ({v_top_start[0]} {v_top_end[0]} {v_top_end[1]} {v_top_start[1]})  // Base top gap region {y_region}\n"

    blockmesh_dict += """        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "blockMeshDict"), 'w') as f:
        f.write(blockmesh_dict)

    print(f"Created blockMeshDict in {system_dir}")
    print(f"  Base: {base_length*1000:.0f}mm x {base_width*1000:.0f}mm x {base_height*1000:.0f}mm")
    print(f"  Fins: {num_fins} fins, each {fin_width*1000:.0f}mm wide, {fin_height*1000:.0f}mm tall")
    print(f"  Y coordinates: {[f'{y*1000:.1f}mm' for y in y_coords]}")
    print(f"  Total blocks: {7 + num_fins} (7 base + {num_fins} fin)")

    return os.path.join(system_dir, "blockMeshDict")


def create_heat_conduction_case(
    case_dir: str,
    heat_flux: float = 10000,      # W/m² heat flux at bottom
    ambient_temp: float = 300,     # K ambient/cooling temperature
    k_aluminum: float = 205        # W/(m·K) thermal conductivity
):
    """
    Create a complete heat conduction case using laplacianFoam.
    """
    os.makedirs(case_dir, exist_ok=True)

    # Create blockMeshDict with 3-fin heatsink
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
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    fieldMinMax
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        fields          (T);
        writeControl    timeStep;
        writeInterval   100;
    }}
}}
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
        tolerance       1e-08;
        relTol          0.001;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 2;
}
"""
    with open(os.path.join(case_dir, "system", "fvSolution"), 'w') as f:
        f.write(fv_solution)

    # constant/transportProperties
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    transport_props = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

// Aluminum thermal diffusivity
// k = 205 W/(m·K), rho = 2700 kg/m³, Cp = 900 J/(kg·K)
// alpha = k / (rho * Cp) = 205 / (2700 * 900) = 8.44e-5 m²/s
DT              [0 2 -1 0 0 0 0] 8.44e-5;
"""
    with open(os.path.join(case_dir, "constant", "transportProperties"), 'w') as f:
        f.write(transport_props)

    # 0/T (temperature field with boundary conditions)
    os.makedirs(os.path.join(case_dir, "0"), exist_ok=True)

    temp_gradient = heat_flux / k_aluminum

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
        // Heat source: q = {heat_flux} W/m², k = {k_aluminum} W/(m·K)
        // gradient = q/k = {temp_gradient:.4f} K/m
        type            fixedGradient;
        gradient        uniform {temp_gradient:.4f};
    }}

    top
    {{
        // Fin tops: Fixed temperature (cooling)
        type            fixedValue;
        value           uniform {ambient_temp};
    }}

    sides
    {{
        // Adiabatic boundaries
        type            zeroGradient;
    }}
}}
"""
    with open(os.path.join(case_dir, "0", "T"), 'w') as f:
        f.write(t_field)

    print(f"\nHeat conduction case created in {case_dir}")
    print(f"  Heat flux: {heat_flux} W/m²")
    print(f"  Cooling temperature: {ambient_temp} K")
    print(f"  Material: Aluminum (k = {k_aluminum} W/(m·K))")
    print(f"  Temperature gradient at bottom: {temp_gradient:.2f} K/m")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create heatsink case with base + 3 fins')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--heat-flux', type=float, default=10000,
                        help='Heat flux at bottom (W/m²), default=10000')
    parser.add_argument('--ambient-temp', type=float, default=300,
                        help='Ambient/cooling temperature (K), default=300')

    args = parser.parse_args()

    create_heat_conduction_case(args.case_dir, args.heat_flux, args.ambient_temp)
