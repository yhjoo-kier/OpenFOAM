#!/usr/bin/env python3
"""
CHT Heatsink Mesh Generation
Creates a combined fluid+solid mesh for chtMultiRegionFoam simulation
Heatsink (solid) with surrounding air channel (fluid)
"""

import os


def create_cht_blockmesh(
    case_dir: str,
    # Heatsink dimensions
    base_length: float = 0.040,    # 40 mm (X)
    base_width: float = 0.030,     # 30 mm (Y)
    base_height: float = 0.003,    # 3 mm base
    fin_width: float = 0.004,      # 4 mm fin thickness
    fin_height: float = 0.015,     # 15 mm fin height
    num_fins: int = 3,
    # Fluid channel dimensions
    channel_height: float = 0.025,  # 25 mm above base
    inlet_length: float = 0.010,    # 10 mm inlet region
    outlet_length: float = 0.020,   # 20 mm outlet region
    # Mesh density
    cells_per_mm: float = 1.0
):
    """
    Create blockMeshDict for CHT simulation.

    Domain structure:
    - Solid: heatsink (base + fins)
    - Fluid: air channel around and above heatsink

    For simplicity, we create a single-block fluid region above the heatsink
    and use topoSet to define solid/fluid zones.
    """
    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)

    # Total dimensions
    total_length = inlet_length + base_length + outlet_length
    total_height = base_height + max(fin_height, channel_height)
    total_width = base_width

    # Heatsink position (centered in X)
    hs_x_start = inlet_length
    hs_x_end = inlet_length + base_length

    # Calculate fin Y positions
    total_fin_width = num_fins * fin_width
    total_gap = base_width - total_fin_width
    gap = total_gap / (num_fins + 1)

    fin_y_coords = []
    for i in range(num_fins):
        y_start = gap + i * (fin_width + gap)
        y_end = y_start + fin_width
        fin_y_coords.append((y_start, y_end))

    # Mesh cells
    nx_inlet = max(5, int(inlet_length * 1000 * cells_per_mm))
    nx_hs = max(20, int(base_length * 1000 * cells_per_mm))
    nx_outlet = max(10, int(outlet_length * 1000 * cells_per_mm))
    ny = max(15, int(base_width * 1000 * cells_per_mm))
    nz_base = max(3, int(base_height * 1000 * cells_per_mm))
    nz_fin = max(10, int(fin_height * 1000 * cells_per_mm))
    nz_fluid = max(10, int((channel_height - fin_height) * 1000 * cells_per_mm)) if channel_height > fin_height else 5

    # For CHT, we'll create a simpler geometry:
    # Single domain with the heatsink as a box with fins
    # Then use topoSet to mark solid region

    blockmesh = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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

// CHT Heatsink Domain
// Solid: Base {base_length*1000:.0f}x{base_width*1000:.0f}x{base_height*1000:.0f}mm + {num_fins} fins
// Fluid: Air channel {total_length*1000:.0f}x{total_width*1000:.0f}x{total_height*1000:.0f}mm

scale 1;

vertices
(
    // Bottom plane (z=0)
    (0 0 0)                                  // 0
    ({total_length} 0 0)                     // 1
    ({total_length} {total_width} 0)         // 2
    (0 {total_width} 0)                      // 3

    // Base top plane (z=base_height)
    (0 0 {base_height})                      // 4
    ({total_length} 0 {base_height})         // 5
    ({total_length} {total_width} {base_height})  // 6
    (0 {total_width} {base_height})          // 7

    // Top plane (z=total_height)
    (0 0 {total_height})                     // 8
    ({total_length} 0 {total_height})        // 9
    ({total_length} {total_width} {total_height}) // 10
    (0 {total_width} {total_height})         // 11
);

blocks
(
    // Bottom layer (includes heatsink base region)
    hex (0 1 2 3 4 5 6 7) ({nx_inlet + nx_hs + nx_outlet} {ny} {nz_base})
    simpleGrading (1 1 1)

    // Top layer (fluid + fin region)
    hex (4 5 6 7 8 9 10 11) ({nx_inlet + nx_hs + nx_outlet} {ny} {nz_fin + nz_fluid})
    simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 3 7 4)
            (4 7 11 8)
        );
    }}

    outlet
    {{
        type patch;
        faces
        (
            (1 5 6 2)
            (5 9 10 6)
        );
    }}

    bottom
    {{
        type wall;
        faces
        (
            (0 1 2 3)
        );
    }}

    top
    {{
        type wall;
        faces
        (
            (8 11 10 9)
        );
    }}

    front
    {{
        type wall;
        faces
        (
            (0 4 5 1)
            (4 8 9 5)
        );
    }}

    back
    {{
        type wall;
        faces
        (
            (3 2 6 7)
            (7 6 10 11)
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "blockMeshDict"), 'w') as f:
        f.write(blockmesh)

    # Create topoSetDict to define solid region (heatsink)
    toposet = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    object      topoSetDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    // Create heatsink base cellSet
    {{
        name    heatsinkBase;
        type    cellSet;
        action  new;
        source  boxToCell;
        sourceInfo
        {{
            box ({hs_x_start} 0 0) ({hs_x_end} {total_width} {base_height});
        }}
    }}

    // Create heatsink fins cellSet
"""

    # Add fin regions
    fin_top = base_height + fin_height
    for i, (y_start, y_end) in enumerate(fin_y_coords):
        action = "new" if i == 0 else "add"
        toposet += f"""
    {{
        name    heatsinkFins;
        type    cellSet;
        action  {action};
        source  boxToCell;
        sourceInfo
        {{
            box ({hs_x_start} {y_start} {base_height}) ({hs_x_end} {y_end} {fin_top});
        }}
    }}
"""

    toposet += f"""
    // Combine base and fins into solid cellSet
    {{
        name    solid;
        type    cellSet;
        action  new;
        source  cellToCell;
        sourceInfo
        {{
            set heatsinkBase;
        }}
    }}

    {{
        name    solid;
        type    cellSet;
        action  add;
        source  cellToCell;
        sourceInfo
        {{
            set heatsinkFins;
        }}
    }}

    // Create solid cellZone
    {{
        name    solid;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        sourceInfo
        {{
            set solid;
        }}
    }}

    // Create fluid cellSet (everything not solid)
    {{
        name    fluid;
        type    cellSet;
        action  new;
        source  cellToCell;
        sourceInfo
        {{
            set solid;
        }}
    }}

    {{
        name    fluid;
        type    cellSet;
        action  invert;
    }}

    // Create fluid cellZone
    {{
        name    fluid;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        sourceInfo
        {{
            set fluid;
        }}
    }}
);

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "topoSetDict"), 'w') as f:
        f.write(toposet)

    print(f"Created CHT mesh files in {system_dir}")
    print(f"  Domain: {total_length*1000:.0f} x {total_width*1000:.0f} x {total_height*1000:.0f} mm")
    print(f"  Heatsink: {base_length*1000:.0f} x {base_width*1000:.0f} mm base + {num_fins} fins")
    print(f"  Fins: {fin_width*1000:.0f}mm wide, {fin_height*1000:.0f}mm tall")

    return {
        'total_length': total_length,
        'total_width': total_width,
        'total_height': total_height,
        'base_height': base_height,
        'fin_height': fin_height,
        'hs_x_start': hs_x_start,
        'hs_x_end': hs_x_end,
        'fin_y_coords': fin_y_coords
    }


def create_cht_case_files(case_dir: str, geometry: dict,
                          inlet_velocity: float = 0.5,
                          inlet_temp: float = 300,
                          heat_flux: float = 5000):
    """Create all case files for CHT simulation."""

    system_dir = os.path.join(case_dir, "system")
    const_dir = os.path.join(case_dir, "constant")
    zero_dir = os.path.join(case_dir, "0")

    for d in [system_dir, const_dir, zero_dir]:
        os.makedirs(d, exist_ok=True)

    # controlDict
    control_dict = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     chtMultiRegionFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         5;
deltaT          0.001;
writeControl    adjustableRunTime;
writeInterval   0.5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
maxCo           0.5;
maxDi           10;
adjustTimeStep  yes;
"""
    with open(os.path.join(system_dir, "controlDict"), 'w') as f:
        f.write(control_dict)

    # fvSchemes (global)
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
    with open(os.path.join(system_dir, "fvSchemes"), 'w') as f:
        f.write(fv_schemes)

    # fvSolution (global)
    fv_solution = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    "rho.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|h|e|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(U|h|e|k|epsilon)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    nOuterCorrectors    1;
    nCorrectors         2;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        "h.*"           1;
        "U.*"           1;
    }
}
"""
    with open(os.path.join(system_dir, "fvSolution"), 'w') as f:
        f.write(fv_solution)

    # regionProperties
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
    with open(os.path.join(const_dir, "regionProperties"), 'w') as f:
        f.write(region_props)

    # decomposeParDict
    decompose = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains 1;
method          simple;

simpleCoeffs
{
    n               (1 1 1);
    delta           0.001;
}
"""
    with open(os.path.join(system_dir, "decomposeParDict"), 'w') as f:
        f.write(decompose)

    # Create fluid region files
    create_fluid_region(case_dir, inlet_velocity, inlet_temp)

    # Create solid region files
    create_solid_region(case_dir, inlet_temp, heat_flux, geometry)

    print(f"\nCHT case files created in {case_dir}")
    print(f"  Inlet velocity: {inlet_velocity} m/s")
    print(f"  Inlet temperature: {inlet_temp} K")
    print(f"  Heat flux: {heat_flux} W/m²")


def create_fluid_region(case_dir: str, inlet_velocity: float, inlet_temp: float):
    """Create fluid region specific files."""

    # system/fluid
    sys_fluid = os.path.join(case_dir, "system", "fluid")
    os.makedirs(sys_fluid, exist_ok=True)

    fv_schemes_fluid = """FoamFile
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
    div(phi,U)      Gauss upwind;
    div(phi,K)      Gauss upwind;
    div(phi,h)      Gauss upwind;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,R)      Gauss upwind;
    div(R)          Gauss linear;
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
    with open(os.path.join(sys_fluid, "fvSchemes"), 'w') as f:
        f.write(fv_schemes_fluid)

    fv_solution_fluid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    "rho.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|h|e|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }

    "(U|h|e|k|epsilon)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    nOuterCorrectors    1;
    nCorrectors         2;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        "h.*"           1;
        "U.*"           1;
    }
}
"""
    with open(os.path.join(sys_fluid, "fvSolution"), 'w') as f:
        f.write(fv_solution_fluid)

    # constant/fluid
    const_fluid = os.path.join(case_dir, "constant", "fluid")
    os.makedirs(const_fluid, exist_ok=True)

    thermo_fluid = """FoamFile
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
        molWeight       28.9;
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
    with open(os.path.join(const_fluid, "thermophysicalProperties"), 'w') as f:
        f.write(thermo_fluid)

    turb_fluid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  laminar;
"""
    with open(os.path.join(const_fluid, "turbulenceProperties"), 'w') as f:
        f.write(turb_fluid)

    rad_fluid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      radiationProperties;
}

radiation       off;
radiationModel  none;
"""
    with open(os.path.join(const_fluid, "radiationProperties"), 'w') as f:
        f.write(rad_fluid)

    g_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}

dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);
"""
    with open(os.path.join(const_fluid, "g"), 'w') as f:
        f.write(g_file)

    # 0/fluid
    zero_fluid = os.path.join(case_dir, "0", "fluid")
    os.makedirs(zero_fluid, exist_ok=True)

    u_file = f"""FoamFile
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
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           $internalField;
    }}
    "(top|bottom|front|back)"
    {{
        type            noSlip;
    }}
    "fluid_to_solid"
    {{
        type            noSlip;
    }}
}}
"""
    with open(os.path.join(zero_fluid, "U"), 'w') as f:
        f.write(u_file)

    p_rgh_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}

dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;

boundaryField
{
    inlet
    {
        type            fixedFluxPressure;
        value           $internalField;
    }
    outlet
    {
        type            fixedValue;
        value           uniform 1e5;
    }
    "(top|bottom|front|back)"
    {
        type            fixedFluxPressure;
        value           $internalField;
    }
    "fluid_to_solid"
    {
        type            fixedFluxPressure;
        value           $internalField;
    }
}
"""
    with open(os.path.join(zero_fluid, "p_rgh"), 'w') as f:
        f.write(p_rgh_file)

    t_file = f"""FoamFile
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
        type            inletOutlet;
        inletValue      uniform {inlet_temp};
        value           $internalField;
    }}
    "(top|bottom|front|back)"
    {{
        type            zeroGradient;
    }}
    "fluid_to_solid"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     fluidThermo;
        value           uniform {inlet_temp};
    }}
}}
"""
    with open(os.path.join(zero_fluid, "T"), 'w') as f:
        f.write(t_file)

    p_file = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}

dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;

boundaryField
{
    ".*"
    {
        type            calculated;
        value           $internalField;
    }
}
"""
    with open(os.path.join(zero_fluid, "p"), 'w') as f:
        f.write(p_file)


def create_solid_region(case_dir: str, initial_temp: float, heat_flux: float, geometry: dict):
    """Create solid region specific files."""

    # system/solid
    sys_solid = os.path.join(case_dir, "system", "solid")
    os.makedirs(sys_solid, exist_ok=True)

    fv_schemes_solid = """FoamFile
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
    with open(os.path.join(sys_solid, "fvSchemes"), 'w') as f:
        f.write(fv_schemes_solid)

    fv_solution_solid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    "h.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0.1;
    }
}

PIMPLE
{
    nNonOrthogonalCorrectors 0;
}
"""
    with open(os.path.join(sys_solid, "fvSolution"), 'w') as f:
        f.write(fv_solution_solid)

    # constant/solid
    const_solid = os.path.join(case_dir, "constant", "solid")
    os.makedirs(const_solid, exist_ok=True)

    thermo_solid = """FoamFile
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
        molWeight       26.98;
    }
    transport
    {
        kappa           205;
    }
    thermodynamics
    {
        Hf              0;
        Cp              900;
    }
    equationOfState
    {
        rho             2700;
    }
}
"""
    with open(os.path.join(const_solid, "thermophysicalProperties"), 'w') as f:
        f.write(thermo_solid)

    rad_solid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      radiationProperties;
}

radiation       off;
radiationModel  none;
"""
    with open(os.path.join(const_solid, "radiationProperties"), 'w') as f:
        f.write(rad_solid)

    # 0/solid
    zero_solid = os.path.join(case_dir, "0", "solid")
    os.makedirs(zero_solid, exist_ok=True)

    # Heat flux boundary condition at bottom of heatsink
    t_file = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}

dimensions      [0 0 0 1 0 0 0];
internalField   uniform {initial_temp};

boundaryField
{{
    bottom
    {{
        type            externalWallHeatFluxTemperature;
        mode            flux;
        q               uniform {heat_flux};
        kappaMethod     solidThermo;
        value           uniform {initial_temp};
    }}
    "(inlet|outlet|top|front|back)"
    {{
        type            zeroGradient;
    }}
    "solid_to_fluid"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     solidThermo;
        value           uniform {initial_temp};
    }}
}}
"""
    with open(os.path.join(zero_solid, "T"), 'w') as f:
        f.write(t_file)


def create_changeDictionaryDict(case_dir: str):
    """Create changeDictionaryDict for fixing BCs after splitMeshRegions."""

    sys_fluid = os.path.join(case_dir, "system", "fluid")
    sys_solid = os.path.join(case_dir, "system", "solid")

    os.makedirs(sys_fluid, exist_ok=True)
    os.makedirs(sys_solid, exist_ok=True)

    # Fluid changeDictionaryDict
    change_fluid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}

T
{
    internalField   uniform 300;

    boundaryField
    {
        ".*"
        {
            type            zeroGradient;
        }

        inlet
        {
            type            fixedValue;
            value           uniform 300;
        }

        outlet
        {
            type            inletOutlet;
            inletValue      uniform 300;
            value           uniform 300;
        }

        "fluid_to_solid"
        {
            type            compressible::turbulentTemperatureCoupledBaffleMixed;
            Tnbr            T;
            kappaMethod     fluidThermo;
            value           uniform 300;
        }
    }
}

U
{
    internalField   uniform (0.5 0 0);

    boundaryField
    {
        ".*"
        {
            type            noSlip;
        }

        inlet
        {
            type            fixedValue;
            value           uniform (0.5 0 0);
        }

        outlet
        {
            type            inletOutlet;
            inletValue      uniform (0 0 0);
            value           uniform (0.5 0 0);
        }

        "fluid_to_solid"
        {
            type            noSlip;
        }
    }
}

p_rgh
{
    internalField   uniform 1e5;

    boundaryField
    {
        ".*"
        {
            type            fixedFluxPressure;
            value           uniform 1e5;
        }

        outlet
        {
            type            fixedValue;
            value           uniform 1e5;
        }
    }
}

p
{
    internalField   uniform 1e5;

    boundaryField
    {
        ".*"
        {
            type            calculated;
            value           uniform 1e5;
        }
    }
}
"""
    with open(os.path.join(sys_fluid, "changeDictionaryDict"), 'w') as f:
        f.write(change_fluid)

    # Solid changeDictionaryDict
    change_solid = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}

T
{
    internalField   uniform 300;

    boundaryField
    {
        ".*"
        {
            type            zeroGradient;
        }

        bottom
        {
            type            externalWallHeatFluxTemperature;
            mode            flux;
            q               uniform 5000;
            kappaMethod     solidThermo;
            value           uniform 300;
        }

        "solid_to_fluid"
        {
            type            compressible::turbulentTemperatureCoupledBaffleMixed;
            Tnbr            T;
            kappaMethod     solidThermo;
            value           uniform 300;
        }
    }
}
"""
    with open(os.path.join(sys_solid, "changeDictionaryDict"), 'w') as f:
        f.write(change_solid)


def setup_full_cht_case(case_dir: str,
                        inlet_velocity: float = 0.5,
                        inlet_temp: float = 300,
                        heat_flux: float = 5000):
    """Set up complete CHT case."""

    print("=" * 60)
    print("Setting up CHT Heatsink Case")
    print("=" * 60)

    # Create mesh
    geometry = create_cht_blockmesh(case_dir)

    # Create case files
    create_cht_case_files(case_dir, geometry, inlet_velocity, inlet_temp, heat_flux)

    # Create changeDictionaryDict
    create_changeDictionaryDict(case_dir)

    print("\n" + "=" * 60)
    print("CHT Case Setup Complete!")
    print("=" * 60)
    print("\nNext steps:")
    print("  1. blockMesh")
    print("  2. topoSet")
    print("  3. splitMeshRegions -cellZones -overwrite")
    print("  4. changeDictionary -region fluid")
    print("  5. changeDictionary -region solid")
    print("  6. chtMultiRegionFoam")

    return geometry


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Setup CHT heatsink case')
    parser.add_argument('case_dir', help='Case directory')
    parser.add_argument('--velocity', type=float, default=0.5, help='Inlet velocity (m/s)')
    parser.add_argument('--temp', type=float, default=300, help='Inlet temperature (K)')
    parser.add_argument('--heat-flux', type=float, default=5000, help='Heat flux (W/m²)')

    args = parser.parse_args()

    setup_full_cht_case(args.case_dir, args.velocity, args.temp, args.heat_flux)
