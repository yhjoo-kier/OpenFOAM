#!/usr/bin/env python3
"""
Setup OpenFOAM case files for heatsink CHT laminar simulation.

This script configures:
- Initial and boundary conditions (0/ directory)
- Thermophysical properties (constant/ directory)
- Region properties for multi-region CHT

Usage:
    # After gmshToFoam and splitMeshRegions:
    python3 setup_openfoam.py --post-split

    # Custom parameters:
    python3 setup_openfoam.py --post-split --Tinlet 293 --heatFlux 60000
"""

import argparse
import os
import glob
import re


# Default parameters
DEFAULTS = {
    "Tinlet": 293.15,     # Inlet temperature (K) = 20°C
    "Tsolid": 293.15,     # Initial solid temperature (K)
    "p_inlet": 10.0,      # Inlet pressure (Pa)
    "p_outlet": 0.0,      # Outlet pressure (Pa)
    "heatFlux": 60000.0,  # Heat flux at base bottom (W/m²)
}


def foam_header(class_type, location, obj_name):
    """Generate OpenFOAM file header."""
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {class_type};
    location    "{location}";
    object      {obj_name};
}}
"""


def write_fluid_T(case_dir, params, solid_regions):
    """Write fluid temperature field."""
    content = foam_header("volScalarField", "0/fluid", "T")
    content += f"""
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {params['Tinlet']};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {params['Tinlet']};
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    fluid_top
    {{
        type            zeroGradient;
    }}
    walls_fluid
    {{
        type            zeroGradient;
    }}
    tube_walls
    {{
        type            zeroGradient;
    }}
"""
    # Add interface patches for each solid region
    for region in solid_regions:
        content += f"""    fluid_to_{region}
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     fluidThermo;
        value           uniform {params['Tinlet']};
    }}
"""
    content += "}\n"

    filepath = os.path.join(case_dir, "0", "fluid", "T")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fluid_U(case_dir, params, solid_regions):
    """Write fluid velocity field."""
    content = foam_header("volVectorField", "0/fluid", "U")
    content += f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
    inlet
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}
    outlet
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}
    fluid_top
    {{
        type            noSlip;
    }}
    walls_fluid
    {{
        type            noSlip;
    }}
    tube_walls
    {{
        type            noSlip;
    }}
"""
    for region in solid_regions:
        content += f"""    fluid_to_{region}
    {{
        type            noSlip;
    }}
"""
    content += "}\n"

    filepath = os.path.join(case_dir, "0", "fluid", "U")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fluid_p(case_dir, params, solid_regions):
    """Write fluid pressure field."""
    content = foam_header("volScalarField", "0/fluid", "p")
    content += f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform {params['p_outlet']};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {params['p_inlet']};
    }}
    outlet
    {{
        type            fixedValue;
        value           uniform {params['p_outlet']};
    }}
    fluid_top
    {{
        type            fixedFluxPressure;
    }}
    walls_fluid
    {{
        type            fixedFluxPressure;
    }}
    tube_walls
    {{
        type            fixedFluxPressure;
    }}
"""
    for region in solid_regions:
        content += f"""    fluid_to_{region}
    {{
        type            fixedFluxPressure;
    }}
"""
    content += "}\n"

    filepath = os.path.join(case_dir, "0", "fluid", "p")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fluid_p_rgh(case_dir, params, solid_regions):
    """Write fluid p_rgh field."""
    content = foam_header("volScalarField", "0/fluid", "p_rgh")
    content += f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform {params['p_outlet']};

boundaryField
{{
    inlet
    {{
        type            fixedValue;
        value           uniform {params['p_inlet']};
    }}
    outlet
    {{
        type            fixedValue;
        value           uniform {params['p_outlet']};
    }}
    fluid_top
    {{
        type            fixedFluxPressure;
    }}
    walls_fluid
    {{
        type            fixedFluxPressure;
    }}
    tube_walls
    {{
        type            fixedFluxPressure;
    }}
"""
    for region in solid_regions:
        content += f"""    fluid_to_{region}
    {{
        type            fixedFluxPressure;
    }}
"""
    content += "}\n"

    filepath = os.path.join(case_dir, "0", "fluid", "p_rgh")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_solid_T(case_dir, params, region_name):
    """Write solid temperature field."""
    content = foam_header("volScalarField", f"0/{region_name}", "T")
    content += f"""
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {params['Tsolid']};

boundaryField
{{
    base_bottom
    {{
        type            externalWallHeatFluxTemperature;
        mode            flux;
        q               uniform {params['heatFlux']};
        kappaMethod     solidThermo;
        value           uniform {params['Tsolid']};
    }}
    walls_solid
    {{
        type            zeroGradient;
    }}
    {region_name}_to_fluid
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     solidThermo;
        value           uniform {params['Tsolid']};
    }}
}}
"""

    filepath = os.path.join(case_dir, "0", region_name, "T")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fluid_thermophysical(case_dir):
    """Write fluid thermophysical properties (water)."""
    content = foam_header("dictionary", "constant/fluid", "thermophysicalProperties")
    content += """
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight   18.015;  // Water
    }
    thermodynamics
    {
        Cp          4182;    // J/kg/K
        Hf          0;
    }
    transport
    {
        mu          0.001;   // Pa.s (water at 20°C)
        Pr          7.0;
        kappa       0.6;     // W/m/K
    }
    equationOfState
    {
        rho         998;     // kg/m³
    }
}
"""

    filepath = os.path.join(case_dir, "constant", "fluid", "thermophysicalProperties")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fluid_turbulence(case_dir):
    """Write fluid turbulence properties (laminar)."""
    content = foam_header("dictionary", "constant/fluid", "turbulenceProperties")
    content += """
simulationType  laminar;
"""

    filepath = os.path.join(case_dir, "constant", "fluid", "turbulenceProperties")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_solid_thermophysical(case_dir, region_name):
    """Write solid thermophysical properties (aluminum)."""
    content = foam_header("dictionary", f"constant/{region_name}", "thermophysicalProperties")
    content += """
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
        molWeight   26.98;   // Aluminum
    }
    thermodynamics
    {
        Cp          900;     // J/kg/K
        Hf          0;
    }
    transport
    {
        kappa       200;     // W/m/K (aluminum)
    }
    equationOfState
    {
        rho         2700;    // kg/m³
    }
}
"""

    filepath = os.path.join(case_dir, "constant", region_name, "thermophysicalProperties")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_solid_turbulence(case_dir, region_name):
    """Write solid turbulence properties."""
    content = foam_header("dictionary", f"constant/{region_name}", "turbulenceProperties")
    content += """
simulationType  laminar;
"""

    filepath = os.path.join(case_dir, "constant", region_name, "turbulenceProperties")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_region_properties(case_dir, solid_regions):
    """Write regionProperties for multi-region case."""
    solid_list = " ".join(solid_regions)
    content = foam_header("dictionary", "constant", "regionProperties")
    content += f"""
regions
(
    fluid       (fluid)
    solid       ({solid_list})
);
"""

    filepath = os.path.join(case_dir, "constant", "regionProperties")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_g(case_dir):
    """Write gravitational acceleration."""
    content = foam_header("uniformDimensionedVectorField", "constant", "g")
    content += """
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);
"""

    filepath = os.path.join(case_dir, "constant", "g")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSchemes_fluid(case_dir):
    """Write fvSchemes for fluid region."""
    content = foam_header("dictionary", "system/fluid", "fvSchemes")
    content += """
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
    grad(p)         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,h)      bounded Gauss linearUpwind default;
    div(phi,K)      bounded Gauss linearUpwind default;
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

    filepath = os.path.join(case_dir, "system", "fluid", "fvSchemes")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSolution_fluid(case_dir):
    """Write fvSolution for fluid region."""
    content = foam_header("dictionary", "system/fluid", "fvSolution")
    content += """
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }

    "(U|h|k|epsilon|omega)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pMinFactor      0.1;
    pMaxFactor      2.0;
    rhoMin          900;
    rhoMax          1100;

    residualControl
    {
        p               1e-5;
        U               1e-5;
        h               1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
        rho             0.5;
    }
    equations
    {
        U               0.7;
        h               0.7;
    }
}
"""

    filepath = os.path.join(case_dir, "system", "fluid", "fvSolution")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSchemes_solid(case_dir, region_name):
    """Write fvSchemes for solid region."""
    content = foam_header("dictionary", f"system/{region_name}", "fvSchemes")
    content += """
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

    filepath = os.path.join(case_dir, "system", region_name, "fvSchemes")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSolution_solid(case_dir, region_name):
    """Write fvSolution for solid region."""
    content = foam_header("dictionary", f"system/{region_name}", "fvSolution")
    content += """
solvers
{
    h
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        h               0.9;
    }
}
"""

    filepath = os.path.join(case_dir, "system", region_name, "fvSolution")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_controlDict(case_dir):
    """Write main controlDict."""
    content = foam_header("dictionary", "system", "controlDict")
    content += """
application     chtMultiRegionSimpleFoam;

startFrom       latestTime;

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
{
    fieldAverage
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        region          fluid;
        fields
        (
            U
            {
                mean        on;
                prime2Mean  off;
                base        time;
            }
            T
            {
                mean        on;
                prime2Mean  off;
                base        time;
            }
        );
    }
}
"""

    filepath = os.path.join(case_dir, "system", "controlDict")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSchemes_main(case_dir):
    """Write main fvSchemes (fallback)."""
    content = foam_header("dictionary", "system", "fvSchemes")
    content += """
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

    filepath = os.path.join(case_dir, "system", "fvSchemes")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_fvSolution_main(case_dir):
    """Write main fvSolution (fallback)."""
    content = foam_header("dictionary", "system", "fvSolution")
    content += """
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
}
"""

    filepath = os.path.join(case_dir, "system", "fvSolution")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def write_decomposeParDict(case_dir):
    """Write decomposeParDict for parallel runs."""
    content = foam_header("dictionary", "system", "decomposeParDict")
    content += """
numberOfSubdomains 4;

method          scotch;
"""

    filepath = os.path.join(case_dir, "system", "decomposeParDict")
    with open(filepath, "w") as f:
        f.write(content)
    print(f"Written: {filepath}")


def detect_solid_regions(case_dir):
    """Detect solid regions from constant/ directory after splitMeshRegions."""
    regions = []
    constant_dir = os.path.join(case_dir, "constant")

    if not os.path.exists(constant_dir):
        return ["solid"]

    for item in os.listdir(constant_dir):
        item_path = os.path.join(constant_dir, item)
        if os.path.isdir(item_path) and item not in ["fluid", "polyMesh"]:
            # Check if it has polyMesh subdirectory
            if os.path.exists(os.path.join(item_path, "polyMesh")):
                regions.append(item)

    if not regions:
        # Check for 'solid' region specifically
        if os.path.exists(os.path.join(constant_dir, "solid", "polyMesh")):
            regions.append("solid")

    return regions if regions else ["solid"]


def fix_boundary_types(case_dir, region):
    """Fix boundary types in polyMesh/boundary file after splitMeshRegions."""
    boundary_file = os.path.join(case_dir, "constant", region, "polyMesh", "boundary")

    if not os.path.exists(boundary_file):
        print(f"Warning: {boundary_file} not found")
        return

    with open(boundary_file, "r") as f:
        content = f.read()

    # Replace patch types with appropriate types
    # inlet and outlet should be 'patch'
    # walls should be 'wall'
    # interface patches should be 'mappedWall'

    # This is a simplified fix - actual implementation may need more sophisticated parsing
    print(f"Boundary file exists: {boundary_file}")


def setup_post_split(case_dir, params):
    """Setup case after splitMeshRegions has been run."""

    # Detect solid regions
    solid_regions = detect_solid_regions(case_dir)
    print(f"Detected solid regions: {solid_regions}")

    # Write all files
    write_region_properties(case_dir, solid_regions)
    write_g(case_dir)

    # Fluid files
    write_fluid_T(case_dir, params, solid_regions)
    write_fluid_U(case_dir, params, solid_regions)
    write_fluid_p(case_dir, params, solid_regions)
    write_fluid_p_rgh(case_dir, params, solid_regions)
    write_fluid_thermophysical(case_dir)
    write_fluid_turbulence(case_dir)
    write_fvSchemes_fluid(case_dir)
    write_fvSolution_fluid(case_dir)

    # Solid files for each region
    for region in solid_regions:
        write_solid_T(case_dir, params, region)
        write_solid_thermophysical(case_dir, region)
        write_solid_turbulence(case_dir, region)
        write_fvSchemes_solid(case_dir, region)
        write_fvSolution_solid(case_dir, region)

    # Main system files
    write_controlDict(case_dir)
    write_fvSchemes_main(case_dir)
    write_fvSolution_main(case_dir)
    write_decomposeParDict(case_dir)

    print("\nSetup complete!")
    print(f"Solid regions: {solid_regions}")
    print(f"Run: chtMultiRegionSimpleFoam")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Setup OpenFOAM case for heatsink CHT laminar simulation"
    )
    parser.add_argument(
        "--post-split",
        action="store_true",
        help="Run setup after splitMeshRegions (required)"
    )
    parser.add_argument(
        "--Tinlet",
        type=float,
        default=DEFAULTS["Tinlet"],
        help=f"Inlet temperature in K (default: {DEFAULTS['Tinlet']})"
    )
    parser.add_argument(
        "--Tsolid",
        type=float,
        default=DEFAULTS["Tsolid"],
        help=f"Initial solid temperature in K (default: {DEFAULTS['Tsolid']})"
    )
    parser.add_argument(
        "--p-inlet",
        type=float,
        default=DEFAULTS["p_inlet"],
        help=f"Inlet pressure in Pa (default: {DEFAULTS['p_inlet']})"
    )
    parser.add_argument(
        "--p-outlet",
        type=float,
        default=DEFAULTS["p_outlet"],
        help=f"Outlet pressure in Pa (default: {DEFAULTS['p_outlet']})"
    )
    parser.add_argument(
        "--heatFlux",
        type=float,
        default=DEFAULTS["heatFlux"],
        help=f"Heat flux at base bottom in W/m² (default: {DEFAULTS['heatFlux']})"
    )
    parser.add_argument(
        "--case-dir",
        type=str,
        default=".",
        help="Case directory (default: current directory)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    params = {
        "Tinlet": args.Tinlet,
        "Tsolid": args.Tsolid,
        "p_inlet": args.p_inlet,
        "p_outlet": args.p_outlet,
        "heatFlux": args.heatFlux,
    }

    case_dir = os.path.abspath(args.case_dir)
    print(f"Case directory: {case_dir}")

    if args.post_split:
        setup_post_split(case_dir, params)
    else:
        print("Please run with --post-split option after running splitMeshRegions")
        print("Usage: python3 setup_openfoam.py --post-split")


if __name__ == "__main__":
    main()
