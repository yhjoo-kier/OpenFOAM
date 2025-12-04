#!/usr/bin/env python3
"""
CHT (Conjugate Heat Transfer) Case Setup for OpenFOAM
Sets up chtMultiRegionFoam case with fluid and solid regions
"""

import os
import shutil


def create_system_files(case_dir: str, end_time: float = 100, delta_t: float = 0.1):
    """Create system directory files for CHT simulation."""

    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)

    # controlDict
    control_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     chtMultiRegionFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         {end_time};

deltaT          {delta_t};

writeControl    adjustableRunTime;

writeInterval   10;

purgeWrite      0;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

maxCo           0.6;

functions
{{
    fieldAverage1
    {{
        type            fieldAverage;
        functionObjectLibs ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        fields
        (
            T
            {{
                mean        on;
                prime2Mean  on;
                base        time;
            }}
        );
    }}
}}

// ************************************************************************* //
"""

    # fvSchemes (main - for region cases)
    fv_schemes = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

// ************************************************************************* //
"""

    # fvSolution (main)
    fv_solution = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    rho
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.1;
    }

    rhoFinal
    {
        $rho;
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

PIMPLE
{
    momentumPredictor   yes;
    nCorrectors         2;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        ".*"            1;
    }
}

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "controlDict"), 'w') as f:
        f.write(control_dict)
    with open(os.path.join(system_dir, "fvSchemes"), 'w') as f:
        f.write(fv_schemes)
    with open(os.path.join(system_dir, "fvSolution"), 'w') as f:
        f.write(fv_solution)


def create_fluid_region_files(case_dir: str, inlet_velocity: float = 1.0,
                               inlet_temp: float = 300):
    """Create files for fluid region."""

    fluid_dir = os.path.join(case_dir, "0", "fluid")
    os.makedirs(fluid_dir, exist_ok=True)

    # U (velocity)
    u_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volVectorField;
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        value           uniform ({inlet_velocity} 0 0);
    }}

    top
    {{
        type            slip;
    }}

    frontAndBack
    {{
        type            slip;
    }}

    bottomFluid
    {{
        type            noSlip;
    }}

    "fluid_to_solid.*"
    {{
        type            noSlip;
    }}
}}

// ************************************************************************* //
"""

    # p_rgh (pressure)
    p_rgh_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 1e5;

boundaryField
{
    inlet
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 1e5;
    }

    top
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
    }

    frontAndBack
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
    }

    bottomFluid
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
    }

    "fluid_to_solid.*"
    {
        type            fixedFluxPressure;
        value           uniform 1e5;
    }
}

// ************************************************************************* //
"""

    # T (temperature)
    t_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volScalarField;
    object      T;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        value           uniform {inlet_temp};
    }}

    top
    {{
        type            zeroGradient;
    }}

    frontAndBack
    {{
        type            zeroGradient;
    }}

    bottomFluid
    {{
        type            zeroGradient;
    }}

    "fluid_to_solid.*"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     fluidThermo;
        value           uniform {inlet_temp};
    }}
}}

// ************************************************************************* //
"""

    with open(os.path.join(fluid_dir, "U"), 'w') as f:
        f.write(u_file)
    with open(os.path.join(fluid_dir, "p_rgh"), 'w') as f:
        f.write(p_rgh_file)
    with open(os.path.join(fluid_dir, "T"), 'w') as f:
        f.write(t_file)

    # Create fluid thermophysicalProperties
    const_fluid = os.path.join(case_dir, "constant", "fluid")
    os.makedirs(const_fluid, exist_ok=True)

    thermo_props = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        molWeight       28.9;  // Air
    }
    thermodynamics
    {
        Cp              1005;  // J/(kg·K)
        Hf              0;
    }
    transport
    {
        mu              1.8e-5;  // Dynamic viscosity Pa·s
        Pr              0.71;    // Prandtl number
    }
}

// ************************************************************************* //
"""

    with open(os.path.join(const_fluid, "thermophysicalProperties"), 'w') as f:
        f.write(thermo_props)

    # g (gravity) for fluid
    g_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);

// ************************************************************************* //
"""
    with open(os.path.join(const_fluid, "g"), 'w') as f:
        f.write(g_file)


def create_solid_region_files(case_dir: str, initial_temp: float = 300,
                               heat_flux: float = 5000):
    """Create files for solid region (heatsink)."""

    solid_dir = os.path.join(case_dir, "0", "solid")
    os.makedirs(solid_dir, exist_ok=True)

    # T (temperature) for solid
    t_file = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       volScalarField;
    object      T;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform {initial_temp};

boundaryField
{{
    heatsinkBottom
    {{
        type            externalWallHeatFluxTemperature;
        mode            flux;
        q               uniform {heat_flux};  // W/m^2
        kappaMethod     solidThermo;
        value           uniform {initial_temp};
    }}

    "solid_to_fluid.*"
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     solidThermo;
        value           uniform {initial_temp};
    }}
}}

// ************************************************************************* //
"""

    with open(os.path.join(solid_dir, "T"), 'w') as f:
        f.write(t_file)

    # Create solid thermophysicalProperties (Aluminum)
    const_solid = os.path.join(case_dir, "constant", "solid")
    os.makedirs(const_solid, exist_ok=True)

    thermo_props = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        molWeight       26.98;  // Aluminum
    }

    transport
    {
        kappa           205;  // Thermal conductivity W/(m·K)
    }

    thermodynamics
    {
        Hf              0;
        Cp              900;  // Specific heat J/(kg·K)
    }

    equationOfState
    {
        rho             2700;  // Density kg/m^3
    }
}

// ************************************************************************* //
"""

    with open(os.path.join(const_solid, "thermophysicalProperties"), 'w') as f:
        f.write(thermo_props)


def create_region_properties(case_dir: str):
    """Create regionProperties file for multi-region simulation."""

    const_dir = os.path.join(case_dir, "constant")
    os.makedirs(const_dir, exist_ok=True)

    region_props = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      regionProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

regions
(
    fluid       (fluid)
    solid       (solid)
);

// ************************************************************************* //
"""

    with open(os.path.join(const_dir, "regionProperties"), 'w') as f:
        f.write(region_props)


def setup_cht_case(case_dir: str, inlet_velocity: float = 1.0,
                   inlet_temp: float = 300, heat_flux: float = 5000,
                   end_time: float = 100):
    """
    Set up complete CHT case for heatsink simulation.

    Args:
        case_dir: OpenFOAM case directory
        inlet_velocity: Air inlet velocity (m/s)
        inlet_temp: Air inlet temperature (K)
        heat_flux: Heat flux at heatsink bottom (W/m^2)
        end_time: Simulation end time (s)
    """
    print(f"Setting up CHT case in: {case_dir}")
    print(f"  Inlet velocity: {inlet_velocity} m/s")
    print(f"  Inlet temperature: {inlet_temp} K")
    print(f"  Heat flux: {heat_flux} W/m²")

    create_system_files(case_dir, end_time=end_time)
    create_region_properties(case_dir)
    create_fluid_region_files(case_dir, inlet_velocity, inlet_temp)
    create_solid_region_files(case_dir, inlet_temp, heat_flux)

    print("CHT case setup complete!")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Setup CHT case for OpenFOAM')
    parser.add_argument('case_dir', help='OpenFOAM case directory')
    parser.add_argument('--velocity', type=float, default=1.0, help='Inlet velocity (m/s)')
    parser.add_argument('--temp', type=float, default=300, help='Inlet temperature (K)')
    parser.add_argument('--heat-flux', type=float, default=5000, help='Heat flux (W/m²)')
    parser.add_argument('--end-time', type=float, default=100, help='End time (s)')

    args = parser.parse_args()

    setup_cht_case(
        args.case_dir,
        inlet_velocity=args.velocity,
        inlet_temp=args.temp,
        heat_flux=args.heat_flux,
        end_time=args.end_time
    )


if __name__ == "__main__":
    main()
