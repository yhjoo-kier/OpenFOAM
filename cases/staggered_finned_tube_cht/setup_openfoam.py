import argparse
import os
from pathlib import Path
from textwrap import dedent

CASE_DIR = Path(__file__).resolve().parent

# Geometric parameters (tube pitches)
SL = 0.060  # longitudinal pitch
ST = 0.060  # transverse pitch
P_FIN = 0.050  # fin center-to-center spacing

# REV domain dimensions for staggered arrangement
# Domain: x in [0, 2*S_L], y in [0, S_T], z in [-P_fin/2, P_fin/2]
DOMAIN_X = 2 * SL      # streamwise length (covers 2 tube rows)
DOMAIN_Y = ST          # transverse height (full pitch, symmetry at y=0 and y=S_T)
DOMAIN_Z = P_FIN       # axial depth (fin spacing)


def write_file(path: Path, contents: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dedent(contents).lstrip("\n"), encoding="utf-8")


def control_dict() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }

    application     chtMultiRegionSimpleFoam;
    startFrom       latestTime;
    startTime       0;
    stopAt          endTime;
    endTime         1000;
    deltaT          1;
    writeControl    timeStep;
    writeInterval   100;
    purgeWrite      3;
    writeFormat     ascii;
    writePrecision  7;
    timeFormat      general;
    timePrecision   6;
    runTimeModifiable true;

    functions
    {
        inletPressureAverage
        {
            type            surfaceFieldValue;
            name            inlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (p_rgh);
            writeFields     false;
            regionType      patch;
        }

        outletPressureAverage
        {
            type            surfaceFieldValue;
            name            outlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (p_rgh);
            writeFields     false;
            regionType      patch;
        }

        inletTemperatureAverage
        {
            type            surfaceFieldValue;
            name            inlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (T);
            writeFields     false;
            regionType      patch;
        }

        outletTemperatureAverage
        {
            type            surfaceFieldValue;
            name            outlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (T);
            writeFields     false;
            regionType      patch;
        }

        interfaceTemperatureAverage
        {
            type            surfaceFieldValue;
            name            solid_to_fluid;
            region          solid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (T);
            writeFields     false;
            regionType      patch;
        }

        wallHeatFlux
        {
            type            wallHeatFlux;
            name            wallHeatFlux;
            libs            ("libfieldFunctionObjects.so");
            region          fluid;
            patches         (fluid_to_solid);
            executeControl  timeStep;
            executeInterval 1;
            writeControl    timeStep;
            writeInterval   100;
        }

        wallHeatFluxIntegral
        {
            type            surfaceFieldValue;
            name            fluid_to_solid;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaIntegrate;
            surfaceFormat   raw;
            fields          (wallHeatFlux);
            writeFields     false;
            regionType      patch;
        }
    }
    """


def fv_schemes() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
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
        grad(U)         Gauss linear;
    }

    divSchemes
    {
        default               none;
        div(phi,U)            Gauss linearUpwind grad(U);
        div(phi,k)            Gauss upwind;
        div(phi,epsilon)      Gauss upwind;
        div(phi,omega)        Gauss upwind;
        div(phi,K)            Gauss upwind;
        div(phi,h)            Gauss linearUpwind grad(h);
        div(phi,T)            Gauss linearUpwind grad(T);
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

    fluxRequired
    {
        default         no;
        p_rgh           ;
        pcorr           ;
    }
    """


def fv_solution() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      fvSolution;
    }

    solvers
    {
        p_rgh
        {
            solver          GAMG;
            tolerance       1e-8;
            relTol          0.05;
            smoother        DIC;
            cacheAgglomeration true;
            agglomerator    faceAreaPair;
            nCellsInCoarsestLevel 10;
            mergeLevels     1;
        }

        pcorr
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-8;
            relTol          0;
        }

        U
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }

        "(k|epsilon|omega)"
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }

        T
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }

        h
        {
            solver          smoothSolver;
            smoother        symGaussSeidel;
            tolerance       1e-8;
            relTol          0.1;
        }
    }

    SIMPLE
    {
        nNonOrthogonalCorrectors 0;
        pRefCell                 0;
        pRefValue                0;
    }

    relaxationFactors
    {
        equations
        {
            "(U|T|h)"            0.3;
            "(k|epsilon|omega)"  0.7;
            "(p|p_rgh)"          0.2;
        }
    }
    """


def fv_options(Ubar: float) -> str:
    return dedent(
        f"""
        /*
        Momentum and energy sources are injected via fvOptions to enforce a prescribed
        bulk velocity and stabilize the periodic thermal field without using fanPressure
        or fixedJump boundaries. meanVelocityForce drives the flow to a target mean
        velocity Ubar. The temperatureRelaxation source may be tuned (scaleMultiplier)
        to maintain a spatially uniform reference mean temperature when needed.
        */

        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       dictionary;
            location    "system";
            object      fvOptions;
        }}

        meanFlowDrive
        {{
            type            meanVelocityForce;
            active          yes;
            selectionMode   all;
            fields          (U);
            Ubar            ({Ubar} 0 0);
            relaxation      0.2;
        }}

        temperatureRelaxation
        {{
            type            scalarSemiImplicitSource;
            active          no; // set to yes to enforce a mean temperature in the fluid
            selectionMode   all;
            volumeMode      specific;
            injectionRateSuSp
            {{
                T           (0 0);
            }}
            // Use scaleMultiplier to tune the thermal source strength while keeping
            // periodic temperature boundary conditions on inlet/outlet.
            scaleMultiplier 0;
        }}
        """
    )


def region_properties() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      regionProperties;
    }

    regions
    (
        fluid (fluid)
        solid (solid)
    );
    """


def turbulence_properties_fluid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/fluid";
        object      turbulenceProperties;
    }

    simulationType  laminar;
    """


def turbulence_properties_solid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/solid";
        object      turbulenceProperties;
    }

    simulationType  laminar;
    """


def thermophysical_fluid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/fluid";
        object      thermophysicalProperties;
    }

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
            molWeight   28.96;
        }
        thermodynamics
        {
            Cp          1005;
            Hf          0;
        }
        transport
        {
            mu          1.8e-5;
            Pr          0.71;
            kappa       0.026;
        }
        equationOfState
        {
            rho         1.18;
        }
    }
    """


def thermophysical_solid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/solid";
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
            molWeight   63.55; // copper
        }
        thermodynamics
        {
            Cp          385;
            Hf          0;
        }
        transport
        {
            kappa       400; // W/m/K
        }
        equationOfState
        {
            rho         8960;
        }
    }
    """


def transport_properties_fluid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/fluid";
        object      transportProperties;
    }

    transportModel  Newtonian;
    nu              1.5e-5;
    """


def U_fluid(Ubar: float) -> str:
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volVectorField;
            location    "0/fluid";
            object      U;
        }}

        dimensions      [0 1 -1 0 0 0 0];
        internalField   uniform ({Ubar} 0 0);

        boundaryField
        {{
            inlet
            {{
                type            cyclicAMI;
                neighbourPatch  outlet;
                transform       translational;
                separationVector ({DOMAIN_X} 0 0);
            }}
            outlet
            {{
                type            cyclicAMI;
                neighbourPatch  inlet;
                transform       translational;
                separationVector (-{DOMAIN_X} 0 0);
            }}
            top
            {{
                type            symmetry;
            }}
            bottom
            {{
                type            symmetry;
            }}
            front
            {{
                type            symmetry;
            }}
            back
            {{
                type            symmetry;
            }}
            fluid_to_solid
            {{
                type            noSlip;
            }}
        }}
        """
    )


def p_rgh_fluid() -> str:
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/fluid";
            object      p_rgh;
        }}

        dimensions      [1 -1 -2 0 0 0 0];
        internalField   uniform 0;

        boundaryField
        {{
            inlet
            {{
                type            cyclicAMI;
                neighbourPatch  outlet;
                transform       translational;
                separationVector ({DOMAIN_X} 0 0);
            }}
            outlet
            {{
                type            cyclicAMI;
                neighbourPatch  inlet;
                transform       translational;
                separationVector (-{DOMAIN_X} 0 0);
            }}
            top
            {{
                type            symmetry;
            }}
            bottom
            {{
                type            symmetry;
            }}
            front
            {{
                type            symmetry;
            }}
            back
            {{
                type            symmetry;
            }}
            fluid_to_solid
            {{
                type            fixedFluxPressure;
            }}
        }}
        """
    )


def p_fluid() -> str:
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/fluid";
            object      p;
        }}

        dimensions      [1 -1 -2 0 0 0 0];
        internalField   uniform 0;

        boundaryField
        {{
            inlet
            {{
                type            cyclicAMI;
                neighbourPatch  outlet;
                transform       translational;
                separationVector ({DOMAIN_X} 0 0);
            }}
            outlet
            {{
                type            cyclicAMI;
                neighbourPatch  inlet;
                transform       translational;
                separationVector (-{DOMAIN_X} 0 0);
            }}
            top
            {{
                type            symmetry;
            }}
            bottom
            {{
                type            symmetry;
            }}
            front
            {{
                type            symmetry;
            }}
            back
            {{
                type            symmetry;
            }}
            fluid_to_solid
            {{
                type            fixedFluxPressure;
            }}
        }}
        """
    )


def T_fluid(initial_T: float) -> str:
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/fluid";
            object      T;
        }}

        dimensions      [0 0 0 1 0 0 0];
        internalField   uniform {initial_T};

        boundaryField
        {{
            inlet
            {{
                type            cyclicAMI;
                neighbourPatch  outlet;
                transform       translational;
                separationVector ({DOMAIN_X} 0 0);
            }}
            outlet
            {{
                type            cyclicAMI;
                neighbourPatch  inlet;
                transform       translational;
                separationVector (-{DOMAIN_X} 0 0);
            }}
            top
            {{
                type            symmetry;
            }}
            bottom
            {{
                type            symmetry;
            }}
            front
            {{
                type            symmetry;
            }}
            back
            {{
                type            symmetry;
            }}
            fluid_to_solid
            {{
                type            compressible::turbulentTemperatureCoupledBaffleMixed;
                Tnbr            T;
                kappaMethod     fluidThermo;
                value           uniform {initial_T};
            }}
        }}
        """
    )


def T_solid(initial_T: float) -> str:
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/solid";
            object      T;
        }}

        dimensions      [0 0 0 1 0 0 0];
        internalField   uniform {initial_T};

        boundaryField
        {{
            solid_to_fluid
            {{
                type            compressible::turbulentTemperatureCoupledBaffleMixed;
                Tnbr            T;
                kappaMethod     solidThermo;
                value           uniform {initial_T};
            }}
            defaultFaces
            {{
                type            zeroGradient;
            }}
        }}
        """
    )


def p_solid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0/solid";
        object      p;
    }

    dimensions      [1 -1 -2 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {
        solid_to_fluid
        {
            type            fixedValue;
            value           uniform 0;
        }
        defaultFaces
        {
            type            fixedValue;
            value           uniform 0;
        }
    }
    """


def U_solid() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volVectorField;
        location    "0/solid";
        object      U;
    }

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform (0 0 0);

    boundaryField
    {
        solid_to_fluid
        {
            type            fixedValue;
            value           uniform (0 0 0);
        }
        defaultFaces
        {
            type            fixedValue;
            value           uniform (0 0 0);
        }
    }
    """


def g_file() -> str:
    return """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       uniformDimensionedVectorField;
        location    "constant";
        object      g;
    }

    dimensions      [0 1 -2 0 0 0 0];
    value           (0 0 0);
    """


def create_case(case_dir: Path, Ubar: float, Tinlet: float, Tsolid: float) -> None:
    write_file(case_dir / "system/controlDict", control_dict())
    write_file(case_dir / "system/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fvSolution", fv_solution())
    write_file(case_dir / "system/fluid/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fluid/fvSolution", fv_solution())
    write_file(case_dir / "system/solid/fvSchemes", fv_schemes())
    write_file(case_dir / "system/solid/fvSolution", fv_solution())
    write_file(case_dir / "system/fvOptions", fv_options(Ubar))
    write_file(case_dir / "constant/regionProperties", region_properties())
    write_file(case_dir / "constant/g", g_file())

    write_file(case_dir / "constant/fluid/turbulenceProperties", turbulence_properties_fluid())
    write_file(case_dir / "constant/fluid/thermophysicalProperties", thermophysical_fluid())
    write_file(case_dir / "constant/fluid/transportProperties", transport_properties_fluid())

    write_file(case_dir / "constant/solid/turbulenceProperties", turbulence_properties_solid())
    write_file(case_dir / "constant/solid/thermophysicalProperties", thermophysical_solid())

    write_file(case_dir / "0/fluid/p", p_fluid())
    write_file(case_dir / "0/fluid/U", U_fluid(Ubar))
    write_file(case_dir / "0/fluid/p_rgh", p_rgh_fluid())
    write_file(case_dir / "0/fluid/T", T_fluid(Tinlet))

    write_file(case_dir / "0/solid/U", U_solid())
    write_file(case_dir / "0/solid/p", p_solid())
    write_file(case_dir / "0/solid/T", T_solid(Tsolid))

    write_file(
        case_dir / "system/decomposeParDict",
        """
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       dictionary;
            object      decomposeParDict;
        }
        numberOfSubdomains 4;
        method          scotch;
        """,
    )

    write_file(
        case_dir / "system/meshQualityDict",
        """
        MeshQualityControls
        {
            maxNonOrtho 70;
        }
        """,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a chtMultiRegionSimpleFoam case for the staggered finned-tube REV.",
    )
    parser.add_argument("--case", default=str(CASE_DIR), help="Case directory to populate")
    parser.add_argument("--Ubar", type=float, default=1.0, help="Target mean streamwise velocity (m/s)")
    parser.add_argument("--Tinlet", type=float, default=300.0, help="Initial/mean fluid temperature (K)")
    parser.add_argument("--Tsolid", type=float, default=320.0, help="Initial solid temperature (K)")
    args = parser.parse_args()

    try:
        create_case(Path(args.case), Ubar=args.Ubar, Tinlet=args.Tinlet, Tsolid=args.Tsolid)
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"Failed to write OpenFOAM dictionaries: {exc}")
    print(f"Case dictionaries written under {Path(args.case)}")
