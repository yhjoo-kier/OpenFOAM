import argparse
import os
import re
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


def detect_solid_regions(case_dir: Path) -> list:
    """Detect all solid regions from constant/ directory."""
    solid_regions = []
    constant_dir = case_dir / "constant"
    if constant_dir.exists():
        for item in constant_dir.iterdir():
            if item.is_dir() and item.name != "fluid" and (item / "polyMesh").exists():
                solid_regions.append(item.name)
    # Default if no regions detected yet
    if not solid_regions:
        solid_regions = ["solid"]
    return sorted(solid_regions)


def detect_fluid_interface_patches(case_dir: Path) -> list:
    """Detect all fluid_to_* interface patches from fluid boundary file."""
    boundary_file = case_dir / "constant/fluid/polyMesh/boundary"
    interface_patches = []
    if boundary_file.exists():
        content = boundary_file.read_text()
        # Find all fluid_to_* patches
        matches = re.findall(r'(fluid_to_\w+)', content)
        interface_patches = list(set(matches))
    # Default if boundary file doesn't exist
    if not interface_patches:
        interface_patches = ["fluid_to_solid"]
    return sorted(interface_patches)


def fix_fluid_boundary_types(case_dir: Path) -> None:
    """Fix boundary types after splitMeshRegions.

    splitMeshRegions converts cyclicAMI and symmetry boundaries to plain 'patch' type.
    This function restores the correct boundary types for the fluid region.
    """
    boundary_file = case_dir / "constant/fluid/polyMesh/boundary"
    if not boundary_file.exists():
        return

    content = boundary_file.read_text()

    # Fix inlet cyclicAMI
    content = re.sub(
        r'(inlet\s*\{[^}]*type\s+)patch;',
        rf'\1cyclicAMI;\n        neighbourPatch  outlet;\n        transform       translational;\n        separationVector ({DOMAIN_X} 0 0);',
        content
    )

    # Fix outlet cyclicAMI
    content = re.sub(
        r'(outlet\s*\{[^}]*type\s+)patch;',
        rf'\1cyclicAMI;\n        neighbourPatch  inlet;\n        transform       translational;\n        separationVector (-{DOMAIN_X} 0 0);',
        content
    )

    # Fix symmetry boundaries
    for patch_name in ['top', 'bottom', 'front', 'back']:
        content = re.sub(
            rf'({patch_name}\s*\{{[^}}]*type\s+)patch;',
            rf'\1symmetry;',
            content
        )

    boundary_file.write_text(content)


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


def region_properties(solid_regions: list) -> str:
    """Generate regionProperties with all detected solid regions."""
    solid_list = " ".join(solid_regions)
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      regionProperties;
    }}

    regions
    (
        fluid (fluid)
        solid ({solid_list})
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


def turbulence_properties_solid(region_name: str = "solid") -> str:
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/{region_name}";
        object      turbulenceProperties;
    }}

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


def thermophysical_solid(region_name: str = "solid") -> str:
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/{region_name}";
        object      thermophysicalProperties;
    }}

    thermoType
    {{
        type            heSolidThermo;
        mixture         pureMixture;
        transport       constIso;
        thermo          hConst;
        equationOfState rhoConst;
        specie          specie;
        energy          sensibleEnthalpy;
    }}

    mixture
    {{
        specie
        {{
            molWeight   63.55; // copper
        }}
        thermodynamics
        {{
            Cp          385;
            Hf          0;
        }}
        transport
        {{
            kappa       400; // W/m/K
        }}
        equationOfState
        {{
            rho         8960;
        }}
    }}
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


def U_fluid(Ubar: float, interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]

    # Generate interface patch entries
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            noSlip;
        }}"""

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
            }}{interface_entries}
        }}
        """
    )


def p_rgh_fluid(interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]

    # Generate interface patch entries
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            fixedFluxPressure;
        }}"""

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
            }}{interface_entries}
        }}
        """
    )


def p_fluid(interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]

    # Generate interface patch entries
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            fixedFluxPressure;
        }}"""

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
            }}{interface_entries}
        }}
        """
    )


def T_fluid(initial_T: float, interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]

    # Generate interface patch entries for CHT coupling
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            compressible::turbulentTemperatureCoupledBaffleMixed;
            Tnbr            T;
            kappaMethod     fluidThermo;
            value           uniform {initial_T};
        }}"""

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
            }}{interface_entries}
        }}
        """
    )


def T_solid(initial_T: float, region_name: str = "solid") -> str:
    # Generate interface patch name (e.g., solid_to_fluid, domain0_to_fluid)
    interface_patch = f"{region_name}_to_fluid"
    return dedent(
        f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/{region_name}";
            object      T;
        }}

        dimensions      [0 0 0 1 0 0 0];
        internalField   uniform {initial_T};

        boundaryField
        {{
            {interface_patch}
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


def p_solid(region_name: str = "solid") -> str:
    interface_patch = f"{region_name}_to_fluid"
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0/{region_name}";
        object      p;
    }}

    dimensions      [1 -1 -2 0 0 0 0];
    internalField   uniform 0;

    boundaryField
    {{
        {interface_patch}
        {{
            type            fixedValue;
            value           uniform 0;
        }}
        defaultFaces
        {{
            type            fixedValue;
            value           uniform 0;
        }}
    }}
    """


def U_solid(region_name: str = "solid") -> str:
    interface_patch = f"{region_name}_to_fluid"
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       volVectorField;
        location    "0/{region_name}";
        object      U;
    }}

    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform (0 0 0);

    boundaryField
    {{
        {interface_patch}
        {{
            type            fixedValue;
            value           uniform (0 0 0);
        }}
        defaultFaces
        {{
            type            fixedValue;
            value           uniform (0 0 0);
        }}
    }}
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
    """Create OpenFOAM case files.

    This function can be called in two modes:
    1. Before mesh generation: Creates initial files with default solid region
    2. After splitMeshRegions: Detects and sets up all solid regions
    """
    # Detect solid regions (returns ["solid"] if mesh not yet split)
    solid_regions = detect_solid_regions(case_dir)

    # Detect fluid interface patches (returns ["fluid_to_solid"] if boundary file doesn't exist)
    interface_patches = detect_fluid_interface_patches(case_dir)

    print(f"Detected solid regions: {solid_regions}")
    print(f"Detected fluid interface patches: {interface_patches}")

    # Write system files
    write_file(case_dir / "system/controlDict", control_dict())
    write_file(case_dir / "system/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fvSolution", fv_solution())
    write_file(case_dir / "system/fluid/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fluid/fvSolution", fv_solution())
    write_file(case_dir / "system/fvOptions", fv_options(Ubar))

    # Write system files for all solid regions
    for region in solid_regions:
        write_file(case_dir / f"system/{region}/fvSchemes", fv_schemes())
        write_file(case_dir / f"system/{region}/fvSolution", fv_solution())

    # Write constant files
    write_file(case_dir / "constant/regionProperties", region_properties(solid_regions))
    write_file(case_dir / "constant/g", g_file())

    # Fluid constant properties
    write_file(case_dir / "constant/fluid/turbulenceProperties", turbulence_properties_fluid())
    write_file(case_dir / "constant/fluid/thermophysicalProperties", thermophysical_fluid())
    write_file(case_dir / "constant/fluid/transportProperties", transport_properties_fluid())

    # Solid constant properties for all detected solid regions
    for region in solid_regions:
        write_file(case_dir / f"constant/{region}/turbulenceProperties", turbulence_properties_solid(region))
        write_file(case_dir / f"constant/{region}/thermophysicalProperties", thermophysical_solid(region))

    # Fluid initial conditions with all interface patches
    write_file(case_dir / "0/fluid/p", p_fluid(interface_patches))
    write_file(case_dir / "0/fluid/U", U_fluid(Ubar, interface_patches))
    write_file(case_dir / "0/fluid/p_rgh", p_rgh_fluid(interface_patches))
    write_file(case_dir / "0/fluid/T", T_fluid(Tinlet, interface_patches))

    # Solid initial conditions for all detected solid regions
    for region in solid_regions:
        write_file(case_dir / f"0/{region}/U", U_solid(region))
        write_file(case_dir / f"0/{region}/p", p_solid(region))
        write_file(case_dir / f"0/{region}/T", T_solid(Tsolid, region))

    # Use simple decomposition method (more compatible with cyclicAMI + multi-region)
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
        method          simple;
        simpleCoeffs
        {
            n           (2 2 1);
            delta       0.001;
        }
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


def setup_after_split_mesh_regions(case_dir: Path, Ubar: float, Tinlet: float, Tsolid: float) -> None:
    """Call this function after running splitMeshRegions to fix boundaries and regenerate files.

    This function:
    1. Fixes boundary types (cyclicAMI, symmetry) that splitMeshRegions converts to 'patch'
    2. Detects all solid regions created by splitMeshRegions
    3. Regenerates all boundary condition files with correct interface patches
    """
    # Fix fluid boundary types first
    fix_fluid_boundary_types(case_dir)

    # Regenerate all case files with detected regions
    create_case(case_dir, Ubar, Tinlet, Tsolid)

    print("Post-splitMeshRegions setup completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a chtMultiRegionSimpleFoam case for the staggered finned-tube REV.",
    )
    parser.add_argument("--case", default=str(CASE_DIR), help="Case directory to populate")
    parser.add_argument("--Ubar", type=float, default=1.0, help="Target mean streamwise velocity (m/s)")
    parser.add_argument("--Tinlet", type=float, default=300.0, help="Initial/mean fluid temperature (K)")
    parser.add_argument("--Tsolid", type=float, default=320.0, help="Initial solid temperature (K)")
    parser.add_argument(
        "--post-split",
        action="store_true",
        help="Run setup after splitMeshRegions (fixes boundaries and regenerates files)",
    )
    args = parser.parse_args()

    try:
        if args.post_split:
            setup_after_split_mesh_regions(
                Path(args.case), Ubar=args.Ubar, Tinlet=args.Tinlet, Tsolid=args.Tsolid
            )
        else:
            create_case(Path(args.case), Ubar=args.Ubar, Tinlet=args.Tinlet, Tsolid=args.Tsolid)
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"Failed to write OpenFOAM dictionaries: {exc}")
    print(f"Case dictionaries written under {Path(args.case)}")
