"""
Setup OpenFOAM case files for full staggered finned-tube heat exchanger.

This case uses real inlet/outlet boundary conditions (not cyclic) to enable
MPI parallelization and capture actual flow development through N REV units.
"""
import argparse
import re
from pathlib import Path
from textwrap import dedent

CASE_DIR = Path(__file__).resolve().parent.parent

# Geometric parameters
SL = 0.060  # longitudinal pitch
ST = 0.060  # transverse pitch
P_FIN = 0.050  # fin spacing
N_DEFAULT = 10  # default number of REV repetitions

# Domain dimensions (computed from N)
def get_domain_x(N):
    return N * 2 * SL


def detect_solid_regions(case_dir: Path) -> list:
    """Detect all solid regions from constant/ directory."""
    solid_regions = []
    constant_dir = case_dir / "constant"
    if constant_dir.exists():
        for item in constant_dir.iterdir():
            if item.is_dir() and item.name != "fluid" and (item / "polyMesh").exists():
                solid_regions.append(item.name)
    if not solid_regions:
        solid_regions = ["solid"]
    return sorted(solid_regions)


def detect_fluid_interface_patches(case_dir: Path) -> list:
    """Detect all fluid_to_* interface patches from fluid boundary file."""
    boundary_file = case_dir / "constant/fluid/polyMesh/boundary"
    interface_patches = []
    if boundary_file.exists():
        content = boundary_file.read_text()
        matches = re.findall(r'(fluid_to_\w+)', content)
        interface_patches = list(set(matches))
    if not interface_patches:
        interface_patches = ["fluid_to_solid"]
    return sorted(interface_patches)


def fix_fluid_boundary_types(case_dir: Path) -> None:
    """Fix boundary types after splitMeshRegions.
    
    For full domain (non-cyclic):
    - inlet: patch (remains as-is)
    - outlet: patch (remains as-is)
    - top, bottom, front, back: symmetry
    """
    boundary_file = case_dir / "constant/fluid/polyMesh/boundary"
    if not boundary_file.exists():
        return
    
    content = boundary_file.read_text()
    
    # Fix symmetry boundaries (NOT inlet/outlet - they stay as patch)
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


def control_dict(N: int) -> str:
    # More iterations for larger domain
    end_time = 500 + N * 50  # Scale with domain size
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }}

    application     chtMultiRegionSimpleFoam;
    startFrom       latestTime;
    startTime       0;
    stopAt          endTime;
    endTime         {end_time};
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
    {{
        inletMassFlow
        {{
            type            surfaceFieldValue;
            name            inlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       sum;
            surfaceFormat   raw;
            fields          (phi);
            writeFields     false;
            regionType      patch;
        }}

        outletMassFlow
        {{
            type            surfaceFieldValue;
            name            outlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       sum;
            surfaceFormat   raw;
            fields          (phi);
            writeFields     false;
            regionType      patch;
        }}

        inletTemperature
        {{
            type            surfaceFieldValue;
            name            inlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (T);
            writeFields     false;
            regionType      patch;
        }}

        outletTemperature
        {{
            type            surfaceFieldValue;
            name            outlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (T);
            writeFields     false;
            regionType      patch;
        }}

        pressureDrop
        {{
            type            pressure;
            name            pressureDrop;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            mode            staticCoeff;
            pRef            0;
            pInf            0;
            UInf            (1 0 0);
            rhoInf          1.18;
            executeControl  timeStep;
            executeInterval 10;
            writeControl    timeStep;
            writeInterval   100;
        }}
    }}
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
            tolerance       1e-7;
            relTol          0.01;
            smoother        DICGaussSeidel;
            cacheAgglomeration true;
            agglomerator    faceAreaPair;
            nCellsInCoarsestLevel 100;
            mergeLevels     1;
        }

        "(U|h|T)"
        {
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-7;
            relTol          0.1;
        }

        "(k|epsilon|omega)"
        {
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-7;
            relTol          0.1;
        }
    }

    SIMPLE
    {
        nNonOrthogonalCorrectors 1;
        pRefCell        0;
        pRefValue       0;
        residualControl
        {
            p_rgh           1e-4;
            U               1e-4;
            h               1e-4;
        }
    }

    relaxationFactors
    {
        equations
        {
            U               0.5;
            "(h|T)"         0.5;
            "(k|epsilon|omega)" 0.7;
        }
        fields
        {
            p_rgh           0.3;
        }
    }
    """


def decompose_par_dict(n_procs: int = 4) -> str:
    """Generate decomposeParDict for MPI parallelization.
    
    Use scotch method for automatic domain decomposition.
    """
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      decomposeParDict;
    }}

    numberOfSubdomains {n_procs};

    method          scotch;

    // Alternative: simple method with explicit decomposition
    // method          simple;
    // simpleCoeffs
    // {{
    //     n           ({n_procs} 1 1);  // Decompose along x (flow direction)
    //     delta       0.001;
    // }}
    """


def region_properties(solid_regions: list) -> str:
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
            molWeight   63.55;
        }}
        thermodynamics
        {{
            Cp          385;
            Hf          0;
        }}
        transport
        {{
            kappa       400;
        }}
        equationOfState
        {{
            rho         8960;
        }}
    }}
    """


def U_fluid(Ubar: float, interface_patches: list = None) -> str:
    """Velocity field with real inlet/outlet BCs."""
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]
    
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            noSlip;
        }}"""
    
    return dedent(f"""
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
                type            fixedValue;
                value           uniform ({Ubar} 0 0);
            }}
            outlet
            {{
                type            zeroGradient;
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
        """)


def p_rgh_fluid(interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]
    
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            fixedFluxPressure;
        }}"""
    
    return dedent(f"""
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
                type            zeroGradient;
            }}
            outlet
            {{
                type            fixedValue;
                value           uniform 0;
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
        """)


def p_fluid(interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]
    
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            calculated;
            value           uniform 0;
        }}"""
    
    return dedent(f"""
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
                type            calculated;
                value           uniform 0;
            }}
            outlet
            {{
                type            calculated;
                value           uniform 0;
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
        """)


def T_fluid(Tinlet: float, interface_patches: list = None) -> str:
    if interface_patches is None:
        interface_patches = ["fluid_to_solid"]
    
    interface_entries = ""
    for patch in interface_patches:
        interface_entries += f"""
        {patch}
        {{
            type            compressible::turbulentTemperatureCoupledBaffleMixed;
            Tnbr            T;
            kappaMethod     fluidThermo;
            value           uniform {Tinlet};
        }}"""
    
    return dedent(f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/fluid";
            object      T;
        }}

        dimensions      [0 0 0 1 0 0 0];
        internalField   uniform {Tinlet};

        boundaryField
        {{
            inlet
            {{
                type            fixedValue;
                value           uniform {Tinlet};
            }}
            outlet
            {{
                type            zeroGradient;
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
        """)


def T_solid(Tsolid: float, region_name: str = "solid") -> str:
    interface_patch = f"{region_name}_to_fluid"
    return dedent(f"""
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0/{region_name}";
            object      T;
        }}

        dimensions      [0 0 0 1 0 0 0];
        internalField   uniform {Tsolid};

        boundaryField
        {{
            {interface_patch}
            {{
                type            compressible::turbulentTemperatureCoupledBaffleMixed;
                Tnbr            T;
                kappaMethod     solidThermo;
                value           uniform {Tsolid};
            }}
            defaultFaces
            {{
                type            zeroGradient;
            }}
        }}
        """)


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


def create_case(case_dir: Path, Ubar: float, Tinlet: float, Tsolid: float, 
                N: int = N_DEFAULT, n_procs: int = 4) -> None:
    """Create OpenFOAM case files for full domain simulation."""
    solid_regions = detect_solid_regions(case_dir)
    interface_patches = detect_fluid_interface_patches(case_dir)
    
    print(f"Detected solid regions: {solid_regions}")
    print(f"Detected fluid interface patches: {interface_patches}")
    
    # System files
    write_file(case_dir / "system/controlDict", control_dict(N))
    write_file(case_dir / "system/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fvSolution", fv_solution())
    write_file(case_dir / "system/decomposeParDict", decompose_par_dict(n_procs))
    write_file(case_dir / "system/fluid/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fluid/fvSolution", fv_solution())
    
    for region in solid_regions:
        write_file(case_dir / f"system/{region}/fvSchemes", fv_schemes())
        write_file(case_dir / f"system/{region}/fvSolution", fv_solution())
    
    # Constant files
    write_file(case_dir / "constant/regionProperties", region_properties(solid_regions))
    write_file(case_dir / "constant/g", g_file())
    write_file(case_dir / "constant/fluid/turbulenceProperties", turbulence_properties_fluid())
    write_file(case_dir / "constant/fluid/thermophysicalProperties", thermophysical_fluid())
    
    for region in solid_regions:
        write_file(case_dir / f"constant/{region}/turbulenceProperties", 
                   turbulence_properties_solid(region))
        write_file(case_dir / f"constant/{region}/thermophysicalProperties", 
                   thermophysical_solid(region))
    
    # Initial conditions
    write_file(case_dir / "0/fluid/p", p_fluid(interface_patches))
    write_file(case_dir / "0/fluid/U", U_fluid(Ubar, interface_patches))
    write_file(case_dir / "0/fluid/p_rgh", p_rgh_fluid(interface_patches))
    write_file(case_dir / "0/fluid/T", T_fluid(Tinlet, interface_patches))
    
    for region in solid_regions:
        write_file(case_dir / f"0/{region}/p", p_solid(region))
        write_file(case_dir / f"0/{region}/T", T_solid(Tsolid, region))


def setup_after_split_mesh_regions(case_dir: Path, Ubar: float, Tinlet: float, 
                                    Tsolid: float, N: int = N_DEFAULT, 
                                    n_procs: int = 4) -> None:
    """Setup after splitMeshRegions - fix boundaries and regenerate files."""
    fix_fluid_boundary_types(case_dir)
    create_case(case_dir, Ubar, Tinlet, Tsolid, N, n_procs)
    print("Post-splitMeshRegions setup completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Setup OpenFOAM case for full staggered finned-tube heat exchanger."
    )
    parser.add_argument("--case", default=str(CASE_DIR), help="Case directory")
    parser.add_argument("--Ubar", type=float, default=1.0, help="Inlet velocity (m/s)")
    parser.add_argument("--Tinlet", type=float, default=300.0, help="Inlet temperature (K)")
    parser.add_argument("--Tsolid", type=float, default=350.0, help="Solid temperature (K)")
    parser.add_argument("-N", "--n-repeat", type=int, default=N_DEFAULT,
                        help="Number of REV repetitions")
    parser.add_argument("--nprocs", type=int, default=4, help="Number of MPI processes")
    parser.add_argument("--post-split", action="store_true",
                        help="Run setup after splitMeshRegions")
    args = parser.parse_args()
    
    case_path = Path(args.case)
    
    try:
        if args.post_split:
            setup_after_split_mesh_regions(
                case_path, args.Ubar, args.Tinlet, args.Tsolid, 
                args.n_repeat, args.nprocs
            )
        else:
            create_case(case_path, args.Ubar, args.Tinlet, args.Tsolid,
                       args.n_repeat, args.nprocs)
    except Exception as exc:
        raise SystemExit(f"Failed to setup case: {exc}")
    
    print(f"Case files written to {case_path}")

