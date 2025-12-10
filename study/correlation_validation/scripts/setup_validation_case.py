"""
Setup OpenFOAM case files for correlation validation.

이 스크립트는 상관식 검증용 CHT 케이스의 OpenFOAM 설정 파일을 생성합니다.
Re 수에 따라 입구 속도를 자동 계산합니다.

무차원수 정의:
- Re_d = ρ * V_max * d_o / μ
- V_max = V_fr / σ (σ = A_min/A_fr)
"""
import argparse
import math
import re
from pathlib import Path
from textwrap import dedent

# 형상 파라미터
D_O = 0.016       # 튜브 외경 [m]
S_T = 0.036       # 가로 피치 [m]
S_L = 0.034       # 세로 피치 [m]
H_FIN = 0.010     # 핀 방사 높이 [m]
T_FIN = 0.0005    # 핀 두께 [m]
S_FIN = 0.004     # clear fin spacing [m]
P_FIN = 0.0045    # 핀 피치 [m]
N_FINS = 5        # 핀 개수
N_REPEAT = 5      # REV 반복 횟수

# 확장 영역
INLET_EXT = 1.0 * S_L
OUTLET_EXT = 3.0 * S_L

# 유체 물성 (공기 300K)
RHO = 1.177       # 밀도 [kg/m³]
MU = 1.846e-5     # 점성 [Pa·s]
K_FLUID = 0.0263  # 열전도도 [W/m-K]
CP = 1005.0       # 비열 [J/kg-K]
PR = 0.707        # 프란틀 수

# 고체 물성 (알루미늄)
RHO_SOLID = 2700.0
K_SOLID = 200.0
CP_SOLID = 900.0


def compute_sigma():
    """
    최소 유동 면적비 σ = A_min/A_fr 계산.
    
    스태거드 배열에서 최소 유동 면적은 가로 갭과 대각선 갭 중 작은 것.
    핀의 영향을 고려하여 유효 간격을 계산합니다.
    """
    # 가로 갭 (튜브 사이)
    g_T = S_T - D_O
    
    # 대각선 거리
    S_d = math.sqrt((S_T / 2.0)**2 + S_L**2)
    g_D = S_d - D_O
    
    # 최소 갭
    g_min = min(g_T, g_D)
    
    # 면적비
    sigma = g_min / S_T
    
    return sigma


def compute_V_fr_from_Re(Re_d):
    """
    목표 Re_d에서 정면 속도 V_fr 계산.
    
    Re_d = ρ * V_max * d_o / μ
    V_max = V_fr / σ
    → V_fr = Re_d * μ * σ / (ρ * d_o)
    """
    sigma = compute_sigma()
    V_max = Re_d * MU / (RHO * D_O)
    V_fr = V_max * sigma
    return V_fr, V_max, sigma


def detect_solid_regions(case_dir: Path) -> list:
    """constant/ 디렉토리에서 솔리드 영역 감지"""
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
    """fluid boundary 파일에서 인터페이스 패치 감지"""
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
    """splitMeshRegions 후 경계 타입 수정"""
    boundary_file = case_dir / "constant/fluid/polyMesh/boundary"
    if not boundary_file.exists():
        return
    
    content = boundary_file.read_text()
    
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
    end_time = 500 + N * 50
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

        inletPressure
        {{
            type            surfaceFieldValue;
            name            inlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (p);
            writeFields     false;
            regionType      patch;
        }}

        outletPressure
        {{
            type            surfaceFieldValue;
            name            outlet;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            operation       areaAverage;
            surfaceFormat   raw;
            fields          (p);
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

        wallHeatFlux
        {{
            type            wallHeatFlux;
            region          fluid;
            libs            ("libfieldFunctionObjects.so");
            patches         ("fluid_to_.*");
            writeFields     true;
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


def fv_solution_fluid() -> str:
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
        nNonOrthogonalCorrectors 2;
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
            U               0.3;
            "(h|T)"         0.3;
            "(k|epsilon|omega)" 0.5;
        }
        fields
        {
            p_rgh           0.2;
        }
    }
    """


def fv_solution_solid() -> str:
    """고체 영역용 fvSolution (DIC preconditioner 사용)"""
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
        "(h|T)"
        {
            solver          PBiCGStab;
            preconditioner  DIC;
            tolerance       1e-7;
            relTol          0.1;
        }
    }

    SIMPLE
    {
        nNonOrthogonalCorrectors 1;
    }

    relaxationFactors
    {
        equations
        {
            "(h|T)"         0.5;
        }
    }
    """


def decompose_par_dict(n_procs: int = 4) -> str:
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
    return f"""
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant/fluid";
        object      thermophysicalProperties;
    }}

    thermoType
    {{
        type            heRhoThermo;
        mixture         pureMixture;
        transport       const;
        thermo          hConst;
        equationOfState rhoConst;
        specie          specie;
        energy          sensibleEnthalpy;
    }}

    mixture
    {{
        specie
        {{
            molWeight   28.96;
        }}
        thermodynamics
        {{
            Cp          {CP};
            Hf          0;
        }}
        transport
        {{
            mu          {MU};
            Pr          {PR};
            kappa       {K_FLUID};
        }}
        equationOfState
        {{
            rho         {RHO};
        }}
    }}
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
            molWeight   26.98;
        }}
        thermodynamics
        {{
            Cp          {CP_SOLID};
            Hf          0;
        }}
        transport
        {{
            kappa       {K_SOLID};
        }}
        equationOfState
        {{
            rho         {RHO_SOLID};
        }}
    }}
    """


def U_fluid(Ubar: float, interface_patches: list = None) -> str:
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


def create_case(case_dir: Path, Re_d: float, Tinlet: float, Tsolid: float,
                N: int = N_REPEAT, n_procs: int = 4) -> None:
    """OpenFOAM 케이스 파일 생성"""
    
    # Re에서 속도 계산
    V_fr, V_max, sigma = compute_V_fr_from_Re(Re_d)
    
    print(f"\n=== 케이스 설정 ===")
    print(f"목표 Re_d: {Re_d:.0f}")
    print(f"σ (A_min/A_fr): {sigma:.4f}")
    print(f"V_max: {V_max:.3f} m/s")
    print(f"V_fr (입구 속도): {V_fr:.3f} m/s")
    
    solid_regions = detect_solid_regions(case_dir)
    interface_patches = detect_fluid_interface_patches(case_dir)
    
    print(f"감지된 솔리드 영역: {solid_regions}")
    print(f"감지된 인터페이스 패치: {interface_patches}")
    
    # System files
    write_file(case_dir / "system/controlDict", control_dict(N))
    write_file(case_dir / "system/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fvSolution", fv_solution_fluid())
    write_file(case_dir / "system/decomposeParDict", decompose_par_dict(n_procs))
    write_file(case_dir / "system/fluid/fvSchemes", fv_schemes())
    write_file(case_dir / "system/fluid/fvSolution", fv_solution_fluid())
    
    for region in solid_regions:
        write_file(case_dir / f"system/{region}/fvSchemes", fv_schemes())
        write_file(case_dir / f"system/{region}/fvSolution", fv_solution_solid())
    
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
    write_file(case_dir / "0/fluid/U", U_fluid(V_fr, interface_patches))
    write_file(case_dir / "0/fluid/p_rgh", p_rgh_fluid(interface_patches))
    write_file(case_dir / "0/fluid/T", T_fluid(Tinlet, interface_patches))
    
    for region in solid_regions:
        write_file(case_dir / f"0/{region}/p", p_solid(region))
        write_file(case_dir / f"0/{region}/T", T_solid(Tsolid, region))
    
    # 케이스 정보 저장
    info_file = case_dir / "case_info.txt"
    with open(info_file, "w") as f:
        f.write(f"Re_d: {Re_d}\n")
        f.write(f"V_fr: {V_fr}\n")
        f.write(f"V_max: {V_max}\n")
        f.write(f"sigma: {sigma}\n")
        f.write(f"d_o: {D_O}\n")
        f.write(f"S_T: {S_T}\n")
        f.write(f"S_L: {S_L}\n")
        f.write(f"h_fin: {H_FIN}\n")
        f.write(f"t_fin: {T_FIN}\n")
        f.write(f"S_fin: {S_FIN}\n")
        f.write(f"rho: {RHO}\n")
        f.write(f"mu: {MU}\n")
        f.write(f"k_fluid: {K_FLUID}\n")
        f.write(f"Pr: {PR}\n")
        f.write(f"T_inlet: {Tinlet}\n")
        f.write(f"T_solid: {Tsolid}\n")


def setup_after_split_mesh_regions(case_dir: Path, Re_d: float, Tinlet: float,
                                    Tsolid: float, N: int = N_REPEAT,
                                    n_procs: int = 4) -> None:
    """splitMeshRegions 후 설정"""
    fix_fluid_boundary_types(case_dir)
    create_case(case_dir, Re_d, Tinlet, Tsolid, N, n_procs)
    print("Post-splitMeshRegions 설정 완료.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Setup OpenFOAM case for correlation validation."
    )
    parser.add_argument("--case", type=str, default=".",
                        help="Case directory")
    parser.add_argument("--Re", type=float, required=True,
                        help="Target Reynolds number (Re_d)")
    parser.add_argument("--Tinlet", type=float, default=300.0,
                        help="Inlet temperature (K)")
    parser.add_argument("--Tsolid", type=float, default=350.0,
                        help="Solid temperature (K)")
    parser.add_argument("-N", "--n-repeat", type=int, default=N_REPEAT,
                        help="Number of REV repetitions")
    parser.add_argument("--nprocs", type=int, default=4,
                        help="Number of MPI processes")
    parser.add_argument("--post-split", action="store_true",
                        help="Run setup after splitMeshRegions")
    args = parser.parse_args()
    
    case_path = Path(args.case)
    
    try:
        if args.post_split:
            setup_after_split_mesh_regions(
                case_path, args.Re, args.Tinlet, args.Tsolid,
                args.n_repeat, args.nprocs
            )
        else:
            create_case(case_path, args.Re, args.Tinlet, args.Tsolid,
                        args.n_repeat, args.nprocs)
    except Exception as exc:
        raise SystemExit(f"케이스 설정 실패: {exc}")
    
    print(f"\n케이스 파일 저장: {case_path}")

