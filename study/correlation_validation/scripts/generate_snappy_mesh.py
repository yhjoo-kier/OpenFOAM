"""
Generate mesh using snappyHexMesh for correlation validation cases.

이 스크립트는 snappyHexMesh를 사용하여 6면체 우세 메시를 생성합니다.
Gmsh 테트라헤드럴 메시보다 CHT 시뮬레이션에 더 안정적입니다.

절차:
1. Gmsh로 STL 형상 파일 생성
2. blockMesh로 배경 메시 생성
3. snappyHexMesh로 형상 메시 생성
4. splitMeshRegions로 영역 분리
"""
import argparse
import os
import subprocess
import math

# 검증용 형상 파라미터 (상관식 적용 범위 중앙값)
PARAMS = {
    "SL": 0.034,           # 34 mm, longitudinal pitch
    "ST": 0.036,           # 36 mm, transverse pitch
    "d_o": 0.016,          # 16 mm, tube outer diameter
    "d_i": 0.012,          # 12 mm, tube inner diameter
    "d_f": 0.036,          # 36 mm, fin outer diameter (h_fin = 10mm)
    "t_fin": 0.0005,       # 0.5 mm, fin thickness
    "P_fin": 0.0045,       # 4.5 mm, fin pitch
    "N_fins": 3,           # 핀 개수 (z방향)
    "N_repeat": 2,         # REV 반복 횟수 (x방향)
    "inlet_ext": 1.0,      # inlet extension in S_L
    "outlet_ext": 2.0,     # outlet extension in S_L
}


def create_stl_geometry(case_dir: str, params: dict) -> str:
    """
    Gmsh를 사용하여 STL 형상 파일 생성.
    유체-고체 인터페이스만 STL로 출력.
    """
    import gmsh
    
    SL = params["SL"]
    ST = params["ST"]
    d_o = params["d_o"]
    d_i = params["d_i"]
    d_f = params["d_f"]
    t_fin = params["t_fin"]
    P_fin = params["P_fin"]
    N_fins = params["N_fins"]
    N = params["N_repeat"]
    
    inlet_ext = params["inlet_ext"] * SL
    outlet_ext = params["outlet_ext"] * SL
    
    REV_length = 2.0 * SL
    core_length = N * REV_length
    domain_x = inlet_ext + core_length + outlet_ext
    domain_y = ST
    domain_z = N_fins * P_fin
    
    R_tube_out = d_o / 2.0
    R_tube_in = d_i / 2.0
    R_fin = d_f / 2.0
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("finned_tube_stl")
    
    occ = gmsh.model.occ
    
    # 모든 솔리드 유닛 생성 (튜브 + 핀)
    all_solids = []
    
    for i in range(N):
        x_offset = inlet_ext + i * REV_length
        
        tube_centers = [
            (x_offset + SL / 2.0, 0.0),
            (x_offset + SL / 2.0, ST),
            (x_offset + 3.0 * SL / 2.0, ST / 2.0),
        ]
        
        for cx, cy in tube_centers:
            # 튜브 외부
            tube_out = occ.addCylinder(cx, cy, 0.0, 0.0, 0.0, domain_z, R_tube_out)
            
            # 핀들
            fin_tags = []
            for j in range(N_fins):
                z_center = (j + 0.5) * P_fin
                z_start = z_center - t_fin / 2.0
                fin = occ.addCylinder(cx, cy, z_start, 0.0, 0.0, t_fin, R_fin)
                fin_tags.append((3, fin))
            
            # 내부 보어
            tube_in = occ.addCylinder(cx, cy, 0.0, 0.0, 0.0, domain_z, R_tube_in)
            
            # 융합
            if fin_tags:
                solid_unit, _ = occ.fuse([(3, tube_out)], fin_tags)
            else:
                solid_unit = [(3, tube_out)]
            solid_tag = solid_unit[0][1]
            
            # 내부 보어 제거
            solid_cut, _ = occ.cut([(3, solid_tag)], [(3, tube_in)])
            solid_tag = solid_cut[0][1]
            
            # 도메인으로 트리밍
            domain_copy = occ.addBox(0.0, 0.0, 0.0, domain_x, domain_y, domain_z)
            trimmed, _ = occ.intersect([(3, solid_tag)], [(3, domain_copy)],
                                       removeObject=True, removeTool=True)
            if trimmed:
                all_solids.extend(trimmed)
    
    occ.synchronize()
    
    # 메시 설정 (STL 출력용 - 표면만)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.0005)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.002)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 24)
    
    # 표면 메시 생성
    gmsh.model.mesh.generate(2)
    
    # STL 저장
    stl_file = os.path.join(case_dir, "constant/triSurface/finnedTube.stl")
    os.makedirs(os.path.dirname(stl_file), exist_ok=True)
    gmsh.write(stl_file)
    
    gmsh.finalize()
    
    print(f"STL 생성: {stl_file}")
    return stl_file, domain_x, domain_y, domain_z


def create_block_mesh_dict(case_dir: str, domain_x: float, domain_y: float, 
                           domain_z: float, cell_size: float = 0.002) -> None:
    """blockMeshDict 생성 (배경 메시)"""
    
    nx = int(math.ceil(domain_x / cell_size))
    ny = int(math.ceil(domain_y / cell_size))
    nz = int(math.ceil(domain_z / cell_size))
    
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale 1;

vertices
(
    (0 0 0)
    ({domain_x} 0 0)
    ({domain_x} {domain_y} 0)
    (0 {domain_y} 0)
    (0 0 {domain_z})
    ({domain_x} 0 {domain_z})
    ({domain_x} {domain_y} {domain_z})
    (0 {domain_y} {domain_z})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
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
    top
    {{
        type symmetry;
        faces
        (
            (3 7 6 2)
        );
    }}
    bottom
    {{
        type symmetry;
        faces
        (
            (0 1 5 4)
        );
    }}
    front
    {{
        type symmetry;
        faces
        (
            (4 5 6 7)
        );
    }}
    back
    {{
        type symmetry;
        faces
        (
            (0 3 2 1)
        );
    }}
);

mergePatchPairs
(
);
"""
    
    dict_path = os.path.join(case_dir, "system/blockMeshDict")
    os.makedirs(os.path.dirname(dict_path), exist_ok=True)
    with open(dict_path, "w") as f:
        f.write(content)
    print(f"blockMeshDict 생성: {dict_path}")


def create_snappy_hex_mesh_dict(case_dir: str, params: dict) -> None:
    """snappyHexMeshDict 생성"""
    
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
    finnedTube.stl
    {
        type triSurfaceMesh;
        name finnedTube;
        regions
        {
            finnedTube
            {
                name finnedTube;
            }
        }
    }
}

castellatedMeshControls
{
    maxLocalCells       500000;
    maxGlobalCells      2000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 3;

    features
    (
    );

    refinementSurfaces
    {
        finnedTube
        {
            level (2 3);
            cellZone solid;
            faceZone solidFaces;
            cellZoneInside inside;
        }
    }

    resolveFeatureAngle 30;

    refinementRegions
    {
    }

    locationInMesh (0.001 0.018 0.005);  // 유체 영역 내부 점
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch    3;
    tolerance       2.0;
    nSolveIter      100;
    nRelaxIter      5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap true;
}

addLayersControls
{
    relativeSizes   true;
    layers
    {
    }
    expansionRatio  1.0;
    finalLayerThickness 0.3;
    minThickness    0.1;
    nGrow           0;
    featureAngle    60;
    nRelaxIter      3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals  3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter      50;
}

meshQualityControls
{
    maxNonOrtho     65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave      80;
    minVol          1e-13;
    minTetQuality   -1e30;
    minArea         -1;
    minTwist        0.01;
    minDeterminant  0.001;
    minFaceWeight   0.05;
    minVolRatio     0.01;
    minTriangleTwist -1;
    nSmoothScale    4;
    errorReduction  0.75;
    relaxed
    {
        maxNonOrtho 75;
    }
}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);

mergeTolerance 1e-6;
"""
    
    dict_path = os.path.join(case_dir, "system/snappyHexMeshDict")
    with open(dict_path, "w") as f:
        f.write(content)
    print(f"snappyHexMeshDict 생성: {dict_path}")


def create_mesh_quality_dict(case_dir: str) -> None:
    """meshQualityDict 생성"""
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      meshQualityDict;
}

maxNonOrtho     65;
maxBoundarySkewness 20;
maxInternalSkewness 4;
maxConcave      80;
minVol          1e-13;
minTetQuality   -1e30;
minArea         -1;
minTwist        0.01;
minDeterminant  0.001;
minFaceWeight   0.05;
minVolRatio     0.01;
minTriangleTwist -1;
nSmoothScale    4;
errorReduction  0.75;
"""
    dict_path = os.path.join(case_dir, "system/meshQualityDict")
    with open(dict_path, "w") as f:
        f.write(content)


def create_surface_feature_extract_dict(case_dir: str) -> None:
    """surfaceFeatureExtractDict 생성"""
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}

finnedTube.stl
{
    extractionMethod    extractFromSurface;
    extractFromSurfaceCoeffs
    {
        includedAngle   150;
    }
    writeObj            yes;
}
"""
    dict_path = os.path.join(case_dir, "system/surfaceFeatureExtractDict")
    with open(dict_path, "w") as f:
        f.write(content)


def create_decompose_par_dict(case_dir: str, n_procs: int = 4) -> None:
    """decomposeParDict 생성"""
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}

numberOfSubdomains {n_procs};
method          scotch;
"""
    dict_path = os.path.join(case_dir, "system/decomposeParDict")
    with open(dict_path, "w") as f:
        f.write(content)


def create_control_dict(case_dir: str) -> None:
    """controlDict 생성"""
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     snappyHexMesh;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         100;
deltaT          1;
writeControl    timeStep;
writeInterval   100;
"""
    dict_path = os.path.join(case_dir, "system/controlDict")
    os.makedirs(os.path.dirname(dict_path), exist_ok=True)
    with open(dict_path, "w") as f:
        f.write(content)


def create_fv_schemes(case_dir: str) -> None:
    """fvSchemes 생성"""
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
    dict_path = os.path.join(case_dir, "system/fvSchemes")
    with open(dict_path, "w") as f:
        f.write(content)


def create_fv_solution(case_dir: str) -> None:
    """fvSolution 생성"""
    content = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers {}
"""
    dict_path = os.path.join(case_dir, "system/fvSolution")
    with open(dict_path, "w") as f:
        f.write(content)


def run_mesh_generation(case_dir: str) -> None:
    """메시 생성 명령 실행"""
    
    commands = [
        "blockMesh",
        "surfaceFeatureExtract",
        "snappyHexMesh -overwrite",
    ]
    
    for cmd in commands:
        print(f"\n실행: {cmd}")
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=case_dir,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            print(f"오류: {result.stderr[-500:] if result.stderr else 'Unknown error'}")
            raise RuntimeError(f"{cmd} 실패")
        print(f"완료: {cmd}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate mesh using snappyHexMesh for validation."
    )
    parser.add_argument("--case", type=str, required=True,
                        help="Case directory name")
    parser.add_argument("-N", "--n-repeat", type=int, default=PARAMS["N_repeat"],
                        help="Number of REV repetitions")
    parser.add_argument("--n-fins", type=int, default=PARAMS["N_fins"],
                        help="Number of fins")
    parser.add_argument("--run", action="store_true",
                        help="Run mesh generation commands")
    args = parser.parse_args()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    study_dir = os.path.dirname(script_dir)
    case_dir = os.path.join(study_dir, args.case)
    
    os.makedirs(case_dir, exist_ok=True)
    
    params = dict(PARAMS)
    params["N_repeat"] = args.n_repeat
    params["N_fins"] = args.n_fins
    
    # 1. 시스템 파일 생성
    create_control_dict(case_dir)
    create_fv_schemes(case_dir)
    create_fv_solution(case_dir)
    create_decompose_par_dict(case_dir)
    create_mesh_quality_dict(case_dir)
    create_surface_feature_extract_dict(case_dir)
    
    # 2. STL 생성
    stl_file, domain_x, domain_y, domain_z = create_stl_geometry(case_dir, params)
    
    # 3. blockMeshDict 생성
    create_block_mesh_dict(case_dir, domain_x, domain_y, domain_z)
    
    # 4. snappyHexMeshDict 생성
    create_snappy_hex_mesh_dict(case_dir, params)
    
    print(f"\n모든 설정 파일 생성 완료: {case_dir}")
    print(f"도메인 크기: {domain_x*1000:.1f} x {domain_y*1000:.1f} x {domain_z*1000:.1f} mm")
    
    if args.run:
        run_mesh_generation(case_dir)
    else:
        print("\n메시 생성 명령:")
        print(f"  cd {case_dir}")
        print("  blockMesh")
        print("  surfaceFeatureExtract")
        print("  snappyHexMesh -overwrite")
        print("  splitMeshRegions -cellZones -overwrite")


if __name__ == "__main__":
    main()

