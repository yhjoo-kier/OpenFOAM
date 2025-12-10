"""
Generate STL geometry for snappyHexMesh.

스태거드 핀-튜브 형상의 STL 파일을 생성합니다.
"""
import numpy as np
import os

# 형상 파라미터
D_O = 0.016       # 튜브 외경 [m]
R_TUBE_OUT = D_O / 2
R_TUBE_IN = 0.006  # 내경
R_FIN = 0.018     # 핀 외경
H_FIN = R_FIN - R_TUBE_OUT  # 10 mm
T_FIN = 0.0005    # 핀 두께
S_FIN = 0.004     # clear fin spacing
P_FIN = T_FIN + S_FIN  # 핀 피치 4.5 mm

S_T = 0.036       # 가로 피치
S_L = 0.034       # 세로 피치

N_FINS = 3        # 핀 개수
N_REPEAT = 2      # REV 반복

INLET_EXT = 1.0 * S_L
OUTLET_EXT = 2.0 * S_L


def create_cylinder_stl(cx, cy, z_start, z_end, radius, n_segments=32):
    """원통 STL 삼각형 생성"""
    triangles = []
    
    # 측면
    for i in range(n_segments):
        theta1 = 2 * np.pi * i / n_segments
        theta2 = 2 * np.pi * (i + 1) / n_segments
        
        x1, y1 = cx + radius * np.cos(theta1), cy + radius * np.sin(theta1)
        x2, y2 = cx + radius * np.cos(theta2), cy + radius * np.sin(theta2)
        
        # 두 개의 삼각형으로 사각형 형성
        triangles.append([(x1, y1, z_start), (x2, y2, z_start), (x1, y1, z_end)])
        triangles.append([(x2, y2, z_start), (x2, y2, z_end), (x1, y1, z_end)])
    
    return triangles


def create_disk_stl(cx, cy, z, r_inner, r_outer, n_segments=32):
    """도넛형 디스크 STL 삼각형 생성"""
    triangles = []
    
    for i in range(n_segments):
        theta1 = 2 * np.pi * i / n_segments
        theta2 = 2 * np.pi * (i + 1) / n_segments
        
        # 내부 점
        xi1, yi1 = cx + r_inner * np.cos(theta1), cy + r_inner * np.sin(theta1)
        xi2, yi2 = cx + r_inner * np.cos(theta2), cy + r_inner * np.sin(theta2)
        
        # 외부 점
        xo1, yo1 = cx + r_outer * np.cos(theta1), cy + r_outer * np.sin(theta1)
        xo2, yo2 = cx + r_outer * np.cos(theta2), cy + r_outer * np.sin(theta2)
        
        # 두 개의 삼각형
        triangles.append([(xi1, yi1, z), (xo1, yo1, z), (xi2, yi2, z)])
        triangles.append([(xi2, yi2, z), (xo1, yo1, z), (xo2, yo2, z)])
    
    return triangles


def write_stl(filename, triangles, solid_name="solid"):
    """STL 파일 쓰기 (ASCII)"""
    with open(filename, 'w') as f:
        f.write(f"solid {solid_name}\n")
        for tri in triangles:
            p1, p2, p3 = tri
            # 법선 벡터 계산
            v1 = np.array(p2) - np.array(p1)
            v2 = np.array(p3) - np.array(p1)
            normal = np.cross(v1, v2)
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm
            else:
                normal = np.array([0, 0, 1])
            
            f.write(f"  facet normal {normal[0]:.6e} {normal[1]:.6e} {normal[2]:.6e}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {p1[0]:.6e} {p1[1]:.6e} {p1[2]:.6e}\n")
            f.write(f"      vertex {p2[0]:.6e} {p2[1]:.6e} {p2[2]:.6e}\n")
            f.write(f"      vertex {p3[0]:.6e} {p3[1]:.6e} {p3[2]:.6e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")


def generate_finned_tube_stl(output_dir):
    """핀-튜브 형상 STL 생성"""
    
    domain_z = N_FINS * P_FIN
    REV_length = 2 * S_L
    core_length = N_REPEAT * REV_length
    domain_x = INLET_EXT + core_length + OUTLET_EXT
    domain_y = S_T
    
    print(f"Domain: {domain_x*1000:.1f} x {domain_y*1000:.1f} x {domain_z*1000:.1f} mm")
    
    all_triangles = []
    
    # 튜브 중심 위치 수집
    tube_centers = []
    for i in range(N_REPEAT):
        x_offset = INLET_EXT + i * REV_length
        
        # 스태거드 배열
        centers = [
            (x_offset + S_L / 2.0, 0.0),           # Row 1: y=0
            (x_offset + S_L / 2.0, S_T),           # Row 1: y=S_T
            (x_offset + 3.0 * S_L / 2.0, S_T / 2.0),  # Row 2: 중앙
        ]
        tube_centers.extend(centers)
    
    # 각 튜브에 대해 형상 생성
    for cx, cy in tube_centers:
        # 튜브 외부 측면
        tube_tris = create_cylinder_stl(cx, cy, 0, domain_z, R_TUBE_OUT, n_segments=24)
        all_triangles.extend(tube_tris)
        
        # 핀들
        for j in range(N_FINS):
            z_center = (j + 0.5) * P_FIN
            z_start = z_center - T_FIN / 2
            z_end = z_center + T_FIN / 2
            
            # 핀 측면 (외부)
            fin_side = create_cylinder_stl(cx, cy, z_start, z_end, R_FIN, n_segments=24)
            all_triangles.extend(fin_side)
            
            # 핀 상면
            fin_top = create_disk_stl(cx, cy, z_end, R_TUBE_OUT, R_FIN, n_segments=24)
            all_triangles.extend(fin_top)
            
            # 핀 하면 (법선 반대)
            fin_bottom = create_disk_stl(cx, cy, z_start, R_TUBE_OUT, R_FIN, n_segments=24)
            # 법선 방향 반전을 위해 점 순서 변경
            fin_bottom = [(t[0], t[2], t[1]) for t in fin_bottom]
            all_triangles.extend(fin_bottom)
    
    # STL 저장
    stl_file = os.path.join(output_dir, "finnedTube.stl")
    write_stl(stl_file, all_triangles, "finnedTube")
    print(f"STL 저장: {stl_file}")
    print(f"삼각형 수: {len(all_triangles)}")
    
    return domain_x, domain_y, domain_z


def generate_blockmesh_dict(output_dir, domain_x, domain_y, domain_z):
    """blockMeshDict 생성"""
    
    # 셀 수 계산 (약 5mm 셀 크기)
    cell_size = 0.003
    nx = max(10, int(domain_x / cell_size))
    ny = max(10, int(domain_y / cell_size))
    nz = max(5, int(domain_z / cell_size))
    
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
    
    filepath = os.path.join(output_dir, "system", "blockMeshDict")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"blockMeshDict 저장: {filepath}")


def generate_snappyhexmesh_dict(output_dir, domain_x, domain_y, domain_z):
    """snappyHexMeshDict 생성"""
    
    # 정제 영역 (입구 확장 후부터 출구 확장 전까지)
    refine_x_min = INLET_EXT * 0.5
    refine_x_max = domain_x - OUTLET_EXT * 0.5
    
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    finnedTube.stl
    {{
        type triSurfaceMesh;
        name finnedTube;
    }}
    
    refinementBox
    {{
        type searchableBox;
        min ({refine_x_min} 0 0);
        max ({refine_x_max} {domain_y} {domain_z});
    }}
}}

castellatedMeshControls
{{
    maxLocalCells 1000000;
    maxGlobalCells 2000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;
    
    features
    (
    );
    
    refinementSurfaces
    {{
        finnedTube
        {{
            level (2 3);
            
            cellZone solid;
            faceZone solid;
            cellZoneInside inside;
        }}
    }}
    
    resolveFeatureAngle 30;
    
    refinementRegions
    {{
        refinementBox
        {{
            mode inside;
            levels ((1E15 1));
        }}
    }}
    
    locationInMesh ({INLET_EXT/2} {domain_y/2} {domain_z/2});
    
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true;
    layers
    {{
    }}
    expansionRatio 1.2;
    finalLayerThickness 0.5;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality -1e30;
    minArea -1;
    minTwist 0.01;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
    relaxed
    {{
        maxNonOrtho 75;
    }}
}}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);

mergeTolerance 1e-6;
"""
    
    filepath = os.path.join(output_dir, "system", "snappyHexMeshDict")
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"snappyHexMeshDict 저장: {filepath}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate STL and snappyHexMesh files")
    parser.add_argument("--case", type=str, required=True, help="Case directory")
    args = parser.parse_args()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    study_dir = os.path.dirname(script_dir)
    case_dir = os.path.join(study_dir, args.case)
    
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant", "triSurface"), exist_ok=True)
    
    # STL 생성
    stl_dir = os.path.join(case_dir, "constant", "triSurface")
    domain_x, domain_y, domain_z = generate_finned_tube_stl(stl_dir)
    
    # blockMeshDict 생성
    generate_blockmesh_dict(case_dir, domain_x, domain_y, domain_z)
    
    # snappyHexMeshDict 생성
    generate_snappyhexmesh_dict(case_dir, domain_x, domain_y, domain_z)
    
    print("\n메시 생성 명령어:")
    print(f"  cd {case_dir}")
    print("  blockMesh")
    print("  snappyHexMesh -overwrite")


if __name__ == "__main__":
    main()

