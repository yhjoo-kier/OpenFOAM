"""
Generate mesh for correlation validation cases.

이 스크립트는 Briggs & Young (1963) 및 ESDU 86022 (1986) 상관식의 적용 범위 내에서
검증용 스태거드 핀-튜브 열교환기 메시를 생성합니다.

기준 형상:
- d_o = 16 mm (튜브 외경)
- S_T = 36 mm (S_T/d_o = 2.25)
- S_L = 34 mm (S_L/d_o = 2.125)
- h_fin = 10 mm (핀 방사 높이)
- t_fin = 0.5 mm (핀 두께)
- S_fin = 4 mm (clear fin spacing)
- 핀 피치 = 4.5 mm
"""
import argparse
import os
import sys

import gmsh


# 검증용 형상 파라미터 (상관식 적용 범위 중앙값)
VALIDATION_PARAMS = {
    "SL": 0.034,           # 34 mm, longitudinal pitch
    "ST": 0.036,           # 36 mm, transverse pitch
    "R_tube_out": 0.008,   # 8 mm (d_o = 16 mm)
    "R_tube_in": 0.006,    # 6 mm
    "R_fin": 0.018,        # 18 mm (d_f = 36 mm, h_fin = 10 mm)
    "t_fin": 0.0005,       # 0.5 mm
    "P_fin": 0.0045,       # 4.5 mm (핀 피치 = t_fin + S_fin)
    "N_fins": 5,           # 핀 개수 (z방향)
    "mesh_size": 0.002,    # 기본 메시 크기
    "interface_size": 0.001,  # 인터페이스 메시 크기
    "N_repeat": 5,         # REV 반복 횟수 (x방향)
    "inlet_extension": 1.0,   # inlet extension in units of S_L
    "outlet_extension": 3.0,  # outlet extension in units of S_L
}


def _surface_sets_for_volume(volume_tag):
    """볼륨의 경계 서피스 집합을 반환"""
    boundary = gmsh.model.getBoundary([(3, volume_tag)], oriented=False, recursive=False)
    return {tag for dim, tag in boundary if dim == 2}


def _select_surfaces_by_coord(surfaces, axis, target, tol=1e-5):
    """특정 좌표에 있는 서피스 선택"""
    selected = []
    for tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        coord_min, coord_max = {
            "x": (xmin, xmax),
            "y": (ymin, ymax),
            "z": (zmin, zmax),
        }[axis]
        if abs(coord_min - target) < tol and abs(coord_max - target) < tol:
            selected.append(tag)
    return selected


def _add_geometry(params):
    """
    스태거드 핀-튜브 형상 생성 (다중 핀 배열).
    
    도메인 레이아웃:
    - x: 유동 방향 (inlet → outlet)
    - y: 가로 방향 (S_T)
    - z: 핀 적층 방향 (P_fin)
    """
    SL = params["SL"]
    ST = params["ST"]
    R_tube_out = params["R_tube_out"]
    R_tube_in = params["R_tube_in"]
    R_fin = params["R_fin"]
    t_fin = params["t_fin"]
    P_fin = params["P_fin"]
    N_fins = params["N_fins"]
    N = params["N_repeat"]
    
    # 확장 영역
    inlet_ext = params["inlet_extension"] * SL
    outlet_ext = params["outlet_extension"] * SL
    
    occ = gmsh.model.occ
    
    # 도메인 크기
    REV_length = 2.0 * SL
    core_length = N * REV_length
    domain_x = inlet_ext + core_length + outlet_ext
    domain_y = ST
    domain_z = N_fins * P_fin
    
    print(f"도메인 레이아웃:")
    print(f"  Inlet extension:  {inlet_ext*1000:.1f} mm")
    print(f"  HX core (N={N}):  {core_length*1000:.1f} mm")
    print(f"  Outlet extension: {outlet_ext*1000:.1f} mm")
    print(f"  Total length:     {domain_x*1000:.1f} mm")
    print(f"  Domain size:      {domain_x*1000:.1f} x {domain_y*1000:.1f} x {domain_z*1000:.1f} mm")
    print(f"  핀 개수:          {N_fins}")
    
    # 메인 유체 박스 생성
    fluid_box = occ.addBox(0.0, 0.0, 0.0, domain_x, domain_y, domain_z)
    
    # 튜브-핀 유닛 생성
    trimmed_solids = []
    
    for i in range(N):
        x_offset = inlet_ext + i * REV_length
        
        # 스태거드 배열의 튜브 중심 위치
        # Row 1: x = x_offset + S_L/2
        # Row 2: x = x_offset + 3*S_L/2
        tube_centers = [
            (x_offset + SL / 2.0, 0.0),           # Row 1: y=0 경계 반쪽 튜브
            (x_offset + SL / 2.0, ST),            # Row 1: y=S_T 경계 반쪽 튜브
            (x_offset + 3.0 * SL / 2.0, ST / 2.0),  # Row 2: 중앙 전체 튜브
        ]
        
        for cx, cy in tube_centers:
            # 튜브 외부 실린더 (전체 z 높이)
            tube_out = occ.addCylinder(cx, cy, 0.0, 0.0, 0.0, domain_z, R_tube_out)
            
            # 다중 핀 디스크 생성
            fin_tags = []
            for j in range(N_fins):
                z_center = (j + 0.5) * P_fin
                z_start = z_center - t_fin / 2.0
                fin = occ.addCylinder(cx, cy, z_start, 0.0, 0.0, t_fin, R_fin)
                fin_tags.append((3, fin))
            
            # 내부 보어
            tube_in = occ.addCylinder(cx, cy, 0.0, 0.0, 0.0, domain_z, R_tube_in)
            
            # 튜브와 핀들을 융합
            if fin_tags:
                solid_unit, _ = occ.fuse([(3, tube_out)], fin_tags)
            else:
                solid_unit = [(3, tube_out)]
            solid_tag = solid_unit[0][1]
            
            # 내부 보어 제거
            solid_cut, _ = occ.cut([(3, solid_tag)], [(3, tube_in)])
            solid_tag = solid_cut[0][1]
            
            # 도메인 박스로 트리밍
            domain_copy = occ.addBox(0.0, 0.0, 0.0, domain_x, domain_y, domain_z)
            trimmed, _ = occ.intersect([(3, solid_tag)], [(3, domain_copy)],
                                       removeObject=True, removeTool=True)
            if trimmed:
                trimmed_solids.extend(trimmed)
    
    if not trimmed_solids:
        raise RuntimeError("솔리드 볼륨이 생성되지 않았습니다.")
    
    # 유체 박스와 솔리드들의 conformal mesh를 위한 fragment
    frag_entities, _ = occ.fragment([(3, fluid_box)], trimmed_solids)
    occ.synchronize()
    
    # 볼륨 분류: 가장 큰 것이 유체
    volumes = gmsh.model.getEntities(dim=3)
    if len(volumes) < 2:
        raise RuntimeError("유체와 솔리드 볼륨이 모두 필요합니다.")
    
    volumes_with_mass = []
    for dim, tag in volumes:
        mass = gmsh.model.occ.getMass(dim, tag)
        volumes_with_mass.append((tag, mass))
    volumes_with_mass.sort(key=lambda item: item[1], reverse=True)
    
    fluid_volume = volumes_with_mass[0][0]
    fluid_surfaces = _surface_sets_for_volume(fluid_volume)
    
    # 유체와 인터페이스를 공유하는 솔리드 볼륨 찾기
    solid_volumes = []
    for dim, tag in volumes:
        if tag == fluid_volume:
            continue
        surfaces = _surface_sets_for_volume(tag)
        if fluid_surfaces.intersection(surfaces):
            solid_volumes.append(tag)
    
    if not solid_volumes:
        raise RuntimeError("유체-솔리드 인터페이스를 찾을 수 없습니다.")
    
    # 불필요한 파편 제거
    for dim, tag in volumes:
        if tag != fluid_volume and tag not in solid_volumes:
            gmsh.model.occ.remove([(dim, tag)], recursive=True)
    gmsh.model.occ.synchronize()
    
    return fluid_volume, solid_volumes, domain_x, domain_y, domain_z


def build_model(params):
    """Gmsh 모델 구축"""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("validation_finned_tube")
    
    fluid_volume, solid_volumes, domain_x, domain_y, domain_z = _add_geometry(params)
    
    fluid_surfaces = _surface_sets_for_volume(fluid_volume)
    
    # 모든 솔리드 서피스 수집
    all_solid_surfaces = set()
    for sv in solid_volumes:
        all_solid_surfaces.update(_surface_sets_for_volume(sv))
    interface_surfaces = fluid_surfaces.intersection(all_solid_surfaces)
    
    # 경계면 식별
    inlet = _select_surfaces_by_coord(fluid_surfaces, "x", 0.0)
    outlet = _select_surfaces_by_coord(fluid_surfaces, "x", domain_x)
    bottom = _select_surfaces_by_coord(fluid_surfaces, "y", 0.0)
    top = _select_surfaces_by_coord(fluid_surfaces, "y", domain_y)
    back = _select_surfaces_by_coord(fluid_surfaces, "z", 0.0)
    front = _select_surfaces_by_coord(fluid_surfaces, "z", domain_z)
    
    if not all([inlet, outlet, bottom, top, back, front]):
        raise RuntimeError("경계면 식별에 실패했습니다.")
    
    # Physical groups
    gmsh.model.addPhysicalGroup(3, [fluid_volume], name="fluid")
    gmsh.model.addPhysicalGroup(3, solid_volumes, name="solid")
    
    gmsh.model.addPhysicalGroup(2, inlet, name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet, name="outlet")
    gmsh.model.addPhysicalGroup(2, top, name="top")
    gmsh.model.addPhysicalGroup(2, bottom, name="bottom")
    gmsh.model.addPhysicalGroup(2, front, name="front")
    gmsh.model.addPhysicalGroup(2, back, name="back")
    gmsh.model.addPhysicalGroup(2, list(interface_surfaces), name="interface_fluid_solid")
    
    gmsh.model.occ.synchronize()
    
    # 메시 설정
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", params["mesh_size"] * 0.3)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", params["mesh_size"])
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT
    
    # 인터페이스 근처 메시 세분화
    if interface_surfaces:
        field_id = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_id, "FacesList", list(interface_surfaces))
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", field_id)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", params["interface_size"])
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", params["mesh_size"])
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.001)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.005)
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
    
    print(f"메시 생성 중...")
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    
    # 메시 통계
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    print(f"총 노드 수: {len(node_tags)}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate validation mesh for correlation verification."
    )
    parser.add_argument("--case", type=str, required=True,
                        help="Case directory name (e.g., case_Re2000)")
    parser.add_argument("--mesh-size", type=float, default=VALIDATION_PARAMS["mesh_size"],
                        help="Target mesh size (m)")
    parser.add_argument("-N", "--n-repeat", type=int, default=VALIDATION_PARAMS["N_repeat"],
                        help="Number of REV repetitions")
    parser.add_argument("--n-fins", type=int, default=VALIDATION_PARAMS["N_fins"],
                        help="Number of fins per tube")
    return parser.parse_args()


def main():
    args = parse_args()
    
    # 스크립트 디렉토리 기준 케이스 경로
    script_dir = os.path.dirname(os.path.abspath(__file__))
    study_dir = os.path.dirname(script_dir)
    case_dir = os.path.join(study_dir, args.case)
    
    if not os.path.exists(case_dir):
        os.makedirs(case_dir, exist_ok=True)
    
    params = dict(VALIDATION_PARAMS)
    params["mesh_size"] = args.mesh_size
    params["N_repeat"] = args.n_repeat
    params["N_fins"] = args.n_fins
    
    output_file = os.path.join(case_dir, "mesh.msh")
    
    try:
        build_model(params)
        gmsh.write(output_file)
    except Exception as exc:
        gmsh.finalize()
        raise SystemExit(f"메시 생성 실패: {exc}")
    
    gmsh.finalize()
    print(f"메시 저장: {output_file}")


if __name__ == "__main__":
    main()

