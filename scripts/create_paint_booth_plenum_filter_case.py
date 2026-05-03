#!/usr/bin/env python3
"""Create a plenum + porous-filter paint-booth case.

This builds on the validated smooth car shell but extends the domain vertically:
work zone z=0..2.95, thin porous filter layer z=2.95..3.05, and upper plenum
z=3.05..3.8. A central top supply opening feeds the plenum; the floor is a full
pressure outlet. The filter is represented as a thin cellZone populated with a
Darcy-Forchheimer fvOptions source.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

from create_paint_booth_baseline import foam_header, generate_car_shell_stl, write_text


def block_mesh_plenum_dict(cell_size: float = 0.25, central_filter_panel: bool = False, filter_z_cells: int = 2) -> str:
    # Split top face to create a central 1.5 m x 1.2 m supply opening.
    # For the central-filter variant, also split at the intended 7.0 x 3.2 m
    # filter-panel edges so the panel and sealed frame zones are cell-aligned.
    if central_filter_panel:
        xs = [-1.5, -1.0, 1.75, 3.25, 6.0, 6.5]
        ys = [-2.0, -1.6, -0.6, 0.6, 1.6, 2.0]
    else:
        xs = [-1.5, 1.75, 3.25, 6.5]
        ys = [-2.0, -0.6, 0.6, 2.0]
    # A 0.10 m thin filter layer centered around z=3.0.
    zs = [0.0, 2.95, 3.05, 3.8]

    def vid(i: int, j: int, k: int) -> int:
        return k * len(xs) * len(ys) + j * len(xs) + i

    vertices = []
    for z in zs:
        for y in ys:
            for x in xs:
                vertices.append((x, y, z))

    blocks = []
    for k in range(len(zs) - 1):
        dz = zs[k + 1] - zs[k]
        nz = max(1, round(dz / cell_size))
        if k == 1:
            nz = max(2, filter_z_cells)  # ensure the porous layer has enough cells in z
        for j in range(len(ys) - 1):
            dy = ys[j + 1] - ys[j]
            ny = max(2, round(dy / cell_size))
            for i in range(len(xs) - 1):
                dx = xs[i + 1] - xs[i]
                nx = max(2, round(dx / cell_size))
                v = [
                    vid(i, j, k),
                    vid(i + 1, j, k),
                    vid(i + 1, j + 1, k),
                    vid(i, j + 1, k),
                    vid(i, j, k + 1),
                    vid(i + 1, j, k + 1),
                    vid(i + 1, j + 1, k + 1),
                    vid(i, j + 1, k + 1),
                ]
                blocks.append((v, nx, ny, nz))

    # Boundary faces. For OpenFOAM outward orientation, reuse validated ordering
    # from the simple blockMeshDict where possible.
    floor_faces = []
    top_supply_faces = []
    top_wall_faces = []
    front_faces = []
    rear_faces = []
    side_faces = []

    k_top = len(zs) - 2
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            # bottom z=0: outward -z
            floor_faces.append((vid(i, j, 0), vid(i, j + 1, 0), vid(i + 1, j + 1, 0), vid(i + 1, j, 0)))
            # top z=max: outward +z
            ftop = (vid(i, j, len(zs) - 1), vid(i + 1, j, len(zs) - 1), vid(i + 1, j + 1, len(zs) - 1), vid(i, j + 1, len(zs) - 1))
            x0, x1 = xs[i], xs[i + 1]
            y0, y1 = ys[j], ys[j + 1]
            if abs(x0 - 1.75) < 1e-9 and abs(x1 - 3.25) < 1e-9 and abs(y0 + 0.6) < 1e-9 and abs(y1 - 0.6) < 1e-9:
                top_supply_faces.append(ftop)
            else:
                top_wall_faces.append(ftop)

    for k in range(len(zs) - 1):
        for j in range(len(ys) - 1):
            # x min, outward -x
            front_faces.append((vid(0, j, k + 1), vid(0, j + 1, k + 1), vid(0, j + 1, k), vid(0, j, k)))
            # x max, outward +x
            rear_faces.append((vid(len(xs) - 1, j, k), vid(len(xs) - 1, j + 1, k), vid(len(xs) - 1, j + 1, k + 1), vid(len(xs) - 1, j, k + 1)))
        for i in range(len(xs) - 1):
            # y min, outward -y
            side_faces.append((vid(i, 0, k), vid(i + 1, 0, k), vid(i + 1, 0, k + 1), vid(i, 0, k + 1)))
            # y max, outward +y
            side_faces.append((vid(i, len(ys) - 1, k + 1), vid(i + 1, len(ys) - 1, k + 1), vid(i + 1, len(ys) - 1, k), vid(i, len(ys) - 1, k)))

    def faces_text(faces):
        return "\n".join("        (" + " ".join(map(str, f)) + ")" for f in faces)

    vertices_text = "\n".join(f"    ({x} {y} {z})" for x, y, z in vertices)
    blocks_text = "\n".join(
        "    hex (" + " ".join(map(str, v)) + f") ({nx} {ny} {nz}) simpleGrading (1 1 1)"
        for v, nx, ny, nz in blocks
    )

    return f"""
{foam_header('blockMeshDict')}

scale 1;

vertices
(
{vertices_text}
);

blocks
(
{blocks_text}
);

edges ();

boundary
(
    supplyInlet
    {{
        type patch;
        faces
        (
{faces_text(top_supply_faces)}
        );
    }}
    plenumTopWall
    {{
        type wall;
        faces
        (
{faces_text(top_wall_faces)}
        );
    }}
    floorExhaust
    {{
        type patch;
        faces
        (
{faces_text(floor_faces)}
        );
    }}
    frontWall
    {{
        type wall;
        faces
        (
{faces_text(front_faces)}
        );
    }}
    rearWall
    {{
        type wall;
        faces
        (
{faces_text(rear_faces)}
        );
    }}
    sideWalls
    {{
        type wall;
        faces
        (
{faces_text(side_faces)}
        );
    }}
);

mergePatchPairs ();
"""


def write_common_system(
    case_dir: Path,
    car_stl_name: str,
    supply_velocity: float,
    filter_forchheimer: float,
    central_filter_panel: bool = False,
    frame_forchheimer: float = 5.0e7,
    cell_size: float = 0.25,
    filter_z_cells: int = 2,
    car_refinement_min: int = 2,
    car_refinement_max: int = 3,
    add_layers: bool = False,
    n_surface_layers: int = 0,
    expansion_ratio: float = 1.2,
    final_layer_thickness: float = 0.3,
    min_layer_thickness: float = 0.05,
    relative_layer_sizes: bool = True,
) -> None:
    write_text(
        case_dir / "system" / "blockMeshDict",
        block_mesh_plenum_dict(cell_size=cell_size, central_filter_panel=central_filter_panel, filter_z_cells=filter_z_cells),
    )
    write_text(
        case_dir / "system" / "snappyHexMeshDict",
        f"""
{foam_header('snappyHexMeshDict')}

castellatedMesh true;
snap            true;
addLayers       {str(add_layers).lower()};

geometry
{{
    {car_stl_name}
    {{
        type triSurfaceMesh;
        name carBody;
    }}
}}

castellatedMeshControls
{{
    maxLocalCells 200000;
    maxGlobalCells 3000000;
    minRefinementCells 0;
    nCellsBetweenLevels 3;
    features ();
    refinementSurfaces
    {{
        carBody
        {{
            level (2 3);
            patchInfo {{ type wall; }}
        }}
    }}
    resolveFeatureAngle 30;
    refinementRegions {{}}
    locationInMesh (2.0 0 2.0);
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
}}

addLayersControls
{{
    relativeSizes {str(relative_layer_sizes).lower()};
    layers
    {{
        carBody
        {{
            nSurfaceLayers {n_surface_layers};
        }}
    }}
    expansionRatio {expansion_ratio:g};
    finalLayerThickness {final_layer_thickness:g};
    minThickness {min_layer_thickness:g};
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
meshQualityControls {{ #include "meshQualityDict" }}
writeFlags (scalarLevels layerSets layerFields);
mergeTolerance 1e-6;
""",
    )
    write_text(
        case_dir / "system" / "meshQualityDict",
        f"""
{foam_header('meshQualityDict')}
maxNonOrtho         65;
maxBoundarySkewness 20;
maxInternalSkewness 4;
maxConcave          80;
minVol              1e-13;
minTetQuality       -1e30;
minArea             -1;
minTwist            0.01;
minDeterminant      0.001;
minFaceWeight       0.05;
minVolRatio         0.01;
minTriangleTwist    -1;
nSmoothScale        4;
errorReduction      0.75;
""",
    )
    write_text(
        case_dir / "system" / "controlDict",
        f"""
{foam_header('controlDict')}
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         200;
deltaT          1;
writeControl    timeStep;
writeInterval   50;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""",
    )
    write_text(
        case_dir / "system" / "fvSchemes",
        f"""
{foam_header('fvSchemes')}
ddtSchemes {{ default steadyState; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default none;
    div(phi,U) bounded Gauss upwind;
    div(phi,k) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
wallDist {{ method meshWave; }}
""",
    )
    write_text(
        case_dir / "system" / "fvSolution",
        f"""
{foam_header('fvSolution')}
solvers
{{
    p {{ solver GAMG; tolerance 1e-7; relTol 0.1; smoother GaussSeidel; }}
    U {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }}
    "(k|omega)" {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }}
}}
SIMPLE
{{
    nNonOrthogonalCorrectors 0;
    consistent yes;
}}
relaxationFactors
{{
    fields {{ p 0.3; }}
    equations {{ U 0.6; k 0.6; omega 0.6; }}
}}
""",
    )
    if central_filter_panel:
        topo_actions = """
    {
        name    filterZone;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (-1.0 -1.6 2.95) (6.0 1.6 3.05);
    }
    {
        name    sealedFrameZone;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (-1.5 -2.0 2.95) (-1.0 2.0 3.05);
    }
    {
        name    sealedFrameZone;
        type    cellZoneSet;
        action  add;
        source  boxToCell;
        box     (6.0 -2.0 2.95) (6.5 2.0 3.05);
    }
    {
        name    sealedFrameZone;
        type    cellZoneSet;
        action  add;
        source  boxToCell;
        box     (-1.0 -2.0 2.95) (6.0 -1.6 3.05);
    }
    {
        name    sealedFrameZone;
        type    cellZoneSet;
        action  add;
        source  boxToCell;
        box     (-1.0 1.6 2.95) (6.0 2.0 3.05);
    }
"""
    else:
        topo_actions = """
    {
        name    filterZone;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (-1.5 -2.0 2.95) (6.5 2.0 3.05);
    }
"""
    write_text(
        case_dir / "system" / "topoSetDict",
        f"""
{foam_header('topoSetDict')}
actions
(
{topo_actions}
);
""",
    )
    frame_source = """
sealedFramePorosity
{
    type            explicitPorositySource;
    active          yes;

    explicitPorositySourceCoeffs
    {
        selectionMode   cellZone;
        cellZone        sealedFrameZone;

        type            DarcyForchheimer;
        DarcyForchheimerCoeffs
        {
            d   (0 0 0);
            f   (FRAME_F FRAME_F FRAME_F);
            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }
            }
        }
    }
}
""".replace("FRAME_F", f"{frame_forchheimer:g}") if central_filter_panel else ""
    write_text(
        case_dir / "system" / "fvOptions",
        f"""
{foam_header('fvOptions')}
filterPorosity
{{
    type            explicitPorositySource;
    active          yes;

    explicitPorositySourceCoeffs
    {{
        selectionMode   cellZone;
        cellZone        filterZone;

        type            DarcyForchheimer;
        DarcyForchheimerCoeffs
        {{
            d   (0 0 0);
            f   ({filter_forchheimer:g} {filter_forchheimer:g} {filter_forchheimer:g});
            coordinateSystem
            {{
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {{
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }}
            }}
        }}
    }}
}}
{frame_source}
""",
    )


def write_fields(case_dir: Path, supply_velocity: float = 4.36) -> None:
    # Supply duct turbulence estimate, I=5%, l=0.15 m.
    k = 1.5 * (supply_velocity * 0.05) ** 2
    omega = math.sqrt(k) / max((0.09 ** 0.25) * 0.15, 1e-9)
    wall_u = """
    plenumTopWall { type noSlip; }
    frontWall { type noSlip; }
    rearWall { type noSlip; }
    sideWalls { type noSlip; }
    carBody { type noSlip; }
"""
    write_text(
        case_dir / "0" / "U",
        f"""
{foam_header('U', 'volVectorField')}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 -0.35);
boundaryField
{{
    supplyInlet {{ type fixedValue; value uniform (0 0 -{supply_velocity}); }}
    floorExhaust {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
{wall_u}
}}
""",
    )
    write_text(
        case_dir / "0" / "p",
        f"""
{foam_header('p', 'volScalarField')}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    supplyInlet {{ type zeroGradient; }}
    floorExhaust {{ type fixedValue; value uniform 0; }}
    plenumTopWall {{ type zeroGradient; }}
    frontWall {{ type zeroGradient; }}
    rearWall {{ type zeroGradient; }}
    sideWalls {{ type zeroGradient; }}
    carBody {{ type zeroGradient; }}
}}
""",
    )
    write_text(
        case_dir / "0" / "k",
        f"""
{foam_header('k', 'volScalarField')}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k:.6g};
boundaryField
{{
    supplyInlet {{ type fixedValue; value uniform {k:.6g}; }}
    floorExhaust {{ type inletOutlet; inletValue uniform {k:.6g}; value uniform {k:.6g}; }}
    plenumTopWall {{ type kqRWallFunction; value uniform {k:.6g}; }}
    frontWall {{ type kqRWallFunction; value uniform {k:.6g}; }}
    rearWall {{ type kqRWallFunction; value uniform {k:.6g}; }}
    sideWalls {{ type kqRWallFunction; value uniform {k:.6g}; }}
    carBody {{ type kqRWallFunction; value uniform {k:.6g}; }}
}}
""",
    )
    write_text(
        case_dir / "0" / "omega",
        f"""
{foam_header('omega', 'volScalarField')}
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega:.6g};
boundaryField
{{
    supplyInlet {{ type fixedValue; value uniform {omega:.6g}; }}
    floorExhaust {{ type inletOutlet; inletValue uniform {omega:.6g}; value uniform {omega:.6g}; }}
    plenumTopWall {{ type omegaWallFunction; value uniform {omega:.6g}; }}
    frontWall {{ type omegaWallFunction; value uniform {omega:.6g}; }}
    rearWall {{ type omegaWallFunction; value uniform {omega:.6g}; }}
    sideWalls {{ type omegaWallFunction; value uniform {omega:.6g}; }}
    carBody {{ type omegaWallFunction; value uniform {omega:.6g}; }}
}}
""",
    )
    write_text(
        case_dir / "0" / "nut",
        f"""
{foam_header('nut', 'volScalarField')}
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    supplyInlet {{ type calculated; value uniform 0; }}
    floorExhaust {{ type calculated; value uniform 0; }}
    plenumTopWall {{ type nutkWallFunction; value uniform 0; }}
    frontWall {{ type nutkWallFunction; value uniform 0; }}
    rearWall {{ type nutkWallFunction; value uniform 0; }}
    sideWalls {{ type nutkWallFunction; value uniform 0; }}
    carBody {{ type nutkWallFunction; value uniform 0; }}
}}
""",
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--case-dir", type=Path, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--supply-velocity", type=float, default=4.36)
    p.add_argument("--filter-forchheimer", type=float, default=6800.0)
    p.add_argument("--central-filter-panel-frame", action="store_true", help="Use a 7.0 x 3.2 m central porous panel and a high-resistance sealed frame in the filter layer.")
    p.add_argument("--frame-forchheimer", type=float, default=5.0e7, help="High Forchheimer coefficient used to numerically seal the outer filter-frame zone.")
    p.add_argument("--cell-size", type=float, default=0.25, help="Target base blockMesh cell size [m].")
    p.add_argument("--filter-z-cells", type=int, default=2, help="Number of base cells through the 0.10 m filter layer.")
    p.add_argument("--car-refinement-min", type=int, default=2, help="snappyHexMesh min refinement level on carBody.")
    p.add_argument("--car-refinement-max", type=int, default=3, help="snappyHexMesh max refinement level on carBody.")
    p.add_argument("--add-layers", action="store_true", help="Enable snappyHexMesh prism layers on carBody.")
    p.add_argument("--n-surface-layers", type=int, default=4, help="Number of prism layers on carBody when --add-layers is used.")
    p.add_argument("--expansion-ratio", type=float, default=1.2, help="Boundary-layer expansion ratio.")
    p.add_argument("--final-layer-thickness", type=float, default=0.3, help="snappy finalLayerThickness. Relative if --absolute-layer-sizes is not set.")
    p.add_argument("--min-layer-thickness", type=float, default=0.05, help="snappy minThickness. Relative if --absolute-layer-sizes is not set.")
    p.add_argument("--absolute-layer-sizes", action="store_true", help="Use absolute layer thickness values instead of relative sizes.")
    args = p.parse_args()

    case_dir = args.case_dir or (
        Path("cases/paint_booth_plenum_filter_panel_frame_v036")
        if args.central_filter_panel_frame
        else Path("cases/paint_booth_plenum_filter_full_floor_v035")
    )
    if case_dir.exists() and args.force:
        shutil.rmtree(case_dir)
    if case_dir.exists() and any(case_dir.iterdir()) and not args.force:
        raise SystemExit(f"Case directory already exists and is non-empty: {case_dir}. Use --force to overwrite.")

    for d in ["0", "constant/triSurface", "system"]:
        (case_dir / d).mkdir(parents=True, exist_ok=True)

    stl_rel = Path("constant/triSurface/simplified_car_shell.stl")
    stats = generate_car_shell_stl(case_dir / stl_rel, window_openings=False)
    write_common_system(
        case_dir,
        stl_rel.name,
        args.supply_velocity,
        args.filter_forchheimer,
        central_filter_panel=args.central_filter_panel_frame,
        frame_forchheimer=args.frame_forchheimer,
        cell_size=args.cell_size,
        filter_z_cells=args.filter_z_cells,
        car_refinement_min=args.car_refinement_min,
        car_refinement_max=args.car_refinement_max,
        add_layers=args.add_layers,
        n_surface_layers=args.n_surface_layers if args.add_layers else 0,
        expansion_ratio=args.expansion_ratio,
        final_layer_thickness=args.final_layer_thickness,
        min_layer_thickness=args.min_layer_thickness,
        relative_layer_sizes=not args.absolute_layer_sizes,
    )
    write_fields(case_dir, args.supply_velocity)
    write_text(
        case_dir / "constant" / "transportProperties",
        f"""
{foam_header('transportProperties')}
transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] 1.5e-05;
""",
    )
    write_text(
        case_dir / "constant" / "turbulenceProperties",
        f"""
{foam_header('turbulenceProperties')}
simulationType RAS;
RAS
{{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}}
""",
    )
    metadata = {
        "case": str(case_dir),
        "model": "plenum_filter_panel_frame_v036" if args.central_filter_panel_frame else "plenum_filter_full_floor_v035",
        "domain_m": {"x": [-1.5, 6.5], "y": [-2.0, 2.0], "z": [0.0, 3.8]},
        "work_zone_z_m": [0.0, 2.95],
        "filter_layer_z_m": [2.95, 3.05],
        "filter_panel_m": {"x": [-1.0, 6.0], "y": [-1.6, 1.6], "area": 22.4} if args.central_filter_panel_frame else "full cross-section",
        "sealed_frame_model": {"type": "high-resistance porous zone", "forchheimer_coeff": args.frame_forchheimer} if args.central_filter_panel_frame else None,
        "plenum_z_m": [3.05, 3.8],
        "supply_opening_m": {"x": [1.75, 3.25], "y": [-0.6, 0.6], "z": 3.8, "area": 1.8},
        "target_filter_velocity_mps": 0.35,
        "supply_velocity_mps": args.supply_velocity,
        "estimated_flow_rate_m3s": 1.8 * args.supply_velocity,
        "filter_forchheimer_coeff": args.filter_forchheimer,
        "mesh_settings": {
            "base_cell_size_m": args.cell_size,
            "filter_z_cells": args.filter_z_cells,
            "car_refinement_level": [args.car_refinement_min, args.car_refinement_max],
            "add_layers": args.add_layers,
            "n_surface_layers": args.n_surface_layers if args.add_layers else 0,
            "expansion_ratio": args.expansion_ratio,
            "final_layer_thickness": args.final_layer_thickness,
            "min_layer_thickness": args.min_layer_thickness,
            "relative_layer_sizes": not args.absolute_layer_sizes,
        },
        "car_geometry": stats,
    }
    write_text(case_dir / "constant" / "geometry_metadata.json", json.dumps(metadata, indent=2))
    print(json.dumps(metadata, indent=2))
    print(f"Created {case_dir}")


if __name__ == "__main__":
    main()
