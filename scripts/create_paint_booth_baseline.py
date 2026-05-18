#!/usr/bin/env python3
"""Create a CFD-friendly baseline paint-booth case with a simplified car shell.

This is intentionally procedural/reproducible: instead of using an over-detailed
scan mesh, it generates a watertight low-detail sedan-like STL and a first-pass
OpenFOAM simpleFoam/snappyHexMesh case skeleton.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path
from typing import Iterable


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def foam_header(obj: str, cls: str = "dictionary") -> str:
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}"""


def superellipse_radius(theta: float, a: float, b: float, n: float = 3.2) -> tuple[float, float]:
    """Return y,z offsets for a superellipse cross-section."""
    c = math.cos(theta)
    s = math.sin(theta)
    y = a * math.copysign(abs(c) ** (2.0 / n), c)
    z = b * math.copysign(abs(s) ** (2.0 / n), s)
    return y, z


def smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def car_profiles(
    x: float,
    length: float,
    width: float,
    body_height: float,
    ground: float,
    *,
    cabin_bulge: bool = True,
) -> tuple[float, float, float]:
    """Return half width, bottom z, top z for a sedan-like body at x."""
    xi = x / length
    # Pointed/rounded ends with broad middle.
    end_taper = math.sin(math.pi * xi) ** 0.32
    half_w = 0.5 * width * max(0.035, end_taper)

    # Hood/trunk lower; optional cabin bulge for the original closed shell.
    base_top = ground + body_height * (0.46 + 0.16 * math.sin(math.pi * xi))
    cabin = smoothstep((xi - 0.30) / 0.18) * (1.0 - smoothstep((xi - 0.70) / 0.18))
    top_z = base_top + (cabin * body_height * 0.44 if cabin_bulge else 0.0)

    # Slight nose/trunk lowering.
    top_z -= body_height * 0.10 * (1.0 - end_taper)
    bottom_z = ground
    return half_w, bottom_z, top_z


def add_box(
    verts: list[tuple[float, float, float]],
    tris: list[tuple[int, int, int]],
    bounds: tuple[float, float, float, float, float, float],
) -> None:
    """Append a closed axis-aligned box with outward triangle winding."""
    x0, x1, y0, y1, z0, z1 = bounds
    base = len(verts)
    verts.extend(
        [
            (x0, y0, z0),
            (x1, y0, z0),
            (x1, y1, z0),
            (x0, y1, z0),
            (x0, y0, z1),
            (x1, y0, z1),
            (x1, y1, z1),
            (x0, y1, z1),
        ]
    )
    faces = [
        (0, 3, 2, 1),  # bottom, -z
        (4, 5, 6, 7),  # top, +z
        (0, 1, 5, 4),  # -y
        (3, 7, 6, 2),  # +y
        (0, 4, 7, 3),  # -x
        (1, 2, 6, 5),  # +x
    ]
    for a, b, c, d in faces:
        tris.append((base + a, base + b, base + c))
        tris.append((base + a, base + c, base + d))


def generate_car_shell_stl(
    out: Path,
    length: float = 4.5,
    width: float = 1.8,
    body_height: float = 1.35,
    ground_clearance: float = 0.18,
    nx: int = 80,
    nt: int = 56,
    x0: float = 0.0,
    window_openings: bool = False,
) -> dict:
    """Generate a CFD-friendly sedan surrogate as ASCII STL.

    Coordinates: x streamwise, y lateral, z vertical. The default is the original
    smooth single-shell car body. A multi-part window-opening prototype remains
    available for experiments but is disabled by default because the boxy pillars
    looked too artificial for the current baseline stage.
    """
    verts: list[tuple[float, float, float]] = []
    rings: list[list[int]] = []
    for i in range(nx):
        x = x0 + length * i / (nx - 1)
        half_w, bottom_z, top_z = car_profiles(
            x - x0,
            length,
            width,
            body_height,
            ground_clearance,
            cabin_bulge=not window_openings,
        )
        zc = 0.5 * (bottom_z + top_z)
        rz = 0.5 * (top_z - bottom_z)
        ring = []
        for j in range(nt):
            theta = 2.0 * math.pi * j / nt
            yoff, zoff = superellipse_radius(theta, half_w, rz)
            verts.append((x, yoff, zc + zoff))
            ring.append(len(verts) - 1)
        rings.append(ring)

    tris: list[tuple[int, int, int]] = []
    # Surface between rings.
    for i in range(nx - 1):
        for j in range(nt):
            a = rings[i][j]
            b = rings[i][(j + 1) % nt]
            c = rings[i + 1][j]
            d = rings[i + 1][(j + 1) % nt]
            # Use outward normals for the tube-like side surface.
            tris.append((a, b, c))
            tris.append((b, d, c))

    # End caps.
    for i, reverse in [(0, True), (nx - 1, False)]:
        ring = rings[i]
        cx = sum(verts[k][0] for k in ring) / nt
        cy = sum(verts[k][1] for k in ring) / nt
        cz = sum(verts[k][2] for k in ring) / nt
        verts.append((cx, cy, cz))
        center = len(verts) - 1
        for j in range(nt):
            a = ring[j]
            b = ring[(j + 1) % nt]
            tris.append((center, b, a) if reverse else (center, a, b))

    opening_parts: list[dict[str, object]] = []
    if window_openings:
        # Simple BIW-like cabin frame. Parts are intentionally low-detail closed
        # solids so that the side windows are real flow gaps, not just textures.
        y_outer = 0.82
        y_inner = 0.60
        belt_z = 0.82
        roof_z0 = 1.42
        roof_z1 = 1.62
        parts = [
            ("roof_rail", (1.18, 3.75, -0.63, 0.63, roof_z0, roof_z1)),
            ("left_A_pillar", (1.12, 1.34, y_inner, y_outer, belt_z, roof_z1)),
            ("left_B_pillar", (2.25, 2.43, y_inner, y_outer, belt_z, roof_z1)),
            ("left_C_pillar", (3.50, 3.76, y_inner, y_outer, belt_z, roof_z1)),
            ("right_A_pillar", (1.12, 1.34, -y_outer, -y_inner, belt_z, roof_z1)),
            ("right_B_pillar", (2.25, 2.43, -y_outer, -y_inner, belt_z, roof_z1)),
            ("right_C_pillar", (3.50, 3.76, -y_outer, -y_inner, belt_z, roof_z1)),
            ("left_side_sill", (1.08, 3.82, y_inner, y_outer, 0.68, 0.84)),
            ("right_side_sill", (1.08, 3.82, -y_outer, -y_inner, 0.68, 0.84)),
            ("front_header", (1.12, 1.42, -0.63, 0.63, 1.26, roof_z1)),
            ("rear_header", (3.44, 3.76, -0.63, 0.63, 1.22, roof_z1)),
        ]
        for name, bounds in parts:
            add_box(verts, tris, bounds)
            opening_parts.append({"name": name, "bounds": bounds})

    def normal(v1, v2, v3):
        ax, ay, az = (v2[k] - v1[k] for k in range(3))
        bx, by, bz = (v3[k] - v1[k] for k in range(3))
        nx_ = ay * bz - az * by
        ny_ = az * bx - ax * bz
        nz_ = ax * by - ay * bx
        mag = math.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_) or 1.0
        return nx_ / mag, ny_ / mag, nz_ / mag

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("solid simplified_car_shell\n")
        for ia, ib, ic in tris:
            v1, v2, v3 = verts[ia], verts[ib], verts[ic]
            nx_, ny_, nz_ = normal(v1, v2, v3)
            f.write(f"  facet normal {nx_:.9g} {ny_:.9g} {nz_:.9g}\n")
            f.write("    outer loop\n")
            for v in (v1, v2, v3):
                f.write(f"      vertex {v[0]:.9g} {v[1]:.9g} {v[2]:.9g}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write("endsolid simplified_car_shell\n")

    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    return {
        "stl": str(out),
        "vertices": len(verts),
        "triangles": len(tris),
        "bounds": {
            "x": [min(xs), max(xs)],
            "y": [min(ys), max(ys)],
            "z": [min(zs), max(zs)],
        },
        "dimensions_m": {"length": length, "width": width, "body_height": body_height, "ground_clearance": ground_clearance},
        "resolution": {"nx": nx, "ntheta": nt},
        "window_openings": window_openings,
        "opening_strategy": "multi_part_watertight_lower_body_roof_ABC_pillars" if window_openings else "single_closed_smooth_shell",
        "opening_parts": opening_parts,
    }


def block_mesh_dict(length=8.0, width=4.0, height=3.0, cell_size=0.25) -> str:
    nx = max(4, round(length / cell_size))
    ny = max(4, round(width / cell_size))
    nz = max(4, round(height / cell_size))
    x0, x1 = -1.5, length - 1.5
    y0, y1 = -width / 2, width / 2
    z0, z1 = 0.0, height
    return f"""
{foam_header('blockMeshDict')}

scale 1;

vertices
(
    ({x0} {y0} {z0}) ({x1} {y0} {z0}) ({x1} {y1} {z0}) ({x0} {y1} {z0})
    ({x0} {y0} {z1}) ({x1} {y0} {z1}) ({x1} {y1} {z1}) ({x0} {y1} {z1})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges ();

boundary
(
    inlet
    {{
        type patch;
        faces ((4 7 3 0));
    }}
    outlet
    {{
        type patch;
        faces ((1 2 6 5));
    }}
    ceiling
    {{
        type wall;
        faces ((4 5 6 7));
    }}
    floorExhaust
    {{
        type wall;
        faces ((0 3 2 1));
    }}
    sideWalls
    {{
        type wall;
        faces ((0 1 5 4) (3 7 6 2));
    }}
);

mergePatchPairs ();
"""


def write_case(case_dir: Path, car_stl_name: str, inlet_velocity: float = 0.35) -> None:
    for d in ["0", "constant/triSurface", "system"]:
        (case_dir / d).mkdir(parents=True, exist_ok=True)

    write_text(case_dir / "system" / "blockMeshDict", block_mesh_dict())
    write_text(
        case_dir / "system" / "snappyHexMeshDict",
        f"""
{foam_header('snappyHexMeshDict')}

castellatedMesh true;
snap            true;
addLayers       false;

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
    maxLocalCells 100000;
    maxGlobalCells 2000000;
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
    locationInMesh (2.0 0 2.5);
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
    relativeSizes true;
    layers {{}}
    expansionRatio 1.2;
    finalLayerThickness 0.3;
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
    #include "meshQualityDict"
}}

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
    equations {{ U 0.7; k 0.7; omega 0.7; }}
}}
""",
    )

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

    # Initial fields. Ceiling/floor are currently walls for cold-flow skeleton; later split floor exhaust/ceiling inlet.
    k = 1.5 * (inlet_velocity * 0.05) ** 2
    omega = math.sqrt(k) / max((0.09 ** 0.25) * 0.3, 1e-9)
    patches_wall = """
    ceiling { type noSlip; }
    floorExhaust { type noSlip; }
    sideWalls { type noSlip; }
    carBody { type noSlip; }
"""
    write_text(
        case_dir / "0" / "U",
        f"""
{foam_header('U', 'volVectorField')}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({inlet_velocity} 0 0);
boundaryField
{{
    inlet {{ type fixedValue; value uniform ({inlet_velocity} 0 0); }}
    outlet {{ type inletOutlet; inletValue uniform ({inlet_velocity} 0 0); value uniform ({inlet_velocity} 0 0); }}
{patches_wall}
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
    inlet {{ type zeroGradient; }}
    outlet {{ type fixedValue; value uniform 0; }}
    ceiling {{ type zeroGradient; }}
    floorExhaust {{ type zeroGradient; }}
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
internalField   uniform {k:.8g};
boundaryField
{{
    inlet {{ type fixedValue; value uniform {k:.8g}; }}
    outlet {{ type inletOutlet; inletValue uniform {k:.8g}; value uniform {k:.8g}; }}
    ceiling {{ type kqRWallFunction; value uniform {k:.8g}; }}
    floorExhaust {{ type kqRWallFunction; value uniform {k:.8g}; }}
    sideWalls {{ type kqRWallFunction; value uniform {k:.8g}; }}
    carBody {{ type kqRWallFunction; value uniform {k:.8g}; }}
}}
""",
    )
    write_text(
        case_dir / "0" / "omega",
        f"""
{foam_header('omega', 'volScalarField')}
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega:.8g};
boundaryField
{{
    inlet {{ type fixedValue; value uniform {omega:.8g}; }}
    outlet {{ type inletOutlet; inletValue uniform {omega:.8g}; value uniform {omega:.8g}; }}
    ceiling {{ type omegaWallFunction; value uniform {omega:.8g}; }}
    floorExhaust {{ type omegaWallFunction; value uniform {omega:.8g}; }}
    sideWalls {{ type omegaWallFunction; value uniform {omega:.8g}; }}
    carBody {{ type omegaWallFunction; value uniform {omega:.8g}; }}
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
    inlet {{ type calculated; value uniform 0; }}
    outlet {{ type calculated; value uniform 0; }}
    ceiling {{ type nutkWallFunction; value uniform 0; }}
    floorExhaust {{ type nutkWallFunction; value uniform 0; }}
    sideWalls {{ type nutkWallFunction; value uniform 0; }}
    carBody {{ type nutkWallFunction; value uniform 0; }}
}}
""",
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--case-dir", type=Path, default=Path("cases/paint_booth_baseline"))
    p.add_argument("--force", action="store_true")
    p.add_argument("--nx", type=int, default=80)
    p.add_argument("--ntheta", type=int, default=56)
    args = p.parse_args()

    case_dir = args.case_dir
    if case_dir.exists() and args.force:
        shutil.rmtree(case_dir)
    if case_dir.exists() and any(case_dir.iterdir()) and not args.force:
        raise SystemExit(f"Case directory already exists and is non-empty: {case_dir}. Use --force to overwrite.")

    stl_rel = Path("constant/triSurface/simplified_car_shell.stl")
    stats = generate_car_shell_stl(case_dir / stl_rel, nx=args.nx, nt=args.ntheta)
    write_case(case_dir, stl_rel.name)
    write_text(case_dir / "constant" / "geometry_metadata.json", json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))
    print(f"Created {case_dir}")


if __name__ == "__main__":
    main()
