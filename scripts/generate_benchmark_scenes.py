#!/usr/bin/env python3
"""Rule-based benchmark scene generator for the image-to-CFD paper.

Generates indoor_cfd_scene_v1 JSON scenes for the 2x2 benchmark matrix:
- A1: rectangular + simple obstacles
- A2: rectangular + complex obstacles
- A3: L-shaped/composite + simple obstacles
- A4: L-shaped/composite + complex obstacles

Design goals:
- reproducible with seeds
- validator-compatible output
- composite rooms represented with exactly two joined blocks
- conservative geometry to maximize meshing stability
- explicit manifest generation for frozen benchmark subsets
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
BENCHMARK_DIR = PROJECT_ROOT / "benchmark"
DEFAULT_OUTDIR = BENCHMARK_DIR / "scenes"
DEFAULT_MANIFEST_DIR = BENCHMARK_DIR / "manifests"
DEFAULT_SCENE_MANIFEST = DEFAULT_MANIFEST_DIR / "scene_manifest.json"

import sys
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from validate_indoor_scene import validate_scene  # noqa: E402

CLEARANCE_TO_WALL = 0.4
CLEARANCE_BETWEEN_OBS = 0.35
OPENING_MARGIN = 0.25
ROOM_HEIGHT_RANGE = (2.6, 3.2)


@dataclass(frozen=True)
class CategorySpec:
    code: str
    room_kind: str  # rectangular | composite
    obstacle_range: tuple[int, int]


CATEGORY_SPECS: dict[str, CategorySpec] = {
    "A1": CategorySpec("A1", "rectangular", (0, 1)),
    "A2": CategorySpec("A2", "rectangular", (2, 3)),
    "A3": CategorySpec("A3", "composite", (0, 1)),
    "A4": CategorySpec("A4", "composite", (3, 4)),
}

SIMPLE_CATEGORIES = {"A1", "A3"}


def randf(rng: random.Random, a: float, b: float, digits: int = 3) -> float:
    return round(rng.uniform(a, b), digits)


def rect_bounds(block: dict[str, Any]) -> tuple[float, float, float, float, float, float]:
    o = block["origin"]
    s = block["size"]
    return o["x"], o["y"], o["z"], o["x"] + s["dx"], o["y"] + s["dy"], o["z"] + s["dz"]


def box_overlap_xy(a: dict[str, Any], b: dict[str, Any], margin: float = 0.0) -> bool:
    ax0, ay0 = a["min"]["x"] - margin, a["min"]["y"] - margin
    ax1 = a["min"]["x"] + a["size"]["dx"] + margin
    ay1 = a["min"]["y"] + a["size"]["dy"] + margin
    bx0, by0 = b["min"]["x"] - margin, b["min"]["y"] - margin
    bx1 = b["min"]["x"] + b["size"]["dx"] + margin
    by1 = b["min"]["y"] + b["size"]["dy"] + margin
    return ax0 < bx1 and ax1 > bx0 and ay0 < by1 and ay1 > by0


def point_in_block_xy(x: float, y: float, block: dict[str, Any], margin: float = 0.0) -> bool:
    ox, oy, _, x1, y1, _ = rect_bounds(block)
    return ox + margin <= x <= x1 - margin and oy + margin <= y <= y1 - margin


def obstacle_inside_any_block(obstacle: dict[str, Any], blocks: list[dict[str, Any]], margin: float = 0.0) -> bool:
    x0 = obstacle["min"]["x"]
    y0 = obstacle["min"]["y"]
    x1 = x0 + obstacle["size"]["dx"]
    y1 = y0 + obstacle["size"]["dy"]
    for block in blocks:
        if point_in_block_xy(x0, y0, block, margin=margin) and point_in_block_xy(x1, y1, block, margin=margin):
            return True
    return False


def overall_room_size(blocks: list[dict[str, Any]]) -> dict[str, float]:
    xmax = max(b["origin"]["x"] + b["size"]["dx"] for b in blocks)
    ymax = max(b["origin"]["y"] + b["size"]["dy"] for b in blocks)
    zmax = max(b["origin"]["z"] + b["size"]["dz"] for b in blocks)
    return {"Lx": xmax, "Ly": ymax, "Lz": zmax}


def make_rectangular_room(rng: random.Random) -> list[dict[str, Any]]:
    lx = randf(rng, 4.8, 7.2)
    ly = randf(rng, 3.8, 5.8)
    lz = randf(rng, *ROOM_HEIGHT_RANGE)
    return [{"origin": {"x": 0.0, "y": 0.0, "z": 0.0}, "size": {"dx": lx, "dy": ly, "dz": lz}}]


def make_composite_room(rng: random.Random) -> list[dict[str, Any]]:
    lz = randf(rng, *ROOM_HEIGHT_RANGE)
    base_dx = randf(rng, 4.5, 6.5)
    base_dy = randf(rng, 3.6, 5.4)
    leg_dx = randf(rng, 1.8, 3.4)
    leg_dy = randf(rng, 1.6, max(1.8, base_dy - 1.0))
    return [
        {"origin": {"x": 0.0, "y": 0.0, "z": 0.0}, "size": {"dx": base_dx, "dy": base_dy, "dz": lz}},
        {"origin": {"x": base_dx, "y": 0.0, "z": 0.0}, "size": {"dx": leg_dx, "dy": leg_dy, "dz": lz}},
    ]


def pick_opening(blocks: list[dict[str, Any]], wall: str, opening_type: str, rng: random.Random) -> dict[str, Any]:
    room = overall_room_size(blocks)
    room_height = room["Lz"]

    if wall in {"west", "east"}:
        candidates = []
        plane_x = 0.0 if wall == "west" else room["Lx"]
        for b in blocks:
            bx0, by0, bz0, bx1, by1, bz1 = rect_bounds(b)
            if abs((bx0 if wall == "west" else bx1) - plane_x) < 1e-9:
                candidates.append((by0, by1, bz0, bz1))
        seg = rng.choice(candidates)
        du = randf(rng, 0.35, min(0.75, seg[1] - seg[0] - 2 * OPENING_MARGIN))
        dv = randf(rng, 0.35, min(0.8, room_height - 0.9))
        u = randf(rng, seg[0] + OPENING_MARGIN + du / 2, seg[1] - OPENING_MARGIN - du / 2)
        v = randf(rng, 0.7 + dv / 2, room_height - 0.2 - dv / 2)
    elif wall in {"south", "north"}:
        candidates = []
        plane_y = 0.0 if wall == "south" else room["Ly"]
        for b in blocks:
            bx0, by0, bz0, bx1, by1, bz1 = rect_bounds(b)
            if abs((by0 if wall == "south" else by1) - plane_y) < 1e-9:
                candidates.append((bx0, bx1, bz0, bz1))
        seg = rng.choice(candidates)
        du = randf(rng, 0.35, min(0.75, seg[1] - seg[0] - 2 * OPENING_MARGIN))
        dv = randf(rng, 0.35, min(0.8, room_height - 0.9))
        u = randf(rng, seg[0] + OPENING_MARGIN + du / 2, seg[1] - OPENING_MARGIN - du / 2)
        v = randf(rng, 0.7 + dv / 2, room_height - 0.2 - dv / 2)
    else:
        raise ValueError(f"Unsupported wall for benchmark opening placement: {wall}")

    return {
        "id": opening_type,
        "type": opening_type,
        "wall": wall,
        "center": {"u": u, "v": v},
        "size": {"du": du, "dv": dv},
    }


def sample_obstacle(idx: int, blocks: list[dict[str, Any]], existing: list[dict[str, Any]], rng: random.Random) -> dict[str, Any] | None:
    room_height = overall_room_size(blocks)["Lz"]
    for _ in range(200):
        block = rng.choice(blocks)
        bx0, by0, _, bx1, by1, _ = rect_bounds(block)
        dx = randf(rng, 0.6, min(1.4, (bx1 - bx0) - 2 * CLEARANCE_TO_WALL))
        dy = randf(rng, 0.8, min(1.8, (by1 - by0) - 2 * CLEARANCE_TO_WALL))
        dz = randf(rng, 1.0, min(2.2, room_height - 0.4))
        if (bx1 - bx0) <= dx + 2 * CLEARANCE_TO_WALL or (by1 - by0) <= dy + 2 * CLEARANCE_TO_WALL:
            continue
        x = randf(rng, bx0 + CLEARANCE_TO_WALL, bx1 - CLEARANCE_TO_WALL - dx)
        y = randf(rng, by0 + CLEARANCE_TO_WALL, by1 - CLEARANCE_TO_WALL - dy)
        obstacle = {
            "id": f"obs_{idx:02d}",
            "type": "box",
            "min": {"x": x, "y": y, "z": 0.0},
            "size": {"dx": dx, "dy": dy, "dz": dz},
        }
        if not obstacle_inside_any_block(obstacle, blocks):
            continue
        if any(box_overlap_xy(obstacle, other, margin=CLEARANCE_BETWEEN_OBS) for other in existing):
            continue
        return obstacle
    return None


def resolved_obstacle_range(category: str, *, min_simple_obstacles: int | None = None) -> tuple[int, int]:
    low, high = CATEGORY_SPECS[category].obstacle_range
    if category in SIMPLE_CATEGORIES and min_simple_obstacles is not None:
        low = max(low, min_simple_obstacles)
        high = max(high, low)
    return low, high


def build_scene(category: str, seed: int, variant_index: int, *, min_simple_obstacles: int | None = None) -> dict[str, Any]:
    spec = CATEGORY_SPECS[category]
    rng = random.Random(seed)
    blocks = make_rectangular_room(rng) if spec.room_kind == "rectangular" else make_composite_room(rng)
    room_size = overall_room_size(blocks)

    obstacles: list[dict[str, Any]] = []
    obstacle_range = resolved_obstacle_range(category, min_simple_obstacles=min_simple_obstacles)
    n_obs = rng.randint(*obstacle_range)
    for idx in range(1, n_obs + 1):
        obstacle = sample_obstacle(idx, blocks, obstacles, rng)
        if obstacle is not None:
            obstacles.append(obstacle)

    wall_pairs = [("west", "east"), ("south", "north"), ("west", "north"), ("south", "east")]
    inlet_wall, outlet_wall = rng.choice(wall_pairs)
    openings = [
        pick_opening(blocks, inlet_wall, "inlet", rng),
        pick_opening(blocks, outlet_wall, "outlet", rng),
    ]

    if spec.room_kind == "rectangular":
        room_obj: dict[str, Any] = {"size": room_size}
    else:
        room_obj = {"blocks": blocks}

    scene = {
        "schema_version": "indoor_cfd_scene_v1",
        "units": "m",
        "coordinate_system": {
            "origin": "room_min_corner",
            "x": "west_to_east",
            "y": "south_to_north",
            "z": "floor_to_ceiling",
        },
        "room": room_obj,
        "obstacles": obstacles,
        "openings": openings,
        "meta": {
            "scene_type": f"benchmark_{category.lower()}",
            "description": (
                f"Rule-based benchmark scene {category}-{variant_index:02d} "
                f"({spec.room_kind}, seed={seed}, obstacles={len(obstacles)})"
            ),
            "benchmark": {
                "category": category,
                "room_kind": spec.room_kind,
                "variant_index": variant_index,
                "seed": seed,
                "requested_min_simple_obstacles": min_simple_obstacles if category in SIMPLE_CATEGORIES else None,
            }
        },
    }
    return scene


def write_scene(scene: dict[str, Any], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(scene, indent=2), encoding="utf-8")


def scene_manifest_entry(scene_path: Path) -> dict[str, Any]:
    scene_path = scene_path.resolve()
    scene = json.loads(scene_path.read_text(encoding="utf-8"))
    benchmark = scene.get("meta", {}).get("benchmark", {})
    room = scene.get("room", {})
    if "blocks" in room:
        room_summary: dict[str, Any] = {
            "blocks": len(room["blocks"]),
            "block_sizes": [b["size"] for b in room["blocks"]],
        }
        room_kind = "composite"
    else:
        room_summary = room.get("size", {})
        room_kind = "rectangular"
    return {
        "scene_file": str(scene_path.relative_to(PROJECT_ROOT)),
        "case_name": f"bench_{scene_path.stem}",
        "category": benchmark.get("category"),
        "variant_index": benchmark.get("variant_index"),
        "seed": benchmark.get("seed"),
        "room_kind": room_kind,
        "obstacle_count": len(scene.get("obstacles", [])),
        "opening_walls": [op.get("wall") for op in scene.get("openings", [])],
        "room": room_summary,
    }


def write_scene_manifest(outdir: Path, manifest_path: Path = DEFAULT_SCENE_MANIFEST) -> list[dict[str, Any]]:
    entries = [scene_manifest_entry(path) for path in sorted(outdir.glob("*.json"))]
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(entries, indent=2), encoding="utf-8")
    return entries


def generate(
    category: str,
    count: int,
    base_seed: int,
    outdir: Path,
    overwrite: bool = False,
    *,
    start_index: int = 1,
    min_simple_obstacles: int | None = None,
) -> list[Path]:
    written: list[Path] = []
    for i in range(count):
        variant_index = start_index + i
        seed = base_seed + (variant_index - 1)
        scene = build_scene(
            category,
            seed=seed,
            variant_index=variant_index,
            min_simple_obstacles=min_simple_obstacles,
        )
        report = validate_scene(scene)
        if not report.ok:
            raise RuntimeError(
                f"Generated invalid scene for {category}-{variant_index:02d} seed={seed}:\n- " + "\n- ".join(report.errors)
            )
        path = outdir / f"{category.lower()}_{variant_index:02d}.json"
        if path.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing file: {path}")
        write_scene(scene, path)
        written.append(path)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate rule-based benchmark scenes for the image-to-CFD paper")
    parser.add_argument("--categories", nargs="+", default=["A1", "A2", "A3", "A4"], help="Subset of categories to generate")
    parser.add_argument("--count", type=int, default=3, help="Scenes per category")
    parser.add_argument("--base-seed", type=int, default=1000, help="Base seed; category offsets are added automatically")
    parser.add_argument("--start-index", type=int, default=1, help="Starting variant index within each category")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory for scene JSONs")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    parser.add_argument(
        "--min-simple-obstacles",
        type=int,
        default=None,
        choices=[0, 1],
        help="Optional lower bound for simple categories (A1/A3) when freezing a subset",
    )
    parser.add_argument(
        "--skip-manifest",
        action="store_true",
        help="Do not rewrite benchmark/manifests/scene_manifest.json after generation",
    )
    args = parser.parse_args()

    category_offsets = {"A1": 0, "A2": 100, "A3": 200, "A4": 300}
    all_written: list[Path] = []
    for category in args.categories:
        category = category.upper()
        if category not in CATEGORY_SPECS:
            raise SystemExit(f"Unsupported category: {category}")
        written = generate(
            category=category,
            count=args.count,
            base_seed=args.base_seed + category_offsets[category],
            outdir=args.outdir,
            overwrite=args.overwrite,
            start_index=args.start_index,
            min_simple_obstacles=args.min_simple_obstacles,
        )
        all_written.extend(written)

    manifest_entries = None if args.skip_manifest else write_scene_manifest(args.outdir)
    print(json.dumps({
        "ok": True,
        "outdir": str(args.outdir),
        "count": len(all_written),
        "files": [str(p) for p in all_written],
        "scene_manifest": str(DEFAULT_SCENE_MANIFEST) if manifest_entries is not None else None,
        "manifest_count": None if manifest_entries is None else len(manifest_entries),
        "min_simple_obstacles": args.min_simple_obstacles,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
