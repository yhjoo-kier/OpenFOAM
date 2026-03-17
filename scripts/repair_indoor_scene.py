#!/usr/bin/env python3
"""Repair/simplify indoor CFD scenes to be more solver-friendly.

Goal: preserve semantic layout while enforcing geometry rules that improve
meshing and CFD stability.
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
import sys
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from validate_indoor_scene import validate_scene  # noqa: E402

EPS = 1e-9


def clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def box_max(obs: dict) -> dict:
    return {
        'x': obs['min']['x'] + obs['size']['dx'],
        'y': obs['min']['y'] + obs['size']['dy'],
        'z': obs['min']['z'] + obs['size']['dz'],
    }


def overlaps(a: dict, b: dict, margin: float = 0.0) -> bool:
    amax = box_max(a)
    bmax = box_max(b)
    return (
        a['min']['x'] < bmax['x'] - margin - EPS and amax['x'] > b['min']['x'] + margin + EPS
        and a['min']['y'] < bmax['y'] - margin - EPS and amax['y'] > b['min']['y'] + margin + EPS
        and a['min']['z'] < bmax['z'] - margin - EPS and amax['z'] > b['min']['z'] + margin + EPS
    )


def enforce_room_bounds(scene: dict, wall_margin: float) -> None:
    room = scene['room']['size']
    for obs in scene['obstacles']:
        obs['size']['dx'] = min(obs['size']['dx'], room['Lx'] - 2 * wall_margin)
        obs['size']['dy'] = min(obs['size']['dy'], room['Ly'] - 2 * wall_margin)
        obs['size']['dz'] = min(obs['size']['dz'], room['Lz'] - wall_margin)
        obs['min']['x'] = clamp(obs['min']['x'], wall_margin, room['Lx'] - wall_margin - obs['size']['dx'])
        obs['min']['y'] = clamp(obs['min']['y'], wall_margin, room['Ly'] - wall_margin - obs['size']['dy'])
        obs['min']['z'] = clamp(obs['min']['z'], 0.0, room['Lz'] - wall_margin - obs['size']['dz'])


def simplify_obstacles(scene: dict, wall_margin: float, max_count: int, max_row_length_ratio: float) -> None:
    room = scene['room']['size']
    obs = sorted(scene['obstacles'], key=lambda o: o['size']['dx'] * o['size']['dy'] * o['size']['dz'], reverse=True)
    obs = obs[:max_count]
    max_len_y = room['Ly'] * max_row_length_ratio
    max_len_x = room['Lx'] * max_row_length_ratio
    for o in obs:
        # Shorten very long obstacles that tend to create unstable channels.
        if o['size']['dy'] > max_len_y:
            o['size']['dy'] = max_len_y
        if o['size']['dx'] > max_len_x:
            o['size']['dx'] = max_len_x
        # Limit tall obstacles slightly below ceiling.
        o['size']['dz'] = min(o['size']['dz'], room['Lz'] - wall_margin)
        # Slightly thicken ultra-thin objects to avoid degenerate passages.
        o['size']['dx'] = max(o['size']['dx'], 0.6)
        o['size']['dy'] = max(o['size']['dy'], 0.8)
        o['size']['dz'] = max(o['size']['dz'], 1.5)
    scene['obstacles'] = obs


def push_apart_obstacles(scene: dict, clearance: float, max_iters: int = 20) -> None:
    room = scene['room']['size']
    obs = scene['obstacles']
    for _ in range(max_iters):
        moved = False
        for i in range(len(obs)):
            for j in range(i + 1, len(obs)):
                a, b = obs[i], obs[j]
                if overlaps(a, b, margin=clearance):
                    # Prefer moving the smaller one.
                    ta, tb = a['size']['dx'] * a['size']['dy'], b['size']['dx'] * b['size']['dy']
                    target = a if ta < tb else b
                    other = b if target is a else a
                    other_max = box_max(other)
                    # Push in +x first, else +y.
                    new_x = other_max['x'] + clearance
                    if new_x + target['size']['dx'] <= room['Lx'] - clearance:
                        target['min']['x'] = new_x
                    else:
                        new_y = other_max['y'] + clearance
                        target['min']['y'] = min(new_y, room['Ly'] - clearance - target['size']['dy'])
                    moved = True
        if not moved:
            break


def enforce_clear_zones(scene: dict, opening_clearance: float) -> None:
    room = scene['room']['size']
    openings = scene.get('openings', [])
    if len(openings) < 2:
        return
    for op in openings:
        if op['wall'] == 'west':
            x_limit = opening_clearance
            for obs in scene['obstacles']:
                if obs['min']['x'] < x_limit:
                    obs['min']['x'] = x_limit
        elif op['wall'] == 'east':
            x_limit = room['Lx'] - opening_clearance
            for obs in scene['obstacles']:
                if box_max(obs)['x'] > x_limit:
                    obs['min']['x'] = x_limit - obs['size']['dx']
        elif op['wall'] == 'south':
            y_limit = opening_clearance
            for obs in scene['obstacles']:
                if obs['min']['y'] < y_limit:
                    obs['min']['y'] = y_limit
        elif op['wall'] == 'north':
            y_limit = room['Ly'] - opening_clearance
            for obs in scene['obstacles']:
                if box_max(obs)['y'] > y_limit:
                    obs['min']['y'] = y_limit - obs['size']['dy']


def normalize_openings(scene: dict, wall: str | None = None, size_du: float | None = None, size_dv: float | None = None) -> None:
    room = scene['room']['size']
    openings = scene.get('openings', [])
    if len(openings) != 2:
        return
    inlet = next((o for o in openings if o['type'] == 'inlet'), openings[0])
    outlet = next((o for o in openings if o['type'] == 'outlet'), openings[-1])
    if wall == 'x':
        inlet['wall'] = 'west'
        outlet['wall'] = 'east'
        center_u = room['Ly'] * 0.5
        center_v = min(1.4, room['Lz'] * 0.5)
        inlet['center'] = {'u': center_u, 'v': center_v}
        outlet['center'] = {'u': center_u, 'v': center_v}
    elif wall == 'y':
        inlet['wall'] = 'south'
        outlet['wall'] = 'north'
        center_u = room['Lx'] * 0.5
        center_v = min(1.4, room['Lz'] * 0.5)
        inlet['center'] = {'u': center_u, 'v': center_v}
        outlet['center'] = {'u': center_u, 'v': center_v}
    if size_du is not None and size_dv is not None:
        inlet['size'] = {'du': size_du, 'dv': size_dv}
        outlet['size'] = {'du': size_du, 'dv': size_dv}


def repair_scene(scene: dict) -> tuple[dict, dict]:
    repaired = deepcopy(scene)
    room = repaired['room']['size']

    wall_margin = max(0.4, 0.05 * min(room['Lx'], room['Ly']))
    obstacle_clearance = max(0.6, 0.08 * min(room['Lx'], room['Ly']))
    opening_clearance = max(1.2, 0.12 * room['Lx'])

    simplify_obstacles(repaired, wall_margin=wall_margin, max_count=3, max_row_length_ratio=0.62)
    normalize_openings(repaired, wall='x', size_du=1.0, size_dv=0.8)
    enforce_room_bounds(repaired, wall_margin=wall_margin)
    enforce_clear_zones(repaired, opening_clearance=opening_clearance)
    push_apart_obstacles(repaired, clearance=obstacle_clearance)
    enforce_room_bounds(repaired, wall_margin=wall_margin)

    repaired.setdefault('meta', {})
    repaired['meta']['repair_applied'] = True
    repaired['meta']['repair_notes'] = (
        'Applied solver-friendly repair: obstacle count capped, long rows shortened, '
        'west-east openings normalized, wall margin and clearance enforced.'
    )

    report = validate_scene(repaired)
    info = {
        'wall_margin': wall_margin,
        'obstacle_clearance': obstacle_clearance,
        'opening_clearance': opening_clearance,
        'valid': report.ok,
        'errors': report.errors,
        'warnings': report.warnings,
    }
    return repaired, info


def main() -> int:
    parser = argparse.ArgumentParser(description='Repair indoor scene JSON for solver stability')
    parser.add_argument('scene_json', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    args = parser.parse_args()

    scene = json.loads(args.scene_json.read_text(encoding='utf-8'))
    repaired, info = repair_scene(scene)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(repaired, indent=2), encoding='utf-8')
    print(json.dumps({'output': str(args.output), 'info': info}, indent=2))
    return 0 if info['valid'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
