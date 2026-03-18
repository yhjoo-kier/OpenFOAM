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
OPENING_MARGIN = 0.25


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


def get_room_blocks(scene: dict) -> list[dict]:
    room = scene['room']
    if 'blocks' in room:
        return room['blocks']
    size = room['size']
    return [{
        'origin': {'x': 0.0, 'y': 0.0, 'z': 0.0},
        'size': {'dx': size['Lx'], 'dy': size['Ly'], 'dz': size['Lz']},
    }]


def overall_room_size(scene: dict) -> dict:
    blocks = get_room_blocks(scene)
    return {
        'Lx': max(block['origin']['x'] + block['size']['dx'] for block in blocks),
        'Ly': max(block['origin']['y'] + block['size']['dy'] for block in blocks),
        'Lz': max(block['origin']['z'] + block['size']['dz'] for block in blocks),
    }


def obstacle_supporting_block(obstacle: dict, blocks: list[dict]) -> dict | None:
    x0, y0 = obstacle['min']['x'], obstacle['min']['y']
    x1, y1 = x0 + obstacle['size']['dx'], y0 + obstacle['size']['dy']
    for block in blocks:
        ox, oy = block['origin']['x'], block['origin']['y']
        dx, dy = block['size']['dx'], block['size']['dy']
        if x0 >= ox - EPS and x1 <= ox + dx + EPS and y0 >= oy - EPS and y1 <= oy + dy + EPS:
            return block
    return None


def exposed_wall_segments(blocks: list[dict], wall: str) -> list[tuple[float, float]]:
    size = overall_room_size({'room': {'blocks': blocks}})
    segments: list[tuple[float, float]] = []
    if wall in {'west', 'east'}:
        plane_x = 0.0 if wall == 'west' else size['Lx']
        for block in blocks:
            bx0 = block['origin']['x']
            bx1 = bx0 + block['size']['dx']
            if abs((bx0 if wall == 'west' else bx1) - plane_x) < 1e-9:
                by0 = block['origin']['y']
                by1 = by0 + block['size']['dy']
                segments.append((by0, by1))
    elif wall in {'south', 'north'}:
        plane_y = 0.0 if wall == 'south' else size['Ly']
        for block in blocks:
            by0 = block['origin']['y']
            by1 = by0 + block['size']['dy']
            if abs((by0 if wall == 'south' else by1) - plane_y) < 1e-9:
                bx0 = block['origin']['x']
                bx1 = bx0 + block['size']['dx']
                segments.append((bx0, bx1))
    return segments


def longest_segment(segments: list[tuple[float, float]]) -> tuple[float, float]:
    if not segments:
        raise ValueError('No exposed wall segment found')
    return max(segments, key=lambda seg: seg[1] - seg[0])


def fit_opening_to_segment(length_target: float, seg: tuple[float, float], margin: float = OPENING_MARGIN) -> tuple[float, float]:
    seg_len = seg[1] - seg[0]
    usable = max(0.35, seg_len - 2 * margin)
    length = min(length_target, usable)
    center = 0.5 * (seg[0] + seg[1])
    return round(center, 3), round(length, 3)


def enforce_room_bounds(scene: dict, wall_margin: float) -> None:
    room = overall_room_size(scene)
    blocks = get_room_blocks(scene)
    for obs in scene['obstacles']:
        support = obstacle_supporting_block(obs, blocks)
        if support is None:
            support = max(blocks, key=lambda b: b['size']['dx'] * b['size']['dy'])
        sx = support['origin']['x']
        sy = support['origin']['y']
        ex = sx + support['size']['dx']
        ey = sy + support['size']['dy']
        obs['size']['dx'] = min(obs['size']['dx'], ex - sx - 2 * wall_margin)
        obs['size']['dy'] = min(obs['size']['dy'], ey - sy - 2 * wall_margin)
        obs['size']['dz'] = min(obs['size']['dz'], room['Lz'] - wall_margin)
        obs['min']['x'] = clamp(obs['min']['x'], sx + wall_margin, ex - wall_margin - obs['size']['dx'])
        obs['min']['y'] = clamp(obs['min']['y'], sy + wall_margin, ey - wall_margin - obs['size']['dy'])
        obs['min']['z'] = clamp(obs['min']['z'], 0.0, room['Lz'] - wall_margin - obs['size']['dz'])


def simplify_obstacles(scene: dict, wall_margin: float, max_count: int, max_row_length_ratio: float) -> None:
    room = overall_room_size(scene)
    obs = sorted(scene['obstacles'], key=lambda o: o['size']['dx'] * o['size']['dy'] * o['size']['dz'], reverse=True)
    obs = obs[:max_count]
    max_len_y = room['Ly'] * max_row_length_ratio
    max_len_x = room['Lx'] * max_row_length_ratio
    for o in obs:
        if o['size']['dy'] > max_len_y:
            o['size']['dy'] = max_len_y
        if o['size']['dx'] > max_len_x:
            o['size']['dx'] = max_len_x
        o['size']['dz'] = min(o['size']['dz'], room['Lz'] - wall_margin)
        o['size']['dx'] = max(o['size']['dx'], 0.6)
        o['size']['dy'] = max(o['size']['dy'], 0.8)
        o['size']['dz'] = max(o['size']['dz'], 1.5)
    scene['obstacles'] = obs


def push_apart_obstacles(scene: dict, clearance: float, max_iters: int = 20) -> None:
    room = overall_room_size(scene)
    blocks = get_room_blocks(scene)
    obs = scene['obstacles']
    for _ in range(max_iters):
        moved = False
        for i in range(len(obs)):
            for j in range(i + 1, len(obs)):
                a, b = obs[i], obs[j]
                if overlaps(a, b, margin=clearance):
                    ta, tb = a['size']['dx'] * a['size']['dy'], b['size']['dx'] * b['size']['dy']
                    target = a if ta < tb else b
                    other = b if target is a else a
                    support = obstacle_supporting_block(target, blocks) or max(blocks, key=lambda blk: blk['size']['dx'] * blk['size']['dy'])
                    other_max = box_max(other)
                    sx = support['origin']['x']
                    sy = support['origin']['y']
                    ex = sx + support['size']['dx']
                    ey = sy + support['size']['dy']
                    new_x = other_max['x'] + clearance
                    if new_x + target['size']['dx'] <= ex - clearance:
                        target['min']['x'] = new_x
                    else:
                        new_y = other_max['y'] + clearance
                        target['min']['y'] = min(max(new_y, sy + clearance), ey - clearance - target['size']['dy'])
                    target['min']['z'] = clamp(target['min']['z'], 0.0, room['Lz'] - clearance - target['size']['dz'])
                    moved = True
        if not moved:
            break


def enforce_clear_zones(scene: dict, opening_clearance: float) -> None:
    room = overall_room_size(scene)
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
    room = overall_room_size(scene)
    blocks = get_room_blocks(scene)
    openings = scene.get('openings', [])
    if len(openings) != 2:
        return
    inlet = next((o for o in openings if o['type'] == 'inlet'), openings[0])
    outlet = next((o for o in openings if o['type'] == 'outlet'), openings[-1])

    if wall == 'x':
        inlet['wall'] = 'west'
        outlet['wall'] = 'east'
        in_seg = longest_segment(exposed_wall_segments(blocks, 'west'))
        out_seg = longest_segment(exposed_wall_segments(blocks, 'east'))
        inlet_center_u, inlet_du = fit_opening_to_segment(size_du or 1.0, in_seg)
        outlet_center_u, outlet_du = fit_opening_to_segment(size_du or 1.0, out_seg)
        center_v = round(min(1.4, room['Lz'] * 0.5), 3)
        dv = round(min(size_dv or 0.8, room['Lz'] - 0.9), 3)
        inlet['center'] = {'u': inlet_center_u, 'v': center_v}
        outlet['center'] = {'u': outlet_center_u, 'v': center_v}
        inlet['size'] = {'du': inlet_du, 'dv': dv}
        outlet['size'] = {'du': outlet_du, 'dv': dv}
    elif wall == 'y':
        inlet['wall'] = 'south'
        outlet['wall'] = 'north'
        in_seg = longest_segment(exposed_wall_segments(blocks, 'south'))
        out_seg = longest_segment(exposed_wall_segments(blocks, 'north'))
        inlet_center_u, inlet_du = fit_opening_to_segment(size_du or 1.0, in_seg)
        outlet_center_u, outlet_du = fit_opening_to_segment(size_du or 1.0, out_seg)
        center_v = round(min(1.4, room['Lz'] * 0.5), 3)
        dv = round(min(size_dv or 0.8, room['Lz'] - 0.9), 3)
        inlet['center'] = {'u': inlet_center_u, 'v': center_v}
        outlet['center'] = {'u': outlet_center_u, 'v': center_v}
        inlet['size'] = {'du': inlet_du, 'dv': dv}
        outlet['size'] = {'du': outlet_du, 'dv': dv}


def repair_scene(scene: dict) -> tuple[dict, dict]:
    repaired = deepcopy(scene)
    room = overall_room_size(repaired)

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
        'Applied solver-friendly repair with composite-room aware opening normalization, '
        'obstacle count cap, clearance enforcement, and wall-margin bounds.'
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
