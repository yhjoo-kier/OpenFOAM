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


def interval_gap(a0: float, a1: float, b0: float, b1: float) -> float:
    if a1 < b0:
        return b0 - a1
    if b1 < a0:
        return a0 - b1
    return 0.0


def interval_overlap(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))


def needs_obstacle_separation(a: dict, b: dict, clearance: float) -> bool:
    amax = box_max(a)
    bmax = box_max(b)
    x_gap = interval_gap(a['min']['x'], amax['x'], b['min']['x'], bmax['x'])
    y_gap = interval_gap(a['min']['y'], amax['y'], b['min']['y'], bmax['y'])
    z_gap = interval_gap(a['min']['z'], amax['z'], b['min']['z'], bmax['z'])
    x_overlap = interval_overlap(a['min']['x'], amax['x'], b['min']['x'], bmax['x'])
    y_overlap = interval_overlap(a['min']['y'], amax['y'], b['min']['y'], bmax['y'])
    z_overlap = interval_overlap(a['min']['z'], amax['z'], b['min']['z'], bmax['z'])
    return (
        (x_gap < clearance - EPS and y_overlap > EPS and z_overlap > EPS)
        or (y_gap < clearance - EPS and x_overlap > EPS and z_overlap > EPS)
        or (z_gap < clearance - EPS and x_overlap > EPS and y_overlap > EPS)
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


def normalize_room_origin(scene: dict) -> dict[str, float]:
    blocks = get_room_blocks(scene)
    min_x = min(block['origin']['x'] for block in blocks)
    min_y = min(block['origin']['y'] for block in blocks)
    min_z = min(block['origin']['z'] for block in blocks)
    shift = {
        'x': -min_x if min_x < 0 else 0.0,
        'y': -min_y if min_y < 0 else 0.0,
        'z': -min_z if min_z < 0 else 0.0,
    }
    if shift['x'] == 0.0 and shift['y'] == 0.0 and shift['z'] == 0.0:
        return shift

    for block in blocks:
        block['origin']['x'] += shift['x']
        block['origin']['y'] += shift['y']
        block['origin']['z'] += shift['z']
    for obs in scene.get('obstacles', []):
        obs['min']['x'] += shift['x']
        obs['min']['y'] += shift['y']
        obs['min']['z'] += shift['z']
    return shift


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


def total_segment_length(segments: list[tuple[float, float]]) -> float:
    return sum(max(0.0, seg[1] - seg[0]) for seg in segments)


def infer_repair_opening_axis(scene: dict) -> str:
    blocks = get_room_blocks(scene)
    openings = scene.get('openings', [])
    inlet = next((o for o in openings if o.get('type') == 'inlet'), openings[0] if openings else None)
    west_segments = exposed_wall_segments(blocks, 'west')
    east_segments = exposed_wall_segments(blocks, 'east')
    south_segments = exposed_wall_segments(blocks, 'south')
    north_segments = exposed_wall_segments(blocks, 'north')
    x_available = bool(west_segments) and bool(east_segments)
    y_available = bool(south_segments) and bool(north_segments)
    axis_counts = {
        'x': sum(1 for op in openings if op.get('wall') in {'west', 'east'}),
        'y': sum(1 for op in openings if op.get('wall') in {'south', 'north'}),
    }

    if axis_counts['x'] > axis_counts['y'] and x_available:
        return 'x'
    if axis_counts['y'] > axis_counts['x'] and y_available:
        return 'y'
    if inlet and inlet.get('wall') in {'west', 'east'} and x_available:
        return 'x'
    if inlet and inlet.get('wall') in {'south', 'north'} and y_available:
        return 'y'

    candidates: list[tuple[str, float]] = []
    if x_available:
        candidates.append(('x', total_segment_length(west_segments) + total_segment_length(east_segments)))
    if y_available:
        candidates.append(('y', total_segment_length(south_segments) + total_segment_length(north_segments)))
    if candidates:
        return max(candidates, key=lambda item: item[1])[0]
    return 'x'


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
                if not needs_obstacle_separation(a, b, clearance):
                    continue
                ta, tb = a['size']['dx'] * a['size']['dy'], b['size']['dx'] * b['size']['dy']
                target = a if ta < tb else b
                other = b if target is a else a
                support = obstacle_supporting_block(target, blocks) or max(blocks, key=lambda blk: blk['size']['dx'] * blk['size']['dy'])
                other_max = box_max(other)
                sx = support['origin']['x']
                sy = support['origin']['y']
                ex = sx + support['size']['dx']
                ey = sy + support['size']['dy']
                candidate_positions: list[tuple[float, float]] = []
                candidate_positions.append((other['min']['x'] - clearance - target['size']['dx'], target['min']['y']))
                candidate_positions.append((other_max['x'] + clearance, target['min']['y']))
                candidate_positions.append((target['min']['x'], other['min']['y'] - clearance - target['size']['dy']))
                candidate_positions.append((target['min']['x'], other_max['y'] + clearance))

                valid_candidates: list[tuple[float, float, float]] = []
                for cand_x, cand_y in candidate_positions:
                    if cand_x < sx + clearance - EPS or cand_x + target['size']['dx'] > ex - clearance + EPS:
                        continue
                    if cand_y < sy + clearance - EPS or cand_y + target['size']['dy'] > ey - clearance + EPS:
                        continue
                    dx_move = cand_x - target['min']['x']
                    dy_move = cand_y - target['min']['y']
                    trial = deepcopy(target)
                    trial['min']['x'] = cand_x
                    trial['min']['y'] = cand_y
                    if needs_obstacle_separation(trial, other, clearance):
                        continue
                    valid_candidates.append((dx_move * dx_move + dy_move * dy_move, cand_x, cand_y))

                if valid_candidates:
                    _, cand_x, cand_y = min(valid_candidates, key=lambda item: item[0])
                    target['min']['x'] = cand_x
                    target['min']['y'] = cand_y
                else:
                    new_x = other['min']['x'] - clearance - target['size']['dx']
                    if new_x >= sx + clearance - EPS:
                        target['min']['x'] = new_x
                    else:
                        new_y = other['min']['y'] - clearance - target['size']['dy']
                        if new_y >= sy + clearance - EPS:
                            target['min']['y'] = new_y
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


def normalize_openings(scene: dict, wall: str | None = None, size_du: float | None = None, size_dv: float | None = None) -> str | None:
    room = overall_room_size(scene)
    blocks = get_room_blocks(scene)
    openings = scene.get('openings', [])
    if len(openings) != 2:
        return None
    inlet = next((o for o in openings if o['type'] == 'inlet'), openings[0])
    outlet = next((o for o in openings if o['type'] == 'outlet'), openings[-1])
    wall = wall or infer_repair_opening_axis(scene)

    if wall == 'x':
        inlet_wall = inlet.get('wall') if inlet.get('wall') in {'west', 'east'} else 'west'
        outlet_wall = 'east' if inlet_wall == 'west' else 'west'
        inlet['wall'] = inlet_wall
        outlet['wall'] = outlet_wall
        in_seg = longest_segment(exposed_wall_segments(blocks, inlet_wall))
        out_seg = longest_segment(exposed_wall_segments(blocks, outlet_wall))
        inlet_center_u, inlet_du = fit_opening_to_segment(size_du or 1.0, in_seg)
        outlet_center_u, outlet_du = fit_opening_to_segment(size_du or 1.0, out_seg)
        center_v = round(min(1.4, room['Lz'] * 0.5), 3)
        dv = round(min(size_dv or 0.8, room['Lz'] - 0.9), 3)
        inlet['center'] = {'u': inlet_center_u, 'v': center_v}
        outlet['center'] = {'u': outlet_center_u, 'v': center_v}
        inlet['size'] = {'du': inlet_du, 'dv': dv}
        outlet['size'] = {'du': outlet_du, 'dv': dv}
        return wall
    if wall == 'y':
        inlet_wall = inlet.get('wall') if inlet.get('wall') in {'south', 'north'} else 'south'
        outlet_wall = 'north' if inlet_wall == 'south' else 'south'
        inlet['wall'] = inlet_wall
        outlet['wall'] = outlet_wall
        in_seg = longest_segment(exposed_wall_segments(blocks, inlet_wall))
        out_seg = longest_segment(exposed_wall_segments(blocks, outlet_wall))
        inlet_center_u, inlet_du = fit_opening_to_segment(size_du or 1.0, in_seg)
        outlet_center_u, outlet_du = fit_opening_to_segment(size_du or 1.0, out_seg)
        center_v = round(min(1.4, room['Lz'] * 0.5), 3)
        dv = round(min(size_dv or 0.8, room['Lz'] - 0.9), 3)
        inlet['center'] = {'u': inlet_center_u, 'v': center_v}
        outlet['center'] = {'u': outlet_center_u, 'v': center_v}
        inlet['size'] = {'du': inlet_du, 'dv': dv}
        outlet['size'] = {'du': outlet_du, 'dv': dv}
        return wall
    return None


def repair_scene(scene: dict) -> tuple[dict, dict]:
    repaired = deepcopy(scene)
    origin_shift = normalize_room_origin(repaired)
    room = overall_room_size(repaired)

    wall_margin = max(0.4, 0.05 * min(room['Lx'], room['Ly']))
    obstacle_clearance = max(0.6, 0.08 * min(room['Lx'], room['Ly']))
    opening_clearance = max(1.2, 0.12 * room['Lx'])

    simplify_obstacles(repaired, wall_margin=wall_margin, max_count=3, max_row_length_ratio=0.62)
    chosen_opening_axis = normalize_openings(repaired, wall=None, size_du=1.0, size_dv=0.8)
    enforce_room_bounds(repaired, wall_margin=wall_margin)
    enforce_clear_zones(repaired, opening_clearance=opening_clearance)
    push_apart_obstacles(repaired, clearance=obstacle_clearance)
    enforce_room_bounds(repaired, wall_margin=wall_margin)

    repaired.setdefault('meta', {})
    repaired['meta']['repair_applied'] = True
    repaired['meta']['repair_notes'] = (
        'Applied solver-friendly repair with room-origin normalization, composite-room aware opening normalization, '
        'obstacle count cap, clearance enforcement, and wall-margin bounds.'
    )
    repaired['meta']['repair_opening_axis'] = chosen_opening_axis
    repaired['meta']['repair_origin_shift'] = origin_shift

    report = validate_scene(repaired)
    info = {
        'wall_margin': wall_margin,
        'obstacle_clearance': obstacle_clearance,
        'opening_clearance': opening_clearance,
        'repair_opening_axis': chosen_opening_axis,
        'origin_shift': origin_shift,
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
