#!/usr/bin/env python3
"""
Post-process OpenFOAM boundary types after gmshToFoam.

gmshToFoam typically imports all surfaces as patch. For this case we want:
  - inlet/outlet : patch
  - topWall/bottomWall/obstacle : wall
  - front/back : symmetryPlane

Usage:
  python scripts/setup_openfoam.py
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent.parent


PATCH_TYPES: dict[str, str] = {
    "inlet": "patch",
    "outlet": "patch",
    "topWall": "wall",
    "bottomWall": "wall",
    "obstacle": "wall",
    "front": "symmetryPlane",
    "back": "symmetryPlane",
}


def _set_patch_type(boundary_text: str, patch_name: str, new_type: str) -> tuple[str, bool]:
    """
    Replace `type ...;` inside a patch block in constant/polyMesh/boundary.
    Returns (updated_text, changed?)
    """
    pat = re.compile(
        rf"(\b{re.escape(patch_name)}\b\s*\{{[^}}]*?\btype\s+)(\w+)\s*;",
        flags=re.MULTILINE,
    )
    m = pat.search(boundary_text)
    if not m:
        return boundary_text, False

    old_type = m.group(2)
    if old_type == new_type:
        return boundary_text, False

    updated = pat.sub(rf"\1{new_type};", boundary_text, count=1)
    return updated, True


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fix patch types in constant/polyMesh/boundary.")
    p.add_argument(
        "--boundary",
        type=str,
        default="",
        help="Path to boundary file (default: case/constant/polyMesh/boundary)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    boundary_path = Path(args.boundary) if args.boundary else (CASE_DIR / "constant" / "polyMesh" / "boundary")

    if not boundary_path.exists():
        print(f"[ERR] boundary file not found: {boundary_path}")
        print("      먼저 다음을 실행하세요:")
        print(f"      - cd {CASE_DIR}")
        print("      - python scripts/generate_mesh.py --output mesh.msh")
        print("      - gmshToFoam mesh.msh")
        return 2

    text = boundary_path.read_text(encoding="utf-8", errors="ignore")
    changed_any = False

    for patch_name, new_type in PATCH_TYPES.items():
        text, changed = _set_patch_type(text, patch_name, new_type)
        if changed:
            changed_any = True
            print(f"[OK] {patch_name}: type -> {new_type}")

    if changed_any:
        boundary_path.write_text(text, encoding="utf-8")
        print(f"[OK] Updated: {boundary_path}")
    else:
        print("[OK] No changes needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


