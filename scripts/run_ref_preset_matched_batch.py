#!/usr/bin/env python3
"""Re-run reference cases with preset matched to their predicted counterpart.

For each (case, view) pair, reads the predicted case's successful preset,
then re-runs the reference case with --force-preset to match inlet velocity.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
SCENES_DIR = PROJECT_ROOT / "benchmark" / "scenes"
RESULTS = PROJECT_ROOT / "results"
CASES = PROJECT_ROOT / "cases"
PYTHON = sys.executable
PIPELINE = PROJECT_ROOT / "scripts" / "run_indoor_stabilized.py"

# 3 solver-failed predicted cases
SOLVER_FAILED = {
    ("bench_a3_04", "perspective"),
    ("bench_a4_02", "perspective"),
    ("bench_a4_02", "wireframe"),
}


def load_json(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def get_pred_preset(case: str, view: str) -> str | None:
    """Get the preset used by the predicted case."""
    pred_dir = RESULTS / f"phase2_pred_{case}_{view}"
    summary = pred_dir / "stabilization_summary.json"
    if summary.exists():
        return load_json(summary).get("successful_preset")
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", "-w", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--mesh-size", type=float, default=0.18)
    parser.add_argument("--solver-timeout", type=int, default=3600)
    parser.add_argument("--import-timeout", type=int, default=600)
    args = parser.parse_args()

    # Collect all (scene_id, preset) pairs needed
    jobs: list[dict] = []
    # Track unique (scene_id, preset) to avoid duplicate ref runs
    ref_runs_needed: dict[tuple[str, str], list[dict]] = defaultdict(list)

    for sp in sorted(EVAL_ROOT.glob("*/*/evaluation_summary.json")):
        d = load_json(sp)
        task = d.get("task", {})
        case = task.get("case_name", "")
        view = task.get("view", "")

        if (case, view) in SOLVER_FAILED:
            continue

        scene_id = case.replace("bench_", "")
        pred_preset = get_pred_preset(case, view)
        if not pred_preset:
            continue

        # Check if ref already ran with this preset
        ref_dir = CASES / f"phase2_ref_{scene_id}_preset_{pred_preset}"
        if ref_dir.exists() and (ref_dir / "VTK").exists():
            continue  # Already done

        ref_runs_needed[(scene_id, pred_preset)].append({"case": case, "view": view})

    print(f"Unique (scene, preset) pairs to run: {len(ref_runs_needed)}")
    for (sid, preset), views in sorted(ref_runs_needed.items()):
        view_list = [v["view"] for v in views]
        print(f"  {sid} with preset={preset}: needed by {len(views)} views ({', '.join(view_list[:3])}...)")

    if args.dry_run:
        return

    # Run each unique (scene_id, preset)
    total = len(ref_runs_needed)
    done = 0
    failed = 0

    for (scene_id, preset), views in sorted(ref_runs_needed.items()):
        scene_json = SCENES_DIR / f"{scene_id}.json"
        case_name = f"phase2_ref_{scene_id}_preset_{preset}"

        print(f"\n[{done+1}/{total}] {case_name} (force-preset={preset})")
        t0 = time.monotonic()

        cmd = [
            PYTHON, str(PIPELINE),
            "--scenario", str(scene_json),
            "--name", case_name,
            "--mesh-size", str(args.mesh_size),
            "--skip-mesh-ladder",
            "--solver-timeout", str(args.solver_timeout),
            "--import-timeout", str(args.import_timeout),
            "--disable-repair",
            "--force-preset", preset,
        ]

        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT),
                                capture_output=True, text=True, timeout=args.solver_timeout + 300)
        elapsed = time.monotonic() - t0

        # Check success
        case_dir = CASES / case_name
        vtk_exists = (case_dir / "VTK").exists() and any((case_dir / "VTK").iterdir()) if (case_dir / "VTK").exists() else False

        if result.returncode == 0 or vtk_exists:
            print(f"  OK ({elapsed:.0f}s)")
            done += 1

            # Run foamToVTK if not already done
            if not vtk_exists:
                subprocess.run([
                    "docker", "run", "--rm",
                    "-v", f"{PROJECT_ROOT}:/app",
                    "-w", f"/app/cases/{case_name}",
                    "openfoam-pipeline-local:latest",
                    "bash", "-lc", "foamToVTK -latestTime"
                ], capture_output=True, timeout=300)
        else:
            print(f"  FAILED ({elapsed:.0f}s) rc={result.returncode}")
            failed += 1

    print(f"\n{'='*50}")
    print(f"Complete: {done}/{total} OK, {failed} failed")


if __name__ == "__main__":
    raise SystemExit(main())
