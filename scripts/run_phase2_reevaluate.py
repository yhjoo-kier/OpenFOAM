#!/usr/bin/env python3
"""Phase 2 re-evaluation: recompute CFD metrics with improved metric script.

Updates symlinks to Phase 2 case dirs, re-runs CFD metrics, updates
evaluation_summary.json, and rebuilds aggregate summaries.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EVAL_ROOT = PROJECT_ROOT / "benchmark" / "evaluations_posthoc_scaled_longest_span"
CASES_DIR = PROJECT_ROOT / "cases"
SCENES_DIR = PROJECT_ROOT / "benchmark" / "scenes"
SCRIPTS = PROJECT_ROOT / "scripts"
PYTHON = sys.executable

# 3 solver-diverged cases to skip (Option B: report as failures)
SOLVER_FAILED = {
    ("bench_a3_04", "perspective"),
    ("bench_a4_02", "perspective"),
    ("bench_a4_02", "wireframe"),
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def update_symlink(target: Path, link: Path) -> None:
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(target)


def find_phase2_ref_case(scene_id: str) -> Path | None:
    """Find Phase 2 reference case dir (e.g., cases/phase2_ref_a1_01)."""
    case_dir = CASES_DIR / f"phase2_ref_{scene_id}"
    if case_dir.exists():
        return case_dir
    return None


def find_phase2_pred_case(case_name: str, view: str) -> Path | None:
    """Find Phase 2 predicted case dir (e.g., cases/phase2_pred_bench_a1_01_floorplan)."""
    case_dir = CASES_DIR / f"phase2_pred_{case_name}_{view}"
    if case_dir.exists():
        return case_dir
    return None


def has_vtk(case_dir: Path) -> bool:
    vtk_dir = case_dir / "VTK"
    return vtk_dir.exists() and any(vtk_dir.iterdir())


def run_cfd_metrics(
    ref_scene: Path, ref_case: Path,
    pred_scene: Path, pred_case: Path,
    output: Path,
) -> dict[str, Any] | None:
    """Run improved CFD metrics script."""
    cmd = [
        PYTHON, str(SCRIPTS / "compute_benchmark_cfd_metrics.py"),
        "--reference-scene", str(ref_scene),
        "--reference-case", str(ref_case),
        "--predicted-scene", str(pred_scene),
        "--predicted-case", str(pred_case),
        "-o", str(output),
    ]
    proc = subprocess.run(cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        print(f"  [WARN] metrics failed: {proc.stderr[-200:]}")
        return None
    if output.exists():
        return load_json(output)
    return None


def update_evaluation_summary(eval_dir: Path, cfd_metrics: dict[str, Any]) -> None:
    """Update evaluation_summary.json with new CFD metrics."""
    summary_path = eval_dir / "evaluation_summary.json"
    if not summary_path.exists():
        return
    summary = load_json(summary_path)
    summary["cfd_summary"] = cfd_metrics
    # Update top-level cfd_score for backward compat
    agg = cfd_metrics.get("aggregate_score", {})
    if "pipeline_summary" in summary:
        summary["pipeline_summary"]["cfd_score"] = agg.get("cfd_score")
    write_json(summary_path, summary)


def main() -> int:
    parser = argparse.ArgumentParser(description="Phase 2 re-evaluation with improved metrics")
    parser.add_argument("--workers", "-w", type=int, default=1,
                        help="Number of parallel workers (default: 1, sequential)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without executing")
    args = parser.parse_args()

    eval_dirs = sorted(EVAL_ROOT.iterdir())
    total = 0
    success = 0
    skipped = 0
    failed = 0
    solver_skipped = 0

    for case_dir in eval_dirs:
        if not case_dir.is_dir() or case_dir.name.endswith(".json"):
            continue
        case_name = case_dir.name  # e.g., bench_a1_01
        scene_id = case_name.replace("bench_", "")  # e.g., a1_01

        for view_dir in sorted(case_dir.iterdir()):
            if not view_dir.is_dir():
                continue
            view = view_dir.name  # e.g., floorplan
            total += 1

            # Skip solver-diverged cases
            if (case_name, view) in SOLVER_FAILED:
                solver_skipped += 1
                print(f"[SKIP-SOLVER] {case_name}/{view} (solver diverged)")
                continue

            # Find Phase 2 cases
            ref_case = find_phase2_ref_case(scene_id)
            pred_case = find_phase2_pred_case(case_name, view)

            if ref_case is None or not has_vtk(ref_case):
                print(f"[SKIP] {case_name}/{view}: no Phase 2 reference VTK")
                skipped += 1
                continue
            if pred_case is None or not has_vtk(pred_case):
                print(f"[SKIP] {case_name}/{view}: no Phase 2 predicted VTK")
                skipped += 1
                continue

            # Scene paths
            ref_scene = SCENES_DIR / f"{scene_id}.json"
            pred_scene = view_dir / "scaled_scene.json"
            if not pred_scene.exists():
                pred_scene = view_dir / "predicted_scene.json"
            if not ref_scene.exists() or not pred_scene.exists():
                print(f"[SKIP] {case_name}/{view}: missing scene JSON")
                skipped += 1
                continue

            if args.dry_run:
                print(f"[DRY] {case_name}/{view}: ref={ref_case.name} pred={pred_case.name}")
                continue

            # Update symlinks
            update_symlink(ref_case, view_dir / "reference_case")
            update_symlink(pred_case, view_dir / "predicted_case")
            ref_results = PROJECT_ROOT / "results" / ref_case.name
            pred_results = PROJECT_ROOT / "results" / pred_case.name
            if ref_results.exists():
                update_symlink(ref_results, view_dir / "reference_results")
            if pred_results.exists():
                update_symlink(pred_results, view_dir / "predicted_results")

            # Run CFD metrics
            cfd_output = view_dir / "cfd_metrics.json"
            print(f"[EVAL] {case_name}/{view} ...", end=" ", flush=True)
            t0 = time.monotonic()
            cfd_metrics = run_cfd_metrics(ref_scene, ref_case, pred_scene, pred_case, cfd_output)
            elapsed = time.monotonic() - t0

            if cfd_metrics and cfd_metrics.get("ok"):
                score = cfd_metrics.get("aggregate_score", {}).get("cfd_agreement_score", "?")
                print(f"OK ({elapsed:.1f}s) agreement={score}")
                update_evaluation_summary(view_dir, cfd_metrics)
                success += 1
            else:
                print(f"FAIL ({elapsed:.1f}s)")
                failed += 1

    print(f"\n{'='*60}")
    print(f"Phase 2 Re-evaluation Complete")
    print(f"  Total:          {total}")
    print(f"  Success:        {success}")
    print(f"  Failed:         {failed}")
    print(f"  Skipped:        {skipped}")
    print(f"  Solver-skipped: {solver_skipped}")
    print(f"{'='*60}")

    # Rebuild aggregate
    if not args.dry_run and success > 0:
        print("\nRebuilding aggregate summaries...")
        agg_cmd = [
            PYTHON, str(SCRIPTS / "build_benchmark_evaluation_aggregate.py"),
            "--evaluation-root", str(EVAL_ROOT),
            "--json-output", str(PROJECT_ROOT / "benchmark" / "manifests" / "evaluation_aggregate_summary_phase2.json"),
            "--markdown-output", str(PROJECT_ROOT / "docs" / "26-03-22_phase2_aggregate_results.md"),
        ]
        proc = subprocess.run(agg_cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True)
        if proc.returncode == 0:
            print("Aggregate rebuild OK")
            # Print headline
            try:
                agg = json.loads(proc.stdout)
                print(f"  Tasks: {agg.get('task_count')}")
                print(f"  Mean structural: {agg.get('mean_structural_score')}")
                print(f"  Mean CFD:        {agg.get('mean_cfd_score')}")
            except (json.JSONDecodeError, KeyError):
                pass
        else:
            print(f"Aggregate rebuild FAILED: {proc.stderr[-300:]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
