#!/usr/bin/env python3
"""Collect grid independence test results into a summary JSON."""

import json
import re
from pathlib import Path

PROJECT_ROOT = Path("/home/yhjoo/projects/OpenFOAM")
CASES_DIR = PROJECT_ROOT / "cases"
RESULTS_DIR = PROJECT_ROOT / "results"
LOG_DIR = PROJECT_ROOT / "benchmark/manifests/grid_test_logs"
OUTPUT = PROJECT_ROOT / "benchmark/manifests/phase1_grid_test_results.json"

CASES = ["a1_01", "a3_03", "a4_03"]
MESHES = [0.35, 0.25, 0.18, 0.10]
MESH_LABELS = ["035", "025", "018", "010"]
TYPES = ["ref", "pred"]


def extract_cell_count(case_dir: Path) -> int | None:
    owner = case_dir / "constant" / "polyMesh" / "owner"
    if not owner.exists():
        return None
    text = owner.read_text(errors="ignore")[:2000]
    m = re.search(r"nCells:\s*(\d+)", text)
    return int(m.group(1)) if m else None


def extract_from_solver_log(log_path: Path) -> dict:
    """Extract solve time, last iteration, and convergence from a solver log."""
    info = {"solve_time": None, "last_iter": None, "reached_end": False}
    if not log_path.exists():
        return info
    text = log_path.read_text(errors="ignore")
    # Last ExecutionTime
    for m in re.finditer(r"ExecutionTime\s*=\s*([0-9.]+)\s*s", text):
        info["solve_time"] = float(m.group(1))
    # Last Time = N
    for m in re.finditer(r"^Time = (\d+)", text, re.MULTILINE):
        info["last_iter"] = int(m.group(1))
    # Check if reached endTime (1000)
    if info["last_iter"] is not None and info["last_iter"] >= 1000:
        info["reached_end"] = True
    lower = text.lower()
    if "\nend\n" in lower or lower.rstrip().endswith("end"):
        info["reached_end"] = True
    return info


def extract_solve_time(results_dir: Path, case_dir: Path) -> tuple[float | None, str | None]:
    """Return (solve_time, source) from stabilization summary or raw log."""
    summary_path = results_dir / "stabilization_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        for attempt in summary.get("attempts", []):
            if attempt.get("reason") == "completed":
                log_path = Path(attempt.get("log", ""))
                if log_path.exists():
                    text = log_path.read_text(errors="ignore")
                    for line in reversed(text.strip().split("\n")[-20:]):
                        m = re.search(r"ExecutionTime\s*=\s*([0-9.]+)\s*s", line)
                        if m:
                            return float(m.group(1)), "stabilization_summary"
    # Fallback: read from case log.simpleFoam
    log_path = case_dir / "log.simpleFoam"
    info = extract_from_solver_log(log_path)
    if info["solve_time"] is not None:
        return info["solve_time"], "raw_log"
    return None, None


def extract_convergence(results_dir: Path, case_dir: Path) -> tuple[bool, str | None, str]:
    """Return (converged, preset, status)."""
    summary_path = results_dir / "stabilization_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        success = summary.get("success", False)
        preset = summary.get("successful_preset")
        return success, preset, "completed" if success else "failed"
    # No summary = pipeline crashed (likely timeout)
    log_path = case_dir / "log.simpleFoam"
    info = extract_from_solver_log(log_path)
    if info["reached_end"]:
        return True, "unknown", "completed_no_summary"
    if info["last_iter"] is not None:
        return False, None, f"solver_timeout_at_iter_{info['last_iter']}"
    return False, None, "no_data"


def main():
    results = {
        "description": "Phase 1 Grid Independence Test Results",
        "cases": CASES,
        "mesh_sizes": MESHES,
        "types": ["reference", "predicted"],
        "solver_timeout": 3600,
        "runs": [],
        "summary_table": [],
    }

    for case_id in CASES:
        for mi, (mesh, ml) in enumerate(zip(MESHES, MESH_LABELS)):
            for run_type in TYPES:
                run_name = f"grid_test_{case_id}_{run_type}_{ml}"
                case_dir = CASES_DIR / run_name
                res_dir = RESULTS_DIR / run_name

                cell_count = extract_cell_count(case_dir)
                solve_time, solve_source = extract_solve_time(res_dir, case_dir)
                converged, preset, status = extract_convergence(res_dir, case_dir)

                # Get last iteration from solver log
                log_path = case_dir / "log.simpleFoam"
                log_info = extract_from_solver_log(log_path)

                # Get elapsed wall time from log
                elapsed_file = LOG_DIR / f"{run_name}.elapsed"
                wall_time = None
                if elapsed_file.exists():
                    try:
                        wall_time = int(elapsed_file.read_text().strip())
                    except ValueError:
                        pass

                rc_file = LOG_DIR / f"{run_name}.rc"
                return_code = None
                if rc_file.exists():
                    try:
                        return_code = int(rc_file.read_text().strip())
                    except ValueError:
                        pass

                entry = {
                    "case": f"bench_{case_id}",
                    "case_id": case_id,
                    "mesh_size": mesh,
                    "mesh_label": ml,
                    "run_type": "reference" if run_type == "ref" else "predicted",
                    "run_name": run_name,
                    "cell_count": cell_count,
                    "converged": converged,
                    "status": status,
                    "successful_preset": preset,
                    "solve_time_s": solve_time,
                    "solve_time_source": solve_source,
                    "last_iteration": log_info["last_iter"],
                    "wall_time_s": wall_time,
                    "pipeline_rc": return_code,
                }
                results["runs"].append(entry)

                results["summary_table"].append({
                    "case": case_id,
                    "type": run_type,
                    "mesh": mesh,
                    "cells": cell_count,
                    "converged": "yes" if converged else "no",
                    "status": status,
                    "last_iter": log_info["last_iter"],
                    "solve_s": solve_time,
                    "wall_s": wall_time,
                })

    # Summary stats
    total = len(results["runs"])
    completed_pipeline = sum(1 for r in results["runs"] if r["pipeline_rc"] == 0)
    timed_out = sum(1 for r in results["runs"] if r["pipeline_rc"] == 1)
    converged_count = sum(1 for r in results["runs"] if r["converged"])
    has_data = sum(1 for r in results["runs"] if r["cell_count"] is not None)
    results["stats"] = {
        "total_runs": total,
        "pipeline_success": completed_pipeline,
        "pipeline_timeout": timed_out,
        "converged": converged_count,
        "has_mesh_data": has_data,
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Written {OUTPUT}")
    print(f"Stats: {completed_pipeline}/{total} pipeline success, {timed_out} timed out, {converged_count} converged")

    # Print table
    print(f"\n{'Case':<10} {'Type':<6} {'Mesh':<6} {'Cells':>8} {'Conv':>5} {'Iter':>6} {'Solve(s)':>10} {'Wall(s)':>8} {'Status'}")
    print("-" * 85)
    for r in results["summary_table"]:
        cells = str(r["cells"]) if r["cells"] else "---"
        solve = f"{r['solve_s']:.1f}" if r["solve_s"] else "---"
        wall = str(r["wall_s"]) if r["wall_s"] else "---"
        last_iter = str(r["last_iter"]) if r["last_iter"] else "---"
        print(f"{r['case']:<10} {r['type']:<6} {r['mesh']:<6} {cells:>8} {r['converged']:>5} {last_iter:>6} {solve:>10} {wall:>8} {r['status']}")


if __name__ == "__main__":
    main()
