#!/usr/bin/env python3
"""Batch runner for scaffolded benchmark image-conditioned evaluation tasks.

Wraps `run_benchmark_evaluation_task.py` so a small filtered subset (for example
one case across all 5 views, or a category subset) can be executed and logged in
one repeatable command.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK = PROJECT_ROOT / "benchmark"
EVALUATIONS = BENCHMARK / "evaluations"
SCRIPTS = PROJECT_ROOT / "scripts"
TASK_RUNNER = SCRIPTS / "run_benchmark_evaluation_task.py"
DEFAULT_VIEWS = ["perspective", "birdseye", "floorplan", "wireframe", "section"]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def selected_tasks(
    *,
    evaluation_root: Path,
    cases: list[str] | None,
    categories: list[str] | None,
    views: list[str] | None,
    statuses: list[str] | None,
    limit: int | None,
) -> list[Path]:
    task_paths = sorted(evaluation_root.glob("*/*/task.json"))
    picked: list[Path] = []
    for task_path in task_paths:
        task = load_json(task_path)
        if cases and task.get("case_name") not in cases:
            continue
        if categories and task.get("category") not in categories:
            continue
        if views and task.get("view") not in views:
            continue
        if statuses and task.get("status") not in statuses:
            continue
        picked.append(task_path)
        if limit is not None and len(picked) >= limit:
            break
    return picked


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a filtered batch of benchmark evaluation tasks")
    parser.add_argument("--evaluation-root", type=Path, default=EVALUATIONS)
    parser.add_argument("--cases", nargs="*", default=None, help="Case names to run, e.g. bench_a1_01")
    parser.add_argument("--categories", nargs="*", default=None, help="Categories to run, e.g. A1 A3")
    parser.add_argument("--views", nargs="*", choices=DEFAULT_VIEWS, default=None, help="View types to run")
    parser.add_argument("--statuses", nargs="*", default=["pending"], help="Only run tasks currently in these statuses")
    parser.add_argument("--backend", choices=["cli", "api"], default="cli")
    parser.add_argument("--model", default="gemini-3-flash-preview")
    parser.add_argument("--mesh-size", type=float, default=0.35)
    parser.add_argument("--end-time", type=int, default=250)
    parser.add_argument("--solver-timeout", type=int, default=600)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--skip-existing-success", action="store_true")
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=BENCHMARK / "manifests" / "evaluation_batch_summary.json",
        help="Batch summary output JSON",
    )
    args = parser.parse_args()

    evaluation_root = args.evaluation_root.expanduser().resolve()
    task_paths = selected_tasks(
        evaluation_root=evaluation_root,
        cases=args.cases,
        categories=args.categories,
        views=args.views,
        statuses=args.statuses,
        limit=args.limit,
    )

    results: list[dict[str, Any]] = []
    for task_path in task_paths:
        task = load_json(task_path)
        cmd = [
            "python3",
            str(TASK_RUNNER),
            "--task",
            str(task_path),
            "--backend",
            args.backend,
            "--model",
            args.model,
            "--mesh-size",
            str(args.mesh_size),
            "--end-time",
            str(args.end_time),
            "--solver-timeout",
            str(args.solver_timeout),
        ]
        if args.skip_existing_success:
            cmd.append("--skip-existing-success")
        proc = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)
        refreshed_task = load_json(task_path)
        summary_path = Path(refreshed_task.get("evaluation_summary")) if refreshed_task.get("evaluation_summary") else None
        evaluation_summary = load_json(summary_path) if summary_path and summary_path.exists() else None
        cfd_score = None
        structural_score = None
        if evaluation_summary:
            structural_score = (((evaluation_summary.get("prediction_summary") or {}).get("structural_score")))
            cfd_score = (((evaluation_summary.get("cfd_summary") or {}).get("aggregate_score") or {}).get("cfd_score"))
        results.append(
            {
                "task": str(task_path),
                "case_name": task.get("case_name"),
                "category": task.get("category"),
                "view": task.get("view"),
                "returncode": proc.returncode,
                "status": refreshed_task.get("status"),
                "last_run_name": refreshed_task.get("last_run_name"),
                "evaluation_summary": str(summary_path) if summary_path else None,
                "structural_score": structural_score,
                "cfd_score": cfd_score,
                "stdout_tail": "\n".join((proc.stdout or "").splitlines()[-30:]),
                "stderr_tail": "\n".join((proc.stderr or "").splitlines()[-30:]),
            }
        )

    payload = {
        "ok": True,
        "evaluation_root": str(evaluation_root),
        "backend": args.backend,
        "model": args.model,
        "filters": {
            "cases": args.cases,
            "categories": args.categories,
            "views": args.views,
            "statuses": args.statuses,
            "limit": args.limit,
        },
        "count": len(results),
        "status_counts": {
            status: sum(1 for row in results if row.get("status") == status)
            for status in sorted({row.get("status") for row in results})
        },
        "results": results,
    }
    write_json(args.summary_out, payload)
    print(json.dumps(payload, indent=2))
    return 0 if all(row.get("status") in {"success", "blocked"} for row in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
