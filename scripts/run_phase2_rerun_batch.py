#!/usr/bin/env python3
"""Phase 2 benchmark re-run: parallel batch orchestration.

Runs reference (20 scenes) and/or predicted (100 evaluation cases) CFD pipelines
in parallel using subprocess workers, with resume support and JSON logging.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_DIR = PROJECT_ROOT / "benchmark"
SCENES_DIR = BENCHMARK_DIR / "scenes"
EVAL_DIR = BENCHMARK_DIR / "evaluations_posthoc_scaled_longest_span"
RUN_SCRIPT = PROJECT_ROOT / "scripts" / "run_indoor_stabilized.py"
LOG_PATH = BENCHMARK_DIR / "manifests" / "phase2_run_log.json"


@dataclass
class Job:
    """A single CFD run job."""
    job_id: str
    phase: str  # "ref" or "pred"
    scenario_path: str
    case_name: str
    # tracking
    status: str = "pending"  # pending, running, completed, failed
    returncode: Optional[int] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    duration_s: Optional[float] = None

    def to_dict(self) -> dict:
        return {
            "job_id": self.job_id,
            "phase": self.phase,
            "scenario_path": self.scenario_path,
            "case_name": self.case_name,
            "status": self.status,
            "returncode": self.returncode,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_s": self.duration_s,
        }


@dataclass
class RunningWorker:
    """Tracks a running subprocess."""
    job: Job
    process: subprocess.Popen
    started: float  # monotonic clock


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_ref_jobs() -> list[Job]:
    """Build job list for reference geometry runs (20 cases)."""
    jobs = []
    for scene_file in sorted(SCENES_DIR.glob("*.json")):
        case_id = scene_file.stem  # e.g. a1_01
        jobs.append(Job(
            job_id=f"ref_{case_id}",
            phase="ref",
            scenario_path=str(scene_file),
            case_name=f"phase2_ref_{case_id}",
        ))
    return jobs


def build_pred_jobs() -> list[Job]:
    """Build job list for predicted case runs (100 cases)."""
    jobs = []
    for case_dir in sorted(EVAL_DIR.iterdir()):
        if not case_dir.is_dir() or case_dir.name.endswith(".json"):
            continue
        case_id = case_dir.name  # e.g. bench_a1_01
        for view_dir in sorted(case_dir.iterdir()):
            if not view_dir.is_dir():
                continue
            view = view_dir.name  # e.g. birdseye
            # Prefer scaled_scene.json, fall back to predicted_scene.json
            scene_json = view_dir / "scaled_scene.json"
            if not scene_json.exists():
                scene_json = view_dir / "predicted_scene.json"
            if not scene_json.exists():
                print(f"[WARN] No scene JSON found in {view_dir}, skipping")
                continue
            jobs.append(Job(
                job_id=f"pred_{case_id}_{view}",
                phase="pred",
                scenario_path=str(scene_json),
                case_name=f"phase2_pred_{case_id}_{view}",
            ))
    return jobs


def is_completed(job: Job) -> bool:
    """Check if a job already has VTK output (for --resume)."""
    case_dir = PROJECT_ROOT / "cases" / job.case_name
    vtk_dir = case_dir / "VTK"
    if vtk_dir.exists() and any(vtk_dir.iterdir()):
        return True
    # Also check results summary
    results_dir = PROJECT_ROOT / "results" / job.case_name
    summary = results_dir / "stabilization_summary.json"
    if summary.exists():
        try:
            data = json.loads(summary.read_text())
            return data.get("success", False)
        except (json.JSONDecodeError, OSError):
            pass
    return False


def launch_worker(job: Job, mesh_size: float, solver_timeout: int) -> RunningWorker:
    """Launch a single CFD worker subprocess."""
    cmd = [
        sys.executable, str(RUN_SCRIPT),
        "--scenario", job.scenario_path,
        "--name", job.case_name,
        "--mesh-size", str(mesh_size),
        "--skip-mesh-ladder",
        "--solver-timeout", str(solver_timeout),
        "--disable-repair",
    ]
    job.status = "running"
    job.start_time = _now_iso()

    log_dir = PROJECT_ROOT / "results" / job.case_name
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_log = log_dir / "batch_stdout.log"
    stderr_log = log_dir / "batch_stderr.log"

    proc = subprocess.Popen(
        cmd,
        cwd=str(PROJECT_ROOT),
        stdout=open(stdout_log, "w"),
        stderr=open(stderr_log, "w"),
    )
    return RunningWorker(job=job, process=proc, started=time.monotonic())


def save_log(jobs: list[Job], phase: str, log_path: Path):
    """Save/update the run log JSON."""
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Load existing log if present
    existing = {}
    if log_path.exists():
        try:
            existing = json.loads(log_path.read_text())
        except (json.JSONDecodeError, OSError):
            existing = {}

    completed = [j for j in jobs if j.status == "completed"]
    failed = [j for j in jobs if j.status == "failed"]
    pending = [j for j in jobs if j.status == "pending"]
    running = [j for j in jobs if j.status == "running"]

    entry = {
        "phase": phase,
        "updated_at": _now_iso(),
        "total": len(jobs),
        "completed": len(completed),
        "failed": len(failed),
        "pending": len(pending),
        "running": len(running),
        "jobs": [j.to_dict() for j in jobs],
    }

    existing[f"phase_{phase}"] = entry
    log_path.write_text(json.dumps(existing, indent=2), encoding="utf-8")


def run_batch(
    jobs: list[Job],
    workers: int,
    mesh_size: float,
    solver_timeout: int,
    resume: bool,
    phase: str,
):
    """Run jobs in parallel with up to N workers."""
    # Filter out completed jobs if resuming
    if resume:
        skipped = 0
        for job in jobs:
            if is_completed(job):
                job.status = "completed"
                job.returncode = 0
                skipped += 1
        print(f"[RESUME] Skipping {skipped}/{len(jobs)} already-completed jobs")

    pending = [j for j in jobs if j.status == "pending"]
    active: list[RunningWorker] = []
    total = len(jobs)
    completed_count = sum(1 for j in jobs if j.status == "completed")
    failed_count = 0

    print(f"\n{'='*60}")
    print(f"Phase 2 Batch: {phase} | Total: {total} | Pending: {len(pending)} | Workers: {workers}")
    print(f"Mesh size: {mesh_size} | Solver timeout: {solver_timeout}s")
    print(f"{'='*60}\n")

    idx = 0  # index into pending list

    while idx < len(pending) or active:
        # Launch new workers up to limit
        while len(active) < workers and idx < len(pending):
            job = pending[idx]
            idx += 1
            print(f"[START] {job.job_id} ({completed_count + failed_count + len(active)}/{total})")
            worker = launch_worker(job, mesh_size, solver_timeout)
            active.append(worker)

        # Poll active workers
        still_active = []
        for w in active:
            retcode = w.process.poll()
            if retcode is not None:
                elapsed = time.monotonic() - w.started
                w.job.end_time = _now_iso()
                w.job.returncode = retcode
                w.job.duration_s = round(elapsed, 1)
                if retcode == 0:
                    w.job.status = "completed"
                    completed_count += 1
                    print(f"[DONE]  {w.job.job_id} ({elapsed:.0f}s)")
                else:
                    w.job.status = "failed"
                    failed_count += 1
                    print(f"[FAIL]  {w.job.job_id} (rc={retcode}, {elapsed:.0f}s)")
                # Close file handles
                if w.process.stdout:
                    w.process.stdout.close()
                if w.process.stderr:
                    w.process.stderr.close()
            else:
                still_active.append(w)
        active = still_active

        # Save progress periodically
        save_log(jobs, phase, LOG_PATH)

        if active:
            time.sleep(5)

    # Final save
    save_log(jobs, phase, LOG_PATH)

    print(f"\n{'='*60}")
    print(f"Phase complete: {phase}")
    print(f"  Completed: {completed_count}/{total}")
    print(f"  Failed:    {failed_count}/{total}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 benchmark batch orchestration"
    )
    parser.add_argument(
        "--workers", "-w", type=int, default=3,
        help="Number of parallel workers (default: 3)",
    )
    parser.add_argument(
        "--mesh-size", type=float, default=0.18,
        help="Mesh size for all runs (default: 0.18)",
    )
    parser.add_argument(
        "--solver-timeout", type=int, default=3600,
        help="Solver timeout in seconds (default: 3600)",
    )
    parser.add_argument(
        "--phase", choices=["ref", "pred", "both"], default="both",
        help="Which phase to run: ref (20 reference), pred (100 predicted), both",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Skip already-completed cases (checks for VTK output)",
    )
    args = parser.parse_args()

    if args.phase in ("ref", "both"):
        ref_jobs = build_ref_jobs()
        if not ref_jobs:
            print("[ERROR] No reference scene files found in", SCENES_DIR)
            return 1
        run_batch(
            ref_jobs,
            workers=args.workers,
            mesh_size=args.mesh_size,
            solver_timeout=args.solver_timeout,
            resume=args.resume,
            phase="ref",
        )

    if args.phase in ("pred", "both"):
        pred_jobs = build_pred_jobs()
        if not pred_jobs:
            print("[ERROR] No predicted evaluation dirs found in", EVAL_DIR)
            return 1
        run_batch(
            pred_jobs,
            workers=args.workers,
            mesh_size=args.mesh_size,
            solver_timeout=args.solver_timeout,
            resume=args.resume,
            phase="pred",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
