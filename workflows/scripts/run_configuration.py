#!/usr/bin/env python3
"""Run parallel mu2e jobs for a selected configuration and summarize outputs."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


@dataclass
class JobResult:
    index: int
    job_dir: Path
    command: list[str]
    returncode: int
    duration_s: float
    log_path: Path


def _run_one_job(index: int, command: list[str], job_dir: Path, env: dict[str, str], dry_run: bool) -> JobResult:
    log_path = job_dir / "job.log"
    start = time.perf_counter()

    if dry_run:
        log_path.write_text("DRY RUN: " + shlex.join(command) + "\n", encoding="utf-8")
        return JobResult(
            index=index,
            job_dir=job_dir,
            command=command,
            returncode=0,
            duration_s=0.0,
            log_path=log_path,
        )

    with log_path.open("w", encoding="utf-8") as log_file:
        process = subprocess.run(
            command,
            cwd=job_dir,
            env=env,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    duration = time.perf_counter() - start
    return JobResult(
        index=index,
        job_dir=job_dir,
        command=command,
        returncode=process.returncode,
        duration_s=duration,
        log_path=log_path,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_version", help="Configuration folder name (for example: config_v06)")
    parser.add_argument("parallel_jobs", type=int, help="Number of jobs to launch")
    parser.add_argument(
        "--events-per-job",
        type=int,
        default=None,
        help="Optional value passed as '-n <events>' to mu2e",
    )
    parser.add_argument(
        "--mu2e-command",
        default="mu2e",
        help="Executable used to run jobs (default: mu2e)",
    )
    parser.add_argument(
        "--seed-start",
        type=int,
        default=1,
        help="First seed value; each job gets seed_start + job_index",
    )
    parser.add_argument(
        "--run-root",
        default=None,
        help="Directory where run outputs are written (default: ../runs from this script)",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Thread workers used to launch jobs (default: parallel_jobs)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print and stage commands without executing mu2e",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.parallel_jobs <= 0:
        raise SystemExit("parallel_jobs must be > 0")
    if args.seed_start <= 0:
        raise SystemExit("seed_start must be > 0")

    script_dir = Path(__file__).resolve().parent
    workflows_dir = script_dir.parent
    fcl_path = workflows_dir / args.config_version / "run1b_beam" / "mubeam.fcl"
    extractor_path = script_dir / "extract_analysis_results.py"

    if not fcl_path.exists():
        raise SystemExit(f"Missing FCL file: {fcl_path}")
    if not extractor_path.exists():
        raise SystemExit(f"Missing extractor script: {extractor_path}")

    run_root = Path(args.run_root).resolve() if args.run_root else workflows_dir / "runs"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = run_root / args.config_version / f"mubeam_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)

    print(f"FCL: {fcl_path}")
    print(f"Run directory: {run_dir}")
    print(f"Launching {args.parallel_jobs} jobs")
    print(f"Seed range: {args.seed_start} to {args.seed_start + args.parallel_jobs - 1}")

    max_workers = args.max_workers if args.max_workers else args.parallel_jobs
    max_workers = max(1, min(max_workers, args.parallel_jobs))

    env = os.environ.copy()
    results: list[JobResult] = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for index in range(args.parallel_jobs):
            job_dir = run_dir / f"job_{index:03d}"
            job_dir.mkdir(parents=True, exist_ok=False)

            seed = args.seed_start + index
            job_fcl = job_dir / "mubeam_job.fcl"
            job_fcl.write_text(
                f"#include \"{fcl_path}\"\n"
                "\n"
                f"services.SeedService.baseSeed : {seed}\n",
                encoding="utf-8",
            )

            command = [args.mu2e_command, "-c", str(job_fcl)]
            if args.events_per_job is not None:
                command.extend(["-n", str(args.events_per_job)])

            # Record command for reproducibility.
            (job_dir / "job_command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
            futures.append(executor.submit(_run_one_job, index, command, job_dir, env, args.dry_run))

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            status = {
                "job_index": result.index,
                "job_dir": str(result.job_dir),
                "command": result.command,
                "returncode": result.returncode,
                "duration_s": round(result.duration_s, 3),
                "log_path": str(result.log_path),
            }
            with (result.job_dir / "job_status.json").open("w", encoding="utf-8") as handle:
                json.dump(status, handle, indent=2, sort_keys=True)

            print(
                f"Job {result.index:03d}: returncode={result.returncode}, "
                f"duration={result.duration_s:.2f}s, log={result.log_path}"
            )

    completed = sum(1 for result in results if result.returncode == 0)
    failed = len(results) - completed
    print(f"Finished jobs: {completed}/{len(results)} successful, {failed} failed")

    summary_path = run_dir / "analysis_summary.json"
    extractor_cmd = [
        sys.executable,
        str(extractor_path),
        "--run-dir",
        str(run_dir),
        "--output",
        str(summary_path),
        "--pretty",
    ]

    print("Running analysis extractor...")
    extractor_run = subprocess.run(extractor_cmd, check=False)
    if extractor_run.returncode != 0:
        print("Extractor failed", file=sys.stderr)
        return extractor_run.returncode

    print(f"Analysis summary: {summary_path}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
