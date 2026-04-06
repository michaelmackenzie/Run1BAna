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


_STAGES = ("mubeam", "mustop")
_MUSTOP_MODES = ("ce", "ce_plus", "flat_gamma")


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
    parser.add_argument(
        "--stage",
        choices=_STAGES,
        default="mubeam",
        help="Processing stage to run (default: mubeam)",
    )
    parser.add_argument("config_version", help="Configuration folder name (for example: config_v06)")
    parser.add_argument(
        "parallel_jobs",
        type=int,
        nargs="?",
        default=None,
        help="Number of jobs to launch (required for stage mubeam)",
    )
    parser.add_argument(
        "--events-per-job",
        type=int,
        required=True,
        help="Number of events per job passed as '-n <events>' to mu2e",
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
    parser.add_argument(
        "--mubeam-run-dir",
        default=None,
        help=(
            "Directory containing mubeam job_* outputs used as TargetStopResampler inputs "
            "for --stage mustop (default: latest mubeam_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--mustop-jobs-per-mode",
        type=int,
        default=1,
        help="Number of jobs to launch for each mustop mode (ce/ce_plus/flat_gamma)",
    )
    return parser.parse_args()


def _find_latest_stage_run(run_root: Path, config_version: str, stage: str) -> Path | None:
    stage_parent = run_root / config_version
    if not stage_parent.exists() or not stage_parent.is_dir():
        return None

    candidates = sorted(
        path for path in stage_parent.glob(f"{stage}_*") if path.is_dir()
    )
    return candidates[-1] if candidates else None


def _collect_target_stop_files(run_dir: Path) -> list[Path]:
    return sorted(run_dir.glob("job_*/sim.mu2e.TargetStops.Run1B.*_*.art"))


def _collect_muminus_stop_files(mubeam_run_dir: Path) -> list[Path]:
    return sorted((mubeam_run_dir / "mu_stops_job").glob("sim.mu2e.MuminusStopsCat.Run1B.*_*.art"))


def _format_fhicl_string_list(paths: list[Path], indent: str = "    ") -> str:
    return ",\n".join(f'{indent}"{path}"' for path in paths)


def _run_mu_stops_job(
    run_dir: Path,
    workflows_dir: Path,
    config_version: str,
    mu2e_command: str,
    env: dict,
    dry_run: bool,
) -> JobResult:
    target_stop_files = _collect_target_stop_files(run_dir)
    if not target_stop_files:
        raise SystemExit(f"No TargetStops files found in {run_dir} for mu_stops job")

    job_dir = run_dir / "mu_stops_job"
    job_dir.mkdir(parents=True, exist_ok=False)

    include_fcl_path = Path("Run1BAna") / "workflows" / config_version / "run1b_beam" / "mu_stops.fcl"
    file_list = _format_fhicl_string_list(target_stop_files)
    job_fcl = job_dir / "mu_stops_job.fcl"
    job_fcl.write_text(
        f"#include \"{include_fcl_path.as_posix()}\"\n"
        "\n"
        f"source.fileNames: [\n"
        f"{file_list}\n"
        "]\n",
        encoding="utf-8",
    )

    command = [mu2e_command, "-c", str(job_fcl)]
    (job_dir / "job_command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")

    print(f"Running mu_stops job with {len(target_stop_files)} TargetStops input files...")
    result = _run_one_job(0, command, job_dir, env, dry_run)

    status = {
        "job_index": 0,
        "job_dir": str(result.job_dir),
        "command": result.command,
        "returncode": result.returncode,
        "duration_s": round(result.duration_s, 3),
        "log_path": str(result.log_path),
    }
    with (job_dir / "job_status.json").open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)

    print(
        f"mu_stops job: returncode={result.returncode}, "
        f"duration={result.duration_s:.2f}s, log={result.log_path}"
    )
    return result


def main() -> int:
    args = parse_args()

    if args.stage == "mubeam" and (args.parallel_jobs is None or args.parallel_jobs <= 0):
        raise SystemExit("parallel_jobs must be > 0 for stage mubeam")
    if args.events_per_job <= 0:
        raise SystemExit("events_per_job must be > 0")
    if args.seed_start <= 0:
        raise SystemExit("seed_start must be > 0")
    if args.mustop_jobs_per_mode <= 0:
        raise SystemExit("mustop_jobs_per_mode must be > 0")

    script_dir = Path(__file__).resolve().parent
    workflows_dir = script_dir.parent
    extractor_path = script_dir / "extract_analysis_results.py"

    run_root = Path(args.run_root).resolve() if args.run_root else workflows_dir / "runs"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = run_root / args.config_version / f"{args.stage}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)

    if args.stage == "mubeam":
        fcl_path = workflows_dir / args.config_version / "run1b_beam" / "mubeam.fcl"
        include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / "run1b_beam" / "mubeam.fcl"
        if not fcl_path.exists():
            raise SystemExit(f"Missing FCL file: {fcl_path}")
        if not extractor_path.exists():
            raise SystemExit(f"Missing extractor script: {extractor_path}")

        job_specs = [
            {
                "index": index,
                "name": "mubeam",
                "job_fcl_name": "mubeam_job.fcl",
                "include_fcl_path": include_fcl_path,
                "fcl_overrides": "",
            }
            for index in range(args.parallel_jobs)
        ]
    else:
        mubeam_run_dir = (
            Path(args.mubeam_run_dir).resolve()
            if args.mubeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mubeam")
        )

        if mubeam_run_dir is None:
            raise SystemExit(
                "Could not find a mubeam run directory. Provide --mubeam-run-dir or run stage mubeam first."
            )
        if not mubeam_run_dir.exists() or not mubeam_run_dir.is_dir():
            raise SystemExit(f"mubeam run directory does not exist: {mubeam_run_dir}")

        muminus_stop_files = _collect_muminus_stop_files(mubeam_run_dir)
        if not muminus_stop_files:
            raise SystemExit(
                f"No MuminusStopsCat files found in {mubeam_run_dir / 'mu_stops_job'}. "
                "Run stage mubeam to completion first."
            )

        print(f"Using mubeam inputs from: {mubeam_run_dir}")
        print(f"MuminusStopsCat input files: {len(muminus_stop_files)}")

        job_specs = []
        for mode in _MUSTOP_MODES:
            mode_fcl_path = workflows_dir / args.config_version / "run1b_mustop" / f"{mode}.fcl"
            include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / "run1b_mustop" / f"{mode}.fcl"
            if not mode_fcl_path.exists():
                raise SystemExit(f"Missing FCL file: {mode_fcl_path}")

            muminus_stop_lines = _format_fhicl_string_list(muminus_stop_files)
            overrides = (
                "physics.filters.TargetStopResampler.fileNames: [\n"
                f"{muminus_stop_lines}\n"
                "]\n"
            )
            for mode_job_index in range(args.mustop_jobs_per_mode):
                index = len(job_specs)
                job_specs.append(
                    {
                        "index": index,
                        "name": f"{mode}_{mode_job_index:03d}",
                        "job_fcl_name": f"{mode}_job_{mode_job_index:03d}.fcl",
                        "include_fcl_path": include_fcl_path,
                        "fcl_overrides": overrides,
                    }
                )

    if args.stage == "mubeam":
        print(f"FCL: {workflows_dir / args.config_version / 'run1b_beam' / 'mubeam.fcl'}")
    else:
        print(
            "FCLs: "
            + ", ".join(
                str(workflows_dir / args.config_version / "run1b_mustop" / f"{mode}.fcl")
                for mode in _MUSTOP_MODES
            )
        )
    print(f"Run directory: {run_dir}")
    print(f"Stage: {args.stage}")
    if args.stage == "mustop":
        print(f"mustop jobs per mode: {args.mustop_jobs_per_mode}")
    print(f"Launching {len(job_specs)} jobs")
    print(f"Seed range: {args.seed_start} to {args.seed_start + len(job_specs) - 1}")

    max_workers = args.max_workers if args.max_workers else len(job_specs)
    max_workers = max(1, min(max_workers, len(job_specs)))

    env = os.environ.copy()

    print("Running getToken to ensure access credentials are ready...")
    token_result = subprocess.run(["getToken"], env=env, check=False)
    if token_result.returncode != 0:
        raise SystemExit(f"getToken failed with exit code {token_result.returncode}")

    results: list[JobResult] = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for spec in job_specs:
            index = spec["index"]
            job_dir = run_dir / f"job_{index:03d}"
            job_dir.mkdir(parents=True, exist_ok=False)

            seed = args.seed_start + index
            job_fcl = job_dir / spec["job_fcl_name"]
            job_fcl.write_text(
                f"#include \"{spec['include_fcl_path'].as_posix()}\"\n"
                "\n"
                f"services.SeedService.baseSeed : {seed}\n"
                f"source.firstSubRun: {index}\n"
                f"{spec['fcl_overrides']}",
                encoding="utf-8",
            )

            command = [args.mu2e_command, "-c", str(job_fcl)]
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

    if args.stage != "mubeam":
        print(f"Stage {args.stage} complete: {run_dir}")
        return 0 if failed == 0 else 1

    if failed > 0:
        print("Skipping mu_stops job due to failed mubeam jobs", file=sys.stderr)
        return 1

    mu_stops_result = _run_mu_stops_job(
        run_dir, workflows_dir, args.config_version,
        args.mu2e_command, env, args.dry_run,
    )
    if mu_stops_result.returncode != 0:
        print("mu_stops job failed", file=sys.stderr)
        return 1

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
