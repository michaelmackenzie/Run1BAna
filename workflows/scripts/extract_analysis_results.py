#!/usr/bin/env python3
"""Summarize outputs from parallel mu2e mubeam jobs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


_OUTPUT_PATTERNS = {
    "EarlyFlashOutput": re.compile(r"EarlyMuBeamFlash", re.IGNORECASE),
    "FlashOutput": re.compile(r"MuBeamFlash", re.IGNORECASE),
    "IPAStopOutput": re.compile(r"IPAStops", re.IGNORECASE),
    "TargetStopOutput": re.compile(r"TargetStops", re.IGNORECASE),
    "TFileService": re.compile(r"mubeam.*\.root$", re.IGNORECASE),
}


_ERROR_PATTERNS = [
    re.compile(r"\bERROR\b"),
    re.compile(r"ArtException"),
    re.compile(r"Fatal Exception"),
]


def _classify_output(path: Path) -> str:
    name = path.name
    for label, pattern in _OUTPUT_PATTERNS.items():
        if pattern.search(name):
            return label
    if name.endswith(".art"):
        return "other_art"
    if name.endswith(".root"):
        return "other_root"
    return "other"


def _scan_logs(log_path: Path) -> dict:
    errors = []
    line_count = 0
    if not log_path.exists():
        return {"line_count": 0, "possible_errors": errors}

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            line_count += 1
            if any(pattern.search(line) for pattern in _ERROR_PATTERNS):
                errors.append({"line": line_number, "text": line.rstrip()[:300]})

    return {"line_count": line_count, "possible_errors": errors}


def build_summary(run_dir: Path) -> dict:
    job_dirs = sorted(path for path in run_dir.glob("job_*") if path.is_dir())

    output_counts: dict[str, int] = {}
    output_bytes: dict[str, int] = {}
    status_rows = []
    warnings = []

    for job_dir in job_dirs:
        status_file = job_dir / "job_status.json"
        if status_file.exists():
            with status_file.open("r", encoding="utf-8") as handle:
                status = json.load(handle)
        else:
            status = {"job_dir": str(job_dir), "returncode": None, "duration_s": None}
            warnings.append(f"Missing status file: {status_file}")

        log_path = job_dir / "job.log"
        log_scan = _scan_logs(log_path)
        status["log"] = {
            "path": str(log_path),
            "line_count": log_scan["line_count"],
            "possible_error_count": len(log_scan["possible_errors"]),
            "possible_errors": log_scan["possible_errors"],
        }

        outputs = []
        for ext in ("*.art", "*.root"):
            for out_file in sorted(job_dir.glob(ext)):
                file_type = _classify_output(out_file)
                size_bytes = out_file.stat().st_size
                output_counts[file_type] = output_counts.get(file_type, 0) + 1
                output_bytes[file_type] = output_bytes.get(file_type, 0) + size_bytes
                outputs.append(
                    {
                        "path": str(out_file),
                        "type": file_type,
                        "size_bytes": size_bytes,
                    }
                )

        status["outputs"] = outputs
        status_rows.append(status)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    return {
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "jobs": status_rows,
        "warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="Directory containing job_* folders")
    parser.add_argument(
        "--output",
        default=None,
        help="Path to JSON summary (default: <run_dir>/analysis_summary.json)",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Also print a compact human-readable summary",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()

    if not run_dir.exists() or not run_dir.is_dir():
        raise SystemExit(f"Run directory does not exist: {run_dir}")

    summary = build_summary(run_dir)
    output_path = Path(args.output).resolve() if args.output else run_dir / "analysis_summary.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(f"Wrote summary: {output_path}")

    if args.pretty:
        print("-----")
        print(f"Jobs: {summary['completed_jobs']}/{summary['total_jobs']} completed")
        print(f"Failed: {summary['failed_jobs']}, Unknown: {summary['unknown_jobs']}")
        print("Output counts:")
        for key, value in sorted(summary["output_counts"].items()):
            print(f"  {key}: {value}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
