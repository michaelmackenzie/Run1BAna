#!/usr/bin/env python3
"""Summarize outputs from parallel mu2e mubeam jobs."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
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


_TS_ROOT_PATTERN = re.compile(r"ts\.mu2e\.mubeam\.Run1B\.001460_(\d+)\.root$", re.IGNORECASE)
_NTS_ROOT_PATTERN = re.compile(r"nts\.mu2e\.mubeam\.Run1B\.001460_(\d+)\.root$", re.IGNORECASE)
_COUNT_EVENTS_PATTERN = re.compile(r"^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$")


def _parse_events_from_command(command: object) -> int | None:
    if not isinstance(command, list):
        return None

    for i, token in enumerate(command):
        if token == "-n" and i + 1 < len(command):
            try:
                return int(command[i + 1])
            except ValueError:
                return None
    return None


def _compute_total_simulated_events(status_rows: list[dict]) -> dict:
    total_events = 0
    jobs_with_known_events = 0

    for row in status_rows:
        events = _parse_events_from_command(row.get("command"))
        if events is None:
            continue
        jobs_with_known_events += 1
        total_events += events

    return {
        "jobs_with_known_events": jobs_with_known_events,
        "total_jobs": len(status_rows),
        "total_events": total_events if jobs_with_known_events == len(status_rows) else None,
    }


def _extract_target_al_entries(run_dir: Path) -> dict:
    try:
        import ROOT  # type: ignore
    except Exception as exc:  # noqa: BLE001
        return {
            "root_available": False,
            "error": f"PyROOT import failed: {exc}",
            "files_analyzed": 0,
            "total_target_al_entries": 0.0,
            "per_file": [],
        }

    root_files = sorted(run_dir.glob("job_*/ts.mu2e.mubeam.Run1B.001460_*.root"))
    if not root_files:
        root_files = sorted(run_dir.glob("job_*/nts.mu2e.mubeam.Run1B.001460_*.root"))

    per_file = []
    total_target_al_entries = 0.0

    for root_path in root_files:
        match = _TS_ROOT_PATTERN.search(root_path.name) or _NTS_ROOT_PATTERN.search(root_path.name)
        subrun = int(match.group(1)) if match else None
        row = {
            "path": str(root_path),
            "subrun": subrun,
            "hist_found": False,
            "bin_found": False,
            "target_al_entries": 0.0,
            "error": None,
        }

        tfile = ROOT.TFile.Open(str(root_path), "READ")
        if not tfile or tfile.IsZombie():
            row["error"] = "Failed to open ROOT file"
            per_file.append(row)
            continue

        hist = tfile.Get("TargetMuonFinder/stopmat")
        if not hist:
            row["error"] = "Histogram TargetMuonFinder/stopmat not found"
            tfile.Close()
            per_file.append(row)
            continue

        row["hist_found"] = True
        xaxis = hist.GetXaxis()
        for bin_idx in range(1, xaxis.GetNbins() + 1):
            if xaxis.GetBinLabel(bin_idx) == "StoppingTarget_Al":
                entries = float(hist.GetBinContent(bin_idx))
                row["bin_found"] = True
                row["target_al_entries"] = entries
                total_target_al_entries += entries
                break

        if not row["bin_found"]:
            row["error"] = "Bin label StoppingTarget_Al not found"

        tfile.Close()
        per_file.append(row)

    return {
        "root_available": True,
        "error": None,
        "files_analyzed": len(root_files),
        "total_target_al_entries": total_target_al_entries,
        "per_file": per_file,
    }


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


def _count_art_file_events(art_path: Path) -> int | None:
    try:
        result = subprocess.run(
            ["count_events", str(art_path)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        for line in result.stdout.split("\n"):
            match = _COUNT_EVENTS_PATTERN.match(line.strip())
            if match:
                return int(match.group(4))
        return None
    except FileNotFoundError:
        return None
    except (subprocess.TimeoutExpired, Exception):
        return None


def _extract_art_event_counts(run_dir: Path) -> dict:
    art_files = sorted(run_dir.glob("job_*/*.art"))
    per_file: list[dict] = []
    total_events_by_type: dict[str, int] = {}
    total_events = 0

    for art_path in art_files:
        file_type = _classify_output(art_path)
        event_count = _count_art_file_events(art_path)

        per_file.append({
            "path": str(art_path),
            "type": file_type,
            "event_count": event_count,
        })

        if event_count is not None:
            total_events_by_type[file_type] = total_events_by_type.get(file_type, 0) + event_count
            total_events += event_count

    return {
        "files_analyzed": len(art_files),
        "total_events": total_events,
        "events_by_type": total_events_by_type,
        "per_file": per_file,
    }


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

    event_stats = _compute_total_simulated_events(status_rows)
    target_al_stats = _extract_target_al_entries(run_dir)
    art_event_stats = _extract_art_event_counts(run_dir)

    if event_stats["total_events"] and event_stats["total_events"] > 0:
        target_al_per_event = target_al_stats["total_target_al_entries"] / event_stats["total_events"]
    else:
        target_al_per_event = None

    art_events_per_simulated = {}
    if event_stats["total_events"] and event_stats["total_events"] > 0:
        for file_type, count in art_event_stats["events_by_type"].items():
            art_events_per_simulated[file_type] = count / event_stats["total_events"]

    return {
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "target_al_analysis": {
            **target_al_stats,
            "target_al_entries_per_simulated_event": target_al_per_event,
        },
        "art_event_analysis": {
            **art_event_stats,
            "events_per_simulated_event": art_events_per_simulated,
        },
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

        target_al = summary["target_al_analysis"]
        sim_events = summary["simulation_events"]
        art_events = summary["art_event_analysis"]

        print("TargetMuonFinder/stopmat summary:")
        print(f"  ROOT files analyzed: {target_al['files_analyzed']}")
        print(f"  Total StoppingTarget_Al entries: {target_al['total_target_al_entries']}")

        if sim_events["total_events"] is None:
            print("  Total simulated events: unknown (not all jobs had '-n' in command)")
            print("  StoppingTarget_Al per simulated event: unavailable")
        else:
            print(f"  Total simulated events: {sim_events['total_events']}")
            print(
                "  StoppingTarget_Al per simulated event: "
                f"{target_al['target_al_entries_per_simulated_event']:.8g}"
            )

        if target_al["error"]:
            print(f"  Note: {target_al['error']}")

        print("\nArt file event counts:")
        print(f"  Files analyzed: {art_events['files_analyzed']}")
        print(f"  Total events: {art_events['total_events']}")
        if art_events["events_by_type"]:
            print("  Events by type:")
            for file_type, count in sorted(art_events["events_by_type"].items()):
                per_sim = art_events["events_per_simulated_event"].get(file_type, 0)
                print(f"    {file_type}: {count} ({per_sim:.8g} per simulated event)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
