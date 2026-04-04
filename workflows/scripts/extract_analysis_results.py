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


_NTS_ROOT_PATTERN = re.compile(r"nts\.mu2e\.mubeam\.Run1B\.(\d+)_(\d+)\.root$", re.IGNORECASE)
_COUNT_EVENTS_PATTERN = re.compile(r"^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$")
_EDEP_SUMMARY_PATTERN = re.compile(
    r"EdepAna summary:\s*Saw\s+(\d+)\s+events,\s*"
    r"average calo energy deposition per event:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"events with Edep > 50 MeV:\s*(\d+)"
)


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

    root_files = sorted(run_dir.glob("job_*/nts.mu2e.mubeam.Run1B.*_*.root"))
    print(f">>> Extracting TargetMuonFinder/stopmat entries from {len(root_files)} ROOT files")

    per_file = []
    total_target_al_entries = 0.0

    for root_path in root_files:
        match = _NTS_ROOT_PATTERN.search(root_path.name)
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
    print(f">>> Counting events in {len(art_files)} art files")

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


def _run_edep_analysis(run_dir: Path) -> dict:
    flash_art_files = sorted(run_dir.glob("job_*/dts.mu2e.MuBeamFlash.Run1B.*_*.art"))
    list_path = run_dir / "mubeam_flash_art_files.txt"
    log_path = run_dir / "edep_analysis.log"

    if not flash_art_files:
        return {
            "ran": False,
            "error": "No MuBeamFlash art files found",
            "input_file_count": 0,
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": None,
            "wrapper_fcl_path": None,
            "nts_output_path": None,
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
        }

    print(f">>> Found {len(flash_art_files)} MuBeamFlash art files for EdepAna processing")
    list_path.write_text("\n".join(str(path) for path in flash_art_files) + "\n", encoding="utf-8")

    workflows_dir = Path(__file__).resolve().parent.parent
    edep_fcl = workflows_dir / "fcl" / "edep.fcl"
    if not edep_fcl.exists():
        return {
            "ran": False,
            "error": f"Missing edep fcl: {edep_fcl}",
            "input_file_count": len(flash_art_files),
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": str(edep_fcl),
            "wrapper_fcl_path": None,
            "nts_output_path": None,
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
        }

    command = ["mu2e", "-c", str(edep_fcl), "-S", str(list_path)]
    try:
        proc = subprocess.run(
            command,
            cwd=run_dir,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        log_path.write_text("mu2e executable not found on PATH\n", encoding="utf-8")
        return {
            "ran": False,
            "error": "mu2e executable not found on PATH",
            "input_file_count": len(flash_art_files),
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": str(edep_fcl),
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
        }

    log_path.write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")

    summary_line = None
    events_seen = None
    average_calo_energy_mev = None
    events_edep_gt_50_mev = None
    for line in (proc.stdout + "\n" + proc.stderr).splitlines():
        match = _EDEP_SUMMARY_PATTERN.search(line)
        if match:
            summary_line = line.strip()
            events_seen = int(match.group(1))
            average_calo_energy_mev = float(match.group(2))
            events_edep_gt_50_mev = int(match.group(3))

    return {
        "ran": proc.returncode == 0,
        "error": None if proc.returncode == 0 else f"mu2e exited with code {proc.returncode}",
        "input_file_count": len(flash_art_files),
        "input_list_path": str(list_path),
        "log_path": str(log_path),
        "fcl_path": str(edep_fcl),
        "mu2e_returncode": proc.returncode,
        "summary_line": summary_line,
        "events_seen": events_seen,
        "average_calo_energy_mev": average_calo_energy_mev,
        "events_edep_gt_50_mev": events_edep_gt_50_mev,
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
    edep_stats = _run_edep_analysis(run_dir)

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
        "edep_analysis": edep_stats,
        "jobs": status_rows,
        "warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        default=None,
        help="Directory containing job_* folders (required unless --from-json is used)",
    )
    parser.add_argument(
        "--from-json",
        default=None,
        help="Load an existing analysis_summary.json and print results without recomputing",
    )
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


def _print_pretty_summary(summary: dict) -> None:
    print("-----")
    print(f"Jobs: {summary['completed_jobs']}/{summary['total_jobs']} completed")
    print(f"Failed: {summary['failed_jobs']}, Unknown: {summary['unknown_jobs']}")
    print("Output counts:")
    nFlash = summary["output_counts"].get("FlashOutput", 0)
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

    edep = summary["edep_analysis"]
    print("\nEdep reprocessing summary:")
    print(f"  Input MuBeamFlash files: {edep['input_file_count']}")
    print(f"  Input list file: {edep['input_list_path']}")
    print(f"  Edep log: {edep['log_path']}")
    if edep["summary_line"]:
        print(f"  {edep['summary_line']}")
        edep_avg = edep["average_calo_energy_mev"]
        ndep_gt_50 = edep["events_edep_gt_50_mev"]
        if edep_avg is not None:
            nsim_events = sim_events["total_events"] if sim_events["total_events"] else 0
            edep_avg_wt = edep_avg * nFlash / nsim_events if nsim_events > 0 else 0
            ndep_gt_50_wt = ndep_gt_50 * nFlash / nsim_events if nsim_events > 0 else 0
            print(f"  Average calo energy deposition per sim event: {edep_avg_wt:.3g} MeV")
            print(f"  Events with Edep > 50 MeV per sim event: {ndep_gt_50_wt:.3g}")
    else:
        print("  EdepAna summary line not found in mu2e output")
    if edep["error"]:
        print(f"  Note: {edep['error']}")


def main() -> int:
    args = parse_args()
    if args.from_json:
        json_path = Path(args.from_json).resolve()
        if not json_path.exists() or not json_path.is_file():
            raise SystemExit(f"Summary JSON does not exist: {json_path}")
        with json_path.open("r", encoding="utf-8") as handle:
            summary = json.load(handle)
        print(f"Loaded summary: {json_path}")
    else:
        if not args.run_dir:
            raise SystemExit("--run-dir is required unless --from-json is used")

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
        _print_pretty_summary(summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
