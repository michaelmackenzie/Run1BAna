#!/usr/bin/env python3
"""Summarize outputs from staged parallel mu2e jobs."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


_OUTPUT_PATTERNS = {
    "EarlyFlashOutput": re.compile(r"EarlyMuBeamFlash", re.IGNORECASE),
    "FlashOutput": re.compile(r"MuBeamFlash", re.IGNORECASE),
    "EleFlashOutput": re.compile(r"EleBeamFlash", re.IGNORECASE),
    "IPAStopOutput": re.compile(r"IPAStops", re.IGNORECASE),
    "TargetStopOutput": re.compile(r"TargetStops", re.IGNORECASE),
    "MuminusStopOutput": re.compile(r"MuminusStopsCat", re.IGNORECASE),
    "MuplusStopOutput": re.compile(r"MuplusStopsCat", re.IGNORECASE),
    "CeEndpointOutput": re.compile(r"CeEndpoint", re.IGNORECASE),
    "CePlusEndpointOutput": re.compile(r"CePlusEndpoint", re.IGNORECASE),
    "FlatGammaOutput": re.compile(r"FlatGamma", re.IGNORECASE),
    "PileupOutput": re.compile(r"Pileup", re.IGNORECASE),
    "TFileService": re.compile(r"mubeam.*\.root$", re.IGNORECASE),
}


_ERROR_PATTERNS = [
    re.compile(r"\bERROR\b"),
    re.compile(r"ArtException"),
    re.compile(r"Fatal Exception"),
]


_COUNT_EVENTS_PATTERN = re.compile(r"^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$")
_EDEP_SUMMARY_PATTERN = re.compile(
    r"EdepAna summary:\s*Saw\s+(\d+)\s+events,\s*"
    r"average calo energy deposition per event:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"events with Edep > 50 MeV:\s*(\d+)"
)
_EDEP_FIT_PATTERN = re.compile(
    r"Primary energy - Edep fit:\s*status\s*=\s*(\d+)\s*,\s*mean\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"FWHM\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV"
)
_EDEP_DISTRIBUTION_PATTERN = re.compile(
    r"Primary energy - Edep distribution:\s*mpv\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"FWHM\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV"
)
_TRACKER_FRONT_FIT_PATTERN = re.compile(
    r"Tracker front - primary energy fit:\s*status\s*=\s*(\d+)\s*,\s*mean\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"FWHM\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV"
)
_PRIMARY_EDEP_MINUS_TRACKER_FRONT_DISTRIBUTION_PATTERN = re.compile(
    r"Primary Edep - tracker front StepPointMC energy distribution:\s*MPV\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"FWHM\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*MeV,\s*"
    r"efficiency\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)
_CALO_STOP_MATERIALS = ("G4_CESIUM_IODIDE", "CarbonFiber", "AluminumHoneycomb")
_STAGES = ("mubeam", "elebeam", "mustop", "mustop_pileup", "run1a_mubeam", "run1a_mustops")
_MUSTOP_EDEP_SAMPLES = {
    "ce": {
        "glob": "job_*/dts.mu2e.CeEndpoint.Run1B.*_*.art",
        "label": "CeEndpoint",
    },
    "ce_plus": {
        "glob": "job_*/dts.mu2e.CePlusEndpoint.Run1B.*_*.art",
        "label": "CePlusEndpoint",
    },
    "flat_gamma": {
        "glob": "job_*/dts.mu2e.FlatGamma.Run1B.*_*.art",
        "label": "FlatGamma",
    },
}
_MUBEAM_INPUT_EFFICIENCY_BY_FCL = {
    "run1b_beam/mubeam.fcl": 0.01278168,
    "run1a_beam/mubeam.fcl": 0.01278168, # 0.00213816
    "run1b_beam/elebeam.fcl": 0.086211877,
    "run1a_beam/elebeam.fcl": 0.01730766,
}


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
    successful_rows = [row for row in status_rows if row.get("returncode") == 0]
    total_events = 0
    jobs_with_known_events = 0

    for row in successful_rows:
        events = _parse_events_from_command(row.get("command"))
        if events is None:
            continue
        jobs_with_known_events += 1
        total_events += events

    return {
        "jobs_with_known_events": jobs_with_known_events,
        "total_jobs": len(successful_rows),
        "total_jobs_all": len(status_rows),
        "total_events": total_events if jobs_with_known_events == len(successful_rows) else None,
    }


def _detect_mubeam_input_efficiency(run_dir: Path) -> dict:
    for job_fcl in sorted(run_dir.glob("job_*/*.fcl")):
        try:
            text = job_fcl.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

        for fcl_fragment, correction in _MUBEAM_INPUT_EFFICIENCY_BY_FCL.items():
            if fcl_fragment in text:
                return {
                    "detected": True,
                    "detected_from": str(job_fcl),
                    "input_fcl": fcl_fragment,
                    "correction_factor": correction,
                }

    return {
        "detected": False,
        "detected_from": None,
        "input_fcl": None,
        "correction_factor": None,
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
            "total_calo_entries": 0.0,
            "per_file": [],
        }

    root_files = sorted(run_dir.glob("job_*/nts.mu2e.mubeam.Run1A.*_*.root"))
    if not root_files:
        root_files = sorted(run_dir.glob("job_*/nts.mu2e.mubeam.Run1B.*_*.root"))
    nts_root_pattern = re.compile(r"nts\.mu2e\.mubeam\.Run1[AB]\.(\d+)_(\d+)\.root$", re.IGNORECASE)
    print(f">>> Extracting TargetMuonFinder/stopmat entries from {len(root_files)} ROOT files")

    per_file = []
    total_target_al_entries = 0.0
    total_calo_entries = 0.0

    for root_path in root_files:
        match = nts_root_pattern.search(root_path.name)
        subrun = int(match.group(1)) if match else None
        row = {
            "path": str(root_path),
            "subrun": subrun,
            "hist_found": False,
            "bin_found": False,
            "target_al_entries": 0.0,
            "calo_entries": 0.0,
            "calo_material_entries": {},
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
            label = xaxis.GetBinLabel(bin_idx)
            entries = float(hist.GetBinContent(bin_idx))

            if label == "StoppingTarget_Al":
                row["bin_found"] = True
                row["target_al_entries"] = entries
                total_target_al_entries += entries
            elif label in _CALO_STOP_MATERIALS:
                row["calo_material_entries"][label] = entries
                row["calo_entries"] += entries

        total_calo_entries += row["calo_entries"]

        errors = []
        missing_calo_labels = [
            material for material in _CALO_STOP_MATERIALS if material not in row["calo_material_entries"]
        ]

        if not row["bin_found"]:
            errors.append("Bin label StoppingTarget_Al not found")
        if missing_calo_labels:
            errors.append("Calo material bin labels not found: " + ", ".join(missing_calo_labels))

        row["error"] = "; ".join(errors) if errors else None

        tfile.Close()
        per_file.append(row)

    return {
        "root_available": True,
        "error": None,
        "files_analyzed": len(root_files),
        "total_target_al_entries": total_target_al_entries,
        "total_calo_entries": total_calo_entries,
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


def _run_edep_analysis_for_files(
    run_dir: Path,
    art_files: list[Path],
    *,
    analysis_name: str,
    input_label: str,
    list_filename: str,
    log_filename: str,
    wrapper_fcl_filename: str,
    nts_output_filename: str,
) -> dict:
    list_path = run_dir / list_filename
    log_path = run_dir / log_filename
    wrapper_fcl_path = run_dir / wrapper_fcl_filename
    nts_output_path = run_dir / nts_output_filename

    if not art_files:
        return {
            "ran": False,
            "error": f"No {input_label} art files found",
            "input_file_count": 0,
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": None,
            "wrapper_fcl_path": str(wrapper_fcl_path),
            "nts_output_path": str(nts_output_path),
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
            "primary_minus_edep_fit_line": None,
            "primary_minus_edep_fit_mean_mev": None,
            "primary_minus_edep_fit_fwhm_mev": None,
            "primary_minus_edep_distribution_line": None,
            "primary_minus_edep_distribution_mean_mev": None,
            "primary_minus_edep_distribution_rms_mev": None,
            "tracker_front_fit_line": None,
            "tracker_front_fit_status": None,
            "tracker_front_fit_mpv_mev": None,
            "tracker_front_fit_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_line": None,
            "primary_edep_minus_tracker_front_distribution_mpv_mev": None,
            "primary_edep_minus_tracker_front_distribution_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_efficiency": None,
        }

    print(f">>> Found {len(art_files)} {input_label} art files for {analysis_name} EdepAna processing")
    list_path.write_text("\n".join(str(path) for path in art_files) + "\n", encoding="utf-8")

    workflows_dir = Path(__file__).resolve().parent.parent
    edep_fcl = workflows_dir / "fcl" / "edep.fcl"
    if not edep_fcl.exists():
        return {
            "ran": False,
            "error": f"Missing edep fcl: {edep_fcl}",
            "input_file_count": len(art_files),
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": str(edep_fcl),
            "wrapper_fcl_path": str(wrapper_fcl_path),
            "nts_output_path": str(nts_output_path),
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
            "primary_minus_edep_fit_line": None,
            "primary_minus_edep_fit_mean_mev": None,
            "primary_minus_edep_fit_fwhm_mev": None,
            "primary_minus_edep_distribution_line": None,
            "primary_minus_edep_distribution_mean_mev": None,
            "primary_minus_edep_distribution_rms_mev": None,
            "tracker_front_fit_line": None,
            "tracker_front_fit_status": None,
            "tracker_front_fit_mpv_mev": None,
            "tracker_front_fit_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_line": None,
            "primary_edep_minus_tracker_front_distribution_mpv_mev": None,
            "primary_edep_minus_tracker_front_distribution_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_efficiency": None,
        }

    wrapper_include_path = Path("Run1BAna") / "workflows" / "fcl" / "edep.fcl"
    wrapper_fcl_path.write_text(
        f"#include \"{wrapper_include_path.as_posix()}\"\n"
        "\n"
        f"services.TFileService.fileName : \"{nts_output_filename}\"\n",
        encoding="utf-8",
    )

    command = ["mu2e", "-c", str(wrapper_fcl_path), "-S", str(list_path)]
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
            "input_file_count": len(art_files),
            "input_list_path": str(list_path),
            "log_path": str(log_path),
            "fcl_path": str(edep_fcl),
            "wrapper_fcl_path": str(wrapper_fcl_path),
            "nts_output_path": str(nts_output_path),
            "mu2e_returncode": None,
            "summary_line": None,
            "events_seen": None,
            "average_calo_energy_mev": None,
            "events_edep_gt_50_mev": None,
            "primary_minus_edep_fit_line": None,
            "primary_minus_edep_fit_mean_mev": None,
            "primary_minus_edep_fit_fwhm_mev": None,
            "primary_minus_edep_distribution_line": None,
            "primary_minus_edep_distribution_mean_mev": None,
            "primary_minus_edep_distribution_rms_mev": None,
            "tracker_front_fit_line": None,
            "tracker_front_fit_status": None,
            "tracker_front_fit_mpv_mev": None,
            "tracker_front_fit_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_line": None,
            "primary_edep_minus_tracker_front_distribution_mpv_mev": None,
            "primary_edep_minus_tracker_front_distribution_fwhm_mev": None,
            "primary_edep_minus_tracker_front_distribution_efficiency": None,
        }

    log_path.write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")

    summary_line = None
    events_seen = None
    average_calo_energy_mev = None
    events_edep_gt_50_mev = None
    primary_minus_edep_fit_line = None
    primary_minus_edep_fit_mean_mev = None
    primary_minus_edep_fit_fwhm_mev = None
    primary_minus_edep_distribution_line = None
    primary_minus_edep_distribution_mean_mev = None
    primary_minus_edep_distribution_rms_mev = None
    tracker_front_fit_line = None
    tracker_front_fit_status = None
    tracker_front_fit_mpv_mev = None
    tracker_front_fit_fwhm_mev = None
    primary_edep_minus_tracker_front_distribution_line = None
    primary_edep_minus_tracker_front_distribution_mpv_mev = None
    primary_edep_minus_tracker_front_distribution_fwhm_mev = None
    primary_edep_minus_tracker_front_distribution_efficiency = None
    for line in (proc.stdout + "\n" + proc.stderr).splitlines():
        match = _EDEP_SUMMARY_PATTERN.search(line)
        if match:
            summary_line = line.strip()
            events_seen = int(match.group(1))
            average_calo_energy_mev = float(match.group(2))
            events_edep_gt_50_mev = int(match.group(3))
        fit_match = _EDEP_FIT_PATTERN.search(line)
        if fit_match:
            primary_minus_edep_fit_line = line.strip()
            primary_minus_edep_fit_mean_mev = float(fit_match.group(2))
            primary_minus_edep_fit_fwhm_mev = float(fit_match.group(3))
        distribution_match = _EDEP_DISTRIBUTION_PATTERN.search(line)
        if distribution_match:
            primary_minus_edep_distribution_line = line.strip()
            primary_minus_edep_distribution_mean_mev = float(distribution_match.group(1))
            primary_minus_edep_distribution_rms_mev = float(distribution_match.group(2))
        tracker_front_fit_match = _TRACKER_FRONT_FIT_PATTERN.search(line)
        if tracker_front_fit_match:
            tracker_front_fit_line = line.strip()
            tracker_front_fit_status = int(tracker_front_fit_match.group(1))
            tracker_front_fit_mpv_mev = float(tracker_front_fit_match.group(2))
            tracker_front_fit_fwhm_mev = float(tracker_front_fit_match.group(3))
        primary_edep_minus_tracker_front_match = _PRIMARY_EDEP_MINUS_TRACKER_FRONT_DISTRIBUTION_PATTERN.search(line)
        if primary_edep_minus_tracker_front_match:
            primary_edep_minus_tracker_front_distribution_line = line.strip()
            primary_edep_minus_tracker_front_distribution_mpv_mev = float(primary_edep_minus_tracker_front_match.group(1))
            primary_edep_minus_tracker_front_distribution_fwhm_mev = float(primary_edep_minus_tracker_front_match.group(2))
            primary_edep_minus_tracker_front_distribution_efficiency = float(primary_edep_minus_tracker_front_match.group(3))

    return {
        "ran": proc.returncode == 0,
        "error": None if proc.returncode == 0 else f"mu2e exited with code {proc.returncode}",
        "input_file_count": len(art_files),
        "input_list_path": str(list_path),
        "log_path": str(log_path),
        "fcl_path": str(edep_fcl),
        "wrapper_fcl_path": str(wrapper_fcl_path),
        "nts_output_path": str(nts_output_path),
        "mu2e_returncode": proc.returncode,
        "summary_line": summary_line,
        "events_seen": events_seen,
        "average_calo_energy_mev": average_calo_energy_mev,
        "events_edep_gt_50_mev": events_edep_gt_50_mev,
        "primary_minus_edep_fit_line": primary_minus_edep_fit_line,
        "primary_minus_edep_fit_mean_mev": primary_minus_edep_fit_mean_mev,
        "primary_minus_edep_fit_fwhm_mev": primary_minus_edep_fit_fwhm_mev,
        "primary_minus_edep_distribution_line": primary_minus_edep_distribution_line,
        "primary_minus_edep_distribution_mean_mev": primary_minus_edep_distribution_mean_mev,
        "primary_minus_edep_distribution_rms_mev": primary_minus_edep_distribution_rms_mev,
        "tracker_front_fit_line": tracker_front_fit_line,
        "tracker_front_fit_status": tracker_front_fit_status,
        "tracker_front_fit_mpv_mev": tracker_front_fit_mpv_mev,
        "tracker_front_fit_fwhm_mev": tracker_front_fit_fwhm_mev,
        "primary_edep_minus_tracker_front_distribution_line": primary_edep_minus_tracker_front_distribution_line,
        "primary_edep_minus_tracker_front_distribution_mpv_mev": primary_edep_minus_tracker_front_distribution_mpv_mev,
        "primary_edep_minus_tracker_front_distribution_fwhm_mev": primary_edep_minus_tracker_front_distribution_fwhm_mev,
        "primary_edep_minus_tracker_front_distribution_efficiency": primary_edep_minus_tracker_front_distribution_efficiency,
    }


def _run_mubeam_edep_analysis(run_dir: Path, run_tag: str = "Run1B") -> dict:
    flash_glob = f"job_*/dts.mu2e.MuBeamFlash.{run_tag}.*_*.art"
    nts_filename = f"nts.owner.edep.mubeam.{run_tag.lower()}.root"
    return _run_edep_analysis_for_files(
        run_dir,
        sorted(run_dir.glob(flash_glob)),
        analysis_name="mubeam",
        input_label="MuBeamFlash",
        list_filename="mubeam_flash_art_files.txt",
        log_filename="edep_analysis.log",
        wrapper_fcl_filename="edep_analysis.fcl",
        nts_output_filename=nts_filename,
    )


def _run_elebeam_edep_analysis(run_dir: Path) -> dict:
    return _run_edep_analysis_for_files(
        run_dir,
        sorted(run_dir.glob("job_*/dts.mu2e.EleBeamFlash.Run1B.*_*.art")),
        analysis_name="elebeam",
        input_label="EleBeamFlash",
        list_filename="elebeam_flash_art_files.txt",
        log_filename="edep_analysis_elebeam.log",
        wrapper_fcl_filename="edep_analysis_elebeam.fcl",
        nts_output_filename="nts.owner.edep.elebeam.root",
    )


def _run_mustop_edep_analyses(run_dir: Path, run_tag: str = "Run1B") -> dict[str, dict]:
    tag_priority = [run_tag]
    for tag in ("Run1A", "Run1B"):
        if tag not in tag_priority:
            tag_priority.append(tag)

    sample_pattern_by_name = {
        "ce": "dts.mu2e.CeEndpoint",
        "ce_plus": "dts.mu2e.CePlusEndpoint",
        "flat_gamma": "dts.mu2e.FlatGamma",
    }

    analyses: dict[str, dict] = {}
    for sample_name, sample_info in _MUSTOP_EDEP_SAMPLES.items():
        selected_files: list[Path] = []
        selected_tag = run_tag
        sample_prefix = sample_pattern_by_name[sample_name]

        for tag in tag_priority:
            candidate_files = sorted(run_dir.glob(f"job_*/{sample_prefix}.{tag}.*_*.art"))
            if candidate_files:
                selected_files = candidate_files
                selected_tag = tag
                break

        analyses[sample_name] = _run_edep_analysis_for_files(
            run_dir,
            selected_files,
            analysis_name=sample_name,
            input_label=sample_info["label"],
            list_filename=f"{sample_name}_art_files.txt",
            log_filename=f"edep_analysis_{sample_name}.log",
            wrapper_fcl_filename=f"edep_analysis_{sample_name}.fcl",
            nts_output_filename=f"nts.owner.edep.{sample_name}.{selected_tag.lower()}.root",
        )
    return analyses


def _run_mustop_pileup_edep_analysis(run_dir: Path) -> dict:
    return _run_edep_analysis_for_files(
        run_dir,
        sorted(run_dir.glob("job_*/dts.mu2e.MuStopPileup.Run1B.*_*.art")),
        analysis_name="mustop_pileup",
        input_label="Pileup",
        list_filename="pileup_art_files.txt",
        log_filename="edep_analysis_pileup.log",
        wrapper_fcl_filename="edep_analysis_pileup.fcl",
        nts_output_filename="nts.owner.edep.pileup.root",
    )


def _collect_job_status(run_dir: Path) -> tuple[dict[str, int], dict[str, int], list[dict], list[str]]:
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

    return output_counts, output_bytes, status_rows, warnings


def _build_mubeam_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    target_al_stats = _extract_target_al_entries(run_dir)
    art_event_stats = _extract_art_event_counts(run_dir)
    edep_stats = _run_mubeam_edep_analysis(run_dir, run_tag="Run1B")
    input_efficiency = _detect_mubeam_input_efficiency(run_dir)

    if event_stats["total_events"] and event_stats["total_events"] > 0:
        target_al_per_event = target_al_stats["total_target_al_entries"] / event_stats["total_events"]
        calo_per_event = target_al_stats["total_calo_entries"] / event_stats["total_events"]
    else:
        target_al_per_event = None
        calo_per_event = None

    art_events_per_simulated = {}
    art_events_absolute_efficiency = {}
    if event_stats["total_events"] and event_stats["total_events"] > 0:
        for file_type, count in art_event_stats["events_by_type"].items():
            rel_eff = count / event_stats["total_events"]
            art_events_per_simulated[file_type] = rel_eff
            if input_efficiency["correction_factor"] is not None:
                art_events_absolute_efficiency[file_type] = rel_eff * input_efficiency["correction_factor"]

    muminus_stops_files = sorted((run_dir / "mu_stops_job").glob("sim.mu2e.MuminusStopsCat.*.art"))
    muminus_stops_events: int | None = None
    if muminus_stops_files:
        counts = [_count_art_file_events(p) for p in muminus_stops_files]
        if all(c is not None for c in counts):
            muminus_stops_events = sum(counts)  # type: ignore[arg-type]

    return {
        "stage": "mubeam",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "input_efficiency": input_efficiency,
        "target_al_analysis": {
            **target_al_stats,
            "target_al_entries_per_simulated_event": target_al_per_event,
            "calo_entries_per_simulated_event": calo_per_event,
            "target_al_entries_absolute_efficiency": (
                target_al_per_event * input_efficiency["correction_factor"]
                if target_al_per_event is not None and input_efficiency["correction_factor"] is not None
                else None
            ),
            "calo_entries_absolute_efficiency": (
                calo_per_event * input_efficiency["correction_factor"]
                if calo_per_event is not None and input_efficiency["correction_factor"] is not None
                else None
            ),
        },
        "art_event_analysis": {
            **art_event_stats,
            "events_per_simulated_event": art_events_per_simulated,
            "absolute_efficiency_by_type": art_events_absolute_efficiency,
        },
        "muminus_stops_events": muminus_stops_events,
        "edep_analysis": edep_stats,
        "jobs": status_rows,
        "warnings": warnings,
    }


def _build_elebeam_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    art_event_stats = _extract_art_event_counts(run_dir)
    edep_stats = _run_elebeam_edep_analysis(run_dir)
    input_efficiency = _detect_mubeam_input_efficiency(run_dir)

    art_events_per_simulated = {}
    art_events_absolute_efficiency = {}
    if event_stats["total_events"] and event_stats["total_events"] > 0:
        for file_type, count in art_event_stats["events_by_type"].items():
            rel_eff = count / event_stats["total_events"]
            art_events_per_simulated[file_type] = rel_eff
            if input_efficiency["correction_factor"] is not None:
                art_events_absolute_efficiency[file_type] = rel_eff * input_efficiency["correction_factor"]

    return {
        "stage": "elebeam",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "input_efficiency": input_efficiency,
        "art_event_analysis": {
            **art_event_stats,
            "events_per_simulated_event": art_events_per_simulated,
            "absolute_efficiency_by_type": art_events_absolute_efficiency,
        },
        "edep_analysis": edep_stats,
        "jobs": status_rows,
        "warnings": warnings,
    }


def _build_mustop_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    art_event_stats = _extract_art_event_counts(run_dir)

    return {
        "stage": "mustop",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "art_event_analysis": art_event_stats,
        "edep_analysis_by_sample": _run_mustop_edep_analyses(run_dir, run_tag="Run1B"),
        "jobs": status_rows,
        "warnings": warnings,
    }


def _build_run1a_mubeam_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    target_al_stats = _extract_target_al_entries(run_dir)
    art_event_stats = _extract_art_event_counts(run_dir)
    edep_stats = _run_mubeam_edep_analysis(run_dir, run_tag="Run1A")
    input_efficiency = _detect_mubeam_input_efficiency(run_dir)

    if event_stats["total_events"] and event_stats["total_events"] > 0:
        target_al_per_event = target_al_stats["total_target_al_entries"] / event_stats["total_events"]
        calo_per_event = target_al_stats["total_calo_entries"] / event_stats["total_events"]
    else:
        target_al_per_event = None
        calo_per_event = None

    art_events_per_simulated = {}
    art_events_absolute_efficiency = {}
    if event_stats["total_events"] and event_stats["total_events"] > 0:
        for file_type, count in art_event_stats["events_by_type"].items():
            rel_eff = count / event_stats["total_events"]
            art_events_per_simulated[file_type] = rel_eff
            if input_efficiency["correction_factor"] is not None:
                art_events_absolute_efficiency[file_type] = rel_eff * input_efficiency["correction_factor"]

    muminus_stops_files = sorted((run_dir / "mu_stops_job").glob("sim.mu2e.MuminusStopsCat.*.art"))
    muminus_stops_events: int | None = None
    if muminus_stops_files:
        counts = [_count_art_file_events(p) for p in muminus_stops_files]
        if all(c is not None for c in counts):
            muminus_stops_events = sum(counts)  # type: ignore[arg-type]

    return {
        "stage": "run1a_mubeam",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "input_efficiency": input_efficiency,
        "target_al_analysis": {
            **target_al_stats,
            "target_al_entries_per_simulated_event": target_al_per_event,
            "calo_entries_per_simulated_event": calo_per_event,
            "target_al_entries_absolute_efficiency": (
                target_al_per_event * input_efficiency["correction_factor"]
                if target_al_per_event is not None and input_efficiency["correction_factor"] is not None
                else None
            ),
            "calo_entries_absolute_efficiency": (
                calo_per_event * input_efficiency["correction_factor"]
                if calo_per_event is not None and input_efficiency["correction_factor"] is not None
                else None
            ),
        },
        "art_event_analysis": {
            **art_event_stats,
            "events_per_simulated_event": art_events_per_simulated,
            "absolute_efficiency_by_type": art_events_absolute_efficiency,
        },
        "muminus_stops_events": muminus_stops_events,
        "edep_analysis": edep_stats,
        "jobs": status_rows,
        "warnings": warnings,
    }


def _build_run1a_mustops_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    art_event_stats = _extract_art_event_counts(run_dir)

    return {
        "stage": "run1a_mustops",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "art_event_analysis": art_event_stats,
        "edep_analysis_by_sample": _run_mustop_edep_analyses(run_dir, run_tag="Run1A"),
        "jobs": status_rows,
        "warnings": warnings,
    }


def _build_mustop_pileup_summary(run_dir: Path) -> dict:
    output_counts, output_bytes, status_rows, warnings = _collect_job_status(run_dir)

    total_jobs = len(status_rows)
    completed_jobs = sum(1 for row in status_rows if row.get("returncode") == 0)
    failed_jobs = sum(1 for row in status_rows if row.get("returncode") not in (0, None))
    unknown_jobs = sum(1 for row in status_rows if row.get("returncode") is None)

    event_stats = _compute_total_simulated_events(status_rows)
    art_event_stats = _extract_art_event_counts(run_dir)

    return {
        "stage": "mustop_pileup",
        "run_dir": str(run_dir),
        "total_jobs": total_jobs,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "unknown_jobs": unknown_jobs,
        "output_counts": output_counts,
        "output_total_bytes": output_bytes,
        "simulation_events": event_stats,
        "art_event_analysis": art_event_stats,
        "edep_analysis": _run_mustop_pileup_edep_analysis(run_dir),
        "jobs": status_rows,
        "warnings": warnings,
    }


def build_summary(run_dir: Path, stage: str) -> dict:
    if stage == "mubeam":
        return _build_mubeam_summary(run_dir)
    if stage == "elebeam":
        return _build_elebeam_summary(run_dir)
    if stage == "mustop":
        return _build_mustop_summary(run_dir)
    if stage == "mustop_pileup":
        return _build_mustop_pileup_summary(run_dir)
    if stage == "run1a_mubeam":
        return _build_run1a_mubeam_summary(run_dir)
    if stage == "run1a_mustops":
        return _build_run1a_mustops_summary(run_dir)
    raise ValueError(f"Unsupported stage: {stage}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=_STAGES,
        default="mubeam",
        help="Workflow stage to analyze (default: mubeam)",
    )
    parser.add_argument(
        "--run-dir",
        default=None,
        help="Directory containing job_* folders (required unless --from-json is used)",
    )
    parser.add_argument(
        "--from-json",
        action="store_true",
        help="Load <run_dir>/analysis_summary.json and print results without recomputing",
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
    stage = summary.get("stage", "mubeam")
    print("-----")
    print(f"Stage: {stage}")
    print(f"Jobs: {summary['completed_jobs']}/{summary['total_jobs']} completed")
    print(f"Failed: {summary['failed_jobs']}, Unknown: {summary['unknown_jobs']}")
    print("Output counts:")
    for key, value in sorted(summary["output_counts"].items()):
        print(f"  {key}: {value}")

    if stage in ("mustop", "run1a_mustops"):
        art_events = summary["art_event_analysis"]

        print("\nArt file event counts:")
        print(f"  Files analyzed: {art_events['files_analyzed']}")
        print(f"  Total events: {art_events['total_events']}")
        if art_events["events_by_type"]:
            print("  Events by type:")
            for file_type, count in sorted(art_events["events_by_type"].items()):
                print(f"    {file_type}: {count}")

        print("\nEdep reprocessing summary:")
        for sample_name, edep in summary["edep_analysis_by_sample"].items():
            print(f"  Sample: {sample_name}")
            print(f"    Input files: {edep['input_file_count']}")
            print(f"    Input list file: {edep['input_list_path']}")
            print(f"    Edep log: {edep['log_path']}")
            if edep["summary_line"]:
                print(f"    {edep['summary_line']}")
            else:
                print("    EdepAna summary line not found in mu2e output")
            if edep.get("primary_minus_edep_fit_line"):
                print(f"    {edep['primary_minus_edep_fit_line']}")
            if edep.get("primary_minus_edep_distribution_line"):
                print(f"    {edep['primary_minus_edep_distribution_line']}")
            if edep.get("tracker_front_fit_line"):
                print(f"    {edep['tracker_front_fit_line']}")
            if edep.get("primary_edep_minus_tracker_front_distribution_line"):
                print(f"    {edep['primary_edep_minus_tracker_front_distribution_line']}")
            if edep["error"]:
                print(f"    Note: {edep['error']}")
        return

    if stage == "mustop_pileup":
        art_events = summary["art_event_analysis"]

        print("\nArt file event counts:")
        print(f"  Files analyzed: {art_events['files_analyzed']}")
        print(f"  Total events: {art_events['total_events']}")
        if art_events["events_by_type"]:
            print("  Events by type:")
            for file_type, count in sorted(art_events["events_by_type"].items()):
                print(f"    {file_type}: {count}")

        edep = summary["edep_analysis"]
        print("\nEdep reprocessing summary:")
        print(f"  Input files: {edep['input_file_count']}")
        print(f"  Input list file: {edep['input_list_path']}")
        print(f"  Edep log: {edep['log_path']}")
        if edep["summary_line"]:
            print(f"  {edep['summary_line']}")
        else:
            print("  EdepAna summary line not found in mu2e output")
        if edep["error"]:
            print(f"  Note: {edep['error']}")
        return

    if stage == "elebeam":
        nFlash = summary["output_counts"].get("EleFlashOutput", 0)
        input_efficiency = summary.get("input_efficiency", {})
        sim_events = summary["simulation_events"]
        art_events = summary["art_event_analysis"]
        edep = summary["edep_analysis"]

        if input_efficiency.get("correction_factor") is not None:
            print(
                "Input efficiency correction: "
                f"{input_efficiency['correction_factor']:.8g} "
                f"(from {input_efficiency.get('input_fcl')})"
            )
        else:
            print("Input efficiency correction: unavailable")

        print("\nArt file event counts:")
        print(f"  Files analyzed: {art_events['files_analyzed']}")
        print(f"  Total events: {art_events['total_events']}")
        if art_events["events_by_type"]:
            print("  Events by type:")
            for file_type, count in sorted(art_events["events_by_type"].items()):
                per_sim = art_events["events_per_simulated_event"].get(file_type, 0)
                abs_eff = art_events.get("absolute_efficiency_by_type", {}).get(file_type)
                if abs_eff is not None:
                    print(f"    {file_type}: {count} ({per_sim:.8g} per simulated event, {abs_eff:.8g} absolute)")
                else:
                    print(f"    {file_type}: {count} ({per_sim:.8g} per simulated event)")

        print("\nEdep reprocessing summary:")
        print(f"  Input EleBeamFlash files: {edep['input_file_count']}")
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
        return

    nFlash = summary["output_counts"].get("FlashOutput", 0)
    input_efficiency = summary.get("input_efficiency", {})

    target_al = summary["target_al_analysis"]
    sim_events = summary["simulation_events"]
    art_events = summary["art_event_analysis"]

    if input_efficiency.get("correction_factor") is not None:
        print(
            "Input efficiency correction: "
            f"{input_efficiency['correction_factor']:.8g} "
            f"(from {input_efficiency.get('input_fcl')})"
        )
    else:
        print("Input efficiency correction: unavailable")

    print("TargetMuonFinder/stopmat summary:")
    print(f"  ROOT files analyzed: {target_al['files_analyzed']}")
    print(f"  Total StoppingTarget_Al entries: {target_al['total_target_al_entries']}")
    print(
        "  Total calo entries "
        "(G4_CESIUM_IODIDE + CarbonFiber + AluminumHoneycomb): "
        f"{target_al.get('total_calo_entries', 0.0)}"
    )

    if sim_events["total_events"] is None:
        print("  Total simulated events: unknown (not all jobs had '-n' in command)")
        print("  StoppingTarget_Al per simulated event: unavailable")
        print("  Calo entries per simulated event: unavailable")
        print("  StoppingTarget_Al absolute efficiency: unavailable")
        print("  Calo entries absolute efficiency: unavailable")
    else:
        print(f"  Total simulated events: {sim_events['total_events']}")
        print(
            "  StoppingTarget_Al per simulated event: "
            f"{target_al['target_al_entries_per_simulated_event']:.8g}"
        )
        print(
            "  Calo entries per simulated event: "
            f"{target_al.get('calo_entries_per_simulated_event', 0.0):.8g}"
        )
        if target_al.get("target_al_entries_absolute_efficiency") is not None:
            print(
                "  StoppingTarget_Al absolute efficiency: "
                f"{target_al['target_al_entries_absolute_efficiency']:.8g}"
            )
        else:
            print("  StoppingTarget_Al absolute efficiency: unavailable")
        if target_al.get("calo_entries_absolute_efficiency") is not None:
            print(
                "  Calo entries absolute efficiency: "
                f"{target_al['calo_entries_absolute_efficiency']:.8g}"
            )
        else:
            print("  Calo entries absolute efficiency: unavailable")

    if target_al["error"]:
        print(f"  Note: {target_al['error']}")

    print("\nArt file event counts:")
    print(f"  Files analyzed: {art_events['files_analyzed']}")
    print(f"  Total events: {art_events['total_events']}")
    if art_events["events_by_type"]:
        print("  Events by type:")
        for file_type, count in sorted(art_events["events_by_type"].items()):
            per_sim = art_events["events_per_simulated_event"].get(file_type, 0)
            abs_eff = art_events.get("absolute_efficiency_by_type", {}).get(file_type)
            if abs_eff is not None:
                print(f"    {file_type}: {count} ({per_sim:.8g} per simulated event, {abs_eff:.8g} absolute)")
            else:
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
        if not args.run_dir:
            raise SystemExit("--run-dir is required when --from-json is used")

        run_dir = Path(args.run_dir).resolve()
        if not run_dir.exists() or not run_dir.is_dir():
            raise SystemExit(f"Run directory does not exist: {run_dir}")

        json_path = run_dir / "analysis_summary.json"
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

        summary = build_summary(run_dir, args.stage)
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
