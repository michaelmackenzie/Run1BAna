#!/usr/bin/env python3
"""Run parallel mu2e jobs for a selected configuration and summarize outputs."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


_STAGES = (
    "mubeam",
    "elebeam",
    "mustop",
    "mustop_pileup",
    "run1a_mubeam",
    "run1a_mustops",
    "final",
    "all",
    "summary",
)
_MUSTOP_MODES = ("ce", "ce_plus", "flat_gamma")
_GEN_RESTRICTION_FACTOR = (1.0 - 0.95) / 2.0  # cos(theta) generation restriction correction
_MUON_STOP_PRESCALE_CORRECTION = 10.0
_DOUBLE_EDEP_PATTERNS = {
    "pot_per_event_average": re.compile(r"N\(POT\)\s*/\s*event average:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"),
    "single_edep_efficiency_per_pot": re.compile(r"Efficiency for single edep / POT:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"),
    "expected_per_event_single_edep": re.compile(r"Expected per event for single edep:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"),
    "expected_per_event_double_edep": re.compile(r"Expected per event for double edep:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"),
    "expected_per_event_triple_edep": re.compile(r"Expected per event for triple edep:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"),
}
_DOUBLE_EDEP_THRESHOLDS_PATTERN = re.compile(
    r"Estimated rates per event:\s*"
    r"E\(50\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"E\(70\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"E\(80\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"E\(90\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)
_ROUGH_SENSITIVITY_PATTERN = re.compile(
    r"Signal MPV\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"FWHM\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"signal rate\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"background rate\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+"
    r"s/sqrt\(b\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)


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
        help=(
            "Number of jobs to launch "
            "(required for stage mubeam/elebeam/mustop_pileup/run1a_mubeam/all)"
        ),
    )
    parser.add_argument(
        "--events-per-job",
        type=int,
        default=0,
        help="Number of events per job passed as '-n <events>' to mu2e (required for mubeam/mustop_pileup/all)",
    )
    parser.add_argument(
        "--run1a-mubeam-events-per-job",
        type=int,
        default=5000,
        help="Number of events per run1a_mubeam job (default: 5000)",
    )
    parser.add_argument(
        "--elebeam-events-per-job",
        type=int,
        default=0,
        help="Number of events per elebeam job (default: 0, uses --events-per-job if not specified)",
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
        "--run1a-mubeam-run-dir",
        default=None,
        help=(
            "Directory containing run1a_mubeam job_* outputs used as TargetStopResampler inputs "
            "for --stage run1a_mustops (default: latest run1a_mubeam_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--elebeam-run-dir",
        default=None,
        help=(
            "Directory containing elebeam job_* outputs for --stage final/summary "
            "(default: latest elebeam_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--mustop-run-dir",
        default=None,
        help=(
            "Directory containing mustop job_* outputs for --stage summary "
            "(default: latest mustop_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--run1a-mustops-run-dir",
        default=None,
        help=(
            "Directory containing run1a_mustops job_* outputs for --stage summary "
            "(default: latest run1a_mustops_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--mustop-pileup-run-dir",
        default=None,
        help=(
            "Directory containing mustop_pileup job_* outputs for --stage summary "
            "(default: latest mustop_pileup_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--final-run-dir",
        default=None,
        help=(
            "Directory containing final stage outputs for --stage summary "
            "(default: latest final_* under run-root/config_version)"
        ),
    )
    parser.add_argument(
        "--mustop-jobs-per-mode",
        type=int,
        default=10,
        help="Number of jobs to launch for each mustop mode (ce/ce_plus/flat_gamma) (default: 10)",
    )
    parser.add_argument(
        "--mustop-events-per-job",
        type=int,
        default=5000,
        help="Number of events per mustop job (default: 5000)",
    )
    parser.add_argument(
        "--mustop-pileup-events-per-job",
        type=int,
        default=100000,
        help="Number of events per mustop_pileup job (default: 100000)",
    )
    return parser.parse_args()


def _find_latest_stage_run(run_root: Path, config_version: str, stage: str) -> Path | None:
    stage_parent = run_root / config_version
    if not stage_parent.exists() or not stage_parent.is_dir():
        return None

    # Match only canonical stage timestamp directories, e.g. "mustop_YYYYMMDD_HHMMSS".
    # This avoids prefix collisions such as "mustop" matching "mustop_pileup_*".
    stage_dir_pattern = re.compile(rf"^{re.escape(stage)}_\d{{8}}_\d{{6}}$")
    candidates = sorted(
        path for path in stage_parent.iterdir()
        if path.is_dir() and stage_dir_pattern.match(path.name)
    )
    return candidates[-1] if candidates else None


def _collect_target_stop_files(run_dir: Path) -> list[Path]:
    files = sorted(run_dir.glob("job_*/sim.mu2e.TargetStops.Run1A.*_*.art"))
    files.extend(sorted(run_dir.glob("job_*/sim.mu2e.TargetStops.Run1B.*_*.art")))
    return sorted(files)


def _collect_muminus_stop_files(mubeam_run_dir: Path) -> list[Path]:
    run_dir = mubeam_run_dir / "mu_stops_job"
    files = sorted(run_dir.glob("sim.mu2e.MuminusStopsCat.Run1A.*_*.art"))
    files.extend(sorted(run_dir.glob("sim.mu2e.MuminusStopsCat.Run1B.*_*.art")))
    return sorted(files)


def _format_fhicl_string_list(paths: list[Path], indent: str = "    ") -> str:
    return ",\n".join(f'{indent}"{path}"' for path in paths)


def _run_mu_stops_job(
    run_dir: Path,
    workflows_dir: Path,
    config_version: str,
    beam_fragment: str,
    mu2e_command: str,
    env: dict,
    dry_run: bool,
) -> JobResult:
    target_stop_files = _collect_target_stop_files(run_dir)
    if not target_stop_files:
        raise SystemExit(f"No TargetStops files found in {run_dir} for mu_stops job")

    job_dir = run_dir / "mu_stops_job"
    job_dir.mkdir(parents=True, exist_ok=False)

    include_fcl_path = Path("Run1BAna") / "workflows" / config_version / beam_fragment / "mu_stops.fcl"
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


def _load_summary(summary_path: Path) -> dict:
    if not summary_path.exists() or not summary_path.is_file():
        raise SystemExit(f"Missing analysis summary: {summary_path}")
    with summary_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _compute_mustop_pileup_absolute_efficiency(mubeam_summary: dict, mustop_pileup_summary: dict) -> float | None:
    input_corr = mubeam_summary.get("input_efficiency", {}).get("correction_factor")
    sim_total = mubeam_summary.get("simulation_events", {}).get("total_events")
    n_muminus_stops = mubeam_summary.get("muminus_stops_events")
    pileup_sim_total = mustop_pileup_summary.get("simulation_events", {}).get("total_events")
    pileup_seen = mustop_pileup_summary.get("edep_analysis", {}).get("events_seen")

    if (
        input_corr is None
        or sim_total in (None, 0)
        or n_muminus_stops is None
        or pileup_sim_total in (None, 0)
        or pileup_seen is None
    ):
        return None

    stopping_factor = n_muminus_stops / sim_total
    return (
        (pileup_seen / pileup_sim_total)
        * input_corr
        * stopping_factor
        * _MUON_STOP_PRESCALE_CORRECTION
    )


def _compute_mustop_sample_absolute_efficiencies(mubeam_summary: dict, mustop_summary: dict) -> dict[str, float | None]:
    input_corr = mubeam_summary.get("input_efficiency", {}).get("correction_factor")
    sim_total = mubeam_summary.get("simulation_events", {}).get("total_events")
    n_muminus_stops = mubeam_summary.get("muminus_stops_events")
    mustop_sim_total = mustop_summary.get("simulation_events", {}).get("total_events")
    mustop_sim_per_mode = (
        mustop_sim_total / len(_MUSTOP_MODES)
        if mustop_sim_total not in (None, 0)
        else None
    )

    if (
        input_corr is None
        or sim_total in (None, 0)
        or n_muminus_stops is None
        or mustop_sim_per_mode in (None, 0)
    ):
        return {mode: None for mode in _MUSTOP_MODES}

    stopping_factor = n_muminus_stops / sim_total
    scale = (
        input_corr
        * stopping_factor
        * _GEN_RESTRICTION_FACTOR
        * _MUON_STOP_PRESCALE_CORRECTION
        / mustop_sim_per_mode
    )

    result: dict[str, float | None] = {}
    for mode in _MUSTOP_MODES:
        sample_gt50 = mustop_summary.get("edep_analysis_by_sample", {}).get(mode, {}).get("events_edep_gt_50_mev")
        result[mode] = sample_gt50 * scale if sample_gt50 is not None else None
    return result


def _find_double_edep_output_path(run_dir: Path) -> Path | None:
    roots = sorted(run_dir.glob("*.root"))
    if not roots:
        return None
    preferred = [path for path in roots if "double" in path.name.lower()]
    return preferred[-1] if preferred else roots[-1]


def _run_rough_sensitivity_analyses(
    run_dir: Path,
    workflows_dir: Path,
    mustop_summary: dict,
    sample_abs_efficiencies: dict[str, float | None],
    double_edep_output_path: Path | None,
    dry_run: bool,
) -> dict[str, dict]:
    analyses: dict[str, dict] = {}

    for sample in _MUSTOP_MODES:
        edep_root_path_str = mustop_summary.get("edep_analysis_by_sample", {}).get(sample, {}).get("nts_output_path")
        abs_eff = sample_abs_efficiencies.get(sample)
        command_path = run_dir / f"rough_sensitivity_{sample}_command.txt"
        log_path = run_dir / f"rough_sensitivity_{sample}.log"

        if not edep_root_path_str:
            analyses[sample] = {
                "ran": False,
                "returncode": None,
                "error": "Missing mustop sample edep output path",
                "command": None,
                "command_path": str(command_path),
                "log_path": str(log_path),
                "edep_root_path": None,
                "absolute_efficiency": abs_eff,
                "double_edep_output_path": str(double_edep_output_path) if double_edep_output_path else None,
                "sensitivity_line": None,
                "signal_mpv": None,
                "signal_fwhm": None,
                "signal_rate": None,
                "background_rate": None,
                "s_over_sqrt_b": None,
            }
            continue

        edep_root_path = Path(edep_root_path_str)
        if abs_eff is None or double_edep_output_path is None or not edep_root_path.exists():
            analyses[sample] = {
                "ran": False,
                "returncode": None,
                "error": "Missing required rough_sensitivity input(s)",
                "command": None,
                "command_path": str(command_path),
                "log_path": str(log_path),
                "edep_root_path": str(edep_root_path),
                "absolute_efficiency": abs_eff,
                "double_edep_output_path": str(double_edep_output_path) if double_edep_output_path else None,
                "sensitivity_line": None,
                "signal_mpv": None,
                "signal_fwhm": None,
                "signal_rate": None,
                "background_rate": None,
                "s_over_sqrt_b": None,
            }
            continue

        macro_arg = (
            f'"{edep_root_path}", {abs_eff:.16g}, '
            f'"{double_edep_output_path}", "{sample}", "{run_dir}"'
        )
        command = ["root", "-q", "-b", f"scripts/rough_sensitivity.C({macro_arg})"]
        command_path.write_text(shlex.join(command) + "\n", encoding="utf-8")

        if dry_run:
            log_path.write_text("DRY RUN\n", encoding="utf-8")
            analyses[sample] = {
                "ran": False,
                "returncode": 0,
                "error": None,
                "command": command,
                "command_path": str(command_path),
                "log_path": str(log_path),
                "edep_root_path": str(edep_root_path),
                "absolute_efficiency": abs_eff,
                "double_edep_output_path": str(double_edep_output_path),
                "sensitivity_line": None,
                "signal_mpv": None,
                "signal_fwhm": None,
                "signal_rate": None,
                "background_rate": None,
                "s_over_sqrt_b": None,
            }
            continue

        try:
            proc = subprocess.run(
                command,
                cwd=workflows_dir,
                capture_output=True,
                text=True,
                check=False,
            )
        except FileNotFoundError:
            log_path.write_text("root executable not found on PATH\n", encoding="utf-8")
            analyses[sample] = {
                "ran": False,
                "returncode": None,
                "error": "root executable not found on PATH",
                "command": command,
                "command_path": str(command_path),
                "log_path": str(log_path),
                "edep_root_path": str(edep_root_path),
                "absolute_efficiency": abs_eff,
                "double_edep_output_path": str(double_edep_output_path),
                "sensitivity_line": None,
                "signal_mpv": None,
                "signal_fwhm": None,
                "signal_rate": None,
                "background_rate": None,
                "s_over_sqrt_b": None,
            }
            continue

        text = proc.stdout + "\n" + proc.stderr
        log_path.write_text(text, encoding="utf-8")

        sensitivity_line = None
        signal_mpv = None
        signal_fwhm = None
        signal_rate = None
        background_rate = None
        s_over_sqrt_b = None
        for line in text.splitlines():
            match = _ROUGH_SENSITIVITY_PATTERN.search(line.strip())
            if match:
                sensitivity_line = line.strip()
                signal_mpv = float(match.group(1))
                signal_fwhm = float(match.group(2))
                signal_rate = float(match.group(3))
                background_rate = float(match.group(4))
                s_over_sqrt_b = float(match.group(5))

        analyses[sample] = {
            "ran": proc.returncode == 0,
            "returncode": proc.returncode,
            "error": None if proc.returncode == 0 else f"root exited with code {proc.returncode}",
            "command": command,
            "command_path": str(command_path),
            "log_path": str(log_path),
            "edep_root_path": str(edep_root_path),
            "absolute_efficiency": abs_eff,
            "double_edep_output_path": str(double_edep_output_path),
            "sensitivity_line": sensitivity_line,
            "signal_mpv": signal_mpv,
            "signal_fwhm": signal_fwhm,
            "signal_rate": signal_rate,
            "background_rate": background_rate,
            "s_over_sqrt_b": s_over_sqrt_b,
        }

    return analyses


def _run_double_edep_analysis(
    run_dir: Path,
    workflows_dir: Path,
    mubeam_edep_root: Path,
    elebeam_edep_root: Path,
    pileup_edep_root: Path,
    mubeam_flash_abs_eff: float,
    elebeam_flash_abs_eff: float,
    pileup_abs_eff: float,
    dry_run: bool,
) -> dict:
    script_path = Path("scripts") / "double_edep.C"
    macro_arg = (
        f'{{"{mubeam_edep_root}", "{elebeam_edep_root}", "{pileup_edep_root}"}}, '
        f'{{{mubeam_flash_abs_eff:.16g}, {elebeam_flash_abs_eff:.16g}, {pileup_abs_eff:.16g}}}, '
        f'{{"MuBeam flash", "EleBeam flash", "MuStop pileup"}}, "{run_dir}"'
    )
    command = ["root", "-q", "-b", f"{script_path.as_posix()}({macro_arg})"]

    command_path = run_dir / "final_command.txt"
    command_path.write_text(shlex.join(command) + "\n", encoding="utf-8")

    log_path = run_dir / "final.log"
    if dry_run:
        log_path.write_text("DRY RUN\n", encoding="utf-8")
        return {
            "ran": False,
            "returncode": 0,
            "command": command,
            "command_path": str(command_path),
            "log_path": str(log_path),
            "metrics": {},
            "error": None,
        }

    try:
        proc = subprocess.run(
            command,
            cwd=workflows_dir,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        log_path.write_text("root executable not found on PATH\n", encoding="utf-8")
        return {
            "ran": False,
            "returncode": None,
            "command": command,
            "command_path": str(command_path),
            "log_path": str(log_path),
            "metrics": {},
            "error": "root executable not found on PATH",
        }

    text = proc.stdout + "\n" + proc.stderr
    log_path.write_text(text, encoding="utf-8")

    metrics: dict[str, float | None] = {key: None for key in _DOUBLE_EDEP_PATTERNS}
    metric_lines: dict[str, str] = {}
    threshold_rates = {
        "estimated_rate_per_event_e50": None,
        "estimated_rate_per_event_e70": None,
        "estimated_rate_per_event_e80": None,
        "estimated_rate_per_event_e90": None,
    }
    for line in text.splitlines():
        stripped = line.strip()
        for key, pattern in _DOUBLE_EDEP_PATTERNS.items():
            match = pattern.search(stripped)
            if match:
                metrics[key] = float(match.group(1))
                metric_lines[key] = stripped
        threshold_match = _DOUBLE_EDEP_THRESHOLDS_PATTERN.search(stripped)
        if threshold_match:
            threshold_rates = {
                "estimated_rate_per_event_e50": float(threshold_match.group(1)),
                "estimated_rate_per_event_e70": float(threshold_match.group(2)),
                "estimated_rate_per_event_e80": float(threshold_match.group(3)),
                "estimated_rate_per_event_e90": float(threshold_match.group(4)),
            }
            metric_lines["estimated_rates_per_event"] = stripped

    metrics.update(threshold_rates)

    return {
        "ran": proc.returncode == 0,
        "returncode": proc.returncode,
        "command": command,
        "command_path": str(command_path),
        "log_path": str(log_path),
        "metrics": metrics,
        "metric_lines": metric_lines,
        "error": None if proc.returncode == 0 else f"root exited with code {proc.returncode}",
    }


def _print_all_stage_compact_summary(
    mubeam_summary: dict,
    elebeam_summary: dict | None,
    mustop_summary: dict,
    mustop_pileup_summary: dict | None = None,
    final_summary: dict | None = None,
    run1a_mubeam_summary: dict | None = None,
    run1a_mustops_summary: dict | None = None,
) -> None:
    input_corr = mubeam_summary.get("input_efficiency", {}).get("correction_factor")
    target_abs = mubeam_summary.get("target_al_analysis", {}).get("target_al_entries_absolute_efficiency")
    calo_abs = mubeam_summary.get("target_al_analysis", {}).get("calo_entries_absolute_efficiency")
    mubeam_flash_abs_eff = (
        mubeam_summary
        .get("art_event_analysis", {})
        .get("absolute_efficiency_by_type", {})
        .get("FlashOutput")
    )
    elebeam_flash_abs_eff = (
        elebeam_summary.get("art_event_analysis", {}).get("absolute_efficiency_by_type", {}).get("EleFlashOutput")
        if elebeam_summary is not None
        else None
    )
    elebeam_edep_avg = (
        elebeam_summary.get("edep_analysis", {}).get("average_calo_energy_mev")
        if elebeam_summary is not None
        else None
    )
    elebeam_input_corr = (
        elebeam_summary.get("input_efficiency", {}).get("correction_factor")
        if elebeam_summary is not None
        else None
    )
    elebeam_sim_total = (
        elebeam_summary.get("simulation_events", {}).get("total_events")
        if elebeam_summary is not None
        else None
    )
    elebeam_edep_gt50 = (
        elebeam_summary.get("edep_analysis", {}).get("events_edep_gt_50_mev")
        if elebeam_summary is not None
        else None
    )

    sim_total = mubeam_summary.get("simulation_events", {}).get("total_events")
    edep = mubeam_summary.get("edep_analysis", {})
    edep_events_seen = edep.get("events_seen")
    edep_avg = edep.get("average_calo_energy_mev")
    edep_gt50 = edep.get("events_edep_gt_50_mev")

    mubeam_gt50_per_sim_abs = None
    if (
        input_corr is not None
        and sim_total not in (None, 0)
        and edep_events_seen is not None
        and edep_gt50 is not None
    ):
        mubeam_gt50_per_sim_abs = (edep_gt50 / sim_total) * input_corr

    elebeam_gt50_per_sim_abs = None
    if (
        elebeam_input_corr is not None
        and elebeam_sim_total not in (None, 0)
        and elebeam_edep_gt50 is not None
    ):
        elebeam_gt50_per_sim_abs = (elebeam_edep_gt50 / elebeam_sim_total) * elebeam_input_corr

    n_muminus_stops = mubeam_summary.get("muminus_stops_events")
    stopping_factor = (
        n_muminus_stops / sim_total
        if n_muminus_stops is not None and sim_total not in (None, 0)
        else None
    )

    mustop_sim_total = mustop_summary.get("simulation_events", {}).get("total_events")
    mustop_sim_per_mode = (
        mustop_sim_total / len(_MUSTOP_MODES)
        if mustop_sim_total not in (None, 0)
        else None
    )

    mubeam_effective_pot = (
        sim_total / input_corr
        if sim_total not in (None, 0) and input_corr not in (None, 0)
        else None
    )
    elebeam_sim_total = elebeam_summary.get("simulation_events", {}).get("total_events") if elebeam_summary else None
    elebeam_input_corr = elebeam_summary.get("input_efficiency", {}).get("correction_factor") if elebeam_summary else None
    elebeam_effective_pot = (
        elebeam_sim_total / elebeam_input_corr
        if elebeam_sim_total not in (None, 0) and elebeam_input_corr not in (None, 0)
        else None
    )
    muon_stop_input_eff = (
        input_corr * stopping_factor * _MUON_STOP_PRESCALE_CORRECTION
        if input_corr is not None and stopping_factor is not None
        else None
    )
    mustop_effective_pot_per_mode = (
        mustop_sim_per_mode / muon_stop_input_eff
        if mustop_sim_per_mode not in (None, 0) and muon_stop_input_eff not in (None, 0)
        else None
    )

    print("-----")
    print("Compact all-stage summary")
    print(f"  Target absolute muon stopping rate: {target_abs:.8g}" if target_abs is not None else "  Target absolute muon stopping rate: unavailable")
    print(f"  Calorimeter absolute muon stopping rate: {calo_abs:.8g}" if calo_abs is not None else "  Calorimeter absolute muon stopping rate: unavailable")
    print(
        f"  mubeam flash output absolute efficiency: {mubeam_flash_abs_eff:.8g}"
        if mubeam_flash_abs_eff is not None
        else "  mubeam flash output absolute efficiency: unavailable"
    )
    print(
        f"  mubeam average Edep / event: {edep_avg:.8g} MeV"
        if edep_avg is not None
        else "  mubeam average Edep / event: unavailable"
    )
    print(
        f"  mubeam events with Edep > 50 per simulated event (absolute): {mubeam_gt50_per_sim_abs:.8g}"
        if mubeam_gt50_per_sim_abs is not None
        else "  mubeam events with Edep > 50 per simulated event (absolute): unavailable"
    )
    print(
        f"  mubeam effective N(POT) simulated: {mubeam_effective_pot:.8g}"
        if mubeam_effective_pot is not None
        else "  mubeam effective N(POT) simulated: unavailable"
    )
    print(
        f"  elebeam flash output absolute efficiency: {elebeam_flash_abs_eff:.8g}"
        if elebeam_flash_abs_eff is not None
        else "  elebeam flash output absolute efficiency: unavailable"
    )
    print(
        f"  elebeam average Edep / event: {elebeam_edep_avg:.8g} MeV"
        if elebeam_edep_avg is not None
        else "  elebeam average Edep / event: unavailable"
    )
    print(
        f"  elebeam events with Edep > 50 per simulated event (absolute): {elebeam_gt50_per_sim_abs:.8g}"
        if elebeam_gt50_per_sim_abs is not None
        else "  elebeam events with Edep > 50 per simulated event (absolute): unavailable"
    )
    print(
        f"  elebeam effective N(POT) simulated: {elebeam_effective_pot:.8g}"
        if elebeam_effective_pot is not None
        else "  elebeam effective N(POT) simulated: unavailable"
    )
    print(
        f"  mustop effective N(POT) simulated (per mode): {mustop_effective_pot_per_mode:.8g}"
        if mustop_effective_pot_per_mode is not None
        else "  mustop effective N(POT) simulated (per mode): unavailable"
    )

    scale = None
    if input_corr is not None and stopping_factor is not None and mustop_sim_per_mode not in (None, 0):
        scale = (
            input_corr
            * stopping_factor
            * _GEN_RESTRICTION_FACTOR
            * _MUON_STOP_PRESCALE_CORRECTION
            / mustop_sim_per_mode
        )

    for sample in _MUSTOP_MODES:
        sample_stats = mustop_summary.get("edep_analysis_by_sample", {}).get(sample, {})
        sample_seen = sample_stats.get("events_seen")
        sample_gt50 = sample_stats.get("events_edep_gt_50_mev")
        sample_mpv = sample_stats.get("primary_minus_edep_fit_mean_mev")
        sample_fwhm = sample_stats.get("primary_minus_edep_fit_fwhm_mev")

        sample_abs_eff_all = sample_seen * scale if scale is not None and sample_seen is not None else None
        sample_abs_eff_gt50 = sample_gt50 * scale if scale is not None and sample_gt50 is not None else None

        eff_all_str = f"{sample_abs_eff_all:.8g}" if sample_abs_eff_all is not None else "unavailable"
        eff_gt50_str = f"{sample_abs_eff_gt50:.8g}" if sample_abs_eff_gt50 is not None else "unavailable"
        mpv_str = f"{sample_mpv:.4g}" if sample_mpv is not None else "unavailable"
        fwhm_str = f"{sample_fwhm:.4g}" if sample_fwhm is not None else "unavailable"

        print(f"  {sample}: abs eff (all)={eff_all_str}, abs eff (Edep>50)={eff_gt50_str}, MPV={mpv_str} MeV, FWHM={fwhm_str} MeV")

    run1a_target_abs = (
        run1a_mubeam_summary.get("target_al_analysis", {}).get("target_al_entries_absolute_efficiency")
        if run1a_mubeam_summary is not None
        else None
    )
    print(
        f"  run1a target absolute muon stopping rate: {run1a_target_abs:.8g}"
        if run1a_target_abs is not None
        else "  run1a target absolute muon stopping rate: unavailable"
    )

    run1a_input_corr = (
        run1a_mubeam_summary.get("input_efficiency", {}).get("correction_factor")
        if run1a_mubeam_summary is not None
        else None
    )
    run1a_sim_total = (
        run1a_mubeam_summary.get("simulation_events", {}).get("total_events")
        if run1a_mubeam_summary is not None
        else None
    )
    run1a_n_muminus_stops = (
        run1a_mubeam_summary.get("muminus_stops_events")
        if run1a_mubeam_summary is not None
        else None
    )
    run1a_stopping_factor = (
        run1a_n_muminus_stops / run1a_sim_total
        if run1a_n_muminus_stops is not None and run1a_sim_total not in (None, 0)
        else None
    )
    run1a_mustops_sim_total = (
        run1a_mustops_summary.get("simulation_events", {}).get("total_events")
        if run1a_mustops_summary is not None
        else None
    )
    run1a_mustops_sim_per_mode = (
        run1a_mustops_sim_total / len(_MUSTOP_MODES)
        if run1a_mustops_sim_total not in (None, 0)
        else None
    )
    run1a_muon_stop_input_eff = (
        run1a_input_corr * run1a_stopping_factor
        if run1a_input_corr is not None and run1a_stopping_factor is not None
        else None
    )
    run1a_mustops_effective_pot_per_mode = (
        run1a_mustops_sim_per_mode / run1a_muon_stop_input_eff
        if run1a_mustops_sim_per_mode not in (None, 0) and run1a_muon_stop_input_eff not in (None, 0)
        else None
    )
    print(
        f"  run1a_mustops effective N(POT) simulated (per mode): {run1a_mustops_effective_pot_per_mode:.8g}"
        if run1a_mustops_effective_pot_per_mode is not None
        else "  run1a_mustops effective N(POT) simulated (per mode): unavailable"
    )

    run1a_scale = None
    if run1a_input_corr is not None and run1a_stopping_factor is not None and run1a_mustops_sim_per_mode not in (None, 0):
        run1a_scale = (
            run1a_input_corr
            * run1a_stopping_factor
            * _GEN_RESTRICTION_FACTOR
            / run1a_mustops_sim_per_mode
        )

    for sample in _MUSTOP_MODES:
        sample_stats = run1a_mustops_summary.get("edep_analysis_by_sample", {}).get(sample, {}) if run1a_mustops_summary else {}
        sample_seen = sample_stats.get("events_seen")
        sample_gt50 = sample_stats.get("events_edep_gt_50_mev")
        sample_mpv = sample_stats.get("primary_minus_edep_fit_mean_mev")
        sample_fwhm = sample_stats.get("primary_minus_edep_fit_fwhm_mev")

        sample_abs_eff_all = sample_seen * run1a_scale if run1a_scale is not None and sample_seen is not None else None
        sample_abs_eff_gt50 = sample_gt50 * run1a_scale if run1a_scale is not None and sample_gt50 is not None else None

        eff_all_str = f"{sample_abs_eff_all:.8g}" if sample_abs_eff_all is not None else "unavailable"
        eff_gt50_str = f"{sample_abs_eff_gt50:.8g}" if sample_abs_eff_gt50 is not None else "unavailable"
        mpv_str = f"{sample_mpv:.4g}" if sample_mpv is not None else "unavailable"
        fwhm_str = f"{sample_fwhm:.4g}" if sample_fwhm is not None else "unavailable"

        print(
            f"  run1a {sample}: abs eff (all)={eff_all_str}, "
            f"abs eff (Edep>50)={eff_gt50_str}, MPV={mpv_str} MeV, FWHM={fwhm_str} MeV"
        )

    if mustop_pileup_summary is None:
        print("  pileup: summary unavailable")
        return

    pileup_sim_total = mustop_pileup_summary.get("simulation_events", {}).get("total_events")
    pileup_stats = mustop_pileup_summary.get("edep_analysis", {})
    pileup_seen = pileup_stats.get("events_seen")
    pileup_gt50 = pileup_stats.get("events_edep_gt_50_mev")
    pileup_avg = pileup_stats.get("average_calo_energy_mev")

    pileup_scale = None
    if input_corr is not None and stopping_factor is not None and pileup_sim_total not in (None, 0):
        pileup_scale = (
            input_corr
            * stopping_factor
            * _MUON_STOP_PRESCALE_CORRECTION
            / pileup_sim_total
        )

    pileup_abs_eff_all = pileup_seen * pileup_scale if pileup_scale is not None and pileup_seen is not None else None
    pileup_abs_eff_gt50 = pileup_gt50 * pileup_scale if pileup_scale is not None and pileup_gt50 is not None else None

    pileup_eff_all_str = f"{pileup_abs_eff_all:.8g}" if pileup_abs_eff_all is not None else "unavailable"
    pileup_eff_gt50_str = f"{pileup_abs_eff_gt50:.8g}" if pileup_abs_eff_gt50 is not None else "unavailable"
    pileup_avg_str = f"{pileup_avg:.8g}" if pileup_avg is not None else "unavailable"
    pileup_effective_pot = (
        pileup_sim_total / muon_stop_input_eff
        if pileup_sim_total not in (None, 0) and muon_stop_input_eff not in (None, 0)
        else None
    )
    print(
        f"  mustop_pileup effective N(POT) simulated: {pileup_effective_pot:.8g}"
        if pileup_effective_pot is not None
        else "  mustop_pileup effective N(POT) simulated: unavailable"
    )
    pileup_effective_pot_str = f"{pileup_effective_pot:.8g}" if pileup_effective_pot is not None else "unavailable"
    print(
        "  pileup: "
        f"abs eff (all)={pileup_eff_all_str}, "
        f"abs eff (Edep>50)={pileup_eff_gt50_str}, "
        f"avg Edep/event={pileup_avg_str} MeV, "
        f"effective N(POT) simulated={pileup_effective_pot_str}"
    )

    if final_summary is None:
        print("  Double-Edep expected/event (single,double,triple): unavailable")
        return

    metrics = final_summary.get("metrics", {})
    single = metrics.get("expected_per_event_single_edep")
    double = metrics.get("expected_per_event_double_edep")
    triple = metrics.get("expected_per_event_triple_edep")
    e50 = metrics.get("estimated_rate_per_event_e50")
    e70 = metrics.get("estimated_rate_per_event_e70")
    e80 = metrics.get("estimated_rate_per_event_e80")
    e90 = metrics.get("estimated_rate_per_event_e90")
    single_str = f"{single:.8g}" if single is not None else "unavailable"
    double_str = f"{double:.8g}" if double is not None else "unavailable"
    triple_str = f"{triple:.8g}" if triple is not None else "unavailable"
    e50_str = f"{e50:.8g}" if e50 is not None else "unavailable"
    e70_str = f"{e70:.8g}" if e70 is not None else "unavailable"
    e80_str = f"{e80:.8g}" if e80 is not None else "unavailable"
    e90_str = f"{e90:.8g}" if e90 is not None else "unavailable"
    print(f"  Double-Edep expected/event (single,double,triple): {single_str}, {double_str}, {triple_str}")
    print(f"  Double-Edep estimated/event E(50,70,80,90): {e50_str}, {e70_str}, {e80_str}, {e90_str}")

    rough = final_summary.get("rough_sensitivity_by_sample", {}) if final_summary else {}
    for sample in _MUSTOP_MODES:
        stats = rough.get(sample, {})
        sens = stats.get("s_over_sqrt_b")
        sig = stats.get("signal_rate")
        bkg = stats.get("background_rate")
        mpv = stats.get("signal_mpv")
        fwhm = stats.get("signal_fwhm")
        sens_str = f"{sens:.8g}" if sens is not None else "unavailable"
        sig_str = f"{sig:.8g}" if sig is not None else "unavailable"
        bkg_str = f"{bkg:.8g}" if bkg is not None else "unavailable"
        mpv_str = f"{mpv:.8g}" if mpv is not None else "unavailable"
        fwhm_str = f"{fwhm:.8g}" if fwhm is not None else "unavailable"
        print(
            f"  Rough sensitivity {sample}: MPV={mpv_str}, FWHM={fwhm_str}, "
            f"signal rate={sig_str}, background rate={bkg_str}, s/sqrt(b)={sens_str}"
        )


def main() -> int:
    args = parse_args()

    # Default elebeam_events_per_job to events_per_job if not specified
    if args.elebeam_events_per_job <= 0:
        args.elebeam_events_per_job = args.events_per_job

    if args.stage in ("mubeam", "elebeam", "mustop_pileup", "run1a_mubeam", "all") and (args.parallel_jobs is None or args.parallel_jobs <= 0):
        raise SystemExit("parallel_jobs must be > 0 for stage mubeam/elebeam/mustop_pileup/run1a_mubeam/all")
    if args.stage != "summary":
        if args.stage not in ("mustop", "mustop_pileup", "final", "run1a_mubeam", "run1a_mustops") and args.events_per_job <= 0:
            raise SystemExit("events_per_job must be > 0")
        if args.stage == "elebeam" and args.elebeam_events_per_job <= 0:
            raise SystemExit("elebeam_events_per_job (or events_per_job) must be > 0")
        if args.stage == "run1a_mubeam" and args.run1a_mubeam_events_per_job <= 0:
            raise SystemExit("run1a_mubeam_events_per_job must be > 0")
        if args.stage in ("mustop", "run1a_mustops", "all") and args.mustop_events_per_job <= 0:
            raise SystemExit("mustop_events_per_job must be > 0")
        if args.stage in ("mustop_pileup", "all") and args.mustop_pileup_events_per_job <= 0:
            raise SystemExit("mustop_pileup_events_per_job must be > 0")
        if args.seed_start <= 0:
            raise SystemExit("seed_start must be > 0")
        if args.stage in ("mustop", "mustop_pileup", "run1a_mustops", "all") and args.mustop_jobs_per_mode <= 0:
            raise SystemExit("mustop_jobs_per_mode must be > 0")

    script_dir = Path(__file__).resolve().parent
    workflows_dir = script_dir.parent
    extractor_path = script_dir / "extract_analysis_results.py"

    run_root = Path(args.run_root).resolve() if args.run_root else workflows_dir / "runs"

    if args.stage == "summary":
        mubeam_run_dir = (
            Path(args.mubeam_run_dir).resolve()
            if args.mubeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mubeam")
        )
        if mubeam_run_dir is None:
            raise SystemExit("Could not locate mubeam run directory; pass --mubeam-run-dir explicitly")
        mustop_run_dir = (
            Path(args.mustop_run_dir).resolve()
            if args.mustop_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mustop")
        )
        if mustop_run_dir is None:
            raise SystemExit("Could not locate mustop run directory; pass --mustop-run-dir explicitly")
        mubeam_summary = _load_summary(mubeam_run_dir / "analysis_summary.json")
        elebeam_run_dir = (
            Path(args.elebeam_run_dir).resolve()
            if args.elebeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "elebeam")
        )
        elebeam_summary = (
            _load_summary(elebeam_run_dir / "analysis_summary.json")
            if elebeam_run_dir is not None
            else None
        )
        mustop_summary = _load_summary(mustop_run_dir / "analysis_summary.json")
        run1a_mubeam_run_dir = (
            Path(args.run1a_mubeam_run_dir).resolve()
            if args.run1a_mubeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "run1a_mubeam")
        )
        run1a_mubeam_summary = (
            _load_summary(run1a_mubeam_run_dir / "analysis_summary.json")
            if run1a_mubeam_run_dir is not None
            else None
        )
        run1a_mustops_run_dir = (
            Path(args.run1a_mustops_run_dir).resolve()
            if args.run1a_mustops_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "run1a_mustops")
        )
        run1a_mustops_summary = (
            _load_summary(run1a_mustops_run_dir / "analysis_summary.json")
            if run1a_mustops_run_dir is not None
            else None
        )
        mustop_pileup_run_dir = (
            Path(args.mustop_pileup_run_dir).resolve()
            if args.mustop_pileup_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mustop_pileup")
        )
        mustop_pileup_summary = (
            _load_summary(mustop_pileup_run_dir / "analysis_summary.json")
            if mustop_pileup_run_dir is not None
            else None
        )
        final_run_dir = (
            Path(args.final_run_dir).resolve()
            if args.final_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "final")
        )
        final_summary = (
            _load_summary(final_run_dir / "analysis_summary.json")
            if final_run_dir is not None
            else None
        )
        print(f"mubeam: {mubeam_run_dir}")
        if elebeam_run_dir is not None:
            print(f"elebeam: {elebeam_run_dir}")
        print(f"mustop: {mustop_run_dir}")
        if run1a_mubeam_run_dir is not None:
            print(f"run1a_mubeam: {run1a_mubeam_run_dir}")
        if run1a_mustops_run_dir is not None:
            print(f"run1a_mustops: {run1a_mustops_run_dir}")
        if mustop_pileup_run_dir is not None:
            print(f"mustop_pileup: {mustop_pileup_run_dir}")
        if final_run_dir is not None:
            print(f"final: {final_run_dir}")
        _print_all_stage_compact_summary(
            mubeam_summary,
            elebeam_summary,
            mustop_summary,
            mustop_pileup_summary,
            final_summary,
            run1a_mubeam_summary,
            run1a_mustops_summary,
        )
        return 0

    if args.stage == "final":
        mubeam_run_dir = (
            Path(args.mubeam_run_dir).resolve()
            if args.mubeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mubeam")
        )
        if mubeam_run_dir is None:
            raise SystemExit("Could not locate mubeam run directory; pass --mubeam-run-dir explicitly")

        mustop_pileup_run_dir = (
            Path(args.mustop_pileup_run_dir).resolve()
            if args.mustop_pileup_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mustop_pileup")
        )
        if mustop_pileup_run_dir is None:
            raise SystemExit("Could not locate mustop_pileup run directory; pass --mustop-pileup-run-dir explicitly")

        mubeam_summary = _load_summary(mubeam_run_dir / "analysis_summary.json")
        elebeam_run_dir = (
            Path(args.elebeam_run_dir).resolve()
            if args.elebeam_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "elebeam")
        )
        if elebeam_run_dir is None:
            raise SystemExit("Could not locate elebeam run directory; pass --elebeam-run-dir explicitly")
        elebeam_summary = _load_summary(elebeam_run_dir / "analysis_summary.json")
        mustop_pileup_summary = _load_summary(mustop_pileup_run_dir / "analysis_summary.json")
        mustop_run_dir = (
            Path(args.mustop_run_dir).resolve()
            if args.mustop_run_dir
            else _find_latest_stage_run(run_root, args.config_version, "mustop")
        )
        if mustop_run_dir is None:
            raise SystemExit("Could not locate mustop run directory; pass --mustop-run-dir explicitly")
        mustop_summary = _load_summary(mustop_run_dir / "analysis_summary.json")

        mubeam_edep_root = Path(mubeam_summary.get("edep_analysis", {}).get("nts_output_path", ""))
        elebeam_edep_root = Path(elebeam_summary.get("edep_analysis", {}).get("nts_output_path", ""))
        pileup_edep_root = Path(mustop_pileup_summary.get("edep_analysis", {}).get("nts_output_path", ""))
        if not mubeam_edep_root.exists():
            raise SystemExit(f"Missing mubeam edep root file: {mubeam_edep_root}")
        if not elebeam_edep_root.exists():
            raise SystemExit(f"Missing elebeam edep root file: {elebeam_edep_root}")
        if not pileup_edep_root.exists():
            raise SystemExit(f"Missing mustop_pileup edep root file: {pileup_edep_root}")

        mubeam_flash_abs_eff = (
            mubeam_summary
            .get("art_event_analysis", {})
            .get("absolute_efficiency_by_type", {})
            .get("FlashOutput")
        )
        if mubeam_flash_abs_eff is None:
            raise SystemExit("Missing mubeam FlashOutput absolute efficiency in mubeam summary")
        elebeam_flash_abs_eff = (
            elebeam_summary
            .get("art_event_analysis", {})
            .get("absolute_efficiency_by_type", {})
            .get("EleFlashOutput")
        )
        if elebeam_flash_abs_eff is None:
            raise SystemExit("Missing elebeam EleFlashOutput absolute efficiency in elebeam summary")

        pileup_abs_eff = _compute_mustop_pileup_absolute_efficiency(mubeam_summary, mustop_pileup_summary)
        if pileup_abs_eff is None:
            raise SystemExit("Could not compute mustop_pileup absolute efficiency from available summaries")

        if args.final_run_dir:
            run_dir = Path(args.final_run_dir).resolve()
            if not run_dir.exists() or not run_dir.is_dir():
                raise SystemExit(f"final run directory does not exist: {run_dir}")
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = run_root / args.config_version / f"final_{timestamp}"
            run_dir.mkdir(parents=True, exist_ok=False)

        final_result = _run_double_edep_analysis(
            run_dir,
            workflows_dir,
            mubeam_edep_root,
            elebeam_edep_root,
            pileup_edep_root,
            mubeam_flash_abs_eff,
            elebeam_flash_abs_eff,
            pileup_abs_eff,
            args.dry_run,
        )
        double_edep_output_path = _find_double_edep_output_path(run_dir)
        sample_abs_efficiencies = _compute_mustop_sample_absolute_efficiencies(mubeam_summary, mustop_summary)
        rough_sensitivity_by_sample = _run_rough_sensitivity_analyses(
            run_dir,
            workflows_dir,
            mustop_summary,
            sample_abs_efficiencies,
            double_edep_output_path,
            args.dry_run,
        )

        summary = {
            "stage": "final",
            "run_dir": str(run_dir),
            "mubeam_run_dir": str(mubeam_run_dir),
            "elebeam_run_dir": str(elebeam_run_dir),
            "mustop_pileup_run_dir": str(mustop_pileup_run_dir),
            "inputs": {
                "mubeam_flash_edep_root": str(mubeam_edep_root),
                "elebeam_flash_edep_root": str(elebeam_edep_root),
                "mustop_pileup_edep_root": str(pileup_edep_root),
                "mubeam_flash_absolute_efficiency": mubeam_flash_abs_eff,
                "elebeam_flash_absolute_efficiency": elebeam_flash_abs_eff,
                "mustop_pileup_absolute_efficiency": pileup_abs_eff,
            },
            "double_edep_analysis": final_result,
            "double_edep_output_path": str(double_edep_output_path) if double_edep_output_path else None,
            "rough_sensitivity_by_sample": rough_sensitivity_by_sample,
            "metrics": final_result.get("metrics", {}),
        }

        summary_path = run_dir / "analysis_summary.json"
        with summary_path.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)

        print(f"Final-stage run directory: {run_dir}")
        print(f"Analysis summary: {summary_path}")
        if final_result.get("error"):
            print("Final-stage analysis failed", file=sys.stderr)
            return 1
        return 0

    if args.stage == "all":
        this_script = Path(__file__).resolve()
        base_cmd = [
            sys.executable,
            str(this_script),
            args.config_version,
            str(args.parallel_jobs),
            "--events-per-job",
            str(args.events_per_job),
            "--run1a-mubeam-events-per-job",
            str(args.run1a_mubeam_events_per_job),
            "--mu2e-command",
            args.mu2e_command,
            "--seed-start",
            str(args.seed_start),
            "--run-root",
            str(run_root),
            "--mustop-jobs-per-mode",
            str(args.mustop_jobs_per_mode),
            "--mustop-events-per-job",
            str(args.mustop_events_per_job),
            "--mustop-pileup-events-per-job",
            str(args.mustop_pileup_events_per_job),
        ]
        if args.max_workers is not None:
            base_cmd.extend(["--max-workers", str(args.max_workers)])
        if args.dry_run:
            base_cmd.append("--dry-run")

        print("Running stage sequence: mubeam -> elebeam -> mustop -> run1a_mubeam -> run1a_mustops -> mustop_pileup -> final")
        mubeam_cmd = base_cmd + ["--stage", "mubeam"]
        mubeam_run = subprocess.run(mubeam_cmd, check=False)
        if mubeam_run.returncode != 0:
            print("mubeam stage failed in all-stage sequence", file=sys.stderr)
            return mubeam_run.returncode

        mubeam_run_dir = _find_latest_stage_run(run_root, args.config_version, "mubeam")
        if mubeam_run_dir is None:
            raise SystemExit("Could not locate mubeam run directory after mubeam stage completion")

        elebeam_cmd = base_cmd + ["--stage", "elebeam"]
        elebeam_run = subprocess.run(elebeam_cmd, check=False)
        if elebeam_run.returncode != 0:
            print("elebeam stage failed in all-stage sequence", file=sys.stderr)
            return elebeam_run.returncode

        elebeam_run_dir = _find_latest_stage_run(run_root, args.config_version, "elebeam")
        if elebeam_run_dir is None:
            raise SystemExit("Could not locate elebeam run directory after elebeam stage completion")

        mustop_cmd = base_cmd + ["--stage", "mustop", "--mubeam-run-dir", str(mubeam_run_dir)]
        mustop_run = subprocess.run(mustop_cmd, check=False)
        if mustop_run.returncode != 0:
            print("mustop stage failed in all-stage sequence", file=sys.stderr)
            return mustop_run.returncode

        run1a_mubeam_cmd = base_cmd + ["--stage", "run1a_mubeam"]
        run1a_mubeam_run = subprocess.run(run1a_mubeam_cmd, check=False)
        if run1a_mubeam_run.returncode != 0:
            print("run1a_mubeam stage failed in all-stage sequence", file=sys.stderr)
            return run1a_mubeam_run.returncode

        run1a_mubeam_run_dir = _find_latest_stage_run(run_root, args.config_version, "run1a_mubeam")
        if run1a_mubeam_run_dir is None:
            raise SystemExit("Could not locate run1a_mubeam run directory after run1a_mubeam stage completion")

        run1a_mustops_cmd = base_cmd + [
            "--stage", "run1a_mustops", "--run1a-mubeam-run-dir", str(run1a_mubeam_run_dir)
        ]
        run1a_mustops_run = subprocess.run(run1a_mustops_cmd, check=False)
        if run1a_mustops_run.returncode != 0:
            print("run1a_mustops stage failed in all-stage sequence", file=sys.stderr)
            return run1a_mustops_run.returncode

        mustop_pileup_cmd = base_cmd + ["--stage", "mustop_pileup", "--mubeam-run-dir", str(mubeam_run_dir)]
        mustop_pileup_run = subprocess.run(mustop_pileup_cmd, check=False)
        if mustop_pileup_run.returncode != 0:
            print("mustop_pileup stage failed in all-stage sequence", file=sys.stderr)
            return mustop_pileup_run.returncode

        latest_mustop_pileup = _find_latest_stage_run(run_root, args.config_version, "mustop_pileup")
        if latest_mustop_pileup is None:
            raise SystemExit("Could not locate mustop_pileup run directory after mustop_pileup stage completion")

        final_cmd = base_cmd + [
            "--stage", "final",
            "--mubeam-run-dir", str(mubeam_run_dir),
            "--elebeam-run-dir", str(elebeam_run_dir),
            "--mustop-pileup-run-dir", str(latest_mustop_pileup),
        ]
        final_run = subprocess.run(final_cmd, check=False)
        if final_run.returncode != 0:
            print("final stage failed in all-stage sequence", file=sys.stderr)
            return final_run.returncode

        print(f"All-stage sequence complete. mubeam run: {mubeam_run_dir}")
        print(f"All-stage sequence complete. elebeam run: {elebeam_run_dir}")
        print(f"All-stage sequence complete. run1a_mubeam run: {run1a_mubeam_run_dir}")
        latest_run1a_mustops = _find_latest_stage_run(run_root, args.config_version, "run1a_mustops")
        if latest_run1a_mustops:
            print(f"All-stage sequence complete. run1a_mustops run: {latest_run1a_mustops}")
        latest_mustop = _find_latest_stage_run(run_root, args.config_version, "mustop")
        if latest_mustop:
            print(f"All-stage sequence complete. mustop run: {latest_mustop}")
            if latest_mustop_pileup:
                print(f"All-stage sequence complete. mustop_pileup run: {latest_mustop_pileup}")
            latest_final = _find_latest_stage_run(run_root, args.config_version, "final")
            if latest_final:
                print(f"All-stage sequence complete. final run: {latest_final}")

            mubeam_summary = _load_summary(mubeam_run_dir / "analysis_summary.json")
            elebeam_summary = _load_summary(elebeam_run_dir / "analysis_summary.json")
            mustop_summary = _load_summary(latest_mustop / "analysis_summary.json")
            run1a_mubeam_summary = _load_summary(run1a_mubeam_run_dir / "analysis_summary.json")
            run1a_mustops_summary = (
                _load_summary(latest_run1a_mustops / "analysis_summary.json")
                if latest_run1a_mustops
                else None
            )
            mustop_pileup_summary = (
                _load_summary(latest_mustop_pileup / "analysis_summary.json")
                if latest_mustop_pileup
                else None
            )
            latest_final = _find_latest_stage_run(run_root, args.config_version, "final")
            final_summary = (
                _load_summary(latest_final / "analysis_summary.json")
                if latest_final
                else None
            )
            _print_all_stage_compact_summary(
                mubeam_summary,
                elebeam_summary,
                mustop_summary,
                mustop_pileup_summary,
                final_summary,
                run1a_mubeam_summary,
                run1a_mustops_summary,
            )
        return 0

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
    elif args.stage == "run1a_mubeam":
        fcl_path = workflows_dir / args.config_version / "run1a_beam" / "mubeam.fcl"
        include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / "run1a_beam" / "mubeam.fcl"
        if not fcl_path.exists():
            raise SystemExit(f"Missing FCL file: {fcl_path}")
        if not extractor_path.exists():
            raise SystemExit(f"Missing extractor script: {extractor_path}")

        job_specs = [
            {
                "index": index,
                "name": "run1a_mubeam",
                "job_fcl_name": "run1a_mubeam_job.fcl",
                "include_fcl_path": include_fcl_path,
                "fcl_overrides": "",
            }
            for index in range(args.parallel_jobs)
        ]
    elif args.stage == "elebeam":
        fcl_path = workflows_dir / args.config_version / "run1b_beam" / "elebeam.fcl"
        include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / "run1b_beam" / "elebeam.fcl"
        if not fcl_path.exists():
            raise SystemExit(f"Missing FCL file: {fcl_path}")
        if not extractor_path.exists():
            raise SystemExit(f"Missing extractor script: {extractor_path}")

        job_specs = [
            {
                "index": index,
                "name": "elebeam",
                "job_fcl_name": "elebeam_job.fcl",
                "include_fcl_path": include_fcl_path,
                "fcl_overrides": "",
            }
            for index in range(args.parallel_jobs)
        ]
    elif args.stage in ("mustop", "run1a_mustops"):
        input_mubeam_stage = "mubeam" if args.stage == "mustop" else "run1a_mubeam"
        input_mubeam_run_dir_arg = args.mubeam_run_dir if args.stage == "mustop" else args.run1a_mubeam_run_dir
        input_mubeam_run_dir = (
            Path(input_mubeam_run_dir_arg).resolve()
            if input_mubeam_run_dir_arg
            else _find_latest_stage_run(run_root, args.config_version, input_mubeam_stage)
        )

        if input_mubeam_run_dir is None:
            raise SystemExit(
                "Could not find a mubeam run directory. "
                "Provide --mubeam-run-dir/--run1a-mubeam-run-dir or run the corresponding mubeam stage first."
            )
        if not input_mubeam_run_dir.exists() or not input_mubeam_run_dir.is_dir():
            raise SystemExit(f"mubeam run directory does not exist: {input_mubeam_run_dir}")

        muminus_stop_files = _collect_muminus_stop_files(input_mubeam_run_dir)
        if not muminus_stop_files:
            raise SystemExit(
                f"No MuminusStopsCat files found in {input_mubeam_run_dir / 'mu_stops_job'}. "
                "Run the corresponding mubeam stage to completion first."
            )

        print(f"Using mubeam inputs from: {input_mubeam_run_dir}")
        print(f"MuminusStopsCat input files: {len(muminus_stop_files)}")

        job_specs = []
        mustop_fragment = "run1b_mustop" if args.stage == "mustop" else "run1a_mustop"
        for mode in _MUSTOP_MODES:
            mode_fcl_path = workflows_dir / args.config_version / mustop_fragment / f"{mode}.fcl"
            include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / mustop_fragment / f"{mode}.fcl"
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
    else: # mustop_pileup
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

        pileup_fcl_path = workflows_dir / args.config_version / "run1b_mustop" / "pileup.fcl"
        include_fcl_path = Path("Run1BAna") / "workflows" / args.config_version / "run1b_mustop" / "pileup.fcl"
        if not pileup_fcl_path.exists():
            raise SystemExit(f"Missing FCL file: {pileup_fcl_path}")

        muminus_stop_lines = _format_fhicl_string_list(muminus_stop_files)
        overrides = (
            "physics.filters.TargetStopResampler.fileNames: [\n"
            f"{muminus_stop_lines}\n"
            "]\n"
        )
        job_specs = [
            {
                "index": index,
                "name": f"pileup_{index:03d}",
                "job_fcl_name": f"pileup_job_{index:03d}.fcl",
                "include_fcl_path": include_fcl_path,
                "fcl_overrides": overrides,
            }
            for index in range(args.parallel_jobs)
        ]

    if args.stage == "mubeam":
        print(f"FCL: {workflows_dir / args.config_version / 'run1b_beam' / 'mubeam.fcl'}")
    elif args.stage == "run1a_mubeam":
        print(f"FCL: {workflows_dir / args.config_version / 'run1a_beam' / 'mubeam.fcl'}")
    elif args.stage == "elebeam":
        print(f"FCL: {workflows_dir / args.config_version / 'run1b_beam' / 'elebeam.fcl'}")
    elif args.stage == "mustop_pileup":
        print(f"FCL: {workflows_dir / args.config_version / 'run1b_mustop' / 'pileup.fcl'}")
    else:
        mustop_fragment = "run1b_mustop" if args.stage == "mustop" else "run1a_mustop"
        print(
            "FCLs: "
            + ", ".join(
                str(workflows_dir / args.config_version / mustop_fragment / f"{mode}.fcl")
                for mode in _MUSTOP_MODES
            )
        )
    print(f"Run directory: {run_dir}")
    print(f"Stage: {args.stage}")
    if args.stage in ("mustop", "run1a_mustops"):
        print(f"mustop jobs per mode: {args.mustop_jobs_per_mode}")
    if args.stage == "mustop_pileup":
        print(f"mustop_pileup jobs: {args.parallel_jobs}")
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
            if args.stage == "mustop":
                n_events = args.mustop_events_per_job
            elif args.stage == "run1a_mubeam":
                n_events = args.run1a_mubeam_events_per_job
            elif args.stage == "elebeam":
                n_events = args.elebeam_events_per_job
            elif args.stage == "mustop_pileup":
                n_events = args.mustop_pileup_events_per_job
            elif args.stage == "run1a_mustops":
                n_events = args.mustop_events_per_job
            else:
                n_events = args.events_per_job
            command.extend(["-n", str(n_events)])

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

    if args.stage in ("mubeam", "run1a_mubeam"):
        if failed > 0:
            print("Warning: proceeding with available successful mubeam outputs despite failed jobs")

        beam_fragment = "run1b_beam" if args.stage == "mubeam" else "run1a_beam"

        mu_stops_result = _run_mu_stops_job(
            run_dir, workflows_dir, args.config_version,
            beam_fragment,
            args.mu2e_command, env, args.dry_run,
        )
        if mu_stops_result.returncode != 0:
            print("mu_stops job failed", file=sys.stderr)
            return 1

    summary_path = run_dir / "analysis_summary.json"
    extractor_cmd = [
        sys.executable,
        str(extractor_path),
        "--stage",
        args.stage,
        "--run-dir",
        str(run_dir),
        "--output",
        str(summary_path),
        "--pretty",
    ]

    print("Running analysis extractor...")
    extractor_command_path = run_dir / "extractor_command.txt"
    extractor_command_path.write_text(shlex.join(extractor_cmd) + "\n", encoding="utf-8")

    extractor_log_path = run_dir / "extractor.log"
    with extractor_log_path.open("w", encoding="utf-8") as log_file:
        extractor_run = subprocess.run(extractor_cmd, stdout=log_file, stderr=subprocess.STDOUT, check=False)

    if extractor_run.returncode != 0:
        print("Extractor failed", file=sys.stderr)
        return extractor_run.returncode

    print(f"Analysis summary: {summary_path}")
    if failed > 0:
        print("Warning: stage completed with failed jobs; efficiencies are evaluated from successful outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
