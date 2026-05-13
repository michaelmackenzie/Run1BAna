#!/usr/bin/env python3
"""Run parallel mu2e jobs for a selected configuration and summarize outputs."""

from __future__ import annotations

import argparse
import json
import math
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

import extract_analysis_results

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_STAGES = (
    "mubeam",
    "elebeam",
    "mustop",
    "mustop_pileup",
    "run1a_mubeam",
    "run1a_mustops",
    "final",
    "all",
    "minimal",
    "summary",
)

_DEFAULT_MUSTOP_MODES = ("ce", "flat_gamma", "flat_electron")
_OPTIONAL_MUSTOP_MODE = "ce_plus"
_ALL_MUSTOP_MODES = _DEFAULT_MUSTOP_MODES + (_OPTIONAL_MUSTOP_MODE,)

# cos(theta) generation restriction correction: (1 - 0.95) / 2
_GEN_RESTRICTION_FACTOR = (1.0 - 0.95) / 2.0
_MUON_STOP_PRESCALE_CORRECTION = 10.0

# Regex for a floating-point number (with optional exponent).
_FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

_DOUBLE_EDEP_PATTERNS: dict[str, re.Pattern[str]] = {
    "pot_per_event_average": re.compile(
        rf"N\(POT\)\s*/\s*event average:\s*({_FLOAT})"
    ),
    "single_edep_efficiency_per_pot": re.compile(
        rf"Efficiency for single edep / POT:\s*({_FLOAT})"
    ),
    "expected_per_event_single_edep": re.compile(
        rf"Expected per event for single edep:\s*({_FLOAT})"
    ),
    "expected_per_event_double_edep": re.compile(
        rf"Expected per event for double edep:\s*({_FLOAT})"
    ),
    "expected_per_event_triple_edep": re.compile(
        rf"Expected per event for triple edep:\s*({_FLOAT})"
    ),
}

_DOUBLE_EDEP_THRESHOLDS_PATTERN = re.compile(
    rf"Estimated rates per event:\s*"
    rf"E\(50\)\s*=\s*({_FLOAT})\s+"
    rf"E\(70\)\s*=\s*({_FLOAT})\s+"
    rf"E\(80\)\s*=\s*({_FLOAT})\s+"
    rf"E\(90\)\s*=\s*({_FLOAT})"
)

_ROUGH_SENSITIVITY_PATTERN = re.compile(
    rf"Signal MPV\s*=\s*({_FLOAT})\s+"
    rf"FWHM\s*=\s*({_FLOAT})\s+"
    rf"signal rate\s*=\s*({_FLOAT})\s+"
    rf"background rate\s*=\s*({_FLOAT})\s+"
    rf"s/sqrt\(b\)\s*=\s*({_FLOAT})"
)

_ROUGH_RUN1A_SENSITIVITY_PATTERN = re.compile(
    rf"Signal box\s*=\s*\[\s*({_FLOAT})\s*,\s*"
    rf"({_FLOAT})\s*\]\s*MeV/c\s*,\s*"
    rf"signal\s*=\s*({_FLOAT})\s*,\s*"
    rf"dio\s*=\s*({_FLOAT})\s*,\s*"
    rf"cosmic\s*=\s*({_FLOAT})\s*-->\s*"
    rf"bkg\s*=\s*({_FLOAT})\s*,\s*"
    rf"S/sqrt\(B\)\s*=\s*({_FLOAT})"
)

# Stage -> (beam_fragment, fcl_name, job_fcl_prefix)
_BEAM_STAGE_CONFIG: dict[str, tuple[str, str, str]] = {
    "mubeam":      ("run1b_beam", "mubeam.fcl",  "mubeam_job.fcl"),
    "run1a_mubeam": ("run1a_beam", "mubeam.fcl", "run1a_mubeam_job.fcl"),
    "elebeam":     ("run1b_beam", "elebeam.fcl",  "elebeam_job.fcl"),
}


# ---------------------------------------------------------------------------
# Mustop mode helpers
# ---------------------------------------------------------------------------

def _selected_mustop_modes(include_ce_plus: bool, ce_only: bool = False) -> tuple[str, ...]:
    """Return the mustop modes to process based on user flags."""
    if ce_only:
        return ("ce",)
    if include_ce_plus:
        return _DEFAULT_MUSTOP_MODES + (_OPTIONAL_MUSTOP_MODE,)
    return _DEFAULT_MUSTOP_MODES


def _mustop_modes_from_summary(
    summary: dict | None, *, include_ce_plus: bool = False
) -> tuple[str, ...]:
    """Read sample names from a summary dict, with fallback to configured modes."""
    fallback = _selected_mustop_modes(include_ce_plus)
    if not summary:
        return fallback
    sample_names = tuple(summary.get("edep_analysis_by_sample", {}).keys())
    return sample_names if sample_names else fallback


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class JobResult:
    """Outcome of a single mu2e job."""
    index: int
    job_dir: Path
    command: list[str]
    returncode: int
    duration_s: float
    log_path: Path


# ---------------------------------------------------------------------------
# Low-level job helpers
# ---------------------------------------------------------------------------

def _run_one_job(
    index: int,
    command: list[str],
    job_dir: Path,
    env: dict[str, str],
    dry_run: bool,
) -> JobResult:
    """Execute a single job (or fake it under *dry_run*)."""
    log_path = job_dir / "job.log"
    start = time.perf_counter()

    if dry_run:
        log_path.write_text(f"DRY RUN: {shlex.join(command)}\n", encoding="utf-8")
        return JobResult(index=index, job_dir=job_dir, command=command,
                         returncode=0, duration_s=0.0, log_path=log_path)

    with log_path.open("w", encoding="utf-8") as fh:
        proc = subprocess.run(
            command, cwd=job_dir, env=env,
            stdout=fh, stderr=subprocess.STDOUT,
            text=True, check=False,
        )

    duration = time.perf_counter() - start
    return JobResult(index=index, job_dir=job_dir, command=command,
                     returncode=proc.returncode, duration_s=duration,
                     log_path=log_path)


def _write_job_status(result: JobResult) -> None:
    """Persist a compact JSON status file next to the job log."""
    status = {
        "job_index": result.index,
        "job_dir": str(result.job_dir),
        "command": result.command,
        "returncode": result.returncode,
        "duration_s": round(result.duration_s, 3),
        "log_path": str(result.log_path),
    }
    with (result.job_dir / "job_status.json").open("w", encoding="utf-8") as fh:
        json.dump(status, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# File-discovery helpers
# ---------------------------------------------------------------------------

def _find_latest_stage_run(
    run_root: Path, config_version: str, stage: str
) -> Path | None:
    """Return the most-recent timestamped directory for *stage*, or ``None``."""
    parent = run_root / config_version
    if not parent.is_dir():
        return None

    # Match only canonical stage timestamp dirs, e.g. "mustop_YYYYMMDD_HHMMSS".
    # Avoids prefix collisions like "mustop" matching "mustop_pileup_*".
    pat = re.compile(rf"^{re.escape(stage)}_\d{{8}}_\d{{6}}$")
    candidates = sorted(p for p in parent.iterdir() if p.is_dir() and pat.match(p.name))
    return candidates[-1] if candidates else None


def _resolve_stage_dir(
    explicit: str | None, run_root: Path, config_version: str, stage: str,
) -> Path | None:
    """Resolve an explicit path or fall back to the latest timestamped run."""
    if explicit:
        return Path(explicit).resolve()
    return _find_latest_stage_run(run_root, config_version, stage)


def _require_stage_dir(
    explicit: str | None, run_root: Path, config_version: str, stage: str,
) -> Path:
    """Like ``_resolve_stage_dir`` but raises ``SystemExit`` on failure."""
    result = _resolve_stage_dir(explicit, run_root, config_version, stage)
    if result is None:
        raise SystemExit(
            f"Could not locate {stage} run directory; "
            f"pass --{stage.replace('_', '-')}-run-dir explicitly"
        )
    return result


def _collect_target_stop_files(run_dir: Path) -> list[Path]:
    """Glob for TargetStops .art files across all job_* directories."""
    files = sorted(run_dir.glob("job_*/sim.mu2e.TargetStops.Run1A.*_*.art"))
    files.extend(sorted(run_dir.glob("job_*/sim.mu2e.TargetStops.Run1B.*_*.art")))
    return sorted(files)


def _collect_muminus_stop_files(mubeam_run_dir: Path) -> list[Path]:
    """Glob for concatenated MuminusStops .art files."""
    mu_dir = mubeam_run_dir / "mu_stops_job"
    files = sorted(mu_dir.glob("sim.mu2e.MuminusStopsCat.Run1A.*_*.art"))
    files.extend(sorted(mu_dir.glob("sim.mu2e.MuminusStopsCat.Run1B.*_*.art")))
    return sorted(files)


def _find_double_edep_output_path(run_dir: Path) -> Path | None:
    """Find the ROOT output from the double-edep analysis in *run_dir*."""
    roots = sorted(run_dir.glob("*.root"))
    if not roots:
        return None
    preferred = [p for p in roots if "double" in p.name.lower()]
    return preferred[-1] if preferred else roots[-1]


# ---------------------------------------------------------------------------
# FHiCL helpers
# ---------------------------------------------------------------------------

def _format_fhicl_string_list(paths: list[Path], indent: str = "    ") -> str:
    """Format a list of paths as a FHiCL string list body."""
    return ",\n".join(f'{indent}"{p}"' for p in paths)


# ---------------------------------------------------------------------------
# Summary I/O
# ---------------------------------------------------------------------------

def _load_summary(summary_path: Path) -> dict:
    """Load a JSON summary; raise ``SystemExit`` if missing."""
    if not summary_path.is_file():
        raise SystemExit(f"Missing analysis summary: {summary_path}")
    with summary_path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _load_summary_optional(summary_path: Path) -> dict | None:
    """Load a JSON summary, returning ``None`` on any failure."""
    if not summary_path.is_file():
        return None
    try:
        with summary_path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


# ---------------------------------------------------------------------------
# Physics calculations
# ---------------------------------------------------------------------------

def _safe_get(d: dict, *keys, default=None):
    """Nested dict get that returns *default* on any missing key."""
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k, default)
        if d is default:
            return default
    return d


def _compute_mustop_pileup_absolute_efficiency(
    mubeam_summary: dict, pileup_summary: dict,
) -> float | None:
    """Absolute efficiency for mustop pileup hits given beam and pileup summaries."""
    input_corr = _safe_get(mubeam_summary, "input_efficiency", "correction_factor")
    sim_total = _safe_get(mubeam_summary, "simulation_events", "total_events")
    n_stops = mubeam_summary.get("muminus_stops_events")
    pileup_sim = _safe_get(pileup_summary, "simulation_events", "total_events")
    pileup_seen = _safe_get(pileup_summary, "edep_analysis", "events_seen")

    if any(v is None for v in (input_corr, n_stops, pileup_seen)):
        return None
    if sim_total in (None, 0) or pileup_sim in (None, 0):
        return None

    stopping_factor = n_stops / sim_total
    return (pileup_seen / pileup_sim) * input_corr * stopping_factor * _MUON_STOP_PRESCALE_CORRECTION


def _compute_mustop_sample_absolute_efficiencies(
    mubeam_summary: dict, mustop_summary: dict,
) -> dict[str, float | None]:
    """Per-mode absolute efficiency for mustop samples (Edep > 50 MeV)."""
    input_corr = _safe_get(mubeam_summary, "input_efficiency", "correction_factor")
    sim_total = _safe_get(mubeam_summary, "simulation_events", "total_events")
    n_stops = mubeam_summary.get("muminus_stops_events")
    modes = _mustop_modes_from_summary(mustop_summary)
    mustop_sim = _safe_get(mustop_summary, "simulation_events", "total_events")
    sim_per_mode = mustop_sim / len(modes) if mustop_sim not in (None, 0) else None

    if (
        input_corr is None or sim_total in (None, 0)
        or n_stops is None or sim_per_mode in (None, 0)
    ):
        return {m: None for m in modes}

    stopping_factor = n_stops / sim_total
    scale = (
        input_corr * stopping_factor * _GEN_RESTRICTION_FACTOR
        * _MUON_STOP_PRESCALE_CORRECTION / sim_per_mode
    )

    result: dict[str, float | None] = {}
    for mode in modes:
        gt50 = _safe_get(
            mustop_summary, "edep_analysis_by_sample", mode, "events_edep_gt_50_mev"
        )
        result[mode] = gt50 * scale if gt50 is not None else None
    return result


def _parse_events_from_command(command: object) -> int | None:
    """Extract the ``-n <events>`` value from a mu2e command list."""
    if not isinstance(command, list):
        return None
    for i, tok in enumerate(command):
        if tok == "-n" and i + 1 < len(command):
            try:
                return int(command[i + 1])
            except ValueError:
                return None
    return None


def _infer_mustop_mode_from_command(command: object) -> str | None:
    """Infer the mustop mode from the ``-c <fcl>`` argument in *command*."""
    if not isinstance(command, list):
        return None
    for i, tok in enumerate(command):
        if tok == "-c" and i + 1 < len(command):
            fcl_name = Path(command[i + 1]).name
            for mode in _ALL_MUSTOP_MODES:
                if fcl_name.startswith(f"{mode}_job_"):
                    return mode
    return None


def _compute_simulated_events_by_mode(summary: dict | None) -> dict[str, int]:
    """Sum events per mustop mode from completed jobs in a run summary."""
    totals: dict[str, int] = {}
    if not summary:
        return totals
    for row in summary.get("jobs", []):
        if row.get("returncode") != 0:
            continue
        mode = _infer_mustop_mode_from_command(row.get("command"))
        events = _parse_events_from_command(row.get("command"))
        if mode is not None and events is not None:
            totals[mode] = totals.get(mode, 0) + events
    return totals


# ---------------------------------------------------------------------------
# Mu-stops concatenation job
# ---------------------------------------------------------------------------

def _run_mu_stops_job(
    run_dir: Path,
    workflows_dir: Path,
    config_version: str,
    beam_fragment: str,
    mu2e_command: str,
    env: dict,
    dry_run: bool,
) -> JobResult:
    """Concatenate TargetStops into a single MuminusStopsCat file."""
    target_stop_files = _collect_target_stop_files(run_dir)
    if not target_stop_files:
        raise SystemExit(f"No TargetStops files found in {run_dir} for mu_stops job")

    job_dir = run_dir / "mu_stops_job"
    job_dir.mkdir(parents=True, exist_ok=False)

    include_fcl = Path("Run1BAna") / "workflows" / config_version / beam_fragment / "mu_stops.fcl"
    file_list = _format_fhicl_string_list(target_stop_files)
    job_fcl = job_dir / "mu_stops_job.fcl"
    job_fcl.write_text(
        f'#include "{include_fcl.as_posix()}"\n'
        "\n"
        f"source.fileNames: [\n{file_list}\n]\n",
        encoding="utf-8",
    )

    command = [mu2e_command, "-c", str(job_fcl)]
    (job_dir / "job_command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")

    print(f"Running mu_stops job with {len(target_stop_files)} TargetStops input files...")
    result = _run_one_job(0, command, job_dir, env, dry_run)
    _write_job_status(result)

    print(
        f"mu_stops job: returncode={result.returncode}, "
        f"duration={result.duration_s:.2f}s, log={result.log_path}"
    )
    return result


# ---------------------------------------------------------------------------
# Double-edep analysis (ROOT macro)
# ---------------------------------------------------------------------------

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
    """Run the double_edep.C ROOT macro and parse its output."""
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
            "ran": False, "returncode": 0,
            "command": command, "command_path": str(command_path),
            "log_path": str(log_path), "metrics": {}, "error": None,
        }

    try:
        proc = subprocess.run(
            command, cwd=workflows_dir,
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError:
        log_path.write_text("root executable not found on PATH\n", encoding="utf-8")
        return {
            "ran": False, "returncode": None,
            "command": command, "command_path": str(command_path),
            "log_path": str(log_path), "metrics": {},
            "error": "root executable not found on PATH",
        }

    text = f"{proc.stdout}\n{proc.stderr}"
    log_path.write_text(text, encoding="utf-8")

    # Parse metrics from output.
    metrics: dict[str, float | None] = {k: None for k in _DOUBLE_EDEP_PATTERNS}
    metric_lines: dict[str, str] = {}
    threshold_rates: dict[str, float | None] = {
        "estimated_rate_per_event_e50": None,
        "estimated_rate_per_event_e70": None,
        "estimated_rate_per_event_e80": None,
        "estimated_rate_per_event_e90": None,
    }

    for line in text.splitlines():
        stripped = line.strip()
        for key, pattern in _DOUBLE_EDEP_PATTERNS.items():
            m = pattern.search(stripped)
            if m:
                metrics[key] = float(m.group(1))
                metric_lines[key] = stripped
        m = _DOUBLE_EDEP_THRESHOLDS_PATTERN.search(stripped)
        if m:
            threshold_rates = {
                "estimated_rate_per_event_e50": float(m.group(1)),
                "estimated_rate_per_event_e70": float(m.group(2)),
                "estimated_rate_per_event_e80": float(m.group(3)),
                "estimated_rate_per_event_e90": float(m.group(4)),
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


# ---------------------------------------------------------------------------
# Summary printing helpers
# ---------------------------------------------------------------------------

def _fmt(value: object, fmt: str = ".8g", na: str = "unavailable") -> str:
    """Format a numeric value, returning *na* when ``None``."""
    if value is None:
        return na
    return f"{value:{fmt}}"


def _fmt_col(value: object, width: int, fmt: str = ".2e", na: str = "N/A") -> str:
    """Format a value right-aligned inside *width* chars, centred *na* on ``None``."""
    if value is None:
        return na.center(width)
    return f"{value:{width}{fmt}}"


def _compute_sample_metrics(
    sample_stats: dict,
    scale: float | None,
) -> dict:
    """Extract per-sample Edep/tracker metrics and absolute efficiencies."""
    seen = sample_stats.get("events_seen")
    gt50 = sample_stats.get("events_edep_gt_50_mev")
    trk_mpv = sample_stats.get("tracker_front_fit_mpv_mev")
    trk_fwhm = sample_stats.get("tracker_front_fit_fwhm_mev")
    delta_mpv = sample_stats.get("primary_edep_minus_tracker_front_distribution_mpv_mev")
    delta_fwhm = sample_stats.get("primary_edep_minus_tracker_front_distribution_fwhm_mev")

    combined_mpv = (
        trk_mpv + delta_mpv
        if trk_mpv is not None and delta_mpv is not None else None
    )
    combined_fwhm = (
        math.sqrt(trk_fwhm ** 2 + delta_fwhm ** 2)
        if trk_fwhm is not None and delta_fwhm is not None else None
    )

    abs_all = seen * scale if scale is not None and seen is not None else None
    abs_gt50 = gt50 * scale if scale is not None and gt50 is not None else None

    return {
        "events_seen": seen,
        "events_gt50": gt50,
        "tracker_mpv": trk_mpv,
        "tracker_fwhm": trk_fwhm,
        "combined_mpv": combined_mpv,
        "combined_fwhm": combined_fwhm,
        "abs_eff_all": abs_all,
        "abs_eff_gt50": abs_gt50,
    }


def _print_compact_summary_table(
    target_abs: float | None,
    calo_abs: float | None,
    total_edep_per_pot: float | None,
    fgam_abs_50: float | None,
    ce_abs_50: float | None,
    ce_mpv: float | None,
    run1a_ce_abs_50: float | None,
    run1a_ce_mpv: float | None,
    run1a_ce_sens: float | None,
    run1a_total_hit_eff_per_pot: float | None,
) -> None:
    """Print compact summary table with given metrics."""
    sep  = "|--------+----------+-----------+------------+----------+----------+--------+----------+------------+-----------+---------------|"
    hdr1 = "| Config | Mu stop  | Calo stop | Pileup dts | RMC eff  |   CE eff | CE TRK |  Run 1A  |   Run 1A   | Run 1A CE | Run 1A pileup |"
    hdr2 = "|        | per POT  |  per POT  |  <edep>    | Edep(50) | Edep(50) |  MPV   |  CE eff  | CE Trk MPV | S/sqrt(B) |    dts eff    |"

    row = "| config | "
    row += _fmt_col(target_abs, 8) + " | "
    row += _fmt_col(calo_abs, 9) + " | "
    row += _fmt_col(total_edep_per_pot, 10) + " | "
    row += _fmt_col(fgam_abs_50, 8) + " | "
    row += _fmt_col(ce_abs_50, 8) + " | "
    row += _fmt_col(ce_mpv, 6, ".1f") + " | "
    row += _fmt_col(run1a_ce_abs_50, 8) + " | "
    row += _fmt_col(run1a_ce_mpv, 10, ".2f") + " | "
    row += _fmt_col(run1a_ce_sens, 9, ".2f") + " | "
    row += _fmt_col(run1a_total_hit_eff_per_pot, 13) + " | "

    print(sep)
    print(hdr1)
    print(hdr2)
    print(sep)
    print(row)
    print(sep)


# ---------------------------------------------------------------------------
# Beam-stage metric extraction helpers (reduce duplication)
# ---------------------------------------------------------------------------

def _extract_beam_metrics(summary: dict, flash_key: str) -> dict:
    """Pull common beam-stage metrics from a summary dict.

    *flash_key* is the art output type name, e.g. ``"FlashOutput"`` or
    ``"EleFlashOutput"``.
    """
    input_corr = _safe_get(summary, "input_efficiency", "correction_factor")
    sim_total = _safe_get(summary, "simulation_events", "total_events")
    flash_eff = _safe_get(summary, "art_event_analysis", "absolute_efficiency_by_type", flash_key)
    edep_avg = _safe_get(summary, "edep_analysis", "average_calo_energy_mev")
    edep_gt50 = _safe_get(summary, "edep_analysis", "events_edep_gt_50_mev")

    gt50_abs = None
    if input_corr is not None and sim_total not in (None, 0) and edep_gt50 is not None:
        gt50_abs = (edep_gt50 / sim_total) * input_corr

    effective_pot = None
    if sim_total not in (None, 0) and input_corr not in (None, 0):
        effective_pot = sim_total / input_corr

    return {
        "input_corr": input_corr,
        "sim_total": sim_total,
        "flash_eff": flash_eff,
        "edep_avg": edep_avg,
        "edep_gt50": edep_gt50,
        "gt50_abs": gt50_abs,
        "effective_pot": effective_pot,
    }


def _compute_muon_stop_scale(
    input_corr: float | None,
    stopping_factor: float | None,
    sim_per_mode: float | None,
) -> float | None:
    """Compute the per-event absolute-efficiency scale factor for mustop samples."""
    if (
        input_corr is None
        or stopping_factor is None
        or sim_per_mode in (None, 0)
    ):
        return None
    return (
        input_corr * stopping_factor
        * _GEN_RESTRICTION_FACTOR * _MUON_STOP_PRESCALE_CORRECTION
        / sim_per_mode
    )


# ---------------------------------------------------------------------------
# Full summary printer
# ---------------------------------------------------------------------------

def _print_all_stage_compact_summary(
    mubeam_summary: dict,
    elebeam_summary: dict | None,
    mustop_summary: dict | None,
    mustop_pileup_summary: dict | None = None,
    final_summary: dict | None = None,
    run1a_mubeam_summary: dict | None = None,
    run1a_mustops_summary: dict | None = None,
) -> None:
    """Print the full multi-stage summary block."""

    # ---- mubeam ----------------------------------------------------------
    mu = _extract_beam_metrics(mubeam_summary, "FlashOutput")
    target_abs = _safe_get(mubeam_summary, "target_al_analysis", "target_al_entries_absolute_efficiency")
    calo_abs = _safe_get(mubeam_summary, "target_al_analysis", "calo_entries_absolute_efficiency")
    n_muminus = mubeam_summary.get("muminus_stops_events")
    sim_total = mu["sim_total"]
    stopping_factor = (
        n_muminus / sim_total
        if n_muminus is not None and sim_total not in (None, 0) else None
    )

    # ---- elebeam ---------------------------------------------------------
    ele = _extract_beam_metrics(elebeam_summary, "EleFlashOutput") if elebeam_summary else {}

    # ---- mustop ----------------------------------------------------------
    mustop_modes = _mustop_modes_from_summary(mustop_summary)
    mustop_sim_total = (
        _safe_get(mustop_summary, "simulation_events", "total_events")
        if mustop_summary else None
    )
    mustop_sim_per_mode = (
        mustop_sim_total / len(mustop_modes)
        if mustop_sim_total not in (None, 0) and mustop_summary is not None
        else None
    )
    muon_stop_input_eff = (
        mu["input_corr"] * stopping_factor * _MUON_STOP_PRESCALE_CORRECTION
        if mu["input_corr"] is not None and stopping_factor is not None
        else None
    )
    mustop_effective_pot = (
        mustop_sim_per_mode / muon_stop_input_eff
        if mustop_sim_per_mode not in (None, 0) and muon_stop_input_eff not in (None, 0)
        else None
    )

    scale = _compute_muon_stop_scale(mu["input_corr"], stopping_factor, mustop_sim_per_mode)

    # ---- Print beam/mustop lines -----------------------------------------
    print("-----")
    print("Compact all-stage summary")
    print(f"  Target absolute muon stopping rate: {_fmt(target_abs)}")
    print(f"  Calorimeter absolute muon stopping rate: {_fmt(calo_abs)}")
    print(f"  mubeam flash output absolute efficiency: {_fmt(mu.get('flash_eff'))}")
    print(f"  mubeam average Edep / event: {_fmt(mu.get('edep_avg'))} MeV")
    print(f"  mubeam events with Edep > 50 per simulated event (absolute): {_fmt(mu.get('gt50_abs'))}")
    print(f"  mubeam effective N(POT) simulated: {_fmt(mu.get('effective_pot'))}")
    print(f"  elebeam flash output absolute efficiency: {_fmt(ele.get('flash_eff'))}")
    print(f"  elebeam average Edep / event: {_fmt(ele.get('edep_avg'))} MeV")
    print(f"  elebeam events with Edep > 50 per simulated event (absolute): {_fmt(ele.get('gt50_abs'))}")
    print(f"  elebeam effective N(POT) simulated: {_fmt(ele.get('effective_pot'))}")
    print(f"  mustop effective N(POT) simulated (per mode): {_fmt(mustop_effective_pot)}")

    # ---- Per-mode mustop efficiencies ------------------------------------
    fgam_abs_50: float | None = None
    ce_abs_50: float | None = None
    ce_mpv: float | None = None
    mustop_abs_eff_all_by_sample: dict[str, float | None] = {}

    if mustop_summary:
        for sample in mustop_modes:
            sample_stats = mustop_summary.get("edep_analysis_by_sample", {}).get(sample, {})
            sm = _compute_sample_metrics(sample_stats, scale)
            mustop_abs_eff_all_by_sample[sample] = sm["abs_eff_all"]

            print(
                f"  {sample}: abs eff (all)={_fmt(sm['abs_eff_all'])}, "
                f"abs eff (Edep>50)={_fmt(sm['abs_eff_gt50'])}, "
                f"Trk MPV={_fmt(sm['tracker_mpv'], '.4g')} MeV, "
                f"Trk FWHM={_fmt(sm['tracker_fwhm'], '.4g')} MeV"
            )

            if sample == "ce":
                ce_abs_50 = sm["abs_eff_gt50"]
                ce_mpv = sm["tracker_mpv"]
            elif sample == "flat_gamma":
                fgam_abs_50 = sm["abs_eff_gt50"]

    # ---- Run 1A ----------------------------------------------------------
    run1a_target_abs = (
        _safe_get(run1a_mubeam_summary, "target_al_analysis", "target_al_entries_absolute_efficiency")
        if run1a_mubeam_summary else None
    )
    print(f"  run1a target absolute muon stopping rate: {_fmt(run1a_target_abs)}")

    run1a_input_corr = _safe_get(run1a_mubeam_summary, "input_efficiency", "correction_factor") if run1a_mubeam_summary else None
    run1a_flash_eff = _safe_get(run1a_mubeam_summary, "art_event_analysis", "absolute_efficiency_by_type", "FlashOutput") if run1a_mubeam_summary else None
    run1a_sim_total = _safe_get(run1a_mubeam_summary, "simulation_events", "total_events") if run1a_mubeam_summary else None
    run1a_n_stops = run1a_mubeam_summary.get("muminus_stops_events") if run1a_mubeam_summary else None
    run1a_stopping_factor = (
        run1a_n_stops / run1a_sim_total
        if run1a_n_stops is not None and run1a_sim_total not in (None, 0) else None
    )

    run1a_mustops_modes = _mustop_modes_from_summary(run1a_mustops_summary)
    run1a_mustops_sim_total = (
        _safe_get(run1a_mustops_summary, "simulation_events", "total_events")
        if run1a_mustops_summary else None
    )
    run1a_mustops_sim_per_mode = (
        run1a_mustops_sim_total / len(run1a_mustops_modes)
        if run1a_mustops_sim_total not in (None, 0) else None
    )
    run1a_sim_by_mode = _compute_simulated_events_by_mode(run1a_mustops_summary)
    run1a_ref_sim = run1a_sim_by_mode.get("ce")
    if run1a_ref_sim in (None, 0):
        run1a_ref_sim = run1a_sim_by_mode.get("ce_plus")
    if run1a_ref_sim in (None, 0):
        run1a_ref_sim = run1a_mustops_sim_per_mode

    run1a_muon_stop_eff = (
        run1a_input_corr * run1a_stopping_factor
        if run1a_input_corr is not None and run1a_stopping_factor is not None else None
    )
    run1a_mustops_effective_pot = (
        run1a_ref_sim / run1a_muon_stop_eff
        if run1a_ref_sim not in (None, 0) and run1a_muon_stop_eff not in (None, 0) else None
    )
    print(f"  run1a_mustops effective N(POT) simulated (per mode): {_fmt(run1a_mustops_effective_pot)}")

    run1a_ce_abs_50: float | None = None
    run1a_ce_mpv: float | None = None
    run1a_ce_sens: float | None = None
    run1a_abs_eff_all_by_sample: dict[str, float | None] = {}

    for sample in run1a_mustops_modes:
        sample_stats = (
            run1a_mustops_summary.get("edep_analysis_by_sample", {}).get(sample, {})
            if run1a_mustops_summary else {}
        )
        seen = sample_stats.get("events_seen")
        gt50 = sample_stats.get("events_edep_gt_50_mev")
        trk_mpv = sample_stats.get("tracker_front_fit_mpv_mev")
        trk_fwhm = sample_stats.get("tracker_front_fit_fwhm_mev")
        delta_mpv = sample_stats.get("primary_edep_minus_tracker_front_distribution_mpv_mev")
        delta_fwhm = sample_stats.get("primary_edep_minus_tracker_front_distribution_fwhm_mev")

        combined_mpv = (
            trk_mpv + delta_mpv
            if trk_mpv is not None and delta_mpv is not None else None
        )
        combined_fwhm = (
            math.sqrt(trk_fwhm ** 2 + delta_fwhm ** 2)
            if trk_fwhm is not None and delta_fwhm is not None else None
        )

        gen_corr = _GEN_RESTRICTION_FACTOR if sample == "flat_gamma" else 1.0
        sample_sim = run1a_sim_by_mode.get(sample)

        # Normalise flat_gamma to same baseline as ce/ce_plus.
        if sample == "flat_gamma" and run1a_ref_sim not in (None, 0):
            sample_sim = run1a_ref_sim
        elif sample_sim in (None, 0):
            sample_sim = run1a_ref_sim

        sample_scale = (
            run1a_input_corr * run1a_stopping_factor / sample_sim
            if run1a_input_corr is not None
            and run1a_stopping_factor is not None
            and sample_sim not in (None, 0)
            else None
        )

        eff_all = seen / sample_sim * gen_corr if seen is not None and sample_sim not in (None, 0) else None
        eff_gt50 = gt50 / sample_sim * gen_corr if gt50 is not None and sample_sim not in (None, 0) else None
        abs_all = seen * sample_scale * gen_corr if sample_scale is not None and seen is not None else None
        abs_gt50 = gt50 * sample_scale * gen_corr if sample_scale is not None and gt50 is not None else None
        run1a_abs_eff_all_by_sample[sample] = abs_all

        print(
            f"  run1a {sample}: eff/mu stop (all)={_fmt(eff_all)}, "
            f"eff/mu stop (Edep>50)={_fmt(eff_gt50)}, "
            f"abs eff (all)={_fmt(abs_all)}, "
            f"abs eff (Edep>50)={_fmt(abs_gt50)}, "
            f"MPV={_fmt(combined_mpv, '.4g')} MeV, FWHM={_fmt(combined_fwhm, '.4g')} MeV, "
            f"trk-front fit MPV={_fmt(trk_mpv, '.4g')} MeV, FWHM={_fmt(trk_fwhm, '.4g')} MeV"
        )

        if sample == "ce":
            run1a_ce_abs_50 = abs_all
            run1a_ce_mpv = trk_mpv

    # Run 1A rough sensitivity
    run1a_sens_line = (
        _safe_get(run1a_mustops_summary, "rough_run1a_sensitivity", "summary_line")
        if run1a_mustops_summary else None
    )
    print(f"  run1a rough sensitivity (ce): {run1a_sens_line or 'unavailable'}")
    if run1a_sens_line:
        run1a_ce_sens = float(run1a_sens_line.split("=")[-1])

    run1a_fe_abs = run1a_abs_eff_all_by_sample.get("flat_electron")
    run1a_total_hit = (
        run1a_flash_eff + run1a_fe_abs
        if run1a_flash_eff is not None and run1a_fe_abs is not None else None
    )
    print(f"  run1a Total hit efficiency per POT (mubeam + flat_electron): {_fmt(run1a_total_hit)}")

    # ---- Pileup ----------------------------------------------------------
    pileup_abs_all: float | None = None
    pileup_avg: float | None = None
    total_edep_per_pot: float | None = None

    if mustop_pileup_summary is None:
        print("  pileup: summary unavailable")
    else:
        pu_sim = _safe_get(mustop_pileup_summary, "simulation_events", "total_events")
        pu_stats = mustop_pileup_summary.get("edep_analysis", {})
        pu_seen = pu_stats.get("events_seen")
        pu_gt50 = pu_stats.get("events_edep_gt_50_mev")
        pileup_avg = pu_stats.get("average_calo_energy_mev")

        pu_scale = None
        if mu["input_corr"] is not None and stopping_factor is not None and pu_sim not in (None, 0):
            pu_scale = mu["input_corr"] * stopping_factor * _MUON_STOP_PRESCALE_CORRECTION / pu_sim

        pileup_abs_all = pu_seen * pu_scale if pu_scale is not None and pu_seen is not None else None
        pileup_abs_gt50 = pu_gt50 * pu_scale if pu_scale is not None and pu_gt50 is not None else None

        pu_effective_pot = (
            pu_sim / muon_stop_input_eff
            if pu_sim not in (None, 0) and muon_stop_input_eff not in (None, 0) else None
        )
        print(f"  mustop_pileup effective N(POT) simulated: {_fmt(pu_effective_pot)}")
        print(
            f"  pileup: abs eff (all)={_fmt(pileup_abs_all)}, "
            f"abs eff (Edep>50)={_fmt(pileup_abs_gt50)}, "
            f"avg Edep/event={_fmt(pileup_avg)} MeV, "
            f"effective N(POT) simulated={_fmt(pu_effective_pot)}"
        )

        if mu["flash_eff"] is not None and ele.get("flash_eff") is not None and pileup_abs_all is not None:
            total_hit = mu["flash_eff"] + ele["flash_eff"] + pileup_abs_all
            print(f"  Total hit efficiency per POT (mubeam + elebeam + pileup): {_fmt(total_hit)}")
            if mu["edep_avg"] is not None and ele.get("edep_avg") is not None and pileup_avg is not None:
                total_edep_per_pot = (
                    mu["flash_eff"] * mu["edep_avg"]
                    + ele["flash_eff"] * ele["edep_avg"]
                    + pileup_abs_all * pileup_avg
                )
        else:
            print("  Total hit efficiency per POT (mubeam + elebeam + pileup): unavailable")

    # ---- Final (double-edep) ---------------------------------------------
    if final_summary is None:
        print("  Double-Edep expected/event (single,double,triple): unavailable")
        _print_compact_summary_table(
            target_abs, calo_abs, total_edep_per_pot,
            fgam_abs_50, ce_abs_50, ce_mpv,
            run1a_ce_abs_50, run1a_ce_mpv, run1a_ce_sens, run1a_total_hit,
        )
        return

    metrics = final_summary.get("metrics", {})
    print(
        f"  Double-Edep expected/event (single,double,triple): "
        f"{_fmt(metrics.get('expected_per_event_single_edep'))}, "
        f"{_fmt(metrics.get('expected_per_event_double_edep'))}, "
        f"{_fmt(metrics.get('expected_per_event_triple_edep'))}"
    )
    print(
        f"  Double-Edep estimated/event E(50,70,80,90): "
        f"{_fmt(metrics.get('estimated_rate_per_event_e50'))}, "
        f"{_fmt(metrics.get('estimated_rate_per_event_e70'))}, "
        f"{_fmt(metrics.get('estimated_rate_per_event_e80'))}, "
        f"{_fmt(metrics.get('estimated_rate_per_event_e90'))}"
    )

    rough = final_summary.get("rough_sensitivity_by_sample", {})
    for sample in sorted(rough):
        s = rough[sample]
        print(
            f"  Rough sensitivity {sample}: "
            f"MPV={_fmt(s.get('signal_mpv'))}, FWHM={_fmt(s.get('signal_fwhm'))}, "
            f"signal rate={_fmt(s.get('signal_rate'))}, "
            f"background rate={_fmt(s.get('background_rate'))}, "
            f"s/sqrt(b)={_fmt(s.get('s_over_sqrt_b'))}"
        )

    _print_compact_summary_table(
        target_abs, calo_abs, total_edep_per_pot,
        fgam_abs_50, ce_abs_50, ce_mpv,
        run1a_ce_abs_50, run1a_ce_mpv, run1a_ce_sens, run1a_total_hit,
    )


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", choices=_STAGES, default="mubeam",
                    help="Processing stage to run (default: mubeam)")
    p.add_argument("config_version", help="Configuration folder name (e.g. config_v06)")
    p.add_argument("parallel_jobs", type=int, nargs="?", default=None,
                    help="Number of jobs to launch (required for mubeam/elebeam/mustop_pileup/run1a_mubeam/all)")

    # Event counts
    p.add_argument("--events-per-job", type=int, default=0,
                    help="Events per job passed as '-n' to mu2e (required for mubeam/mustop_pileup/all)")
    p.add_argument("--run1a-mubeam-events-per-job", type=int, default=5000,
                    help="Events per run1a_mubeam job (default: 5000)")
    p.add_argument("--elebeam-events-per-job", type=int, default=0,
                    help="Events per elebeam job (default: uses --events-per-job)")
    p.add_argument("--mustop-events-per-job", type=int, default=10000,
                    help="Events per mustop job (default: 10000)")
    p.add_argument("--mustop-pileup-events-per-job", type=int, default=100000,
                    help="Events per mustop_pileup job (default: 100000)")

    # Execution
    p.add_argument("--mu2e-command", default="mu2e",
                    help="Executable used to run jobs (default: mu2e)")
    p.add_argument("--seed-start", type=int, default=1,
                    help="First seed value; each job gets seed_start + job_index")
    p.add_argument("--run-root", default=None,
                    help="Directory where run outputs are written (default: ../runs from this script)")
    p.add_argument("--max-workers", type=int, default=None,
                    help="Thread workers used to launch jobs (default: parallel_jobs)")
    p.add_argument("--dry-run", action="store_true",
                    help="Print and stage commands without executing mu2e")

    # Run-directory overrides
    for stage_name in ("mubeam", "run1a-mubeam", "elebeam", "mustop",
                       "run1a-mustops", "mustop-pileup", "final"):
        p.add_argument(
            f"--{stage_name}-run-dir", default=None,
            help=f"Directory containing {stage_name} outputs",
        )

    # Mustop mode options
    p.add_argument("--mustop-jobs-per-mode", type=int, default=10,
                    help="Jobs to launch for each mustop mode (default: 10)")
    p.add_argument("--include-ce-plus", action="store_true",
                    help="Include ce_plus in mustop/run1a_mustops modes")
    p.add_argument("--mustop-ce-only", action="store_true",
                    help="Run only ce mode for mustop/run1a_mustops")

    # Cleanup
    p.add_argument("--clean", action="store_true",
                    help="Delete intermediate sim.*.art and dts.*.art files (only for --stage all)")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

def _validate_args(args: argparse.Namespace) -> None:
    """Validate argument combinations; raise ``SystemExit`` on errors."""
    # Default elebeam events to general events if unset.
    if args.elebeam_events_per_job <= 0:
        args.elebeam_events_per_job = args.events_per_job

    needs_parallel = {"mubeam", "elebeam", "mustop_pileup", "run1a_mubeam", "all", "minimal"}
    if args.stage in needs_parallel and (args.parallel_jobs is None or args.parallel_jobs <= 0):
        raise SystemExit(f"parallel_jobs must be > 0 for stage {args.stage}")

    if args.stage == "summary":
        return

    no_general_events = {"mustop", "mustop_pileup", "final", "run1a_mubeam", "run1a_mustops", "minimal"}
    if args.stage not in no_general_events and args.events_per_job <= 0:
        raise SystemExit("events_per_job must be > 0")

    if args.stage == "elebeam" and args.elebeam_events_per_job <= 0:
        raise SystemExit("elebeam_events_per_job (or events_per_job) must be > 0")
    if args.stage in ("run1a_mubeam", "minimal") and args.run1a_mubeam_events_per_job <= 0:
        raise SystemExit("run1a_mubeam_events_per_job must be > 0")
    if args.stage in ("mustop", "run1a_mustops", "all") and args.mustop_events_per_job <= 0:
        raise SystemExit("mustop_events_per_job must be > 0")
    if args.stage in ("mustop_pileup", "all") and args.mustop_pileup_events_per_job <= 0:
        raise SystemExit("mustop_pileup_events_per_job must be > 0")
    if args.seed_start <= 0:
        raise SystemExit("seed_start must be > 0")
    if args.stage in ("mustop", "mustop_pileup", "run1a_mustops", "all") and args.mustop_jobs_per_mode <= 0:
        raise SystemExit("mustop_jobs_per_mode must be > 0")


# ---------------------------------------------------------------------------
# Job-spec builders
# ---------------------------------------------------------------------------

def _build_beam_job_specs(
    args: argparse.Namespace,
    workflows_dir: Path,
    extractor_path: Path,
) -> list[dict]:
    """Build job specs for beam stages (mubeam / run1a_mubeam / elebeam)."""
    beam_frag, fcl_name, job_fcl_name = _BEAM_STAGE_CONFIG[args.stage]
    fcl_path = workflows_dir / args.config_version / beam_frag / fcl_name
    include_fcl = Path("Run1BAna") / "workflows" / args.config_version / beam_frag / fcl_name

    if not fcl_path.exists():
        raise SystemExit(f"Missing FCL file: {fcl_path}")
    if not extractor_path.exists():
        raise SystemExit(f"Missing extractor script: {extractor_path}")

    return [
        {
            "index": i,
            "name": args.stage,
            "job_fcl_name": job_fcl_name,
            "include_fcl_path": include_fcl,
            "fcl_overrides": "",
        }
        for i in range(args.parallel_jobs)
    ]


def _build_mustop_job_specs(
    args: argparse.Namespace,
    workflows_dir: Path,
    run_root: Path,
) -> list[dict]:
    """Build job specs for mustop / run1a_mustops stages."""
    input_stage = "mubeam" if args.stage == "mustop" else "run1a_mubeam"
    explicit_dir = args.mubeam_run_dir if args.stage == "mustop" else args.run1a_mubeam_run_dir
    input_dir = _require_stage_dir(explicit_dir, run_root, args.config_version, input_stage)

    if not input_dir.is_dir():
        raise SystemExit(f"mubeam run directory does not exist: {input_dir}")

    muminus_files = _collect_muminus_stop_files(input_dir)
    if not muminus_files:
        raise SystemExit(
            f"No MuminusStopsCat files found in {input_dir / 'mu_stops_job'}. "
            "Run the corresponding mubeam stage to completion first."
        )
    print(f"Using mubeam inputs from: {input_dir}")
    print(f"MuminusStopsCat input files: {len(muminus_files)}")

    fragment = "run1b_mustop" if args.stage == "mustop" else "run1a_mustop"
    modes = _selected_mustop_modes(args.include_ce_plus, args.mustop_ce_only)
    muminus_lines = _format_fhicl_string_list(muminus_files)
    overrides = (
        f"physics.filters.TargetStopResampler.fileNames: [\n{muminus_lines}\n]\n"
    )

    specs: list[dict] = []
    for mode in modes:
        fcl = workflows_dir / args.config_version / fragment / f"{mode}.fcl"
        include = Path("Run1BAna") / "workflows" / args.config_version / fragment / f"{mode}.fcl"
        if not fcl.exists():
            raise SystemExit(f"Missing FCL file: {fcl}")
        for j in range(args.mustop_jobs_per_mode):
            specs.append({
                "index": len(specs),
                "name": f"{mode}_{j:03d}",
                "job_fcl_name": f"{mode}_job_{j:03d}.fcl",
                "include_fcl_path": include,
                "fcl_overrides": overrides,
            })
    return specs


def _build_pileup_job_specs(
    args: argparse.Namespace,
    workflows_dir: Path,
    run_root: Path,
) -> list[dict]:
    """Build job specs for the mustop_pileup stage."""
    mubeam_dir = _require_stage_dir(args.mubeam_run_dir, run_root, args.config_version, "mubeam")
    if not mubeam_dir.is_dir():
        raise SystemExit(f"mubeam run directory does not exist: {mubeam_dir}")

    muminus_files = _collect_muminus_stop_files(mubeam_dir)
    if not muminus_files:
        raise SystemExit(
            f"No MuminusStopsCat files found in {mubeam_dir / 'mu_stops_job'}. "
            "Run stage mubeam to completion first."
        )
    print(f"Using mubeam inputs from: {mubeam_dir}")
    print(f"MuminusStopsCat input files: {len(muminus_files)}")

    fcl = workflows_dir / args.config_version / "run1b_mustop" / "pileup.fcl"
    include = Path("Run1BAna") / "workflows" / args.config_version / "run1b_mustop" / "pileup.fcl"
    if not fcl.exists():
        raise SystemExit(f"Missing FCL file: {fcl}")

    muminus_lines = _format_fhicl_string_list(muminus_files)
    overrides = f"physics.filters.TargetStopResampler.fileNames: [\n{muminus_lines}\n]\n"

    return [
        {
            "index": i,
            "name": f"pileup_{i:03d}",
            "job_fcl_name": f"pileup_job_{i:03d}.fcl",
            "include_fcl_path": include,
            "fcl_overrides": overrides,
        }
        for i in range(args.parallel_jobs)
    ]


def _events_per_job_for_stage(args: argparse.Namespace) -> int:
    """Return the correct per-job event count for the current stage."""
    return {
        "mustop": args.mustop_events_per_job,
        "run1a_mubeam": args.run1a_mubeam_events_per_job,
        "elebeam": args.elebeam_events_per_job,
        "mustop_pileup": args.mustop_pileup_events_per_job,
        "run1a_mustops": args.mustop_events_per_job,
    }.get(args.stage, args.events_per_job)


# ---------------------------------------------------------------------------
# Parallel job execution
# ---------------------------------------------------------------------------

def _launch_parallel_jobs(
    args: argparse.Namespace,
    job_specs: list[dict],
    run_dir: Path,
    env: dict[str, str],
) -> list[JobResult]:
    """Set up FCL files, launch jobs in parallel, and return results."""
    n_events = _events_per_job_for_stage(args)
    max_workers = args.max_workers or len(job_specs)
    max_workers = max(1, min(max_workers, len(job_specs)))

    results: list[JobResult] = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for spec in job_specs:
            idx = spec["index"]
            job_dir = run_dir / f"job_{idx:03d}"
            job_dir.mkdir(parents=True, exist_ok=False)

            seed = args.seed_start + idx
            job_fcl = job_dir / spec["job_fcl_name"]
            job_fcl.write_text(
                f'#include "{spec["include_fcl_path"].as_posix()}"\n'
                "\n"
                f"services.SeedService.baseSeed : {seed}\n"
                f"source.firstSubRun: {idx}\n"
                f"{spec['fcl_overrides']}",
                encoding="utf-8",
            )

            command = [args.mu2e_command, "-c", str(job_fcl), "-n", str(n_events)]
            (job_dir / "job_command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
            futures.append(executor.submit(_run_one_job, idx, command, job_dir, env, args.dry_run))

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            _write_job_status(result)
            print(
                f"Job {result.index:03d}: returncode={result.returncode}, "
                f"duration={result.duration_s:.2f}s, log={result.log_path}"
            )

    return results


# ---------------------------------------------------------------------------
# Art file cleanup
# ---------------------------------------------------------------------------

def _clean_intermediate_art_files(stage_dirs: list[Path | None]) -> None:
    """Delete sim.*.art and dts.*.art files from the given stage directories."""
    print("\nCleaning up intermediate .art files from this run...")
    deleted = 0
    for d in stage_dirs:
        if d is None:
            continue
        for pattern in ("job_*/sim.*.art", "job_*/dts.*.art"):
            for f in d.glob(pattern):
                try:
                    f.unlink()
                    deleted += 1
                except OSError as e:
                    print(f"Warning: Could not delete {f}: {e}", file=sys.stderr)
    print(f"Deleted {deleted} intermediate .art files")


# ---------------------------------------------------------------------------
# Stage handlers
# ---------------------------------------------------------------------------

def _run_stage_summary(args: argparse.Namespace, run_root: Path) -> int:
    """Handle ``--stage summary``."""
    mubeam_dir = _require_stage_dir(args.mubeam_run_dir, run_root, args.config_version, "mubeam")
    mubeam_summary = _load_summary(mubeam_dir / "analysis_summary.json")

    # Load optional summaries.
    stage_dirs: dict[str, Path | None] = {}
    summaries: dict[str, dict | None] = {}
    for name in ("elebeam", "mustop", "run1a_mubeam", "run1a_mustops", "mustop_pileup", "final"):
        attr = getattr(args, f"{name.replace('-', '_')}_run_dir", None)
        d = _resolve_stage_dir(attr, run_root, args.config_version, name)
        stage_dirs[name] = d
        summaries[name] = _load_summary(d / "analysis_summary.json") if d else None

    print(f"mubeam: {mubeam_dir}")
    for name, d in stage_dirs.items():
        if d is not None:
            print(f"{name}: {d}")

    _print_all_stage_compact_summary(
        mubeam_summary,
        summaries["elebeam"],
        summaries["mustop"],
        summaries["mustop_pileup"],
        summaries["final"],
        summaries["run1a_mubeam"],
        summaries["run1a_mustops"],
    )
    return 0


def _run_stage_final(
    args: argparse.Namespace,
    run_root: Path,
    workflows_dir: Path,
) -> int:
    """Handle ``--stage final``."""
    mubeam_dir = _require_stage_dir(args.mubeam_run_dir, run_root, args.config_version, "mubeam")
    elebeam_dir = _require_stage_dir(args.elebeam_run_dir, run_root, args.config_version, "elebeam")
    pileup_dir = _require_stage_dir(args.mustop_pileup_run_dir, run_root, args.config_version, "mustop_pileup")
    mustop_dir = _require_stage_dir(args.mustop_run_dir, run_root, args.config_version, "mustop")

    mubeam_summary = _load_summary(mubeam_dir / "analysis_summary.json")
    elebeam_summary = _load_summary(elebeam_dir / "analysis_summary.json")
    pileup_summary = _load_summary(pileup_dir / "analysis_summary.json")
    mustop_summary = _load_summary(mustop_dir / "analysis_summary.json")

    # Resolve ROOT input files.
    mubeam_edep = Path(mubeam_summary.get("edep_analysis", {}).get("nts_output_path", ""))
    elebeam_edep = Path(elebeam_summary.get("edep_analysis", {}).get("nts_output_path", ""))
    pileup_edep = Path(pileup_summary.get("edep_analysis", {}).get("nts_output_path", ""))
    for label, path in [("mubeam", mubeam_edep), ("elebeam", elebeam_edep), ("pileup", pileup_edep)]:
        if not path.exists():
            raise SystemExit(f"Missing {label} edep root file: {path}")

    # Absolute efficiencies.
    mubeam_flash = _safe_get(mubeam_summary, "art_event_analysis", "absolute_efficiency_by_type", "FlashOutput")
    elebeam_flash = _safe_get(elebeam_summary, "art_event_analysis", "absolute_efficiency_by_type", "EleFlashOutput")
    if mubeam_flash is None:
        raise SystemExit("Missing mubeam FlashOutput absolute efficiency")
    if elebeam_flash is None:
        raise SystemExit("Missing elebeam EleFlashOutput absolute efficiency")
    pileup_eff = _compute_mustop_pileup_absolute_efficiency(mubeam_summary, pileup_summary)
    if pileup_eff is None:
        raise SystemExit("Could not compute mustop_pileup absolute efficiency")

    # Run directory.
    if args.final_run_dir:
        run_dir = Path(args.final_run_dir).resolve()
        if not run_dir.is_dir():
            raise SystemExit(f"final run directory does not exist: {run_dir}")
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = run_root / args.config_version / f"final_{ts}"
        run_dir.mkdir(parents=True, exist_ok=False)

    final_result = _run_double_edep_analysis(
        run_dir, workflows_dir,
        mubeam_edep, elebeam_edep, pileup_edep,
        mubeam_flash, elebeam_flash, pileup_eff,
        args.dry_run,
    )

    dedep_root = _find_double_edep_output_path(run_dir)
    sample_eff = _compute_mustop_sample_absolute_efficiencies(mubeam_summary, mustop_summary)
    rough_sens = extract_analysis_results.run_rough_sensitivity_analyses(
        run_dir, workflows_dir, mustop_summary, sample_eff, dedep_root, args.dry_run,
    )

    summary = {
        "stage": "final",
        "run_dir": str(run_dir),
        "mubeam_run_dir": str(mubeam_dir),
        "elebeam_run_dir": str(elebeam_dir),
        "mustop_pileup_run_dir": str(pileup_dir),
        "inputs": {
            "mubeam_flash_edep_root": str(mubeam_edep),
            "elebeam_flash_edep_root": str(elebeam_edep),
            "mustop_pileup_edep_root": str(pileup_edep),
            "mubeam_flash_absolute_efficiency": mubeam_flash,
            "elebeam_flash_absolute_efficiency": elebeam_flash,
            "mustop_pileup_absolute_efficiency": pileup_eff,
        },
        "double_edep_analysis": final_result,
        "double_edep_output_path": str(dedep_root) if dedep_root else None,
        "rough_sensitivity_by_sample": rough_sens,
        "metrics": final_result.get("metrics", {}),
    }

    summary_path = run_dir / "analysis_summary.json"
    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    print(f"Final-stage run directory: {run_dir}")
    print(f"Analysis summary: {summary_path}")
    if final_result.get("error"):
        print("Final-stage analysis failed", file=sys.stderr)
        return 1
    return 0


def _build_all_stage_base_cmd(args: argparse.Namespace, run_root: Path) -> list[str]:
    """Build the common command prefix for ``--stage all`` / ``--stage minimal``."""
    cmd = [
        sys.executable, str(Path(__file__).resolve()),
        args.config_version, str(args.parallel_jobs),
        "--events-per-job", str(args.events_per_job),
        "--run1a-mubeam-events-per-job", str(args.run1a_mubeam_events_per_job),
        "--mu2e-command", args.mu2e_command,
        "--seed-start", str(args.seed_start),
        "--run-root", str(run_root),
        "--mustop-jobs-per-mode", str(args.mustop_jobs_per_mode),
        "--mustop-events-per-job", str(args.mustop_events_per_job),
        "--mustop-pileup-events-per-job", str(args.mustop_pileup_events_per_job),
        "--elebeam-events-per-job", str(args.elebeam_events_per_job),
    ]
    if args.include_ce_plus:
        cmd.append("--include-ce-plus")
    if args.max_workers is not None:
        cmd.extend(["--max-workers", str(args.max_workers)])
    if args.dry_run:
        cmd.append("--dry-run")
    if args.clean:
        cmd.append("--clean")
    return cmd


def _run_substage(base_cmd: list[str], extra_args: list[str], label: str) -> int:
    """Run a single sub-stage in the all/minimal pipeline."""
    rc = subprocess.run(base_cmd + extra_args, check=False).returncode
    if rc != 0:
        print(f"{label} stage failed", file=sys.stderr)
    return rc


def _run_stage_all(args: argparse.Namespace, run_root: Path) -> int:
    """Handle ``--stage all``."""
    base = _build_all_stage_base_cmd(args, run_root)
    cv = args.config_version

    print("Running stage sequence: mubeam -> elebeam -> mustop -> run1a_mubeam -> run1a_mustops -> mustop_pileup -> final")

    # mubeam
    rc = _run_substage(base, ["--stage", "mubeam"], "mubeam")
    if rc != 0:
        return rc
    mubeam_dir = _require_stage_dir(None, run_root, cv, "mubeam")

    # elebeam
    rc = _run_substage(base, ["--stage", "elebeam"], "elebeam")
    if rc != 0:
        return rc
    elebeam_dir = _require_stage_dir(None, run_root, cv, "elebeam")

    # mustop
    rc = _run_substage(base, ["--stage", "mustop", "--mubeam-run-dir", str(mubeam_dir)], "mustop")
    if rc != 0:
        return rc

    # run1a_mubeam
    rc = _run_substage(base, ["--stage", "run1a_mubeam"], "run1a_mubeam")
    if rc != 0:
        return rc
    r1a_mubeam_dir = _require_stage_dir(None, run_root, cv, "run1a_mubeam")

    # run1a_mustops
    rc = _run_substage(
        base,
        ["--stage", "run1a_mustops", "--run1a-mubeam-run-dir", str(r1a_mubeam_dir), "--mustop-ce-only"],
        "run1a_mustops",
    )
    if rc != 0:
        return rc

    # mustop_pileup
    rc = _run_substage(base, ["--stage", "mustop_pileup", "--mubeam-run-dir", str(mubeam_dir)], "mustop_pileup")
    if rc != 0:
        return rc
    pileup_dir = _require_stage_dir(None, run_root, cv, "mustop_pileup")

    # final
    rc = _run_substage(
        base,
        ["--stage", "final", "--mubeam-run-dir", str(mubeam_dir),
         "--elebeam-run-dir", str(elebeam_dir), "--mustop-pileup-run-dir", str(pileup_dir)],
        "final",
    )
    if rc != 0:
        return rc

    # Print summary
    r1a_mustops_dir = _find_latest_stage_run(run_root, cv, "run1a_mustops")
    mustop_dir = _find_latest_stage_run(run_root, cv, "mustop")
    final_dir = _find_latest_stage_run(run_root, cv, "final")

    print(f"All-stage sequence complete. mubeam run: {mubeam_dir}")
    print(f"All-stage sequence complete. elebeam run: {elebeam_dir}")
    print(f"All-stage sequence complete. run1a_mubeam run: {r1a_mubeam_dir}")
    if r1a_mustops_dir:
        print(f"All-stage sequence complete. run1a_mustops run: {r1a_mustops_dir}")
    if mustop_dir:
        print(f"All-stage sequence complete. mustop run: {mustop_dir}")
    if pileup_dir:
        print(f"All-stage sequence complete. mustop_pileup run: {pileup_dir}")
    if final_dir:
        print(f"All-stage sequence complete. final run: {final_dir}")

    if mustop_dir:
        _print_all_stage_compact_summary(
            _load_summary(mubeam_dir / "analysis_summary.json"),
            _load_summary(elebeam_dir / "analysis_summary.json"),
            _load_summary(mustop_dir / "analysis_summary.json"),
            _load_summary(pileup_dir / "analysis_summary.json") if pileup_dir else None,
            _load_summary(final_dir / "analysis_summary.json") if final_dir else None,
            _load_summary(r1a_mubeam_dir / "analysis_summary.json"),
            _load_summary(r1a_mustops_dir / "analysis_summary.json") if r1a_mustops_dir else None,
        )

    if args.clean:
        _clean_intermediate_art_files([
            mubeam_dir, elebeam_dir, mustop_dir,
            r1a_mustops_dir, r1a_mubeam_dir, pileup_dir, final_dir,
        ])

    return 0


def _run_stage_minimal(args: argparse.Namespace, run_root: Path) -> int:
    """Handle ``--stage minimal``."""
    base = _build_all_stage_base_cmd(args, run_root)
    cv = args.config_version

    print("Running minimal stage sequence: mubeam -> run1a_mubeam -> run1a_mustops")

    rc = _run_substage(base, ["--stage", "mubeam"], "mubeam")
    if rc != 0:
        return rc
    mubeam_dir = _require_stage_dir(None, run_root, cv, "mubeam")

    rc = _run_substage(base, ["--stage", "run1a_mubeam"], "run1a_mubeam")
    if rc != 0:
        return rc
    r1a_mubeam_dir = _require_stage_dir(None, run_root, cv, "run1a_mubeam")

    rc = _run_substage(
        base,
        ["--stage", "run1a_mustops", "--run1a-mubeam-run-dir", str(r1a_mubeam_dir), "--mustop-ce-only"],
        "run1a_mustops",
    )
    if rc != 0:
        return rc
    r1a_mustops_dir = _require_stage_dir(None, run_root, cv, "run1a_mustops")

    print(f"Minimal-stage sequence complete. mubeam run: {mubeam_dir}")
    print(f"Minimal-stage sequence complete. run1a_mubeam run: {r1a_mubeam_dir}")
    print(f"Minimal-stage sequence complete. run1a_mustops run: {r1a_mustops_dir}")

    # Load and print summary information.
    mu_sum = _load_summary_optional(mubeam_dir / "analysis_summary.json")
    r1a_mu_sum = _load_summary_optional(r1a_mubeam_dir / "analysis_summary.json")
    r1a_ms_sum = _load_summary_optional(r1a_mustops_dir / "analysis_summary.json")

    if mu_sum or r1a_mu_sum or r1a_ms_sum:
        print("\n--- Summary Information ---")

        target_abs = _safe_get(mu_sum, "target_al_analysis", "target_al_entries_absolute_efficiency") if mu_sum else None
        calo_abs = _safe_get(mu_sum, "target_al_analysis", "calo_entries_absolute_efficiency") if mu_sum else None

        if mu_sum:
            print(f"mubeam total events simulated: {_safe_get(mu_sum, 'simulation_events', 'total_events', default='unavailable')}")
            print(f"mubeam input efficiency correction: {_fmt(_safe_get(mu_sum, 'input_efficiency', 'correction_factor'))}")
        else:
            print("mubeam summary: unavailable")

        if r1a_mu_sum:
            print(f"run1a_mubeam total events simulated: {_safe_get(r1a_mu_sum, 'simulation_events', 'total_events', default='unavailable')}")
            print(f"run1a_mubeam input efficiency correction: {_fmt(_safe_get(r1a_mu_sum, 'input_efficiency', 'correction_factor'))}")
        else:
            print("run1a_mubeam summary: unavailable")

        if r1a_ms_sum:
            print(f"run1a_mustops total events simulated: {_safe_get(r1a_ms_sum, 'simulation_events', 'total_events', default='unavailable')}")
            modes = _mustop_modes_from_summary(r1a_ms_sum)
            print(f"run1a_mustops modes: {', '.join(modes)}")
        else:
            print("run1a_mustops summary: unavailable")

        # Extract run1a metrics for table.
        run1a_ce_abs_50: float | None = None
        run1a_ce_mpv: float | None = None
        run1a_ce_sens: float | None = None
        run1a_total_hit: float | None = None

        if r1a_ms_sum:
            sim_by_mode = _compute_simulated_events_by_mode(r1a_ms_sum)
            ref_sim = sim_by_mode.get("ce")
            if ref_sim in (None, 0):
                ref_sim = sim_by_mode.get("ce_plus")

            r1a_ic = _safe_get(r1a_mu_sum, "input_efficiency", "correction_factor") if r1a_mu_sum else None
            r1a_stops = r1a_mu_sum.get("muminus_stops_events") if r1a_mu_sum else None
            r1a_sim = _safe_get(r1a_mu_sum, "simulation_events", "total_events") if r1a_mu_sum else None
            r1a_sf = r1a_stops / r1a_sim if r1a_stops is not None and r1a_sim not in (None, 0) else None

            sc = (
                r1a_ic * r1a_sf / ref_sim
                if r1a_ic is not None and r1a_sf is not None and ref_sim not in (None, 0) else None
            )

            ce = r1a_ms_sum.get("edep_analysis_by_sample", {}).get("ce", {})
            ce_gt50 = ce.get("events_edep_gt_50_mev")
            run1a_ce_abs_50 = ce_gt50 * sc if sc is not None and ce_gt50 is not None else None
            run1a_ce_mpv = ce.get("tracker_front_fit_mpv_mev")

            sens_line = _safe_get(r1a_ms_sum, "rough_run1a_sensitivity", "summary_line")
            if sens_line:
                run1a_ce_sens = float(sens_line.split("=")[-1])

            r1a_flash = _safe_get(r1a_mu_sum, "art_event_analysis", "absolute_efficiency_by_type", "FlashOutput") if r1a_mu_sum else None
            fe = r1a_ms_sum.get("edep_analysis_by_sample", {}).get("flat_electron", {})
            fe_seen = fe.get("events_seen")
            fe_abs = fe_seen * sc if sc is not None and fe_seen is not None else None

            run1a_total_hit = (
                r1a_flash + fe_abs
                if r1a_flash is not None and fe_abs is not None else None
            )

        _print_compact_summary_table(
            target_abs, calo_abs, None, None, None, None,
            run1a_ce_abs_50, run1a_ce_mpv, run1a_ce_sens, run1a_total_hit,
        )

    if args.clean:
        _clean_intermediate_art_files([mubeam_dir, r1a_mubeam_dir, r1a_mustops_dir])

    return 0


def _run_single_stage(
    args: argparse.Namespace,
    run_root: Path,
    workflows_dir: Path,
    extractor_path: Path,
) -> int:
    """Handle a single simulation stage (mubeam/elebeam/mustop/run1a_*/mustop_pileup)."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = run_root / args.config_version / f"{args.stage}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)

    # Build job specs.
    if args.stage in _BEAM_STAGE_CONFIG:
        job_specs = _build_beam_job_specs(args, workflows_dir, extractor_path)
    elif args.stage in ("mustop", "run1a_mustops"):
        job_specs = _build_mustop_job_specs(args, workflows_dir, run_root)
    else:  # mustop_pileup
        job_specs = _build_pileup_job_specs(args, workflows_dir, run_root)

    # Print configuration.
    if args.stage in _BEAM_STAGE_CONFIG:
        frag, fcl_name, _ = _BEAM_STAGE_CONFIG[args.stage]
        print(f"FCL: {workflows_dir / args.config_version / frag / fcl_name}")
    elif args.stage == "mustop_pileup":
        print(f"FCL: {workflows_dir / args.config_version / 'run1b_mustop' / 'pileup.fcl'}")
    else:
        fragment = "run1b_mustop" if args.stage == "mustop" else "run1a_mustop"
        modes = _selected_mustop_modes(args.include_ce_plus, args.mustop_ce_only)
        print("FCLs: " + ", ".join(
            str(workflows_dir / args.config_version / fragment / f"{m}.fcl") for m in modes
        ))

    print(f"Run directory: {run_dir}")
    print(f"Stage: {args.stage}")
    if args.stage in ("mustop", "run1a_mustops"):
        print(f"mustop jobs per mode: {args.mustop_jobs_per_mode}")
    if args.stage == "mustop_pileup":
        print(f"mustop_pileup jobs: {args.parallel_jobs}")
    print(f"Launching {len(job_specs)} jobs")
    print(f"Seed range: {args.seed_start} to {args.seed_start + len(job_specs) - 1}")

    # Ensure credentials.
    env = os.environ.copy()
    print("Running getToken to ensure access credentials are ready...")
    token_rc = subprocess.run(["getToken"], env=env, check=False).returncode
    if token_rc != 0:
        raise SystemExit(f"getToken failed with exit code {token_rc}")

    # Run jobs.
    results = _launch_parallel_jobs(args, job_specs, run_dir, env)
    completed = sum(1 for r in results if r.returncode == 0)
    failed = len(results) - completed
    print(f"Finished jobs: {completed}/{len(results)} successful, {failed} failed")

    # Post-processing: mu_stops concatenation for beam stages.
    if args.stage in ("mubeam", "run1a_mubeam"):
        if failed > 0:
            print("Warning: proceeding with available successful mubeam outputs despite failed jobs")
        beam_frag = "run1b_beam" if args.stage == "mubeam" else "run1a_beam"
        mu_result = _run_mu_stops_job(
            run_dir, workflows_dir, args.config_version,
            beam_frag, args.mu2e_command, env, args.dry_run,
        )
        if mu_result.returncode != 0:
            print("mu_stops job failed", file=sys.stderr)
            return 1

    # Run analysis extractor.
    summary_path = run_dir / "analysis_summary.json"
    extractor_cmd = [
        sys.executable, str(extractor_path),
        "--stage", args.stage,
        "--run-dir", str(run_dir),
        "--output", str(summary_path),
        "--pretty",
    ]
    if args.stage == "run1a_mustops" and args.run1a_mubeam_run_dir:
        extractor_cmd.extend(["--run1a-mubeam-run-dir", str(Path(args.run1a_mubeam_run_dir).resolve())])
    if args.stage in ("mustop", "run1a_mustops") and args.include_ce_plus:
        extractor_cmd.append("--include-ce-plus")

    print("Running analysis extractor...")
    (run_dir / "extractor_command.txt").write_text(shlex.join(extractor_cmd) + "\n", encoding="utf-8")

    with (run_dir / "extractor.log").open("w", encoding="utf-8") as log:
        ext_rc = subprocess.run(extractor_cmd, stdout=log, stderr=subprocess.STDOUT, check=False).returncode

    if ext_rc != 0:
        print("Extractor failed", file=sys.stderr)
        return ext_rc

    print(f"Analysis summary: {summary_path}")
    if failed > 0:
        print("Warning: stage completed with failed jobs; efficiencies are evaluated from successful outputs")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    _validate_args(args)

    script_dir = Path(__file__).resolve().parent
    workflows_dir = script_dir.parent
    extractor_path = script_dir / "extract_analysis_results.py"
    run_root = Path(args.run_root).resolve() if args.run_root else workflows_dir / "runs"

    if args.stage == "summary":
        return _run_stage_summary(args, run_root)
    if args.stage == "final":
        return _run_stage_final(args, run_root, workflows_dir)
    if args.stage == "all":
        return _run_stage_all(args, run_root)
    if args.stage == "minimal":
        return _run_stage_minimal(args, run_root)

    # Single simulation stage.
    return _run_single_stage(args, run_root, workflows_dir, extractor_path)


if __name__ == "__main__":
    raise SystemExit(main())
