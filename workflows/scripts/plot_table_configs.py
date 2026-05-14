#!/usr/bin/env python3
"""Parse table.org configuration summaries and plot numeric fields by version."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import ROOT


DESCRIPTION_RE = re.compile(r"^config_(v\d+):\s*(.+?)\s*$")
VERSION_RE = re.compile(r"^v(\d+)$")

CLASS_NOMINAL = "Nominal Run 1A"
CLASS_DEGRADER = "Movable target"
CLASS_SOLID = "Plate target"
CLASS_HOLE = "Plate + hole"
CLASS_PLUG = "Plate + hole and offset plug"
CLASS_FILL = "Plate + hole and thin fill"
CLASS_FOILS = "Plate + hole + foils"
CLASS_HELICAL = "Plate + helical plug"
CLASS_UNKNOWN = "Unknown"

CLASS_ORDER = [
    CLASS_NOMINAL,
    CLASS_DEGRADER,
    CLASS_SOLID,
    CLASS_HOLE,
    CLASS_FILL,
    CLASS_PLUG,
    CLASS_FOILS,
    CLASS_HELICAL,
    CLASS_UNKNOWN,
]

CLASS_COLORS = {
    CLASS_NOMINAL: 1,
    CLASS_DEGRADER: 861,
    CLASS_SOLID: 632,
    CLASS_HOLE: 600,
    CLASS_FILL: 797,
    CLASS_PLUG: 418,
    CLASS_FOILS: ROOT.kMagenta - 9,
    CLASS_HELICAL: 880,
    CLASS_UNKNOWN: 920,
}


CLASS_MARKERS = {
    CLASS_NOMINAL: ROOT.kFullStar,
    CLASS_DEGRADER: ROOT.kFullCircle,
    CLASS_SOLID: ROOT.kFullCircle,
    CLASS_HOLE: ROOT.kFullCircle,
    CLASS_FILL: ROOT.kFullCircle,
    CLASS_PLUG: ROOT.kFullCircle,
    CLASS_FOILS: ROOT.kFullCircle,
    CLASS_HELICAL: ROOT.kFullCircle,
    CLASS_UNKNOWN: ROOT.kFullCircle,
}


@dataclass(frozen=True)
class ConfigRecord:
    version: str
    version_number: int
    description: str
    config_class: str
    fields: dict[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "table_path",
        nargs="?",
        default=Path(__file__).resolve().parent.parent / "table.org",
        type=Path,
        help="Path to the org table file (default: ../table.org)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "figures" / "table_config_plots",
        help="Directory for plot output files",
    )
    parser.add_argument(
        "--root-file",
        default="table_config_plots.root",
        help="Output ROOT file name",
    )
    return parser.parse_args()


def sanitize_field_name(name: str) -> str:
    sanitized = re.sub(r"[^0-9A-Za-z]+", "_", name.strip().lower()).strip("_")
    return sanitized or "field"


def combine_header_rows(row_one: list[str], row_two: list[str]) -> list[str]:
    headers: list[str] = []
    used: dict[str, int] = {}
    for top, bottom in zip(row_one, row_two):
        parts = [part.strip() for part in (top, bottom) if part.strip()]
        combined = " ".join(parts)
        sanitized = sanitize_field_name(combined)
        count = used.get(sanitized, 0)
        used[sanitized] = count + 1
        headers.append(sanitized if count == 0 else f"{sanitized}_{count + 1}")
    return headers


def parse_org_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def parse_float(value: str) -> float | None:
    text = value.strip()
    if not text or text.upper() == "N/A":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def classify_description(description: str) -> str:
    lowered = description.lower()
    if "nominal run 1a" in lowered:
        return CLASS_NOMINAL
    if "degrader" in lowered:
        return CLASS_DEGRADER
    if "plug offset" in lowered or " plug " in f" {lowered} ":
        return CLASS_PLUG
    if "fill" in lowered:
        return CLASS_FILL
    if "foils" in lowered:
        return CLASS_FOILS
    if "helical" in lowered:
        return CLASS_HELICAL
    if "hole" in lowered:
        return CLASS_HOLE
    if "plate target" in lowered:
        return CLASS_SOLID

    return CLASS_UNKNOWN


def parse_table_file(path: Path) -> tuple[list[ConfigRecord], list[str]]:
    descriptions: dict[str, str] = {}
    table_lines: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            match = DESCRIPTION_RE.match(raw_line.rstrip("\n"))
            if match:
                descriptions[match.group(1)] = match.group(2).strip()
            if raw_line.lstrip().startswith("|"):
                table_lines.append(raw_line.rstrip("\n"))

    data_lines = [line for line in table_lines if not line.startswith("|-") and "<l>" not in line and "<c>" not in line]
    if len(data_lines) < 3:
        raise ValueError(f"Could not find a complete table in {path}")

    headers = combine_header_rows(parse_org_row(data_lines[0]), parse_org_row(data_lines[1]))
    records: list[ConfigRecord] = []
    warnings: list[str] = []

    for line in data_lines[2:]:
        cells = parse_org_row(line)
        if not cells or not cells[0].startswith("v"):
            continue

        version = cells[0]
        version_match = VERSION_RE.match(version)
        if version_match is None:
            warnings.append(f"Skipping malformed version entry: {version}")
            continue

        description = descriptions.get(version, "")
        config_class = classify_description(description)
        if config_class == CLASS_UNKNOWN:
            warnings.append(f"Warning: {version} has unknown class: {description or '<missing description>'}")

        row_values = cells + [""] * (len(headers) - len(cells))
        fields = dict(zip(headers, row_values))
        records.append(
            ConfigRecord(
                version=version,
                version_number=int(version_match.group(1)),
                description=description,
                config_class=config_class,
                fields=fields,
            )
        )

    records.sort(key=lambda record: record.version_number)
    return records, warnings


def numeric_fields(records: Iterable[ConfigRecord]) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for record in records:
        for name, value in record.fields.items():
            if name == "config" or name in seen:
                continue
            if parse_float(value) is not None:
                seen.add(name)
                names.append(name)
    return names


def import_root():
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise SystemExit("PyROOT is required to run this script") from exc
    return ROOT


def make_graph(ROOT, configs: list[ConfigRecord], field_name: str):
    graph = ROOT.TGraph(len(configs))
    graph.SetName(f"{field_name}_{sanitize_field_name(configs[0].config_class) if configs else 'graph'}")
    for index, config in enumerate(configs):
        graph.SetPoint(index, config.version_number, parse_float(config.fields[field_name]))
    return graph


def make_xy_graph(ROOT, configs: list[ConfigRecord], x_field: str, y_field: str, class_name: str):
    graph = ROOT.TGraph(len(configs))
    graph.SetName(f"{sanitize_field_name(y_field)}_vs_{sanitize_field_name(x_field)}_{sanitize_field_name(class_name)}")
    for index, config in enumerate(configs):
        graph.SetPoint(index, parse_float(config.fields[x_field]), parse_float(config.fields[y_field]))
    return graph


def plot_fields(records: list[ConfigRecord], output_dir: Path, root_file_name: str) -> None:
    ROOT = import_root()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    output_dir.mkdir(parents=True, exist_ok=True)
    root_output_path = output_dir / root_file_name
    root_file = ROOT.TFile(str(root_output_path), "RECREATE")

    field_names = numeric_fields(records)
    grouped_records = {
        class_name: [record for record in records if record.config_class == class_name]
        for class_name in CLASS_ORDER
    }

    for field_name in field_names:
        multigraph = ROOT.TMultiGraph()
        legend = ROOT.TLegend(0.12, 0.80, 0.88, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetNColumns(3)

        canvas = ROOT.TCanvas(f"canvas_{field_name}", field_name, 1200, 700)
        plotted_graphs = []

        for class_name in CLASS_ORDER:
            class_records = [record for record in grouped_records[class_name] if parse_float(record.fields.get(field_name, "")) is not None]
            if not class_records:
                continue

            graph = make_graph(ROOT, class_records, field_name)
            color = CLASS_COLORS[class_name]
            marker = CLASS_MARKERS[class_name]
            graph.SetTitle(class_name)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.2)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetLineWidth(2)
            multigraph.Add(graph, "P")
            legend.AddEntry(graph, class_name, "p")
            plotted_graphs.append(graph)

        if not plotted_graphs:
            canvas.Close()
            continue

        multigraph.SetTitle(f"{field_name} vs version;version;{field_name}")
        multigraph.Draw("A")
        multigraph.GetXaxis().SetLimits(-1.0, max(record.version_number for record in records) + 1.0)
        multigraph.GetXaxis().SetNdivisions(510)
        canvas.SetGrid()
        legend.Draw()
        canvas.Modified()
        canvas.Update()
        canvas.Write()
        canvas.SaveAs(str(output_dir / f"{field_name}.png"))

    scatter_x_field = "run_1a_ce_s_sqrt_b"
    scatter_y_field = "calo_stop_per_pot"
    nominal_x_values = [
        parse_float(record.fields.get(scatter_x_field, ""))
        for record in grouped_records.get(CLASS_NOMINAL, [])
    ]
    x_scale = next((value for value in nominal_x_values if value not in (None, 0.0)), 1.0)
    if scatter_x_field in field_names and scatter_y_field in field_names:
        scatter_canvas = ROOT.TCanvas("canvas_calo_stop_vs_run1a_s_sqrt_b", "calo_stop_vs_run1a_s_sqrt_b", 1000, 800)
        scatter_multigraph = ROOT.TMultiGraph()
        scatter_legend = ROOT.TLegend(0.10, 0.80, 0.88, 0.9)
        scatter_legend.SetBorderSize(0)
        scatter_legend.SetFillStyle(0)
        scatter_legend.SetNColumns(3)
        scatter_legend.SetTextSize(0.028)
        plotted_scatter_graphs = []

        for class_name in CLASS_ORDER:
            class_records = []
            for record in grouped_records[class_name]:
                raw_x_value = parse_float(record.fields.get(scatter_x_field, ""))
                y_value = parse_float(record.fields.get(scatter_y_field, ""))
                if raw_x_value is None or y_value is None:
                    continue
                x_value = raw_x_value / x_scale
                record = ConfigRecord(
                    version=record.version,
                    version_number=record.version_number,
                    description=record.description,
                    config_class=record.config_class,
                    fields={**record.fields, scatter_x_field: str(x_value)},
                )
                class_records.append(record)

            if not class_records:
                continue

            graph = make_xy_graph(ROOT, class_records, scatter_x_field, scatter_y_field, class_name)
            color = CLASS_COLORS[class_name]
            marker = CLASS_MARKERS[class_name]
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.2 if marker != ROOT.kFullStar else 2.0)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            scatter_multigraph.Add(graph, "P")
            scatter_legend.AddEntry(graph, class_name, "p")
            plotted_scatter_graphs.append(graph)

        if plotted_scatter_graphs:
            scatter_multigraph.SetTitle("Calo muon stops per POT vs. Relative Run 1A S/#sqrt{B};Relative Run 1A S/#sqrt{B};Calo muon stops per POT")
            scatter_multigraph.Draw("A")
            scatter_canvas.SetGrid()
            scatter_legend.Draw()
            scatter_canvas.Modified()
            scatter_canvas.Update()
            scatter_canvas.Write()
            scatter_multigraph.GetYaxis().SetRangeUser(1.e-8, 1.e-4)
            scatter_canvas.SetLogy()
            scatter_canvas.SaveAs(str(output_dir / "calo_stop_vs_run1a_s_sqrt_b.png"))
        else:
            scatter_canvas.Close()

    root_file.Close()


def print_summary(records: list[ConfigRecord], warnings: list[str]) -> None:
    print(f"Parsed {len(records)} configurations")
    for warning in warnings:
        print(warning)

    for record in records:
        print(f"{record.version}: {record.config_class}: {record.description}")


def main() -> int:
    args = parse_args()
    records, warnings = parse_table_file(args.table_path.resolve())
    print_summary(records, warnings)
    plot_fields(records, args.output_dir.resolve(), args.root_file)
    print(f"Wrote plots to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
