#!/usr/bin/env python3
"""Standalone graphical results viewer for ribofootPrinter2 GUI metadata.

Usage:
    python ribofootprinter2_viewer.py [run.metadata.json]

If a metadata path is omitted, a native file chooser is shown.  The viewer is
also launched by ribofootprinter2_gui.py when its post-run viewer switch is on.
It uses matplotlib (already a ribofootPrinter2 dependency) and the Python
standard library; it never modifies analysis output files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
except ImportError as exc:  # A clear startup error is friendlier than a traceback.
    raise SystemExit(
        "The results viewer requires matplotlib. Run it with the same Python "
        "environment used for ribofootPrinter2.\n" + str(exc)
    )


APP_NAME = "ribofootPrinter2 Results Viewer"
SCHEMA = "ribofootprinter2-gui-metadata-v1"
BLUE = "#52699a"
LIGHT_BLUE = "#9cb8e5"
GREEN = "#a8c994"
DARK_GREEN = "#77ad63"
GREY = "#999999"
PINK = "#d99a9e"
RED = "#b95151"
REGION_COLORS = {"UTR5": GREEN, "start": DARK_GREEN, "CDS": GREY, "stop": "#ff8f7d", "UTR3": PINK}


def read_csv(path: Path) -> list[list[str]]:
    """Read a CSV while preserving the irregular tables used by the tools."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [row for row in csv.reader(handle)]


def number(value, default=math.nan):
    """Convert a CSV cell to float without failing an entire visualization."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def finite(values):
    return [v for v in values if isinstance(v, (int, float)) and math.isfinite(v)]


def output_root(script: str, saved: dict) -> Path | None:
    """Get the output-root parameter despite its two upstream spellings."""
    params = saved.get("parameters", {})
    key = "outfile_path" if script == "metagene_3D" else "outfile"
    value = params.get(key)
    return Path(os.path.expandvars(os.path.expanduser(value))) if value else None


def existing(paths):
    """Return unique existing paths in the order supplied."""
    seen, answer = set(), []
    for path in paths:
        path = Path(path)
        if path.is_file() and path not in seen:
            seen.add(path)
            answer.append(path)
    return answer


def discover_outputs(script: str, saved: dict) -> list[Path]:
    """Reproduce each upstream script's output-name rules from metadata."""
    params = saved.get("parameters", {})
    root = output_root(script, saved)
    if script == "builddense" and root:
        return existing([str(root) + ".rocc"])
    if script in {"writegene2", "genelist", "metagene", "smorflist", "region_size_and_abundance"} and root:
        candidates = [str(root) + ".csv"]
        if script == "genelist":
            candidates.append(str(root) + "_frame.csv")
        return existing(candidates)
    if script == "posavg" and root:
        return existing([str(root) + "_avgdata.csv", str(root) + "_score.csv"])
    if script == "posstats" and root:
        # posstats appends each ROCC sample name to the configured root.
        return sorted(root.parent.glob(root.name + "_*.csv"))
    if script == "metagene_3D" and root:
        return existing([
            str(root) + "_1Dmetas.csv", str(root) + "_3Dmeta_end5.csv", str(root) + "_3Dmeta_end3.csv"
        ])
    if script == "metagene_3D_plot":
        value = params.get("csv_in")
        return existing([value]) if value else []
    return []


def transposed_records(path: Path) -> list[dict[str, str]]:
    """Undo tools.transposecsv's column-oriented record representation."""
    rows = read_csv(path)
    if not rows or len(rows[0]) < 2:
        return []
    records = []
    for column in range(1, len(rows[0])):
        record = {"id": rows[0][column]}
        for row in rows[1:]:
            if row and row[0] and column < len(row):
                record[row[0]] = row[column]
        records.append(record)
    return records


def transpose_rows(rows: list[list[str]]) -> list[list[str]]:
    """Transpose an irregular CSV table, padding missing cells with blanks."""
    width = max((len(row) for row in rows), default=0)
    return [[row[column] if column < len(row) else "" for row in rows] for column in range(width)]


def transposed_series(path: Path) -> list[tuple[str, list[float]]]:
    """Read a transposed file whose columns are labeled numeric series."""
    rows = read_csv(path)
    if not rows:
        return []
    series = []
    for column, label in enumerate(rows[0]):
        values = [number(row[column]) for row in rows[1:] if column < len(row)]
        series.append((label or f"series {column + 1}", finite(values)))
    return series


def row_series(path: Path) -> list[tuple[str, list[float]]]:
    """Read a regular file in which every row starts with a series label."""
    answer = []
    for index, row in enumerate(read_csv(path)):
        if row:
            values = finite([number(cell) for cell in row[1:]])
            if values:
                answer.append((row[0] or f"series {index + 1}", values))
    return answer


def style_axis(axis, title, xlabel="", ylabel=""):
    """Apply a clean style similar to the example figures on GitHub."""
    axis.set_title(title, loc="left", fontsize=11, fontweight="bold")
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#e8e8e8", linewidth=0.7, zorder=0)


def message_figure(title: str, message: str, files=()) -> Figure:
    figure = Figure(figsize=(9, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    axis.axis("off")
    file_text = "\n\nOutput files:\n" + "\n".join(f"• {path}" for path in files) if files else ""
    axis.text(0.04, 0.94, title, va="top", fontsize=18, fontweight="bold", color=BLUE)
    axis.text(0.04, 0.84, message + file_text, va="top", fontsize=11, linespacing=1.5, wrap=True)
    return figure


def plot_writegene2(saved, files):
    if not files:
        return message_figure("writegene2", "No output CSV was found.")
    data = transposed_series(files[0])
    figure = Figure(figsize=(10, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    for label, values in data:
        axis.plot(range(len(values)), values, linewidth=0.8, label=label, color=BLUE if len(data) == 1 else None)
    style_axis(axis, "Read occupancy by transcript position", "Transcript position (nt)", "Reads (RPM)")
    if len(data) <= 12:
        axis.legend(frameon=False, fontsize=8)
    return figure


def _find_field(record, prefix, exclude="raw"):
    return next((key for key in record if key.startswith(prefix) and exclude not in key.lower()), None)


def plot_genelist(saved, files):
    if not files:
        return message_figure("genelist", "No output CSV was found.")
    # Unlike writegene2/metagene, genelist keeps ordinary row-oriented CSVs.
    records = normal_records(files[0])
    if not records:
        return message_figure("genelist", "The output CSV is empty or not readable.", files)
    keys = {prefix: _find_field(records[0], prefix) for prefix in ("UTR5_", "CDS_", "UTR3_")}
    cds = [number(r.get(keys["CDS_"], "")) for r in records]
    utr5 = [number(r.get(keys["UTR5_"], "")) for r in records]
    utr3 = [number(r.get(keys["UTR3_"], "")) for r in records]
    triples = [(x, y, z) for x, y, z in zip(cds, utr5, utr3) if all(math.isfinite(v) and v >= 0 for v in (x, y, z))]
    if not triples:
        return message_figure("genelist", "No numeric region counts were found.", files)
    frame_records = normal_records(files[1]) if len(files) > 1 and files[1].name.endswith("_frame.csv") else []
    panels = 3 if frame_records else 2
    figure = Figure(figsize=(11, 5.5), dpi=100, constrained_layout=True)
    scatter = figure.add_subplot(1, panels, 1)
    for idx, (label, color) in enumerate((("5′-UTR", GREEN), ("3′-UTR", PINK)), start=1):
        x = [max(v[0], 1e-4) for v in triples]
        y = [max(v[idx], 1e-4) for v in triples]
        scatter.scatter(x, y, s=9, alpha=0.45, color=color, label=label, edgecolors="none")
    scatter.set_xscale("log"); scatter.set_yscale("log")
    style_axis(scatter, "Regional reads per transcript", "CDS reads (RPKM)", "UTR reads (RPKM)")
    scatter.legend(frameon=False)
    bars = figure.add_subplot(1, panels, 2)
    totals = [sum(v[i] for v in triples) for i in range(3)]
    denom = sum(totals) or 1
    bars.bar(["5′-UTR", "CDS", "3′-UTR"], [v / denom for v in totals], color=[GREEN, GREY, PINK])
    style_axis(bars, "Fraction of mapped signal", "Transcript region", "Fraction of RPKM")
    if frame_records:
        frame_axis = figure.add_subplot(1, 3, 3)
        frame_keys = [_find_field(frame_records[0], f"CDS{i}_") for i in range(3)]
        frame_totals = [sum(number(record.get(key, ""), 0) for record in frame_records) for key in frame_keys]
        frame_denom = sum(frame_totals) or 1
        frame_axis.bar(["0", "+1", "−1"], [v / frame_denom for v in frame_totals], color=[GREY, "#bbbbbb", "#dddddd"])
        style_axis(frame_axis, "Main ORF frame", "Frame", "Fraction of RPKM")
    return figure


def plot_posavg(saved, files):
    avg = next((path for path in files if path.name.endswith("_avgdata.csv")), None)
    score = next((path for path in files if path.name.endswith("_score.csv")), None)
    if not avg:
        return message_figure("posavg", "No average-data CSV was found.", files)
    series = row_series(avg)
    figure = Figure(figsize=(11, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(121 if score else 111)
    half = int(saved.get("parameters", {}).get("bkndwindow", 0) or 0)
    for index, (label, values) in enumerate(series):
        x = list(range(-half, -half + len(values))) if half else list(range(len(values)))
        axis.plot(x, values, linewidth=1.6, color=BLUE if index == 0 else None, label=label)
    axis.axvline(0, color="#222222", linewidth=0.8, linestyle="--")
    style_axis(axis, "Average signal around sequence motif", "Distance from motif (nt)", "Average reads")
    if len(series) <= 12:
        axis.legend(frameon=False, fontsize=8)
    if score:
        rows = read_csv(score)
        labels = rows[1] if len(rows) > 1 else []
        values = [number(v) for v in (rows[2] if len(rows) > 2 else [])]
        pairs = [(label, value) for label, value in zip(labels, values) if math.isfinite(value)]
        score_axis = figure.add_subplot(122)
        if pairs:
            score_axis.bar([p[0] for p in pairs], [p[1] for p in pairs], color=BLUE)
            score_axis.tick_params(axis="x", labelrotation=70, labelsize=7)
        style_axis(score_axis, "Motif pause scores", "Motif", "Pause score")
    return figure


def plot_metagene(saved, files):
    if not files:
        return message_figure("metagene", "No output CSV was found.")
    series = transposed_series(files[0])
    params = saved.get("parameters", {})
    left = int(params.get("range5", 0) or 0)
    figure = Figure(figsize=(10, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    for index, (label, values) in enumerate(series):
        axis.plot(range(-left, -left + len(values)), values, linewidth=1.2, color=BLUE if index == 0 else None, label=label)
    axis.axvline(0, color="#222222", linewidth=0.8, linestyle="--")
    anchor = "start" if str(params.get("kind")) == "1" else "stop"
    style_axis(axis, f"Metagene around {anchor} codon", f"Distance from {anchor} codon (nt)", "Average reads")
    if len(series) <= 12:
        axis.legend(frameon=False, fontsize=8)
    return figure


def normal_records(path: Path):
    rows = read_csv(path)
    if len(rows) < 2:
        return []
    width = len(rows[0])
    return [dict(zip(rows[0], row + [""] * (width - len(row)))) for row in rows[1:]]


def plot_smorflist(saved, files):
    if not files:
        return message_figure("smorflist", "No output CSV was found.")
    records = normal_records(files[0])
    if not records:
        return message_figure("smorflist", "No small ORFs were reported.", files)
    sm_key = next((k for k in records[0] if k.startswith("sm_orf_")), None)
    cds_key = next((k for k in records[0] if k.startswith("CDS_")), None)
    points = [(number(r.get(cds_key)), number(r.get(sm_key))) for r in records]
    points = [(max(x, 1e-4), max(y, 1e-4)) for x, y in points if math.isfinite(x) and math.isfinite(y) and x >= 0 and y >= 0]
    if not points:
        return message_figure("smorflist", "No numeric small-ORF/CDS counts were found.", files)
    figure = Figure(figsize=(8, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    axis.scatter([p[0] for p in points], [p[1] for p in points], s=12, color=GREY, alpha=0.45, edgecolors="none")
    axis.set_xscale("log"); axis.set_yscale("log")
    style_axis(axis, "Small ORF abundance", "CDS reads (RPKM)", "Small ORF reads (RPKM)")
    return figure


def plot_posstats(saved, files):
    groups, labels = [], []
    for path in files:
        records = normal_records(path)
        scores = [number(r.get("pausescore")) for r in records]
        scores = [math.log10(v) for v in scores if math.isfinite(v) and v > 0]
        if scores:
            groups.append(scores)
            labels.append(path.stem.rsplit("_", 1)[-1])
    if not groups:
        return message_figure("posstats", "No positive pause scores were found.", files)
    figure = Figure(figsize=(9, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    # Set tick labels separately for compatibility across matplotlib versions.
    box = axis.boxplot(groups, notch=True, patch_artist=True)
    axis.set_xticks(range(1, len(labels) + 1))
    axis.set_xticklabels(labels)
    for index, patch in enumerate(box["boxes"]):
        patch.set_facecolor(RED if index % 2 else "white")
        patch.set_edgecolor("#555555")
    style_axis(axis, "Pause-score distributions", "Sample", "log10 (pause score)")
    return figure


def plot_region_size(saved, files):
    if not files:
        return message_figure("region_size_and_abundance", "No output CSV was found.")
    # region_size_and_abundance calls tools.transposecsv before it exits.
    rows = transpose_rows(read_csv(files[0]))
    if len(rows) < 7:
        return message_figure("region_size_and_abundance", "The output table is incomplete.", files)
    lengths = [number(v) for v in rows[1][1:-1]]  # Exclude all_lengths from line plots.
    counts = {row[0]: [number(v, 0) for v in row[1:1 + len(lengths)]] for row in rows[2:7]}
    normalizers = {row[0]: [number(v, 0) for v in row[1:7]] for row in rows[9:] if row}
    modes = [("probability_density", "Fraction of reads"), ("rpm_density", "RPM"), ("rpm_length_density", "RPM per nt")]
    figure = Figure(figsize=(11, 5.5), dpi=100, constrained_layout=True)
    for panel, (mode, ylabel) in enumerate(modes, start=1):
        axis = figure.add_subplot(1, 3, panel)
        factors = normalizers.get(mode, [1] * 6)
        for region_index, region in enumerate(("UTR5", "start", "CDS", "stop", "UTR3")):
            factor = factors[region_index] if region_index < len(factors) else 1
            y = [value * factor for value in counts.get(region, [])]
            axis.plot(lengths, y, linewidth=1.7, color=REGION_COLORS[region], label=region.replace("UTR5", "5′-UTR").replace("UTR3", "3′-UTR"))
        style_axis(axis, mode.replace("_", " ").title(), "Read length (nt)", ylabel)
        if panel == 3:
            axis.legend(frameon=False, fontsize=8)
    return figure


def plot_metagene_3d(saved, files):
    """Show a lightweight 1D preview, but direct users to the preferred tool."""
    one_d = next((path for path in files if path.name.endswith("_1Dmetas.csv")), None)
    note = "For full 3D heatmaps and interactive contrast control, metagene_3D_plot.py is the preferred way to view these outputs (not this tool)."
    if not one_d:
        return message_figure("metagene_3D", note, files)
    rows = read_csv(one_d)
    params = saved.get("parameters", {})
    left = int(params.get("window_left", 0) or 0)
    end5 = [number(row[0]) for row in rows[1:] if row]
    end3 = [number(row[1]) for row in rows[1:] if len(row) > 1]
    figure = Figure(figsize=(10, 6), dpi=100, constrained_layout=True)
    axis = figure.add_subplot(111)
    axis.plot(range(-left, -left + len(end5)), end5, color=BLUE, label="5′-assigned reads")
    axis.plot(range(-left, -left + len(end3)), end3, color=PINK, label="3′-assigned reads")
    axis.axvline(0, color="#222222", linewidth=0.8, linestyle="--")
    style_axis(axis, "metagene_3D — 1D preview", "Distance from anchor codon (nt)", "Average reads (RPM)")
    axis.legend(frameon=False)
    axis.text(0.01, 0.99, note, transform=axis.transAxes, va="top", fontsize=9,
              bbox={"facecolor": "#fff4cf", "edgecolor": "#d8b34d", "pad": 6})
    return figure


def plot_3d_plot_input(saved, files):
    return message_figure(
        "metagene_3D_plot",
        "This script is itself the preferred interactive viewer for metagene_3D heatmaps. "
        "Its selected input is listed below.", files
    )


PLOTTERS = {
    "writegene2": plot_writegene2,
    "genelist": plot_genelist,
    "posavg": plot_posavg,
    "metagene": plot_metagene,
    "smorflist": plot_smorflist,
    "posstats": plot_posstats,
    "region_size_and_abundance": plot_region_size,
    "metagene_3D": plot_metagene_3d,
    "metagene_3D_plot": plot_3d_plot_input,
}


class ResultsViewer(tk.Tk):
    """Resizable Tk shell around matplotlib's interactive canvas and toolbar."""

    def __init__(self, metadata_path: str | None = None):
        super().__init__()
        self.title(APP_NAME)
        self.geometry("1250x820")
        self.minsize(850, 560)
        self.metadata_path = None
        self.metadata = None
        self.entries = []
        self.canvas = None
        self.toolbar = None

        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)
        top = ttk.Frame(self, padding=8)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(1, weight=1)
        ttk.Button(top, text="Open metadata…", command=self.choose_metadata).grid(row=0, column=0, padx=(0, 8))
        self.path_label = ttk.Label(top, text="No metadata loaded", foreground="#555555")
        self.path_label.grid(row=0, column=1, sticky="w")

        pane = ttk.Panedwindow(self, orient="horizontal")
        pane.grid(row=1, column=0, padx=8, pady=(0, 8), sticky="nsew")
        left = ttk.LabelFrame(pane, text="Available results", padding=6, width=260)
        left.grid_propagate(False)
        plot_area = ttk.Frame(pane)
        pane.add(left, weight=0); pane.add(plot_area, weight=1)
        left.rowconfigure(0, weight=1); left.columnconfigure(0, weight=1)
        self.listbox = tk.Listbox(left, exportselection=False)
        self.listbox.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(left, orient="vertical", command=self.listbox.yview)
        scrollbar.grid(row=0, column=1, sticky="ns"); self.listbox.configure(yscrollcommand=scrollbar.set)
        self.listbox.bind("<<ListboxSelect>>", self.show_selected)
        plot_area.rowconfigure(0, weight=1); plot_area.columnconfigure(0, weight=1)
        self.plot_host = ttk.Frame(plot_area)
        self.plot_host.grid(row=0, column=0, sticky="nsew")
        self.status = tk.StringVar(value="Open a metadata file to begin.")
        ttk.Label(self, textvariable=self.status, padding=(10, 0, 10, 8)).grid(row=2, column=0, sticky="ew")

        if metadata_path:
            self.after(0, lambda: self.load_metadata(Path(metadata_path)))
        else:
            self.after(100, self.choose_metadata)

    def choose_metadata(self):
        path = filedialog.askopenfilename(
            title="Open ribofootPrinter2 run metadata", initialdir=str(self.metadata_path.parent) if self.metadata_path else str(Path.cwd()),
            filetypes=[("ribofootPrinter2 metadata", "*.json"), ("All files", "*")]
        )
        if path:
            self.load_metadata(Path(path))

    def load_metadata(self, path: Path):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if data.get("schema") != SCHEMA:
                raise ValueError("The file is not ribofootPrinter2 GUI v1 metadata.")
        except (OSError, json.JSONDecodeError, ValueError, TypeError) as exc:
            messagebox.showerror(APP_NAME, f"Could not open metadata:\n{exc}")
            return
        self.metadata_path, self.metadata = path.resolve(), data
        self.path_label.configure(text=str(self.metadata_path))
        self.entries = []
        self.listbox.delete(0, "end")
        for script, saved in data.get("scripts", {}).items():
            if not saved.get("selected"):
                continue
            files = discover_outputs(script, saved)
            self.entries.append((script, saved, files))
            suffix = f" ({len(files)} file{'s' if len(files) != 1 else ''})" if files else " (no output found)"
            self.listbox.insert("end", script + suffix)
        if not self.entries:
            self.entries = [("run summary", {}, [])]
            self.listbox.insert("end", "Run summary")
        self.listbox.selection_set(0)
        self.show_selected()
        self.status.set(f"Run status: {data.get('status', 'unknown')}  •  {len(self.entries)} result group(s)")

    def show_selected(self, _event=None):
        selection = self.listbox.curselection()
        if not selection:
            return
        script, saved, files = self.entries[selection[0]]
        try:
            if script == "run summary":
                figure = message_figure("Run summary", "No selected scripts were recorded in this metadata.")
            elif script == "builddense":
                figure = message_figure(
                    "builddense",
                    "builddense creates a binary ROCC occupancy file rather than a directly plotted table. "
                    "Use the downstream ribofootPrinter2 analyses to visualize its contents.", files
                )
            else:
                figure = PLOTTERS.get(script, lambda _s, f: message_figure(script, "No specialized preview is available.", f))(saved, files)
        except Exception as exc:
            figure = message_figure(script, f"The output could not be plotted:\n{exc}", files)
        self._install_figure(figure)

    def _install_figure(self, figure):
        for child in self.plot_host.winfo_children():
            child.destroy()
        self.canvas = FigureCanvasTkAgg(figure, master=self.plot_host)
        self.canvas.draw()
        widget = self.canvas.get_tk_widget()
        widget.pack(side="top", fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_host, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="View outputs recorded by ribofootPrinter2 GUI metadata.")
    parser.add_argument("metadata", nargs="?", help="Path to a *.metadata.json file; omit to choose interactively.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    app = ResultsViewer(args.metadata)
    app.mainloop()


if __name__ == "__main__":
    main()
