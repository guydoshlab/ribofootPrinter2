#!/usr/bin/env python3
"""A single-file Tkinter GUI for the ribofootPrinter2 command-line tools.

The GUI deliberately uses only Python's standard library.  ribofootPrinter2 and
its scientific dependencies are needed only when an analysis is actually run.
Place this file beside ribofootPrinter2's ``code`` directory (or beside the
individual scripts), point "Code/package location" at a checkout, or install a
package that exposes the scripts as importable modules.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import queue
import re
import shlex
import subprocess
import sys
import threading
import webbrowser
from datetime import datetime
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk


APP_NAME = "ribofootPrinter2 GUI"
APP_VERSION = "1.0"
SCRIPT_DIR = Path(__file__).resolve().parent
RUN_SUFFIX = re.compile(r"_\d{8}_\d{6}$")


# Parameter dictionaries are kept here, rather than inferred by importing the
# upstream scripts, so opening the GUI does not require BioPython/pandas/etc.
# The order of each tuple is the positional command-line argument order.
SCRIPTS = [
    {
        "name": "builddense",
        "summary": "Convert an aligned SAM file into a normalized ROCC file.",
        "params": [
            ("fasta_in", "Longnames FASTA", "file", "", "Annotated longnames transcriptome FASTA."),
            ("sam_in", "Aligned SAM", "file", "", "Bowtie-aligned SAM input file."),
            ("outfile", "Output root", "output", "builddense_output", "ROCC output path without an extension."),
            ("normalize", "Normalization", "text", "-1", "-1 uses all mapped reads for RPM; 1000000 disables normalization; or enter a divisor."),
            ("smallsize", "Minimum read length", "int", "25", "Smallest included read length (inclusive)."),
            ("largesize", "Maximum read length", "int", "34", "Largest included read length (inclusive)."),
            ("endmode", "Read assignment", "choice", "1", "0 = coverage, 1 = 5′ end, -1 = 3′ end.", ["1", "0", "-1"]),
        ],
    },
    {
        "name": "writegene2",
        "summary": "Export positional read values for one or more genes.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("list_of_genes", "Gene name(s)", "text", "ACTB", "Aliases or ENSG names; use the format accepted by ribofootPrinter2."),
            ("outfile", "Output root", "output", "writegene2_output", "CSV output path without an extension."),
        ],
    },
    {
        "name": "genelist",
        "summary": "Count reads in transcript regions and optionally by CDS frame.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("shift", "Position shift", "int", "12", "Nucleotide shift, commonly 12 for a ribosomal P-site."),
            ("doextra", "Frame output", "choice", "1", "1 writes the extra CDS-frame CSV; 0 does not.", ["1", "0"]),
            ("outfile", "Output root", "output", "genelist_output", "CSV output path without an extension."),
        ],
    },
    {
        "name": "posavg",
        "summary": "Average reads around motifs or calculate motif pause scores.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("motif", "Motif", "text", "all_1", "Nucleotide/amino-acid motif, comma list, or all_N."),
            ("kind", "Motif type", "choice", "1", "0 = nucleotide; 1 = amino acid.", ["1", "0"]),
            ("frame", "Frame", "choice", "0", "0, 1, or 2 for one frame; 3 for all frames.", ["0", "1", "2", "3"]),
            ("bkndwindowthresh", "Background threshold", "int", "0", "Minimum RPKM in the background window."),
            ("bkndwindow", "Background half-window", "int", "30", "Bases on either side of the motif."),
            ("ORFnorm", "ORF normalization", "float", "0", "0 disables it; a positive value is the ORF RPKM threshold."),
            ("shift", "Position shift", "int", "12", "Nucleotide shift to the site of interest."),
            ("UTRmode", "Transcript region", "choice", "1", "0 = 5′ UTR; 1 = CDS; 2 = 3′ UTR.", ["1", "0", "2"]),
            ("subsetlist", "Gene subset", "optional_file", "none", "Excel gene list, or 'none' for all genes."),
            ("outfile", "Output root", "output", "posavg_output", "Output path without an extension."),
        ],
    },
    {
        "name": "metagene",
        "summary": "Average ROCC read counts around start or stop codons.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("kind", "Anchor codon", "choice", "1", "1 = start codon; 2 = stop codon.", ["1", "2"]),
            ("weighting", "Gene weighting", "choice", "1", "1 weights genes equally; 0 averages RPM values directly.", ["1", "0"]),
            ("genethresh", "Gene RPKM threshold", "int", "5", "Minimum whole-gene RPKM for inclusion."),
            ("range5", "5′ window", "int", "50", "Window size on the 5′ side."),
            ("range3", "3′ window", "int", "300", "Window size on the 3′ side."),
            ("subsetlist", "Gene subset", "optional_file", "none", "Excel gene list, or 'none' for all genes."),
            ("outfile", "Output root", "output", "metagene_output", "CSV output path without an extension."),
        ],
    },
    {
        "name": "smorflist",
        "summary": "Find and quantify small ORFs in 5′ or 3′ UTRs.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("lengththresh", "Length threshold (aa)", "int", "4", "ORFs must be longer than this; a negative value selects only that length."),
            ("shift", "Position shift", "int", "12", "Nucleotide shift, commonly 12 for a P-site."),
            ("smallest", "ORF choice", "choice", "0", "1 keeps the smallest ORF; 0 keeps the longest.", ["0", "1"]),
            ("mismatches", "Start mismatches", "choice", "0", "Allowed start-codon mismatches (0 or 1).", ["0", "1"]),
            ("UTR", "UTR region", "choice", "5", "5 = 5′ UTR (uORFs); 3 = 3′ UTR (dORFs).", ["5", "3"]),
            ("outfile", "Output root", "output", "smorflist_output", "CSV output path without an extension."),
        ],
    },
    {
        "name": "posstats",
        "summary": "Report pause scores for every occurrence of a motif.",
        "params": [
            ("inputfiles", "ROCC input(s)", "files", "", "One or more ROCC files; multiple paths are comma-separated."),
            ("motif", "Motif", "text", "PPG", "Nucleotide or amino-acid sequence to score."),
            ("kind", "Motif type", "choice", "1", "0 = nucleotide; 1 = amino acid.", ["1", "0"]),
            ("frame", "Frame", "choice", "0", "Frame in which to find motifs.", ["0", "1", "2"]),
            ("genethresh", "Region RPKM threshold", "int", "10", "Minimum RPKM in the region used as denominator."),
            ("pkwindow", "Peak half-window", "int", "1", "Bases on each side of the peak in the numerator."),
            ("shift", "Position shift", "int", "12", "Positive for 5′-assigned reads; negative for 3′-assigned reads."),
            ("UTRmode", "Transcript region", "choice", "1", "0 = 5′ UTR; 1 = CDS; 2 = 3′ UTR.", ["1", "0", "2"]),
            ("outfile", "Output root", "output", "posstats_output", "Output path without an extension."),
        ],
    },
    {
        "name": "region_size_and_abundance",
        "summary": "Summarize read lengths and abundance by transcript region.",
        "params": [
            ("fasta_in", "Longnames FASTA", "file", "", "Annotated longnames transcriptome FASTA."),
            ("sam_in", "Aligned SAM", "file", "", "Bowtie-aligned SAM input file."),
            ("outfile", "Output root", "output", "region_size_and_abundance_output", "CSV output path without an extension."),
            ("smallsize", "Minimum read length", "int", "25", "Smallest included read length (inclusive)."),
            ("largesize", "Maximum read length", "int", "34", "Largest included read length (inclusive)."),
            ("window", "Codon half-window", "int", "4", "Required overlap window around start/stop codons."),
            ("subset_list", "Gene subset", "optional_file", "none", "Excel gene list, or 'none' for all genes."),
        ],
    },
    {
        "name": "metagene_3D",
        "summary": "Average SAM reads around codons as a function of read length.",
        "warning": "For full output viewing, metagene_3D_plot.py is preferred over the general results viewer.",
        "params": [
            ("fasta_in", "Longnames FASTA", "file", "", "Annotated longnames transcriptome FASTA."),
            ("sam_in", "Aligned SAM", "file", "", "Bowtie-aligned SAM input file."),
            ("outfile_path", "Output root", "output", "metagene_3D_output", "Root used for several output files; no extension."),
            ("subset_list", "Gene subset", "optional_file", "none", "Excel gene list, or 'none' for all genes."),
            ("smallsize", "Minimum read length", "int", "25", "Smallest included read length (inclusive)."),
            ("largesize", "Maximum read length", "int", "34", "Largest included read length (inclusive; must be below 101)."),
            ("window_left", "Upstream window", "int", "30", "Window before the selected codon."),
            ("window_right", "Downstream window", "int", "75", "Window after the selected codon."),
            ("metagene", "Anchor codon", "choice", "1", "1 = start codon; 2 = stop codon.", ["1", "2"]),
        ],
    },
    {
        "name": "metagene_3D_plot",
        "summary": "Open an interactive heatmap from a metagene_3D CSV file.",
        "params": [
            ("csv_in", "3D metagene CSV", "file", "", "A *_3Dmeta_end3.csv or *_3Dmeta_end5.csv file."),
        ],
    },
]


def param_parts(param):
    """Return a uniform 6-tuple for parameter definitions with optional choices."""
    return (*param, None) if len(param) == 5 else param


def find_runner(script_name: str, source: str, python: str) -> tuple[list[str] | None, str]:
    """Resolve a script without importing ribofootPrinter2 or its dependencies."""
    roots = []
    if source.strip():
        roots.append(Path(os.path.expandvars(os.path.expanduser(source.strip()))))
    roots.append(SCRIPT_DIR)
    for root in roots:
        for candidate in (root / f"{script_name}.py", root / "code" / f"{script_name}.py"):
            if candidate.is_file():
                return [python, str(candidate.resolve())], str(candidate.resolve())

    # Common package spellings are supported if a distributor packages the tools.
    modules = [
        f"ribofootprinter2.{script_name}",
        f"ribofootprinter2.code.{script_name}",
        f"ribofootPrinter2.{script_name}",
        f"ribofootPrinter2.code.{script_name}",
        script_name,  # Some installers expose the original scripts as top-level modules.
    ]
    for module in modules:
        try:
            if importlib.util.find_spec(module) is not None:
                return [python, "-m", module], module
        except (ImportError, ModuleNotFoundError, AttributeError):
            pass
    return None, "not found"


def shell_join(parts: list[str]) -> str:
    """Render an exact, copyable command preview."""
    return shlex.join([str(part) for part in parts])


def timestamped_output_root(value: str, run_id: str) -> str:
    """Add or replace the GUI's timestamp suffix on an extension-free output root."""
    path = Path(value)
    base_name = RUN_SUFFIX.sub("", path.name)
    return str(path.with_name(f"{base_name}_{run_id}"))


class ScrollableFrame(ttk.Frame):
    """A vertically scrollable frame that tracks horizontal resizing."""

    def __init__(self, master):
        super().__init__(master)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        bar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.body = ttk.Frame(self.canvas, padding=(8, 4, 12, 8))
        self.window = self.canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.canvas.configure(yscrollcommand=bar.set)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        bar.grid(row=0, column=1, sticky="ns")
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.body.bind("<Configure>", self._sync_scroll)
        self.canvas.bind("<Configure>", self._sync_width)

    def _sync_scroll(self, _event=None):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _sync_width(self, event):
        self.canvas.itemconfigure(self.window, width=event.width)


class RibofootPrinterGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"{APP_NAME} {APP_VERSION}")
        self.geometry("1260x850")
        self.minsize(900, 620)
        self.option_add("*tearOff", False)

        self.specs = {spec["name"]: spec for spec in SCRIPTS}
        self.selected = {name: tk.BooleanVar(value=False) for name in self.specs}
        self.values = {}
        self.widgets = {}
        self.script_buttons = {}
        self.current_name = SCRIPTS[0]["name"]
        self.process = None
        self.worker = None
        self.events = queue.Queue()
        self.stop_requested = False

        default_out = SCRIPT_DIR / "ribofootprinter_output"
        self.output_dir = tk.StringVar(value=str(default_out))
        self.record_dir = tk.StringVar(value=str(default_out / "run_records"))
        self.source_dir = tk.StringVar(value="")
        self.python_exe = tk.StringVar(value=sys.executable)
        self.open_results_viewer = tk.BooleanVar(value=False)
        self.status = tk.StringVar(value="Ready")
        self.default_output_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        self._make_variables()
        self._build_ui()
        self._show_script(self.current_name)
        self._update_preview()
        self.after(100, self._drain_events)
        self.protocol("WM_DELETE_WINDOW", self._close)

    def _make_variables(self):
        for spec in SCRIPTS:
            name = spec["name"]
            self.values[name] = {}
            for param in spec["params"]:
                key, _label, kind, default, _help, _choices = param_parts(param)
                if kind == "output":
                    root = str(Path(self.output_dir.get()) / name / default)
                    value = timestamped_output_root(root, self.default_output_stamp)
                else:
                    value = default
                var = tk.StringVar(value=value)
                var.trace_add("write", lambda *_args: self._update_preview())
                self.values[name][key] = var

    def _build_ui(self):
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)

        paths = ttk.LabelFrame(self, text="Run locations", padding=8)
        paths.grid(row=0, column=0, padx=10, pady=(8, 4), sticky="ew")
        paths.columnconfigure(1, weight=1)
        self._path_row(paths, 0, "Default output folder", self.output_dir, self._choose_output_dir)
        self._path_row(paths, 1, "Metadata / log folder", self.record_dir, self._choose_record_dir)
        self._path_row(paths, 2, "Code/package location", self.source_dir, self._choose_source_dir)
        ttk.Label(paths, text="Python executable").grid(row=3, column=0, padx=(0, 8), pady=3, sticky="w")
        py = ttk.Entry(paths, textvariable=self.python_exe)
        py.grid(row=3, column=1, pady=3, sticky="ew")
        ttk.Button(paths, text="Browse…", command=self._choose_python).grid(row=3, column=2, padx=(8, 0), pady=3)
        ttk.Label(paths, text="Leave code location blank for automatic local/package discovery.", foreground="#555555").grid(
            row=4, column=1, sticky="w", pady=(0, 5)
        )
        ttk.Separator(paths, orient="horizontal").grid(row=5, column=0, columnspan=3, pady=(5, 6), sticky="ew")
        ttk.Checkbutton(
            paths, text="Open the graphical results viewer after the run finishes",
            variable=self.open_results_viewer
        ).grid(row=6, column=1, sticky="w", pady=(0, 2))
        for var in (self.source_dir, self.python_exe):
            var.trace_add("write", lambda *_args: self._update_preview())

        main = ttk.Panedwindow(self, orient="horizontal")
        main.grid(row=1, column=0, padx=10, pady=4, sticky="nsew")
        left = ttk.LabelFrame(main, text="Scripts", padding=6)
        right = ttk.Frame(main)
        main.add(left, weight=0)
        main.add(right, weight=1)
        left.configure(width=315)
        left.grid_propagate(False)
        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)

        list_actions = ttk.Frame(left)
        list_actions.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        ttk.Button(list_actions, text="All", width=7, command=lambda: self._set_all(True)).pack(side="left")
        ttk.Button(list_actions, text="None", width=7, command=lambda: self._set_all(False)).pack(side="left", padx=4)
        script_scroll = ScrollableFrame(left)
        script_scroll.grid(row=1, column=0, sticky="nsew")
        script_scroll.body.columnconfigure(1, weight=1)
        # A dedicated style makes the script currently open in the editor clear.
        self.current_button_font = tkfont.nametofont("TkDefaultFont").copy()
        self.current_button_font.configure(weight="bold")
        self.link_font = tkfont.nametofont("TkDefaultFont").copy()
        self.link_font.configure(underline=True)
        ttk.Style(self).configure("CurrentScript.TButton", font=self.current_button_font)
        for row, spec in enumerate(SCRIPTS):
            name = spec["name"]
            cb = ttk.Checkbutton(script_scroll.body, variable=self.selected[name], command=self._update_preview)
            cb.grid(row=row, column=0, padx=(0, 3), pady=2, sticky="w")
            button = ttk.Button(
                script_scroll.body, text=name, width=31,
                command=lambda n=name: self._show_script(n)
            )
            button.grid(row=row, column=1, pady=2, sticky="ew")
            self.script_buttons[name] = button
        ttk.Separator(left, orient="horizontal").grid(row=2, column=0, pady=(8, 7), sticky="ew")
        ttk.Button(left, text="About", command=self._show_about).grid(row=3, column=0, sticky="ew")

        right.rowconfigure(1, weight=3)
        right.rowconfigure(3, weight=1)
        right.columnconfigure(0, weight=1)
        heading = ttk.Frame(right)
        heading.grid(row=0, column=0, sticky="ew")
        heading.columnconfigure(1, weight=1)
        self.edit_selected = ttk.Checkbutton(heading, text="Run this script", variable=self.selected[self.current_name], command=self._update_preview)
        self.edit_selected.grid(row=0, column=0, padx=(6, 10), pady=5)
        self.title_label = ttk.Label(heading, text="", font=("TkDefaultFont", 13, "bold"))
        self.title_label.grid(row=0, column=1, pady=5, sticky="w")
        self.form = ScrollableFrame(right)
        self.form.grid(row=1, column=0, sticky="nsew")

        preview_box = ttk.LabelFrame(right, text="Live command preview", padding=5)
        preview_box.grid(row=2, column=0, pady=(5, 3), sticky="ew")
        preview_box.columnconfigure(0, weight=1)
        self.preview = tk.Text(preview_box, height=5, wrap="word", state="disabled", font="TkFixedFont")
        self.preview.grid(row=0, column=0, sticky="ew")
        preview_bar = ttk.Scrollbar(preview_box, orient="vertical", command=self.preview.yview)
        preview_bar.grid(row=0, column=1, sticky="ns")
        self.preview.configure(yscrollcommand=preview_bar.set)

        log_box = ttk.LabelFrame(right, text="Run log", padding=5)
        log_box.grid(row=3, column=0, pady=(3, 0), sticky="nsew")
        log_box.rowconfigure(0, weight=1)
        log_box.columnconfigure(0, weight=1)
        self.log = tk.Text(log_box, height=9, wrap="word", state="disabled", font="TkFixedFont")
        self.log.grid(row=0, column=0, sticky="nsew")
        log_bar = ttk.Scrollbar(log_box, orient="vertical", command=self.log.yview)
        log_bar.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=log_bar.set)

        actions = ttk.Frame(self, padding=(10, 4, 10, 9))
        actions.grid(row=2, column=0, sticky="ew")
        actions.columnconfigure(0, weight=1)
        ttk.Label(actions, textvariable=self.status).grid(row=0, column=0, sticky="w")
        ttk.Button(actions, text="Load metadata…", command=self._load_metadata).grid(row=0, column=1, padx=3)
        ttk.Button(actions, text="Save metadata…", command=self._save_metadata_dialog).grid(row=0, column=2, padx=3)
        self.stop_button = ttk.Button(actions, text="Stop", command=self._stop, state="disabled")
        self.stop_button.grid(row=0, column=3, padx=3)
        self.run_button = ttk.Button(actions, text="Run selected", command=self._start_run)
        self.run_button.grid(row=0, column=4, padx=(3, 0))

    def _path_row(self, parent, row, label, variable, callback):
        ttk.Label(parent, text=label).grid(row=row, column=0, padx=(0, 8), pady=3, sticky="w")
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, pady=3, sticky="ew")
        ttk.Button(parent, text="Browse…", command=callback).grid(row=row, column=2, padx=(8, 0), pady=3)

    def _show_about(self):
        """Show project and citation links in a compact modal window."""
        popup = tk.Toplevel(self)
        popup.title(f"About {APP_NAME}")
        popup.transient(self)
        popup.resizable(False, False)
        content = ttk.Frame(popup, padding=18)
        content.grid(row=0, column=0, sticky="nsew")

        items = (
            ("For more information:", "https://github.com/guydoshlab/ribofootPrinter2"),
            ("When using ribofootPrinter code, use this citation:", "https://www.biorxiv.org/content/10.1101/2021.07.04.451082"),
        )
        for row, (description, url) in enumerate(items):
            item = ttk.Frame(content)
            item.grid(row=row, column=0, pady=(0, 14 if row == 0 else 4), sticky="w")
            ttk.Label(item, text=description).pack(anchor="w")
            link = ttk.Label(item, text=url, foreground="#1a5fb4", cursor="hand2", font=self.link_font)
            link.pack(anchor="w", pady=(2, 0))
            link.bind("<Button-1>", lambda _event, address=url: webbrowser.open_new_tab(address))

        ttk.Button(content, text="Close", command=popup.destroy).grid(row=len(items), column=0, pady=(8, 0), sticky="e")
        popup.update_idletasks()
        x = self.winfo_rootx() + max(20, (self.winfo_width() - popup.winfo_width()) // 2)
        y = self.winfo_rooty() + max(20, (self.winfo_height() - popup.winfo_height()) // 2)
        popup.geometry(f"+{x}+{y}")
        popup.grab_set()
        popup.focus_set()

    def _show_script(self, name):
        self.current_name = name
        spec = self.specs[name]
        for script_name, button in self.script_buttons.items():
            button.configure(style="CurrentScript.TButton" if script_name == name else "TButton")
        self.title_label.configure(text=name)
        self.edit_selected.configure(variable=self.selected[name])
        for child in self.form.body.winfo_children():
            child.destroy()
        self.widgets[name] = {}
        self.form.body.columnconfigure(1, weight=1)
        ttk.Label(self.form.body, text=spec["summary"], wraplength=500, foreground="#444444").grid(
            row=0, column=0, columnspan=3, padx=3, pady=(3, 10), sticky="w"
        )
        first_parameter_row = 1
        if spec.get("warning"):
            ttk.Label(
                self.form.body, text=spec["warning"], wraplength=500,
                foreground="#8a5a00", font=self.current_button_font
            ).grid(row=1, column=0, columnspan=3, padx=3, pady=(0, 10), sticky="w")
            first_parameter_row = 2
        for index, param in enumerate(spec["params"]):
            row = first_parameter_row + index * 2
            key, label, kind, _default, help_text, choices = param_parts(param)
            ttk.Label(self.form.body, text=label).grid(row=row, column=0, padx=(3, 10), pady=5, sticky="nw")
            var = self.values[name][key]
            if kind == "choice":
                widget = ttk.Combobox(self.form.body, textvariable=var, values=choices, state="readonly", width=16)
            else:
                widget = ttk.Entry(self.form.body, textvariable=var)
            widget.grid(row=row, column=1, pady=5, sticky="ew")
            self.widgets[name][key] = widget
            if kind in {"file", "files", "optional_file", "output"}:
                ttk.Button(
                    self.form.body,
                    text="Browse…",
                    command=lambda n=name, p=param: self._browse_parameter(n, p),
                ).grid(row=row, column=2, padx=(8, 3), pady=5, sticky="n")
            ttk.Label(self.form.body, text=help_text, wraplength=500, foreground="#555555").grid(
                row=row + 1, column=1, columnspan=2, pady=(0, 3), sticky="w"
            )
        self.form._sync_scroll()
        self._update_preview()

    def _browse_parameter(self, script, param):
        key, label, kind, _default, _help, _choices = param_parts(param)
        current = self.values[script][key].get()
        initial = self.output_dir.get() if kind == "output" else (str(Path(current).parent) if current and current != "none" else str(SCRIPT_DIR))
        if kind == "files":
            paths = filedialog.askopenfilenames(title=f"Choose {label}", initialdir=initial)
            if paths:
                self.values[script][key].set(",".join(paths))
        elif kind == "output":
            path = filedialog.asksaveasfilename(title=f"Choose {label}", initialdir=initial, initialfile=Path(current).name)
            if path:
                self.values[script][key].set(path)
        else:
            path = filedialog.askopenfilename(title=f"Choose {label}", initialdir=initial)
            if path:
                self.values[script][key].set(path)

    def _choose_output_dir(self):
        path = filedialog.askdirectory(title="Choose default output folder", initialdir=self.output_dir.get())
        if not path:
            return
        old = Path(self.output_dir.get())
        self.output_dir.set(path)
        # Relocate only output roots that still live beneath the previous default.
        for spec in SCRIPTS:
            for param in spec["params"]:
                key, _label, kind, _default, _help, _choices = param_parts(param)
                if kind == "output":
                    value = Path(self.values[spec["name"]][key].get())
                    try:
                        rel = value.relative_to(old)
                    except ValueError:
                        continue
                    self.values[spec["name"]][key].set(str(Path(path) / rel))

    def _choose_record_dir(self):
        path = filedialog.askdirectory(title="Choose metadata and log folder", initialdir=self.record_dir.get())
        if path:
            self.record_dir.set(path)

    def _choose_source_dir(self):
        path = filedialog.askdirectory(title="Choose ribofootPrinter2 folder", initialdir=self.source_dir.get() or str(SCRIPT_DIR))
        if path:
            self.source_dir.set(path)

    def _choose_python(self):
        path = filedialog.askopenfilename(title="Choose Python executable", initialdir=str(Path(self.python_exe.get()).parent))
        if path:
            self.python_exe.set(path)

    def _set_all(self, state):
        for var in self.selected.values():
            var.set(state)
        self._update_preview()

    def _timestamp_output_roots(self, names, run_id):
        """Give every selected output root the timestamp shared by this run."""
        for name in names:
            for param in self.specs[name]["params"]:
                key, _label, kind, _default, _help, _choices = param_parts(param)
                if kind == "output":
                    value = self.values[name][key].get().strip()
                    if value:
                        self.values[name][key].set(timestamped_output_root(value, run_id))

    def _command_for(self, name):
        runner, source = find_runner(name, self.source_dir.get(), self.python_exe.get().strip() or sys.executable)
        values = [self.values[name][param_parts(p)[0]].get().strip() for p in self.specs[name]["params"]]
        return (runner + values if runner else None), source

    def _update_preview(self):
        if not hasattr(self, "preview"):
            return
        # Always show the pane being edited, followed by any other selected jobs.
        names = [self.current_name] + [
            n for n, var in self.selected.items() if var.get() and n != self.current_name
        ]
        lines = []
        for name in names:
            command, source = self._command_for(name)
            if command:
                lines.append(shell_join(command))
            else:
                values = [self.values[name][param_parts(p)[0]].get().strip() for p in self.specs[name]["params"]]
                lines.append(f"# {name}.py not found — choose a code location\n" + shell_join([self.python_exe.get(), f"<{name}.py>", *values]))
        self.preview.configure(state="normal")
        self.preview.delete("1.0", "end")
        self.preview.insert("1.0", "\n".join(lines))
        self.preview.configure(state="disabled")

    def _validate(self, names):
        errors = []
        commands = []
        for name in names:
            command, source = self._command_for(name)
            if command is None:
                errors.append(f"{name}: script was not found (searched selected location and installed packages)")
            for param in self.specs[name]["params"]:
                key, label, kind, _default, _help, choices = param_parts(param)
                value = self.values[name][key].get().strip()
                if not value:
                    errors.append(f"{name} — {label}: value is required")
                    continue
                if kind == "int":
                    try:
                        int(value)
                    except ValueError:
                        errors.append(f"{name} — {label}: expected an integer")
                elif kind == "float":
                    try:
                        float(value)
                    except ValueError:
                        errors.append(f"{name} — {label}: expected a number")
                elif kind == "choice" and value not in choices:
                    errors.append(f"{name} — {label}: choose one of {', '.join(choices)}")
                elif kind == "file" and not Path(os.path.expanduser(value)).is_file():
                    errors.append(f"{name} — {label}: file does not exist: {value}")
                elif kind == "files":
                    for path in value.split(","):
                        if not Path(os.path.expanduser(path.strip())).is_file():
                            errors.append(f"{name} — {label}: file does not exist: {path.strip()}")
                elif kind == "optional_file" and value.lower() != "none" and not Path(os.path.expanduser(value)).is_file():
                    errors.append(f"{name} — {label}: file does not exist: {value}")
            if command:
                commands.append((name, command, source))
            # Validate relationships that single-field type checks cannot catch.
            if name in {"builddense", "region_size_and_abundance", "metagene_3D"}:
                try:
                    small = int(self.values[name]["smallsize"].get())
                    large = int(self.values[name]["largesize"].get())
                    if small > large:
                        errors.append(f"{name}: minimum read length must not exceed maximum read length")
                    if name == "metagene_3D" and large > 100:
                        errors.append("metagene_3D: maximum read length must be 100 or less")
                except ValueError:
                    pass  # The field-specific messages above already explain these.
        return errors, commands

    def _metadata(self, commands=None, run_id=None, status="configured"):
        command_map = {name: shell_join(command) for name, command, _source in (commands or [])}
        source_map = {name: source for name, _command, source in (commands or [])}
        scripts = {}
        for spec in SCRIPTS:
            name = spec["name"]
            scripts[name] = {
                "selected": self.selected[name].get(),
                "parameters": {param_parts(p)[0]: self.values[name][param_parts(p)[0]].get() for p in spec["params"]},
            }
            if name in command_map:
                scripts[name]["command"] = command_map[name]
                scripts[name]["resolved_source"] = source_map[name]
        return {
            "schema": "ribofootprinter2-gui-metadata-v1",
            "gui_version": APP_VERSION,
            "saved_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "run_id": run_id,
            "status": status,
            "settings": {
                "default_output_folder": self.output_dir.get(),
                "record_folder": self.record_dir.get(),
                "code_or_package_location": self.source_dir.get(),
                "python_executable": self.python_exe.get(),
                "open_results_viewer_after_run": self.open_results_viewer.get(),
            },
            "scripts": scripts,
        }

    def _save_json(self, path, data):
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)

    def _save_metadata_dialog(self):
        path = filedialog.asksaveasfilename(
            title="Save metadata", initialdir=self.record_dir.get(), initialfile="ribofootprinter_settings.metadata.json",
            defaultextension=".json", filetypes=[("JSON metadata", "*.json"), ("All files", "*")]
        )
        if not path:
            return
        try:
            names = [n for n, var in self.selected.items() if var.get()]
            _errors, commands = self._validate(names)
            self._save_json(Path(path), self._metadata(commands))
            self.status.set(f"Metadata saved: {path}")
        except OSError as exc:
            messagebox.showerror(APP_NAME, f"Could not save metadata:\n{exc}")

    def _load_metadata(self):
        path = filedialog.askopenfilename(
            title="Load metadata", initialdir=self.record_dir.get(),
            filetypes=[("JSON metadata", "*.json"), ("All files", "*")]
        )
        if not path:
            return
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
            if data.get("schema") != "ribofootprinter2-gui-metadata-v1":
                raise ValueError("This is not ribofootPrinter2 GUI v1 metadata.")
            settings = data.get("settings", {})
            for key, var in (
                ("default_output_folder", self.output_dir), ("record_folder", self.record_dir),
                ("code_or_package_location", self.source_dir), ("python_executable", self.python_exe),
            ):
                if key in settings:
                    var.set(str(settings[key]))
            if "open_results_viewer_after_run" in settings:
                self.open_results_viewer.set(bool(settings["open_results_viewer_after_run"]))
            for name, saved in data.get("scripts", {}).items():
                if name not in self.specs:
                    continue
                self.selected[name].set(bool(saved.get("selected", False)))
                for key, value in saved.get("parameters", {}).items():
                    if key in self.values[name]:
                        self.values[name][key].set(str(value))
            self._show_script(self.current_name)
            self.status.set(f"Loaded metadata: {path}")
        except (OSError, json.JSONDecodeError, ValueError, TypeError) as exc:
            messagebox.showerror(APP_NAME, f"Could not load metadata:\n{exc}")

    def _start_run(self):
        names = [n for n, var in self.selected.items() if var.get()]
        if not names:
            messagebox.showwarning(APP_NAME, "Select at least one script to run.")
            return
        errors, _commands = self._validate(names)
        if errors:
            messagebox.showerror(APP_NAME, "Please fix these inputs:\n\n" + "\n".join(f"• {e}" for e in errors[:20]))
            return
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        # Update the visible panes before creating commands and metadata so the
        # preview, actual command, metadata, and resulting filenames all agree.
        self._timestamp_output_roots(names, run_id)
        errors, commands = self._validate(names)
        if errors:  # Defensive: timestamping should never invalidate a root.
            messagebox.showerror(APP_NAME, "Could not prepare timestamped outputs:\n\n" + "\n".join(errors[:20]))
            return
        record_dir = Path(os.path.expandvars(os.path.expanduser(self.record_dir.get())))
        metadata_path = record_dir / f"ribofootprinter_run_{run_id}.metadata.json"
        log_path = record_dir / f"ribofootprinter_run_{run_id}.log"
        try:
            record_dir.mkdir(parents=True, exist_ok=True)
            # Create output parent folders before launching upstream scripts.
            for name in names:
                for param in self.specs[name]["params"]:
                    key, _label, kind, _default, _help, _choices = param_parts(param)
                    if kind == "output":
                        Path(os.path.expanduser(self.values[name][key].get())).parent.mkdir(parents=True, exist_ok=True)
            self._save_json(metadata_path, self._metadata(commands, run_id, "running"))
        except OSError as exc:
            messagebox.showerror(APP_NAME, f"Could not prepare output/record folders:\n{exc}")
            return
        self.stop_requested = False
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status.set(f"Running {len(commands)} script(s)…")
        self._append_log(f"Run {run_id}\nMetadata: {metadata_path}\nLog: {log_path}\n")
        self.worker = threading.Thread(
            target=self._run_worker, args=(commands, metadata_path, log_path, run_id), daemon=True
        )
        self.worker.start()

    def _run_worker(self, commands, metadata_path, log_path, run_id):
        results = []
        final_status = "completed"
        try:
            with log_path.open("w", encoding="utf-8", buffering=1) as log_file:
                log_file.write(f"{APP_NAME} {APP_VERSION}\nRun ID: {run_id}\nStarted: {datetime.now().astimezone().isoformat()}\n\n")
                for name, command, _source in commands:
                    if self.stop_requested:
                        final_status = "stopped"
                        break
                    rendered = shell_join(command)
                    banner = f"\n===== {name} =====\n$ {rendered}\n"
                    log_file.write(banner)
                    self.events.put(("log", banner))
                    try:
                        self.process = subprocess.Popen(
                            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            bufsize=1, cwd=str(SCRIPT_DIR)
                        )
                        assert self.process.stdout is not None
                        for line in self.process.stdout:
                            log_file.write(line)
                            self.events.put(("log", line))
                        return_code = self.process.wait()
                    except OSError as exc:
                        return_code = -1
                        line = f"Could not start process: {exc}\n"
                        log_file.write(line)
                        self.events.put(("log", line))
                    finally:
                        self.process = None
                    results.append({"script": name, "return_code": return_code})
                    if self.stop_requested:
                        final_status = "stopped"
                        break
                    if return_code != 0:
                        final_status = "failed"
                        break
                log_file.write(f"\nFinished: {datetime.now().astimezone().isoformat()}\nStatus: {final_status}\n")
            data = json.loads(metadata_path.read_text(encoding="utf-8"))
            data["status"] = final_status
            data["finished_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
            data["results"] = results
            self._save_json(metadata_path, data)
        except Exception as exc:  # Keep unexpected worker failures visible and recorded.
            final_status = "failed"
            self.events.put(("log", f"\nGUI worker error: {exc}\n"))
        self.events.put(("done", final_status, str(metadata_path), str(log_path)))

    def _stop(self):
        self.stop_requested = True
        process = self.process
        if process is not None and process.poll() is None:
            process.terminate()
        self.status.set("Stopping after the current process exits…")
        self.stop_button.configure(state="disabled")

    def _append_log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text)
        self.log.see("end")
        self.log.configure(state="disabled")

    def _drain_events(self):
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "log":
                    self._append_log(event[1])
                elif event[0] == "done":
                    _kind, status, metadata_path, log_path = event
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    self.status.set(f"Run {status}. Metadata: {metadata_path} | Log: {log_path}")
                    if self.open_results_viewer.get() and status != "stopped":
                        self._launch_viewer(metadata_path)
        except queue.Empty:
            pass
        self.after(100, self._drain_events)

    def _launch_viewer(self, metadata_path):
        """Launch the separate viewer and pass it this run's final metadata."""
        viewer = SCRIPT_DIR / "ribofootprinter2_viewer.py"
        if not viewer.is_file():
            self._append_log(f"\nResults viewer not found: {viewer}\n")
            messagebox.showerror(APP_NAME, f"Results viewer not found:\n{viewer}")
            return
        python = self.python_exe.get().strip() or sys.executable
        command = [python, str(viewer), str(metadata_path)]
        try:
            subprocess.Popen(command, cwd=str(SCRIPT_DIR))
            self._append_log("\nOpened results viewer:\n" + shell_join(command) + "\n")
        except OSError as exc:
            self._append_log(f"\nCould not open results viewer: {exc}\n")
            messagebox.showerror(APP_NAME, f"Could not open results viewer:\n{exc}")

    def _close(self):
        if self.worker and self.worker.is_alive():
            if not messagebox.askyesno(APP_NAME, "A run is active. Stop it and close the GUI?"):
                return
            self._stop()
        self.destroy()


def main():
    app = RibofootPrinterGUI()
    app.mainloop()


if __name__ == "__main__":
    main()
