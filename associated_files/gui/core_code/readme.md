# ribofootPrinter2 GUI: Quick Start

## What does it do?

This graphical interface helps you run the command-line tools in [ribofootPrinter2](https://github.com/guydoshlab/ribofootPrinter2).

It lets you:

- Choose which ribofootPrinter2 scripts to run.
- Enter each script's input settings in a form.
- Preview the exact Python commands before running them.
- Run several selected scripts in order.
- Save a log and a JSON metadata file describing the run.
- Reload metadata from an earlier run to restore its settings.
- View common output files as graphs.

## Files

The application has two Python files:

- `ribofootprinter2_gui.py` — the main interface used to configure and run analyses.
- `ribofootprinter2_viewer.py` — an optional graphical viewer for completed outputs.

Keep these two files together in the same folder.

## Requirements

You need:

- Python 3.
- Tkinter, which is included with many Python installations.
- A working ribofootPrinter2 environment with its required scientific packages.
- Matplotlib to use the graphical results viewer.

Use the same Python environment that you normally use to run ribofootPrinter2.

## Finding the ribofootPrinter2 code

The GUI can find the analysis scripts in any of these ways:

1. Place the GUI files in the ribofootPrinter2 folder, beside its `code` folder.
2. Select the downloaded ribofootPrinter2 folder in **Code/package location**.
3. Install the ribofootPrinter2 scripts as Python modules in the selected Python environment.

If the GUI is already beside the ribofootPrinter2 code, you can leave **Code/package location** blank.

## Start the GUI

Open a terminal in the folder containing the GUI and run:

```bash
python3 ribofootprinter2_gui.py
```

If your ribofootPrinter2 environment uses another Python executable, launch the GUI with that Python instead.

## Run an analysis

1. Choose a **Default output folder**.
2. Choose a **Metadata / log folder**.
3. Confirm the **Python executable** points to your ribofootPrinter2 environment.
4. If necessary, select the ribofootPrinter2 folder under **Code/package location**.
5. Click a script name in the list on the left.
6. Enter its input files and settings.
7. Select the checkbox beside every script you want to run.
8. Review the **Live command preview**.
9. Click **Run selected**.

Selected scripts run one at a time in the order shown on the left. Progress and messages appear in the **Run log** pane.

Note that the ROCC file output of `builddense.py` is not directly funneled into other scripts (a new run selecting the new ROCC file is needed).

Output roots receive a timestamp so a new run does not normally overwrite an earlier run. Note that timestamps are added at runtime to whatever output root name you create.

## Metadata and logs

Each run creates:

- A `.metadata.json` file containing the selected scripts, inputs, settings, and commands.
- A `.log` file containing output and error messages from the analysis scripts.

Use **Load metadata…** to restore the settings from an earlier run. Loading metadata does not immediately rerun the analysis.

## View results

Select **Open the graphical results viewer after the run finishes** if you want the viewer to open automatically.

You can also start the viewer by itself:

```bash
python3 ribofootprinter2_viewer.py
```

Choose a metadata JSON file when prompted.

You may also pass the metadata path directly:

```bash
python3 ribofootprinter2_viewer.py /path/to/ribofootprinter_run_TIMESTAMP.metadata.json
```

For full `metagene_3D` heatmaps, use ribofootPrinter2's `metagene_3D_plot.py`. It is the preferred viewer for those outputs. In addition, there will be nothing shown for the output of `builddense.py`.

## Help and citation

Click **About** in the lower-left corner of the GUI for links to:

- The ribofootPrinter2 GitHub repository.
- The publication to cite when using ribofootPrinter code.

For detailed descriptions of each analysis and its inputs, see the main [ribofootPrinter2 documentation](https://github.com/guydoshlab/ribofootPrinter2).
