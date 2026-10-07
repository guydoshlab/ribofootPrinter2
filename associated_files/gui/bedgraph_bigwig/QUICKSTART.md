# SAM to IGV Tracks

This Python GUI converts SAM alignment files into tracks that can be viewed in IGV.

It can create:

- 5′-assigned BedGraph files
- 3′-assigned BedGraph files
- Both 5′ and 3′ BedGraph files in one run
- BigWig coverage files

The 5′ and 3′ BedGraphs can have separate shift values. Multiple SAM files can be processed together. Temporary files are deleted automatically, and the commands used for each SAM file are saved with the results.

## Requirements

- Python 3 with Tkinter
- samtools
- bedtools
- deepTools (`bamCoverage`)

These command-line tools can be installed with Conda:

```bash
conda create -n sam-to-igv -c conda-forge -c bioconda python samtools bedtools deeptools
conda activate sam-to-igv
```

## Run the program

Download `igv_track_gui.py`, activate the environment containing the required tools, and run:

```bash
python igv_track_gui.py
```

## Using the GUI

1. Click **Choose** and select one or more SAM files.
2. Choose an output folder.
3. Select 5′ BedGraph, 3′ BedGraph, BigWig, or any combination.
4. Enter separate shift values for the selected BedGraph outputs. Use `0` for no shift.
5. Review the full output paths.
6. Click **Create tracks**.

Output filenames use the basename of each input SAM file. A command script and run log are also saved for every input file.

For more information, see [ribofootPrinter2](https://github.com/guydoshlab/ribofootPrinter2).
