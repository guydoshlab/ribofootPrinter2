# This page is currently under construction

# ribofootPrinter2.0-beta
ribofootPrinter2.0 is a package of Python tools for analyzing data from ribosome profiling and other RNA-based sequencing experiments. This is a new version of the original ribofootPrinter that adds support for version 1.4 of the MANE transcriptome. It updates the functionality of all the previous core functionalities and adds several new capabilities. For guidance on use of the package, please see the associated manuscript on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.07.04.451082


# Download ribofootPrinter from Github
We suggest to download the entire *ribofootPrinter* folder from Github. 

## Download
You can either download ribofootPrinter manually from the [main Github page](https://github.com/guydoshlab/ribofootPrinter2.0-beta) or copy the entire folder as follows:
```unix
git clone https://github.com/guydoshlab/ribofootPrinter2.0-beta/
```

Following download, navigate to the folder:
```unix
cd ./ribofootPrinter2.0-beta/
```


# Requirements for running ribofootPrinter
## Environment
We run the Python toolbox ribofootPrinter locally inside a conda generated virtual environment using the terminal which has following dependencies installed:
python v3.9.23

biopython v1.85

matplotlib v3.9.4

pandas v2.3.2

numpy v2.0.2

openpyxl v3.1.5

## Alignment files
*ribofootPrinter* is compatible with bowtie1 aligned SAM files. We use the MANEv1.4_shortnames.fasta reduced transcriptome file for alignments, found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_Preparation/output_files).

An overview on how to run the different *ribofootPrinter* scripts is listed below. This includes example code which can be used on the datasets provided [here](MENDELEY LINK). 

Create a new directory in the ribofootPrinter folder and download the SAM files here:
```unix
mkdir -p ./sam_files/
```

More details on preparation of the MANE transcriptome for alignment (shortnames.FASTA) and ribofootPrinter (longnames.FASTA) can be found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_Preparation/output_files).

More details how to view your aligned reads in IGV can be found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_IGV).

# Running the *ribofootPrinter* package
All scripts are run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is provided in table format for each script.
Most scripts require a ROCC file as input file. The ROCC file is generated from a SAM file using the *builddense* script.
1. *builddense*

ROCC file dependent scripts:
2. *writegene2*
3. *metagene*
4. *genelist*
5. *smorflist*
6. *posavg*
7. *posstats*

The remaining scripts use SAM files as input files (ROCC file generation not required):

8. *region_size_and_abundance*
9. *metagene_3D*


![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/ribofootprinter.png)


## Generate folders for output files
```unix
mkdir -p ./rocc_files/
mkdir -p ./output_files/ 
```

## Create variables for different paths
Provide full path to create variables. The full path can be obtained through: 
```unix
pwd
```
```unix
CODE=/full_path/ribofootPrinter2.0-beta/code
FASTA=/full_path/ribofootPrinter2.0-beta/preparation/MANE_v1.4_Preparation/output_files
SAM=/full_path/ribofootPrinter2.0-beta/sam_files
ROCC=/full_path/ribofootPrinter2.0-beta/rocc_files
DATA=/full_path/ribofootPrinter2.0-beta/output_files
```
## Unzip files
```unix
for file in $FASTA/*.zip; do unzip $file; done
```


## SCRIPT 1. Run *builddense* to generate ROCC files from SAM files
*builddense* requires a SAM file containing bowtie1 aligned reads located in ./ribofootPrinter2.0-beta/sam_files (see above). 
This example code will generate 5'-end mapped ROCC files from single-end reads from ribosome profiling data. 

```unix
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $ROCC/80S -1 25 34 1
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_builddense.png)

## SCRIPT 2. Run *writegene2* from ROCC files
This Python script takes a ROCC file and determines the abundance of 5'- or 3'-end mapped reads at a certain location within the transcript.



## SCRIPT 3. Run *metagene* from ROCC files




## SCRIPT 4. Run *genelist* from ROCC files





## SCRIPT 5. Run *smorflist* from ROCC files





## SCRIPT 6. Run *posavg* from ROCC files





## SCRIPT 7. Run *posstats* from ROCC files




## SCRIPT 8. Run *region_size_and_abundance* from SAM files




## SCRIPT 9. Run * metagene_3D* from SAM files









