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
*ribofootPrinter* is compatible with bowtie1 aligned SAM files. We use the bowtie indexed version of the MANEv1.4_shortnames.fasta reduced transcriptome file for alignments (EBWT file extension), found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_Preparation/output_files) (unzip before use).

An overview on how to run the different *ribofootPrinter* scripts is listed below. This includes example code which can be used directly on the datasets provided [here](MENDELEY LINK). 

Create a new directory in the ribofootPrinter folder and download the SAM files here:
```unix
mkdir -p ./sam_files/
```

More details on preparation of the MANE transcriptome for alignment (shortnames.FASTA) and ribofootPrinter (longnames.FASTA) can be found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_Preparation/).

More details how to view your aligned reads in IGV can be found [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_IGV).

# Running the *ribofootPrinter* package
All scripts are run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is provided in table format for each script.
Most scripts require a ROCC file as input file. The ROCC file is generated from a SAM file using the *builddense* script.
1. *builddense*

ROCC file dependent scripts:

2. *writegene2*
3. *genelist*
4. *posavg*
5. *metagene*
6. *smorflist*
7. *posstats*

The remaining scripts use SAM files as input files (ROCC file generation not required):

8. *region_size_and_abundance*
9. *metagene_3D*


![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/ribofootprinter.png)


## Generate folders for data files
```unix
mkdir -p ./rocc_files/
mkdir -p ./data_files/ 
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
DATA=/full_path/ribofootPrinter2.0-beta/data_files
```

```unix
mkdir -p ./data_files/writegene2 ./data_files/metagene ./data_files/genelist ./data_files/smorflist ./data_files/posavg ./data_files/posstats ./data_files/region_size_and_abundance ./data_files/metagene_3D
```

## Unzip files
```unix
for file in $FASTA/*.zip; do unzip $file; done
```
```unix
mv *.fasta $FASTA
mv *.ebwt $FASTA  
```






## SCRIPT 1. Run *builddense* to generate ROCC files from SAM files
This script requires a SAM file containing bowtie1 aligned reads located in ./ribofootPrinter2.0-beta/sam_files (see above). 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/builddense.png)

The example code generates 5'-end mapped ROCC files from single-end reads from example datasets (ribosome profiling and 40S-seq):
```unix
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $ROCC/80S -1 25 34 1 > $ROCC/80S_bd_metadata.txt 
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $ROCC/40S -1 20 80 1 > $ROCC/40S_bd_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_builddense.png)







## SCRIPT 2. Run *writegene2* from ROCC files
This Python script takes a ROCC file and determines the abundance of 5'- or 3'-end mapped reads at a certain location within the transcript. The script tolerates gene accession numbers (e.g. ENSG00000075624.17) or gene names (e.g. ACTB) as inputs. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/writegene2.png)

The example code outputs data for ACTB and EIF4G2 from ribosome profiling data:
```unix
python $CODE/writegene2.py $ROCC/80S.rocc ACTB $DATA/writegene2/80S_wg2_ACTB > $DATA/writegene2/80S_wg2_ACTB_metadata.txt
python $CODE/writegene2.py $ROCC/80S.rocc EIF4G2 $DATA/writegene2/80S_wg2_EIF4G2 > $DATA/writegene2/80S_wg2_EIF4G2_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_writegene2.png)







## SCRIPT 3. Run *genelist* from ROCC files
This Python script takes a ROCC file and counts reads that map to the main ORF and UTRs of gene models. It also calculates open reading frames on individual transcripts.
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/genelist.png)

The 5'-end aligned data is shifted to accommodate the P-site of the ribosome for example. We used a shift of 12 for our riboseq in the example code. Note that genes with 5'-UTRs shorter than the shift value are excluded from the analysis.

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/shift.png)

```unix
python $CODE/genelist.py $ROCC/80S.rocc 12 1 $DATA/genelist/80S_gl > $DATA/genelist/80S_gl_metadata.txt 
python $CODE/genelist.py $ROCC/40S.rocc 12 1 $DATA/genelist/40S_gl > $DATA/genelist/40S_gl_metadata.txt 
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_genelist.png)







## SCRIPT 4. Run *posavg* from ROCC files
This Python script takes a ROCC file and averages ribosome profiling data around any sequence feature of interest or computes a pause score for every occurrence of all 61 codons or 20 amino acids.
python. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posavg1.png)

The example code calculates pause scores for all amino acids in the main frame for the CDS:
```unix
python $CODE/posavg.py $ROCC/80S.rocc all_1 1 0 0 30 0 12 1 none $DATA/posavg/80S_pa_all_1_aa_CDS_frame0 > $DATA/posavg/80S_pa_all_1_aa_CDS_frame0_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posavg2.png)

The example code outputs data for uORFs with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm0 > $DATA/posavg/80S_pa_uORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm1 > $DATA/posavg/80S_pa_uORF_norm1_metadata.txt 
```

The example code outputs data for dORFs with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm0 > $DATA/posavg/80S_pa_dORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm1 > $DATA/posavg/80S_pa_dORF_norm1_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_posavg.png)






## SCRIPT 5. Run *metagene* from ROCC files
This Python script takes a ROCC file and calculates the average around start or stop codons. Note that if the UTRs are shorter than the 5' or 3' range settings, these transcripts will be excluded from the analysis. Information of number of transcripts included in the analysis can be found as an output.

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/metagene2.png)

How the code works: 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/metagene1.png)


The example code generates a start and stop codon metagene for all transcripts:
```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 none $DATA/metagene/80S_mg_start > $DATA/metagene/80S_mg_start_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 none $DATA/metagene/80S_mg_stop > $DATA/metagene/80S_mg_stop_metadata.txt
```

The example code generates a start and stop codon metagene for a subset of the dataset containing predicted uORF transcripts:
```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 ./associated_files/subset_list.xlsx $DATA/metagene/80S_mg_start_subset > $DATA/metagene/80S_mg_start_subset_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 ./associated_files/subset_list.xlsx $DATA/metagene/80S_mg_stop_subset > $DATA/metagene/80S_mg_stop_subset_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_metagene.png)






## SCRIPT 6. Run *smorflist* from ROCC files
This Python script takes a ROCC file and counts reads that map in frame to uORF or dORFs in the 5’ or 3’ UTR, respectively.
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/smorflist.png)

The example code provides a list of predicted uORF containing transcripts with a length of 4 amino acids strictly starting with ATG:
```unix
python $CODE/smorflist.py $ROCC/80S.rocc 4 12 0 0 5 $DATA/smorflist/80S_uorflist > $DATA/smorflist/80S_uorflist_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_smorflist.png)



## SCRIPT 7. Run *posstats* from ROCC files
This Python script takes a ROCC file and computes pause scores for every individual occurrence of a sequence feature of interest. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posstats.png)
The example code calculates pause scores for PPG and AAA tri-amino acid motifs in our ribosome profiling dataset:
```unix
python $CODE/posstats.py $ROCC/80S.rocc PPG 1 0 10 1 12 1 $DATA/posstats/80S_ps_PPG >  $DATA/posstats/80S_ps_PPG_metadata.txt
python $CODE/posstats.py $ROCC/80S.rocc AAA 1 0 10 1 12 1 $DATA/posstats/80S_ps_AAA >  $DATA/posstats/80S_ps_AAA_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_posstats.png)



## SCRIPT 8. Run *region_size_and_abundance* from SAM files
This Python script takes a bowtie aligned SAM file and outputs the number of reads that aligned to each transcript feature (i.e. UTRs, ORF, start and stop codon) for different read lengths. The saved output consists of a csv file that includes separate sections for read length distributions (left; raw counts for each of the 5 transcript regions) and multiple kinds of abundance measurement for each transcript region (right; can be used to normalize the read length distributions by multiplying the raw count values).
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/rsa2.png)

The window parameter can be adjusted. This setting determines the regions around the start and stop codon in which a mapped read needs to overlap:
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/rsa1.png)

The example code outputs data for ribosome profiling data for all transcripts: 
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/region_size_and_abundance/80S_rsa 25 34 4 none > $DATA/region_size_and_abundance/80S_rsa_metadata.txt
```

The example code outputs data for ribosome profiling data for *EIF4G2*: 
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/region_size_and_abundance/80S_rsa_subset 25 34 4 ./associated_files/subset_list_EIF4G2.xlsx > $DATA/region_size_and_abundance/80S_rsa_subset_metadata.txt
```

The example code outputs data for 40S-seq data for all transcripts: 
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $DATA/region_size_and_abundance/40S_rsa 20 80 4 none > $DATA/region_size_and_abundance/40S_rsa_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_rsa.png)





## SCRIPT 9. Run *metagene_3D* from SAM files
This script takes a bowtie aligned SAM file (single-end reads) and averages reads around start or stop codons as a function of read length. The output consists of 2 csv files with 3-D metagene data (5’- and 3’-assigned reads) and 1 csv file with conventional metagene data (both 5’- and 3’-assigned included). 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/3D_metagene.png)

The example code outputs 3D metagene data for ribosome profiling data on start and stop codons:
```unix
python $CODE/3D_metagene.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_start none 25 34 30 75 1 > $DATA/metagene_3D/80S_3dmg_start_metadata.txt
python $CODE/3D_metagene.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_stop none 25 34 75 30 2 > $DATA/metagene_3D/80S_3dmg_stop_metadata.txt
```

The example code outputs 3D metagene data for 40S-seq data on start codons:
```unix
python $CODE/3D_metagene.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $DATA/metagene_3D/40S_3dmg_start none 20 80 50 50 1 > $DATA/metagene_3D/40S_3dmg_start_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_3D_metagene.png)

Plots can be generated from the 3D output csv files using the script below that generates a matplotlib window with controls to view/save the metagene image. The read abundance is shown by a color scheme. 
Ribosome profiling plots:
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/80S_3dmg_start_3Dmeta_end3.csv 
```
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/80S_3dmg_start_3Dmeta_end5.csv 
```
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/80S_3dmg_stop_3Dmeta_end3.csv 
```
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/80S_3dmg_stop_3Dmeta_end5.csv 
```


40S plots:
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/40S_3dmg_start_3Dmeta_end3.csv 
```
```unix
python $CODE/3D_metagene_plot.py $DATA/metagene_3D/40S_3dmg_start_3Dmeta_end5.csv 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_3D_metagene_plot.png)

