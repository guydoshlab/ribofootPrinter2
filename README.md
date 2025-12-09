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
FILES=/full_path/ribofootPrinter2.0-beta/associated_files
```

```unix
mkdir -p ./data_files/writegene2 ./data_files/metagene ./data_files/genelist ./data_files/smorflist ./data_files/posavg ./data_files/posstats ./data_files/region_size_and_abundance ./data_files/metagene_3D
```

## Unzip files
```unix
for file in $FASTA/*.zip; do unzip $file; done
for file in $FILES/*.zip; do unzip $file; done
```
```unix
mv *.fasta $FASTA
mv *.ebwt $FASTA  
mv *.xlsx $FILES 
```






## SCRIPT 1. Run *builddense* to generate ROCC files from SAM files
*builddense* converts SAM files into the ROCC files which are required for multiple packages, including *writegene2*, *metagene*, *genelist*, *smorflist*, *posavg* and *posstats*. Standard use of these packages assumes RPM normalization option has been selected. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/builddense.png)

We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript to generate 5’-mapped ROCC files with data normalized in units of RPM:
```unix
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $ROCC/80S -1 25 34 1 > $ROCC/80S_bd_metadata.txt 
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $ROCC/40S -1 20 80 1 > $ROCC/40S_bd_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_builddense.png)







## SCRIPT 2. Run *writegene2* from ROCC files
*writegene2* outputs RPM values derived from reads that aligned to a particular gene model.

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/writegene2.png)

 We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript, generating data on *ACTB* and *EIF4G2 *from ROCC files:
```unix
python $CODE/writegene2.py $ROCC/80S.rocc ACTB $DATA/writegene2/80S_wg2_ACTB > $DATA/writegene2/80S_wg2_ACTB_metadata.txt
python $CODE/writegene2.py $ROCC/80S.rocc EIF4G2 $DATA/writegene2/80S_wg2_EIF4G2 > $DATA/writegene2/80S_wg2_EIF4G2_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_writegene2.png)







## SCRIPT 3. Run *genelist* from ROCC files
*genelist* counts reads that map to CDS and UTR regions of gene models, and individual frames of the CDS. The outputs are provided in RPKM and raw counts units. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/genelist.png)

We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript:
```unix
python $CODE/genelist.py $ROCC/80S.rocc 12 1 $DATA/genelist/80S_gl > $DATA/genelist/80S_gl_metadata.txt 
python $CODE/genelist.py $ROCC/40S.rocc 12 1 $DATA/genelist/40S_gl > $DATA/genelist/40S_gl_metadata.txt 
```
The 5'-end aligned data is shifted to accommodate the P-site of the ribosome by using a shift of 12. Note that genes with 5'-UTRs shorter than the shift value are excluded from the analysis.

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/shift.png)

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_genelist.png)






## SCRIPT 4. Run *posavg* from ROCC files
*posavg* averages ribosome profiling data around any sequence feature of interest or computes a pause score for every occurrence of all 61 codons or 20 amino acids. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posavg1.png)

We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript.
Calculation of pause scores for all amino acids in the main frame for the CDS:
```unix
python $CODE/posavg.py $ROCC/80S.rocc all_1 1 0 0 30 0 12 1 none $DATA/posavg/80S_pa_all_1_aa_CDS_frame0 > $DATA/posavg/80S_pa_all_1_aa_CDS_frame0_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posavg2.png)

uORF analysis with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm0 > $DATA/posavg/80S_pa_uORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm1 > $DATA/posavg/80S_pa_uORF_norm1_metadata.txt 
```

dORF analysis with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm0 > $DATA/posavg/80S_pa_dORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm1 > $DATA/posavg/80S_pa_dORF_norm1_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_posavg.png)






## SCRIPT 5. Run *metagene* from ROCC files
*metagene* averages read counts over transcripts for a windowed region around start or stop codons.
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/metagene2.png)

How the code works: 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/metagene1.png)


We used the following settings, and redirect (via > operator) for saving of metadata, to generate a start and stop codon metagene in this study:
```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 none $DATA/metagene/80S_mg_start > $DATA/metagene/80S_mg_start_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 none $DATA/metagene/80S_mg_stop > $DATA/metagene/80S_mg_stop_metadata.txt
```

The subset excel file to limit analysis to genes with predicted uORFs. Setting to generate uORF-containing metagenes include:

```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 $FILES/subset_list.xlsx $DATA/metagene/80S_mg_start_subset > $DATA/metagene/80S_mg_start_subset_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 $FILES/subset_list.xlsx $DATA/metagene/80S_mg_stop_subset > $DATA/metagene/80S_mg_stop_subset_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_metagene.png)






## SCRIPT 6. Run *smorflist* from ROCC files
*smorflist* counts reads that map in frame to uORFs or dORFs in the 5’ or 3’-UTR, respectively. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/smorflist.png)

We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript identifying uORFs within the 5’-UTR of ribosome profiling data:
```unix
python $CODE/smorflist.py $ROCC/80S.rocc 4 12 0 0 5 $DATA/smorflist/80S_uorflist > $DATA/smorflist/80S_uorflist_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_smorflist.png)






## SCRIPT 7. Run *posstats* from ROCC files
posstats computes and reports pause scores for every occurrence of a sequence feature of interest. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/posstats.png)
We used the following settings, and redirect (via > operator) for saving of metadata, for data analysis in this manuscript to calculate the distribution of pause scores for PPG and AAA tri-amino acid motifs in our ribosome profiling dataset:
```unix
python $CODE/posstats.py $ROCC/80S.rocc PPG 1 0 10 1 12 1 $DATA/posstats/80S_ps_PPG >  $DATA/posstats/80S_ps_PPG_metadata.txt
python $CODE/posstats.py $ROCC/80S.rocc AAA 1 0 10 1 12 1 $DATA/posstats/80S_ps_AAA >  $DATA/posstats/80S_ps_AAA_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_posstats.png)



## SCRIPT 8. Run *region_size_and_abundance* from SAM files
*region_size_and_abundance* outputs the number of reads that aligned to each transcript feature (i.e. UTRs, ORF, start and stop codon) for different read lengths. The saved output consists of a csv file that includes separate sections for read length distributions (left; raw counts for each of the 5 transcript regions) and multiple kinds of abundance measurement for each transcript region (right; can be used to normalize the read length distributions by multiplying the raw count values). 
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/rsa2.png)

The window parameter can be adjusted. This setting determines the regions around the start and stop codon in which a mapped read needs to overlap:
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/rsa1.png)

The example code outputs data for ribosome profiling data for all transcripts: 
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/region_size_and_abundance/80S_rsa 25 34 4 none > $DATA/region_size_and_abundance/80S_rsa_metadata.txt
```

The subset excel file is used to limit analysis to *EIF4G2*:
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/region_size_and_abundance/80S_rsa_subset 25 34 4 $FILES/subset_list_EIF4G2.xlsx > $DATA/region_size_and_abundance/80S_rsa_subset_metadata.txt
```

The example code outputs data for 40S-seq data for all transcripts: 
```unix
python $CODE/region_size_and_abundance.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $DATA/region_size_and_abundance/40S_rsa 20 80 4 none > $DATA/region_size_and_abundance/40S_rsa_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_rsa.png)





## SCRIPT 9. Run *metagene_3D* from SAM files
*metagene_3D* averages reads around start or stop codons as a function of read length. The output consists of 2 csv files with 3-D metagene data (5’- and 3’-assigned reads) and 1 csv file with conventional metagene data (both 5’- and 3’-assigned included). The example code below generates plots used in this study (around start codons). Plots are created for both 5’- and 3’-assigned reads. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/3D_metagene.png)

Ribosome profiling data on start and stop codons:
```unix
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_start none 25 34 30 75 1 > $DATA/metagene_3D/80S_3dmg_start_metadata.txt
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_stop none 25 34 75 30 2 > $DATA/metagene_3D/80S_3dmg_stop_metadata.txt
```

40S-seq data on start codons:
```unix
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $DATA/metagene_3D/40S_3dmg_start none 20 80 50 50 1 > $DATA/metagene_3D/40S_3dmg_start_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_3D_metagene.png)

Plots can be generated from the 3D output csv files using the script below that generates a matplotlib window with controls to view/save the metagene image. The read abundance is shown by a color scheme. 
Ribosome profiling plots:
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/80S_3dmg_start_3Dmeta_end3.csv 
```
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/80S_3dmg_start_3Dmeta_end5.csv 
```
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/80S_3dmg_stop_3Dmeta_end3.csv 
```
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/80S_3dmg_stop_3Dmeta_end5.csv 
```


40S plots:
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/40S_3dmg_start_3Dmeta_end3.csv 
```
```unix
python $CODE/metagene_3D_plot.py $DATA/metagene_3D/40S_3dmg_start_3Dmeta_end5.csv 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_tables/settings_3D_metagene_plot.png)

