# ribofootPrinter2.1
ribofootPrinter2 is a package of Python 3 tools for analyzing data from ribosome profiling and other RNA-based sequencing experiments. This is a new version of the [first ribofootPrinter](https://github.com/guydoshlab/ribofootPrinter) package, which was originally based on the [ReducedTranscriptome](https://github.com/guydoshlab/ReducedTranscriptome) package. 

ribofootPrinter2 adds support for version 1.4 of the MANE transcriptome. It updates the functionality of all the previous core tools and adds several new capabilities. For guidance on use of the package, please see the associated manuscript on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.07.04.451082

We reccomend going to the "releases" page and downloading a zip file of the last stable release for use. Releases are also archived on Zenodo. If there are code updates since the last formal release, they will be available here on the main page as a beta version (not fully tested).

# Download ribofootPrinter from Github
To begin, we suggest starting with the entire *ribofootPrinter* folder that is created after downloading and unzipping the compressed file available on the releases page. 

Move the folder to a convenient place to work with it (such as your home directory) and navigate to the folder from the command line using the *cd* comand. Test this by typing the following:
```unix
pwd
```
You should see the prompt return the name of the folder: ribofootPrinter2.1

# Requirements for running ribofootPrinter
## Environment
We run the Python toolbox ribofootPrinter locally inside a virtual environment. 

A virtual environment can be created using the Python *venv* command (i.e. *python -m venv venv*). The .py files can then be copied to the site-packages folder associated with the virtual environment to facilitate use. Other needed packages (below) can be installed with *pip*.

Alternatively, the conda package manager can be used to handle the virtual environment.

The ribofootPrinter package has the following dependencies:

python v3.9.23\
biopython v1.85\
matplotlib v3.9.4\
pandas v2.3.2\
numpy v2.0.2\
openpyxl v3.1.5

## Alignment files
*ribofootPrinter* is compatible with bowtie1 aligned SAM files. We use the bowtie indexed version of the MANEv1.4_shortnames.fasta reduced transcriptome file for alignments (EBWT file extension), found [here](https://github.com/guydoshlab/ribofootPrinter2/tree/main/preparation/MANE_v1.4_Preparation/output_files) (unzip before use).

An overview on how to run the different *ribofootPrinter* scripts is listed below. This includes example code which can be used directly on the data files (SAM and ROCC) provided [here](https://zenodo.org/records/17917807). 

More details on preparation of the MANE transcriptome for alignment (shortnames FASTA file) and ribofootPrinter (longnames FASTA file) can be found [here](https://github.com/guydoshlab/ribofootPrinter2/tree/main/preparation/MANE_v1.4_Preparation/).

More details how to view your aligned reads in IGV can be found [here](https://github.com/guydoshlab/ribofootPrinter2/tree/main/preparation/MANE_v1.4_IGV).

# Running the *ribofootPrinter* package
All scripts are run in the terminal by calling Python and providing the Python script and settings as described below. An overview of the different settings is provided in table format for each script.
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


![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/ribofootprinter.png)


What follows below is example code for running each script on the provided test data.


## Generate folders for data files
```unix
mkdir -p ./rocc_files/
mkdir -p ./data_files/ 
mkdir -p ./sam_files/
```


## Create variables for different paths
Provide full path to create variables. The full path can be obtained through: 
```unix
FILEPATH=$(pwd)
```

```unix
CODE=/$FILEPATH/code
FASTA=/$FILEPATH/preparation/MANE_v1.4_Preparation/output_files
SAM=/$FILEPATH/sam_files
ROCC=/$FILEPATH/rocc_files
DATA=/$FILEPATH/data_files
FILES=/$FILEPATH/associated_files
```

```unix
mkdir -p ./data_files/writegene2 ./data_files/metagene ./data_files/genelist ./data_files/smorflist ./data_files/posavg ./data_files/posstats ./data_files/region_size_and_abundance ./data_files/metagene_3D
```
## Download SAM files:
```unix
wget https://zenodo.org/records/17917807/files/40S.SAM.zip
wget https://zenodo.org/records/17917807/files/80S.SAM.zip
```


## Unzip files
```unix
for file in $FASTA/*.zip; do unzip $file; done
for file in $FILES/*.zip; do unzip $file; done
for file in $SAM/*.zip; do unzip $file; done
```
```unix
mv *.fasta $FASTA
mv *.ebwt $FASTA  
mv *.xlsx $FILES 
mv *.SAM $SAM
```






## SCRIPT 1. Run *builddense* to generate ROCC files from SAM files
*builddense* converts SAM files into the ROCC files that are required for multiple packages, including *writegene2*, *metagene*, *genelist*, *smorflist*, *posavg* and *posstats*. Standard use of these packages assumes RPM normalization option has been selected. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/builddense.png)

We used the following settings, and redirect (via the > operator) to save metadata, for data analysis to generate 5’-mapped ROCC files with data normalized in units of RPM:
```unix
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $ROCC/80S -1 25 34 1 > $ROCC/80S_bd_metadata.txt 
python $CODE/builddense.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $ROCC/40S -1 20 80 1 > $ROCC/40S_bd_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_builddense.png)







## SCRIPT 2. Run *writegene2* from ROCC files
*writegene2* outputs RPM values derived from reads that aligned to a particular gene model.

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/writegene2.png)

 We used the following settings, and redirect (via the > operator) to save metadata, for data analysis to generate csv files with data for *ACTB* and *EIF4G2* from ROCC files:
```unix
python $CODE/writegene2.py $ROCC/80S.rocc ACTB $DATA/writegene2/80S_wg2_ACTB > $DATA/writegene2/80S_wg2_ACTB_metadata.txt
python $CODE/writegene2.py $ROCC/80S.rocc EIF4G2 $DATA/writegene2/80S_wg2_EIF4G2 > $DATA/writegene2/80S_wg2_EIF4G2_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_writegene2.png)







## SCRIPT 3. Run *genelist* from ROCC files
*genelist* counts reads that map to CDS and UTR regions of gene models, and individual frames of the CDS. The outputs are provided in RPKM and raw counts units. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/genelist.png)

We used the following settings, and redirect (via the > operator) to save metadata, for data analysis:
```unix
python $CODE/genelist.py $ROCC/80S.rocc 12 1 $DATA/genelist/80S_gl > $DATA/genelist/80S_gl_metadata.txt 
python $CODE/genelist.py $ROCC/40S.rocc 12 1 $DATA/genelist/40S_gl > $DATA/genelist/40S_gl_metadata.txt 
```
For this analysis, the 5'-end aligned data are shifted to accommodate the P-site of the ribosome by using a shift of 12. Note that genes with 5'-UTRs shorter than the shift value are excluded from the analysis.

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/shift.png)

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_genelist.png)






## SCRIPT 4. Run *posavg* from ROCC files
*posavg* averages ribosome profiling data around any sequence feature of interest or computes a pause score for every nucleotide or amino-acid combination for a sequence of given length. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/posavg1.png)

We used the following settings, and redirect (via the > operator) to save metadata, for data analysis to calculate pause scores for all amino acids in the ribosome P-site in the main frame of the CDS:
```unix
python $CODE/posavg.py $ROCC/80S.rocc all_1 1 0 0 30 0 12 1 none $DATA/posavg/80S_pa_all_1_aa_CDS_frame0 > $DATA/posavg/80S_pa_all_1_aa_CDS_frame0_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/posavg2.png)

ATG triplets in 5'-UTRs (uORF analysis) with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm0 > $DATA/posavg/80S_pa_uORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 0 none $DATA/posavg/80S_pa_ATG_UTR5_allframe_norm1 > $DATA/posavg/80S_pa_uORF_norm1_metadata.txt 
```

ATG triplets in 3'-UTRs (dORF analysis) with different normalization settings:
```unix
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 0 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm0 > $DATA/posavg/80S_pa_dORF_norm0_metadata.txt 
python $CODE/posavg.py $ROCC/80S.rocc ATG 0 3 0 30 1 12 2 none $DATA/posavg/80S_pa_ATG_UTR3_allframe_norm1 > $DATA/posavg/80S_pa_dORF_norm1_metadata.txt 
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_posavg.png)






## SCRIPT 5. Run *metagene* from ROCC files
*metagene* averages read counts over transcripts for a windowed region around start or stop codons.
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/metagene2.png)

How the code works: 
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/metagene1.png)


We used the following settings, and redirect (via the > operator) to save metadata, to generate a start and stop codon metagene plots:
```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 none $DATA/metagene/80S_mg_start > $DATA/metagene/80S_mg_start_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 none $DATA/metagene/80S_mg_stop > $DATA/metagene/80S_mg_stop_metadata.txt
```

The subset excel file limits analysis to genes with predicted uORFs. Settings to generate metagene plots from predicted uORF-containing genes:

```unix
python $CODE/metagene.py $ROCC/80S.rocc 1 1 5 50 300 $FILES/subset_list.xlsx $DATA/metagene/80S_mg_start_subset > $DATA/metagene/80S_mg_start_subset_metadata.txt
python $CODE/metagene.py $ROCC/80S.rocc 2 1 5 300 50 $FILES/subset_list.xlsx $DATA/metagene/80S_mg_stop_subset > $DATA/metagene/80S_mg_stop_subset_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_metagene.png)






## SCRIPT 6. Run *smorflist* from ROCC files
*smorflist* counts reads that map in frame to uORFs or dORFs in the 5’ or 3’-UTR, respectively. 
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/smorflist.png)

We used the following settings, and redirect (via the > operator) to save metadata, for data analysis to identify uORFs within the 5’-UTR of ribosome profiling data:
```unix
python $CODE/smorflist.py $ROCC/80S.rocc 4 12 0 0 5 $DATA/smorflist/80S_uorflist > $DATA/smorflist/80S_uorflist_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_smorflist.png)






## SCRIPT 7. Run *posstats* from ROCC files
posstats computes and reports pause scores for every occurrence of a sequence feature of interest. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/posstats.png)
We used the following settings, and redirect (via the > operator) to save metadata, to calculate the distribution of pause scores for PPG and AAA tri-amino acid motifs:
```unix
python $CODE/posstats.py $ROCC/80S.rocc PPG 1 0 10 1 12 1 $DATA/posstats/80S_ps_PPG >  $DATA/posstats/80S_ps_PPG_metadata.txt
python $CODE/posstats.py $ROCC/80S.rocc AAA 1 0 10 1 12 1 $DATA/posstats/80S_ps_AAA >  $DATA/posstats/80S_ps_AAA_metadata.txt
```

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_posstats.png)



## SCRIPT 8. Run *region_size_and_abundance* from SAM files
*region_size_and_abundance* outputs the number of reads that aligned to each transcript feature (i.e. 5'-UTR, 3'-UTR, main ORF, start codon, and stop codon) as a function of read length. The saved output consists of a csv file that includes separate sections for read length distributions (left; raw counts for each of the 5 transcript regions) and multiple kinds of abundance measurement for each transcript region (right; can be used to normalize the read length distributions by multiplying the raw count values). 
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/rsa2.png)

The window parameter can be adjusted. This setting determines the regions around the start and stop codon in which a mapped read needs to overlap:
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/rsa1.png)

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

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_rsa.png)





## SCRIPT 9. Run *metagene_3D* from SAM files
*metagene_3D* averages reads around start or stop codons as a function of read length. The output consists of 2 csv files with 3-D metagene data (one file each for 5’- and 3’-assigned reads) and 1 csv file with conventional metagene data (both 5’- and 3’-assigned included in the file). The example code below generates plots around start codons. Plots are created for both 5’- and 3’-assigned reads. 

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_figures/3D_metagene.png)

Ribosome profiling data on start and stop codons:
```unix
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_start none 25 34 30 75 1 > $DATA/metagene_3D/80S_3dmg_start_metadata.txt
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/80S.SAM $DATA/metagene_3D/80S_3dmg_stop none 25 34 75 30 2 > $DATA/metagene_3D/80S_3dmg_stop_metadata.txt
```

40S-seq data on start codons:
```unix
python $CODE/metagene_3D.py $FASTA/MANEv1.4_longnames.fasta $SAM/40S.SAM $DATA/metagene_3D/40S_3dmg_start none 20 80 50 50 1 > $DATA/metagene_3D/40S_3dmg_start_metadata.txt
```
![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_3D_metagene.png)

Plots can be generated from the 3D output csv files using the script below that generates a matplotlib window with controls to view/save the metagene image. The read abundance is shown by a color scheme. 
Ribosome profiling (80S) plots:
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

![alt text](https://github.com/guydoshlab/ribofootPrinter2/blob/main/Github_tables/settings_3D_metagene_plot.png)

