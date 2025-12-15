# Step by step guide for viewing MANEv1.4 transcriptome aligned ribosome footprints in IGV
This code is written in bash and should be executed in the terminal.
IGV can be downloaded here:
https://igv.org/doc/desktop/#DownloadPage/

# General outline
1. Obtain the 3 neccesary input files:

      a. Download transcriptome FASTA shortnames file (see MANE_v1.4_Preparation Github).

      b. Get SAM file(s) from your alignment to the shortnames FASTA file using software such as bowtie (see MANE_v1.4_Preparation Github for instructions). 

      c. Create or download CDS_information_refseq.txt file as described for MANE annotation (see MANE_v1.4_Preparation Github) (only needed for GTF file generation).
   
      d. Download MANE.GRCh38.v1.4.summary.txt from MANE website (only needed for GTF file generation).

3. Prepare alignment files (SAM -> BAM -> BEDGRAPH or BW).
4. Prepare GTF file from FASTA and CDS_information_refseq text file (or download here from output_files folder) to annotate the CDS of transcripts using Python script MANE_GTF_formatter.py.
5. Load outputs with IGV.
   
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/MANE_IGV_viewing.png)

Note: It is important to that the FASTQ footprints files are aligned against the same transcriptome used for IGV. In this example we use the MANEv1.4_shornames.FASTQ reduced transcriptome for mapping and IGV. 

# Setting up the environment:
SAM files should be converted to sorted BAM files using samtools. See the samtools website for download and install instructions. This can be done locally using a virtual environment, such as Anaconda, but not required. 
Reference for Anaconda: https://github.com/conda/conda?tab=readme-ov-file.
Details for Anaconda follow.
## Generate conda environment (only once)
```unix
conda create --name samtools
```
## Activate the environment (each time)
```unix
conda activate samtools
```
## Install seqtk, samtools, bedtools (only once)
https://anaconda.org/bioconda/samtools

https://anaconda.org/bioconda/bedtools

https://anaconda.org/bioconda/deeptools

```unix
conda install bioconda::samtools
```
```unix
conda install bioconda::bedtools
```
```unix
conda install bioconda::deeptools
```

# 1. Obtaining input files
Create MANE_IGV folder, then navigate to folder:
```unix
mkdir -p MANEv1.4_IGV
cd ./MANEv1.4_IGV
```
   a. Download transcriptome FASTA shortnames file to this folder (see MANE_v1.4_Preparation Github; MANEv1.4_shortnames.FASTA).

   b. Place  SAM file(s) from your alignment in this folder.
   
   c. Download CDS_information_refseq.txt file to this folder (see MANE_v1.4_Preparation Github) (only needed for GTF file generation).

   d. Download and unzip MANE.GRCh38.v1.4.summary.txt (only needed for GTF file generation).
   
```unix
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz
```
```unix
gunzip MANE.GRCh38.v1.4.summary.txt.gz
```

# 2. Convert alignment files to IGV compatible files

## Convert and sort SAM file into BAM file
Use samtools to convert SAM files into BAM files as follows:
```unix
samtools sort -o 80S.bam 80S.SAM
```
## Index BAM file to enable IGV viewing
This will index all bam files within the current directory. 
```unix
find *.bam -exec echo samtools index {} \; | sh
```

## Generate BEDGRAPH or BIGWIG files from BAM files
BEDGRAPH files contain reads that have been assigned a single count at either the 5' or 3' end. These end mapped reads can then be converted to have the single count at for example the P-site of ribosome footprint. 

The code below determines the normalization factor for each BAM file and generates BEDGRAPH files that contain a single count at the 5' end. If 3'-end assigned reads are desired, use -bg -3 (instead of -bg -5).
```unix
for file in *.bam
do
	echo $file
	count=$(samtools view -F 4 $file | wc -l | xargs)
	echo $count
	scalefactor=$(echo "1000000 / $count" | bc -l)
	echo $scalefactor
	bedtools genomecov -ibam $file -scale $scalefactor -bg -5 > ${file%.bam}.bedgraph	
done
```

## Shifting BEDGRAPH files 
If desired, the optional code below will shift the assigned 5' read positions by +12 which generally aligns with the P-site in riboseq footprints.
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/shift.png)
```unix
for i in *.bedgraph;
do awk '{print $1, $2+12, $3+12, $4}' $i > ${i%.bedgraph}_shiftadd12.bedgraph;
done
```

## Coverage (instead of end mapped reads)
BAM files can also be converted to generate BIGWIG or BW files which have uniform coverage across the read. Coverage plots are generated using Deeptools using code below.

Note: this code runs for a long time.

```unix
for file in *.bam
do
	echo $file	
	bamCoverage -b ./$file -o ./${file%.bam}.bw -bs 1 --normalizeUsing BPM
done
```


# 3. Prepare GTF file using Python script
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/MANE_transcriptome_GTF.png)
The GTF file functions as a lookup table and contains information on CDS boundries important to identify footprints outside the coding region. 
Navigate to folder:
```unix
cd ./MANEv1.4_IGV/
```
Downloaded files inside MANEv1.4_IGV folder:

1. MANE.GRCh38.v1.4.summary.txt (downloaded above)
2. CDS_information_refseq.txt (see MANE_v1.4_Preparation Github) 

Run Python script MANE_GTF_formatter.py

Note: This script has dependencies on Pandas. 

As inputs, it will ask for "workingfolder" to match the location of the MANEv1.4 folder (e.g. "/Users/yourusername/Desktop/MANEv1.4") where it can find MANE.GRCh38.v1.4.summary.txt and CDS_information_refseq.txt.


Once the GTF file is obtained, run following script to remove unwanted characters.
```unix
sed -i -e 's/"""/""/g' ./MANEv1.4_CDS.gtf
sed -i -e 's/""/"/g' ./MANEv1.4_CDS.gtf
sed -i -e 's/"gene_id/gene_id/g' ./MANEv1.4_CDS.gtf
sed '1d' ./MANEv1.4_CDS.gtf > MANEv1.4_CDS_tmp.gtf
mv MANEv1.4_CDS_tmp.gtf MANEv1.4_CDS.gtf
rm MANEv1.4_CDS.gtf-e
```

# 4. Load files in IGV
## Load transcriptome sequence (FASTA)
Load the MANE transcriptome using following steps:

Genomes -> Load Genome from File... -> MANEv1.4_shortnames.fasta

This should automatically create the MANEv1.4_shortnames.fasta.fai file required by IGV.

## Load alignment files (BAM, BEDGRAPH or BW)
Drag in the files or upload using following steps:

File -> Load from File... -> example.bam

File -> Load from File... -> example.bedgraph

File -> Load from File... -> example.bw


When uploading a bam file make sure the .bam.bai index is located in the same folder.

## Load CDS annotations (GTF)
Drag in the files or upload using following steps:

File -> Load from File... -> MANEv1.4_CDS.gtf


