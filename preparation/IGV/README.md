# Step by step guide for viewing transcriptome aligned ribosome footprints in IGV
This code is written in bash and should be executed in the terminal.
IGV can be downloaded here:
https://igv.org/doc/desktop/#DownloadPage/

# General outline
1. Obtain the neccesary files:

      a. Obtain transcriptome FASTA shortnames and longnames files (for example, see MANE_v1.4_Preparation Github).

      b. Get SAM file(s) from your alignment to the shortnames FASTA file using software such as bowtie (see MANE_v1.4_Preparation Github for instructions).

2. Prepare alignment files (SAM -> BAM -> BEDGRAPH or BW).
3. Prepare GTF file from longnames FASTA to annotate the CDS of transcripts using Python script longnames_to_GTF.py.
4. Load outputs with IGV.
   
![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/MANE_IGV_viewing.png)

Note: It is important to that the FASTQ footprints files are aligned against the same transcriptome used for IGV. For example, we use the MANEv1.4_shornames.FASTQ reduced transcriptome for mapping and IGV. 

# Setting up the environment:
SAM files should be converted to sorted BAM files using samtools. Then, bedtools and deeptools are used to generate bedgraph or bigwig files, respectively. See the respective websites for download and install instructions. This can also be handled using package managers and virtual environments. 

# 1. Obtaining input files
Create transcriptome_IGV folder, then navigate to folder:
```unix
mkdir -p transcriptome_IGV
cd ./transcriptome_IGV
```
   a. Download transcriptome FASTA shortnames and longnames files to this folder (i.e. see MANE_v1.4_Preparation Github).

   b. Place  SAM file(s) from your alignment in this folder (i.e. 80S.SAM).

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
Note that the 28 nt arrow depicted corresponds to the 28 bonds within a read made up of 29 nt.
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
The GTF file functions as a lookup table and contains information on CDS boundries important to identify footprints outside the coding region. 
Navigate to folder:
```unix
cd ./transcriptome_IGV/
```
Run Python script longnames_to_GTF.py

# 4. Load files in IGV
## Load transcriptome sequence (FASTA)
Load, for example, the MANE transcriptome using following steps:

Genomes -> Load Genome from File... -> MANEv1.4_shortnames.fasta

This should automatically create the MANEv1.4_shortnames.fasta.fai file required by IGV.

## Load alignment files (BAM, BEDGRAPH or BW)
Drag in the files or upload using following steps:

File -> Load from File... -> 80S.bam

File -> Load from File... -> 80S.bedgraph

File -> Load from File... -> 80S.bw


When uploading a bam file make sure the .bam.bai index is located in the same folder.

## Load CDS annotations (GTF)
Drag in the files or upload using following steps:

File -> Load from File... -> MANEv1.4_CDS.gtf


