# Step by step guide for MANEv1.4 transcriptome preparation for *ribofootprinter*
This code is written in bash and should be executed in the terminal. In addition, a Python script is used to generate the transcriptome. 
MS Excel is also used for formatting.

Note: The longnames and shortnames FASTA files are uploaded in the output_files folder as zip compressed version and can be downloaded [here](https://github.com/guydoshlab/ribofootPrinter2.0-beta/tree/main/preparation/MANE_v1.4_Preparation/output_files).

# General outline
1. Download transcriptome files from MANE website (TXT, GBFF, FASTA).
2. Obtain CDS and transcript length information for each transcript from GBFF file. The output file is a TXT file.
3. Run Python script "MANE_longnames_formatter.py" to generate longnames transcriptome file which contains ENSEMBL and REFSEQ GeneIDs, GeneNames, 5'-UTR, CDS and 3'-UTR start and end positions, gene length and sequence.
4. Convert MANE_longnames (used with Ribofootprinter) into MANE_shortnames (used for alignment).

![alt text](https://github.com/guydoshlab/ribofootPrinter2.0-beta/blob/main/Github_figures/MANE_preparation.png)


# 1. Download files from MANE transcriptome website
MANE transcriptome website:
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/

More information on the MANE transcriptome:
https://www.ncbi.nlm.nih.gov/refseq/MANE/

Create MANE folder and download 3 files listed below here:
```unix
mkdir -p MANEv1.4
cd MANEv1.4
```

### FASTA - Transcript sequences with ENSEMBL identifiers in FASTA format (.fna.gz)
This file contains 19404 transcripts in FASTA format. 
```unix
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_rna.fna.gz
```

### GBFF - Transcript sequences with REFSEQ identifiers in NCBI/nucleotide format.
This file is important to obtain CDS information for each gene.
```unix
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_rna.gbff.gz
```

### TXT - Overview file.
This file contains additional information for each transcript. 
Note the presence of 66 MANE Plus Clinical transcripts. These transcripts cause duplicates in GeneID and will cause problems in downstream steps. These transcripts will be removed from the FASTA transcriptome later on.
```unix
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz  
```

### Unzip all files: 
```unix
gunzip *.gz
```

# 2. CDS and transcript length information
This script obtains transcriptID information, transcripts length and CDS coordinates.
```unix
cat MANE.GRCh38.v1.4.refseq_rna.gbff | grep '^LOCUS \|CDS    ' > CDS_information_refseq.txt
```

Open this txt file in Excel in order to modify formatting. Instructions below texted on Mac version of Excel.
1) Open the CDS_information_refseq.txt file in Excel using the text import wizard.
2) Choose Delimited
3) Check boxes for:
   
   Space
   
   Other, using "." as the delimiter.
   
   (Note that box for treat consecutive delimiters as one should be checked)
5) Finish import and resave the file over the old name (tab-delimited text).

# 3. Run Python script "MANE_longnames_formatter.py". 
This script takes three files as input and outputs the longnames_transcriptome FASTA file. 

Note: This script has dependencies on Biopython and Pandas.

Note: Update "workingfolder" to match the location of the MANEv1.4 folder (e.g. workingfolder = "/Users/yourusername/Desktop/MANEv1.4")

Note: This script filters out MANE Plus Clinical and non-coding transcripts.

The Python script reads in the "MANE.GRCh38.v1.4.summary.TXT" file and reformats to obtain
1. ENSEMBL_GeneID (e.g. ENSG00000121410.12)
2. ENSEMBL_TranscriptID (e.g. ENST00000263100.8)
3. ENSEMBL_ProteinID (e.g. ENSP00000263100.2)
4. REFSEQ_TranscriptID (e.g. NM_130786.4)
5. REFSEQ_ProteinID (e.g. NP_570602.2)
6. GeneName / Alias (e.g. A1BG)

The Python script reads in the "CDS_information_refseq.TXT" file and reformats to obtain CDS and transcript length information. 
5'-UTR and 3'-UTR start and end coordinates are calculated based on CDS and transcript lenght information.

7. Transcriptlength
8. UTR5
9. CDS
10. UTR3

The Python scripts reads in the FASTA (or FNA) file and extracts the sequences for each transcript. 

11. "nextline"
12. Sequence

It merges all these items together into the expected longnames_transcriptome format needed for Ribofootprinter. 

Once the file is obtained, run following bash script to remove unwanted characters.
```unix
sed -e 's/,//g' ./MANEv1.4_longnames.fasta > MANEv1.4_longnames_tmp.fasta
sed -e 's/"//g' ./MANEv1.4_longnames_tmp.fasta > MANEv1.4_longnames_tmp2.fasta
sed '1d' ./MANEv1.4_longnames_tmp2.fasta > MANEv1.4_longnames.fasta
rm MANEv1.4_longnames_tmp*.fasta
```
    
# 4. Convert MANE_longnames into MANE_shortnames
MANE_longnames is used for Ribofootprinter, while MANE_shortnames is used for mapping with alignement software such as Bowtie.
```unix
sed -r 's/\|.+//' MANEv1.4_longnames.fasta > MANEv1.4_shortnames.fasta
```

For alignments, bowtie executable can be run from the command line following download. Alternatively, it can be installed with Anaconda: https://anaconda.org/bioconda/bowtie

Bowtie requires a Bowtie index for sequence alignments. This will generate 6 files with suffixes .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt. These are the only files required for alignments. The original fasta file will not be used.

```unix
bowtie-build MANEv1.4_shortnames.fasta MANEv1.4
```

The fastq files can now be aligned against the MANEv1.4 shortnames transcriptome using Bowtie. For more information, see the Bowtie manual: https://bowtie-bio.sourceforge.net/manual.shtml. 

