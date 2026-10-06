import os
import pandas as pd
from Bio import SeqIO

workingfolder = input("Enter folder path where MANE files are located: ")

## Navigate to MANE folder containing all the downloaded files
os.chdir(workingfolder)

df_names = pd.read_csv('MANE.GRCh38.v1.4.summary.txt', sep="\t", header=0)  # Read in txt file which contains genename and geneID information
df_CDS = pd.read_csv('CDS_information_refseq.txt', sep="\t", header=None)   # Read in CDS information txt file (obtained from .GBFF file using bash script)

with open('MANE.GRCh38.v1.4.ensembl_rna.fna') as fasta_file:  # Will close handle cleanly
    identifiers = []
    sequence = []
    for record in SeqIO.parse(fasta_file, 'fasta'):  # Read in fasta file and store genename and sequence as lists.
        identifiers.append(record.id)
        sequence.append(record.seq)   

sequence_df = pd.DataFrame({"genename": identifiers, "sequence": sequence}) 
sequence_df = sequence_df.set_index("genename") # set NM transcriptIDs as index for sorting
sequence_df = sequence_df.sort_index()

print("The three input files have been uploaded and contain " + str(len(df_names)) + " entries")


## Formatting file containing transcript names 
df_names = df_names[df_names["MANE_status"] == 'MANE Select'] # Filter out MANE Plus Clinical transcripts, only keep MANE select transcripts to avoid isoform duplicate genenames
print("MANE Plus Clinical transcripts have been filtered out, resulting in " + str(len(df_names)) + " entries")

transcriptnames = df_names[["Ensembl_Gene","Ensembl_nuc","Ensembl_nuc", "Ensembl_prot","RefSeq_nuc","RefSeq_nuc","RefSeq_prot","symbol"]] # Filter columns to retain geneID and genename information
transcriptnames.columns = ["Ensembl_Gene","Ensembl_nuc","Ensembl_nuc2","Ensembl_prot","RefSeq_nuc_short","RefSeq_nuc","RefSeq_prot","symbol"]
transcriptnames['RefSeq_nuc_short'] = transcriptnames['RefSeq_nuc_short'].str.split('.').str[0] # Removes everything after the refseq names, e.g. NM_130786.4 becomes NM_130786

transcriptnames = transcriptnames.set_index("Ensembl_nuc2") # set enseble transcriptIDs as index for merging with sequence file.
transcriptnames = transcriptnames.sort_index()

transcriptnames = pd.merge(transcriptnames, sequence_df, left_index=True, right_index=True) # Add transcript sequences to dataframe.
transcriptnames["sequence"] = transcriptnames["sequence"].astype(str)   # Convert sequences from object to string.

transcriptnames = transcriptnames.set_index("RefSeq_nuc_short") # set short NM transcriptIDs as index for sorting
transcriptnames = transcriptnames.sort_index()

transcriptnames = transcriptnames[transcriptnames.index.str.contains('NM_')] # Removes non-coding RNAs that start with NR_ refseq names.

print("Non-coding RNAs have been filtered out, resulting in " + str(len(transcriptnames)) + (" entries"))
#[19288 rows x 6 columns]


## Formatting file containing CDS information 
df_CDS_filtered = df_CDS[df_CDS[4] != "RNA"]    # Filter out any non-coding RNAs. Retains mRNA (first line) and CDS information (line below first line) only.
df_CDS_filtered_1 = df_CDS_filtered[df_CDS_filtered.index % 2 == 0] # Retains rows with transcriptID and transcriptlength information (i.e. first line)
df_CDS_filtered_2 = df_CDS_filtered[df_CDS_filtered.index % 2 != 0] # Retains rows with CDS coordinates information (i.e. line below first line)

CDSinformation_1 = df_CDS_filtered_1.iloc[:, [1,2]]     # Subsets file to only contain refseq name and transcript length
CDSinformation_1.columns = ["ID", "transcriptlength"]   # Adds column titles
CDSinformation_1.index = range(len(CDSinformation_1.index)) # Resets index

CDSinformation_2 = df_CDS_filtered_2.iloc[:, [2,3]]     # Subsets file to only contain CDS start and CDS end information.
CDSinformation_2.columns = ["CDS_start", "CDS_end"]     # Adds column titles
CDSinformation_2.index = range(len(CDSinformation_2.index)) # Resets index

CDSinformation = pd.concat([CDSinformation_1, CDSinformation_2],ignore_index=False, sort=False, axis=1) # merges transcriptID, transcript length and CDS information into one file. axis = 1 to concatenate by column

CDSinformation = CDSinformation.set_index("ID") # set NM transcriptIDs as index for sorting
CDSinformation = CDSinformation.sort_index()

CDSinformation = CDSinformation[CDSinformation.index.str.contains('NM_')]

print(CDSinformation)
print(transcriptnames)

longnames = pd.merge(transcriptnames, CDSinformation, left_index=True, right_index=True)
#print(longnames)



## Calculate UTR5, CDS, UTR3 variable to store in file
name1 = longnames.iloc[:, 0].to_list() 
name2 = longnames.iloc[:, 1].to_list() 
name3 = longnames.iloc[:, 2].to_list() 
name4 = longnames.iloc[:, 3].to_list() 
name5 = longnames.iloc[:, 4].to_list() 
name6 = longnames.iloc[:, 5].to_list() 
seq = longnames.iloc[:, 6].to_list() 
transcriptlength = longnames.iloc[:, 7].to_list() 

CDS_start = longnames.iloc[:, 8].to_list()
CDS_end = longnames.iloc[:, 9].to_list()

UTR5_end_df = longnames["CDS_start"].astype(int) - 1
UTR3_start_df = longnames["CDS_end"].astype(int) + 1
UTR5_end = UTR5_end_df.iloc[:, ].to_list()
UTR3_start = UTR3_start_df.iloc[:, ].to_list()

listlength = len(CDS_start)

UTR5_string = ["UTR5:1-"] * listlength
UTR5 = list(zip(UTR5_string, UTR5_end))
UTR5 = [''.join(map(str, item)) for item in UTR5]

CDS_string_1 = ["CDS:"] * listlength
CDS_string_2 = ["-"] * listlength
CDS = list(zip(CDS_string_1, CDS_start, CDS_string_2, CDS_end))
CDS = [''.join(map(str, item)) for item in CDS]

UTR3_string_1 = ["UTR3:"] * listlength
UTR3_string_2 = ["-"] * listlength
UTR3 = list(zip(UTR3_string_1, UTR3_start, UTR3_string_2,transcriptlength))
UTR3 = [''.join(map(str, item)) for item in UTR3]

pipe = ["|"] * listlength
fastasymbol = [">"] * listlength
newline = ['\n'] * listlength

longnames_full = pd.DataFrame(list(zip(
    fastasymbol,
    name1, pipe,
    name2, pipe,
    name3, pipe,
    name4, pipe,
    name5, pipe,
    name6, pipe,
    transcriptlength, pipe,
    UTR5, pipe,
    CDS, pipe,
    UTR3, newline,
    seq
    )))

# Joining all genenames, UTR and CDS information into one column 
longnames_full["longname"] = longnames_full[longnames_full.columns[0:20]].apply(
    lambda x: ''.join(x.dropna().astype(str)),
    axis=1
)

longnames_full = longnames_full.iloc[:, [22,20,21]]     # Subsets file 
#print(longnames_full)


df = pd.DataFrame(longnames_full)
df.columns = [''] * len(df.columns)
df.to_csv('MANEv1.4_longnames.fasta', index=False)
