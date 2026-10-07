import os
import pandas as pd

origworkingfolder = os.getcwd()
#workingfolder = input("Enter folder path where MANE files are located: ")
workingfolder = ("/Users/kerkhofsk2/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/3_Papers/2025-09_Kerkhofs_Guydosh_Ribofootprinter/3.NARGAB_reviews/MANEv1.4_IGV")

## Navigate to MANE folder containing all the downloaded files
os.chdir(workingfolder)

df_names = pd.read_csv('MANE.GRCh38.v1.4.summary.txt', sep="\t", header=0)  # Read in txt file which contains genename and geneID information
df_CDS = pd.read_csv('CDS_information_refseq.txt', sep="\t", header=None)   # Read in CDS information txt file (obtained from .GBFF file using bash script)

## Formatting file containing transcript names 
#df_names = df_names[df_names["MANE_status"] == 'MANE Select'] # Filter out MANE Plus Clinical transcripts, only keep MANE select transcripts to avoid isoform duplicate genenames
transcriptnames = df_names[["Ensembl_Gene","Ensembl_nuc", "Ensembl_prot","RefSeq_nuc","RefSeq_nuc","RefSeq_prot","symbol"]] # Filter columns to retain geneID and genename information
transcriptnames.columns = ["Ensembl_Gene","Ensembl_nuc","Ensembl_prot","RefSeq_nuc_short","RefSeq_nuc","RefSeq_prot","symbol"]
transcriptnames['RefSeq_nuc_short'] = transcriptnames['RefSeq_nuc_short'].str.split('.').str[0] # Removes everything after the refseq names, e.g. NM_130786.4 becomes NM_130786

transcriptnames = transcriptnames.set_index("RefSeq_nuc_short") # set short NM transcriptIDs as index for sorting
transcriptnames = transcriptnames.sort_index()

transcriptnames = transcriptnames[transcriptnames.index.str.contains('NM_')] # Removes non-coding RNAs that start with NR_ refseq names.
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

#print(CDSinformation)
#print(transcriptnames)

longnames = pd.merge(transcriptnames, CDSinformation, left_index=True, right_index=True)
print(longnames)



## Calculate UTR5, CDS, UTR3 variable to store in file
name1 = longnames.iloc[:, 0].to_list() 
name2 = longnames.iloc[:, 1].to_list() 
name3 = longnames.iloc[:, 2].to_list() 
name4 = longnames.iloc[:, 3].to_list() 
name5 = longnames.iloc[:, 4].to_list() 
name6 = longnames.iloc[:, 5].to_list() 
#transcriptlength = longnames.iloc[:, 6].to_list() 

CDS_start = longnames.iloc[:, 7].to_list()
CDS_end = longnames.iloc[:, 8].to_list()

CDS_end_GTF_df = longnames["CDS_end"].astype(int) - 3   # GTF files require CDS end at the start of stop codon
CDS_end_GTF = CDS_end_GTF_df.iloc[:, ].to_list()

listlength = len(CDS_start)

quotation = list(['\"'] * listlength)

name1_string = ["gene_id"] * listlength
name1_info = list(zip(name1_string, name1))
name1_info = [' \"'.join(map(str, item)) for item in name1_info]
name1_info = list(zip(name1_info, ['"'] * listlength))
name1_info = [''.join(map(str, item)) for item in name1_info]

name2_string = ["transcript_id"] * listlength
name2_info = list(zip(name2_string, name2))
name2_info = [' \"'.join(map(str, item)) for item in name2_info]
name2_info = list(zip(name2_info, ['"'] * listlength))
name2_info = [''.join(map(str, item)) for item in name2_info]

name7_info = ["gene_type \"protein_coding\""] * listlength

name6_string  = ["gene_name"] * listlength
name6_info = list(zip(name6_string, name6))
name6_info = [' \"'.join(map(str, item)) for item in name6_info]
name6_info = list(zip(name6_info, ['"'] * listlength))
name6_info = [''.join(map(str, item)) for item in name6_info]

name8_info = ["transcript_type \"protein_coding\""] * listlength

name3_string = ["protein_id"] * listlength
name3_info = list(zip(name3_string, name1))
name3_info = [' \"'.join(map(str, item)) for item in name3_info]
name3_info = list(zip(name3_info, ['"'] * listlength))
name3_info = [''.join(map(str, item)) for item in name3_info]

name4_string = ["db_xref RefSeq"] * listlength
name4_info = list(zip(name4_string, name4))
name4_info = [' \"'.join(map(str, item)) for item in name4_info]
name4_info = list(zip(name4_info, ['"'] * listlength))
name4_info = [''.join(map(str, item)) for item in name4_info]

geneinfo = longnames_full = pd.DataFrame(list(zip(
    name1_info, name2_info, name7_info, name6_info, name8_info, name3_info, name4_info
    )))

# Joining all geneinformation into one column 
geneinfo["longname"] = geneinfo[geneinfo.columns[0:6]].apply(
    lambda x: '; '.join(x.dropna().astype(str)),
    axis=1
)

source = ["MANE_v1.4"] * listlength
feature = ["CDS"] * listlength
dot = ["."] * listlength
plus = ["+"] * listlength
zero = ["0"] * listlength


gtf_file = pd.DataFrame(list(zip(
    name2, source, feature, CDS_start, CDS_end,     # Change name1 to name2 to switch between ENSG and ENST names, respectively
    dot, plus, zero,
    geneinfo["longname"]
    )))

#gtf_file.columns = ["geneID", "source", "feature", "start","end","6","7","8","geneinfo"]   # Adds column titles

#print(longnames_full)


df = pd.DataFrame(gtf_file)
#df.columns = [''] * len(df.columns)
os.chdir(origworkingfolder)
df.to_csv('MANEv1.4_CDS_transcriptID.gtf', sep='\t', index=False)




