# Transcriptome for MANEv1.4 including MANE Select Plus Clinical transcripts
The FASTA files provided here can be used to run *ribofootPrinter2* on the MANEv1.4 transcriptome which includes MANE Clinical isoforms. The longnames FASTA file was generated using the Python script included here.

This version of the MANE transcriptome uses transcriptIDs (ENST) as lookup values, as opposed to the ENSG nomenclature used in the original 
MANEv1.4 transcriptome which would result in multiple identical entries caused by isoforms (same geneID, different transcriptID).

The matching GTF is also deposited in this folder. The longnames_to_GTF.py file can be modified to generate this GTF with this alternate structure. 

