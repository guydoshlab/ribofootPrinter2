# This folder contains additional files that are useful for running ribofootPrinter

## Subset lists
Subset lists are MS Excel (.xlsx) files that limit a particular script to operating on a subset of the transcriptome.

subset_list.xlsx - this is an example of a subset gene list of genes with predicted 5'-UTRs.\
subset_list_EIF4G.xlsx - this is a subset list that includes a single gene, EIF4G2.

## caller_scripts_and_input_files
To facilitate input data organization, it is recommended to use some kind caller (i.e. wrapper) functions to run the code. The text files and the .py caller functions here can optionally facilitate this. As an example, the writegene2 function is run this way:
```unix
nohup python writegene2_caller.py writegene2_EIF4G2.txt > 80S_wg2_EIF4G2_metadata.txt 2>&1 &
```
The *nohup* command (and related character string at the end) allows it to run in the background. 

## Parallelization
To parallelize running of ribofootPrinter, bash scripts can be created to include multiple python commands.
