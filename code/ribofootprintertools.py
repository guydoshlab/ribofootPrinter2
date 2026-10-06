import csv
import os
import re
import pickle
import gzip
import pandas as pd

# This is a collection of helper functions for the ribofootprinter package.

# Transposes a csv file by reading it in, transposing, writing out, deleting the original, and renaming the transposed version.
def transposecsv(csvfile):
	writerfile=open(csvfile+"_transposed.csv", "w")
	writer = csv.writer(writerfile,delimiter=',')
	f=open(csvfile+".csv")
	readercsv=csv.reader(f)
	maxrow=0
	for row in readercsv:
		if len(row)>maxrow:
			maxrow=len(row)
	
	columnnum=maxrow
	for column in range(columnnum):
		f.close() 
		f=open(csvfile+".csv")
		readercsv=csv.reader(f)
		col=extractcolumn(readercsv,column)
		writer.writerow(col)
	f.close()
	writerfile.close()
	os.remove(csvfile+".csv")
	os.rename(csvfile+"_transposed.csv",csvfile+".csv")
	
# Extract the given column from a csv file.				
def extractcolumn(rowgen,colnum):
	column=[]
	for row in rowgen:
		if len(row)>colnum:
			column.append(row[colnum])
		else:
			#column.append(float('nan'))
			column.append("")	# Using this since the option above does not result in NaN in Excel (only in R). This gives NA in R and nothing in Excel.
	return column

# Function that takes in a list of proper genenames and converts any aliases (gene names) to proper IDs. Needs a roccfile to use for the conversions.
def alias_converter(genenames,roccfile):
	
	genelist=list(genenames.split(","))
	
	# Handle Alias conversion:
	aliasdict={}	# Dictionary of all gene names, keyed by formal ENSG gene names.
	rocc_load=roccfile_loader(roccfile)
	footprints=rocc_load[0]
	for gene in footprints.keys():
		aliasdict[footprints[gene][0]]=gene
	for gnnum in range(len(genelist)):
		if genelist[gnnum][0:4]!="ENSG":	# Not the ENSG name, so need to convert.
			alias=genelist[gnnum]	# Name of the alternative name in question.
			if alias in aliasdict:	# Look for this alternative name in the alias dictionary.
				genelist[gnnum]=aliasdict[alias]
			else:
				print("Error - alias "+alias+" does not exist. Try the ENSG... format name.")
				exit()
	return genelist

# Subsetter takes a subset list (path to excel file) and eliminates genes from a footprints data structure.
def subsetter(subsetlist,footprints):
	gene_list = None
	if subsetlist != "none":
		selection = pd.read_excel(subsetlist) # Provide list of genenames to subset dataset
		gene_list = set(selection["geneid"].tolist())	# Take the column with "genename" as a header and make a list from it.
		print ("Subset list was loaded containing " + str(len(gene_list)) + " genes.")
	
	if gene_list is not None:
		footprints_new={}	# New footprints that will lack excluded genes.
		for gene in footprints.keys():
			if gene not in gene_list:  # Filter out any genes that are not in the subset list
				continue
			else:
				footprints_new[gene]=footprints[gene]
		return footprints_new			
	else:
		return footprints

# This loads roccfiles and returns a tuple of the footprint data structure (dictionary), samplename, endmode, and mapped reads.
# It also removes the metadata key.
def roccfile_loader(roccfile):
	samplename=re.split("[./]",roccfile)[-2]
	
	f=gzip.open(roccfile,"rb")
	footprints=pickle.load(f)
	print("\nData loaded for "+samplename)
	f.close()
	
	# Remove and extract metadata dictionary. Provide defaults if not found.
	metadata=footprints.pop("info",{"endmode":["all_5","all_3"],"mappreads":0})
	
	if len(metadata["endmode"])!=1:
		print("Warning: more than one endmode detected in rocc file.")
	endmode=metadata["endmode"][0]
	print("Using endmode = "+endmode)
	
	mappedreads=metadata["mappedreads"]
	
	return (footprints,samplename,endmode,mappedreads)