import csv
import sys
import tools # This is the accompanying tools file with helper functions.

#### This writegene2 script writes out reads on genes to a csv file to be used, for example, to make figures.
# Inputs:
# list of genes to be output (can be alias or formal ENSG... name)
# inputfiles - input rocc files of ribosome footprints
# outfile - name of the csv file root (no extension) that will be the output of the reads.
# endmode supported are: "all_3" or "all_5" or "cov"
def main(inputfiles,list_of_genes,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))

	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.
	print("\nLoading first roccfile for alias conversion.")
	genelist = tools.alias_converter(list_of_genes, inputfiles[0])	# Convert names to a python list. Also convert any aliases to proper names using first inputfile for lookup.
	print("\nLoading roccfiles for output of read data.")
	writegene2(genelist,  		# genename or alias 
				inputfiles, 	# rocc files
				outfile)	# output file name    				

# This is the function that loops through the genes and counts reads.
def writegene2(
	genenames,
	roccfiles,
	outfile
):
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)

	for roccfile in roccfiles:
		# Load the roccfile.
		rocc_load=tools.roccfile_loader(roccfile)
		footprints=rocc_load[0]
		samplename=rocc_load[1]
		endmode=rocc_load[2]
	
		for gene in genenames:
			footprintcounts=list(footprints[gene][2][endmode])
			alias=footprints[gene][0]
			writer.writerow([samplename+"_"+alias+"_"+gene]+footprintcounts)
	
	writerfile.close()	
	tools.transposecsv(outfile)	
	


# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":	
	if len(sys.argv)!=4:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3])
