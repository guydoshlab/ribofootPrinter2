import csv
import sys
import re
import tools

### This metagene script counts reads on genes.
# Inputs:
# inputfiles - input rocc files of ribosome footprints
# kind - whether to do average around start codons (1) or stop codons (2)
# weighting - whether to equally weight all genes (1) or to make an unweighted average based on rpm values (0)
# genethresh - the minimal counts (rpkm) in the entire gene to be included in the analysis.
# range5 - the 5' window of the metagene.
# range3 - the 3' window of the metagene.
# subsetlist - excel file with column of genes you want in the average; if no filtering is required put "none".
# outfile - name of the csv file root (no extension) that will be the output of the reads.
# Note this function is currently not compatible with 3' end aligned reads since those require negative shifts, so only "all_5" or "cov"
def main(inputfiles,kind,weighting,genethresh,range5,range3,subsetlist,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	
	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	
	for roccfile in inputfiles:
		samplename=re.split("[./]",roccfile)[-2]
		metagene_values=metagene(roccfile,
									kind,
									weighting,
									genethresh,
									range5,
									range3,
									subsetlist,
									outfile)
		writer.writerow([samplename]+metagene_values)

	writerfile.close()	
	tools.transposecsv(outfile)	


# This is the function that loops through the genes and averages reads.
def metagene(
	roccfile,
	kind,
	weighting,
	genethresh,
	range5,
	range3,
	subsetlist,
	outfile_path
):
	# Load the roccfile.
	rocc_load=tools.roccfile_loader(roccfile)
	footprints=rocc_load[0]
	samplename=rocc_load[1]
	endmode=rocc_load[2]
	if endmode=="all_3":
		print("ERROR - right now the genelist code doesn't support anything except end5 or cov.")
		exit()
	if endmode=="cov" and shift!=0: 
		print("Warning, coverage being used without a shift of 0.")
	
	# kind is either 1 or 2 (1 for start codons, 2 for stop codons; other features can come later).
	kind=int(kind)
	weighting=int(weighting)
	genethresh=int(genethresh)
	range5=int(range5)
	range3=int(range3)
	
	growingmetagene=[0 for x in range(range5+range3)]
	genecount=0

	# Take care of any subsetting.
	footprints_new=tools.subsetter(subsetlist,footprints)		
	footprints=footprints_new		
	
	for gene in footprints.keys():			
		tempmetagene=[0 for x in range(range5+range3)]
		
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		
		# For thresholding and equalweighting, need an rpkm value of the gene.
		# This calculation assumes no end effects on ends of genes (can be done to control for stop/stop peaks) and a shift of 12.
		if (UTR3start-ORFstart)<=0 or ORFstart<12 or (len(footprints[gene][2][endmode])-UTR3start)<12:	# Eliminate genes with no ORF or a UTR so short that reads would not align at ends of ORF.
			continue
		else:
			ORFcounts=(sum(footprints[gene][2][endmode][ORFstart+0-12:UTR3start-0-12]))/((UTR3start-ORFstart-0)/1000)	# This is the rpkm value of the gene.
			if ORFcounts<genethresh:
				continue
		
		if kind==1:
			if ORFstart<range5 or (UTR3start-ORFstart)<range3:	# Check that there is room for ranges.
				continue
			else:
				tempmetagene=footprints[gene][2][endmode][ORFstart-range5:ORFstart+range3]
				genecount+=1
			
		elif kind==2:
			if (len(footprints[gene][2][endmode])-UTR3start)<range3 or (UTR3start-ORFstart)<range5:	# Check that there is room for ranges.
				continue
			else:
				tempmetagene=footprints[gene][2][endmode][UTR3start-range5:UTR3start+range3]
				genecount+=1
		
		else:
			print("error kind")
			exit()
			
		for position in range(range5+range3):
			if weighting==0:
				growingmetagene[position]+=tempmetagene[position]
			else:
				growingmetagene[position]+=(tempmetagene[position]/(ORFcounts/1000))	# Normalize after converting rpkm to rpm of the gene.
			
	for position in range(range5+range3):
		growingmetagene[position]/=genecount
		
	print("Genes included = "+str(genecount))	
	
	return growingmetagene


# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":
	if len(sys.argv)!=9:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])