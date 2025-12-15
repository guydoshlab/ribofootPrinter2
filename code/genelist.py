import csv
import sys
import tools

# This script counts reads on genes.
# Inputs:
# inputfiles - input rocc files of ribosome footprints
# shift - amount to shift data before counting. Typically do for P sites.
# doextra - Usually set to 0. Setting this variable to 1 will output an extra file with CDS reads broken down by frame. Other extras could be added in future. 
# outfile - name of the csv file root (no extension) that will be the output of the reads.
# Note the dictionaries for gene and frame and pause are called "list" but they are dictionaries.
# Note that counts can only computed when the gene has both UTRs longer than the shift value. Otherwise, genes are left out of the output.
# Note this function is currently not compatible with 3' end aligned reads since those require negative shifts.
# endmode currently supported are: "all_5" or "cov"
def main(inputfiles,shift,doextra,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	
	firstwrite=0	
	doextra=int(doextra)
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	if doextra==1:
		writerfile2=open(outfile+"_frame.csv", "w")	
		writer2 = csv.writer(writerfile2)
	
	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.	
	for roccfile in inputfiles:
		outdictionaries=genelist(roccfile,
				shift)
		
		genedict=outdictionaries[0]
		if doextra==1:
			framedict=outdictionaries[1]
		else:
			framedict={}		# Just a placeholder.
		
		if firstwrite==0:
			growinggenedict=genedict
			growingframedict=framedict
			firstwrite=1
			continue
		else:
			for key in genedict.keys():
				growinggenedict[key]+=genedict[key][6:]
			for key in framedict.keys():
				growingframedict[key]+=framedict[key][1:]
	
	for key in genedict.keys():
		writer.writerow([key]+growinggenedict[key])					
	for key in framedict.keys():
		writer2.writerow([key]+growingframedict[key])		
	
	writerfile.close()		
	if doextra==1:
		writerfile2.close()
		
# This is the function that loops through the genes and counts reads.
def genelist(
	roccfile,
	shift
):
	shift=int(shift)
	# Load the roccfile.
	rocc_load=tools.roccfile_loader(roccfile)
	footprints=rocc_load[0]
	samplename=rocc_load[1]
	endmode=rocc_load[2]
	mappedreads=rocc_load[3]
	if endmode=="all_3":
		print("ERROR - right now the genelist code doesn't support anything except end5 or cov.")
		exit()
	if endmode=="cov" and shift!=0: 
		print("Warning, coverage being used without a shift of 0.")
	
	genecount=0
	tooshort=0
	genelist={}
	genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","seq","UTR5_"+samplename,"CDS_"+samplename,"UTR3_"+samplename,"UTR5raw_"+samplename,"CDSraw_"+samplename,"UTR3raw_"+samplename]
	framelist={}
	framelist["headers"]=["alias","CDS0_"+samplename,"CDS1_"+samplename,"CDS2_"+samplename,"CDS0raw_"+samplename,"CDS1raw_"+samplename,"CDS2raw_"+samplename]
			
	print("To compute raw reads, using total mapped reads = "+str(mappedreads))
	
	for gene in footprints.keys():
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])

		if ((ORFstart-shift)<=0) or (len(footprints[gene][2][endmode])-UTR3start)<=0:		# Either UTR is too short to be compatible with the shift value.
			CDScount=float('nan')
			UTR5count=float('nan')
			UTR3count=float('nan')
			CDSraw=float('nan')
			UTR5raw=float('nan')
			UTR3raw=float('nan')
			CDScount_0=float('nan')
			CDScount_1=float('nan')
			CDScount_2=float('nan')
			CDScount_0raw=float('nan')
			CDScount_1raw=float('nan')
			CDScount_2raw=float('nan')
			tooshort+=1
			continue	# This continue can be deactivated if preferred to keep the genes in the list. Doing so would facilitate comparison of different transcriptomes and make the output independent of shift. 
			
		else:
			CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift])/((UTR3start-ORFstart)/1000)
			UTR5count=sum(footprints[gene][2][endmode][0:ORFstart-shift])/((ORFstart-shift)/1000)
			UTR3count=sum(footprints[gene][2][endmode][UTR3start-shift:-shift])/((len(footprints[gene][2][endmode])-UTR3start)/1000)
			CDSraw=int(round((CDScount*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
			UTR5raw=int(round((UTR5count*mappedreads/1E6)*((ORFstart-shift)/1000)))
			UTR3raw=int(round((UTR3count*mappedreads/1E6)*((len(footprints[gene][2][endmode])-UTR3start)/1000)))
			
			CDScount_0=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
			CDScount_1=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][1::3])/((UTR3start-ORFstart)/1000)
			CDScount_2=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][2::3])/((UTR3start-ORFstart)/1000)
			CDScount_0raw=int(round((CDScount_0*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
			CDScount_1raw=int(round((CDScount_1*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
			CDScount_2raw=int(round((CDScount_2*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
			
		genecount+=1
	
		genelist[gene]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),footprints[gene][1],0,0,0,0,0,0]
		genelist[gene][6]=UTR5count
		genelist[gene][7]=CDScount
		genelist[gene][8]=UTR3count
		genelist[gene][9]=UTR5raw
		genelist[gene][10]=CDSraw
		genelist[gene][11]=UTR3raw
		
		if genelist[gene][4]>32767:
			genelist[gene][5]="sequence too long for Excel - lookup in fasta instead"
		
		# Output frame information		
		framelist[gene]=[footprints[gene][0],CDScount_0,CDScount_1,CDScount_2,CDScount_0raw,CDScount_1raw,CDScount_2raw]
	
	print("Genes included = "+str(genecount))	
	print("Genes dropped due to short UTR = "+str(tooshort))	
	return [genelist,framelist]
	
# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":
	if len(sys.argv)!=5:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

