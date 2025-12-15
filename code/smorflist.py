import csv
import sys
import tools
from Bio import Seq


### Script to find small ORFs in UTRs.
# Inputs:
# inputfiles - input rocc files of ribosome footprints
# lengththresh - this is the threshold of amino acids. Must be longer than this. Negative is only this.
# shift - amount to shift data before counting. Typically do for P sites.
# smallest - put 1 to keep smallest ORF (last start codon). Put 0 for longest.
# mismatches - this is how many mismatches to allow in start codon - 1 or 0. For finding near cognates.
# UTR - this is whether to do 5'UTR or 3'UTR (5 or 3).
# outfile - name of the csv file root (no extension) that will be the output of the reads.
# Note this function is currently not compatible with 3' end aligned reads since those require negative shifts, so only "all_5" or "cov"
def main(inputfiles,lengththresh,shift,smallest,mismatches,UTR,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	firstwrite=0
	smallest=int(smallest)
	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.
	for roccfile in inputfiles:
		dorfdict=smorflist(roccfile,
							shift,
							lengththresh,
							mismatches,
							UTR)		
		
		if firstwrite==0:
			growinggenedict=dorfdict
			firstwrite=1
			continue
		else:
			for key in dorfdict.keys():
				growinggenedict[key]+=dorfdict[key][9:]
		
	tossedgenes=0	
	startkey="startkey"	
	towrite="towrite"
	for key in growinggenedict.keys():
		tossedgenescounter=0
		# 0s for all dorfs are gone.
		for i in range(len(inputfiles)):
			if growinggenedict[key][9+i*2]!=0:
				tossedgenescounter=1
				break
		if tossedgenescounter==0:	
			tossedgenes+=1
			continue
		
		
		# CODE TO TAKE LONGEST ISOFORM ONLY. 
		keyrootnum=key.find("_smorf")
		testkey=key[0:keyrootnum]
		if testkey!=startkey:		# New gene (new key).
			startkey=testkey
			currentstop={}
			if towrite!="towrite":
				writer.writerow(towrite)		# Write out old stop codon.
			
			currentstop[growinggenedict[key][8]]=0		# 0 is just a placeholder.
			towrite=[key]+growinggenedict[key]		
			continue				
		if growinggenedict[key][8] in currentstop:	# Check if we are on a new stop codon. If yes, check if length is better and save it.
			if smallest==1:		# Check if we are looking for smaller or longer ORFs. 
				if growinggenedict[key][7] > towrite[8]:		# > signals shortest
					towrite=[key]+growinggenedict[key]
				elif growinggenedict[key][7] == towrite[8]:	#Same length should not occur (not possible).
					print("Error")
					exit()
			else:
				if growinggenedict[key][7] < towrite[8]:		# < signals longest
					towrite=[key]+growinggenedict[key]	
				elif growinggenedict[key][7] == towrite[8]: #Same length should not occur (not possible).
					print("Error")
					exit()
	
		else:	
			writer.writerow(towrite)					# New stop codon, get ready to write out the best option and reset.
			currentstop[growinggenedict[key][8]]=0	
			towrite=[key]+growinggenedict[key]
			
	writer.writerow(towrite)							# Write out very last one.



	writerfile.close()
	print("Genes tossed since no reads in any sample.")
	print(tossedgenes)
	
	

def smorflist(
	roccfile,
	shift,
	lengththresh,
	mismatches,
	UTR
):

	if int(lengththresh)<0:		# Allow a neg length to be used absolutely.
		abslengththresh=int(lengththresh)*-1
		lengththresh=0
	
	# Load the roccfile.
	rocc_load=tools.roccfile_loader(roccfile)
	footprints=rocc_load[0]
	samplename=rocc_load[1]
	endmode=rocc_load[2]
	if endmode=="all_3":
		print("ERROR - right now the genelist code doesn't support anything except end5 or cov.")
		exit()

	shift=int(shift)
	genelist={}
	######
	genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","smallorfseq","smallorfaa","sm_orfstartinUTR","sm_orfstopinUTR","sm_orf_"+samplename,"CDS_"+samplename]

	dorfcounttotal=0
	genecount=0
	for gene in footprints.keys():		
		dorfcount=0
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		
		
		if (ORFstart-shift)<0:	# Check for negative indexes:
				continue
		if (UTR3start-shift)<0:	# Check for negative indexes:
				continue
				
		if int(UTR)==3:
			
			counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
			genesequence=footprints[gene][1][UTR3start:]
		else:	
			counts=footprints[gene][2][endmode][0:ORFstart-shift]
			genesequence=footprints[gene][1][shift:ORFstart]
		
		
		
		genelen=len(counts)

		i=0
		dorfcount=0		# Counting number of dorfs per gene.
		while (i+3)<genelen: 
			testmotif=str(genesequence[i:i+3])
			score=0
			if testmotif[0]!="A":
				score+=1
			if testmotif[1]!="T":
				score+=1
			if testmotif[2]!="G":
				score+=1
			if score>int(mismatches):		# How many mismatches to tolerate.
				i+=1
				continue
				
			# Find in frame stop codon.
				
			j=i+3
				
			while (j+3)<genelen:
				testmotif2=str(genesequence[j:j+3])
				if testmotif2=="TAA" or testmotif2=="TAG" or testmotif2=="TGA":
					endpos=j
					break
				j+=3
			else:
				i+=1
				continue # NO in frame stop.
				
			dorfseq=genesequence[i:j+3]
			dorfaa=Seq.Seq(dorfseq[:-3]).translate()
			dorflen=len(dorfaa)
			if dorflen<int(lengththresh) or (lengththresh==0 and dorflen!=abslengththresh):
				i+=1
				continue
				
			if ((ORFstart-shift)<=0):
				CDScount=-10
			else:
				CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
				
			# Note we don't shift since this is already shifted.
			dorfreads=sum(counts[i:j+3][0::3])/(((dorflen+1)*3)/1000)
						
						
						
			dorfkey=gene+"_smorf_"+str(dorfcount)
		
			
			genelist[dorfkey]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),dorfseq,dorfaa,i,j,0,0]
			genelist[dorfkey][9]=dorfreads
			genelist[dorfkey][10]=CDScount
			dorfcount+=1	
			dorfcounttotal+=1
			i+=1
	
		genecount+=1
		if genecount%1000==0:
			print("Genes finished so far="+str(genecount))		
	print("Total smorfs reported = "+str(dorfcounttotal))	
	return genelist

# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":	
	if len(sys.argv)!=8:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

