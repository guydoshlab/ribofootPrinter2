import csv
import tools
import sys
import re
from Bio import Seq

### Codon analysis posstats script.
# Inputs:
# inputfiles - input rocc file of ribosome footprints
# motif - aa or nt sequence to compute scores at
# kind is 0 for nt or 1 for aa.
# frame is 0, 1, or 2 to be looking for motifs in.
# genethresh - this is the minimal counts (rpkm) in the region of interest (pause score denominator) to be included in the analysis.
# pkwindow - how much on either side of peak to use in numerator of pause score.
# shift is a positive value for 5' end aligned, negative for 3' end aligned. It is amount to shift data before computing score (P site, A site, etc).
# UTRmode is whether doing UTR5,CDS,or UTR3 = 0,1,2, respectively. Note that the denominator of the pause score is computed from total reads in this region of interest.
# outfile is the name of the csv file for output.
# endmode supported are: "all_3" or "all_5" or "cov"
def main(inputfiles,motif,kind,frame,genethresh,pkwindow,shift,UTRmode,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.
	
	for roccfile in inputfiles:
		outfiletemp=outfile+"_"+re.split("[./]",roccfile)[-2]
		posstats(roccfile,
						motif,
						kind,
						frame,
						genethresh,
						pkwindow,
						shift,
						UTRmode,
						outfiletemp)

def posstats(
	roccfile,
	motif,
	kind,
	frame,
	genethresh,
	pkwindow,
	shift,
	UTRmode,
	outfile
):
	kind=int(kind)
	frame=int(frame)
	genethresh=int(genethresh)
	pkwindow=int(pkwindow)
	shift=int(shift)	
	UTRmode=int(UTRmode)
	
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
	
	i=0
	motiflen=len(motif)
	writerfile=open(outfile+".csv", "w")
	writer = csv.writer(writerfile)
	
	outputsite=["headers","alias","mRNA_position","localsequence_nt","localsequence_aa","numerator","denominator","pausescore"]
	writer.writerow(outputsite)
	for gene in footprints.keys():
		
		alias=footprints[gene][0]
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		

		# For positive shifts:			#Note 0 shifts won't work.
		if shift>=0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][0:ORFstart-shift]
				genesequence=footprints[gene][1][shift:ORFstart]
				if (ORFstart-shift)<0:		# Check for negative indexes: 
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)<0 or (UTR3start-shift)<0:	# Check for negative indexes:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
				genesequence=footprints[gene][1][UTR3start:]
				if (UTR3start-shift)<0:	# Check for negative indexes:
					continue
		
		# FOR negative SHIFTS (3' end aligned):
		if shift<0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
				genesequence=footprints[gene][1][0:ORFstart]
				if (ORFstart-shift)>len(counts):		# Check for indexes off end:
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)>len(counts) or (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:]
				genesequence=footprints[gene][1][UTR3start:shift]
				if (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
					
					
					
		
		genelen=len(genesequence)
		# Find motifs of interest
		i=frame
		
		if kind==0:
			while (i+motiflen)<genelen: 

				if motif==str(genesequence[i:i+motiflen]):
					if (i+pkwindow)<genelen:
						numerator=sum(counts[i-pkwindow:i+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)

						
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1

					else:
						i+=motiflen
						continue
					if (i-9)>=0 and (i+motiflen+9)<genelen:
						localsequence=genesequence[i-9:i+motiflen+9]
						localsequenceaa=Seq.Seq(genesequence[i-9:i+motiflen+9-motiflen%3]).translate()
					else:
						localsequence="none_tooclose"
						localsequenceaa="none_tooclose"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue	
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen
		
		elif kind==1:
			i=0
			counts=counts[frame:]
			genesequencent=genesequence[frame:]
			frameadj=(len(genesequence[frame:]))%3
			templen=len(genesequence)
			genesequence=Seq.Seq(genesequence[frame:templen-frameadj]).translate()
			genelen_seq=len(genesequence)
			genelen=len(counts)
	
			while (i+motiflen)<genelen_seq: 
				if motif==str(genesequence[i:i+motiflen]):
					
					if (i*3+pkwindow)<genelen:
						numerator=sum(counts[i*3-pkwindow:i*3+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1
					
					else:
						i+=motiflen
						continue
						
					if (i-9)>=0 and (i+motiflen+9)<genelen_seq:
						localsequenceaa=genesequence[i-9:i+motiflen+9]
						localsequence=genesequencent[i*3-18:i*3+motiflen*3+18]
												
						
					else:
						localsequenceaa="none_tooclose"
						localsequence="none_tooclose"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue
					
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen

# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":	
	if len(sys.argv)!=10:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9])
