import csv
import sys
import re
import tools
from Bio import Seq

#### Codon analysis position average analysis (metacodon analysis) - it finds motifs of interest in gene sequences and averages the reads around them.
# inputfiles - input rocc files of ribosome footprints
# outfile - the name of the output csv file (do not include .csv extension)
# motif - the nt of aa of the motif that is being averaged. Can be a comma separated list. Note you can specify the first or last occurence only be included, as noted in the code comments with a "f" or "l" preceding the motif. It allows use of "all_N" to look at every amino acid or nt combination (N is number of nt or aa); in this case pause scores computed from the average data are also written out (-1 for motifs that don't exist in the transcriptome).
# kind - whether the motif is an amino acid (1) or nucleotide sequence (0)
# frame - can be 0, 1, or 2. It's the frame to be looking in with respect to the first base in the region of interest (5'UTR, main ORF, or 3'UTR). Put 3 for all 3 frames.
# bkndwindowthresh - the minimal rpkm counts within the bkndwindow to be included in the average.
# bkndwindow - this is the size of the half-window around the codon of interest to be included in the average.
# ORFnorm - normalize to the ORF if doing UTRs (UTRmode 0 or 2). 0 does not normalize. A positive value will ORF normalize and specify the rpkm ORF threshold for genes to include. This uses a hard-coded shift of +/-12 to compute the main ORF rpkm.
# shift - nt to shift data for site of interest (P site, A site, etc).
# UTRmode - over what region motifs are found and averaged. UTR5,CDS,or UTR3 = 0,1,2, respectively.
# subsetlist - excel file with column of genes you want in the average; if no filtering is required put "none".
# Pause score is computed from the average plot, uses a hardcoded window of +/- 1 nt around the peak.
# There is a hard-coded requirement that if using ORFnorm>0 (UTR regions only), the UTR region background cannot be >10x the ORF (these are cases of a bad annotation where UTR is actually main ORF). 
# endmode supported are: "all_3" or "all_5" or "cov"
def main(inputfiles,motif,kind,frame,bkndwindowthresh,bkndwindow,ORFnorm,shift,UTRmode,subsetlist,outfile):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	inputfiles=inputfiles.split(",")	# Turn inputfiles into a python list.
	kind=int(kind)
	frame=int(frame)
	bkndwindowthresh=int(bkndwindowthresh)
	bkndwindow=int(bkndwindow)
	ORFnorm=float(ORFnorm)
	shift=int(shift)	
	UTRmode=int(UTRmode)
	
	if len(motif)>4 and motif[0:3]=="all":			# This will do all aas or triplets, and compute pause scores in addition to averages.
		allmotiflen=int(motif[4:])
		if (allmotiflen>3 and kind==1) or (allmotiflen>6 and kind==0):
			print("Current limit is 3 for amino acid and 6 for nucleotide motifs. Input argument of "+str(allmotiflen)+" not allowed. Exiting.")
			exit()
		motif="all"
		pausescores=[]		
		filenames=[]
		motifused=[]
		if kind==0:
			motifs_loop=["A","G","C","T"]
		elif kind==1:
			motifs_loop=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] # No stop codons since few cases.
		else:
			exit()
		motifs=[""]
		for ii in range(allmotiflen):
			motifs_old=list(motifs)
			motifs=[]
			for jj in motifs_loop:
				for kk in motifs_old:
					motifs.append(kk+jj)
					
	else:
		motifs=list(motif.split(","))
		
	writerfile=open(outfile+"_avgdata.csv", "w")	
	writer = csv.writer(writerfile)
	
	for loc_motif in motifs:
		for roccfile in inputfiles:
			avggene=posavg(roccfile,
					loc_motif,
					kind,
					frame,
					bkndwindowthresh,
					bkndwindow,
					ORFnorm,
					shift,
					UTRmode,
					subsetlist)

			#Compute pause score
			if motif=="all":	
				num=sum(avggene[bkndwindow-1:bkndwindow+2])/3		# assume a peak of 3 nt.
				denom=sum(avggene)/len(avggene)	
				if denom!=0:
					pausescores.append(num/denom)
				else:
					pausescores.append(-1)	# Case where there were no cases of the motif. Set to -1 so it is clear there is no pause score.
				filenames.append(re.split("[./]",roccfile)[-2])
				motifused.append(loc_motif)
			
			writer.writerow([re.split("[./]",roccfile)[-2]+"_"+loc_motif]+avggene)	# Write avggene.
		
	writerfile.close()
	tools.transposecsv(outfile+"_avgdata")		
	
	# Write out scores. This is only done when all motifs are used.
	if motif=="all":		
		writerfile=open(outfile+"_score.csv", "w")	
		writer = csv.writer(writerfile)
		writer.writerow(filenames)
		writer.writerow(motifused)
		writer.writerow(pausescores)
		writerfile.close()
		tools.transposecsv(outfile+"_score")	


def posavg(
	roccfile,
	motif,
	kind,
	frame,
	bkndwindowthresh,
	bkndwindow,
	ORFnorm,
	shift,
	UTRmode,
	subsetlist
):
	# Load the roccfile.
	rocc_load=tools.roccfile_loader(roccfile)
	footprints=rocc_load[0]
	samplename=rocc_load[1]
	endmode=rocc_load[2]
	
	# Take care of any subsetting.
	footprints_new=tools.subsetter(subsetlist,footprints)		
	footprints=footprints_new	
	
	if ORFnorm>0 and UTRmode==1:
		print("Error: ORFnorm can only be used when doing a UTR (UTRmode 0 or 2).")
		exit()

	firstlast=0		# Set to 1 if doing first only; 2 if last only. Note this will check all frames.
	if motif[0]=="f":
		motif=motif[1:]
		firstlast=1
	elif motif[0]=="l":
		motif=motif[1:]
		firstlast=2
	if firstlast>0 and kind!=0:
		print("Error. Use of first and last for a motif is valid for nucleotide only.")
		exit()

	count=0
	i=0
	motiflen=len(motif)
	averagegene=[0 for x in range(2*bkndwindow)]

	if frame==3:
		framelist=[0,1,2]
	else:
		framelist=[frame]
			
	for frame in framelist:
		genecount=0
		print("Frame = "+str(frame)+". Motif = "+str(motif)+".")

		for gene in footprints.keys():
			#if genecount%9000==0:
			#	print("Genes finished so far="+str(genecount))
			genecount+=1
			ORFstart=int(footprints[gene][3])
			UTR3start=int(footprints[gene][4])
		
			# For positive shifts:			
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
					if shift==0:
						counts=footprints[gene][2][endmode][UTR3start-shift:]	# Can't slice to a negative 0.
					else:
						counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
					genesequence=footprints[gene][1][UTR3start:]
					if (UTR3start-shift)<0:	# Check for negative indexes:
						continue
		
			# FOR negative SHIFTS (3' end aligned):
			if shift<0:
				if UTRmode==0:
					counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
					genesequence=footprints[gene][1][0:ORFstart]
					if (ORFstart-shift)>len(footprints[gene][2][endmode]):		# Check for indexes off end:
						continue
				elif UTRmode==1:
					counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
					genesequence=footprints[gene][1][ORFstart:UTR3start]
					if (ORFstart-shift)>len(footprints[gene][2][endmode]) or (UTR3start-shift)>len(footprints[gene][2][endmode]):	# Check for indexes off end:
						continue
				elif UTRmode==2:
					counts=footprints[gene][2][endmode][UTR3start-shift:]
					genesequence=footprints[gene][1][UTR3start:shift]
					if (UTR3start-shift)>len(footprints[gene][2][endmode]):	# Check for indexes off end:
						continue
			
			if len(counts)==0:
				#Skip gene because shift isn't allowing us to get the full region.
				continue
			
			# Set shift_loc which will be used for estimating the rpkm value of the main ORF.
			if shift<0:
				shift_loc=-12
			else:
				shift_loc=12
			ORFlevel=(sum(footprints[gene][2][endmode][ORFstart-shift_loc:UTR3start-shift_loc]))/(UTR3start-ORFstart) #avg rpm per len level						
		
			# Set a value that helps to eliminate poorly annotated UTRs (so UTR includes some main ORF reads).
			# Set the level at which ORFnorm no longer includes a 3'UTR.
			badUTRcheck=(10*ORFlevel)
			#badUTRcheck=1000000000000000000   # Essentially skips this.
		
			# Take out genes if ORFnorm is being used and you're reading off ends of genes due to short UTR, or if the ORFnorm threshold is not reached:
			if ORFnorm>0:
				if shift<0:
					if (ORFstart-shift_loc)>len(counts) or (UTR3start-shift_loc)>len(counts):
						continue
				if shift>=0:
					if (ORFstart-shift_loc)<0 or (UTR3start-shift_loc)<0:	# Check for negative indexes:
						continue
				if ((ORFlevel*1000)<ORFnorm): # rpkm conversion and check that ORF is over ORF thresh.
					continue
		
			genelen=len(genesequence)

			
			### Find motifs of interest ###
			i=frame
			
			if kind==0:
				while (i+motiflen)<genelen: 
					if motif==str(genesequence[i:i+motiflen]):
						
						if firstlast==1:	# Check if another instance before this starting just past first codon.
							if motif in str(genesequence[1:i]):
								i+=motiflen
								continue
						
						
						if firstlast==2:	# CHeck if another instance of motif is coming in any frame.
							if motif in str(genesequence[i+1:]):
								i+=motiflen
								continue
					
						if (i+bkndwindow)<genelen and (i-bkndwindow)>0:
							bknd=counts[i-bkndwindow:i+bkndwindow]
							
							if ((sum(bknd))/(len(bknd)/1000))<bkndwindowthresh:
								i+=motiflen
								continue
								
							if ORFnorm==0:
								normfactor=sum(bknd)/(len(bknd))	#rpm per length units.
							else:
								normfactor=ORFlevel
								if (sum(bknd))/(len(bknd))>=badUTRcheck:	# Check for bad UTR annotation
									i+=motiflen
								#	print(gene) # For debugging.
									continue

							if normfactor>0:		# Note this will allow 0-read genes when ORFnorm used but not when ORFnorm=0.
								for j in range(len(averagegene)):
									averagegene[j]+=((counts[i-bkndwindow+j])/normfactor)	
								count+=1
								# For debugging
								#outputsite=[gene]+[i]+averagegene
								#writer.writerow(outputsite)
							else:
								i+=motiflen
								continue
						else:
							i+=motiflen
							continue
					else:
						i+=motiflen
						continue	
					i+=motiflen

			elif kind==1:
				i=0
				counts=counts[frame:]		# Note that this loss here of a base or two on the end eliminates around 100 genes. 
				genesequence=genesequence[frame:]
				genelen_seq=len(genesequence)
				trimlen=genelen_seq%3
				if trimlen!=0:
					genesequence=genesequence[0:-trimlen]
				genesequence=Seq.Seq(genesequence).translate()
				genelen_seq=len(genesequence)
				genelen=len(counts)
	
				while (i+motiflen)<genelen_seq: 
					if motif==str(genesequence[i:i+motiflen]):
					
						if (i*3+bkndwindow)<genelen and (i*3-bkndwindow)>0:
							bknd=counts[i*3-bkndwindow:i*3+bkndwindow]
						
							if ((sum(bknd))/(len(bknd)/1000))<bkndwindowthresh:		
								i+=motiflen
								continue
						
							if ORFnorm==0:
								normfactor=sum(bknd)/(len(bknd))			#rpm per length units.
							else:
								normfactor=ORFlevel
								if (sum(bknd))/(len(bknd))>=badUTRcheck:	# Check for bad UTR annotation
									i+=motiflen
									continue
								
							if normfactor>0:	# Note this will allow 0-read genes when ORFnorm used but not when ORFnorm=0.
						
								for j in range(len(averagegene)):
									averagegene[j]+=((counts[i*3-bkndwindow+j])/normfactor)
								count+=1
							else:
								i+=motiflen
								continue

						else:
							i+=motiflen
							continue
					else:
						i+=motiflen
						continue	
					i+=motiflen
		

	#if count==0:
		#print("Error 0 count")
		#quit()

	#Normalize and deal with motifs that are not in genome:
	for i in range(len(averagegene)):
		if count==0:
			averagegene[i]=0
		else:
			averagegene[i]/=count

		
	print("Number of positions in average")
	print(count)
	return(averagegene)

# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can parse a txt input file and call main. 
if __name__ == "__main__":	
	if len(sys.argv)!=12:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11])