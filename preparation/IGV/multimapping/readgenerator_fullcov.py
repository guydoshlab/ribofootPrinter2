import sys
from Bio import SeqIO

# This script will divide up the FASTA longnames MANE v1.4 transcriptome into predefined short reads. 
# It will make sure to have a read at every nucleotide (unlike bbmap, which generates random reads).
# FASTQ files are generated.
def main(fasta_in, fastq_out, readlength):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))
	
	# Add in sequences etc. 
	with open (fastq_out, "w") as fastq_out:
		readlength=int(readlength)
		outputdata={} # Data structure for holding information.
		counter={}
		counter["all"]=0
		
		qual_score = "E"*readlength
		
		# Go through the fasta file and populate outputdata. Similar to the rocc file but delete count data entry.
		for record in SeqIO.parse(fasta_in, "fasta"):
			gene=record.id.split("|")[0]
			outputdata[gene]=["alias","sequence","ORFstart","utr3start","transcriptlen"]
			outputdata[gene][0]=record.id.split("|")[5] # alias
			outputdata[gene][1]=str(record.seq) # Sequence of transcript
			outputdata[gene][2]=(record.id.split("|")[7]).split("-")[-1]	# ORF start
			outputdata[gene][3]=(record.id.split("|")[8]).split("-")[-1]	# UTR3 start
			outputdata[gene][4]=record.id.split("|")[6] # Transcript length
			
			genelength=int(outputdata[gene][4])
			seq = str(outputdata[gene][1])

			for position in range(genelength):	# this will loop through the entire length of each transcript
				if position < (genelength - readlength): # only if the end of the transcript is not reached yet. This will depend on the lenght of the readlenght that will be generated.
					read = seq[position:(position+readlength)] # stores each read with fixed length into read variable
					header = str("@SYN_+_") + str(position+1) + str("_") + str(gene)
					fastq_out.write(header + "\n")
					fastq_out.write(read + "\n")
					fastq_out.write("+\n")
					fastq_out.write(qual_score + "\n")
					counter["all"]+=1
					
					if counter["all"] <= 10:	# writes out the first 10 entries. You can check here if the sequence is moving with +1 in each loop.
						print("Line 1:", header)
						print("Line 2:	", read)
						print("Line 3:	", qual_score)
						print("-" * 20)

		print(str(counter["all"])+" reads generated from fasta file that are " + str(readlength) + " nt long")

# Here is where commands are input.			  
if __name__ == "__main__":	
	if len(sys.argv)!=4:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3])