from Bio import SeqIO
import csv
import sys

### This is the function for converting a longnames fasta file to a GTF for IGV.

# fasta_in is the longnames fasta (described elsewhere) with transcriptome information and annotation.
# outfile is the name of the GTF outfile.
# transcriptome is the name of the transcriptome (i.e. MANE_v1.4 or Scer_R64-1-1)
def main(fasta_in,outfile,transcriptome):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals))+"\n")
	
	counter=0
	f = open(outfile, "w")
	
	# Loop through fasta file.
	for record in SeqIO.parse(fasta_in, "fasta"):       # Reading in a fasta file using the SeqIO package.
		gene=record.id.split("|")[0]                    # Obtains genename from fasta file (ENSG000000001.12).
		transcript=record.id.split("|")[1]  			# Transcript name in current convention.
		protein=record.id.split("|")[2] 				# Protein name, fixed from original bug.
		# In future, could add alternative IDs here.
		
		startcodon=1+int((record.id.split("|")[7]).split("-")[-1])
		cdsend=record.id.split("|")[8].split("-")[-1]
		outputline=[]
		outputline.append(gene)
		outputline.append(transcriptome)
		outputline.append("CDS")
		outputline.append(startcodon)	#main ORF start codon
		outputline.append(cdsend)	#UTR3 start
		outputline.append(".")
		outputline.append("+")
		outputline.append("0")
		outputline.append('gene_id "'+gene+'";'+' transcript_id "'+transcript+'";'+' gene_type "protein_coding";'+' gene_name "'+record.id.split("|")[5]+'";'+' transcript_type "protein_coding";'+' protein_id "'+protein+'"')
		f.write("\t".join(map(str, outputline)) + "\n")
		counter+=1
		
	print("Completed.")
	print("Number of genes = "+str(counter))
	
	f.close()


# For running from the command line, the code below pulls in the variables and calls main.
# Alternatively, a separate python caller script can call main. 
if __name__ == "__main__":
	if len(sys.argv)!=4:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1],sys.argv[2],sys.argv[3])

