import matplotlib.pyplot as plt # to plot the data
import pandas as pd
import os
import sys

# This function plots 3D metagene csv files. 
# To adjust contrast, use buttons on plot over contrast bar, or manually change below.
# The zoom and panning can be adjusted on the plot using the buttons.
# To adjust labels or manually zoom plot, adjust below.
def main(csv_in):
	print("\nName of python script:",(__file__.split("/")[-1]))
	print("Total arguments passed:", len(locals()))
	print("Argument names: "+("; ".join(list(locals().keys()))))
	argkeys=list(locals().keys())
	argvals=[]
	for key in argkeys:
		argvals.append((locals()[key]))
	print("Argument values: "+("; ".join(argvals)))   
	
	df = pd.read_table(csv_in,sep=",")

	# Set x and y axis 
	y = df["footprint size (nt)"].tolist()
	#x = list(range(1, df.shape[1]))  
	x = df.columns[1:]
	array = df.iloc[:, 1:].values

	# Plot heatmap
	fig, ax = plt.subplots()
	im = ax.imshow(array, aspect='auto',cmap='plasma')
	#im = ax.imshow(array, aspect='auto', vmin=0, vmax=200)	# Change setting to set initial contrast rather than auto-adjustable default.

	# Set ticks for low footprint range, for example 25-34
	yticks = list(range(len(y)))

	# Can adjust here to add more digits on Y axis.
	ax.set_yticks(yticks)
	labels=["" for i in yticks]
	#for i in range(0,len(y)):	# Use this for statement instead for max y-axis labels numbers.
	for i in range(0,len(y), 5):	
		labels[i]=y[i] 
	ax.set_yticklabels(labels)

	numxlabels=6		### Set this to number of label numbers on the x-axis. Currently set for 6 with extra ticks in between.
	multfactor=int((len(x)/((numxlabels-1)*2)))
	#multfactor=int((len(x))/10)
	xticks=list(range(len(x)))[::multfactor]
	ax.set_xticks(xticks)
	labelsx=["" for i in xticks]
	for i in list(range(0,int(len(xticks)/2+1))):	
		labelsx[i*2]=str(x[i*2*multfactor])
	ax.set_xticklabels(labelsx)
	
	#manual_xticks = [0, 100, 200, 300, 400]  # Manual x axis option
	#ax.set_xticks(manual_xticks)  

	# Formatting of the graph
	ax.tick_params(axis='x', labelsize=10)
	ax.tick_params(axis='y', labelsize=10)
	ax.invert_yaxis()
	ax.set_title('3D metagene for: '+csv_in, fontsize=10)
	plt.colorbar(im, ax=ax, label='Average reads (rpm). Pan or zoom on colorbar to adjust contrast')
	plt.ylabel("Readlength (nt)")
	plt.xlabel("Distance from read end to start or stop codon (nt)")
	plt.show()


# For running from the command line, the code below pulls in the variables and calls main.
if __name__ == "__main__":	
	if len(sys.argv)!=2:
		print("Wrong number of inputs.")
		exit()
	main(sys.argv[1])
       