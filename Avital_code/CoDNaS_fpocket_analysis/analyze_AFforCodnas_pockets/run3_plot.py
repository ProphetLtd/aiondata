import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot(input_csv,xaxis_nm,yaxis_nm,title):
	# Load the CSV file into a DataFrame
	df = pd.read_csv(input_csv, sep=',')
	# Plot the distribution
	bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	hist, bin_edges = np.histogram(df[df.columns[4]], bins=bins)
	plt.bar(bin_edges[:-1], hist, width=0.1,edgecolor='black', linewidth=1.2, align='edge')
#	plt.hist(df[df.columns[4]], bins=[0,0.00001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], edgecolor='black', linewidth=1.2)  # Adjust the number of bins as needed
	# Highlight the bin for zero values
	zero_bin_count = (df[df.columns[4]] == 0).sum()
#	plt.bar(0, zero_bin_count, width=0.1, color='red', edgecolor='black', linewidth=1.2, align='edge', label='0 values')
	# Highlight the 0 values
#	plt.ylim(0, 70)
	plt.xlabel(xaxis_nm)
	plt.ylabel(yaxis_nm)
	plt.title(title )
	plt.show()


plot('codnas_fpockets_compared_to_AFhighest.csv','Coverage of codnas fpocket by AF','# Instances','All pockets with AF model')
plot('codnas_fpockets_compared_to_AFhighest_biolip.csv','Coverage of codnas fpocket by AF','# Instances','HOLO (Pockets with bound ligand in BioLip)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip.csv','Coverage of codnas fpocket by AF','# Instances','APO (Pockets without bound ligand)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_wHoloConf','Coverage of codnas fpocket by AF','# Instances','APO , is HOLO in another confromation')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_woutHoloConf','Coverage of codnas fpocket by AF','# Instances','APO , no bound ligand in any known conformation')

