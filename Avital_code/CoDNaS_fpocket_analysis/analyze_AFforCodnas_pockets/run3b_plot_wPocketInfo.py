import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot(input_csv,xaxis_nm,yaxis_nm,avg_col_name,title):
	# Load the CSV file into a DataFrame
	df = pd.read_csv(input_csv, sep=',')
	# Plot the distribution
	bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	hist, bin_edges = np.histogram(df[df.columns[4]], bins=bins)

	# df of codnas pockets
	df_cod_pockets=pd.read_csv('../../codnas/codnas_Signif_reduced_fpockets',index_col=False)
	df_cod_pockets['pdb_chain_pocket'] = df_cod_pockets['pdb_chain'] + '_'+df_cod_pockets['pocket']
	df_new = df.merge(df_cod_pockets[['pdb_chain_pocket','Pscore','Dscore','Aspheres','Vol']], on='pdb_chain_pocket', how='right')

	hist, bin_edges = np.histogram(df_new[df_new.columns[4]], bins=bins)
	df_new['bin'] = pd.cut(df_new[df_new.columns[4]], bins=bins,right=False)
	df_new.to_csv('try')
	bin_averages = df_new.groupby('bin')[avg_col_name].mean()
	print (bin_edges[:-1], hist)
	print (bin_averages)
	plt.bar(bin_edges[:-1], hist, width=0.1,edgecolor='black', linewidth=1.2, align='edge')
	# Plot the averages as a line plot on the same axis
	#plt.plot(bin_edges[:-1] + 0.05, bin_averages.to_numpy(), color='red', marker='o', label=f'Avg {avg_col_name}')
	for i, avg in enumerate(bin_averages):
		plt.text(bin_edges[i] + 0.05, hist[i] + 0.5, f'{avg:.2f}', ha='center', color='blue')   
		print (i, bin_edges[i], hist[i], avg) 
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


plot('codnas_fpockets_compared_to_AFhighest.csv','Coverage of codnas fpocket by AF','# Instances','Vol','All pockets with AF model')
plot('codnas_fpockets_compared_to_AFhighest_biolip.csv','Coverage of codnas fpocket by AF','# Instances','Vol','HOLO (Pockets with bound ligand in BioLip)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip.csv','Coverage of codnas fpocket by AF','# Instances','Vol','APO (Pockets without bound ligand)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_wHoloConf','Coverage of codnas fpocket by AF','# Instances','Vol','APO , is HOLO in another confromation')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_woutHoloConf','Coverage of codnas fpocket by AF','# Instances','Vol','APO , no bound ligand in any known conformation')

plot('codnas_fpockets_compared_to_AFhighest.csv','Coverage of codnas fpocket by AF','# Instances','Pscore','All pockets with AF model')
plot('codnas_fpockets_compared_to_AFhighest_biolip.csv','Coverage of codnas fpocket by AF','# Instances','Pscore','HOLO (Pockets with bound ligand in BioLip)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip.csv','Coverage of codnas fpocket by AF','# Instances','Pscore','APO (Pockets without bound ligand)')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_wHoloConf','Coverage of codnas fpocket by AF','# Instances','Pscore','APO , is HOLO in another confromation')
plot('codnas_fpockets_compared_to_AFhighest_Nobiolip_woutHoloConf','Coverage of codnas fpocket by AF','# Instances','Pscore','APO , no bound ligand in any known conformation')
