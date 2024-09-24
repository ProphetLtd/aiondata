import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot(df,column_nm,x,y,title,wcolor):
	# Plot the distribution
	# remove lines with nan
	df = df.dropna(subset=[column_nm])
	# turn to numeric
	df[column_nm] = pd.to_numeric(df[column_nm], errors='coerce')
	plt.hist(df[column_nm].values ,bins=20,edgecolor='black', linewidth=1.2,color=wcolor)  # Adjust the number of bins as needed
	plt.xlabel(x)
	plt.ylabel(y)
	plt.title(title)
	plt.show()




	#plt.xlabel(f'Ave codnas-protein recall (RMSD) ')
	#plt.ylabel('Frequency')
	#if df['pockets_covered'][0]==True:
#		cov=''
#	else:
#		cov='not ' 
#	plt.title(f'Recall for codnas proteins {cov} covered by AF ')
