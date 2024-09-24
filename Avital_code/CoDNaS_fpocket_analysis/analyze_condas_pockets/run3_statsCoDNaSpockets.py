import pandas as pd 
import numpy as np
from plot_scatter_column_df import plot
import matplotlib.pyplot as plt
import seaborn as sns

df=pd.read_csv('codnas_fpockets',index_col=False)


# Count unique values in 'cluster_nm '
count_prots=df['cluster_nm'].nunique()
print ('number of pdbs:   ',count_prots)
# Count unique values in 'pdb_chain' which is  model !
count_confs = df['pdb_chain'].nunique()
print ('number of conformations:   ' ,count_confs)
# For each unique value in pdb model: 'pdb_chan', count unique values in 'pocket'
pockets_per_conf = df.groupby('pdb_chain')['pocket'].nunique()
print ('Average, min, max number of pockets per real conf',round(np.mean(pockets_per_conf),2),min(pockets_per_conf),max(pockets_per_conf))

# Per predicted Conformation
#----------------------
	# max pocket score per protein conformation prediction , pdb_chain= model, conformation
max_pocket_score_per_conf=df.groupby('pdb_chain')['Pscore'].max()
# Convert Series to DataFrame
df_max_pocket_score_per_conf = max_pocket_score_per_conf.reset_index()
# Rename columns for clarity
df_max_pocket_score_per_conf.columns = ['pdb_chain', 'Pscore']
mean_of_max_pocket_score_per_conf=df_max_pocket_score_per_conf['Pscore'].mean()
print('mean of max pocket_score per conformation ',mean_of_max_pocket_score_per_conf)
plot(df_max_pocket_score_per_conf,'Pscore','max pocket score per conformation', 'freq','distribution of highest pocket score per conformation ','darkblue')
	# max Dscore per conf
max_dpocket_score_per_conf=df.groupby('pdb_chain')['Dscore'].max()
# Convert Series to DataFrame
df_max_dpocket_score_per_conf = max_dpocket_score_per_conf.reset_index()
# Rename columns for clarity
df_max_dpocket_score_per_conf.columns = ['pdb_chain', 'Dscore']
mean_of_max_dpocket_score_per_conf=df_max_dpocket_score_per_conf['Dscore'].mean()
print('mean of max dpocket_score per conformation ',mean_of_max_dpocket_score_per_conf)
plot(df_max_dpocket_score_per_conf,'Dscore','max drug pocket score per conformation ', 'freq','distribution of highest drug-pocket score per conformation ','blue')
	# mean of Pscore per conformaiton
mean_pocket_score_per_conf=df.groupby('pdb_chain')['Pscore'].mean()
# Convert Series to DataFrame
df_mean_pocket_score_per_conf = mean_pocket_score_per_conf.reset_index()
# Rename columns for clarity
df_mean_pocket_score_per_conf.columns = ['pdb_chain', 'Pscore']
mean_of_mean_pocket_score_per_conf=df_mean_pocket_score_per_conf['Pscore'].mean()
print ('mean of mean pocket_score per conformation' , mean_of_mean_pocket_score_per_conf)
#plot(df_mean_pocket_score_per_conf,'Pscore','mean pocket score of conf', 'freq','distribution of highest pocket score per conformation ','darkblue')
	# mean of Dscore per conf
mean_dpocket_score_per_conf=df.groupby('pdb_chain')['Dscore'].mean()
# Convert Series to DataFrame
df_mean_dpocket_score_per_conf = mean_dpocket_score_per_conf.reset_index()
# Rename columns for clarity
df_mean_dpocket_score_per_conf.columns = ['pdb_chain', 'Dscore']
mean_of_mean_dpocket_score_per_conf=df_mean_dpocket_score_per_conf['Dscore'].mean()
print ('mean of mean dpocket_score per conformation' , mean_of_mean_dpocket_score_per_conf)
#plot(df_mean_dpocket_score_per_conf,'Dscore','mean drug pocket score of conf', 'freq','distribution of highest pocket score per conformation ','blue')

# Per protein!
#------------------
	# max pocket score per protein for all its predictions
max_pocket_score_per_prot=df.groupby('cluster_nm')['Pscore'].max()
# Convert Series to DataFrame
df_max_pocket_score_per_prot = max_pocket_score_per_prot.reset_index()
# Rename columns for clarity
df_max_pocket_score_per_prot.columns = ['cluster_nm', 'Pscore']
mean_of_max_pocket_score_per_prot=df_max_pocket_score_per_prot['Pscore'].mean()
print ('mean of max pocket score per prot',mean_of_max_pocket_score_per_prot)
plot(df_max_pocket_score_per_prot,'Pscore','max pocket score of prot', 'freq','distribution of highest pocket score per protein ','darkblue')
	# max Drug pocket score per protein for all its predictions
max_dpocket_score_per_prot=df.groupby('cluster_nm')['Dscore'].max()
# Convert Series to DataFrame
df_max_dpocket_score_per_prot = max_dpocket_score_per_prot.reset_index()
# Rename columns for clarity
df_max_dpocket_score_per_prot.columns = ['cluster_nm', 'Dscore']
mean_of_max_dpocket_score_per_prot=df_max_dpocket_score_per_prot['Dscore'].mean()
print ('mean of max drug-pocket score per prot',mean_of_max_dpocket_score_per_prot)
plot(df_max_dpocket_score_per_prot,'Dscore','max drug pocket score of prot', 'freq','distribution of highest drug-pocket score per protein ','darkblue')
	# mean pocket score per protein for all its predictions
mean_pocket_score_per_prot=df.groupby('cluster_nm')['Pscore'].mean()
# Convert Series to DataFrame
df_mean_pocket_score_per_prot = mean_pocket_score_per_prot.reset_index()
# Rename columns for clarity
df_mean_pocket_score_per_prot.columns = ['cluster_nm', 'Pscore']
mean_of_mean_pocket_score_per_prot=df_mean_pocket_score_per_prot['Pscore'].mean()
print ('mean of mean pocket score per prot',mean_of_mean_pocket_score_per_prot)
#plot(df_mean_pocket_score_per_prot,'Pscore','mean pocket score of prot', 'freq','distribution of highest pocket score per protein ','darkblue')
	# mean Drug pocket score per protein for all its predictions
mean_dpocket_score_per_prot=df.groupby('cluster_nm')['Dscore'].mean()
# Convert Series to DataFrame
df_mean_dpocket_score_per_prot = mean_dpocket_score_per_prot.reset_index()
# Rename columns for clarity
df_mean_dpocket_score_per_prot.columns = ['cluster_nm', 'Dscore']
mean_of_mean_dpocket_score_per_prot=df_mean_dpocket_score_per_prot['Dscore'].mean()
print ('mean of mean drug-pocket score per prot',mean_of_mean_dpocket_score_per_prot)
#plot(df_mean_dpocket_score_per_prot,'Dscore','mean drug-pocket score of prot', 'freq','distribution of highest drug-pocket score per protein ','darkblue')
