from collect_pocket import pocket_overlap 
import pandas as pd
from align_uniprot_res_to_pdb import renumber_uniprot_residues
# read AF pockets into df
df_AFpockets=pd.read_csv('AF_fpockets')
# read codnas pockets into df
df_codnaSpockets=pd.read_csv('../analyze_condas_pockets/codnas_Signif_reduced_fpockets')
# read df of sifts to renumber uniprot pocket residues to match the predicted pdb
df_sifts = pd.read_csv('/data/sifts/pdb_chain_uniprot.csv',skiprows=1)
# create a df_codnas_fpockets_compared_to_AF
df_codnas_fpockets_compared_to_AF=pd.DataFrame(columns=['cluster_nm', 'pdb_chain_pocket', 'AF_model', 'AF_model_pocket', 'pocket_coverage', 'AF_coverage', 'pockets_covered'])
# go over all codnas pockets and see if they are modeled here:
# Function to apply to each row
def process_row(row):
	global df_codnas_fpockets_compared_to_AF
	cluster_nm,pdb_chain,pocket,p_residues=row['cluster_nm'],row['pdb_chain'],row['pocket'],row['residues']
	# reduce df_sifts only to the line of pdb_chain:
	pdb,chain=pdb_chain.split('_')[0].split('-')[0].lower(),pdb_chain.split('_')[1]
	line_sifts=df_sifts[(df_sifts['PDB'] == pdb) & (df_sifts['CHAIN'] == chain) ].iloc[0]
	df_line_sifts = pd.DataFrame(line_sifts).transpose()

	p_residues=p_residues.replace("'","").replace("[","").replace("]","").split(',')
	# Filter the DataFrame
	filtered_df_AFpockets = df_AFpockets[df_AFpockets['cluster_nm'] == cluster_nm]
	residues_AF=filtered_df_AFpockets['residues']
	
	for af_p_residues in residues_AF.values:
		AF_model,AF_model_pocket=filtered_df_AFpockets[filtered_df_AFpockets['residues']==af_p_residues]['pdb_chain'].iloc[0],filtered_df_AFpockets[filtered_df_AFpockets['residues']==af_p_residues]['pocket'].iloc[0]
	## here - add a function that aligns af_p_residues
		af_p_residues=af_p_residues.replace("[","").replace("]","").split(',')
		renumbered_af_p_residues=renumber_uniprot_residues(df_line_sifts,AF_model.split('_')[0],pdb_chain,af_p_residues)
		overlap1,overlap2=pocket_overlap(p_residues, renumbered_af_p_residues)
		#print (cluster_nm,pdb_chain,AF_model,AF_model_pocket,overlap1,overlap2,p_residues,renumbered_af_p_residues)	
		# append the results to the df 
		new_row_dict={'cluster_nm': cluster_nm,'pdb_chain_pocket': pdb_chain+'_'+pocket,'AF_model':AF_model,'AF_model_pocket':AF_model_pocket,'pocket_coverage':overlap1,'AF_coverage':overlap2,'pockets_covered':overlap1>0.6}
		new_row_df = pd.DataFrame([new_row_dict])
		df_codnas_fpockets_compared_to_AF = pd.concat([df_codnas_fpockets_compared_to_AF, new_row_df], ignore_index=True)
# Apply the function to each row
df_codnaSpockets.apply(process_row, axis=1)

df_codnas_fpockets_compared_to_AF.to_csv('codnas_fpockets_compared_to_AF.csv',index=False)
print (df_codnas_fpockets_compared_to_AF)

# create a new df with only the highest coverage AF pocket for each codnas pocket
# Group by the 'pdb_chain_pocket' column and get the index of the row with the highest 'pockets_covered' in each group
idx = df_codnas_fpockets_compared_to_AF.groupby('pdb_chain_pocket')['pockets_covered'].idxmax()
# Use these indices to create the new DataFrame
df_codnas_fpockets_compared_to_AF_highest= df_codnas_fpockets_compared_to_AF.loc[idx].reset_index(drop=True)

df_codnas_fpockets_compared_to_AF_highest.to_csv('codnas_fpockets_compared_to_AFhighest.csv',index=False)
