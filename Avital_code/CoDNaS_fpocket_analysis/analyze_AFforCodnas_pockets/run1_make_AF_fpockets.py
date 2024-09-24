import os
from collect_pocket import collect_pocket
import pandas as pd

def read_codnas():
	df_codnas=pd.read_csv('/data/codnas_dataset/codnas_dataset.csv')
	return df_codnas

###! define pred_fpocket_df and fill a new row with pdb and its pockets for each codnas entry
pred_fpocket_df = pd.DataFrame(columns=['cluster_nm','pdb_chain','pocket','Pscore','Dscore','residues'])
pred_signif_fpocket_df = pd.DataFrame(columns=['cluster_nm','pdb_chain','pocket','Pscore','Dscore','residues'])
##read codnas ensemble and loop over its pdb_chains, to create and rows to the pockts df
df_codnas=read_codnas()
for cluster_nm in df_codnas['Pool_repr'].values:
	pdb=str(df_codnas[df_codnas['Pool_repr']==cluster_nm]['uniprot'].iloc[0])
	# go over all the avialble models of the pdb and their pockets to add to the fpocket df. find who are they from curr dir
	# remove the '_out' from the name, to stay with pdb_chain alone , chain here is 'model_X'
	pdb_chains = [f[:-4] for f in os.listdir('.') if f.startswith(pdb) and f.endswith('.pdb')]
	for pdb_chain in pdb_chains:
		pockets=collect_pocket(pdb_chain)
		for pocket_num in pockets:
			Pscore,Dscore,residues=pockets[pocket_num]
			residues ='[' + ','.join(residues) + ']' # residues is a list 1_ASN etc
			pdb_chain_dict={'cluster_nm': cluster_nm, 'pdb_chain': pdb_chain, 'pocket': pocket_num,'Pscore': Pscore,'Dscore': Dscore, 'residues': residues}
			new_row_df = pd.DataFrame([pdb_chain_dict])
			pred_fpocket_df = pd.concat([pred_fpocket_df, new_row_df], ignore_index=True)
			# only when pockets are significant:
			if float(Pscore)>0.2 and float(Dscore)>0.2:
				pred_signif_fpocket_df= pd.concat([pred_signif_fpocket_df, new_row_df], ignore_index=True)
pred_fpocket_df.to_csv('AF_fpockets', index=False)
pred_signif_fpocket_df.to_csv('AF_Signif_fpockets', index=False)	


