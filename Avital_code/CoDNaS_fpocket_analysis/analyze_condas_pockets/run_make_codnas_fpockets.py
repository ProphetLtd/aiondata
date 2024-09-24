import os, re
import pandas as pd
from collect_pocket import collect_pocket,pocket_overlap
from itertools import combinations
from functions import read_BioLip,read_AA_3_to_1
 
def read_codnas():
	df_codnas=pd.read_csv('/data/codnas_dataset/codnas_dataset.csv')	
	return df_codnas

# read BioLip into format similar to codnas:
df_biolip=read_BioLip()
df_biolip['pdb_chain']=df_biolip['pdb_id'].str.upper()+'_'+df_biolip['chain']
# rm lines where binding residues contain 'X' residue
df_biolip = df_biolip[~df_biolip['BSres_frompdb'].str.contains('X')]

residue_dict=read_AA_3_to_1()
	# convert residues gets the list of residues from biolip and change to codnas format:
def convert_residues(res_string):
	residues,residues_new=res_string.split(),[]
	for res in residues:	
		let,num=res[0],res[1:]
		resnew=num+'_'+residue_dict[let]
		residues_new.append(resnew)
	return residues_new
# convert biolip binding site residues to same format codnas will have 
df_biolip['residues'] = df_biolip['BSres_frompdb'].apply(convert_residues)

###! define ensemble_fpocket_df and fill a new row with pdb and its pockets for each codnas entry
#----------------------------------------
ensemble_fpocket_df = pd.DataFrame(columns=['cluster_nm','pdb_chain','pocket','Pscore','Dscore','residues','Aspheres','Vol'])
ensemble_signif_fpocket_df = pd.DataFrame(columns=['cluster_nm','pdb_chain','pocket','Pscore','Dscore','residues','Aspheres','Vol'])
ensemble_BioLipLigBound_fpocket_df = pd.DataFrame(columns=['cluster_nm','pdb_chain','pocket','Pscore','Dscore','residues','Aspheres','Vol'])
##read codnas ensemble and loop over its pdb_chains, to create and rows to the pockts df
df_codnas=read_codnas()
for cluster_nm in df_codnas['Pool_repr'].values:
	for pdb_chain in df_codnas[df_codnas['Pool_repr']==cluster_nm]['Cluster'].iloc[0].replace("'","").replace("[","").replace("]","").split(','):
		pdb_chain=pdb_chain.replace(" ","")
		pockets=collect_pocket(pdb_chain)
		for pocket_num in pockets:
			Pscore,Dscore,residues,Aspheres,Vol=pockets[pocket_num]
			residues ='[' + ','.join(residues) + ']' 
			pdb_chain_dict={'cluster_nm': cluster_nm, 'pdb_chain': pdb_chain, 'pocket': pocket_num,'Pscore': Pscore,'Dscore': Dscore, 'residues': residues, 'Aspheres': Aspheres, 'Vol':Vol}
			new_row_df = pd.DataFrame([pdb_chain_dict])		
			ensemble_fpocket_df = pd.concat([ensemble_fpocket_df, new_row_df], ignore_index=True) 
		# check if the pocket is an actual solved ligand bound pocket
			# for that look at the pdb chain of nmr cases without the nmr model, because there is no such in biolip
			if '-' in pdb_chain:
				pdb_chain_noModel=re.sub(r'-\d+', '', pdb_chain)
			else:
				pdb_chain_noModel=pdb_chain
			df_biolip_pdb_chain = df_biolip[df_biolip['pdb_chain']==pdb_chain_noModel]
			for index, row in df_biolip_pdb_chain.iterrows():
				residues_cod=residues.replace("[","").replace("]","").split(',') 
				overlap1,overlap2=pocket_overlap(residues_cod,row['residues'])
				if float(overlap2)>0.7:
					ensemble_BioLipLigBound_fpocket_df= pd.concat([ensemble_BioLipLigBound_fpocket_df, new_row_df], ignore_index=True)	# this is not subset of significant pockets. some biolip pockets are with <0.3 Pscore
			if float(Pscore)>0.3 and float(Dscore)>0.3:	
				ensemble_signif_fpocket_df= pd.concat([ensemble_signif_fpocket_df, new_row_df], ignore_index=True)
ensemble_fpocket_df.to_csv('codnas_fpockets', index=False)	
ensemble_signif_fpocket_df.to_csv('codnas_Signif_fpockets', index=False)	
ensemble_BioLipLigBound_fpocket_df.to_csv('codnas_BioLipLigBound_fpockets', index=False)
# reduce pockets which are actually the same one (in different proteins, but represent same pocket) 
# Function to filter rows within each group
def filter_group(group_df):
	to_drop = set()
	for (i1, row1), (i2, row2) in combinations(group_df.iterrows(), 2):
		if i1 not in to_drop and i2 not in to_drop:
			p1,p2=row1['residues'].replace("'","").replace("[","").replace("]","").split(','),row2['residues'].replace("'","").replace("[","").replace("]","").split(',')
			overlap1,overlap2 = pocket_overlap(p1,p2)
			alpha1,alpha2=float(row1['Aspheres']),float(row2['Aspheres'])
			print (row1,row2,alpha1,alpha2)
			if len(p1)>=len(p2) and overlap2 > 0.95 and abs(alpha2-alpha1)<3.0: # essentially 100%
				to_drop.add(i2)
			elif len(p2) >= len(p1) and overlap1>0.95 and abs(alpha2-alpha1)<3.0: # essentially 100%
				to_drop.add(i1)
	return group_df.drop(to_drop)

# Apply filtering function to each group
result_df = ensemble_signif_fpocket_df.groupby('cluster_nm').apply(filter_group).reset_index(drop=True)
result_df.to_csv('codnas_Signif_reduced_fpockets', index=False)
