from collect_pocket import pocket_overlap
import pandas as pd                                                 
import ast 

#1) read codnas pockets
#  read codnas significant & reduced pockets
df_codnaSpockets=pd.read_csv('codnas_Signif_reduced_fpockets')
# read codnas BioLip with ligand pockets (these are taken from the larger non reduced and non-signif codnas )
df_codnaSpockets_BioLip=pd.read_csv('codnas_BioLipLigBound_fpockets')

#2) convert first the resiudes pockets columns to actual lists of residues !!
# Function to convert the custom string format to actual list
def convert_string_to_list(s):
    # Remove brackets and split by commas
    return s.strip('[]').split(',')
# Apply the function to the column
df_codnaSpockets['residues'] = df_codnaSpockets['residues'].apply(convert_string_to_list)
df_codnaSpockets_BioLip['residues'] = df_codnaSpockets_BioLip['residues'].apply(convert_string_to_list)

#3) go over pockets and find all the apo pockets to work on them
# Extract first three columns from both DataFrames
df_codnaSpockets_pocket = df_codnaSpockets[['cluster_nm', 'pdb_chain', 'pocket']]
df_codnaSpockets_BioLip_pocket = df_codnaSpockets_BioLip[['cluster_nm', 'pdb_chain', 'pocket']]
# Filter dfpockets by checking if the first three items are not in BioLip-holo 
df_codnaSpockets_apo = df_codnaSpockets[~df_codnaSpockets_pocket.apply(tuple, axis=1).isin(df_codnaSpockets_BioLip_pocket.apply(tuple, axis=1))]

print (df_codnaSpockets_apo)

##4) check if remaining apo pockets are actually such that have a holo form as well, by overlapping pocket residues in different conformation
# Define a function to check if 70% overlap occurs for a pocket with with any row in holo pockets 
def has_overlap(list1):
    return any(pocket_overlap(list1, list2)[0] >= 0.7 for list2 in df_codnaSpockets_BioLip['residues']) 
# Filter apo pockets into two types based on the 70% overlap with holo criteria
apo_with_holo = df_codnaSpockets_apo[df_codnaSpockets_apo['residues'].apply(has_overlap)]
apo_wout_holo = df_codnaSpockets_apo[~df_codnaSpockets_apo['residues'].apply(has_overlap)]

apo_with_holo.to_csv('apo_with_holo')
apo_wout_holo.to_csv('apo_wout_holo')
