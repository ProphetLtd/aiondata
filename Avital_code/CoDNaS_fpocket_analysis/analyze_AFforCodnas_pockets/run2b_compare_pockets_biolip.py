from collect_pocket import pocket_overlap 
import pandas as pd

# read covered pockets into df
df_covered_pockets=pd.read_csv('codnas_fpockets_compared_to_AFhighest.csv')
# read codnas pockets that are bound according to BioLip into df
df_codnas_biolip=pd.read_csv('../analyze_condas_pockets/codnas_BioLipLigBound_fpockets')
df_codnas_biolip['pdb_chain_pocket']=df_codnas_biolip['pdb_chain']+'_'+df_codnas_biolip['pocket']
#read apo pockets which are found as holo in other conformation and those that are not
df_apo_w_holo=pd.read_csv('../analyze_condas_pockets/apo_with_holo')
df_apo_wout_holo=pd.read_csv('../analyze_condas_pockets/apo_wout_holo')
df_apo_w_holo['pdb_chain_pocket']= df_apo_w_holo['pdb_chain']+'_'+df_apo_w_holo['pocket']
df_apo_wout_holo['pdb_chain_pocket']= df_apo_wout_holo['pdb_chain']+'_'+df_apo_wout_holo['pocket']

# which from the covered pockets is a ligand bound biolip pocket
df_covered_biolip = df_covered_pockets[df_covered_pockets['pdb_chain_pocket'].isin(df_codnas_biolip['pdb_chain_pocket'])]
df_covered_biolip.to_csv('codnas_fpockets_compared_to_AFhighest_biolip.csv',index=False)

# which from the NoBioLip=apo pockets are actually holo in other conformations and which are not
df_covered_Nobiolip = df_covered_pockets[~df_covered_pockets['pdb_chain_pocket'].isin(df_codnas_biolip['pdb_chain_pocket'])]
df_covered_Nobiolip.to_csv('codnas_fpockets_compared_to_AFhighest_Nobiolip.csv',index=False)

df_covered_Nobiolip_w_holoConf = df_covered_Nobiolip[df_covered_Nobiolip['pdb_chain_pocket'].isin(df_apo_w_holo['pdb_chain_pocket'])]
df_covered_Nobiolip_wout_holoConf = df_covered_Nobiolip[df_covered_Nobiolip['pdb_chain_pocket'].isin(df_apo_wout_holo['pdb_chain_pocket'])]

df_covered_Nobiolip_w_holoConf.to_csv('codnas_fpockets_compared_to_AFhighest_Nobiolip_wHoloConf',index=False)
df_covered_Nobiolip_wout_holoConf.to_csv('codnas_fpockets_compared_to_AFhighest_Nobiolip_woutHoloConf',index=False)
