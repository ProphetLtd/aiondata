from Bio.PDB import PDBParser, PDBIO
import pandas as pd

def split_pdb_by_chain(input_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_file)

    io = PDBIO()

    for model in structure:
        for chain in model:
            chain_id = chain.id
            output_file = f"{input_file[:-4]}_{chain_id}.pdb"

            io.set_structure(chain)
            io.save(output_file)
            print(f"Chain {chain_id} saved to {output_file}")

def read_BioLip():
	columns=['pdb_id', 'chain', 'Res', 'BS', 'ligCCD', 'lig_chain', 'Lig_ser_num', 'BSres_frompdb', 'BSres_fromscratch', 'Catres_frompdb', 'Catres_fromscratch', 'ECnum', 'GOterms', 'affinity', 'affinity_MOAD', 'aff_CNPDBbind', 'affinity_BDB', 'UnirpotID', 'PubmedID', 'lig_res_seq', 'Receptor_seq']
	df=pd.read_csv('/home/tali/Prophet/BioLip/BioLiP.txt.gz',sep='\t',names=columns, header=None)
	return df
def read_AA_3_to_1():
	# Initialize an empty dictionary
	residue_dict = {}
	with open('AA_1_to_3.txt', 'r') as file:
		for line in file:
			oneRes, threeRes = line.split()[0],line.split()[1]
			residue_dict[oneRes] = threeRes
	return residue_dict		
