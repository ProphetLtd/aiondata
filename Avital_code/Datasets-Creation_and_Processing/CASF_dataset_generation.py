# read all folder or file names in a directory
import aiondata
import os
from Bio.PDB import PDBParser, is_aa
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
from rdkit import Chem
import json
from tqdm import tqdm
import requests
import pandas as pd
from bs4 import BeautifulSoup
from aiondata.raw import protein_structure
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.warning")
pdbh = protein_structure.PDBHandler()

def get_files_in_dir(directory):
	files = []
	try:
		for file in os.listdir(directory):
			if os.path.isfile(directory + file):
				files.append(file)
	except Exception as e:
		print(f"Error reading files in {directory}: {e}")
	return files


def get_folder_names(directory):
    folder_names = []
    try:
        for root, dirs, files in os.walk(directory):
            for name in dirs:
                folder_names.append(name)
    except Exception as e:
        print(f"Error reading folders in {directory}: {e}")
    return folder_names

# Function to extract the sequence
def extract_sequence_from_pdb(pdb_file):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        # Dictionary for 3-letter to 1-letter amino acid conversion
        aa_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'  }
        sequence = ""
        for model in structure:
                for chain in model:
                        for residue in chain:
                                res_name = residue.get_resname()
                                # Convert 3-letter code to 1-letter code
                                if res_name in aa_dict:
                                        sequence += aa_dict[res_name]
        return str(SeqRecord(Seq(sequence), id="PDB_sequence").seq)

def extract_uniprot_seq(pdb):
	try:
		uniprot_name = pdbh.fetch_PDB_uniprot_accession(pdb)[0]
	except:
		uniprot_sequence = None
		uniprot_name = None
		print (pdb, uniprot_name)
	if uniprot_name:
		uniprot_sequence = pdbh.fetch_uniprot_sequence(uniprot_name)
	else:
		uniprot_sequence = None
			
	return uniprot_sequence

def extract_smiles_from_mol2(file_path):
	"""Extracts isomeric SMILES from a .mol2 file."""
	try:
		# First attempt with sanitization
		mol = Chem.MolFromMol2File(file_path, removeHs=True)
		if mol is None:
			print("Issue with sanitization. Attempting without sanitization.")
			mol = Chem.MolFromMol2File(file_path, removeHs=True, sanitize=False)
		if mol is not None:
			return Chem.MolToSmiles(mol, isomericSmiles=True)
		else:
			print(f"Failed to read molecule from {file_path}")
			return None
	except Exception as e:
		print(f"An error occurred: {str(e)}")
		return None


#make a table
s=extract_sequence_from_pdb('1a30/1a30_protein.pdb')
print (s)
mol=extract_smiles_from_mol2('1a30/1a30_ligand.mol2')
print (str(mol))
extract_uniprot_seq('3gnw')

#-----

CASF_df = pd.DataFrame(columns=['pdb','uniprot_seq','pdb_seq','smiles' ])

folder_names=get_folder_names('./')
folder_names= [f for f in folder_names if f!= '__pycache__' and f!= 'obsolete' ]
for pdb in folder_names:
	print (pdb)
	pdb_file=pdb+'/'+pdb+'_protein.pdb'
	lig_file=pdb+'/'+pdb+'_ligand.mol2'
	pdb_dict={'pdb': pdb, 'uniprot_seq': extract_uniprot_seq(pdb), 'pdb_seq': extract_sequence_from_pdb(pdb_file) ,'smiles': extract_smiles_from_mol2(lig_file) }
	new_row_df = pd.DataFrame([pdb_dict])
	CASF_df = pd.concat([CASF_df, new_row_df], ignore_index=True)
CASF_df.to_csv('CASF_core_dataset', index=False)
