# read all folder or file names in a directory
import os
from Bio.PDB import PDBParser, is_aa
from Bio.Data import IUPACData
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


def get_sequence_with_resseq(pdb_file):
    parser = PDBParser(QUIET=True)  # Suppress warnings
    sequence = ""
    resseq_list = []
    try:
        structure = parser.get_structure("protein", pdb_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue, standard=True):
                        resname = residue.get_resname().capitalize()  # 3-letter code
                        if resname in IUPACData.protein_letters_3to1:
                            one_letter_code = IUPACData.protein_letters_3to1[resname]
                            sequence += one_letter_code
                            resseq = residue.id[1]  # Get resseq number
                            resseq_list.append(resseq)
    except Exception as e:
        print(f"Error processing PDB file {pdb_file}: {e}")
    return sequence, resseq_list


# def extract_smiles_from_sdf_gz(file_path):
#     smiles_list = []
#     try:
#         with gzip.open(file_path, "rb") as gz_file:
#             suppl = Chem.ForwardSDMolSupplier(gz_file,removeHs=True)
#             for mol in suppl:
#                 if mol is not None:
#                     smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
#                     smiles_list.append(smiles)
#     except Exception as e:
#         print(f"Error extracting SMILES from {file_path}: {e}")
#     return smiles_list


# def extract_smiles_from_ism(file_path):
#     smiles_list = []
#     try:
#         with open(file_path, "r") as file:
#             for line in file.readlines():
#                 smiles_list.append(line.split()[0])
#     except Exception as e:
#         print(f"Error extracting SMILES from {file_path}: {e}")
#     return smiles_list


##new version
def extract_smiles_from_sdf_gz(file_path):
    """Extracts SMILES from a .sdf.gz file."""
    smiles_list = []
    with gzip.open(file_path, "rb") as gz_file:
        suppl = Chem.ForwardSDMolSupplier(gz_file, removeHs=True)
        for mol in suppl:
            if mol is not None:  # Check if the molecule is valid
                if not mol.GetNumConformers():  # Check if it lacks 3D coordinates
                    AllChem.EmbedMolecule(mol)  # Embed 3D coordinates if needed
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                smiles_list.append(smiles)
                smiles_list = list(set(smiles_list))
    return smiles_list


def extract_smiles_from_ism(file_path):
    """Extracts SMILES from a .ism file."""
    smiles_list = []

    with open(file_path, "r") as file:
        for line in file.readlines():
            smiles = line.split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:  # Check if the molecule is valid
                mol = Chem.RemoveHs(mol)  # Remove hydrogen atoms
                if not mol.GetNumConformers():  # Check if it lacks 3D coordinates
                    AllChem.EmbedMolecule(mol)  # Embed 3D coordinates if needed
                smiles = Chem.MolToSmiles(
                    mol, isomericSmiles=True
                )  # Regenerate SMILES after embedding
                smiles_list.append(smiles)
                smiles_list = list(set(smiles_list))
    return smiles_list


def get_smiles_list(sample_dir, file_prefix):
    sdf_path = os.path.join(sample_dir, f"{file_prefix}.sdf.gz")
    mol2_path = os.path.join(sample_dir, f"{file_prefix}.mol2.gz")
    ism_path = os.path.join(sample_dir, f"{file_prefix}.ism")

    try:
        if os.path.exists(sdf_path):
            return extract_smiles_from_sdf_gz(sdf_path)
        elif os.path.exists(mol2_path):
            print(f"{file_prefix}.mol2.gz found, processing...")
            return extract_smiles_from_sdf_gz(mol2_path)
        elif os.path.exists(ism_path):
            print(f"{file_prefix}.ism found, processing...")
            return extract_smiles_from_ism(ism_path)
        else:
            print(f"No {file_prefix} files found for {sample_dir}")
    except Exception as e:
        print(f"Error getting SMILES for {sample_dir}: {e}")
    return []


def get_table(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check if the request was successful
        soup = BeautifulSoup(response.text, "html.parser")
        table = soup.find("table")
        df = pd.read_html(str(table))[0]
    except Exception as e:
        print(f"Error fetching table from {url}: {e}")
        return pd.DataFrame()  # Return an empty DataFrame if there's an error
    return df


def get_pdb_from_table(target_name):
    try:
        for i in range(len(table)):
            if str(table["Target Name"][i]).upper() == target_name.upper():
                return table["PDB"][i]
    except Exception as e:
        print(f"Error finding PDB for target {target_name}: {e}")
    return None


# Path to the DUDE dataset directory
DUDE_base = "local data/DUDE/"
DUDE_dir = DUDE_base + "DUDE samples/"

# Get the table from the URL
url = "http://dude.docking.org/targets"
table = get_table(url)

info = {}
folders = get_folder_names(DUDE_dir)

for sample in tqdm(folders, desc="Processing Samples"):
    print(sample)
    sample_dir = DUDE_dir + sample + "/"
    sample_files = get_files_in_dir(sample_dir)
    info[sample] = {}
    info[sample]["sample_dir"] = sample_dir

    try:
        info[sample]["pdb"] = get_pdb_from_table(sample)

        # Fetch the sequence from a pdb file
        try:
            PDB_seq = get_sequence_with_resseq(sample_dir + "receptor.pdb")
            info[sample]["PDB_seq"] = PDB_seq[0]
        except Exception as e:
            print(f"Failed to fetch PDB sequence for {sample}: {e}")
            info[sample]["PDB_seq"] = None

        # Fetching the sequence from Uniprot
        try:
            uniprot_name = pdbh.fetch_PDB_uniprot_accession(info[sample]["pdb"])[0]
            uniprot_seq = pdbh.fetch_uniprot_sequence(uniprot_name)
            info[sample]["uniprot_name"] = uniprot_name
            info[sample]["uniprot_seq"] = uniprot_seq
        except Exception as e:
            print(f"Failed to fetch UniProt data for {sample}: {e}")
            info[sample]["uniprot_name"] = None
            info[sample]["uniprot_seq"] = None

        # Extracting SMILES strings active
        try:
            info[sample]["actives_smiles"] = get_smiles_list(
                sample_dir, "actives_final"
            )
        except Exception as e:
            print(f"Error processing {sample}: {e}")
            info[sample]["actives_smiles"] = []
        try:
            info[sample]["decoys_smiles"] = get_smiles_list(sample_dir, "decoys_final")
        except Exception as e:
            print(f"Error processing {sample}: {e}")
            info[sample]["decoys_smiles"] = []

    except Exception as e:
        print(f"Error processing {sample}: {e}")

    print(sample + " done")

# Write info to a json file
try:
    with open(DUDE_base + "DUDE_processed.json", "w") as f:
        json.dump(info, f, indent=4)
except Exception as e:
    print(f"Error writing to info.json: {e}")
