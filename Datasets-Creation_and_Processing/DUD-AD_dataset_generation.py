import os
import urllib.request
import requests
import pandas as pd
from bs4 import BeautifulSoup
from Bio.PDB import PDBParser, is_aa
from Bio.Data import IUPACData
import os
import gzip
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from aiondata.raw import protein_structure
from tqdm import tqdm
from aiondata.raw import protein_structure
import json

pdbh = protein_structure.PDBHandler()
RDLogger.DisableLog("rdApp.warning")


def get_table(url):
    # Send a GET request to fetch the page content
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful

    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, "html.parser")

    # Find the table in the HTML content
    table = soup.find("table")

    # Read the table into a DataFrame
    df = pd.read_html(str(table))[0]

    # Show the DataFrame
    return df


# find the target name in the table and return its PDB from the table
def get_pdb_from_table(target_name):
    for i in range(len(table)):
        if str(table["Target Name"][i]).upper() == target_name.upper():
            return table["PDB"][i]
    return None


def get_files_in_dir(directory):
    files = []
    for file in os.listdir(directory):
        if os.path.isfile(directory + file):
            files.append(file)
    return files


def get_sequence_with_resseq(pdb_file):
    """
    Reads a PDB file and returns the amino acid sequence and resseq numbers as a list of tuples.

    Parameters:
    pdb_file (str): Path to the PDB file.

    Returns:
    tuple: (str, list) Amino acid sequence and list of (residue sequence number, 1-letter code).
    """
    parser = PDBParser(QUIET=True)  # Suppress warnings
    structure = parser.get_structure("protein", pdb_file)
    sequence = ""
    resseq_list = []

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

    return sequence, resseq_list


# exaple
# PDB_seq=get_sequence_with_resseq(f"{DUDE_original_samples}{sample}/receptor.pdb")
# info[sample]["PDB_seq"]=PDB_seq[0]


def extract_smiles_from_sdf_gz(file_path):
    """Extracts SMILES from a .sdf.gz file."""
    smiles_list = []

    suppl = Chem.ForwardSDMolSupplier(file_path, removeHs=True)
    for mol in suppl:
        if mol is not None:  # Check if the molecule is valid
            if not mol.GetNumConformers():
                # Embed the molecule to generate 3D coordinates if needed
                AllChem.EmbedMolecule(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            smiles_list.append(smiles)
            # unique list
            smiles_list = list(set(smiles_list))

    return smiles_list


#### Main code ####
DUDE_AD_base = "local data/DUDE/"
DUDE_AD_dir = DUDE_AD_base + "102_AD_dataset/"
DUDE_original_samples = DUDE_AD_base + "DUDE samples/"


## Download the DUDE AD dataset ##
file_path = DUDE_AD_base + "102_targets_AD_dataset.tar"
url = "http://www.lehman.edu/faculty/tkurtzman/files/102_targets_AD_dataset.tar"

# Check if the file exists
if not os.path.exists(file_path):
    print(f"{file_path} does not exist. Downloading...")
    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Download the file
    urllib.request.urlretrieve(url, file_path)
    print("Download complete.")
else:
    print(f"{file_path} already exists.")

# Extract the tar file
if os.path.exists(file_path):
    print(f"Extracting {file_path}...")
    # Ensure to handle spaces by adding quotes around the path
    os.system(f'tar -xf "{file_path}" -C "{DUDE_AD_base}"')
    print("Extraction complete.")


## Fetch the protein names from the DUDE AD site ##
# URL of the site
url = "https://dude.docking.org/targets"
# Get the table from the URL
table = get_table(url)


## Get the unique protein names from the files in the directory ##
files = get_files_in_dir(DUDE_AD_dir)
unique_proteins = list(set([file.split("_")[0] for file in files]))


## Process the samples ##
info = {}

# Wrap the loop with tqdm to show progress
for sample in tqdm(unique_proteins, desc="Processing samples", unit="sample"):
    print(f"Processing {sample}...")
    info[sample] = {}
    active_dir = f"{DUDE_AD_dir}{sample}_actives.sdf"
    decoy_dir = f"{DUDE_AD_dir}{sample}_AD.sdf"

    # Get the PDB ID from the site
    try:
        info[sample]["pdb"] = get_pdb_from_table(sample)
    except Exception as e:
        print(f"Error fetching PDB ID for {sample}: {e}")
        info[sample]["pdb"] = None  # Assign None or handle accordingly

    if info[sample]["pdb"] is not None:
        try:
            uniprot_name = pdbh.fetch_PDB_uniprot_accession(info[sample]["pdb"])[0]
            uniprot_seq = pdbh.fetch_uniprot_sequence(uniprot_name)
            info[sample]["uniprot_name"] = uniprot_name
            info[sample]["uniprot_seq"] = uniprot_seq
        except Exception as e:
            print(f"Error fetching UniProt info for {info[sample]['pdb']}: {e}")
            info[sample]["uniprot_name"] = None
            info[sample]["uniprot_seq"] = None

    # Extract SMILES from active and decoy SDF files
    try:
        info[sample]["actives_smiles"] = extract_smiles_from_sdf_gz(active_dir)
    except Exception as e:
        print(f"Error extracting SMILES from {active_dir}: {e}")
        info[sample]["actives_smiles"] = None

    try:
        info[sample]["decoys_smiles"] = extract_smiles_from_sdf_gz(decoy_dir)
    except Exception as e:
        print(f"Error extracting SMILES from {decoy_dir}: {e}")
        info[sample]["decoys_smiles"] = None

    # Get sequence from PDB file
    try:
        PDB_seq = get_sequence_with_resseq(
            f"{DUDE_original_samples}{sample}/receptor.pdb"
        )
        info[sample]["PDB_seq"] = PDB_seq[0]
    except Exception as e:
        print(f"Error fetching PDB sequence for {sample}: {e}")
        info[sample]["PDB_seq"] = None

# Save the info dictionary to a file
info_file = DUDE_AD_base + "AD_dataset_info.json"
with open(info_file, "w") as f:
    json.dump(info, f)
print(f"Info saved to {info_file}")
