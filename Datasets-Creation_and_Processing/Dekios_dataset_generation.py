import os
import urllib.request
import zipfile
import pandas as pd
import os
from aiondata.raw import protein_structure
import os
import gzip
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from tqdm import tqdm

RDLogger.DisableLog("rdApp.warning")
pdbh = protein_structure.PDBHandler()


def download_and_extract_decoys(base_dir, extraction_folder, urls):
    """
    Downloads and extracts decoy files from given URLs.

    Parameters:
    base_dir (str): Base directory for downloads.
    extraction_folder (str): Directory where extracted files will be stored.
    urls (list): List of URLs to download.
    """
    # Define the full path for the extraction folder
    extraction_path = os.path.join(base_dir, extraction_folder)

    # Create the extraction folder if it doesn't exist
    os.makedirs(extraction_path, exist_ok=True)

    # Download and extract each part into the same folder
    for url in urls:
        # Define the file path for the downloaded zip
        file_path = os.path.join(base_dir, os.path.basename(url))

        # Check if the file exists
        if not os.path.exists(file_path):
            print(f"{file_path} does not exist. Downloading...")

            # Download the file
            urllib.request.urlretrieve(url, file_path)
            print("Download complete.")
        else:
            print(f"{file_path} already exists.")

        # Extract the zip file
        print(f"Extracting {file_path}...")
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall(
                extraction_path
            )  # Extract to the specified extraction path
        print("Extraction complete.")

    print("All files downloaded and extracted.")


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

    ### Main code ###

    ## Download and extract the Dekois dataset ##


# Define download URLs
urls = [
    "https://uni-tuebingen.de/securedl/sdl-eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpYXQiOjE3MjcxNTYxNTYsImV4cCI6MTcyNzI0NjE1NSwidXNlciI6MCwiZ3JvdXBzIjpbMCwtMV0sImZpbGUiOiJmaWxlYWRtaW5cL1VuaV9UdWViaW5nZW5cL0Zha3VsdGFldGVuXC9DaGVtaWVQaGFybWFcL0luc3RpdHV0ZVwvUGhhcm0uX0luc3RpdHV0XC9QaGFybS5fQ2hlbWllXC9Nb2xEZXNpZ25cL2Rla29pc1wvZGVjb3lzX3BhcnRfMS56aXAiLCJwYWdlIjoyNjkwNzN9.J7ZvGcoyHv7hmdguy-Q_QxwYdt8Fb4TmhchDz5k58ss/decoys_part_1.zip",
    "https://uni-tuebingen.de/securedl/sdl-eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpYXQiOjE3MjcxNTYxNTYsImV4cCI6MTcyNzI0NjE1NSwidXNlciI6MCwiZ3JvdXBzIjpbMCwtMV0sImZpbGUiOiJmaWxlYWRtaW5cL1VuaV9UdWViaW5nZW5cL0Zha3VsdGFldGVuXC9DaGVtaWVQaGFybWFcL0luc3RpdHV0ZVwvUGhhcm0uX0luc3RpdHV0XC9QaGFybS5fQ2hlbWllXC9Nb2xEZXNpZ25cL2Rla29pc1wvZGVjb3lzX3BhcnRfMi56aXAiLCJwYWdlIjoyNjkwNzN9.KNp1TtFpxYXV2JgJeeYbK6wAd27_MUniZ_BRahnprK0/decoys_part_2.zip",
    "https://uni-tuebingen.de/securedl/sdl-eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpYXQiOjE3MjcxNTYxNTYsImV4cCI6MTcyNzI0NjE1NSwidXNlciI6MCwiZ3JvdXBzIjpbMCwtMV0sImZpbGUiOiJmaWxlYWRtaW5cL1VuaV9UdWViaW5nZW5cL0Zha3VsdGFldGVuXC9DaGVtaWVQaGFybWFcL0luc3RpdHV0ZVwvUGhhcm0uX0luc3RpdHV0XC9QaGFybS5fQ2hlbWllXC9Nb2xEZXNpZ25cL2Rla29pc1wvZGVjb3lzX3BhcnRfMy56aXAiLCJwYWdlIjoyNjkwNzN9.QBg4Hh5YYsYd29nzphYkl1qr44ZLzwqZmsZEDQRHOos/decoys_part_3.zip",
]
Dekios_base = "local data/Dekois/"
Dekios_decoys = os.path.join(Dekios_base, "Decoys")
Dekios_actives = os.path.join(Dekios_base, "Actives")

# Call the function with the base directory, extraction folder, and URLs
download_and_extract_decoys(Dekios_base, "Decoys", urls)
# Define download URLs
urls = [
    "https://uni-tuebingen.de/securedl/sdl-eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpYXQiOjE3MjcxNTYxNTYsImV4cCI6MTcyNzI0NjE1NSwidXNlciI6MCwiZ3JvdXBzIjpbMCwtMV0sImZpbGUiOiJmaWxlYWRtaW5cL1VuaV9UdWViaW5nZW5cL0Zha3VsdGFldGVuXC9DaGVtaWVQaGFybWFcL0luc3RpdHV0ZVwvUGhhcm0uX0luc3RpdHV0XC9QaGFybS5fQ2hlbWllXC9Nb2xEZXNpZ25cL2Rla29pc1wvbGlnYW5kcy56aXAiLCJwYWdlIjoyNjkwNzN9.XJF85c9VekiRX_E_I5dNbjPglWRpMaXClQrIjS7VjlY/ligands.zip"
]
# Call the function with the base directory, extraction folder, and URLs
download_and_extract_decoys(Dekios_base, "Actives", urls)


## Process the Dekois supplentarry table (use pdb ids) ##
df = pd.read_excel(f"{Dekios_base}/tables2.xlsx", header=None, keep_default_na=False)


## Process Dekois ##
# Run per row in the df with tqdm for progress and error handling
info = {}
for i in tqdm(range(df.shape[0]), desc="Processing samples"):
    try:
        sample = df.iloc[i, 0]
        pdb = df.iloc[i, 1]
        info[sample] = {}
        info[sample]["pdb"] = pdb
        print(f"Processing sample {sample}")
        # Fetch UniProt data
        try:
            uniprot_name = pdbh.fetch_PDB_uniprot_accession(info[sample]["pdb"])[0]
            uniprot_seq = pdbh.fetch_uniprot_sequence(uniprot_name)
            info[sample]["uniprot_name"] = uniprot_name
            info[sample]["uniprot_seq"] = uniprot_seq
        except Exception as e:
            print(f"Error fetching UniProt data for sample {sample}, PDB {pdb}: {e}")
            info[sample]["uniprot_name"] = None
            info[sample]["uniprot_seq"] = None

        # Process active and decoy directories
        # sample2 = sample.replace("-", "")
        decoy_dir = f"{Dekios_decoys}/{sample}_Celling-v1.12_decoyset.sdf.gz"
        active_dir = f"{Dekios_actives}/ligands/{sample}.sdf.gz"

        # Extract SMILES with error handling
        try:
            info[sample]["actives_smiles"] = extract_smiles_from_sdf_gz(active_dir)
        except Exception as e:
            print(f"Error extracting active SMILES for sample {sample}: {e}")
            info[sample]["actives_smiles"] = None

        try:
            info[sample]["decoys_smiles"] = extract_smiles_from_sdf_gz(decoy_dir)
        except Exception as e:
            print(f"Error extracting decoy SMILES for sample {sample}: {e}")
            info[sample]["decoys_smiles"] = None

    except Exception as e:
        print(f"Error processing sample {i}: {e}")

    ## save the info ##
import json

with open(f"{Dekios_base}/Dekios_processed.json", "w") as f:
    json.dump(info, f, indent=2)
