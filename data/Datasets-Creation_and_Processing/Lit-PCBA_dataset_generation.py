# read all folder or file names in a directory
import sys
import os
from Bio.Data import IUPACData
from rdkit import Chem
import tempfile
from rdkit import Chem
from rdkit import RDLogger
from aiondata.raw import protein_structure
import requests
from bs4 import BeautifulSoup
import pandas as pd
from Bio import Entrez, SeqIO
import concurrent.futures
from tqdm import tqdm

pdbh = protein_structure.PDBHandler()


# Suppress RDKit warnings
RDLogger.DisableLog("rdApp.*")
RDLogger.DisableLog("rdApp.warning")


def get_files_in_dir(directory):
    files = []
    for file in os.listdir(directory):
        if os.path.isfile(directory + file):
            files.append(file)
    return files


def get_folder_names(directory):
    folder_names = []
    for root, dirs, files in os.walk(directory):
        for name in dirs:
            folder_names.append(name)
    return folder_names


def get_sequence_with_mol2(mol2_file):
    sequence = []
    previous_residue = None
    processing = False  # Flag to indicate if we're processing atom lines

    with open(mol2_file, "r") as mol2:
        for line in mol2:
            # Check for the start of atom section
            if line.startswith("@<TRIPOS>ATOM"):
                processing = True
                continue  # Skip the line that starts the atom section

            # Check for the next section
            if line.startswith("@"):
                processing = False
                continue  # Stop processing when a new section starts

            if processing:
                # Split the line into columns
                columns = line.split()
                if len(columns) >= 9:
                    # Extract residue name and number
                    residue_name = columns[7][
                        :3
                    ]  # First 3 letters are residue code (e.g., 'THR')
                    residue_number = columns[6]  # 7th column is the residue number

                    # Only append if it's a new residue
                    if residue_number != previous_residue:
                        residue_name2 = residue_name.lower().capitalize()
                        if residue_name2 in IUPACData.protein_letters_3to1:
                            sequence.append(
                                IUPACData.protein_letters_3to1[residue_name2]
                            )
                        previous_residue = residue_number

    # Join the sequence list into a single string
    return "".join(sequence)


def fix_cationic_carbons(mol):
    """Fix cationic carbons by changing them to regular carbons."""
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 1
        ):  # Carbon with +1 charge
            atom.SetFormalCharge(0)
    return mol


def extract_smiles_from_mol2_CarbonicFix(file_path):
    """Extracts isomeric SMILES from a .mol2 file."""
    try:
        # Read the mol2 file as text
        with open(file_path, "r") as f:
            mol2_text = f.read()

        # Replace C.cat with C
        mol2_text = mol2_text.replace("C.cat", "C")

        # Create a temporary file with the modified content
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".mol2", delete=False
        ) as temp_file:
            temp_file.write(mol2_text)
            temp_file_path = temp_file.name

        # First attempt with sanitization
        mol = Chem.MolFromMol2File(temp_file_path, removeHs=True)

        if mol is None:
            print("Issue with sanitization. Attempting without sanitization.")
            mol = Chem.MolFromMol2File(temp_file_path, removeHs=True, sanitize=False)
            return mol
        if mol is not None:
            return mol
        else:
            print(f"Failed to read molecule from {file_path}")
            return None
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


# # Example usage:
# file_path = "../../local data/LIT-PCBA/PDBs/ALDH1/4wp7_ligand.mol2"
# smiles = extract_smiles_from_mol2_CarbonicFix(file_path)
# if smiles:
#     print(f"SMILES: {smiles}")
# else:
#     print("Failed to extract SMILES.")
# from rdkit import Chem


def extract_smiles_from_mol2(file_path):
    """Extracts isomeric SMILES from a .mol2 file."""
    try:
        # First attempt with sanitization
        mol = Chem.MolFromMol2File(file_path, removeHs=True)

        if mol is None:
            print("Issue with sanitization. Attempting without sanitization.")
            mol = Chem.MolFromMol2File(file_path, removeHs=True, sanitize=False)
            if mol is None:
                # Try to fix cationic carbons
                print("Attempting to fix cationic carbons.")
                mol = extract_smiles_from_mol2_CarbonicFix(file_path)
                if mol is None:
                    print("Failed to fix cationic carbons.")
                    return None
        if mol is not None:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        else:
            print(f"Failed to read molecule from {file_path}")
            return None
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


# # Example usage:
# file_path = "../../local data/LIT-PCBA/PDBs/ALDH1/4wp7_ligand.mol2"
# smiles = extract_smiles_from_mol2(file_path)
# if smiles:
#     print(f"SMILES: {smiles}")
# else:
#     print("Failed to extract SMILES.")


def fix_cationic_carbons(mol):
    """Fixes cationic carbons by resetting their formal charge."""
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 1
        ):  # Carbon with +1 charge
            atom.SetFormalCharge(0)
    return mol


def extract_smiles_from_mol2_CarbonicFix(file_path):
    """Attempts to extract molecule from a .mol2 file after fixing cationic carbons."""
    try:
        # Read and modify the mol2 file content
        with open(file_path, "r") as f:
            mol2_text = f.read()

        # Replace 'C.cat' with 'C' to correct cationic carbons
        mol2_text = mol2_text.replace("C.cat", "C")

        # Write the modified content to a temporary file
        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".mol2", delete=False
        ) as temp_file:
            temp_file.write(mol2_text)
            temp_file_path = temp_file.name

        # Try loading the modified molecule file
        mol = Chem.MolFromMol2File(temp_file_path, removeHs=True)
        if mol is None:
            print("Issue with sanitization. Trying without sanitization.")
            mol = Chem.MolFromMol2File(temp_file_path, removeHs=True, sanitize=False)

        if mol is None:
            print(f"Failed to read molecule from {temp_file_path}")
        return mol

    except Exception as e:
        print(f"An error occurred during CarbonicFix: {str(e)}")
        return None


def extract_smiles_from_mol2(file_path):
    """Extracts isomeric SMILES from a .mol2 file, fixing cationic carbons if necessary."""
    try:
        # Try to read the molecule with sanitization
        mol = Chem.MolFromMol2File(file_path, removeHs=True)

        if mol is None:
            print("Sanitization issue. Trying without sanitization.")
            mol = Chem.MolFromMol2File(file_path, removeHs=True, sanitize=False)

            if mol is None:
                print("Sanitization failed. Attempting to fix cationic carbons.")
                mol = extract_smiles_from_mol2_CarbonicFix(file_path)

        # Return the SMILES if the molecule was successfully loaded
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=True)

        print(f"Failed to read molecule from {file_path}")
        return None

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


# # Example usage
# file_path = "../../local data/LIT-PCBA/PDBs/ALDH1/4wp7_ligand.mol2"
# smiles = extract_smiles_from_mol2(file_path)
# if smiles:
#     print(f"SMILES: {smiles}")
# else:
#     print("Failed to extract SMILES.")

# read a txt file and make a list per row.
# also take only the first column


def get_smiles_file(file_path):
    """Reads a text file and returns a list of lines."""
    with open(file_path, "r") as file:
        lines = file.readlines()
    return [line.split()[0] for line in lines]


def process_smiles(smiles):
    """
    Processes a single SMILES string, embedding if necessary and converting to canonical SMILES.

    Args:
        smiles (str): The SMILES string to process.

    Returns:
        str: A canonical SMILES string or None if invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    return smiles


def extract_smiles_from_ism(file_path):
    """
    Extracts and processes SMILES strings from an .ism file,
    removing duplicates and ensuring valid molecular structures, using parallel processing.

    Args:
        file_path (str): The path to the .ism file.

    Returns:
        list: A list of unique and valid SMILES strings.
    """

    with open(file_path, "r") as file:
        smiles_list = [line.split()[0] for line in file]

    smiles_set = set()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = []
        for result in tqdm(
            executor.map(process_smiles, smiles_list),
            total=len(smiles_list),
            desc="Processing SMILES",
        ):
            results.append(result)
            sys.stdout.flush()  # Force output flush

    for smiles in results:
        if smiles:
            smiles_set.add(smiles)

    return list(smiles_set)


# # Example usage (Make sure to replace with your actual path)
# ism_path = "../../local data/LIT-PCBA/PDBs/ADRB2/inactives.smi"
# inactives = extract_smiles_from_ism(ism_path)
# print(inactives)


def extract_table_to_pandas(url, header=True):
    """
    Extracts the table from the given URL and returns a Pandas DataFrame.

    Args:
      url: The URL of the webpage containing the table.
      header: A boolean value indicating whether the table has a header row.
        Defaults to True. If False, no header will be used.

    Returns:
      A Pandas DataFrame containing the table data.
    """

    try:
        # Fetch the webpage content
        response = requests.get(url)
        response.raise_for_status()

        # Parse the HTML content
        soup = BeautifulSoup(response.content, "html.parser")

        # Find the table (assuming it's the only table on the page)
        table = soup.find("table")

        # Extract the table data using Pandas' read_html with header option
        df = pd.read_html(str(table), header=0 if header else None)[0]

        return df

    except Exception as e:
        print(f"Error: Could not extract table data - {e}")
        return None


# # Example usage with header
# url = "https://drugdesign.unistra.fr/LIT-PCBA/index.html"
# df_with_header = extract_table_to_pandas(url, header=True)


# df_with_header


def aid_to_seq(aid):
    """
    Retrieves the protein sequence associated with a given assay ID (aid) from PubChem.

    Args:
      aid: The PubChem assay ID.

    Returns:
      The protein sequence as a string, or None if not found.
    """
    try:
        # Construct the API request URL for the assay description
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/description/JSON"

        # Send the request and raise an exception for bad status codes
        response = requests.get(url)
        response.raise_for_status()

        # Parse the JSON response
        data = response.json()

        # Extract the protein accession
        protein_accession = data["PC_AssayContainer"][0]["assay"]["descr"]["target"][0][
            "mol_id"
        ]["protein_accession"]

        # Set your email address for Entrez
        Entrez.email = "your_email@example.com"

        # Fetch the protein record from NCBI
        handle = Entrez.efetch(
            db="protein", id=protein_accession, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Return the protein sequence
        return str(record.seq)

    except Exception as e:
        print(f"Error: Could not retrieve sequence - {e}")
        return None


# # Example usage
# aid = 2101
# sequence = aid_to_seq(aid)

# if sequence:
#   print(f"Protein Sequence for aid {aid}:\n{sequence}")

#### Main ####


LIT_base = "local data/LIT-PCBA/"
LIT_dir = LIT_base + "PDBs/"
info = {}
folders = get_folder_names(LIT_dir)
url = "https://drugdesign.unistra.fr/LIT-PCBA/index.html"
LIT_site_df = extract_table_to_pandas(url, header=True)

# esr gene are not written as ESR1 in the folder names - fix
if LIT_site_df["Set"].str.contains("ESR").any():
    # add the 1 to the ESR - ESR1
    LIT_site_df["Set"] = LIT_site_df["Set"].str.replace("ESR", "ESR1")


for sample in folders:

    print("starting with sample: ", sample)
    sample_dir = LIT_dir + sample + "/"
    sample_files = get_files_in_dir(sample_dir)
    info[sample] = {}
    info[sample]["sample_dir"] = sample_dir
    info[sample]["pdbs"] = {}
    saved_uniprot_names = {}

    result_string = LIT_site_df[LIT_site_df["Set"].str.contains(sample)]["AID"].values[
        0
    ]
    print("Gene name: ", result_string)
    info[sample]["Gene_seq"] = aid_to_seq(result_string)

    # find mol2 files that contain the word "protein"
    for file in sample_files:
        if "protein" in file:
            # take the first 4 characters of the file name
            pdb_id = file[:4]

            if pdb_id not in info[sample]["pdbs"]:
                info[sample]["pdbs"][pdb_id] = {}
            PDB_seq = get_sequence_with_mol2(sample_dir + file)
            info[sample]["pdbs"][pdb_id]["PDB_seq"] = PDB_seq

            # ### Fetching the sequence from Uniprot ###
            # uniprot_name = pdbh.fetch_PDB_uniprot_accession(pdb_id)[0]
            # # if uniprot_name has a value, then fetch the sequence
            # if uniprot_name:
            #     # if the uniprot name is not in the saved_uniprot_names dictionary, then fetch the sequence
            #     if uniprot_name not in saved_uniprot_names:
            #         sequence = pdbh.fetch_uniprot_sequence(uniprot_name)
            #     else:
            #         sequence = saved_uniprot_names[uniprot_name]
            # else:
            #     sequence = None
            # info[sample]["pdbs"][pdb_id]["uniprot_seq"] = sequence

        if "ligand" in file:
            pdb_id = file[:4]
            print(pdb_id)
            # if the pdb_id is not in the dictionary, then add it
            if pdb_id not in info[sample]["pdbs"]:
                info[sample]["pdbs"][pdb_id] = {}
            info[sample]["pdbs"][pdb_id]["active_smiles"] = extract_smiles_from_mol2(
                sample_dir + file
            )

    info[sample]["decoys_smiles"] = extract_smiles_from_ism(sample_dir + "actives.smi")
    info[sample]["active_smiles"] = extract_smiles_from_ism(
        sample_dir + "inactives.smi"
    )
    print("finished with sample: ", sample)


# create json file for dictionary
import json

with open(LIT_base + "LIT_processed.json", "w") as fp:
    json.dump(info, fp, indent=4)
