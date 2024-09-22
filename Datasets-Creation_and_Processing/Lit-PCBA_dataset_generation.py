# read all folder or file names in a directory
import os
from Bio.Data import IUPACData
from rdkit import Chem
import tempfile
from rdkit import Chem
from rdkit import RDLogger
from aiondata.raw import protein_structure

pdbh = protein_structure.PDBHandler()


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


# Suppress RDKit warnings
RDLogger.DisableLog("rdApp.*")


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


def get_inactive_file(file_path):
    """Reads a text file and returns a list of lines."""
    with open(file_path, "r") as file:
        lines = file.readlines()
    return [line.split()[0] for line in lines]


#### Main ####


LIT_base = "local data/LIT-PCBA/"
LIT_dir = LIT_base + "PDBs/"
info = {}
folders = get_folder_names(LIT_dir)


for sample in folders:

    print("starting with sample: ", sample)
    sample_dir = LIT_dir + sample + "/"
    sample_files = get_files_in_dir(sample_dir)
    info[sample] = {}
    info[sample]["sample_dir"] = sample_dir
    info[sample]["pdbs"] = {}
    saved_uniprot_names = {}

    # find mol2 files that contain the word "protein"
    for file in sample_files:
        if "protein" in file:
            # take the first 4 characters of the file name
            pdb_id = file[:4]

            if pdb_id not in info[sample]["pdbs"]:
                info[sample]["pdbs"][pdb_id] = {}
            PDB_seq = get_sequence_with_mol2(sample_dir + file)
            info[sample]["pdbs"][pdb_id]["PDB_seq"] = PDB_seq

            ### Fetching the sequence from Uniprot ###
            uniprot_name = pdbh.fetch_PDB_uniprot_accession(pdb_id)[0]
            # if uniprot_name has a value, then fetch the sequence
            if uniprot_name:
                # if the uniprot name is not in the saved_uniprot_names dictionary, then fetch the sequence
                if uniprot_name not in saved_uniprot_names:
                    sequence = pdbh.fetch_uniprot_sequence(uniprot_name)
                else:
                    sequence = saved_uniprot_names[uniprot_name]
            else:
                sequence = None
            info[sample]["pdbs"][pdb_id]["uniprot_seq"] = sequence

        if "ligand" in file:
            pdb_id = file[:4]
            print(pdb_id)
            # if the pdb_id is not in the dictionary, then add it
            if pdb_id not in info[sample]["pdbs"]:
                info[sample]["pdbs"][pdb_id] = {}
            info[sample]["pdbs"][pdb_id]["active_smiles"] = extract_smiles_from_mol2(
                sample_dir + file
            )

    info[sample]["decoys_smiles"] = get_inactive_file(sample_dir + "inactives.smi")
    print("finished with sample: ", sample)


# create json file for dictionary
import json

with open(LIT_base + "LIT_processed.json", "w") as fp:
    json.dump(info, fp, indent=4)
