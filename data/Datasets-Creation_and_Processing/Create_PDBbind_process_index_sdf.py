import pandas as pd
import requests
import tarfile
import io
import os
from rdkit import Chem
from aiondata.raw import protein_structure
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def download_extract(link, path):
    response = requests.get(link)
    tar = tarfile.open(fileobj=io.BytesIO(response.content))
    tar.extractall(path=path)
    tar.close()


def sdf_to_smiles(sdf_path):
    """Convert a sdf file to a SMILES string."""
    # Load the molecule from an SDF file
    supplier = Chem.SDMolSupplier(sdf_path)
    molecules = [mol for mol in supplier if mol is not None]
    if molecules != []:
        for mol in molecules:
            isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
            canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
            # print(mol)
    else:
        isomeric = None
        canonical = None
    return isomeric, canonical


def process_folder_to_csv(folder_path, output_csv_path):
    """Convert all sdf files in a folder to a CSV with Isomeric and Canonical SMILES."""
    data = []

    for filename in os.listdir(folder_path):
        if filename.endswith(".sdf"):
            sdf_path = os.path.join(folder_path, filename)
            Isomeric, canonical = sdf_to_smiles(sdf_path)
            if Isomeric or canonical:
                # mol = Chem.MolFromSmiles(smiles)
                # if mol:
                # canonical_smiles = Chem.MolToSmiles(sdf_path, isomericSmiles=False)
                pdb_name = filename.split("_")[0]
                data.append(
                    {
                        "filename": filename,
                        "Isomeric SMILES": Isomeric,
                        "Canonical SMILES": canonical,
                        "PDB": pdb_name,
                    }
                )

    df = pd.DataFrame(data)
    df.to_csv(output_csv_path, index=False)


def read_and_clean_df(file_path):

    # Read the file content from line 7 onwards
    with open(file_path, "r") as file:
        lines = file.readlines()[6:]  # Skip the first six lines

    # split each line into columns
    data = [line.split() for line in lines if line.strip()]

    # Create a DataFrame
    df = pd.DataFrame(data)

    # remove column 5 and 8
    df = df.drop([5, 8], axis=1)

    column_names = [
        "PDB code",
        "resolution",
        "release year",
        "-logKd/Ki",
        "Kd/Ki",
        "reference",
        "ligand name",
    ]
    df.columns = column_names

    # remove the "(" and ")"" from the column "ligand name"
    df["ligand name"] = df["ligand name"].str.replace("(", "").str.replace(")", "")

    return df


def merge_and_clean(ps_df, df):
    # Merge DataFrames based on 'ligand name'
    merged_df = df.merge(ps_df, how="left", left_on="PDB code", right_on="PDB")

    # count the number of missing values in the column "SMILES"
    print(
        "Number of missing values in the column 'Isomeric SMILES': ",
        merged_df["Isomeric SMILES"].isnull().sum(),
    )

    # remove the rows with missing values in the column "SMILES"
    merged_df = merged_df.dropna(subset=["Isomeric SMILES"])

    # remove the column "Ligand HET ID in PDB"
    merged_df = merged_df.drop("PDB", axis=1)
    return merged_df


def fetch_sequence_for_pdb(pdb_code, pdbh):
    try:
        uni = pdbh.fetch_PDB_uniprot_accession(pdb_code)[0]
        if uni:
            seq = pdbh.fetch_uniprot_sequence(uni)
            return seq if seq else None
        else:
            return None
    except Exception as e:
        print(f"Error processing PDB ID {pdb_code}: {e}")
        return None


def add_pdb_based_sequence(df, pdb_column_name):
    pdbh = protein_structure.PDBHandler()
    sequences = []

    # Use ThreadPoolExecutor for parallelization
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Submit all tasks to the executor
        futures = {
            executor.submit(
                fetch_sequence_for_pdb, df[pdb_column_name].iloc[i], pdbh
            ): i
            for i in range(len(df))
        }

        # Use tqdm to track progress
        for future in tqdm(
            as_completed(futures), total=len(df), desc="Fetching Sequences"
        ):
            idx = futures[future]
            try:
                sequences.append(future.result())
            except Exception as e:
                print(f"Error at index {idx}: {e}")
                sequences.append(None)
    df["sequence"] = sequences
    return df


if __name__ == "__main__":
    # print current working directory
    print(os.getcwd())

    local_path = "../local data/PDBbind/"

    # Create folder called PDBbind in the local data folder

    # Check if the PDB index file exists
    if os.path.exists(local_path + "index/INDEX_general_PL_data.2020"):
        print("The PDB index file exists.")
    else:
        print("The PDB index file doesnt exist. Downloading...")
        PDB_index_link = (
            "http://www.pdbbind.org.cn/download/PDBbind_v2020_plain_text_index.tar.gz"
        )
        # Download and extract the PDB index file
        download_extract(PDB_index_link, local_path)

    # Check if the PDB sdf files exist
    if os.path.exists(local_path + "sdf/"):
        print("The PDB sdf files exist.")
    else:
        print("The PDB sdf files dont exist. Downloading...")
        PDB_sdf_link = "http://www.pdbbind.org.cn/download/PDBbind_v2020_sdf.tar.gz"

        # Download and extract the PDB sdf files
        download_extract(PDB_sdf_link)

    # Read clean and merge the data

    # Define the file path
    file_path = local_path + "index/INDEX_general_PL_data.2020"
    df = read_and_clean_df(file_path)

    ### Process the sdf files to a CSV file ###
    print("Processing sdf files")
    folder_path = "sdf"
    output_csv_path = local_path + "processed_smiles.csv"
    # process_folder_to_csv(local_path + folder_path, output_csv_path)

    # Read the processed SMILES CSV file
    ps_df = pd.read_csv(output_csv_path)

    # Merge the DataFrames
    merged_df = merge_and_clean(ps_df, df)

    # add sequence to the dataframe
    merged_df = add_pdb_based_sequence(merged_df, "PDB code")

    # save the merged DataFrame to a CSV file
    merged_df.to_csv(local_path + "PDBbind_SmilesAddition.csv", index=False)

    print("PDBbind edited has been saved to a CSV file.")
