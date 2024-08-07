import pandas as pd
import requests
import tarfile
import io

local_path="./"
#create folder called PDBbind in the local data folder


PDB_index_link="http://www.pdbbind.org.cn/download/PDBbind_v2020_plain_text_index.tar.gz"

def download_extract(link, path): 
    response = requests.get(link)
    tar = tarfile.open(fileobj=io.BytesIO(response.content))
    tar.extractall(path=path)
    tar.close()

# Download and extract the PDB index file
download_extract(PDB_index_link, local_path)


# Define the file path
file_path = local_path+"index/INDEX_general_PL_data.2020"

# Read the file content from line 7 onwards
with open(file_path, 'r') as file:
    lines = file.readlines()[6:]  # Skip the first six lines

#split each line into columns
data = [line.split() for line in lines if line.strip()]

# Create a DataFrame
df = pd.DataFrame(data)

#remove column 5 and 8
df = df.drop([5, 8], axis=1)

column_names=["PDB code", "resolution", "release year", "-logKd/Ki", "Kd/Ki", "reference", "ligand name"]
df.columns = column_names

#remove the "(" and ")"" from the column "ligand name"
df["ligand name"] = df["ligand name"].str.replace("(", "").str.replace(")", "")



# # Display the DataFrame
df

# Download sdf_link and open the tar gz file
import requests
import tarfile
import io

sdf_link="https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_sdf.tar.gz"
download_extract(sdf_link, local_path)
sdffolder=local_path+"sdf/"



import os
import pandas as pd
from rdkit import Chem

def sdf_to_smiles(sdf_path):
    """Convert a sdf file to a SMILES string."""
    # Load the molecule from an SDF file
    supplier = Chem.SDMolSupplier(sdf_path)
    molecules = [mol for mol in supplier if mol is not None]
    if molecules !=[]:
    	for mol in molecules:
    		sm=Chem.MolToSmiles(mol,isomericSmiles=True)
    		print (mol)
    else:
    	sm=None
#    if mol is None:
 #   	return None
    return sm

def process_folder_to_csv(folder_path, output_csv_path):
    """Convert all sdf files in a folder to a CSV with Isomeric and Canonical SMILES."""
    data = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            sdf_path = os.path.join(folder_path, filename)
            smiles = sdf_to_smiles(sdf_path)
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                    pdb_name = filename.split('_')[0]
                    data.append({
                        'filename': filename,
                        'Isomeric SMILES': smiles,
                        'Canonical SMILES': canonical_smiles,
                        'PDB': pdb_name
                    })
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv_path, index=False)

# Usage
folder_path = sdffolder
output_csv_path = local_path+'processed_smiles.csv'
process_folder_to_csv(folder_path, output_csv_path)

# Read the processed SMILES CSV file
ps_df = pd.read_csv(output_csv_path)

# Merge DataFrames based on 'ligand name'
merged_df = df.merge(ps_df, how='left', left_on='PDB code', right_on='PDB')

# Display the merged DataFrame
merged_df

#count the number of missing values in the column "SMILES"
print(merged_df["Isomeric SMILES"].isnull().sum())

#remove the rows with missing values in the column "SMILES"
merged_df = merged_df.dropna(subset=["Isomeric SMILES"])

#remove the column "Ligand HET ID in PDB"
merged_df = merged_df.drop("PDB", axis=1)

#save the merged DataFrame to a CSV file
merged_df.to_csv(local_path+"PDBbind_SmilesAddition.csv", index=False)

