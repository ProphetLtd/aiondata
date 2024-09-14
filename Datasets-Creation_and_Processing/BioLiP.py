import pandas as pd
import requests
import os

# File paths
txt_file = "../../local data/Biolip/biolip.txt"

# CSV file paths
full_csv_file = "../../local data/Biolip/biolip.csv"
clean_csv_file = "../../local data/Biolip/biolip_clean.csv"


#### Download the data from the URL if not found locally ####

# URL to download the txt file if not found locally
url = "https://zhanggroup.org/BioLiP/qsearch.cgi?outfmt=txt&order=pdbid"

# Check if the txt file exists locally
if os.path.exists(txt_file):
    print(f"Loading data from local file: {txt_file}")
    with open(txt_file, "r") as file:
        data = file.read()
else:
    print(f"Downloading data from {url}")
    # Fetch the content from the URL
    response = requests.get(url)
    data = response.text
    # Save the downloaded text to a file for future use
    with open(txt_file, "w") as file:
        file.write(data)
    print(f"Data saved to local file: {txt_file}")


#### Parse the data and create a DataFrame ####


# Convert the text data to a list of lines
lines = data.splitlines()

# Define the column headers as specified
columns = [
    "PDB ID",
    "Receptor chain",
    "Resolution",
    "Binding site number code",
    "Ligand ID",
    "Ligand chain",
    "Ligand serial number",
    "Binding site residues (PDB numbering)",
    "Binding site residues (renumbered)",
    "Catalytic site residues (PDB numbering)",
    "Catalytic site residues (renumbered)",
    "EC number",
    "GO terms",
    "Binding affinity (manual)",
    "Binding affinity (MOAD)",
    "Binding affinity (PDBbind-CN)",
    "Binding affinity (BindingDB)",
    "UniProt ID",
    "PubMed ID",
    "Ligand residue sequence number",
    "Receptor sequence",
]

# Initialize a list to store data
data_list = []

# Iterate over the lines and split each line by tab
for line in lines:
    split_line = line.split("\t")
    if len(split_line) == len(columns):  # Ensure the row matches expected column count
        data_list.append(split_line)

# Create a DataFrame from the data list
df = pd.DataFrame(data_list, columns=columns)

# Filter the DataFrame to only include rows with binding affinity data
df_ba = df[
    (df["Binding affinity (manual)"] != "")
    | (df["Binding affinity (MOAD)"] != "")
    | (df["Binding affinity (PDBbind-CN)"] != "")
    | (df["Binding affinity (BindingDB)"] != "")
]


#### Save the DataFrame to a CSV file####

# Save the DataFrame to a CSV file
df.to_csv(full_csv_file, index=False)

print(f"CSV file has been created: {full_csv_file}")

# Save the DataFrame to a CSV file
df_ba.to_csv(clean_csv_file, index=False)

print(f"CSV file has been created: {clean_csv_file}")
