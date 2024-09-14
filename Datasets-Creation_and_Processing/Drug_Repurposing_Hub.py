import pandas as pd

url = "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt"
rep_df = pd.read_csv(url, sep="\t", header=9)


url = "https://storage.googleapis.com/cdot-general-storage/repurposing_samples_20240610.txt"
sample_df = pd.read_csv(url, sep="\t", header=9)

# merge the two dataframes on the "pert_iname" column
merged_df = pd.merge(rep_df, sample_df, on="pert_iname")

# explode the "target" column to get a row for each target
merged_df2 = merged_df.assign(target=merged_df["target"].str.split("|")).explode(
    "target"
)

# save the dataframe to a csv file
merged_df2.to_csv(
    "../local data/Drug_Repurposing_Hub/Drug_Repurposing_Hub_databased_edited.csv",
    index=False,
)
