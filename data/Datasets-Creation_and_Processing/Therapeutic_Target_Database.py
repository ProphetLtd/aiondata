import requests
import pandas as pd


def getfile(url):
    response = requests.get(url)
    response.raise_for_status()  # Raise an exception for HTTP errors
    return response.text.splitlines(keepends=True)


def parse_lines_to_dataframe(lines, columns):
    # Initialize an empty list to collect data
    data_list = []

    # Initialize an empty dictionary to hold data for the current TARGETID
    current_data = {}

    # Iterate through each line in the list
    for line in lines:
        # Check for the end of the current TARGETID block (empty line)
        if line == "\r\n":
            # If there's data collected, add it to the list and reset the dictionary
            if current_data:
                data_list.append(current_data)
                current_data = {}
        else:
            # Split the line by tab characters and strip newline characters
            line = line.strip().split("\t")

            # Ensure the line has at least 3 elements
            if len(line) >= 3:
                # Extract the TARGETID, column name, and value
                target_id = line[0]
                column_name = line[1]
                value = line[2]

                # If it's a new TARGETID, set it
                if not current_data:
                    current_data["TARGETID"] = target_id

                # Add the value to the corresponding column in the dictionary
                current_data[column_name] = value

    # Add any remaining data (in case the last block doesn't end with a newline)
    if current_data:
        data_list.append(current_data)

    # Create a DataFrame from the collected data
    df = pd.DataFrame(data_list, columns=columns)

    return df


# Download the data files
# target URL
url = "https://idrblab.net/ttd/sites/default/files/ttd_database/P1-01-TTD_target_download.txt"
TargetLines = getfile(url)
# remove header lines and lines with DRUGINFO
TargetLines = [line for line in TargetLines[32:] if "DRUGINFO" not in line]

# drug URL
url = "https://idrblab.net/ttd/sites/default/files/ttd_database/P1-02-TTD_drug_download.txt"
DrugLines = getfile(url)[28:]

# target-compound URL
url = "https://idrblab.net/ttd/sites/default/files/ttd_database/P1-09-Target_compound_activity.txt"
Target_compound_df = pd.read_csv(url, sep="\t", header=0)


columns = [
    "TARGETID",
    "FORMERID",
    "UNIPROID",
    "TARGNAME",
    "GENENAME",
    "TARGTYPE",
    "SYNONYMS",
    "FUNCTION",
    "PDBSTRUC",
    "BIOCLASS",
    "ECNUMBER",
    "SEQUENCE",
]
Target_df = parse_lines_to_dataframe(TargetLines, columns)

columns = [
    "DRUG__ID",
    "TRADNAME",
    "DRUGCOMP",
    "THERCLAS",
    "DRUGTYPE",
    "DRUGINCH",
    "DRUGINKE",
    "DRUGSMIL",
    "HIGHSTAT",
    "COMPCLAS",
]
Drug_df = parse_lines_to_dataframe(DrugLines, columns)
# merge the two dataframes using "TTD Target ID" in the Target_compound_df and "TARGETID" in the Target_df
merged_df1 = pd.merge(
    Target_df, Target_compound_df, left_on="TARGETID", right_on="TTD Target ID"
)

# merge the merged_df1 with Drug_df using "DRUG__ID" in the Drug_df and "TTD Drug/Compound ID" in the merged_df1
merged_df2 = pd.merge(
    merged_df1, Drug_df, left_on="TTD Drug/Compound ID", right_on="DRUG__ID"
)

# save the merged dataframe to a csv file
merged_df2.to_csv("../local data/TTD/TTD_edited_data.csv", index=False)
