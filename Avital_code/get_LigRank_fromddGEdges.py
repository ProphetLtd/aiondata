import pandas as pd
from collections import Counter


def rank_from_edges(input_edges_csv,out_rank):
	df = pd.read_csv(input_edges_csv, header=None, skiprows=1) # assuming file has the columns Lig1,Lig2,ddG. ignore current header.
	df.columns=['Ligand 1' , 'Ligand 2', 'ddG (kcal/mol)'] # and replace with htis header
	### choose the reference ligand !
	# Count occurrences of each ligand
	ligand1_counts = df['Ligand 1'].value_counts()
	ligand2_counts = df['Ligand 2'].value_counts()
	# Combine counts from both columns
	ligand_counts = ligand1_counts.add(ligand2_counts, fill_value=0)
	# Choose the ligand with the highest count as the reference
	reference_ligand = ligand_counts.idxmax()
	print(f"Reference Ligand: {reference_ligand}")

	# Initialize dictionaries to store energy differences
	energy_dict = {}
	visited = set()

	# Function to compute energy differences
	def compute_energy_differences(current_ligand):
		if current_ligand in visited:
			return
		visited.add(current_ligand)
		for index, row in df.iterrows():
			ligand1 = row['Ligand 1']
			ligand2 = row['Ligand 2']
			ddG = float(row['ddG (kcal/mol)'])
			if ligand1 == current_ligand:
				if ligand2 not in energy_dict:
					energy_dict[ligand2] = energy_dict[current_ligand] + ddG
				else:
					energy_dict[ligand2] = min(energy_dict[ligand2], energy_dict[current_ligand] + ddG)
				compute_energy_differences(ligand2)
			elif ligand2 == current_ligand:
				if ligand1 not in energy_dict:
					energy_dict[ligand1] = energy_dict[current_ligand] - ddG
				else:
					energy_dict[ligand1] = min(energy_dict[ligand1], energy_dict[current_ligand] - ddG)
				compute_energy_differences(ligand1)

	# Initialize reference ligand
	energy_dict[reference_ligand] = 0
	compute_energy_differences(reference_ligand)

	# Handle unvisited ligands
	unvisited_ligands = set(df['Ligand 1']).union(df['Ligand 2']) - visited
	for ligand in unvisited_ligands:
    		energy_dict[ligand] = float('inf')  # Assign a high value to unvisited ligands

	# Create output file with ranked ligands
	with open(out_rank, 'w') as f:
		sorted_ligands = sorted(energy_dict.items(), key=lambda x: x[1])
		for ligand, score in sorted_ligands:
			score=round(score,3)
			f.write(f"{ligand}, {score}\n")

	print(f"Ranking complete. Results saved to {out_rank }")
