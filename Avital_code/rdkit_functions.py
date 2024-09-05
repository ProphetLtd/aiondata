from rdkit import Chem
from rdkit.Chem import AllChem


def SDFfile_to_Smiles(sdf_file,outf):
	outf=open(outf,'w')
	outf.write(f'ligand\t Canonical SMILES \t Isomeric Smiles\t iso==canonic?\n')
	supplier = Chem.SDMolSupplier(sdf_file)
	for mol in supplier:
		if mol is not None:
			try:
				# Sanitize the molecule (checks for structural issues)
				Chem.SanitizeMol(mol)
				# Generate canonical and isomeric SMILES
				canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
				isomeric_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
				# Get the molecule name or ID (assuming the name is in the SDF's properties)
				name = mol.GetProp('_Name') if mol.HasProp('_Name') else 'Unknown'
				# Write the molecule info
				outf.write(f"{name}\t {canonical_smiles}\t {isomeric_smiles}\t" + str(canonical_smiles == isomeric_smiles)+"\n")

			except Exception as e:
				outf.write(f"Error processing molecule: {e}\n")
	# Close the output file
	outf.close()

SDFfile_to_Smiles('pfkfb3_automap_ligands.sdf','out_try')
