from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def fasta_from_PDB(pdbf,out_fasta): #pdbf is a pdb in current dir
	# Initialize the PDB parser
	parser = PDB.PDBParser(QUIET=True)
	structure = parser.get_structure('protein', pdbf)
	# Extract the sequence from the first model
	ppb = PDB.PPBuilder()
	sequence = ""
	for pp in ppb.build_peptides(structure):
		sequence += str(pp.get_sequence())
	# Create a SeqRecord to write in FASTA format
	record = SeqRecord(Seq(sequence), id="Protein_from_PDB", description="FASTA sequence derived from PDB")
	# Save to FASTA format
	SeqIO.write(record, out_fasta, "fasta")
	print(f"FASTA sequence saved to {out_fasta}")


fasta_from_PDB('cmet_protein.pdb','cmet_fasta')
