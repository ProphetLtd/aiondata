import os
import pandas as pd

def collect_pocket(pdb_chain):
	pockets={}
	pockets_path=os.path.join(pdb_chain+'_out','pockets')
	# make a list of the files with pocket info in the aove dir
	if os.path.exists(pockets_path):
		pocket_files= [file for file in os.listdir(pockets_path)]
		for filename in pocket_files:
			if filename.endswith('.pdb'):
				file_path = os.path.join(pockets_path, filename)	
		# read the lines
				lines=open(file_path).readlines()		
				resisudes=[]
		# find Pocketscore Drug score and residues that compose the pocket from the pocket file
				Dscore= float([line for line in lines if 'Drug Score' in line][0].split(':')[1])
				Pscore= float([line for line in lines if 'Pocket Score' in line][0].split(':')[1])	
				Aspheres = float([line for line in lines if 'Number of alpha spheres' in line][0].split(':')[1])
				Vol = float([line for line in lines if 'Pocket volume (convex hull)' in line][0].split(':')[1])
			# Define a function to extract the leading number from a string
				def extract_number(s):
					return int(''.join(filter(str.isdigit, s.split('_')[0])))
			# use extract_num to sort the residues and also remove duplicates:
				residues = list(dict.fromkeys(sorted([line.split()[5]+'_'+line.split()[3] for line in lines if 'ATOM' in line],key=extract_number)))
				pocket_num= filename.split('_')[0]
				pockets[pocket_num]=Pscore,Dscore,residues,Aspheres,Vol
	else:
		pocket_num,Pscore,Dscore,residues,Aspheres,Vol=None,None,None,[],None,None
	return  pockets

# funct to calculate the overlap between two pockets in fpocket csv file
def pocket_overlap(list_pocket1,list_pocket2):
	# Convert lists to sets and compute the intersection
	overlap = set(list_pocket1) & set(list_pocket2)
	# Convert the result back to a list if needed
	overlap1,overlap2 = len(list(overlap))/len(list_pocket1),len(list(overlap))/len(list_pocket2)
	return round(overlap1,2),round(overlap2,2)
