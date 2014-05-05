#!/usr/bin/python

# Isidro Cortes Ciriano.	6/8/2013
# Institut Pasteur
# isidrolauscher@gmail.com

# Import modules
import argparse
import numpy as np
import os,sys
# Arguments passed to the scripts
parser = argparse.ArgumentParser(prog='PROG',description='Get Morgan Fingerprints for compounds codified in either SMILES or SDF format using RDkit. Isidro Cortes Ciriano. August/September 2013')
parser.add_argument('--bits', required='TRUE',type=int, help="Size of the hashed Morgan Fingerprints (binary and with counts)")
parser.add_argument('--rad', required='TRUE', type=int, help="Maximum radius. Deafault is two, equivalent to ECFP4 from PipelinePilot")
parser.add_argument('--f', required='TRUE', type=str, help="Format of the input file")
parser.add_argument('--mols', type=str,help="File containing the molecules {.smi|.smiles|.sdf}. If the format is smiles, each line should contain the smiles and the name separated by a comma (in this order)")
parser.add_argument('--image', action='store_true', help="Write --image if you want the images of the substructures")
parser.add_argument('--unhashed', action='store_true', help="Write --unhashed if you want the unhashed fingerprints")
parser.add_argument('--v', action='store_true', help="Verbose")
parser.add_argument('--extF',type=str, help="Type -extF followed by the format {.smi|.smiles|.sdf} of the external file for which you want to calculate HASHED circular fingerprints")
parser.add_argument('--molsEXT',type=str,help="External file")
parser.add_argument('--unhashedEXT', action='store_true', help="Write --unhashedEXT if you want the unhashed fingerprints for the external file. The substructures of the molecules in the external file will be compared to the pool of substructures contained in the molecules of the main file")
parser.add_argument('--RDkitPath', required='TRUE', type=str, help="Path to the directory where the RDkit files are")
parser.add_argument('--output', required='TRUE', type=str, help="Name of the output files")
args = vars(parser.parse_args())

image=args['image']
unhashed=args['unhashed']
verbose=args['v']
formatFile=args['f']
fileMols=str(args['mols'])
nbBits=int(args['bits'])
fp_diam=int(args['rad'])
# External file.
formatFileEXT=args['extF']
fileMolsEXT=str(args['molsEXT'])
unhashedEXT=args['unhashedEXT']
RDkitPath=args['f']
outname=args['output']
sys.path.append(RDkitPath)

if verbose:
	if image:
		print "\nCalculation of Morgan Fingerprints with diameter %d hashed into a fingerprint size equal to %d.\nMolecules file: %s.\nImages for the chemical substructures will be created.\n" %(args['rad'],args['bits'],args['mols'])
	else :
		print "\nCalculation of Morgan Fingerprints with diameter %d hashed into a fingerprint size equal to %d.\nMolecules file: %s.\n NO Images for the chemical substructures will be created.\n" %(args['rad'],args['bits'],args['mols'])

#####################################
# Import Modules
#####################################
import gzip
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.rdBase
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem import AllChem
from rdkit import DataStructs 
from rdkit.DataStructs import BitVectToText
from rdkit.Chem import Draw

#####################################
# Define Functions
#####################################
# To search within sublists:
def insublist(item, list):
    for l in list:
        if np.array_equal(item,l):
            return True
    return False

# Define a function to create matrices of empty strings
def nans(shape, dtype=str):
    a = np.empty(shape, dtype)
    a[:]=""
    return a

### Read Mol2 files
def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
	import rdkit.Chem
	"""generator which retrieves one mol2 block at a time
	"""
	mol2 = []
	for line in fileLikeObject:
		if line.startswith(delimiter) and mol2:
			yield "".join(mol2)
			mol2 = []
		mol2.append(line)
	if mol2:
		yield "".join(mol2)


#####################################
# Read Molecules
#####################################
if formatFile == 'smi' or formatFile == 'smiles':
	if verbose:
		print "Format of the main file = SMILES"
	suppl = Chem.SmilesMolSupplier(fileMols,smilesColumn=0,nameColumn=1,delimiter=',',titleLine=False)
	mols=[]
	molserr=[]
	for i,m in enumerate(suppl):
		if m is not None:
			mols.append(m)
		else:
			molserr.append(i)
	nbMols=len(mols)


elif formatFile == 'mol2':
	molss=[]

	with open(fileMols) as fi:
		for mol2 in RetrieveMol2Block(fi):
			rdkMolecule = rdkit.Chem.MolFromMol2Block(mol2)
			molss.append(rdkMolecule)
	
	molserr=[]
	mols=[]
	for i,m in enumerate(molss):
		if m is not None:
			mols.append(m)
		else:
			molserr.append(i)
			mols.append(m)   #### XXXX 
	nbMols=len(mols)
	print nbMols


else:
	if verbose:
		print "Format of the main file = SDF"
	suppl = Chem.SDMolSupplier(fileMols)
	mols=[]
	molserr=[]
	for i,m in enumerate(suppl):
		if m is not None:
			mols.append(m)
		else:
			#molserr.append(m.GetProp('_Name'))
			molserr.append(i)
	nbMols=len(mols)



#mols2 = RetrieveMol2Block(fileMols)
#mols3 = [ x for x in mols2]
#print len(mols3)




if verbose: 
	if len(molserr) !=0:
		print "The following %d molecules (starting at zero) could not be processed:\n"%(len(molserr))
		for x in molserr: print x
		print "NOTE: the indexes of the molecules start at zero. Thus the first molecule is molecule 0."
		errfile="incorrect_molecules_"+outname+".csv"
		print "This information has been saved in the following file: %s\n"%(errfile)
		# Save the information about which molecules could not be processed correctly.
		np.savetxt(errfile,molserr,fmt="%d")
		del errfile
	else:
		print "All molecules in the input file were processed correctly"


	if verbose:
		print 'Your molecules file has %d molecules\n' % (len(mols))
		if formatFileEXT:
			print 'Your external file contains %d molecules\n' % (len(molsEXT))



# Define the list that will contain all the submolecules
subm_all=[]
# Define the list where the erroneous molecules will be saved.
err_mols=[]

# Define a progess bar
from progressbar import *
widgets = ['Progression: ', Percentage(), ' ', Bar(marker='.',left='[',right=']'),' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

pbar = ProgressBar(widgets=widgets, maxval=nbMols)
pbar.start()


# Get the total number of features:

nbFeatTot=0
ids_all_feat=[]
smiles_subs_kept =[]
kept_mols=[]
name_mols=[]
#Loop over the molecules
for molecule_nb,m in enumerate(mols):
	info2={}
	if m is None:
			print "Erroneous input at molecule: %d" %(molecule_nb)
			err_mols.append(molecule_nb)
			#print m.GetProp('_Name')
	else:
		kept_mols.append(m)
		fp2 = AllChem.GetMorganFingerprint(m,fp_diam,bitInfo=info2)
		ids_now=np.asarray([info2.items()[i][0] for i in range(0,len(info2.items()))])
		for i in range(0,len(info2.items())):
			radius=info2.items()[i][1][0][1]
			atom=info2.items()[i][1][0][0]
			env=Chem.FindAtomEnvironmentOfRadiusN(m,radius,atom)
			amap={}
			submol=Chem.PathToSubmol(m,env,atomMap=amap)
			if submol.GetNumAtoms() >1: # This means that the radius is zero, so the feature is a single atom
				ids_all_feat.append(ids_now[i])
				smiles_subs_kept.append(Chem.MolToSmiles(submol))
#				print 'pep', Chem.MolToSmiles(submol)
#			else:
#				print 'pep', Chem.MolToSmiles(submol)

#print 'smiles_feat'
#print smiles_subs_kept
#print 'feat ids\n'
#print ids_all_feat

#ids_all_feat =  np.unique(np.array(ids_all_feat))
ids_all_feat=np.array(ids_all_feat)

fp_matrix = np.zeros((nbMols,ids_all_feat.shape[0]))
#print 'shape fps matrix',fp_matrix.shape
coordinates = np.zeros((nbMols,ids_all_feat.shape[0]*3))
#print 'shape Matrix', coordinates.shape

# now for each submol, we have to get the coordinates..
# and also the unhashed fps.. 

# we are going to keep the data in different lists:
# compound ids:
comp_ids=[]
# feat_id
feat_ids_expanded=[]
#xcoordinates
x_coords=[]
# ycoordinates
y_coords=[]
# zcoordinates
z_coords=[]
smiles_all=[]
smiles_whole_mol=[]

nbFeatTot=0
#Loop over the molecules
for molecule_nb,m in enumerate(kept_mols):
	info={}; info2={}
	if m is None:
			print "Erroneous input at molecule: %d" %(molecule_nb)
			err_mols.append(molecule_nb)
	else:

		fp2 = AllChem.GetMorganFingerprint(m,fp_diam,bitInfo=info2)

		for i in range(0,len(info2.items())): # iterate over the substructures
			info3 = [ (key, val) for key, value in info2.iteritems() for val in value ]
			radius=info3[i][1][1]
			atom=info3[i][1][0]
			feat=info3[i][0] # substructure ID
			if feat in ids_all_feat:
				smiles_whole_mol.append(Chem.MolToSmiles(m))
				name_mols.append(m.GetProp('_Name'))
				comp_ids.append(molecule_nb)
				feat_ids_expanded.append(feat)
				env=Chem.FindAtomEnvironmentOfRadiusN(m,radius,atom)
				amap={}
				submol=Chem.PathToSubmol(m,env,atomMap=amap)
				image_name="./images/%s_Feature_%s_%d.pdf"%(outname,m.GetProp('_Name'),feat)
				Draw.MolToFile(m,image_name,size=(300,300),highlightAtoms=amap.keys())
				smiles_all.append(Chem.MolToSmiles(submol))
#				print 'ID',feat,'  smiles ',Chem.MolToSmiles(submol)
				atoms=amap.keys()
				xcoords = [m.GetConformer().GetAtomPosition(x).x for x in atoms]
				ycoords = [m.GetConformer().GetAtomPosition(x).y for x in atoms]
				zcoords = [m.GetConformer().GetAtomPosition(x).z for x in atoms]
				x_coords.append(sum(xcoords) / len(xcoords))
				y_coords.append(sum(ycoords) / len(ycoords))
				z_coords.append(sum(zcoords) / len(zcoords))


print 'Comps',len(comp_ids), 'feats',len(feat_ids_expanded),'\n' 
##len(xcoords), len(ycoords), len(zcoords)
# write the results to a file
filename = outname+"_results.txt"
f = open(filename,'w')
data='Compound_name,Compound_nb,Feature_ID,X_Coord,Y_Coord,Z_Coord,Smiles_Feat,Smiles_Molecule\n'
f.write(data)
for i in range(0,len(smiles_all)):
	data='%s,%d,%d,%f,%f,%f,%s,%s\n'%(name_mols[i],comp_ids[i],feat_ids_expanded[i],x_coords[i],y_coords[i],z_coords[i],smiles_all[i],smiles_whole_mol[i])
	f.write(data)
f.close()

print 'total number of features',len(smiles_all)
print 'total nb of molecules',len(kept_mols)
