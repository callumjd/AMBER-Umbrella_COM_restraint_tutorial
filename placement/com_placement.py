#!/usr/bin/env python
import sys
import math
import os.path
import argparse
import numpy as np

################################################################################
# com_placement. Run as:
# ./com_placement.py -i bilayer.pdb -d drug_to_place.pdb -z point_along_z 
#
# IMPORTANT: Membrane bilayer center of mass is defined by N31 head group atoms.
# If your membrane does not contain these atom types, or you want to define
# using other/all atoms of the membrane, you will have to update the code.
#
################################################################################

# x,y,z cartesian coordinate object
class Coord(object):
	def __init__(self,x,y,z):
		self.x=x
		self.y=y
		self.z=z

# If your system has unusual atom types you may need to add the atomic
# mass into the dictionary below 
def atom_mass(str):
	weight={'C':'12.01','H':'1.008','P':'30.97','N':'14.01','O':'16.0','Cl':'35.45','F':'19.0','S':'32.065'}
	if str[0].isalpha() == True:
		if str[1].isalpha() == True:
			try:
				return weight[str[0:2]]
			except:
				return weight[str[0]]
		else:
			return weight[str[0]] 
	elif str[0].isalpha() == False:
		if str[1].isalpha() == True and str[2].isalpha() == True:
			try:
				return weight[str[1:3]]
			except:
				return weight[str[1]]
		else:
			return weight[str[1]]

################################################################################
#### Read input from command line ####
################################################################################
position=None
input_file_name=None
drug_file=None

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Membrane bilayer PDB file")
parser.add_argument("-d", type=str, help="Drug to position PDB file")
parser.add_argument("-z", type=float, help="Placement position along z-axis")

# Check command line arguments exist 
if len(sys.argv)<=4:
	parser.print_help()
	sys.exit(1)
elif len(sys.argv)>4:
	args = parser.parse_args()
	if (args.i != None and args.d):
		if (os.path.isfile(args.i) and os.path.isfile(args.d)):
			input_file_name=args.i
			drug_file=args.d
		elif not os.path.isfile(args.d):
			print('Cannot find file: %s' % (args.d))
			sys.exit(1)
		elif not os.path.isfile(args.i):
			print('Cannot find file: %s' % (args.i))
			sys.exit(1)
		elif (not os.path.isfile(args.i) and not os.path.isfile(args.d)):
			print('Cannot find inputs')
			sys.exit(1)
	else:
		parser.print_help()
		sys.exit(1)

	if args.z != None:
		position=float(args.z)
	else:
		#z-placement position was not specified, set to z=0
		position=0.0

if input_file_name==None or drug_file==None or position==None:
	print('Error: input options not set')
	sys.exit(1)

################################################################################
#### Membrane centre of mass ####
################################################################################
# Get number of atoms defining bilayer center-of-mass
length=0
with open(input_file_name, 'r') as f:
	for line in f:
		if line[0:3] != "TER" and line[0:3] != "END" and line[0:6] != "REMARK" and line[17:20] != "WAT" and line[17:20] != "Na+" and line[17:19] != "K+" and line[17:20] != "Cl-":
			if line[0:4] == "ATOM" or line[0:6] == "HETATM":
				if line[13:16] == "N31": #or line[13:16] == "P31":
					length+=1

# Check that bilayer PDB does contain N31 atom types as set above
if length==0:
	print('Error: Cannot set bilayer center-of-mass')
	sys.exit(1)

# Now set bilayer center-of-mass
i=0
bilayer=np.zeros((length,4))
with open(input_file_name, 'r') as f:
	for line in f:
		if line[0:3] != "TER" and line[0:3] != "END" and line[0:6] != "REMARK" and line[17:20] != "WAT" and line[17:20] != "Na+" and line[17:19] != "K+" and line[17:20] != "Cl-":
			if line[0:4] == "ATOM" or line[0:6] == "HETATM":
				if line[13:16] == "N31": #or line[13:16] == "P31":
					bilayer[i][0]=float(atom_mass(str(line[12:16])))
					bilayer[i][1]=str(line[30:38])
					bilayer[i][2]=str(line[38:46])
					bilayer[i][3]=str(line[46:54])
					i+=1

total_bilayer_weight=0.0
for i in range(0,length):
	total_bilayer_weight+=bilayer[i][0]

bilayer_com=Coord(0.000,   0.000,  0.000)
for i in range(0,length):
	bilayer_com.x+=bilayer[i][1]*bilayer[i][0]
	bilayer_com.y+=bilayer[i][2]*bilayer[i][0]
	bilayer_com.z+=bilayer[i][3]*bilayer[i][0]

bilayer_com.x=bilayer_com.x/total_bilayer_weight
bilayer_com.y=bilayer_com.y/total_bilayer_weight
bilayer_com.z=bilayer_com.z/total_bilayer_weight

################################################################################
#### Drug centre of mass ####
################################################################################
drug_length=0
with open(drug_file,'r') as f:
	for line in f:
		if line[0:3] != "TER" and line[0:3] != "END" and line[0:6] != "REMARK" and line[17:20]:
			if line[0:4] == "ATOM" or line[0:6] == "HETATM":
				drug_length+=1

i=0
drug=np.zeros((drug_length,4))
with open(drug_file,'r') as f:
	for line in f:
		if line[0:3] != "TER" and line[0:3] != "END" and line[0:6] != "REMARK" and line[17:20]:
			if line[0:4] == "ATOM" or line[0:6] == "HETATM":
				drug[i][0]=float(atom_mass(str(line[12:16])))
				drug[i][1]=str(line[30:38])
				drug[i][2]=str(line[38:46])
				drug[i][3]=str(line[46:54])
				i+=1

total_drug_weight=0.0
for i in range(0,drug_length):
	total_drug_weight+=drug[i][0]

drug_com=Coord(0.000,   0.000,  0.000)
for i in range(0,drug_length):
        drug_com.x+=drug[i][1]*drug[i][0]
        drug_com.y+=drug[i][2]*drug[i][0]
        drug_com.z+=drug[i][3]*drug[i][0]

drug_com.x=drug_com.x/total_drug_weight
drug_com.y=drug_com.y/total_drug_weight
drug_com.z=drug_com.z/total_drug_weight

move=Coord((bilayer_com.x-drug_com.x), (bilayer_com.y-drug_com.y), (bilayer_com.z-drug_com.z)+position)

################################################################################
# Print output to std out
################################################################################
with open(drug_file,'r') as f:
	for line in f:
	# Move
		if line[0:3] == "TER":
			print("TER")
		elif line[0:3] == "END":
                	print("END")
		elif (line[0:4] == "ATOM" or line[0:6] == "HETATM"):
			print("%s   %7.3f %7.3f %7.3f" % (line[0:28],float(line[30:38])+move.x,float(line[38:46])+move.y,float(line[46:54])+move.z))

	print("TER")
