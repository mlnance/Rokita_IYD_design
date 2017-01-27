#!/usr/bin/python
__author__="morganlnance"


import argparse
parser = argparse.ArgumentParser(description="Use PyRosetta to make low-energy mutants of pre-determined IYD design residues")
parser.add_argument("pdb_file", type=str, help="the path to the relevant PDB file")
input_args = parser.parse_args()


from rosetta import *
from rosetta.protocols.simple_moves import RotamerTrialsMover, \
    MinMover
from pyrosetta import *
from util import *
init()


# load in params to enable processing of IYD with cofactor and substrate
params = [ "/Users/Research/pyrosetta4/Rokita_IYD_design/params/2_iodophenol.params", "/Users/Research/pyrosetta4/Rokita_IYD_design/params/flavin_mononucleotide.params" ]
pose = Pose()
print "\ngenerate_nonstandard_residue_set"
nonstandard_res_set = generate_nonstandard_residue_set( pose, params )
print "\npose_from_file"
try:
    pose_from_file( pose, nonstandard_res_set, input_args.pdb_file )
except:
    print "\nThere was some error loading your PDB. Is this a valid PDB file?: %s\n" %input_args.pdb_file
    sys.exit()
pose.pdb_info().name( "IYD" )
# create a PyMOLMover
pmm = PyMOLMover()

# generate a score function
print "\nget_fa_scorefxn"
sf = get_fa_scorefxn()
