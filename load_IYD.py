#!/usr/bin/python
__author__="morganlnance"


import argparse
parser = argparse.ArgumentParser(description="Use PyRosetta to make low-energy mutants of pre-determined IYD design residues")
parser.add_argument("pdb_file", type=str, help="the path to the relevant PDB file")
parser.add_argument("params_dir", type=str, help="the path to the directory holding the relevant parameter files")
input_args = parser.parse_args()


from rosetta import *
from rosetta.protocols.simple_moves import RotamerTrialsMover, \
    MinMover
from pyrosetta import *
from util import *
import sys, os
init()


# load in params to enable processing of IYD with cofactor and substrate
# organize the parameter directory and pull param files from there
if not input_args.params_dir.endswith( '/' ):
    params_dir = input_args.params_dir + '/'
else:
    params_dir = input_args.params_dir
params = [ os.path.join( params_dir, param ) for param in os.listdir( params_dir ) ]
print "\ngenerate_nonstandard_residue_set"
in_pose = Pose()
try:
    nonstandard_res_set = generate_nonstandard_residue_set( in_pose, params )
except:
    print "\nThere is something wrong with your parameter files. Did you give me the proper params directory?\n"
    sys.exit()
print "\npose_from_file"
try:
    pose_from_file( in_pose, nonstandard_res_set, input_args.pdb_file )
except:
    print "\nThere was some error loading your PDB. Is this a valid PDB file?: %s\n" %input_args.pdb_file
    sys.exit()
in_pose.pdb_info().name( "IYD" )
# create a PyMOLMover
pmm = PyMOLMover()

# generate a score function
print "\nget_fa_scorefxn"
sf = get_fa_scorefxn()
