#!/usr/bin/python
__author__="morganlnance"


import argparse
parser = argparse.ArgumentParser(description="Use PyRosetta to make a low-energy base structure of bact_IYD_chainA.pdb")
parser.add_argument("rounds", type=int, help="the number of base poses you want to generate to find the lowest energy PDB")
input_args = parser.parse_args()


# imports
from rosetta import *
from rosetta.protocols.simple_moves import RotamerTrialsMover, \
    MinMover
from pyrosetta import *
from util import *
import pandas as pd
init()



#################################
#### LOAD IN PDB WITH PARAMS ####
#################################
# load in params to enable processing of IYD with cofactor and substrate
params = [ "/Users/Research/pyrosetta4/Rokita_IYD_design/params/2_iodophenol.params", "/Users/Research/pyrosetta4/Rokita_IYD_design/params/flavin_mononucleotide.params" ]
orig_pose = Pose()
print "\ngenerate_nonstandard_residue_set"
nonstandard_res_set = generate_nonstandard_residue_set( orig_pose, params )
print "\npose_from_file"
#pose_from_file( orig_pose, nonstandard_res_set, "/Users/Research/pyrosetta4/Rokita_IYD_design/bact_IYD_chainA.pdb" )
pose_from_file( orig_pose, nonstandard_res_set, "/Users/Research/pyrosetta4/Rokita_IYD_design/test.pdb" )
orig_pose.pdb_info().name( "IYD" )

# generate a score function
print "\nget_fa_scorefxn"
sf = get_fa_scorefxn()
orig_E = sf( orig_pose )

# create a Pandas DataFrame and data holders
df = pd.DataFrame()
orig_E_vals = []
base_E_vals = []
dG_vals = []
poses = []


###########################
#### PACK AND MINIMIZE ####
###########################
# create a packer task
print "\nstandard_packer_task"
task = standard_packer_task( orig_pose )
task.or_include_current( True )
task.restrict_to_repacking()

# create a MoveMap used for minimization
print "\nMoveMap"
mm = MoveMap()
mm.set_bb( True )
mm.set_chi( True )
mm.set_jump( True )

# for as many decoy poses as specified
for ii in range( 1, input_args.rounds + 1 ):
    # get a fresh copy of the original pose
    pose = orig_pose.clone()
    pose.pdb_info().name( "IYD_" + str( ii ) )

    # create a RotamerTrialsMover
    rtm = RotamerTrialsMover( sf, task )

    # create a MinMover using the MoveMap
    minmover = MinMover( movemap_in=mm, 
                         scorefxn_in=sf, 
                         min_type_in="lbfgs_armijo_nonmonotone", 
                         tolerance_in=0.01, 
                         use_nb_list_in=True, 
                         deriv_check_in=False, 
                         deriv_check_verbose_in=False )

    # run multiple RotamerTrialsMover and MinMover rounds
    # 1
    print "\nRotamerTrialsMover"
    rtm.apply( pose )
    print "\nMinMover"
    minmover.apply( pose )
    # 2
    print "\nRotamerTrialsMover"
    rtm.apply( pose )
    print "\nMinMover"
    minmover.apply( pose )
    # 3
    print "\nRotamerTrialsMover"
    rtm.apply( pose )
    print "\nMinMover"
    minmover.apply( pose )

    # relay energy of pose
    base_E = sf( pose )
    print "\n\norig energy:", orig_E
    print "current energy:", base_E

    # store pose and energies in DataFrame lists
    orig_E_vals.append( orig_E )
    base_E_vals.append( base_E )
    dG_vals.append( base_E - orig_E )
    poses.append( pose )

# add the data to the DataFrame
df["orig_E"] = orig_E_vals
df["base_E"] = base_E_vals
df["dG"] = dG_vals
df["poses"] = poses
df = df.sort( "dG" )
