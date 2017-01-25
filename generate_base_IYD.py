#!/usr/bin/python
__author__="morganlnance"

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
pose_from_file( pose, nonstandard_res_set, "/Users/Research/pyrosetta4/Rokita_IYD_design/bact_IYD_chainA.pdb" )
pose.pdb_info().name( "IYD" )

# copy the input pose
orig_pose = pose.clone()
orig_pose.pdb_info().name( "orig_IYD" )

# generate a score function
print "\nget_fa_scorefxn"
sf = get_fa_scorefxn()

# create a packer task
print "\nstandard_packer_task"
task = standard_packer_task( pose )
task.or_include_current( True )
task.restrict_to_repacking()

# create and apply a RotamerTrialsMover
print "\nRotamerTrialsMover"
rtm = RotamerTrialsMover( sf, task )
rtm.apply( pose )

# create a MoveMap used for minimization
print "\nMoveMap"
mm = MoveMap()
mm.set_bb( True )
mm.set_chi( True )
mm.set_jump( True )

# create and apply a MinMover using the MoveMap
print "\nMinMover"
minmover = MinMover( movemap_in=mm, 
                     scorefxn_in=sf, 
                     min_type_in="lbfgs_armijo_nonmonotone", 
                     tolerance_in=0.01, 
                     use_nb_list_in=True, 
                     deriv_check_in=False, 
                     deriv_check_verbose_in=False )
minmover.apply( pose )

# run more RotamerTrialsMover and MinMover rounds
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

print "\n\norig energy:", sf( orig_pose )
print "current energy:", sf( pose )
