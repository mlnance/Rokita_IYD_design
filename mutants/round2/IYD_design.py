#!/usr/bin/python
__author__="morganlnance"


import argparse
parser = argparse.ArgumentParser(description="Use PyRosetta to make low-energy mutants of pre-determined IYD design residues")
parser.add_argument("pdb_file", type=str, help="the path to the relevant PDB file")
parser.add_argument("params_dir", type=str, help="the path to the directory holding the relevant parameter files")
parser.add_argument("mutation_file", type=str, help="the path to the file containing the mutations to be made (PDB_num chain)")
input_args = parser.parse_args()


from rosetta import *
from rosetta.protocols.simple_moves import RotamerTrialsMover, \
    MinMover
from pyrosetta import *
from util import *
import pandas as pd
import sys, os
init()


# relevant mutation locations
with open( input_args.mutation_file, "rb" ) as fh:
    mutations = fh.readlines()
    # mutations = [ [ 101, A ], [ 105, B ] ]
    mutation_locations = [ mut.strip().split( ' ' ) for mut in mutations ]
# from Rokita email 9/29/16
#mutation_locations = [ 91, 95, 99, 103, 107, 112, 113, 116, 172]
# from Zuodong's PDB and finding residues within 5 Ang of 2IP
# this includes both chain A and B. Not doing symmetric mutations
#mutation_locations_A = [ 128, 129, 130, 131, 211, 212 ]
#mutation_locations_B = [ 104, 157, 161, 165, 169, 173, 176, 178, 184, 239 ]
#mutation_locations = []
#mutation_locations.extend( mutation_locations_A )
#mutation_locations.extend( mutation_locations_B )
# tester
#mutation_locations = [ 104 ]
#AA_name1_list = [ 'A' ]


# load in params to enable processing of IYD with cofactor and substrate
# organize the parameter directory and pull param files from there
if not input_args.params_dir.endswith( '/' ):
    params_dir = input_args.params_dir + '/'
else:
    params_dir = input_args.params_dir
params = [ os.path.join( params_dir, param ) for param in os.listdir( params_dir ) ]
print "\ngenerate_nonstandard_residue_set"
pose = Pose()
try:
    nonstandard_res_set = generate_nonstandard_residue_set( pose, params )
except:
    print "\nThere is something wrong with your parameter files. Did you give me the proper params directory?\n"
    sys.exit()
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

# store energy of original pose
orig_E = sf( pose )

# set up the lists that will hold relevant data
pdb_mutant_locations = []
pose_mutant_locations = []
chains = []
current_AAs = []
mutant_AAs = []
delta_energies = []
poses = []

# create a Pandas DataFrame
df = pd.DataFrame()


###################
#### MUTATIONS ####
###################
# for each relevant mutation location (given in PDB number)
for mutation in mutation_locations:
    # mutation_locations = [ [ 101, A ], [ 105, B ] ]
    # parse mutation file accordingly
    pdb_mutant_residue = int( mutation[0] )
    chain = mutation[1]
    # get the original residue and the Pose number
    pose_mutant_residue = pose.pdb_info().pdb2pose( chain, pdb_mutant_residue )
    orig_AA = pose.residue( pose_mutant_residue ).name1()
    # sample each single point mutant
    for mut_AA in AA_name1_list:
        # but pick the lowest E of three clone mutants
        best_mut_pose = None
        for ii in range( 3 ):
            print "\nmutate_residue"
            # we're giving a pose number
            mut_pose = Pose()
            mut_pose.assign( mutate_residue( pose_mutant_residue, 
                                             mut_AA, 
                                             pose, 
                                             sf, 
                                             pdb_num=False ) )
            # rename this mutant using the mutation information
            mut_pose.pdb_info().name( ''.join( [ orig_AA, str( pdb_mutant_residue ), mut_AA ] ) )

            # create a packer task to pack only the mutation site
            print "\nstandard_packer_task"
            task = standard_packer_task( mut_pose )
            task.or_include_current( True )
            task.restrict_to_repacking()
            [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in range( 1, mut_pose.size() + 1 ) if res_num != pose_mutant_residue ]
            print "\nRotamerTrialsMover"
            rtm = RotamerTrialsMover( sf, task )
            rtm.apply( mut_pose )

            # create a MoveMap to use in a MinMover for only the mutation site
            print "\nMoveMap"
            mm = MoveMap()
            mm.set_bb_true_range( pose_mutant_residue, pose_mutant_residue )
            mm.set_chi_true_range( pose_mutant_residue, pose_mutant_residue )
            print "\nMinMover"
            minmover = MinMover( movemap_in=mm,
                                 scorefxn_in=sf,
                                 min_type_in="lbfgs_armijo_nonmonotone",
                                 tolerance_in=0.01,
                                 use_nb_list_in=True,
                                 deriv_check_in=False,
                                 deriv_check_verbose_in=False )
            minmover.apply( mut_pose )

            # create a packer task to pack the whole protein
            print "\nstandard_packer_task"
            task = standard_packer_task( mut_pose )
            task.or_include_current( True )
            task.restrict_to_repacking()
            print "\nRotamerTrialsMover"
            rtm = RotamerTrialsMover( sf, task )
            rtm.apply( mut_pose )

            # create a MoveMap to use in a MinMover
            print "\nMoveMap"
            mm = MoveMap()
            mm.set_bb( True )
            mm.set_chi( True )
            mm.set_jump( True )
            print "\nMinMover"
            minmover = MinMover( movemap_in=mm,
                                 scorefxn_in=sf,
                                 min_type_in="lbfgs_armijo_nonmonotone",
                                 tolerance_in=0.01,
                                 use_nb_list_in=True,
                                 deriv_check_in=False,
                                 deriv_check_verbose_in=False )
            minmover.apply( mut_pose )

            # check if this is the best mutant of the three clones made
            # or, if this is the first mutant made of the three clones
            if best_mut_pose is None or ( sf( mut_pose ) < sf( best_mut_pose ) ):
                best_mut_pose = mut_pose.clone()

        # store the best mutation information in the Pandas DataFrame
        pdb_mutant_locations.append( pdb_mutant_residue )
        pose_mutant_locations.append( pose_mutant_residue )
        chains.append( chain )
        current_AAs.append( orig_AA )
        mutant_AAs.append( mut_AA )
        delta_energies.append( sf( best_mut_pose ) - orig_E )
        poses.append( best_mut_pose )


# append all mutation information into the Pandas Dataframe
df["pose_num"] = pose_mutant_locations
df["pdb_num"] = pdb_mutant_locations
df["chain"] = chains
df["orig_res"] = current_AAs
df["mutation"] = mutant_AAs
df["dG"] = delta_energies
df["poses"] = poses
df = df.sort( "dG" )
