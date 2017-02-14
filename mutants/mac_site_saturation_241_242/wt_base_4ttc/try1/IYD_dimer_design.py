#!/usr/bin/python
__author__="morganlnance"


import argparse
parser = argparse.ArgumentParser(description="Use PyRosetta to make low-energy mutants of pre-determined IYD design residues")
parser.add_argument("pdb_file", type=str, help="the path to the relevant PDB file")
parser.add_argument("params_dir", type=str, help="the path to the directory holding the relevant parameter files")
parser.add_argument("mutation_file", type=str, help="the path to the file containing the mutations to be made (PDBnum1_PDBnum2_etc)")
input_args = parser.parse_args()


from rosetta import *
from rosetta.protocols.simple_moves import RotamerTrialsMover, \
    MinMover
from pyrosetta import *
from util import *
import pandas as pd
import sys, os
from itertools import product
init()


# relevant mutation locations
with open( input_args.mutation_file, "rb" ) as fh:
    mutations = fh.readlines()
    mutation_set_locations = [ mut.strip().split( '_' ) for mut in mutations ]
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
mut_residues = []
mutations_made = []
delta_energies = []
poses = []

# set up how many rounds each mutation set should be made
# in order to find the lowest E mutant pose
rounds = 3

# create a Pandas DataFrame
df = pd.DataFrame()


###################
#### MUTATIONS ####
###################
# for each relevant set of mutation locations (given in PDB number)
# mutation_set_locations = [ [ 101, 102, 103 ], [ 50, 51 ] ]
for mutation_set in mutation_set_locations:
    # mutation_set = [ 101, 102, 103 ]
    # generate all mutation combinations using the Cartesian product
    # http://stackoverflow.com/questions/3099987/generating-permutations-with-repetitions-in-python
    mutation_set_mutations = [ list( p ) for p in product( AA_name1_list, repeat=len( mutation_set ) ) ]
    # for each mutation in each mutation set
    for mutations in mutation_set_mutations:
        # pick the lowest E of three clone mutants
        best_mut_pose = None
        # go one-by-one and mutate all residues in the mutation_set
        # to the appropriate mutation in mutation_set_mutations
        # three cycles so we can pick the best mutant conformation
        for round in range( 1, rounds + 1 ):
            # get a fresh copy of the input pose
            mut_pose = pose.clone()
            # for pose naming
            orig_AAs = []
            # make each individual mutation
            for ii in range( len( mutation_set ) ):
                # get the mutant residue and according mutation
                pdb_mutant_residue = int( mutation_set[ ii ] )
                mut_AA = mutations[ ii ]
                # get the original residue and the Pose number
                pose_mutant_residueA = pose.pdb_info().pdb2pose( 'A', pdb_mutant_residue )
                pose_mutant_residueB = pose.pdb_info().pdb2pose( 'B', pdb_mutant_residue )
                # for pose naming
                # this is a symmetrical dimer so we only need
                # the name of residue A (same for chain A and B)
                orig_AA = pose.residue( pose_mutant_residueA ).name1()
                orig_AAs.append( orig_AA )
                print "\nmutate_residue"
                print "round %s of %s" %( round, rounds )
                print "mutation %s of %s" %( ii + 1, 
                                             len( mutation_set ) )
                print "mutation set %s of %s" %( '_'.join( mutations ), 
                                                 '_'.join( mutation_set ) )
                # we're giving a pose number
                mut_pose.assign( mutate_residue( pose_mutant_residueA, 
                                                 mut_AA, 
                                                 mut_pose, 
                                                 sf, 
                                                 pdb_num=False ) )
                mut_pose.assign( mutate_residue( pose_mutant_residueB, 
                                                 mut_AA, 
                                                 mut_pose, 
                                                 sf, 
                                                 pdb_num=False ) )

                # create a packer task to pack only the mutation site
                print "\nstandard_packer_task"
                task = standard_packer_task( mut_pose )
                task.or_include_current( True )
                task.restrict_to_repacking()
                [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in range( 1, mut_pose.size() + 1 ) if res_num != pose_mutant_residueA and res_num != pose_mutant_residueB ]
                print "\nRotamerTrialsMover"
                rtm = RotamerTrialsMover( sf, task )
                rtm.apply( mut_pose )

                # create a MoveMap to use in a MinMover for only the mutation site
                print "\nMoveMap"
                mm = MoveMap()
                mm.set_bb_true_range( pose_mutant_residueA, pose_mutant_residueA )
                mm.set_bb_true_range( pose_mutant_residueB, pose_mutant_residueB )
                mm.set_chi_true_range( pose_mutant_residueA, pose_mutant_residueA )
                mm.set_chi_true_range( pose_mutant_residueB, pose_mutant_residueB )
                print "\nMinMover"
                minmover = MinMover( movemap_in=mm,
                                     scorefxn_in=sf,
                                     min_type_in="lbfgs_armijo_nonmonotone",
                                     tolerance_in=0.01,
                                     use_nb_list_in=True,
                                     deriv_check_in=False,
                                     deriv_check_verbose_in=False )
                minmover.apply( mut_pose )

            # by now, all mutations have been made for the set
            # rename this mutant using the mutation information
            mutation_made = '_'.join( [ orig_AAs[ii] + mutation_set[ii] + mutations[ii] for ii in range( len( mutation_set ) ) ] )
            mut_pose.pdb_info().name( mutation_made )

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
            print "round %s of %s" %( round, rounds )
            print "current mutant: %s" %sf( mut_pose )
            print "best mutant: %s" %sf( best_mut_pose )

        # store the best mutation information in the Pandas DataFrame
        mut_residues.append( [ '_'.join( mut_set ) for mut_set in mutation_set_locations ] )
        mutations_made.append( mutation_made )
        delta_energies.append( sf( best_mut_pose ) - orig_E )
        poses.append( best_mut_pose )


# append all mutation information into the Pandas Dataframe
df["mutations_made"] = mutations_made
df["dG"] = delta_energies
df["poses"] = poses
df = df.sort( "dG" )
