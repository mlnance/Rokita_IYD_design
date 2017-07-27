#!/usr/bin/python
__author__ = "morganlnance"

'''
PYRSETTA INITIALIZATION MUST CHANGE FOR EACH USE OF SCRIPT

Use PyRosetta4 to analyze a directory of structures
using an input native structure, or, if not available,
the lowest-scoring structure as the native
'''

###################
# GENERAL IMPORTS #
###################
import sys
import os
import argparse
try:
    import pandas as pd
    pandas_on = True
except ImportError:
    import csv
    pandas_on = False
# helper functions for analysis
from analyze_pdbs_util import get_sequence, get_atom_pair_distance, write_fasta_file


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description="Use PyRosetta to analyze a directory of PDBs.")
parser.add_argument("pdb_directory", type=str, help="/path/to/PDB directory")
parser.add_argument("--native_pdb_file", type=str, help="/path/to/native PDB, if available")
parser.add_argument("--dump_dir", type=str, help="/path/to/where you want to dump your data (a .csv file). \
                                                 Default is the current directory")
input_args = parser.parse_args()

# ensure input arguments are valid
# input pdb directory
if not os.path.isdir(input_args.pdb_directory):
    print("\nYou did not give me a valid pdb_directory argument.\n")
    sys.exit()
# if it is a valid directory, get its absolute path and add a trailing '/'
else:
    pdb_directory = os.path.abspath(input_args.pdb_directory) + '/'
# input dump directory, if given
if input_args.dump_dir is not None:
    if not os.path.isdir(input_args.dump_dir):
        print("\nYou did not give me a valid dump_dir argument.\n")
        sys.exit()
    # if it is a valid directory, get its absolute path and add a trailing '/'
    else:
        dump_dir = os.path.abspath(input_args.dump_dir) + '/'
# if no dump directory was given, store it as an empty variable
# this way it can be added to a path without changing anything
else:
    dump_dir = ''


##############
# INITIALIZE #
##############
# import and initialize pyrosetta
import rosetta
import pyrosetta

# IYD param files for FMN and 2IP molecules
extra_params = ' '.join(["params/2-iodophenol.params  params/flavin_mononucleotide.params"])
# initialize with flags for IYD design
# ref2015 (beta_nov15) sf was used for fast relax, so use corrections flag in init
pyrosetta.init(extra_options="-constraints:cst_fa_file FMN_2IP.cst -constraints:cst_fa_weight 10 \
                              -corrections:beta_nov15 -in:file:extra_res_fa %s" % extra_params)

# ScoreFunction should be ref2015 (beta_nov15) as defined in the init
sf = pyrosetta.get_score_function()
# test pose
# test = pyrosetta.pose_from_pdb("out_pdb/bact_IYD_00001.pdb")


###########
# ANALYZE #
###########
# holders for data
pdb_names = []
pdb_numbers = []
full_sequences = []
lid_sequences = []
scores = []
# analyzing a directory of PDBs requires a loop
# load in the PDBs one by one from the pdb_directory
for pdb in os.listdir(pdb_directory):
    # ensure we just grab .pdb files
    if pdb.endswith(".pdb") or pdb.endswith(".gz"):
        ########
        # LOAD #
        ########
        # get the absolute path to the pdb by adding in the pdb_directory
        pdb = pdb_directory + pdb
        # add this pdb name to the pdb_names holder
        pdb_name = pdb.split('/')[-1].split(".gz")[0].split(".pdb")[0]
        # get the pdb number
        # a typical pdb for the IYD fast relax looks like
        # bact_IYD_*.pdb (ex. bact_IYD_00117.pdb)
        try:
            pdb_num = int(pdb.split('/')[-1].split('.')[0].split('_')[-1])
        except:
            # if something doesn't work with the above, just ignore and set number as its name
            pdb_num = pdb_name
        # load in the pdb as a pose
        try:
            pose = pyrosetta.pose_from_pdb(pdb)
        # if the pose can't be loaded in, move on to the next one
        except:
            continue

        ############
        # SEQUENCE #
        ############
        # get the sequence of the pose
        try:
            full_seq = get_sequence(pose)
        except:
            continue
        # get the sequence of the lid region that was mutated
        # pose numbers of lid (within 5A of catalytic site)
        lid_seq = [86, 90, 94, 98, 102, 107, 108, 111, 167]
        try:
            lid_seq = get_sequence(pose, lid_seq)
        except:
            continue
        # compare to native sequence?

        ####################
        # HYDROGEN BONDING #
        ####################
        '''
        # calculate the distance between atoms of particular interest for IYD
        # catalytic Thr (pose resi 168) HG1 (atom num 11) to FMN (pose resi 219) N5 (atom num 11)
        thr_HG1_fmn_N5 = get_atom_pair_distance(pose, 168, 11, 219, 11)
        # FMN (pose resi 219) HO4' (atom num 50) to 2IP (pose resi 220) to O1 (atom num 8)
        fmn_HO4_2ip_O1 = get_atom_pair_distance(pose, 219, 50, 220, 8)
        # ALA 64 (pose resi 59) H (atom num 6) to 2IP (pose resi 220) to O1 (atom num 8)
        ala64_H_2ip_O1 = get_atom_pair_distance(pose, 59, 6, 220, 8)
        print thr_HG1_fmn_N5, fmn_HO4_2ip_O1, ala64_H_2ip_O1
        '''

        ###########
        # SCORING #
        ###########
        # re-score for now because I can
        # the sf and the pose need to know the constraints you put in
        try:
            rosetta.core.scoring.constraints.add_fa_constraints_from_cmdline(pose, sf)
            score = sf(pose)
        except:
            continue


        ################
        # COLLECT INFO #
        ################
        # a pose was able to be loaded in and all the information was acquired
        # it is likely safe to say that all data can be stored now
        pdb_names.append(pdb_name)
        pdb_numbers.append(pdb_num)
        full_sequences.append(full_seq)
        lid_sequences.append(lid_seq)
        scores.append(score)


##############
# WRITE DATA #
##############
# specifically for IYD, write out the lid sequences
with open(dump_dir + "bact_IYD_lid_sequences.txt","w") as fh:
    for line in lid_sequences:
        fh.write( line + "\n" )

# write out a fasta file of the pdb names and sequences
write_fasta_file(pdb_names=pdb_names,
                 pdb_sequences=full_sequences,
                 filename=dump_dir+"bact_IYD_fasta.txt")

# use pandas if it is available
if pandas_on:
    # store collected data in a pandas dataframe
    df = pd.DataFrame()
    df["pdb_name"] = pdb_names
    df["pdb_num"] = pdb_numbers
    df["lid_seq"] = lid_sequences
    df["score"] = scores
    df["full_seq"] = full_sequences

    print(df)
    df.to_csv(dump_dir + "analyze_pdbs.csv")
# otherwise, use python's csv module
else:
    # zip the data together to write row-by-row
    data = zip(pdb_names, pdb_numbers, full_sequences, lid_sequences, scores)
    with open(dump_dir + "analyze_pdbs.csv", 'w') as fh:
        writer = csv.writer(fh)
        header = ["pdb_name", "pdb_num", "lid_seq", "score", "full_seq"]
        print header
        writer.writerow(header)
        for line in data:
            print line
            writer.writerow(line)
