#!/usr/bin/python
__author__ = "morganlnance"

'''
Use Python to read PDB files and return the bact_IYD lid sequences in a .csv file
'''

#################
# BACT_IYD INFO #
#################
# SPECIFICALLY FOR BACT_IYD PROJECT
# PDB RESIDUE NUMBERS FOR THE BACT_IYD LID SEQUENCE
# COLLECT ONLY FOR CHAIN A SINCE IT IS A SYMMETRIC PROTEIN
bact_IYD_lid = [91, 95, 99, 103, 107, 112, 113, 116, 172]
bact_IYD_chain = 'A'
# dictionary for connecting residue name to residue letter
AA_dict = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
           "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
           "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"}


###########
# IMPORTS #
###########
import sys
import os
try:
    import pandas as pd
    pandas_on = True
except ImportError:
    import csv
    pandas_on = False
import argparse


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description="Use Python to read a list of PDB files and \
                                             return a .csv file of the bact_IYD lid sequences. \
                                             This program can handle gzipped PDBs.")
parser.add_argument("pdb_dir", type=str, help="/path/to/pdb directory")
parser.add_argument("--dump_dir", type=str, help="/path/to/where you want to dump your data (a .csv file). \
                                                 Default is the current directory")
input_args = parser.parse_args()


###################
# CHECK ARGUMENTS #
###################
# check the pdb_dir
if not os.path.isdir(input_args.pdb_dir):
    print "\nYou did not give me a valid pdb_dir argument.\n"
    sys.exit()
# if it is a valid directory, get its absolute path and add a trailing '/'
else:
    pdb_dir = os.path.abspath(input_args.pdb_dir) + '/'
# check the dump_dir, if given
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


############
# GET PDBs #
############
# get a list of pdb files from the pdb_dir
# use the pdb_dir path to get the absolute path to all pdbs
pdb_files = [pdb_dir + pdb for pdb in os.listdir(pdb_dir) if pdb.endswith(".pdb") or pdb.endswith(".gz")]


#############
# READ PDBs #
#############
# for each pdb in the directory, unzip (if needed) and read it
# looking for the bact_IYD lid sequence information set at top of file
pdb_names = []
lid_sequences = []
#full_sequences = []
for pdb in pdb_files:
    # grab the pdb name, which should end with .pdb, if not also .gz
    pdb_name = pdb.split('/')[-1].split(".gz")[0].split(".pdb")[0]

    #########
    # UNZIP #
    #########
    # check if pdb needs to be unzipped
    rezip_pdb = False
    if pdb.endswith(".gz"):
        # unzip using terminal
        try:
            os.popen("gunzip %s" % pdb)
        except:
            continue
        # now pdb name is the same minus the ending ".gz", so update pdb name for rezipping
        pdb = pdb.split(".gz")[0]
        # update the bool saying that this pdb needs to be rezipped
        rezip_pdb = True

    ###########
    # READ IN #
    ###########
    # try to open the pdb, skip if it doesn't work for any reason
    try:
        with open(pdb, 'r') as fh:
            pdb_lines = fh.readlines()
    except:
        continue

    ############
    # READ PDB #
    ############
    # make local copy of bact_IYD_lid
    # as each residue is noted, remove it from list and
    # start looking for the next residue
    local_bact_IYD_lid = [ii for ii in bact_IYD_lid]
    # lid sequence will be stored as a string
    lid_sequence = ''
    # for each line, check for chain A and residue number in local_bact_IYD_lid
    for line in pdb_lines:
        # if all lid residues have been checked, break and continue to next pdb
        # lid residues are removed one-by-one as they are seen to save time
        if len(local_bact_IYD_lid) == 0:
            break
        # otherwise, look for bact_IYD lid residues and get its sequence
        else:
            # check ATOM lines
            if line.startswith("ATOM"):
                ################
                # LID SEQUENCE #
                ################
                # check residue number, which is always in the same spot in a PDB
                # remove whitespace and make it an integer
                res_num = int(line[22:26].replace(' ', ''))
                # check only chain A for the lid sequence, for clarity
                chain = line[21:22]
                # if this residue is a lid residue and on the chain of interest (A)
                if res_num in local_bact_IYD_lid and chain == bact_IYD_chain:
                    # check the amino acid sequence of this residue
                    # should be three letters, but replace whitespace anyway
                    res_name = line[17:20].replace(' ', '')
                    # use dictionary to get one-letter amino acid code
                    try:
                        res_letter = AA_dict[res_name]
                    # if key error, store res_letter as X
                    except KeyError:
                        res_letter = 'X'
                    # add the res_letter to the growing lid_sequence string
                    lid_sequence += res_letter
                    # this lid residue has contributed its sequence
                    # remove from list and continue
                    local_bact_IYD_lid.remove(res_num)

    #########
    # REZIP #
    #########
    # rezip the pdb, if needed
    if rezip_pdb:
        try:
            os.popen("gzip %s" % pdb)
        except:
            continue

    ################
    # COLLECT DATA #
    ################
    pdb_names.append(pdb_name)
    lid_sequences.append(lid_sequence)


##############
# WRITE DATA #
##############
# specifically for IYD, write out the lid sequences
with open(dump_dir + "bact_IYD_lid_sequences.txt","w") as fh:
    for line in lid_sequences:
        fh.write( line + "\n" )

# use pandas if it is available
if pandas_on:
    # store collected data in a pandas dataframe
    df = pd.DataFrame()
    df["pdb_name"] = pdb_names
    df["lid_seq"] = lid_sequences
    print(df)
    df.to_csv(dump_dir + "bact_IYD_lid_sequences.csv")
# otherwise, use python's csv module
else:
    # zip the data together to write row-by-row
    data = zip(pdb_names, lid_sequences)
    with open(dump_dir + "bact_IYD_lid_sequences.csv", 'w') as fh:
        writer = csv.writer(fh)
        header = ["pdb_name", "lid_seq"]
        print header
        writer.writerow(header)
        for line in data:
            print line
            writer.writerow(line)
