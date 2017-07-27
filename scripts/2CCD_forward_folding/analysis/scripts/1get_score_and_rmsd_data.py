#!/usr/bin/python

'''
THIS SCRIPT USES PYMOL FOR RMSD CALCULATION

Example: pymol -qc script.py -- /path/native.pdb /path/data_filename.csv /path/decoy_filenames_list.txt
    if you want pymol to do its job but not actually open a window
    otherwise, remove -qc
    the -- tells pymol to stop interpreting input and take rest as args
    can't do this interactively, hence asking for a filename for dumping the dataframe

Python functions to be used by PyMOL to get the rmsd of certain selections 
of decoys vs a given native structure
PyMOL commands at the bottom of all the function declarations
Requires that you specify:
    1) native pdb file path
    2) decoys (specified at the end like /path/to/decoys/*pdb)
    3) filename for the .csv file to be made with the data
        a) using Pandas dataframe, so best to dump it each run
    4) which residues and atom names you want to get the rmsd of
        a) split up into CA, BB (CA+C+N), and ALL_ATOM
'''

###########
# IMPORTS #
###########
import pandas as pd
import argparse
import sys
import os
import gzip


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description='Use PyMOL to calculate RMSDs.')
parser.add_argument('native', type=str, help='/path/to/native PDB.')
parser.add_argument('data_filename', type=str, help='/path/to/filename where you want your .csv file dumped.')
parser.add_argument('decoy_list', type=str, help='A list of decoy files to analyze.')
input_args = parser.parse_args()

# read and interpret the list of decoys to analyze
with open(input_args.decoy_list, 'r') as fh:
    decoy_filenames = fh.readlines()
# newline character throws stuff off
decoy_filenames = [f.strip() for f in decoy_filenames]
# check that are filenames are valid
for f in decoy_filenames:
    if not os.path.isfile(f):
        print '\n%s is not a valid file.\n' %f
        sys.exit()


#############
# FUNCTIONS #
#############
def get_rmsd(native, decoy, rmsd_type):
    '''
    Get the rmsd (of a hardcoded selection) of <native> to <decoy> using <rmsd_type>
    Assumes <native> and <decoy> have already been loaded using cmd.load()
    Also assumes that <native> and <decoy> have the same sequence and residue/atom order
    <rmsd_types> can be 1) 'CA', 2) 'BB', or 3) 'ALL_ATOM'
    :param native: str(native pdb file path)
    :param decoy: str(decoy pdb file path)
    :param rmsd_type: str(CA, BB, or ALL_ATOM)
    :return: float(rmsd)
    '''
    ## define the residues and atoms of interest
    # will be prefaced by the name of the native and decoy
    # CA atom rmsd
    CA_rmsd_residues = ' and resi 89-113 and chain A and name CA'
    # BB backbone atom rmsd
    BB_rmsd_residues = ' and resi 89-113 and chain A and name CA+C+N'
    # ALL_ATOM rmsd
    ALL_ATOM_rmsd_residues = ' and resi 89-113 and chain A'

    # the first result from pymol's align command is the rmsd
    # this is now a pymol command, so use cmd.
    # cycles=0 means don't throw out bad atoms
    # transform=0 means don't actually move the structure during alignment
    # CA rmsd
    if rmsd_type == 'CA':
        rmsd = cmd.align(decoy + CA_rmsd_residues, native + CA_rmsd_residues, cycles=0, transform=0)[0]
    elif rmsd_type == 'BB':
        rmsd = cmd.align(decoy + BB_rmsd_residues, native + BB_rmsd_residues, cycles=0, transform=0)[0]
    elif rmsd_type == 'ALL_ATOM':
        rmsd = cmd.align(decoy + ALL_ATOM_rmsd_residues, native + ALL_ATOM_rmsd_residues, cycles=0, transform=0)[0]
    else:
        print "\nYou didn't give me a proper <rmsd_type>. Use 'CA', 'BB', or 'ALL_ATOM'.\n"
        return None

    return rmsd

def get_total_score_from_pdb_file(pdb_filepath):
    '''
    Open the <pdb_filepath> and find the 'pose' line of the energy table.
    The last entry for that line should be the total score
    :param pdb_filepath: str(/path/to/pdb to get energy)
    :return: float(total score)
    '''
    # ensure the filepath exists
    if not os.path.isfile(pdb_filepath):
        return None

    # open the file according to whether it is compressed or not
    if pdb_filepath.endswith('.gz'):
        # open and get the lines of the pdb
        try:
            with gzip.open(pdb_filepath, 'r') as fh:
                pdb_lines = [line for line in fh]
        # but sometimes the cake is a lie apparently and a .gz file isn't actually compressed...
        except IOError:
            with open(pdb_filepath, 'r') as fh:
                pdb_lines = fh.readlines()
    else:
        # otherwise it's not compressed, so open normally
        with open(pdb_filepath, 'r') as fh:
            pdb_lines = fh.readlines()

    # find the 'pose' line
    # or return None if that isn't present
    for line in pdb_lines:
        if line.startswith('pose'):
            # return None if something is wrong with this line
            try:
                return float(line.split()[-1])
            except ValueError:
                return None
    
    # if nothing was found, return None
    return None
    


##################
# PYMOL COMMANDS #
##################
# load the native and get its name
native_name = input_args.native.split('/')[-1].split('.')[0]
cmd.load(input_args.native)

## for data collection
pdb_names = []
CA_rmsds = []
BB_rmsds = []
ALL_ATOM_rmsds = []
scores = []

# for each decoy given
# load it
# get the rmsd to native
# add rmsd to dataframe
for decoy in decoy_filenames:
    # get the decoy name by splitting on '/' and '.'
    decoy_name = decoy.split('/')[-1].split('.')[0]
    print decoy_name

    # load the decoy
    cmd.load(decoy)

    # get the rmsd of the decoy to native
    CA_rmsd = get_rmsd(native_name, decoy_name, 'CA')
    BB_rmsd = get_rmsd(native_name, decoy_name, 'BB')
    ALL_ATOM_rmsd = get_rmsd(native_name, decoy_name, 'ALL_ATOM')

    # remove the decoy as to not overload poor pymol
    cmd.remove(decoy_name)
    cmd.delete(decoy_name)

    # get the score of this structure from its pdb file
    # function returns None if it can't find a score
    score = get_total_score_from_pdb_file(decoy)

    ## collect data in dataframe
    pdb_names.append(decoy_name)
    CA_rmsds.append(CA_rmsd)
    BB_rmsds.append(BB_rmsd)
    ALL_ATOM_rmsds.append(ALL_ATOM_rmsd)
    scores.append(score)

########
# DATA #
########
# create the dataframe
df = pd.DataFrame()
df['name'] = pdb_names
df['CA_RMSD'] = CA_rmsds
df['BB_RMSD'] = BB_rmsds
df['ALL_ATOM_RMSD'] = ALL_ATOM_rmsds
df['score'] = scores
# dump
if not input_args.data_filename.endswith('.csv'):
    df.to_csv(input_args.data_filename + '.csv')
else:
    df.to_csv(input_args.data_filename)
