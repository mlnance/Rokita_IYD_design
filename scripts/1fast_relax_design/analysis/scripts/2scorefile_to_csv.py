#!/usr/bin/python
__author__ = "morganlnance"

'''
Use Python to read a score.sc file and return a .csv file of the information
Allows for overwriting multiple entries of the same decoy in the score.sc file
where only the most recent (last entry) decoy in the file is saved
i.e. handles when a score.sc file has multiple scores for the same decoy
'''


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
parser = argparse.ArgumentParser(description="Use Python to read a score.sc file and return a .csv file of the info.")
parser.add_argument("score_file", type=str, help="/path/to/score.sc file")
parser.add_argument("--dump_dir", type=str, help="/path/to/where you want to dump your data (a .csv file). \
                                                 Default is the current directory")
input_args = parser.parse_args()


#############
# READ FILE #
#############
# try to open and read the file
score_file = input_args.score_file
try:
    with open(score_file, "r") as fh:
        lines = fh.readlines()
except:
    print "\nSomething was wrong with your score file: %s\n" % score_file
    sys.exit()
# get the score file's name
score_file_name = score_file.split('/')[-1].split(".sc")[0]


##################
# CHECK DUMP DIR #
##################
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


##################
# INTERPRET FILE #
##################
# get the header information (the entry information for each score value in the file)
# this should always be line two (line 1 is "SEQUENCE: ")
header = lines[1].split()
# skip the first entry in the header because that's just "SCORE: "
header = header[1:]
# now pull out the data, also skipping the first "SCORE: " entry in each line
data = [line.split()[1:] for line in lines[2:]]

# since there will be repeating data, it will be best to create a dictionary
# that way, entries that repeated will be overwritten and the last pdb info will be saved
data_dict = {}
# use the "description" header entry to get the pdb names
# a pdb name will be used as the dictionary's key
description_index = header.index("description")
# the dictionary's value will be set as the corresponding data line
for line in data:
    # use the pdb name as the key and the rest of the line as the value
    # this should allow overwrite to where the last version of the pdb is the one that's kept
    data_dict[line[description_index]] = line

# now consolidate the data again from the dictionary back into the data object
# all the values of the dictionary are the score lines corresponding to the last pdb of each entry
# this was done to avoid having multiple entries for the same pdb
data = data_dict.values()
# convert all relevant entries to floats (description column will not work)
# for each line of data
for ii in range(len(data)):
    # for each entry in that line
    line = data[ii]
    for jj in range(len(line)):
        # try to turn this entry in this line to a float
        try:
            data[ii][jj] = float(data[ii][jj])
        # if it doesn't work, it's a string, so keep it what it is
        except ValueError:
            data[ii][jj] = data[ii][jj]
# zip the data together so that all types of entries are grouped together accordingly
data = zip(*data)

# use the description index to pull out the pdb names
#description_index = header.index("description")
#pdb_names = data[description_index]
# split on '_' and use the last entry of that list as the pdb number
# ex) bact_IYD_00117
#pdb_nums = [int(name.split('_')[-1]) for name in pdb_names]
# one-line version to avoid having to make a copy of "description" data
try:
    pdb_nums = [int(name.split('_')[-1]) for name in data[header.index("description")]]
except:
    # if something didn't work, oh well
    pdb_nums = None


##############
# WRITE DATA #
##############
# if pandas was loaded in, use that module
if pandas_on:
    # connect the header with the data in a pandas DataFrame
    df = pd.DataFrame()
    for ii in range(len(header)):
        # pull out the corresponding header entry
        entry = header[ii]
        # connect this header with its data
        df[entry] = data[ii]
    # if pdb numbers were acquired, add that in
    if pdb_nums is not None:
        df["pdb_num"] = pdb_nums
        df = df.sort_values("pdb_num")
    # output the dataframe
    print df
    df.to_csv(dump_dir + score_file_name + ".csv")
# otherwise, pandas is not accessible, so use python's csv module
# this won't write the index column, but oh well
else:
    # add the pdb_num header and data, if collected
    if pdb_nums is not None:
        header.append("pdb_num")
        data.append(pdb_nums)
    # write the csv
    with open(dump_dir + score_file_name + ".csv", 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for line in zip(*data):
            writer.writerow(line)
