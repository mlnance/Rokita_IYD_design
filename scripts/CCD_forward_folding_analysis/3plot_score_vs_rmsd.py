#!/usr/bin/python

'''
Using the .tsv file (or .csv) that the score_vs_rmsd_full.py script outputs
Use pandas and matplot lib to plot score vs rmsd
Manually change columns to get different plots
CA_RMSD
BB_RMSD
ALL_ATOM_RMSD
'''

import sys
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )
import pandas as pd


# import the .tsv file
try:
    data_file = sys.argv[1]
except IndexError:
    print "\nGive me a .tsv file from score_vs_rmsd.py output\n"
    sys.exit()

# load the .tsv file
# or it could be a .csv file
try:
    if data_file.endswith('.tsv'):
        data_df = pd.DataFrame().from_csv(data_file, sep='\t', index_col=None)
    elif data_file.endswith('.csv'):
        data_df = pd.DataFrame().from_csv(data_file)
        #data_df = pd.DataFrame().from_csv(data_file, index_col=None)
    else:
        print "\nWhat kind of data file did you give me?\n"
        sys.exit()
except:
    print "\nDid you give me an appropriate .tsv file?"
    print "I could not use Pandas to read the file.\n"
    sys.exit()

# get rmsd and score data
x = data_df.CA_RMSD.tolist()
#x = data_df.BB_RMSD.tolist()
#x = data_df.ALL_ATOM_RMSD.tolist()
y = data_df.total.tolist()

# native data
# bact_IYD_0291
native_score = -1161.86

# plot
plt.scatter(x, y, zorder=3)
plt.grid(zorder=0)
# native point
plt.scatter(0, native_score, c='red', 
            clip_on=False, s=22, zorder=100)
plt.xlabel("CA_RMSD")
#plt.xlabel("BB_RMSD")
#plt.xlabel("ALL_ATOM_RMSD")
plt.ylabel("REU")
plt.xlim([0, 8])
plot_title = "CCD forward folding test bact_IYD_02961"
plt.title(plot_title, y=1.02)
plt.show(block=False)
