#!/usr/bin/python

'''
Using the .tsv file (or .csv) that the score_vs_rmsd_full.py script outputs
Use pandas and matplot lib to plot score vs rmsd
Manually change columns to get different plots
CA_RMSD
BB_RMSD
ALL_ATOM_RMSD
'''

###########
# IMPORTS #
###########
import sys
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )
import pandas as pd
import argparse


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description='Use matplotlib to plot RMSD vs score from data from a .csv file.')
parser.add_argument('data_file', type=str, help='/path/to/data.csv file.')
input_args = parser.parse_args()


###########
# LOADING #
###########
# load the .csv file
if input_args.data_file.endswith('.csv'):
    try:
        data_df = pd.DataFrame().from_csv(input_args.data_file)
    except:
        print "\nDid you give me an appropriate .csv file?"
        print "I could not use Pandas to read the file.\n"
        sys.exit()
else:
    print "\nWhat kind of data file did you give me? It wasn't .csv\n"
    sys.exit()


###################
# WARNING MESSAGE #
###################
print '\nDid you look in the script and update the values needed?'


##############################
# GET RMSD AND SCORE FROM DF #
##############################
x = data_df.CA_RMSD.tolist()
xlabel = 'CA_RMSD'
#x = data_df.BB_RMSD.tolist()
#xlabel = 'BB_RMSD'
#x = data_df.ALL_ATOM_RMSD.tolist()
#xlabel = 'ALL_ATOM_RMSD'
y = data_df.score.tolist()


####################
# CHANGE PER DECOY #
####################
### bact_IYD_native
#native_score = -1056.110
#plot_title = "CCD forward folding test bact_IYD_native"
#ylim = [-1080, -1000]
### bact_IYD_02961
native_score = -1161.86
lot_title = "CCD forward folding test bact_IYD_02961"
ylim = [-1170, -1100]
### bact_IYD_01372
#native_score = -1177.22
#ylim = [-1180, -1080]
#plot_title = "CCD forward folding test bact_IYD_01372"


# plot
plt.scatter(x, y, zorder=3)
plt.grid(zorder=0)
# native point
plt.scatter(0, native_score, c='red', 
            clip_on=False, s=22, zorder=100)
plt.xlabel(xlabel)
plt.ylabel('REU')
plt.xlim([0, 10])
plt.ylim(ylim)
plt.title(plot_title, y=1.02)
plt.show(block=False)
# save plot
#plt.savefig('CCD_bact_IYD_02961_rmsd_vs_score.png', transparent=True, dpi=1080)
#plt.savefig('CCD_bact_IYD_01372_rmsd_vs_score.png', transparent=True, dpi=1080)
