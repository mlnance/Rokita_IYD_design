#!/usr/bin/python

'''
Example execution:
for f in in_pdb/*pdb
do
  ./0prepare_score_vs_rmsd_full.py $f
done

non-gzipped copies of all PDBs to be analyzed
should be available since this script
will overwrite the input PDB
for bact_IYD, need to replace
flavin_mononucleotide with FMN and
2-iodophenol with 2IP
in the decoy scripts
'''

import sys
import os

pdb_name = sys.argv[1]

short_pdb_name = pdb_name.split( ".pdb" )[0]

sed_cmd = "sed 's/flavin_mononucleotide/FMN/g' %s > %s.temp" %( pdb_name, short_pdb_name )
os.popen( sed_cmd )
sed_cmd = "sed 's/2-iodophenol/2IP/g' %s.temp > %s.temp2" %( short_pdb_name, short_pdb_name )
os.popen( sed_cmd )

# remove the .temp file
rm_cmd = "rm %s.temp" % short_pdb_name
os.popen( rm_cmd )

# overwrite the input .pdb file with the .temp2 file
rename_cmd = "mv %s.temp2 %s" %( short_pdb_name, pdb_name )
os.popen( rename_cmd )
