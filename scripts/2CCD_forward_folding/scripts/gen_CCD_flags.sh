#!/bin/bash

: '
THIS ALSO MAKES THE CORRESPONDING DIRECTORY IN out_pdb

Provide the name of the PDB file and generate an appropriate CCD lid forward folding flags file

Assuming specific paths and file name trailings, so check for that in
(1) in file flag
(2) Robetta fragment file names (.frag3 and .frag9)
(3) score file name
(4) out pdb file name

Example Usage:
./gen_CCD_flags.sh bact_IYD_02961
where $1 = bact_IYD_02961 gives:
(1) -in:file:s in_pdb/$1.pdb
(2) -loops:frag_files Robetta_fragments/$1/$1.frag9 Robetta_fragments/$1/$1.frag3
(3) -out:file:scorefile CCD_$1_score.sc
(4) -out:path:pdb out_pdb/$1

'
# check that an argument was given
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Give me the name of a PDB with which to generate a flags file"
    echo "Note that this script will make a corresponding directory in out_pdb/"
    exit 1
fi

# for the out_pdb directory, if not already there
mkdir -p "out_pdb/$1"

# for the flag file
echo "# Executable: /Rosetta/main/source/bin/loopmodel.<system>

-in:file:s in_pdb/$1.pdb
# params for 2IP and FMN are (or should	be) in the Rosetta database
# update both fa_standard and centroid params
# be sure to add I atom to centroid atom properties
# and add the param paths to the respective residue_types.txt file

# loop modeling flag
-loops:remodel quick_ccd
# unknown
-loops:refine refine_ccd
-loops:loop_file CCD_bact_IYD_lid.loop
# force extended on loops, independent of loop input file
-loops:extended true

-loops:frag_sizes 9 3
-loops:frag_files Robetta_fragments/$1/$1.frag9 Robetta_fragments/$1/$1.frag3

-nstruct 5000

-ex1
-ex2

-use_input_sc
-flip_HNQ
-no_optH false

-out:file:scorefile	CCD_$1_score.sc
-out:path:pdb		out_pdb/$1
-out:pdb_gz

#-mpi_tracer_to_file	outerr/relax_tracer.out
# Jazz
-multiple_processes_writing_to_one_directory True

#-show_simulation_in_pymol 3
#-keep_pymol_simulation_history"
