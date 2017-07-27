#!/bin/bash

: '
THIS ALSO MAKES THE CORRESPONDING DIRECTORY IN out_pdb/multiple_ligand_docking

Provide the name of the PDB file and generate an appropriate ligand docking flags file
Assuming use of a specific PDB for ligand docking
ligands/2IP.pdb ligands/FMN.pdb

Assuming specific paths and file name trailings, so check for that in
(1) in file flag
(2) protocol name
(3) score file name
(4) out pdb file name

Example Usage:
./gen_lig_dock_mld_flags.sh bact_IYD_02961
where $1 = bact_IYD_02961 gives:
(1) -in:file:s in_pdb/multiple_ligand_docking/$1_docking.pdb ligands/2IP.pdb
(2) -parser:protocol 1$1_dock.xml
(3) -out:file:scorefile $1_multiple_ligand_docking_score.sc
(4) -out:path:pdb out_pdb/multiple_ligand_docking/$1

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
mkdir -p "out_pdb/multiple_ligand_docking/$1"

# for the flag file
echo "# Pound signs indicate comments 

# -in:file:s option imports the protein and ligand PDB structures
# -in:file:extra_res_fa option imports the parameters for the ligand
-in:file:s 'in_pdb/multiple_ligand_docking/$1_docking.pdb ligands/2IP.pdb ligands/FMN.pdb'
-extra_res_fa params/2-iodophenol.params params/flavin_mononucleotide.params


# the packing options allow Rosetta to sample additional rotamers for
# protein sidechain angles chi 1 (ex1) and chi 2 (ex2) 
# no_optH false tells Rosetta to optimize hydrogen placements
# flip_HNQ tells Rosetta to consider HIS,ASN,GLN hydrogen flips
# ignore_ligand_chi prevents Roseta from adding additional ligand rotamer
-packing
	-ex1
	-ex2
	-no_optH false
	-flip_HNQ true
# warning is given that this flag is not used (MLN)
	-ignore_ligand_chi true


# number of decoy structures
-nstruct 10000


# parser:protocol locates the XML file for RosettaScripts
-parser
	-protocol 1$1_mld_dock.xml


# overwrite allows Rosetta to write over previous structures and scores
# not using overwrite as using the multiple processes flag (MLN)
#-overwrite


# Ligand docking is not yet benchmarked with the updated scoring function
# This flag restores certain parameters to previously published values
-mistakes
	-restore_pre_talaris_2013_behavior true 


# out files and score file
-out:file:scorefile     $1_multiple_ligand_docking_score.sc
-out:path:pdb   out_pdb/multiple_ligand_docking/$1
-out:pdb_gz

# for visualizing in pymol
#-show_simulation_in_pymol 0.2
#-keep_pymol_simulation_history


# for running on Jazz
-multiple_processes_writing_to_one_directory True"
