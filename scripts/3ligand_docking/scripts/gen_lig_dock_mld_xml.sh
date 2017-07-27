#!/bin/bash

: '
Provide the name of the PDB file and generate an appropriate ligand docking xml file

Assuming specific paths and file name trailings, so check for that in

Example Usage:
./gen_lig_dock_mld_xml.sh bact_IYD_02961
where $1 = bact_IYD_02961 gives:
(1) native=native_pdb/multiple_ligand_docking/$1.pdb (in InterfaceScoreCalculator)
'

# check that an argument was given
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Give me the name of a PDB with which to generate an XML file"
    exit 1
fi

# for the xml file
echo "<ROSETTASCRIPTS>
		Reweight lines taken as per Rocco's suggestion from 
		https://github.com/RosettaCommons/main/blob/53f5472af0a6c12a26bcebbc3cf1f2f4864006e0/source/src/apps/benchmark/performance/ligand_dock/ligand_dock_script.xml
		<SCOREFXNS>
			<ScoreFunction name=\"ligand_soft_rep\" weights=\"ligand_soft_rep\">
				<Reweight scoretype=\"fa_elec\" weight=\"0.42\"/>
				<Reweight scoretype=\"hbond_bb_sc\" weight=\"1.3\"/>
				<Reweight scoretype=\"hbond_sc\" weight=\"1.3\"/>
				<Reweight scoretype=\"rama\" weight=\"0.2\"/>
			</ScoreFunction>
			<ScoreFunction name=\"hard_rep\" weights=\"ligand\">
				<Reweight scoretype=\"fa_intra_rep\" weight=\"0.004\"/>
				<Reweight scoretype=\"fa_elec\" weight=\"0.42\"/>
				<Reweight scoretype=\"hbond_bb_sc\" weight=\"1.3\"/>
				<Reweight scoretype=\"hbond_sc\" weight=\"1.3\"/>
				<Reweight scoretype=\"rama\" weight=\"0.2\"/>
			</ScoreFunction>
		</SCOREFXNS>

		<LIGAND_AREAS>
		  	LigandArea describes parameters specific for each ligand
			Cutoff is the distance in A from the ligand an amino acid's C-beta atom 
			can be and that residue still be part of the interface
			<LigandArea name=\"inhibitor_dock_sc\" chain=\"X\" cutoff=\"6.0\" add_nbr_radius=\"true\" all_atom_mode=\"false\"/>
			<LigandArea name=\"FMN_dock_sc\" chain=\"Z\" cutoff=\"6.0\" add_nbr_radius=\"true\" all_atom_mode=\"false\"/>
			<LigandArea name=\"inhibitor_final_sc\" chain=\"X\" cutoff=\"6.0\" add_nbr_radius=\"true\" all_atom_mode=\"false\"/>
			<LigandArea name=\"FMN_final_sc\" chain=\"Z\" cutoff=\"6.0\" add_nbr_radius=\"true\" all_atom_mode=\"false\"/>
			<LigandArea name=\"inhibitor_final_bb\" chain=\"X\" cutoff=\"7.0\" add_nbr_radius=\"false\" all_atom_mode=\"true\" Calpha_restraints=\"0.3\"/>
			<LigandArea name=\"FMN_final_bb\" chain=\"Z\" cutoff=\"7.0\" add_nbr_radius=\"false\" all_atom_mode=\"true\" Calpha_restraints=\"0.3\"/>
		</LIGAND_AREAS>

		<INTERFACE_BUILDERS>
		  	InterfaceBuilder describes how to choose residues that will be 
			part of a protein-ligand interface. Residiues chosen are for 
			repacking, rotamer trails, and bb min during ligand docking
			InterfaceBuilder builds off of LigandArea, where there is
			a LigandArea for each ligand specified to be docked
			<InterfaceBuilder name=\"side_chain_for_docking\" ligand_areas=\"inhibitor_dock_sc,FMN_dock_sc\"/>
			<InterfaceBuilder name=\"side_chain_for_final\" ligand_areas=\"inhibitor_final_sc,FMN_final_sc\"/>
			<InterfaceBuilder name=\"backbone\" ligand_areas=\"inhibitor_final_bb,FMN_final_bb\" extension_window=\"3\"/>
		</INTERFACE_BUILDERS>

		<MOVEMAP_BUILDERS>
		  	MoveMapBuilder builds off of InterfaceBuilder which was
			built off of LigandArea. So each MoveMapBuilder is built
			for each ligand specified in InterfaceBuilder that had
			a LigandArea also made for it. i.e. one MoveMapBuilder
			for all ligands for each type of movemap 
			<MoveMapBuilder name=\"docking\" sc_interface=\"side_chain_for_docking\" minimize_water=\"false\"/>
			<MoveMapBuilder name=\"final\" sc_interface=\"side_chain_for_final\" bb_interface=\"backbone\" minimize_water=\"false\"/>
		</MOVEMAP_BUILDERS>

		<SCORINGGRIDS ligand_chain=\"X\" width=\"15\">
			<ClassicGrid grid_name=\"classic\" weight=\"1.0\"/>
		</SCORINGGRIDS>

		<MOVERS>
			added initial_perturb flag as per Rocco's suggestion to get 
			additional movement from initial starting position as to 
			not bias the simulation too much with the known answer
			Rocco suggest initial_perturb to be between 3-5 Angstroms
			Rocco also said typically 5 Ang is chosen for box_size, so going with that
			Chain X is 2IP, Rocco suggested that FMN (chain Z) won't need the Transform
			because it likely won't rearrange much in the pocket
			<Transform name=\"transform\" chain=\"X\" box_size=\"5.0\" move_distance=\"0.2\" angle=\"20\" cycles=\"500\" repeats=\"1\" temperature=\"5\" initial_perturb=\"4\"/>
			increased the number of cycles to increase sampling 
			but with also getting a decoy about every minute or two 
			(cycles=100 more or less arbitrary choice though)
			also increase number of times a repack occurs (can be adjusted more)
			Additionally, Rocco said that there was not much benefit in 
			increasing the number of cycles in the HighResRocker 
			(aggressive sampling didn't improve results that much)
			<HighResDocker name=\"high_res_docker\" cycles=\"20\" repack_every_Nth=\"1\" scorefxn=\"ligand_soft_rep\" movemap_builder=\"docking\"/>
			<FinalMinimizer name=\"final\" scorefxn=\"hard_rep\" movemap_builder=\"final\"/>
			<InterfaceScoreCalculator name=\"add_scores\" chains=\"X,Z\" scorefxn=\"hard_rep\" native=\"native_pdb/multiple_ligand_docking/bact_IYD_02961.pdb\"/> 
		</MOVERS>

		<PROTOCOLS>
			<Add mover_name=\"transform\"/>
			<Add mover_name=\"high_res_docker\"/>
			<Add mover_name=\"final\"/>
			<Add mover_name=\"add_scores\"/>
		</PROTOCOLS>

</ROSETTASCRIPTS>"
