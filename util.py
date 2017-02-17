2#!/usr/bin/python
__author__ = 'morganlnance'


# global variables
AA_name1_list = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
AA_name3_list = [ "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR" ]
AA_name1_to_name3 = { 'A':"ALA", 'C':"CYS", 'D':"ASP", 'E':"GLU", 'F':"PHE", 'G':"GLY", 'H':"HIS", 'I':"ILE", 'K':"LYS", 'L':"LEU", 'M':"MET", 'N':"ASN", 'P':"PRO", 'Q':"GLN", 'R':"ARG", 'S':"SER", 'T':"THR", 'V':"VAL", 'W':"TRP", 'Y':"TYR" }
AA_name3_to_name1 = { "ALA":'A', "CYS":'C', "ASP":'D', "GLU":'E', "PHE":'F', "GLY":'G', "HIS":'H', "ILE":'I', "LYS":'K', "LEU":'L', "MET":'M', "ASN":'N', "PRO":'P', "GLN":'Q', "ARG":'R', "SER":'S', "THR":'T', "VAL":'V', "TRP":'W', "TYR":'Y' }



def mutate_residue( pose_num, new_res_name, input_pose, sf, pdb_num = False, pdb_chain = None ):
    '''
    Mutate residue at position <pose_num> to <new_res_name>
    <new_res_name> can be a single-letter or three-letter residue code
    If you are giving a pdb number, set <pdb_num> to True AND give me a <pdb_chain> letter id
    :param pose_num: int( Pose number for residue )
    :param new_res_name: str( one- or three-letter code for the new amino acid. Example 'A' or "THR" )
    :param input_pose: Pose
    :param sf: ScoreFunction ( used for packing )
    :param pdb_num: bool( did you give me a PDB number instead? Set to True if so. Give me a <pdb_chain> too then ) Default = False (Pose number)
    :param pdb_chain: str( PDB chain id such as 'A' or 'X'. Must have set <pdb_num> to True as well
    :return: mutated Pose
    '''
    # imports
    from pyrosetta import Pose, pose_from_sequence
    from rosetta.core.conformation import ResidueFactory


    # copy over the input pose
    pose = input_pose.clone()

    # check if <pdb_chain> was given if <pdb_num> is True
    if pdb_num == True:
        if pdb_chain is None:
            print "\nYou told me you gave me a PDB number, but you did not give me a PDB chain id. Set <pdb_chain> to the appropriate chain id. Returning the original pose."
            return pose

    # check the logic of the input arguments
    # ensure the <new_res_name> is a valid residue
    if len( new_res_name ) != 1 and len( new_res_name ) != 3:
        print "\nYou did not give me a single- or three-letter amino acid code. '%s' did not work. Returning the original pose." %new_res_name
        return pose
    # if it is a single-letter code
    if len( new_res_name ) == 1 and new_res_name.upper() not in AA_name1_list:
        print "\nIt appears that '%s' is not a valid single-letter amino acid code. Returning the original pose." %new_res_name
        return pose
    # if it is a three-letter code
    elif len( new_res_name ) == 3 and new_res_name.upper() not in AA_name3_list:
        print "\nIt appears that '%s' is not a valid three-letter amino acid code. Returning the original pose." %new_res_name
        return pose
   # otherwise, use the <new_res_name> argument to get the appropriate three-letter amino acid code
    if len( new_res_name ) == 1:
        single_new_res_name = new_res_name.upper()
        new_res_name = AA_name1_to_name3[ single_new_res_name ]
    else:
        new_res_name = new_res_name.upper()
        single_new_res_name = AA_name3_to_name1[ new_res_name ]

    # ensure <pose_num> (and <pdb_chain>) exists in the pose
    if not pdb_num:
        if not 1 <= pose_num <= pose.size():
            print "\nYou appear to have given me an invalid Pose residue number. Ensure residue number %s exists in your Pose. Returning the original pose." %pose_num
            return pose
    # if it's a PDB number, check it exists as well using the <pdb_chain> too
    else:
        # get the actual pose number
        pose_num = pose.pdb_info().pdb2pose( pdb_chain, pose_num )
        if pose_num == 0:
            print "\nYour PDB number and chain ( %s chain %s ) don't seem to exist in the pose. Check your input. Returning the original pose." %( pose_num, pdb_chain )
            return pose

    # move on to the mutation
    # instantiate a ResidueFactory
    res_factory = ResidueFactory()

    # create a three-mer of the <new_res_name> desired
    # want a three-mer because it's easier to deal with a new amino acid that does not have a special end VariantType
    threemer = pose_from_sequence( single_new_res_name * 3 )

    # get the ResidueType from the middle <new_res_name> in the threemer
    res_type = threemer.conformation().residue_type( 2 )

    # build the new residue and preserve the CB information from the original pose
    new_residue = res_factory.create_residue( res_type, 
                                              current_rsd = pose.residue( pose_num ), 
                                              conformation = pose.conformation(), 
                                              preserve_c_beta = False )

    # replace the old residue in the pose
    pose.replace_residue( pose_num, new_residue, orient_backbone = True )

    return pose


def get_sum_hbond_E( sf, pose ):
    '''
    Get the sum of hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc
    :param sf: ScoreFunction
    :param pose: Pose
    :return: float( sum of hbond energies )
    '''
    # imports
    from rosetta.core.scoring import score_type_from_name
    
    # list of hbond energies to grab
    hbond_E_names = [ "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc" ]

    return sum( [ sf.score_by_scoretype( pose, score_type_from_name( st ) ) for st in hbond_E_names ] )
    #return sum( [ pose.energies().total_energies().get( score_type_from_name( n ) ) for n in hbond_E_names ] )


def show_score_breakdown( sf, pose ):
    '''
    Shows the breakdown of the <pose>'s total score by printing the score of each nonzero weighted ScoreType in <sf>
    :param sf: ScoreFunction
    :param pose: Pose
    '''
    # print out each score
    print "\n".join( [ "%s: %s" %( score_type, round( sf.score_by_scoretype( pose, score_type ), 3 ) ) for score_type in sf.get_nonzero_weighted_scoretypes() ] )
    print


def get_res_nums_within_radius( res_num_in, input_pose, radius, include_res_num = False ):
    """
    Use the nbr_atom_xyz to find residue numbers within <radius> of <pose_num> in <pose>
    The nbr_atom seems to be C4 on carbohydrates
    :param res_num_in: int( Pose residue number )
    :param input_pose: Pose
    :param radius: int or float( radius around <pose_num> to use to select resiudes )
    :param include_res_num: bool( do you want to include <res_num> in the return list? ) Default = False
    :return: list( Pose residue numbers within <radius> of <pose_num>
    """
    # clone the <input_pose>
    pose = input_pose.clone()

    # container for the centers of each residue in pose
    centers_of_res = []

    # fill up the centers container
    for res_num in range( 1, pose.size() + 1 ):
        center = pose.residue( res_num ).nbr_atom_xyz()
        centers_of_res.append( center )

    # container for residues inside the <radius>
    res_nums_in_radius = []

    # nbr_xyz of the residue of interest
    res_num_xyz = pose.residue( res_num_in ).nbr_atom_xyz()

    for res_num in range( 1, pose.size() + 1 ):
        # this will get the xyz of the residue of interest, but it will be removed from the final list if desired
        # (since it will be added as 0 will always be less than <radius>)
        # get the center of the residue
        center = pose.residue( res_num ).nbr_atom_xyz()

        # keep the residue number if the nbr_atom_xyz is less than <radius>
        if center.distance( res_num_xyz ) <= radius:
            res_nums_in_radius.append( res_num )

    # if the user didn't want the residue of interest in the return list, remove it
    if not include_res_num:
        res_nums_in_radius.remove( res_num_in )

    return res_nums_in_radius


def get_res_nums_within_radius_of_residue_list( residues, input_pose, radius, include_res_nums = False ):
    """
    Find all residue numbers around the list of <residues> given in <input_pose> within <radius> Angstroms.
    Set <include_residues> if you want to include the list of passed <residues> in the return list of residue numbers.
    Uses the nbr_atom to calculate distance. The nbr_atom seems to be C4 on carbohydrates
    :param residues: list( Pose residue numbers )
    :param input_pose: Pose
    :param radius: int() or float( radius in Angstroms )
    :param include_res_nums: bool( do you want to include the passed <residues> in the return list of resiude numbers? ) Default = False
    :return: list( residues around passed <residues> list within <radius> Angstroms
    """
    # argument check: ensure passed <residues> argument is a list
    if type( residues ) != list:
        print "\nArgument error. You're supposed to past me a list of residue numbers for the <residues> argument. Returning None."
        return None

    # use get_res_nums_within_radius to get all residue numbers
    residues_within_radius = []
    for res_num in residues:
        residues_within_radius.extend( get_res_nums_within_radius( res_num, input_pose, radius, include_res_num = include_res_nums ) )

    # get the set of the list and sort the residue numbers
    set_of_residues_within_radius = [ res for res in set( residues_within_radius ) ]

    # it is possible that there are still residues from <residues> in the list, so remove them one by one if desired
    if not include_res_nums:
        for res in residues:
            try:
                set_of_residues_within_radius.remove( res )
            except ValueError:
                pass

    # sort
    set_of_residues_within_radius.sort()

    return set_of_residues_within_radius
