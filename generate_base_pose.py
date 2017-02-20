#!/usr/bin/python
__author__="morganlnance"



#################################
#### LOAD IN PDB WITH PARAMS ####
#################################
def generate_base_pose( pdb_file, rounds, params_dir=None ):
    '''
    Takes an input pdb_file or a Pose object and generates a base conformation by running four rounds of packing and minimization
    Can be run as an independent script or imported as a function
    :param pdb_file: str( /path/to/pdb_file if script ) or Pose() if function
    :param rounds: int( number of structures to generate in order to find the best one
    :param params_dir: str( /path/to/dir with parameter files needed, if run as script
    :return: Pose( low_E_pose ), Pandas.DataFrame() containing data
    '''
    # imports
    from rosetta.protocols.simple_moves import \
        RotamerTrialsMover, MinMover
    from rosetta.core.scoring import fa_atr, fa_rep
    from pyrosetta import PyMOLMover, get_fa_scorefxn, \
        standard_packer_task, MoveMap, Pose, \
        generate_nonstandard_residue_set, \
        pose_from_file
    import util
    import pandas as pd
    import sys, os


    # if params_dir is given
    if params_dir is not None:
        # load in params to enable processing of IYD with cofactor and substrate
        # organize the parameter directory and pull param files from there
        if not params_dir.endswith( '/' ):
            params_dir = params_dir + '/'
        else:
            params_dir = params_dir
        params = [ os.path.join( params_dir, param ) for param in os.listdir( params_dir ) ]
        print "\ngenerate_nonstandard_residue_set"
        orig_pose = Pose()
        try:
            nonstandard_res_set = generate_nonstandard_residue_set( orig_pose, params )
        except:
            print "\nThere is something wrong with your parameter files. Did you give me the proper params directory?"
            print "*** %s\n" %params_dir
            sys.exit()
        print "\npose_from_file"
        try:
            pose_from_file( orig_pose, nonstandard_res_set, pdb_file )
        except:
            print "\nThere was some error loading your PDB. Is this a valid PDB file?: %s\n" %pdb_file
            sys.exit()

        # checks pass, load pose
        orig_pose.pdb_info().name( "orig_pose" )

    # otherwise, a pose should have been given
    else:
        try:
            orig_pose = pdb_file.clone()
        except:
            "\nYou ran this as a function, so I expected a Pose object. Check your input\n"
            sys.exit()

    # create a PyMOLMover
    pmm = PyMOLMover()

    # generate a score function
    print "\nget_fa_scorefxn"
    sf = get_fa_scorefxn()
    orig_E = sf( orig_pose )

    # get orig weights of fa_atr and fa_rep
    # these weights will be adjusted during pack/min
    orig_fa_atr = sf.get_weight( fa_atr )
    orig_fa_rep = sf.get_weight( fa_rep )

    # create a Pandas DataFrame and data holders
    df = pd.DataFrame()
    orig_E_vals = []
    base_E_vals = []
    dG_vals = []
    poses = []


    ###########################
    #### PACK AND MINIMIZE ####
    ###########################
    # create a packer task
    print "\nstandard_packer_task"
    task = standard_packer_task( orig_pose )
    task.or_include_current( True )
    task.restrict_to_repacking()

    # create a MoveMap used for minimization
    print "\nMoveMap"
    mm = MoveMap()
    mm.set_bb( True )
    mm.set_chi( True )
    mm.set_jump( True )

    # set up for storing the best E pose
    best_pose = None
    times_changed = 0

    # for as many decoy poses as specified
    for round in range( 1, rounds + 1 ):
        # get a fresh copy of the original pose
        pose = orig_pose.clone()
        pose.pdb_info().name( "pose_%s" %str( round ) )

        # create a RotamerTrialsMover
        rtm = RotamerTrialsMover( sf, task )

        # create a MinMover using the MoveMap
        minmover = MinMover( movemap_in=mm, 
                             scorefxn_in=sf, 
                             min_type_in="lbfgs_armijo_nonmonotone", 
                             tolerance_in=0.01, 
                             use_nb_list_in=True, 
                             deriv_check_in=False, 
                             deriv_check_verbose_in=False )

        # run multiple RotamerTrialsMover and MinMover rounds
        # 1
        print "\nRotamerTrialsMover"
        # multiply up fa_atr and down fa_rep
        sf.set_weight( fa_atr, orig_fa_atr * 2.0 )
        sf.set_weight( fa_rep, orig_fa_rep * 0.5 )
        rtm.score_function( sf )
        rtm.apply( pose )
        print "\nMinMover"
        minmover.score_function( sf )
        minmover.apply( pose )
        # 2
        # then ramp fa_atr back down and fa_rep back up
        # over each pack/min trial
        sf.set_weight( fa_atr, orig_fa_atr * 1.67 )
        sf.set_weight( fa_rep, orig_fa_rep * 0.67 )
        print "\nRotamerTrialsMover"
        rtm.score_function( sf )
        rtm.apply( pose )
        print "\nMinMover"
        minmover.score_function( sf )
        minmover.apply( pose )
        # 3
        sf.set_weight( fa_atr, orig_fa_atr * 1.33 )
        sf.set_weight( fa_rep, orig_fa_rep * 0.83 )
        print "\nRotamerTrialsMover"
        rtm.score_function( sf )
        rtm.apply( pose )
        print "\nMinMover"
        minmover.score_function( sf )
        minmover.apply( pose )
        # 4
        sf.set_weight( fa_atr, orig_fa_atr )
        sf.set_weight( fa_rep, orig_fa_rep )
        print "\nRotamerTrialsMover"
        rtm.score_function( sf )
        rtm.apply( pose )
        print "\nMinMover"
        minmover.score_function( sf )
        minmover.apply( pose )

        # is this the first pose?
        # does this pose have the lowest E seen?
        if best_pose is None or ( sf( pose ) < sf( best_pose ) ):
            best_pose = pose.clone()
            best_pose_E = sf( best_pose )
            times_changed += 1

        # relay energy of pose
        base_E = sf( pose )
        print "\n\nround %s of %s" %( round, rounds )
        print "orig energy:", orig_E
        print "current energy:", base_E
        print "best seen energy:", best_pose_E
        print "times best_pose changed:", times_changed

        # store pose and energies in DataFrame lists
        orig_E_vals.append( orig_E )
        base_E_vals.append( base_E )
        dG_vals.append( base_E - orig_E )
        poses.append( pose )

    # add the data to the DataFrame
    df["orig_E"] = orig_E_vals
    df["base_E"] = base_E_vals
    df["dG"] = dG_vals
    df["poses"] = poses
    df = df.sort( "dG" )

    return best_pose, df



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Use PyRosetta to make a low-energy base structure of IYD")
    parser.add_argument("pdb_file", type=str, help="the path to the PDB file of which you want to create a low-E version")
    parser.add_argument("params_dir", type=str, help="the path to the directory holding the relevant parameter files")
    parser.add_argument("rounds", type=int, help="the number of base poses you want to generate to find the lowest energy PDB")
    input_args = parser.parse_args()

    from pyrosetta import init
    init()
    generate_base_pose( input_args.pdb_file, input_args.rounds, params_dir=input_args.params_dir )
