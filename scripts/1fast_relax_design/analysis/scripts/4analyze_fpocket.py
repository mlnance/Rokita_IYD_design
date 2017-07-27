#!/usr/bin/python
__author__ = "morganlnance"

'''
Analysis of fpocket pdb data using python
MUST MANUALLY DEFINE convex_hull_residue_numbers, which is a set of PDB residue numbers
.csv data gets dumped to current working directory, or the optional --dump_dir location
.pdb files of all STP spheres within the convex hull for each pdb can be dumped, if desired
they would be dumped in the same directory that the pdb file from which the pdb file gets read
'''

###################
# GENERAL IMPORTS #
###################
import sys
import os
import argparse
try:
    import pandas as pd
    pandas_on = True
except ImportError:
    import csv
    pandas_on = False
# helper functions for analysis
from analyze_fpocket_util import pull_out_STP_lines, get_xyz_coords, \
    define_plane_three_points, distance_from_point_to_plane, \
    pdb_residues_to_convex_hull, get_lines_from_xyz_coords, \
    get_xyz_coords_to_lines_dict, write_pdb_file


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description="Use python to analyze fpocket data.")
parser.add_argument("pdb_list", type=str, help="/path/to/file containing paths to fpocket PDB files")
parser.add_argument("--dump_dir", type=str, help="/path/to/where you want to dump your data. \
                                                 Default is the current directory \
                                                 a .pdb file gets dumped for each input pdb \
                                                 and one .csv file for all data.")
parser.add_argument("--dump_STP_pdb", default=False, action="store_true",
                                    help="Do you want to dump a .pdb file \
                                    of all STP spheres within the convex hull \
                                    for each pdb in pdb_list? Default = False")
parser.add_argument("--dump_convex_hull_pdb", default=False, action="store_true",
                    help="Do you want to dump a .pdb file of the normal PDB with the STP spheres that are within the convex hull? Default = False")
input_args = parser.parse_args()

# ensure input arguments are valid
# make sure pdb_list is an actual file
if not os.path.isfile(input_args.pdb_list):
    print "\nYour pdb_list argument is invalid.\n"
    sys.exit()
else:
    pdb_list = input_args.pdb_list
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


##############
# READ FILES #
##############
# read the pdb_list file
with open(pdb_list, 'r') as fh:
    pdb_files = fh.readlines()


#############
# READ PDBs #
#############
pdb_names = []
all_num_stp_spheres_in_convex_hull = []
# for each pdb file in the list
for pdb_file in pdb_files:
    # read in the pdb file
    # skip the pdb if for any reason opening and reading doesn't work
    try:
        with open(pdb_file.strip(), 'r') as fh:
            pdb_lines = fh.readlines()
    except:
        continue
    # add this pdb name to the pdb_names holder
    pdb_name = pdb_file.split('/')[-1].split(".gz")[0].split(".pdb")[0].split("_out")[0]
    # store the directory that this pdb is in for use later
    pdb_dir = '/'.join(pdb_file.split('/')[:-1]) + '/'

    ######################
    # DEFINE CONVEX HULL #
    ######################
    # predetermined residue numbers (in PDB numbering) are used to
    # create a convex hull using the coordinates of each residue's CA
    # or, for the case of FMN and 2IP (401 and 402), atoms C7M and C4, respectively
    convex_hull_residue_numbers = [[90, 98, 114], [90, 98, 173],
                                   [90, 114, 173], [98, 114, 402],
                                   [98, 173, 401], [98, 401, 402],
                                   [114, 173, 401], [114, 401, 402]]
    planes = pdb_residues_to_convex_hull(pdb_lines=pdb_lines,
                                         residues=convex_hull_residue_numbers)

    ###################
    # GET STP SPHERES #
    ###################
    stp_sphere_lines = pull_out_STP_lines(pdb_lines)
    # also keep all lines that ARENT STP sphere lines
    # not using a set because that puts the lines out of order
    # keeping these lines separate in order to make a PDB of the structure
    # with the STP spheres and lines between the convex hull residues
    normal_pdb_lines = [line for line in pdb_lines if line not in stp_sphere_lines]

    ##################
    # GET STP COORDS #
    ##################
    # this function returns a dictionary
    # :return: dict[(x, y, z)] = pdb_line
    stp_sphere_xyz_coords_to_lines_dict = get_xyz_coords_to_lines_dict(pdb_lines=stp_sphere_lines)

    ##########
    # PLANES #
    ##########
    # for each stp point
    # for each plane
    # for each vertex of every other plane
    # ensure same direction of vertex as stp point
    # repeat
    # for all STP water spheres, determine their direction with respect to each plane
    stp_spheres_inside_convex_hull = []
    # the stp_sphere_xyz_coords_to_lines_dict has as its keys the x,y,z coords of the STP spheres
    for stp_point in stp_sphere_xyz_coords_to_lines_dict.keys():
        # for each ii_plane, compare the stp_direction to each vertex
        # of each jj_plane compared to ii_plane
        ii_to_jj_directions = []
        for ii in range(len(planes)):
            ############
            # ii PLANE #
            ############
            # each coordinate set is a CA atom of a residue
            # [(x1, y1, y2), (x2, y2, z2), (x3, y3, z3)]
            ii_plane = planes[ii]

            ###################
            # ii PLANE POINTS #
            ###################
            # three 3D points define the plane of ii_plane
            # p1: (x1, y1, y2), p2: (x2, y2, z2), p3: (x3, y3, z3)
            ii_p1, ii_p2, ii_p3 = ii_plane

            #########################
            # ii PLANE COEFFICIENTS #
            #########################
            # use the plane points to get the variables that describe the plane's equation
            # ax + by + cz + d = 0
            a, b, c, d = define_plane_three_points(p1=ii_p1, p2=ii_p2, p3=ii_p3)

            #######################
            # STP POINT DIRECTION #
            #######################
            # get the direction of an STP sphere from the plane
            # ax + by + cz + d = a "direction" (positive or negative)
            # use the x,y,z coordinates of the STP sphere
            x, y, z = stp_point
            stp_direction = a*x + b*y + c*z + d
            # we only care about if the direction is positive or negative
            # True == positive, False == negative
            stp_direction = True if stp_direction > 0 else False

            #####################################
            # COMPARE DIRECTION TO OTHER PLANES #
            #####################################
            # the stp_direction should match the direction of every other vertex
            # of every other plane from the current plane
            # every vertex on every plane that isn't ii_plane must be compared
            # in terms of direction with respect to ii_plane
            # so, every point defining all jj_planes are compared to the ii_plane
            # equation to result in a direction. each direction must be the same as
            # the stp_sphere's direction with respect to ii_plane. if all directions
            # (leeway of 2) are the same for all planes, then that STP sphere
            # is within the convex hull
            for jj in range(len(planes)):
                ############
                # jj PLANE #
                ############
                # don't compare the same plane vertices
                # though other planes will share the same vertices
                # this case will be handled below
                if ii != jj:
                    # extract the jjth plane
                    # for every vertex point of every plane
                    # [(x1, y1, y2), (x2, y2, z2), (x3, y3, z3)]
                    jj_plane = planes[jj]

                    #############################
                    # jj PLANE POINTS DIRECTION #
                    #############################
                    # for each vertex of jj_plane
                    for jj_point in jj_plane:
                        # only check direction for vertices that are not also in ii_plane
                        # ii_plane and jj_plane are not the same plane, but they may share vertices
                        # meaning a vertex in ii_plane may also be a vertex in jj_plane
                        # so check that the point being considered in jj_plane is not in ii_plane
                        if jj_point != ii_p1 and jj_point != ii_p2 and jj_point != ii_p3:
                            # pull out x, y, z from the a vertex of jj_plane
                            jj_x, jj_y, jj_z = jj_point
                            # get the direction of jj_plane_point using ii_plane's plane equation coefficients
                            jj_direction = a*jj_x + b*jj_y + c*jj_z + d
                            # we only care about if the direction is positive or negative
                            # True == positive, False == negative
                            jj_direction = True if jj_direction > 0 else False

                            # compare the jj_plane direction with respect to ii_plane to
                            # the direction of the stp_point with respect to the ii_plane
                            ii_to_jj_directions.append(stp_direction == jj_direction)

        ###############################
        # DETERMINE IF IN CONVEX HULL #
        ###############################
        # use the set of comparisons of directions to determine
        # if the STP sphere is in the convex hull defined by planes
        # since the convex hull chosen is not perfectly convex,
        # allow a max of two falses (meaning jj_direction != stp_direction)
        # if 2 or fewer falses in comparison, STP sphere is in the convex hull
        if ii_to_jj_directions.count(False) <= 7:
            line = get_lines_from_xyz_coords(pdb_lines, [stp_point])[0]
            stp_spheres_inside_convex_hull.append(line)

    ################
    # COLLECT DATA #
    ################
    # we only want the unique STP spheres that were inside the convex hull
    stp_spheres_inside_convex_hull = list(set(stp_spheres_inside_convex_hull))
    num_stp_spheres_in_convex_hull = len(stp_spheres_inside_convex_hull)

    # add to lists
    pdb_names.append(pdb_name)
    all_num_stp_spheres_in_convex_hull.append(num_stp_spheres_in_convex_hull)

    # write out a pdb file of the unique STP spheres, if desired
    if input_args.dump_STP_pdb:
        write_pdb_file(pdb_lines = stp_spheres_inside_convex_hull,
                       filename = pdb_dir + pdb_name + "_STP_spheres_in_convex_hull.pdb")
    # write out a pdb file of the normal PDB with the STP spheres in the convex hull, if desired
    if input_args.dump_convex_hull_pdb:
        file_lines = normal_pdb_lines + stp_spheres_inside_convex_hull
        write_pdb_file(pdb_lines = file_lines, 
                       filename = pdb_dir + pdb_name + '_convex_hull.pdb')


##############
# WRITE DATA #
##############
# all unique STP spheres within the defined convex hull have been found
# specifically only looking at chain A, though, but bact_IYD is symmetrical
# compile the data into a dataframe
# use pandas if it is available
if pandas_on:
    # store collected data in a pandas dataframe
    df = pd.DataFrame()
    df["pdb_name"] = pdb_names
    df["num_STP_spheres"] = all_num_stp_spheres_in_convex_hull

    print(df)
    df.to_csv(dump_dir + "fpocket_analysis.csv")
# otherwise, use python's csv module
else:
    # zip the data together to write row-by-row
    data = zip(pdb_names, all_num_stp_spheres_in_convex_hull)
    with open(dump_dir + "fpocket_analysis.csv", 'w') as fh:
        writer = csv.writer(fh)
        header = ["pdb_name", "num_STP_spheres"]
        print header
        writer.writerow(header)
        for line in data:
            print line
            writer.writerow(line)
