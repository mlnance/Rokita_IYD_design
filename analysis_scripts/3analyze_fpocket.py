#!/usr/bin/python
__author__ = "morganlnance"

'''
Analysis of fpocket pdb data using python
MUST MANUALLY DEFINE convex_hull_residue_numbers, which is a set of PDB residue numbers
'''

###################
# GENERAL IMPORTS #
###################
import sys
import os
import argparse
# helper functions for analysis
from analyze_fpocket_util import pull_out_STP_lines, get_xyz_coords, \
    define_plane_three_points, distance_from_point_to_plane, \
    pdb_residues_to_convex_hull


#############
# ARGUMENTS #
#############
parser = argparse.ArgumentParser(description="Use python to analyze fpocket data.")
parser.add_argument("pdb_list", type=str, help="/path/to/file containing paths to fpocket PDB files")
input_args = parser.parse_args()

# ensure input arguments are valid
# make sure pdb_list is an actual file
if not os.path.isfile(input_args.pdb_list):
    print "\nYour pdb_list argument is invalid.\n"
    sys.exit()
else:
    pdb_list = input_args.pdb_list


##############
# READ FILES #
##############
# read the pdb_list file
with open(pdb_list, 'r') as fh:
    pdb_files = fh.readlines()


#############
# READ PDBs #
#############
# for each pdb file in the list
for pdb_file in pdb_files:
    # read in the pdb file
    # skip the pdb if for any reason opening and reading doesn't work
    try:
        with open(pdb_file.strip(), 'r') as fh:
            pdb_lines = fh.readlines()
    except:
        continue

    ######################
    # DEFINE CONVEX HULL #
    ######################
    # predetermined residue numbers (in PDB numbering) are used to
    # create a convex hull using the coordinates of each residue's CA
    convex_hull_residue_numbers = [[95, 99, 103],
                                   [95, 99, 112]]
    planes = pdb_residues_to_convex_hull(pdb_lines=pdb_lines,
                                         residues=convex_hull_residue_numbers)

    ###################
    # GET STP SPHERES #
    ###################
    stp_sphere_lines = pull_out_STP_lines(pdb_lines)

    ##################
    # GET STP COORDS #
    ##################
    stp_sphere_xyz_coords = get_xyz_coords(stp_sphere_lines)
    stp_sphere_xyz_coords = stp_sphere_xyz_coords[:2]

    ##########
    # PLANES #
    ##########
    # for each set of points that define a plane of interest
    #all_plane_points = convex_hull_planes_df.values.tolist()
    # for each stp point
    # for each plane
    # for each vertex of every other plane
    # ensure same direction of vertex as stp point
    # repeat
    # for all STP water spheres, determine their direction with respect to each plane
    for stp_point in stp_sphere_xyz_coords:
        # for each plane
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
            # not considering cases where a point lies on the plane, for clarity
            # True == positive, False == negative
            stp_direction = True if stp_direction > 0 else False

            #####################################
            # COMPARE DIRECTION TO OTHER PLANES #
            #####################################
            # the stp_direction should match the direction of every other vertex
            # of every other plane from the current plane
            for jj in range(len(planes)):
                ############
                # jj PLANE #
                ############
                # don't compare the same plane vertices
                if ii != jj:
                    # extract the jjth plane
                    # for every vertex point of every plane
                    # [(x1, y1, y2), (x2, y2, z2), (x3, y3, z3)]
                    jj_plane = planes[jj]

                    #############################
                    # jj PLANE POINTS DIRECTION #
                    #############################
                    # get direction of each vertex point of jj_plane with respect to ii_plane
                    # each jj_plane_point is (x, y, z) of a vertex on jj_plane
                    for jj_plane_point in jj_plane:
                        # extract x, y, z
                        jj_x, jj_y, jj_z = jj_plane_point
                        # get the direction of jj_plane_point using ii_plane's plane equation coefficients
                        jj_direction = a*jj_x + b*jj_y + c*jj_z + d
                        # we only care about if the direction is positive or negative
                        # not considering cases whejre a point lies on the plane, for clarity
                        # True == positive, False == negative
                        jj_direction = True if jj_direction > 0 else False

                        print stp_direction == jj_direction
