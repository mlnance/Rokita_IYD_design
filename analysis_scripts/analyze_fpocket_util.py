#!/usr/bin/python
__author__ = "morganlnance"

'''
Analysis functions using python
'''


def pull_out_STP_lines(pdb_lines):
    '''
    Using residue names, use <pdb_lines> to find all STP residues
    :param pdb_lines: list(lines from pdb file)
    :return: list(lines of STP residues)
    '''
    # go through each line, keeping lines of STP residues
    stp_lines = []
    for line in pdb_lines:
        # residue names are at this position in the PDB
        if line[17:20] == "STP":
            stp_lines.append(line)

    return stp_lines


def get_xyz_coords(pdb_lines):
    '''
    Using lines from the pdb, pull out the x, y, z coordinates into a list
    :param pdb_lines: list(lines from pdb file)
    :return: list( (x, y, z) )
    '''
    # x, y, z coordinates are in the same location for each pdb file line
    # if a list of lines is given
    if type(pdb_lines) is list:
        xyz_coords = []
        for line in pdb_lines:
            line_xyz = (float(line[30:38].replace(' ', '')),
                        float(line[38:46].replace(' ', '')),
                        float(line[46:54].replace(' ', '')))
            xyz_coords.append(line_xyz)
    # if a single line is given
    if type(pdb_lines) is str:
        line = pdb_lines
        xyz_coords = (float(line[30:38].replace(' ', '')),
                      float(line[38:46].replace(' ', '')),
                      float(line[46:54].replace(' ', '')))

    return xyz_coords


def define_plane_three_points(p1, p2, p3):
    '''
    Use three 3D points <x>, <y>, <z> to define a plane equation
    :param p1: tuple(x1, y1, z1)
    :param p2: tuple(x2, y2, z2)
    :param p3: tuple(x3, y3, z3)
    :return: floats(a, b, c, d defining a plane as ax + by + cz - d)
    '''
    # information source
    # steps/code: http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
    # math info: https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FEquationofaPlane3Points
    # imports
    import numpy as np

    # three points
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1

    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp

    # This evaluates a * x3 + b * y3 + c * z3 which equals d
    # p3 is arbitrary
    # -1 * d because the person who wrote this code has d on the right of the =
    d = -1 * np.dot(cp, p3)

    #print('The equation is {0}x + {1}y + {2}z + {3} = 0'.format(a, b, c, d))
    return (a, b, c, d)


def distance_from_point_to_plane(p1, a, b, c, d):
    '''
    Using the plane parameters <a>, <b>, <c>, <d>, calculate the distance from the plane to point <p1>
    :param p1: tuple(x1, y1, z1)
    :param a: float
    :param b: float
    :param c: float
    :param d: float
    :return: float(distance from p1 to plane ax1 + by1 + cz1 + d = 0)
    '''
    # information source
    # math info: http://www.songho.ca/math/plane/plane.html
    # imports
    import numpy as np
    from math import sqrt

    #                   numerator                   denominator
    # distance = ( ax1 + by1 + cz1 + d ) / ( srt( a**2 + b**2 + c**2 ) )
    numerator = (a*p1[0]) + (b*p1[1]) + (c*p1[2]) + d
    denominator = sqrt(a**2 + b**2 + c**2)

    return numerator / denominator


def pdb_residues_to_convex_hull(pdb_lines, residues):
    '''
    For a pdb described in <pdb_lines>, use <residues> to pull out CA atoms that will describe a set of planes
    Example: residues = [ [1, 2, 3], [3, 5, 6], [6, 8, 9] ]
    Where each set of three residues numbers defines one plane using their CA coordinate
    :param pdb_lines: list(lines in the pdb file)
    :param residues: list(sets of residue numbers that define each plane)
    :return:
    '''
    # get a unique list of all residue numbers whose information is needed
    all_res_nums = []
    for l in residues:
        for ii in l:
            all_res_nums.append(ii)
    all_res_nums = list(set(all_res_nums))

    # construct a dictionary of residue number to CA coordinate
    CA_dict = {}
    for line in pdb_lines:
        if line.startswith("ATOM"):
            # check if it's chain A
            chain = line[21:22]
            if chain == 'A':
                # check if it's a CA atom of a residue
                atom_name = line[12:16].replace(' ', '')
                if atom_name == "CA":
                    # check the residue number and if its a CA atom and chain A
                    res_num = int(line[22:26].replace(' ', ''))
                    if res_num in all_res_nums:
                        # add to dictionary
                        CA_dict[res_num] = get_xyz_coords(line)

    # now that the CA coordinate dictionary has been filled
    # with the CA coordinates of only the residues of interest
    # construct the coordinates of each specified plane
    planes = []
    for plane_residues in residues:
        plane = [ CA_dict[res] for res in plane_residues ]
        planes.append(plane)

    return planes