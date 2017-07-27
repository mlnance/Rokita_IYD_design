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


def write_pdb_file(pdb_lines, filename):
    '''
    Write <pdb_lines> into a file named <filename>
    Include a path in <filename> if you don't want it in the current directory
    :param pdb_lines: list(pdb_lines)
    :param filename: str(/path/to/filename)
    :return: bool(True if successful)
    '''
    # the input pdb_lines may or may not end with a newline character
    # strip each line and add it, for safety
    pdb_lines = [line.strip() + "\n" for line in pdb_lines]

    # write the file
    try:
        with open(filename, 'w') as fh:
            fh.writelines(pdb_lines)
        return True
    except:
        return False


def get_xyz_coords_to_lines_dict(pdb_lines):
    '''
    Using lines from the pdb, pull out the x, y, z coordinates into a dictionary
    Coordinates for the atom of that line will be the key, and the corresponding line is the value
    :param pdb_lines: list(lines from pdb file)
    :return: dict[(x, y, z)] = pdb_line
    '''
    # set up the dictionary
    coords_to_pdb_line = {}
    # x, y, z coordinates are in the same location for each pdb file line
    # if a list of lines is given
    if type(pdb_lines) is list:
        for line in pdb_lines:
            line_xyz = (float(line[30:38].replace(' ', '')),
                        float(line[38:46].replace(' ', '')),
                        float(line[46:54].replace(' ', '')))
            coords_to_pdb_line[line_xyz] = line
    # if a single line is given
    if type(pdb_lines) is str:
        line = pdb_lines
        line_xyz = (float(line[30:38].replace(' ', '')),
                    float(line[38:46].replace(' ', '')),
                    float(line[46:54].replace(' ', '')))
        coords_to_pdb_line[line_xyz] = line

    return coords_to_pdb_line


def get_xyz_coords(pdb_lines):
    '''
    Using lines from the pdb, pull out the x, y, z coordinates into a dictionary
    Coordinates for the atom of that line will be the key, and the corresponding line is the value
    :param pdb_lines: list(lines from pdb file)
    :return: dict[(x, y, z)] = pdb_line
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
        line_xyz = (float(line[30:38].replace(' ', '')),
                    float(line[38:46].replace(' ', '')),
                    float(line[46:54].replace(' ', '')))
        xyz_coords = line_xyz

    return xyz_coords


def get_lines_from_xyz_coords(pdb_lines, xyz_coords):
    '''
    From a list of x, y, z coords in <xyz_coords>, find the corresponding line in <pdb_lines>
    :param pdb_lines: list(pdb_lines)
    :param xyz_coords: list([x, y, z])
    :return: list(corresponding pdb_lines)
    '''
    return_lines = []
    # for each pdb line
    for line in pdb_lines:
        # only look at ATOM and HETATM
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # for each set of xyz coordinates
            for xyz_set in xyz_coords:
                # pull out x, y, z
                x, y, z = xyz_set
                # pull out the line's x, y, z coords
                line_x, line_y, line_z = [float(line[30:38].replace(' ', '')),
                                          float(line[38:46].replace(' ', '')),
                                          float(line[46:54].replace(' ', ''))]
                # compare
                if x == line_x and y == line_y and z == line_z:
                    return_lines.append(line)

    return return_lines


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
    NOTE: Function also can handle FMN and 2IP, which don't have the CA atoms
    NOTE: If residue FMN (#401) and 2IP (#402) are passed, different atoms are used
    NOTE: For FMN, atom C7M is used. For 2IP, atom C4 is used
    :param pdb_lines: list(lines in the pdb file)
    :param residues: list(sets of residue numbers that define each plane)
    :return:
    '''
    # get a unique list of all residue numbers whose information is needed
    all_res_nums = []
    for l in residues:
        for r in l:
            all_res_nums.append(r)
    all_res_nums = list(set(all_res_nums))

    # construct a dictionary of residue number to CA coordinate
    CA_dict = {}
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # check if it's chain A
            chain = line[21:22]
            if chain == 'A':
                # check the residue number and if its a CA atom and chain A
                # also, if it's residue 401 or 402, use their corresponding atoms
                # residue 401 (FMN) atom C7M, residue 402 (2IP) atom C4
                res_num = int(line[22:26].replace(' ', ''))
                if res_num in all_res_nums:
                    # pull out the atom name
                    atom_name = line[12:16].replace(' ', '')
                    # check if it's residue 401 (FMN)
                    # and atom_name C7M
                    if res_num == 401 and atom_name == "C7M":
                        # add to dictionary
                        CA_dict[res_num] = get_xyz_coords(line)
                    # check if it's residue 402 (2IP)
                    # and atom_name C4
                    elif res_num == 402 and atom_name == "C4":
                        # add to dictionary
                        CA_dict[res_num] = get_xyz_coords(line)
                    # now that residues 401 and 402 have been handled
                    # check if it's any other protein residue in the list
                    # and check if it's a CA atom of a residue
                    elif res_num in all_res_nums and atom_name == "CA":
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
