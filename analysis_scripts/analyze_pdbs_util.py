#!/usr/bin/python
__author__ = "morganlnance"

'''
Analysis functions using PyRosetta4
'''


def get_sequence(pose, res_nums=None):
    # type: (Pose, list) -> str
    """
    Return the sequence of the <pose>, or, return the sequence listed in <res_nums>
    :param pose: Pose
    :param res_nums: list() of Pose residue numbers
    :return: str(Pose sequence)
    """
    # if no res_nums were given, return the pose's sequence
    if res_nums is None:
        return str(pose.sequence())
    # else, return the sequence of the specified res_nums
    else:
        return str(''.join([pose.residue(r).name1() for r in res_nums]))


def get_atom_pair_distance(pose, res1, atom1, res2, atom2):
    """
    Get the xyz distance between <atom1> of <res1> to <atom2> in <res2> in the <pose>
    :param pose: Pose
    :param res1: int(residue number)
    :param atom1: int(atom number)
    :param res2: int(residue number)
    :param atom2: int(atom number)
    :return: float(xyz distance)
    """
    # pull out the atom objects from pose
    atom1 = pose.residue(res1).atom(atom1)
    atom2 = pose.residue(res2).atom(atom2)

    # calculate and return the distance between atom1 and atom2
    return float(atom1.xyz().distance(atom2.xyz()))


def write_fasta_file(pdb_names, pdb_sequences, filename, dump_dir=''):
    """
    Use a list of <pdb_names> and their corresponding <pdb_sequences> to write out a FASTA formatted file
    Need a <filename> to work with. Include a path to a dump directory, if desired
    :param pdb_names: list(pdb names)
    :param pdb_sequences: list(pdb sequences)
    :param filename: str(filename)
    :return: Bool
    """
    # ensure that the pdb_names and pdb_sequences lists are the same length
    if len(pdb_names) != len(pdb_sequences):
        return False

    # add .txt to the filename, if needed
    if not filename.endswith(".txt"):
        filename += ".txt"

    # write out the fasta file
    with open(filename, 'w') as fh:
        for pdb_name, pdb_seq in zip(pdb_names, pdb_sequences):
            fh.write(">%s\n%s\n") %(pdb_name, pdb_seq)
