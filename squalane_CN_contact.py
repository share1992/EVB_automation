import openbabel
import numpy as np
import pandas as pd
import os
from scipy.spatial.distance import cdist
import MDAnalysis as md

# Dictionary of atom INDEXES from VMD. C: (H, (H, (H)))
all_atoms = {0: (30, 31, 32),
             2: (34, 35, 36),
             7: (44, 45, 46),
             12: (54, 55, 56),
             18: (66, 67, 68),
             23: (76, 77, 78),
             28: (86, 87, 88),
             29: (89, 90, 91),
             3: (37, 38),
             4: (39, 40),
             5: (41, 42),
             8: (47, 48),
             9: (49, 50),
             10: (51, 52),
             13: (57, 58),
             14: (59, 60),
             15: (61, 62),
             16: (63, 64),
             19: (69, 70),
             20: (71, 72),
             21: (73, 74),
             24: (79, 80),
             25: (81, 82),
             26: (83, 84),
             1: 33,
             6: 43,
             11: 53,
             17: 65,
             22: 75,
             27: 85}

def read_pdb(path):
    """ Reads in an pdb file from path as a DataFrame. This DataFrame is then turned into a 3D array such that the
    dimensions are (number of points) X (number of atoms) X 3 (Cartesian coordinates). The system name (based on the
    filename), list of atoms in the system, and Cartesian coordinates are output.
    :param path: path to pdb file to be read
    :return frames_split: numpy array composed of individual frames as pandas DataFrames
    """
    data = pd.read_csv(path, skiprows=1, header=None, delim_whitespace=True)
    data = data[data[0] != 'END']
    data.columns = ['type', 'global index', 'local index', 'fragment', 'N', 'fragment index', 'X', 'Y', 'Z', 'Nan',
                    'Nan', 'Nan']
    n_atoms = int(data['global index'].max())
    n_frames = int(len(data)/n_atoms)

    frames_split = []
    for i in range(int(n_frames)):
        start = n_atoms * i + i
        end = start + n_atoms - 1
        single_frame_df = data.loc[start:end]
        frames_split.append(single_frame_df)

    return frames_split


def get_atoms_within_threshold_distance(traj_frames_df, threshold_distance):

    C_index = traj_frames_df[0].index[(traj_frames_df[0]['fragment'] == 'CYA') & (traj_frames_df[0]['local index'] == 'C1')][0]
    N_index = traj_frames_df[0].index[(traj_frames_df[0]['fragment'] == 'CYA') & (traj_frames_df[0]['local index'] == 'N2')][0]

    atom_indexes = []
    for i in range(len(traj_frames_df)):
        # Distances between CN C and all other atom_indexes for first frame
        distances = cdist(traj_frames_df[i].iloc[C_index, 6:9][None, :], traj_frames_df[i].iloc[:, 6:9])
        indexes_within_threshold_distance = np.argwhere(distances < threshold_distance)[:, 1]

        for index in indexes_within_threshold_distance:
            if int(index) != C_index and int(index) != N_index:
                if index not in atom_indexes:
                    atom_indexes.append(index)

    return atom_indexes


def get_relevant_indices_H_atoms_only(traj_frames_df, atom_indexes):

    frame_1 = traj_frames_df[0]
    local_indexes = []
    fragments = []

    for atom_index in atom_indexes:
        local_index = frame_1.loc[atom_index]['local index']
        local_index_split = list(local_index)

        if local_index_split[0] == 'H':
            fragment_index = int(frame_1.loc[atom_index]['fragment index'])
            local_index_number = "".join([local_index_split[1], local_index_split[2]])

            local_indexes.append(int(local_index_number))
            fragments.append(int(fragment_index))

    d = {'local index': local_indexes, 'fragment index': fragments}
    relevant_indices_and_fragments_df = pd.DataFrame(d)

    return relevant_indices_and_fragments_df


def generate_patch_defs(H_atoms, atom_dictionary, squalane_mol2_file, skeleton_file='squalane-patches_skeleton.txt'):
    """
    Print necessary patch definitions for EVB parameters
    :param H_atoms: List of H atoms that get within threshold distance of CN carbon
    :return:
    """
    with open(skeleton_file) as skeletonfile, open('squalane_patch_defs_all.txt', 'w') as outputfile:
        outputfile.writelines(skeletonfile)

    outputfile = open('squalane_patch_defs_all.txt', 'a')

    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("mol2", "mol2")

    mol = openbabel.OBMol()

    obconversion.ReadFile(mol, squalane_mol2_file)

    atom_table = {}
    for atom in openbabel.OBMolAtomBFSIter(mol):
        atom_table[atom.GetIndex()] = atom.GetType()

    for H_atom in H_atoms:
        # To get numbering to start at 0
        H_atom = H_atom - 1
        if list(atom_table[H_atom])[0] != 'H':
            print("ERROR: Atom %s is not an H atom." % H_atom)
            break
        else:
            outputfile.write("PRES RH%s \t ! PATCH CN + squalane -> HCN + squalyl (missing H%s)" % ((H_atom+1), (H_atom+1)))
            for C_atom, connected_H_atoms in sorted(atom_dictionary.items()):
                if type(connected_H_atoms) == tuple and H_atom in connected_H_atoms:
                    outputfile.write("\nATOM 1C%s\t CSP2\t 0.00" % (C_atom+1))
                elif type(connected_H_atoms) == int and H_atom == connected_H_atoms:
                    outputfile.write("\nATOM 1C%s\t CSP2\t 0.00" % (C_atom+1))
                else:
                    outputfile.write("\nATOM 1C%s\t CR\t 0.00" % (C_atom+1))
            for index, atom_type in atom_table.items():
                if 'H' == atom_type:
                    outputfile.write("\nATOM 1H%s\t HC\t 0.00" % (int(index)+1))

            outputfile.write("\nATOM 2C1\t CSP\t 0.00")
            outputfile.write("\nATOM 2N2\t NSP\t 0.00")

            for C_atom, connected_H_atoms in atom_dictionary.items():
                if type(connected_H_atoms) == tuple and H_atom in connected_H_atoms:
                    outputfile.write("\nDELETE BOND 1C%s 1H%s" % ((C_atom+1), (H_atom+1)))
                    outputfile.write("\nBOND 2C1 1H%s" % (H_atom+1))
                elif type(connected_H_atoms) == int and H_atom == connected_H_atoms:
                    outputfile.write("\nDELETE BOND 1C%s 1H%s" % ((C_atom+1), (H_atom+1)))
                    outputfile.write("\nBOND 2C1 1H%s" % (H_atom+1))

            dihedrals_to_delete = []
            for C_atom, connected_H_atoms in atom_dictionary.items():
                for obtorsion in openbabel.OBMolTorsionIter(mol):
                    if H_atom in obtorsion and C_atom in obtorsion and obtorsion not in dihedrals_to_delete:
                        dihedrals_to_delete.append(obtorsion)
                        outputfile.write("\nDELETE DIHE 1%s%s 1%s%s 1%s%s 1%s%s" % (list(atom_table[obtorsion[0]])[0],
                                                                     (obtorsion[0] + 1),
                                                                     list(atom_table[obtorsion[1]])[0],
                                                                     (obtorsion[1] + 1),
                                                                     list(atom_table[obtorsion[2]])[0],
                                                                     (obtorsion[2] + 1),
                                                                     list(atom_table[obtorsion[3]])[0],
                                                                     (obtorsion[3] + 1)))
            outputfile.write("\n\n")


def print_input_file_lines(relevant_indices_and_fragments):
    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        index_1 = relevant_indices_and_fragments['local index'][i]
        index_2 = relevant_indices_and_fragments['fragment index'][i]

        print("if @NODE .eq. %s PATCH RH%s A %s A 101 setup" % ((i + 1), index_1, index_2))

        # Maximum of 16 nodes possible
        if i >= 14:
            break

    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        print("if @NODE .eq. %s AUTO ANGL DIH" % (i + 1))

        # Maximum of 16 nodes possible
        if i >= 14:
            break

    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        print("if @NODE .eq. %s print psf" % (i + 1))

        # Maximum of 16 nodes possible
        if i >= 14:
            break


def print_patch_lines(relevant_indices_and_fragments):

    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        index_1 = relevant_indices_and_fragments['local index'][i]
        index_2 = relevant_indices_and_fragments['fragment index'][i]

        print("if @NODE .eq. %s PATCH RH%s A %s A 101 setup" % ((i + 1), index_1, index_2))

        # Maximum of 16 nodes possible
        if i >= 14:
            break


def print_autogeneration_lines(relevant_indices_and_fragments):

    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        print("if @NODE .eq. %s AUTO ANGL DIH" % (i + 1))

        # Maximum of 16 nodes possible
        if i >= 14:
            break


def print_psf_lines(relevant_indices_and_fragments):
    print("\n")
    for i in range(len(relevant_indices_and_fragments)):
        print("if @NODE .eq. %s print psf" % (i + 1))

        # Maximum of 16 nodes possible
        if i >= 14:
            break


def print_patch_lines_to_input_file(relevant_indices_and_fragments, template_file_path):

    with open(template_file_path) as template_file, open('test_input.inp', 'w') as output_file:
        output_file.writelines(template_file)

    with open('test_input.inp') as f:
        content = f.readlines()

    patch_line = content.index("INSERT NECESSARY PATCHES HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        index_1 = relevant_indices_and_fragments['local index'][i]
        index_2 = relevant_indices_and_fragments['fragment index'][i]

        line_to_write = "if @NODE .eq. %s PATCH RH%s A %s A 101 setup" % ((i + 1), index_1, index_2)
        line_number = patch_line + i + 2

        replace_command = 'sed -i "%si\%s" test_input.inp' % (line_number, line_to_write)

        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd test_input.inp' % (patch_line + 1)
            os.system(delete_command)
            break


def print_auto_lines_to_input_file(relevant_indices_and_fragments):

    with open('test_input.inp') as f:
        content = f.readlines()

    auto_line = content.index("INSERT NECESSARY AUTOGENERATIONS HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "if @NODE .eq. %s AUTO ANGL DIH" % (i + 1)
        line_number = auto_line + i + 2

        replace_command = 'sed -i "%si\%s" test_input.inp' % (line_number, line_to_write)

        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd test_input.inp' % (auto_line + 1)
            os.system(delete_command)
            break


def print_psf_lines_to_input_file(relevant_indices_and_fragments):

    with open('test_input.inp') as f:
        content = f.readlines()

    psf_line = content.index("INSERT NECESSARY PSF LINES HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "if @NODE .eq. %s print psf" % (i + 1)
        line_number = psf_line + i + 2

        replace_command = 'sed -i "%si\%s" test_input.inp' % (line_number, line_to_write)
        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd test_input.inp' % (psf_line + 1)
            os.system(delete_command)
            break


def print_shft_and_coup_lines_to_input_file(relevant_indices_and_fragments, shift_parameter=-25.0,
                                            coupling_parameter=105.0):

    with open('test_input.inp') as f:
        content = f.readlines()

    shift_line = content.index("         SHFT 0   0.0 \n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "         SHFT %s   %s" % ((i + 1), shift_parameter)
        line_number = shift_line + i + 2

        replace_command = 'sed -i "%si\%s" test_input.inp' % (line_number, line_to_write)
        os.system(replace_command)

        if i >= 14:
            break

    blank_line_command = 'sed -i "%si\ " test_input.inp' % (line_number + 1)
    os.system(blank_line_command)

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "         COUP 0 %s CONST %s" % ((i + 1), coupling_parameter)
        new_line_number = line_number + i + 2

        replace_command = 'sed -i "%si\%s" test_input.inp' % (new_line_number, line_to_write)
        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            break

def generate_probe_trajectory(squalane_simulation_dcd, squalane_simulation_crd, CN_simulation_dcd, CN_simulation_crd):

    return squalane_simulation_crd, CN_simulation_crd


def determine_reactive_CH_bonds(probe_trajectory_path, threshold_distance=5.0):

    #TODO: Add line to convert DCD to pdb (or function to use DCD)

    frames = read_pdb(probe_trajectory_path)
    indexes = get_atoms_within_threshold_distance(frames, threshold_distance)
    H_atoms = get_relevant_indices_H_atoms_only(frames, indexes)

    return H_atoms


def generate_production_trajectory(H_atoms, template_input_path):

    #TODO: Add functions that put probe crd and rst files with production input

    print_patch_lines_to_input_file(H_atoms, template_input_path)
    print_auto_lines_to_input_file(H_atoms)
    print_psf_lines_to_input_file(H_atoms)
    print_shft_and_coup_lines_to_input_file(H_atoms)


