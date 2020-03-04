import openbabel
import numpy as np
import pandas as pd
import os
from scipy.spatial.distance import cdist
import MDAnalysis as md
import shutil
import subprocess

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
    """
    This will currently only work for a cyano radical (fragment name = CYA). If it is to be extended to other systems,
    this will need to be modified.
    :param traj_frames_df:
    :param threshold_distance:
    :return:
    """
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
    """
    This function takes the list of atom indexes that come within the threshold distance of the C of the CN radical
    and (1) determines which indexes are hydrogens and (2) what the local index and fragment index are for that H.
    :param traj_frames_df:
    :param atom_indexes:
    :return:
    """
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
            # NOTE: these fragment indices start at 1. In order to correctly visualize them in VMD, you have to subtract
            #  1 from the fragment index listed here.
            fragments.append(int(fragment_index))
            # fragments.append(int(fragment_index)-1)

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


def print_patch_lines_to_input_file(relevant_indices_and_fragments, template_file_path, number=1):

    with open(template_file_path) as template_file, open('reactive%s.inp' % number, 'w') as output_file:
        output_file.writelines(template_file)

    with open('reactive%s.inp' % number) as f:
        content = f.readlines()

    patch_line = content.index("INSERT NECESSARY PATCHES HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        index_1 = relevant_indices_and_fragments['local index'][i]
        index_2 = relevant_indices_and_fragments['fragment index'][i]

        line_to_write = "if @NODE .eq. %s PATCH RH%s A %s A 101 setup" % ((i + 1), index_1, index_2)
        line_number = patch_line + i + 2

        replace_command = 'sed -i "%si\%s" reactive%s.inp' % (line_number, line_to_write, number)

        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd reactive%s.inp' % ((patch_line + 1), number)
            os.system(delete_command)
            break

    #TODO: Have input be how many processors are available
    #TODO: Once threshold number of processors is hit, make NEW input file with the rest of the CH bonds


def print_auto_lines_to_input_file(relevant_indices_and_fragments, number=1):

    with open('reactive%s.inp' % number) as f:
        content = f.readlines()

    auto_line = content.index("INSERT NECESSARY AUTOGENERATIONS HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "if @NODE .eq. %s AUTO ANGL DIH" % (i + 1)
        line_number = auto_line + i + 2

        replace_command = 'sed -i "%si\%s" reactive%s.inp' % (line_number, line_to_write, number)

        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd reactive%s.inp ' % ((auto_line + 1), number)
            os.system(delete_command)
            break

        # TODO: Have input be how many processors are available
        # TODO: Once threshold number of processors is hit, make NEW input file with the rest of the CH bonds


def print_psf_lines_to_input_file(relevant_indices_and_fragments, number=1):

    with open('reactive%s.inp' % number) as f:
        content = f.readlines()

    psf_line = content.index("INSERT NECESSARY PSF LINES HERE\n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "if @NODE .eq. %s print psf" % (i + 1)
        line_number = psf_line + i + 2

        replace_command = 'sed -i "%si\%s" reactive%s.inp' % (line_number, line_to_write, number)
        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            delete_command = 'sed -i %sd reactive%s.inp' % ((psf_line + 1), number)
            os.system(delete_command)
            break

        # TODO: Have input be how many processors are available
        # TODO: Once threshold number of processors is hit, make NEW input file with the rest of the CH bonds


def print_shift_and_coupling_lines_to_input_file(relevant_indices_and_fragments, number=1, shift_parameter=-25.0,
                                            coupling_parameter=105.0):

    with open('reactive%s.inp' % number) as f:
        content = f.readlines()

    shift_line = content.index("         SHFT 0   0.0 \n")

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "         SHFT %s   %s" % ((i + 1), shift_parameter)
        line_number = shift_line + i + 2

        replace_command = 'sed -i "%si\%s" reactive%s.inp' % (line_number, line_to_write, number)
        os.system(replace_command)

        if i >= 14:
            break

    blank_line_command = 'sed -i "%si\ " reactive%s.inp' % ((line_number + 1), number)
    os.system(blank_line_command)

    for i in range(len(relevant_indices_and_fragments)):
        line_to_write = "         COUP 0 %s CONST %s" % ((i + 1), coupling_parameter)
        new_line_number = line_number + i + 2

        replace_command = 'sed -i "%si\%s" reactive%s.inp' % (new_line_number, line_to_write, number)
        os.system(replace_command)

        # Maximum of 16 nodes possible
        if i >= 14:
            break

    #TODO: Have input be how many processors are available
    #TODO: Once threshold number of processors is hit, make NEW input file with the rest of the CH bonds


def generate_probe_trajectory_starting_coordinates(squalane_coordinates_dcd, squalane_crd, squalane_snapshot_number=0):

    squalane_u = md.Universe(squalane_crd, squalane_coordinates_dcd)

    # Get coordinates from random snapshot from squalane trajectory
    squalane_coordinates = squalane_u.trajectory[squalane_snapshot_number]

    # Put CN normal to surface, C facing down
    CN_coordinates_shifted = [[0, 0, 75], [0, 0, 76]]
    # Stick the two snapshots together (CN X A from origin)
    squalane_and_CN_coordinates = np.append(squalane_coordinates, CN_coordinates_shifted, axis=0)

    return squalane_and_CN_coordinates


def generate_probe_trajectory_starting_velocities(squalane_velocities_dcd, squalane_crd, squalane_snapshot_number=0, cn_velocity=1800):

    squalane_u = md.Universe(squalane_crd, squalane_velocities_dcd, velocities=True)

    # Get coordinates from random snapshot from squalane trajectory
    squalane_velocities = squalane_u.trajectory[squalane_snapshot_number]

    # CHARMM units use AKMA system. AKMA unit of time is 4.888821E-14 seconds and velocity is in Angstroms/AKMA time
    # units. Need to convert cn_velocity in m/s to AKMA units.
    cn_velocity_akma = cn_velocity * (10 ** 10) * (4.888821 * 10 ** (-14))

    squalane_and_CN_velocities = squalane_velocities.append(cn_velocity_akma)

    return squalane_and_CN_velocities


def generate_probe_trajectory_crd_file(squalane_crd_file_path, cn_atom_coordinates=[[0, 0, 75], [0, 0, 76.172]]):
    # Copy and read in crd file
    file_copy_path = squalane_crd_file_path + '_with_cn'

    for i in range(len(cn_atom_coordinates)):
        for j in range(len(cn_atom_coordinates[i])):
            cn_atom_coordinates[i][j] = "{:.5f}".format(float(cn_atom_coordinates[i][j]))

    with open(squalane_crd_file_path, 'r') as f:
        i = 0
        line_no = 0
        for line in f:
            line_no += 1
            if not line.lstrip().startswith('*'):
                if i == 0:
                    sqa_natoms = int(line)
                    sqa_cn_natoms = sqa_natoms + 2
                    natom_replacement_command = "sed -e '0,/%s/ s/%s/%s/' %s > %s" % (sqa_natoms, sqa_natoms, sqa_cn_natoms,
                                                                                    "'" + squalane_crd_file_path + "'",
                                                                                    "'" + file_copy_path + "'")
                    os.system(natom_replacement_command)
                i += 1
                if i > 1 and line.split()[0] == str(sqa_natoms):
                    crd_line = line_no + 1
                    frag_num = int(line.split()[1]) + 1

    column_values_a = [sqa_natoms + 1, frag_num, 'CYA', 'C1', cn_atom_coordinates[0][0], cn_atom_coordinates[0][1],
                       cn_atom_coordinates[0][2], 'NODI', frag_num, '0.00000']
    column_values_b = [sqa_natoms + 2, frag_num, 'CYA', 'N2', cn_atom_coordinates[1][0], cn_atom_coordinates[1][1],
                       cn_atom_coordinates[1][2], 'NODI', frag_num, '0.00000']
    column_formats = ["{:>5}", "{:>4}", "{:<4}", "{:<4}", "{:>9}", "{:>9}", "{:>9}", "{:>4}", "{:<6}", "{:>7}"]

    for i in range(len(column_values_a)):
        column_values_a[i] = list(str(column_values_a[i]))
        column_values_a[i] = ("%s" % column_formats[i]).format("".join(column_values_a[i]))
        column_values_b[i] = list(str(column_values_b[i]))
        column_values_b[i] = ("%s" % column_formats[i]).format("".join(column_values_b[i]))

    column_values_a = ' '.join(column_values_a)
    column_values_b = ' '.join(column_values_b)

    # Insert blank line at end of file to be able to insert new lines there
    blank_line_command = 'sed -i -e "\$a\  " %s' % ("'" + file_copy_path + "'")
    os.system(blank_line_command)

    # Insert CN coordinates at point n
    insert_command_a = 'sed -i -e \'%si %s\' %s' % (crd_line, "\\" + column_values_a, "'" + file_copy_path + "'")
    insert_command_b = 'sed -i -e \'%si %s\' %s' % (crd_line + 1, "\\" + column_values_b, "'" + file_copy_path + "'")

    os.system(insert_command_a)
    os.system(insert_command_b)

    f.close()


def generate_probe_trajectory_rst_file(squalane_rst_file_path, cn_atom_coordinates=[[0, 0, 75], [0, 0, 76.172]],
                                       cn_velocity_ms=-1800, timestep=0.5*10**(-15)):

    # CHARMM units use AKMA system. AKMA unit of time is 4.888821E-14 seconds and velocity is in Angstroms/AKMA time
    # units. Need to convert cn_velocity in m/s to AKMA units.
    cn_velocity_akma = cn_velocity_ms * (10**10) * (4.888821 * 10**(-14))
    cn_velocity_As = cn_velocity_ms * (10**10)

    # Create name for copy of rst file to pipe to
    file_copy_path = squalane_rst_file_path + "_with_cn"

    #TODO: Fix this so "with_cn_" is at the beginning of the file name (not the path)

    with open(squalane_rst_file_path, 'r') as f:
        i = 1
        for line in f:
            i += 1
            if line.split(',')[0] == " !NATOM":
                sqa_natoms = int(next(f).split()[0])
                sqa_cn_natoms = sqa_natoms + 2
                natom_replacement_command = "sed -e '/!NATOM/!b' -e ':a' -e '"'s/%s/%s/;t trail'"' -e 'n;ba' -e " \
                                            "':trail' -e 'n;btrail' %s > %s" % (sqa_natoms, sqa_cn_natoms,
                                                                                "'" + squalane_rst_file_path + "'",
                                                                                "'" + file_copy_path + "'")
                os.system(natom_replacement_command)
            elif line.split(',')[0] == " !XOLD":
                xold_line = i
            elif line.split(',')[0] == " !VX":
                vx_line = i
            elif line.split(',')[0] == " !X":
                x_line = i

    cn_atom_coordinates_old = list(np.array(cn_atom_coordinates) - np.array([0, 0, cn_velocity_As*timestep]))
    for i in range(len(cn_atom_coordinates_old)):
        cn_atom_coordinates_old[i] = list(cn_atom_coordinates_old[i])

    cn_atom_velocities = [[0, 0, cn_velocity_akma], [0, 0, cn_velocity_akma]]

    cn_atom_delta_coordinates = list(np.array(cn_atom_coordinates) - cn_atom_coordinates_old)
    for i in range(len(cn_atom_delta_coordinates)):
        cn_atom_delta_coordinates[i] = list(cn_atom_delta_coordinates[i])

    for atom_value in [cn_atom_coordinates_old, cn_atom_velocities, cn_atom_delta_coordinates]:
        for i in range(len(atom_value)):
            new = [0, 0, 0]
            for j in range(len(atom_value[i])):
                new[j] = "{:.15e}".format(atom_value[i][j])
                new[j] = list(new[j])
                new[j] = ['D' if x == 'e' else x for x in new[j]]
                new[j] = "{:>22}".format("".join(new[j]))

            new = ''.join(new)
            atom_value[i] = new

    # Insert blank line at end of file to be able to insert new lines there
    blank_line_command = 'sed -i -e "\$a\  " %s' % ("'" + file_copy_path + "'")
    os.system(blank_line_command)

    # Insert CN coordinates at point n-1
    insert_command_1a = 'sed -i -e \'%si %s\' %s' % (xold_line + sqa_natoms + 1, "\\" + cn_atom_coordinates_old[0], "'" + file_copy_path + "'")
    insert_command_1b = 'sed -i -e \'%si %s\' %s' % (xold_line + sqa_natoms + 2, "\\" + cn_atom_coordinates_old[1], "'" + file_copy_path + "'")
    os.system(insert_command_1a)
    os.system(insert_command_1b)

    # Insert CN velocities at point n
    insert_command_2a = 'sed -i -e \'%si %s\' %s' % (vx_line + sqa_natoms + 3, "\\" + cn_atom_velocities[0], "'" + file_copy_path + "'")
    insert_command_2b = 'sed -i -e \'%si %s\' %s' % (vx_line + sqa_natoms + 4, "\\" + cn_atom_velocities[1], "'" + file_copy_path + "'")
    os.system(insert_command_2a)
    os.system(insert_command_2b)

    # Insert CN coordinates at point n
    insert_command_3a = 'sed -i -e \'%si %s\' %s' % (x_line + sqa_natoms + 5, "\\" + cn_atom_delta_coordinates[0], "'" + file_copy_path + "'")
    insert_command_3b = 'sed -i -e \'%si %s\' %s' % (x_line + sqa_natoms + 6, "\\" + cn_atom_delta_coordinates[1], "'" + file_copy_path + "'")

    os.system(insert_command_3a)
    os.system(insert_command_3b)

    f.close()


def determine_reactive_CH_bonds(probe_trajectory_dcd_path, probe_trajectory_crd_path, threshold_distance=5.0):
    #TODO: Make sure this is deprecated
    u = md.Universe(probe_trajectory_crd_path, probe_trajectory_dcd_path)

    with md.Writer("probe.pdb", u.atoms) as W:
        for ts in u.trajectory:
            W.write(u)

    frames = read_pdb('probe.pdb')
    indexes = get_atoms_within_threshold_distance(frames, threshold_distance)
    H_atoms = get_relevant_indices_H_atoms_only(frames, indexes)

    return H_atoms


def generate_production_trajectory(relevant_indices_and_fragments, template_input_path):

    #TODO: Add functions that put probe crd and rst files with production input
    #TODO: Have input be how many processors are available
    #TODO: Once threshold number of processors is hit, make NEW input file with the rest of the CH bonds

    print_patch_lines_to_input_file(relevant_indices_and_fragments, template_input_path)
    print_auto_lines_to_input_file(relevant_indices_and_fragments)
    print_psf_lines_to_input_file(relevant_indices_and_fragments)
    print_shift_and_coupling_lines_to_input_file(relevant_indices_and_fragments)


