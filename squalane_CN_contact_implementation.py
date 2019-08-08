from squalane_CN_contact import *

# If not already complete, generate the patch definitions file (squalane_patch_defs_all.txt) by uncommenting the two
# lines below
# squalane_mol2_file = "squalane_monomer-a.mol2"
# generate_patch_defs(range(31, 93), all_atoms, squalane_mol2_file)

# Determine relevant indexes of H atoms that come within a threshold distance of the CN C and the squalane fragments
# they are associated with
probe_path = '/Users/ec18006/OneDrive - University of Bristol/CHAMPS/Research_Topics/Squalane_project/Probe_and_Reactive_Dynamics_300719/input_file_generation/test_probe_traj.pdb'
distance_threshold = 5.0

probe_dfs = read_pdb(probe_path)
H_atoms = get_atoms_within_threshold_distance(probe_dfs, distance_threshold)
relevant_indices_and_fragments = get_relevant_indices_H_atoms_only(probe_dfs, H_atoms)

# Print lines associated with these H atoms that need to be added to CHARMM input file for reactive trajectory
print_input_file_lines(relevant_indices_and_fragments)

# Print lines associated with these H atoms directly to the CHARMM input file for a reactive trajectory

