import os
from squalane_cn import *

dirname = os.path.dirname(__file__)
probe_dirname = os.path.join(dirname, 'test_input_files_probe')
reactive_dirname = os.path.join(dirname, 'test_input_files_reactive')
template_dirname = os.path.join(dirname, 'template_and_patch_files')

# If not already complete, generate the patch definitions file (squalane_patch_defs_all.txt) by uncommenting the two
# lines below
# squalane_mol2_file = "squalane_monomer-a.mol2"
# generate_patch_defs(range(31, 93), all_atoms, squalane_mol2_file)

# Generate probe trajectory
# Path to restart file being altered (a copy will be made)
sqa_rst = os.path.join(probe_dirname, 'dyn100.res')
# sqa_rst = os.path.join(probe_dirname, 'dyn101.res')
# sqa_rst = os.path.join(probe_dirname, 'dyn102.res')

# Path to crd file being altered (a copy will be made)
sqa_crd = os.path.join(probe_dirname, 'dyn100.crd')
# sqa_crd = os.path.join(probe_dirname, 'dyn101.crd')
# sqa_crd = os.path.join(probe_dirname, 'dyn102.crd')

# CN velocity in m/s in the z-direction (default=-1800, or 1800 m/s in the -z-direction)
cn_velocity_ms = -1800
# Be sure to include correct time step of trajectory used to generate rst file (used to calculate "X/Y/ZOLD" of CN) in
# units of s (default=0.5*10**(-15), or 0.5 fs)
timestep = 0.5*10**(-15)

generate_probe_trajectory_crd_file(sqa_crd)
generate_probe_trajectory_rst_file(sqa_rst, cn_velocity_ms=cn_velocity_ms, timestep=timestep)

# Determine relevant indexes of H atoms that come within a threshold distance of the CN C and the squalane fragments
# they are associated with
# probe_path = os.path.join(reactive_dirname, 'dyn100_with_cn_probe_test_1.pdb')
# template_input_path = os.path.join(template_dirname, 'template_input.inp')

# distance_threshold = 5.0

# probe_dfs = read_pdb(probe_path)
# all_atoms_within_threshold = get_atoms_within_threshold_distance(probe_dfs, distance_threshold)
# relevant_indices_and_fragments = get_relevant_indices_H_atoms_only(probe_dfs, all_atoms_within_threshold)

# Print lines associated with these H atoms directly to the CHARMM input file for a reactive trajectory
# generate_production_trajectory(relevant_indices_and_fragments, template_input_path)

