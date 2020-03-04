# EVB Trajectory Automation for CN-Squalane System

This repository consists of a set of scripts that can be used to automate generation of Empirical Valence Bond (EVB) trajectories for the simulation of a CN radical being shot into a bath of squalane molecules using CHARMM. 

## Overview of Process

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**
2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 
3. **Sample snapshots from sampling trajectory and add CN radical with chosen coordinates and velocity to snapshots.** 
4. **Run these snapshots as NVE "probe" trajectories and save the outputs in the PDB format.**
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**
6. **Run "reactive" EVB trajectories.**

## Dependencies

Dependencies can be found in requirements.txt. To install all dependencies, use `pip install` in Terminal:

`pip install -r requirements.txt`

&nbsp;

## Details for Conducting Each Step

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**
   
Equilibration can be done in Gromacs, CHARMM, or any other MD program. All that is needed from this step is the coordinates generated at the end of the equilibration step (a single frame) in .crd format. 
   
**Output**
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) equilibration_frame.crd
   
&nbsp;

2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 

**Inputs** 
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) equilibration_frame.crd
- dynstrt.inp
- dyn.inp
- squalanemonomermmff.rtf
- dyn.csh
- lobos.csh
- cryst.str
- psfcrd.str
- rtfprm.str
- dyn.rea

This set of scripts basically "daisy-chains" several CHARMM trajectories together in order to generate several restart and coordinate files, which are usually only generated at the end of a trajectory. dynstrt.inp is the CHARMM input file for the first simulation and dyn.inp is the CHARMM input file for all subsequent simulations. dyn.csh is the main workhorse script and might need to be altered depending on the queueing system of the cluster being used. lobos.csh might also need to be changed (ever so slightly) depending on the queueing system. lobos.csh is the script that is submitted to the cluster queue and will continually submit CHARMM input files.

We still need to figure out how to limit the amount of times the dyn.csh script is submitted, because as of now, lobos.csh will run forever.

**Outputs**
- ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) dyn[#].res (many)
- ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) dyn[#].crd (many)

&nbsp;
 
3. **Extract snapshots from sampling trajectory and add CN radical with chosen coordinates and velocities.** 

**Inputs**
- ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) restart file (.res or .rst) path
- ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) coordinate file (.crd) path
- CN radical velocity
- timestep (of probe trajectory)

Uncomment [relevant lines] of generate_probe_and_reactive_trajectories.py to alter restart and coordinate files from sampling trajectory. 

**Outputs**
- ![#9dc010](https://placehold.it/15/9dc010/000000?text=+) with_cn_dyn[#].res (many)
- ![#106dc0](https://placehold.it/15/106dc0/000000?text=+) with_cn_dyn[#].crd (many)

&nbsp;
 
4. **Run these snapshots as NVE "probe" trajectories and save the output trajectory in the PDB format.**

**Inputs**
- ![#9dc010](https://placehold.it/15/9dc010/000000?text=+) restart file (with_cn_dyn[#].res)
- ![#106dc0](https://placehold.it/15/106dc0/000000?text=+) coordinate file (with_cn_dyn[#].crd) 
- CHARMM input (.inp) file
- squalane_cyanide_system_mmff.rtf

Submit NVE CHARMM input file using altered restart and coordinate files extracted from sampling trajectory. Save the output as a PDB file, for example by loading the original .crd file into VMD and loading the output .dcd file into it, and saving as a .pdb.

**Output**
- ![#cc99ff](https://placehold.it/15/cc99ff/000000?text=+) probe[#].pdb 

&nbsp;
 
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**

**Inputs**
- ![#cc99ff](https://placehold.it/15/cc99ff/000000?text=+) "Probe" trajectory PDB file (.pdb) path
- template CHARMM input file (.inp) path
- distance cutoff (in Angstroms)

Uncomment [relevant lines] of generate_probe_and_reactive_trajectories.py to determine C-H bonds to make reactive and generate a CHARMM input file that will allow these bonds to break. 

**Outputs**
- ![#ffc921](https://placehold.it/15/ffc921/000000?text=+) CHARMM input file (reactive[#].inp)

&nbsp;
 
6. **Run "reactive" EVB trajectories.**

**Inputs**
- ![#ffc921](https://placehold.it/15/ffc921/000000?text=+) CHARMM input file (reactive[#].inp)
- ![#9dc010](https://placehold.it/15/9dc010/000000?text=+) restart files *from probe trajectories* (.res or .rst)
- ![#106dc0](https://placehold.it/15/106dc0/000000?text=+) coordinate files *from probe trajectories* (.crd)
- squalane_patch_defs_all.txt

Run the reactive CHARMM input filesto generate your final EVB trajectories. 

**Outputs**
- CHARMM output files (.out)
- CHARMM DCD files (.dcd)

&nbsp;

