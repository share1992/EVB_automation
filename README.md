# EVB Trajectory Automation for CN-Squalane System

This repository consists of a set of scripts that can be used to automate generation of Empirical Valence Bond (EVB) trajectories for the simulation of a CN radical being shot into a bath of squalane molecules using CHARMM. 

## Overview of Process

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**
2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 
3. **Sample snapshots from sampling trajectory and add CN radical with chosen coordinates and velocity to snapshots.** 
4. **Run these snapshots as NVE "probe" trajectories and save the outputs in the PDB format.**
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**
6. **Run "reactive" EVB trajectories.**

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

   **Outputs**
- probe[#].pdb 

&nbsp;
 
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**

   **Inputs**
- "Probe" trajectory PDB file (.pdb) path
- template CHARMM input file (.inp) path
- distance cutoff (in Angstroms)

   **Outputs**
- CHARMM input file (reactive[#].inp)

&nbsp;
 
6. **Run "reactive" EVB trajectories.**

   **Inputs**
- CHARMM input file (reactive[#].inp)
- ![#9dc010](https://placehold.it/15/9dc010/000000?text=+) restart files *from probe trajectories* (.res or .rst)
- ![#106dc0](https://placehold.it/15/106dc0/000000?text=+) coordinate files *from probe trajectories* (.crd)
- squalane_patch_defs_all.txt

   **Outputs**
- CHARMM output files (.out)
- CHARMM DCD files (.dcd)

&nbsp;

