## EVB Trajectory Automation for CN-Squalane System

This repository consists of a set of scripts that can be used to automate generation of Empirical Valence Bond (EVB) trajectories for the simulation of a CN radical being shot into a bath of squalane molecules using CHARMM. 

The major steps are as follows:

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**
2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 
3. **Sample snapshots from sampling trajectory and add CN radical with chosen coordinates and velocity to snapshots.** 
4. **Run these snapshots as NVE "probe" trajectories and save the outputs in the PDB format.**
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**
6. **Run "reactive" EVB trajectories.**

&nbsp;

And here are the details for conducting each step:

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**

&nbsp;

2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 

   **Inputs:** 
- dynstrt.inp
- dyn.inp
- [input].crd
- squalanemonomermmff.rtf
- dyn.csh
- lobos.csh

   **Outputs:**
- dyn[#].res (many)
- dyn[#].crd (many)

&nbsp;
 
3. **Extract snapshots from sampling trajectory and add CN radical with chosen coordinates and velocities.** 

   **Inputs:**
- restart file (.res or .rst) path
- coordinate file (.crd) path
- CN radical velocity
- timestep (of probe trajectory)

   **Outputs:**
- with_cn_dyn[#].res (many)
- with_cn_dyn[#].crd (many)

&nbsp;
 
4. **Run these snapshots as NVE "probe" trajectories and save the outputs in the PDB format.**

   **Inputs:**
- restart file (.res or .rst)
- coordinate file (.crd) 
- CHARMM input (.inp) file
- squalane_cyanide_system_mmff.rtf

   **Outputs:**
- probe[#].pdb 

&nbsp;
 
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**

   **Inputs:**
- "Probe" trajectory PDB file (.pdb) path
- template CHARMM input file (.inp) path
- distance cutoff (in Angstroms)

   **Outputs:**
- reactive[#].inp 

&nbsp;
 
6. **Run "reactive" EVB trajectories.**

   **Inputs:**
- CHARMM input file (reactive[#].inp)
- restart files (.res or .rst)
- coordinate files (.crd)
- squalane_patch_defs_all.txt

&nbsp;

