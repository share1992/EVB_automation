## EVB Trajectory Automation for CN-Squalane System

This repository consists of a set of scripts that can be used to automate generation of Empirical Valence Bond (EVB) trajectories for the simulation of a CN radical being shot into a bath of squalane molecules using CHARMM. The major steps are as follows, with the details of how to conduct each step outlined below:

1. **Equilibrate slab of squalane molecules. _This step should be done prior to using this repository._**
2. **Run an NVT sampling trajectory of the equilibrated slab of squalane molecules.** 
3. **Sample snapshots from sampling trajectory and add CN radical with chosen coordinates and velocity to snapshots.** 
4. **Run these snapshots as NVE "probe" trajectories and save the outputs in the PDB format.**
5. **Determine C-H bonds that come close to CN radical and generate relevant "reactive" trajectory input files.**
6. **Run "reactive" EVB trajectories.**
