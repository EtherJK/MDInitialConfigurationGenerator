# MD Initial Configuration Generator

## Description

When you want to simulate the behavior of linear polymers with the freely jointed chain model, this code will help to generate a random initial configuration for software LAMMPS and DL_MESO. The format of input file follows the example file.

The current version must be used with the periodic boundary condition.

## Usage

Set the information of system and write it in `polymer_info.txt`, run `initial_generate_read_LAMMPS.py` or `initial_generate_read_DL_Meso.py`, and the `polymer.in` will be generated for LAMMPS, or `FIELD` be generated for DL_MESO.

For `initial_generate_read_LAMMPS.py`, the `polymer.in` includes all the polymers with random configurations and positions, and the interaction should be further given somewhere else in other input files. The beads may overlap with each other. So, a pre-equilibrium simulation should run in order to make sure that any two beads are not too close. In the pre-equilibrium simulation, the potential between beads should be soft.

For `initial_generate_read_DL_Meso.py`, the `FIELD` includes one random configuration for one species, and the interaction parameter has been included. The current version only applies to DPD potential and harmony bond potential.



