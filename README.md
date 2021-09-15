# ASOHF
Adaptive Spherical Overdensity Halo Finder (Planelles &amp; Quilis, 2010. A&amp;A, 519:A94)

## Short documentation (to be completed)

#### Code structure

The routines in this programme are distributed in the following files:

- ``asohf.f``: main program, contains the main execution workflow.
- ``grids.f``: build the base and AMR grids
- ``haloes_grids.f``: routines for halo treatment over the AMR grids
- ``haloes_particles.f``: routines for halo treatment using DM/stellar particles
- ``merger_tree.f``: build the merger tree of the halo catalogue between each pair of code outputs (reduced or complete)
- ``nr.f``: includes several routines found in 'Numerical Recipes in Fortran90', Press, Teukoslky et al.
- ``reader.f``: read the input data (may need to be adapted to read data from other simulation codes
- ``nomfile.f``: generates the filenames for I/O
