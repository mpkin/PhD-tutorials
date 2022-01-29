# PAMR

The code in this directory illustrates basic functionality of the PAMR (Parallel Adaptive Mesh Refinement) and AMRD (Adaptive Mesh Refinement Driver) libraries. This directory contains several examples:
1. `wave-local` -- solution of a basic wave equation in 1D, 2D, or 3D (based on Matt's `pamr-wave`, but modified to run locally)
2. `compactified` -- solution of a basic wave equation in axisymmetry using compactified coordinates (also uses `RNPL`)
3. `integration` -- an illustration of how to integrate quantities over arbitrary AMR levels (also uses `RNPL`)
4. `regrid_scripts.txt` -- a tutorial on using regrid scripts

In `pamr-installguide.md` and `pamr-installguide-CC.md` I provide instructions on how to install the PAMR libraries on Ubuntu and on the ComputeCanada systems.

Further reading:
```
J. Oliger and M. J. Berger. Adaptive mesh refinement for hyperbolic partial differential equations. J. Comput. Phys., 53:484–512, 1984.
F. Pretorius and M. W. Choptuik. Adaptive Mesh Refinement for Coupled Elliptic-Hyperbolic Systems. J. Comput. Phys., 218:246–274, 2006.
F. Pretorius. Numerical Simulations of Gravitational Collapse. PhD thesis, The University of British Columbia, 2002.
```
