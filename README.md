# SUNSET code

The **S**calable, **U**nstructured **N**ode-**SET** code, for high-fidelity mesh-free simulations of viscoelastic flows.

Developed by Dr Jack King, University of Manchester.

## Overview and features

- Numerically solves the isothermally compressible Navier-Stokes equations.
- Along with a conformation tensor evolution equation
- Can model Oldroyd B and sPTT fluids (for now).
- Domain discretised with unstructured node-set.
- Spatial discretisation between 4th and 10th order using **LABFM**.
- Temporal discretisation 3rd order explicit Runge-Kutta.
- Characteristic based boundary conditions
   + Walls (can be curved)
   + Inflow (non-reflecting or hard) (TBC)
   + Outflow (non-reflecting) (TBC)
   + Symmetry
   + Periodic
- Parallelised with OpenMP + MPI.

## Documentation

Documentation in docs folder is from the combustion variant (sunset-flames). I haven't yet produced a version of the documentation for this viscoelastic variant. For now, see ref 3 below, and the combustion documentation. See comments in the Makefile for compiler flag options.

- make newt=1    - Newtonian simulations
- make ceform=X  - Viscoelastic simulations (X=0 Direct integration, X=1 Cholesky, X=2 Log-Cholesky, X=3 Log-conformation)

## References

- LABFM fundamentals: King, Lind & Nasar (2020) JCP 415:109549 https://doi.org/10.1016/j.jcp.2020.109549
- LABFM fundamentals + BCs: King & Lind (2022) JCP 449:110760 https://doi.org/10.1016/j.jcp.2021.110760
- LABFM for viscoelastic flows: King & Lind (2024) JNNFM 330:105278 https://doi.org/10.1016/j.jnnfm.2024.105278


