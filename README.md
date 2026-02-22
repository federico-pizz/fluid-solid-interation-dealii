# FSI Solver

A monolithic fluid-structure interaction solver built on [deal.II](https://www.dealii.org/), coupling Stokes flow with linear elasticity on a shared domain.

## What it does

Solves the steady FSI problem on $\Omega = \Omega_f \cup \Omega_s$: Stokes equations in the fluid, linear elasticity in the solid, stress continuity at the interface $\Sigma$. Works in 2D and 3D.

## Project structure

```
src/
├── FSIParallel.hpp / .cpp   # parallel solver (MPI + OpenMP)
├── FSI.hpp / .cpp           # sequential solver
├── simulate_parallel.cpp    # parallel entry point
└── simulate_serial.cpp      # serial entry point
```

## Key features

- Taylor-Hood elements (Q2/Q1) for the fluid, satisfying the inf-sup condition
- Block-triangular preconditioner with AMG subsolves (GMRES outer solver)
- Adaptive mesh refinement via Kelly error estimator
- Hybrid MPI + OpenMP parallelism — tested up to 16 threads on HPC hardware

## Dependencies

- deal.II ≥ 9.3
- Trilinos (for AMG and parallel linear algebra)
- Boost ≥ 1.72 (`filesystem`, `iostreams`, `serialization`)
- MPI + OpenMP

## Build & run

```bash
cmake -DDEAL_II_DIR=/path/to/dealii . && make
mpirun -n $(nproc) ./fsi_test_parallel
```

> To switch between 2D and 3D, set `static constexpr unsigned int dim` in `FSIParallel.hpp` before building.

Results are exported in VTK format and can be visualized with ParaView.
