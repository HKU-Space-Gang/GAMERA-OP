# GAMERA-OP (Orthogonal-Plus)

## Publication

The algorithms, numerical methods, and verification tests underlying GAMERA-OP are described in the following open-access paper:

> H. Luo, B. Zhang, J. Tian, J. Cai, J. Chen, E. Feng, Z. Zheng, S. Xi, and J. G. Lyon, “GAMERA-OP: A Three-dimensional Finite-volume Magnetohydrodynamic Solver for Orthogonal Curvilinear Geometries,” *The Astrophysical Journal Supplement Series*, **285**(1), 15 (2026).
>
> [Open-access article](https://doi.org/10.3847/1538-4365/ae7344)

If you use GAMERA-OP in your research, please cite this paper.

<details>
<summary>BibTeX</summary>

```bibtex
@article{Luo2026GAMERAOP,
  author    = {Luo, Hongyang and Zhang, Binzheng and Tian, Jiaxing and Cai, Jinshu and Chen, Junjie and Feng, Enhao and Zheng, Zhiqi and Xi, Sheng and Lyon, John G.},
  title     = {{GAMERA-OP}: A Three-dimensional Finite-volume Magnetohydrodynamic Solver for Orthogonal Curvilinear Geometries},
  journal   = {The Astrophysical Journal Supplement Series},
  year      = {2026},
  month     = jun,
  volume    = {285},
  number    = {1},
  pages     = {15},
  publisher = {The American Astronomical Society},
  doi       = {10.3847/1538-4365/ae7344},
  url       = {https://doi.org/10.3847/1538-4365/ae7344}
}
```

</details>

## Introduction

GAMERA-OP (Orthogonal Plus) is a three-dimensional finite-volume magnetohydrodynamics (MHD) solver designed for orthogonal curvilinear geometries. 
It serves as the successor to the LFM (Lyon–Fedder–Mobarry) and GAMERA codes, re-designed in the C language to emphasize improved numerics, 
modularity, high performance, and scalability for applications in space physics and astrophysical research.

This public release supports Cartesian and Spherical coordinates and includes a suite of standard benchmarks, 
such as the Orszag–Tang vortex, magnetic field loop advection, and blast wave tests. Its modular architecture allows 
researchers to easily implement new coordinates and physical configurations by leveraging the provided template structures.

## Non-orthogonal Earth production solver

The existing orthogonal solver above remains unchanged. The frozen Yin-Yang
global-magnetosphere production solver, including solar-wind driving,
magnetosphere-ionosphere coupling, electron precipitation, and conductance, is
maintained in [`nonorthogonal/`](nonorthogonal/README.md). New users should
start with its local two-rank smoke example before using the 256-rank Sugon
example.

Student entry points:

- [Production-model overview](nonorthogonal/README.md)
- [Local installation and first-run guide](nonorthogonal/docs/LOCAL_SETUP.md)
- [Grid sizes, MPI decompositions and intended uses](nonorthogonal/docs/GRID_OPTIONS.md)
- [Solar-wind input for idealized schedules and observed events](nonorthogonal/docs/SOLAR_WIND.md)
- [Local two-rank example and scripts](nonorthogonal/examples/earth_magnetosphere/local_smoke/)
- [Sugon setup and production-run guide](nonorthogonal/docs/SUGON.md)
- [Sugon 256-rank example and Slurm scripts](nonorthogonal/examples/earth_magnetosphere/sugon_128/)
- [Standard XY/XZ, FAC/potential and MI diagnostics](nonorthogonal/docs/DIAGNOSTICS.md)
- [Yin-Yang and schema-3 output data model](nonorthogonal/docs/OUTPUT_DATA_MODEL.md)
- [Accepted production provenance](nonorthogonal/docs/PRODUCTION_PROVENANCE.md)

## Requirements

GCC, MPI, OpenMP, HDF5, and CMake are required.

## Compilation

To compile the project, use the following commands:

```bash
mkdir build && cd build

# -DCOORD_TYPE=0 : Cartesian geometry
# -DCOORD_TYPE=2 : Spherical geometry

# -DPROBLEM=1 : Orszag–Tang test (Cartesian)
# -DPROBLEM=2 : Field loop test (Cartesian)
# -DPROBLEM=3 : Field loop test (Spherical)
# -DPROBLEM=4 : Blast wave test (Cartesian)
# -DPROBLEM=5 : Blast wave test (Spherical)

cmake -DCMAKE_BUILD_TYPE=Release -DCOORD_TYPE=2 -DPROBLEM=3 ..
make
```

## Execution
Create a script named run.sh with the following content, and execute it using ./run.sh.

```bash
#!/bin/bash

# Set the number of OpenMP threads
export OMP_NUM_THREADS=1

# Set the number of MPI processes
# Make sure the number of MPI processes equals the total number of subdomains defined in the config.yaml 
# (e.g., proc_dims_i x proc_dims_j x proc_dims_k = MPI processes)
mpirun -n 4 ./run_mhd config.yaml
```

A sample configuration file named `config.yaml` is provided in the `scripts/` directory.

## Data Processing

GAMERA-OP outputs data in a rank-specific and time-step-specific manner. Each MPI rank writes its own subdomain data to a separate HDF5 file. The naming convention follows the pattern:
```
mhd_XX-YY-ZZ_NNNNNN.h5
```

where:
- `XX-YY-ZZ` : MPI rank index in the (x, y, z) decomposition (e.g., `00-00-00` for rank (0,0,0))
- `NNNNNN`   : Output step number, zero-padded to 6 digits (e.g., `000000`, `000001`, `000002`, ...)

**Examples:**
- `mhd_00-00-00_000000.h5` – Rank (0,0,0), first output (step 0)
- `mhd_00-00-00_000001.h5` – Rank (0,0,0), second output (step 1)
- `mhd_01-00-00_000000.h5` – Rank (1,0,0), first output (step 0)

This per-rank, per-step output strategy facilitates parallel I/O and simplifies post-processing for distributed memory simulations. 
However, it also requires post-processing tools to combine or analyze the partitioned data.

A MATLAB script is provided in the `scripts/` directory to handle the rank-split output files. 
The script reads all HDF5 files for a given output step, combines them into a full 3D field, and saves the results in a more convenient hdf5 format for analysis and visualization.
