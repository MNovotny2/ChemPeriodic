# ChemPeriodic

A simple tool for generating symmetry-inequivalent surface structures. 

## Required

- **Fortran Compiler**: GFortran 7.0+ or Intel ifort 2021+
- **MPI Library**: OpenMPI 3.0+ or MPICH 3.3+
- **Build Tool**: Make

## Installation

### Ubuntu/Debian
```bash
sudo apt-get install gfortran libopenmpi-dev make
```
Compile using the makefile

```bash
make
```

## Usage Example

```bash
mpirun -np 4 ChemPeriodic
```
## Implemented symmetry operations
Currently, the code can perform only symmetry operations on patterns in a square type lattice: translation along the lattice vector and rotation of 90&deg. 
