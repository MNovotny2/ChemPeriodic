# ChemPeriodic

A simple tool for generating symmetry-inequivalent surface structures.

## Requirements

- **Fortran Compiler**: GFortran 7.0+ or Intel ifort 2021+
- **MPI Library**: OpenMPI 3.0+ or MPICH 3.3+
- **Build Tool**: Make

## Installation

### Ubuntu/Debian

```bash
sudo apt-get install gfortran libopenmpi-dev make
```

Compile using the Makefile:

```bash
make
```

## Usage

### Example run command

```bash
mpirun -np 4 ChemPeriodic < example.inp
```

### Process

The program prompts for the number of substituents (**n**) and substitution sites (**k**), from which it generates and encodes all possible variations (**n<sup>k</sup>**). These are stored in `variations.dat`. Variations with the same substituent frequencies are grouped and stored in `frequency.dat`.

The program then filters out all symmetry-equivalent variants and outputs only the unique ones, together with the number of symmetry-equivalent variations represented by each pattern. These are stored in `variations_unique.dat`.

## Implemented Symmetry Operations

Currently, the code supports symmetry operations only for patterns on a square lattice:

- Translation along the lattice vectors
- Rotation by 90°
