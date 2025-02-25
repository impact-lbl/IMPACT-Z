IMPACT-Z is a parallel+serial particle-in-cell code whose primary purpose is to
model the dynamics of multiple charged particle beams in linear and ring acceler
ators. The code uses longitudinal position (z) as independent variable and 
includes the effects of externally applied fields from magnets and accelerating 
cavities as well as the effect of self-fields (space charge fields). Mathematically,
the code solves the Vlasov/Poisson equations using a particle-based technique. 
The code, which is written in Fortran90 with MPI, runs on both single-processor 
and multi-processor systems. It has been applied to studies of halo formation and
coupling resonance in high intensity beams, microbunching instability in high 
brightness electron linac, beam dynamics in SNS linac, JARPC linac, RIA driver 
linac, CERN superconducting linac, LEDA halo experiment, Proton Synchrotron at 
CERN, etc.


The ImpactZexeMac, ImpactZexeUbuntu, and ImpactZexeWin.exe are old executables.

Main contact: Ji Qiang (jqiang@lbl.gov), Lawrence Berkeley National Laboratory

Citation: 
J. Qiang, R. Ryne, S. Habib, V. Decyk, "An Object-Oriented Parallel Particle-In-Cell Code for Beam Dynamics Simulation in Linear Accelerators," 
J. Comp. Phys. vol. 163, 434, (2000).

# Installation using Anaconda

Information about Anaconda, including install instructions, can be found on the [Conda Reference](https://docs.conda.io/projects/conda/en/latest/) website.

IMPACT-Z is available through conda-forge and can be installed via:
```bash
conda create -n impact
source activate impact # or conda activate impact
# For non-MPI
conda install -c conda-forge impact-z

# For OpenMPI
conda install -c conda-forge impact-z=*=mpi_openmpi*

# For MPICH
conda install -c conda-forge impact-z=*=mpi_mpich*
```
After these steps, the IMPACT-Z executable `ImpactZexe` or `ImpactZexe-mpi`, respectively, will be in your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)) environment variable and is thus ready to use like any regular command-line command.

# Compiling The Code

If you are new to CMake, [this short tutorial](https://hsf-training.github.io/hsf-training-cmake-webpage/) from the HEP Software foundation is the perfect place to get started with it.

If you just want to use CMake to build the project, jump into sections *1. Introduction*, *2. Building with CMake* and *9. Finding Packages*.

## Unix

### Single Processor Code:

```shell script
# inside the IMPACT-Z src/ directory:
cmake -S . -B build
cmake --build build
# the executable in now in build/

# this command needs sudo if you install into system paths:
cmake --build build --target install
```

### Multi Processor Code:

```shell script
# inside the IMPACT-Z src/ directory:
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
cmake --build build --target install
