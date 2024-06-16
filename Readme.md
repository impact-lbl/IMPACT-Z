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

To compile the ImpactZ code, one can follow the same procedure as described in 
the Impact-T github Readme file.

The ImpactZexeMac, ImpactZexeUbuntu, and ImpactZexeWin.exe are old executables.

Main contact: Ji Qiang (jqiang@lbl.gov), Lawrence Berkeley National Laboratory

# Compiling The Code

If you are new to CMake, [this short tutorial](https://hsf-training.github.io/hsf-training-cmake-webpage/) from the HEP Software foundation is the perfect place to get started with it.

If you just want to use CMake to build the project, jump into sections *1. Introduction*, *2. Building with CMake* and *9. Finding Packages*.

## Unix

### Single Processor Code:

```shell script
# inside the IMPACT-T src/ directory:
cmake -S . -B build
cmake --build build
# the executable in now in build/bin/

# this command needs sudo if you install into system paths:
cmake --build build --target install
```

