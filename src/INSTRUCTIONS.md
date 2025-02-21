# Using Anaconda

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

## Using conda to compile IMPACT-Z

`conda-forge` has all of the necessary compilers and dependencies to build IMPACT-Z from source.

Create a build environment like so:

```bash
conda create -n impactz-build -c conda-forge compilers cmake openmpi
conda activate impactz-build
```

Then to build the non-parallel version:

```bash
cmake -S src/ -B build-single
cmake --build build-single -j 4
ls build-single/ImpactZexe
```

And to build the MPI-parallelized version with OpenMPI:

```bash
cmake -S src/ -B build-mpi -DUSE_MPI=ON
cmake --build build-mpi -j 4
ls build-mpi/ImpactZexe-mpi
```

## Unix

### Single Processor Code:

```shell script
# inside the IMPACT-Z src/ directory:
cmake -S . -B build
cmake --build build
# the executable in now in build/bin/

# this command needs sudo if you install into system paths:
cmake --build build --target install
```
If you like to install IMPACT-Z into another directory than the default, pass to the `cmake -S . -B build` line the additional argument `-DCMAKE_INSTALL_PREFIX=/your/custom/install/path`.

### Multi Processor Code:

```shell script
# inside the IMPACT-Z src/ directory:
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
cmake --build build --target install
```

## Windows

For Windows it will be necessary to use `NMake` to read and execute the generated makefiles.

`NMake` is command-line tool included with Visual Studio that builds projects based on commands that are contained in a description file.

More information on `NMake` can be found on the [NMAKE Reference](https://docs.microsoft.com/en-us/cpp/build/reference/nmake-reference?view=msvc-160) website.

### Single Processor Code:

```shell script
cmake -S . -B build -G "NMake Makefiles"
cmake --build build
cmake --build build --target install
cmake --install
```

### Multi Processor Code:

**Not Tested**


## Testing

```shell script
cd examples
cd Example1
# Execute non-MPI
ImpactZexe

# Execute MPI
mpirun -n <cores> ImpactZexe-mpi
```

# Using at NERSC

If you decide to use IMPACT-Z at NERSC you can do so via Anaconda or building from source.
The instructions for Anaconda are the same as above.

## Compiling the code

### For Haswell
```bash
module load openmpi # if using OpenMPI otherwise skip for MPICH
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

### For KNL
```bash
module swap craype-haswell craype-mic-knl
module load openmpi # if using OpenMPI otherwise skip for MPICH
cmake -S . -B build -DUSE_MPI=ON
cmake --build build
# find the executable in build/bin/
```

## Running at NERSC

There is one caveat with the conda-forge installed MPI version of the code.
Instead of running with `srun -n <cores> ` you must use `mpirun -n <cores>`.
