cmake_minimum_required(VERSION 3.1)

project(ImpactZ Fortran)

# Location of fortran modules:
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/include/fortran/impact CACHE PATH "Single Directory for all fortran modules."
)
#gfortran -c -O3 mpistub.f90
add_library(impact_mpi mpistub.f90)

add_subdirectory(DataStruct)
add_subdirectory(Func)
add_subdirectory(Appl)
add_subdirectory(Contrl)

