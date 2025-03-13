# set minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Uncomment for VERBOSE Makefiles. Useful for debugging.
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Uncomment for Debug Build
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# project name and language
project(ImpactZ LANGUAGES Fortran)

option(USE_MPI "Activate MPI build" OFF)


list(APPEND _sources
    DataStruct/NumConst.f90
    DataStruct/Pgrid.f90
    DataStruct/Data.f90
    DataStruct/PhysConst.f90
    Func/FFT.f90
    Func/Timer.f90
    Func/Fldmger.f90
    Func/Transpose.f90
    Func/Ptclmger.f90
    Appl/BPM.f90
    Appl/CCDTL.f90
    Appl/ConstFoc.f90
    Appl/Distribution.f90
    Appl/Field.f90
    Appl/SC.f90
    Appl/TWS.f90
    Appl/Wiggler.f90
    Appl/BeamBunch.f90
    Appl/CCL.f90
    Appl/DTL.f90
    Appl/DriftTube.f90
    Appl/Multipole.f90
    Appl/Sol.f90
    Appl/BeamLineElem.f90
    Appl/CompDom.f90
    Appl/Dipole.f90
    Appl/EMfld.f90
    Appl/Quadrupole.f90
    Appl/SolRF.f90
    Contrl/Output.f90
    Contrl/AccSimulator.f90
    Contrl/Input.f90
    Contrl/main.f90
)

set(OUTPUT_NAME "ImpactZexe")
if(USE_MPI)
    set(OUTPUT_NAME "ImpactZexe-mpi")
endif()
add_executable(${OUTPUT_NAME})

if(USE_MPI)
    find_package(MPI REQUIRED)
    target_compile_definitions(${OUTPUT_NAME} PUBLIC USE_MPI)
    target_link_libraries(${OUTPUT_NAME} PUBLIC MPI::MPI_Fortran)
    include_directories(${MPI_Fortran_INCLUDE_PATH})

    # Allow type mismatch for GCC 10
    # (Adapted from from: https://github.com/Unidata/netcdf-fortran/blob/main/CMakeLists.txt)
    # Check to see if compiler supports
    # -fallow-argument-mismatch-flag introduced in gcc 10.
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag("-fallow-argument-mismatch" COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)
    if(COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)
        target_compile_options(${OUTPUT_NAME} PRIVATE "-fallow-argument-mismatch")
    endif(COMPILER_HAS_ALLOW_ARGUMENT_MISMATCH)
else()
    # Add mpistub to source list
    list(APPEND _sources mpistub.f90)

    # Create mpif.h from mpif_single.h
    configure_file(mpif_single.h ${ImpactZ_BINARY_DIR}/include/mpif.h COPYONLY)
    # Add include dir above into search path
    target_include_directories(${OUTPUT_NAME} PRIVATE $<BUILD_INTERFACE:${ImpactZ_BINARY_DIR}/include>)
endif()

# Associate source files with target
target_sources(${OUTPUT_NAME} PRIVATE ${_sources})

# Enables preprocessing
target_compile_options(${OUTPUT_NAME} PRIVATE -cpp)

enable_testing()
add_test(NAME example1 COMMAND ${OUTPUT_NAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example1/)
add_test(NAME example2 COMMAND ${OUTPUT_NAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example2/)
add_test(NAME example3 COMMAND ${OUTPUT_NAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example3/)

install(TARGETS ${OUTPUT_NAME} RUNTIME DESTINATION bin)
