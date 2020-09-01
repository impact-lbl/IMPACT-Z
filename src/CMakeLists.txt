cmake_minimum_required(VERSION 3.1)

project(ImpactZ Fortran)

# Location of fortran modules:
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/include/fortran/impact CACHE PATH "Single Directory for all fortran modules."
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -fcheck=bounds ${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_FLAGS_RELEASE   "-O3 ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-Og -g -traceback -check bounds ${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_FLAGS_RELEASE   "-O3 ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

add_library(impact_mpi mpistub.f90)

option(USE_MPI "Activate MPI build" OFF)

if(USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    link_directories(${MPI_Fortran_LIBRARIES})
    target_link_libraries(impact_mpi ${MPI_Fortran_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(DataStruct)
add_subdirectory(Func)
add_subdirectory(Appl)
add_subdirectory(Contrl)

enable_testing()
add_test(NAME example1 COMMAND ImpactZ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example1/)
add_test(NAME example2 COMMAND ImpactZ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example2/)
add_test(NAME example3 COMMAND ImpactZ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Example3/)
