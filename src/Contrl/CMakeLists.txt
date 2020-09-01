add_library(impact Input.f90 Output.f90 AccSimulator.f90)
target_link_libraries(impact impact_datastruct impact_appl)

add_executable(ImpactZ main.f90)
target_link_libraries(ImpactZ impact)

if(USE_MPI)
    target_link_libraries(impact ${MPI_Fortran_LIBRARIES})
    target_link_libraries(ImpactZ ${MPI_Fortran_LIBRARIES})
endif()
