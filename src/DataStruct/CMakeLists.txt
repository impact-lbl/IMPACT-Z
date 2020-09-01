add_library(impact_datastruct NumConst.f90 PhysConst.f90 Pgrid.f90 Data.f90)

target_link_libraries(impact_datastruct impact_mpi)
if(USE_MPI)
    target_link_libraries(impact_datastruct ${MPI_Fortran_LIBRARIES})
endif()

