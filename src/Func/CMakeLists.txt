add_library(impact_func Timer.f90 Transpose.f90 Fldmger.f90 Ptclmger.f90 FFT.f90)

target_link_libraries(impact_func impact_datastruct impact_mpi)
if(USE_MPI)
    target_link_libraries(impact_func ${MPI_Fortran_LIBRARIES})
endif()
