
#gfortran -c -O3 DataStruct/NumConst.f90
#gfortran -c -O3 DataStruct/PhysConst.f90
#gfortran -c -O3 DataStruct/Pgrid.f90
#gfortran -c -O3 DataStruct/Data.f90

add_library(impact_datastruct NumConst.f90 PhysConst.f90 Pgrid.f90 Data.f90)

target_link_libraries(impact_datastruct impact_mpi)
