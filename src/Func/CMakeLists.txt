
#gfortran -c -O3 Func/Timer.f90
#gfortran -c -O3 Func/Transpose.f90
#gfortran -c -O3 Func/Fldmger.f90
#gfortran -c -O3 Func/Ptclmger.f90
#gfortran -c -O3 Func/FFT.f90

add_library(impact_func Timer.f90 Transpose.f90 Fldmger.f90 Ptclmger.f90 FFT.f90)

target_link_libraries(impact_func impact_mpi)
