
#gfortran -c -O3 Contrl/Input.f90
#gfortran -c -O3 Contrl/Output.f90
#gfortran -c -O3 Contrl/AccSimulator.f90
#gfortran -c -O3 Contrl/main.f90

add_library(impact Input.f90 Output.f90 AccSimulator.f90)
target_link_libraries(impact impact_datastruct impact_appl)

add_executable(ImpactZ main.f90)
target_link_libraries(ImpactZ impact)
