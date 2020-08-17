
#gfortran -c -O3 Appl/BPM.f90
#gfortran -c -O3 Appl/CCL.f90
#gfortran -c -O3 Appl/CCDTL.f90
#gfortran -c -O3 Appl/DTL.f90
#gfortran -c -O3 Appl/SC.f90
#gfortran -c -O3 Appl/DriftTube.f90
#gfortran -c -O3 Appl/Quadrupole.f90
#gfortran -c -O3 Appl/ConstFoc.f90
#gfortran -c -O3 Appl/SolRF.f90
#gfortran -c -O3 Appl/Sol.f90
#gfortran -c -O3 Appl/Dipole.f90
#gfortran -c -O3 Appl/Multipole.f90
#gfortran -c -O3 Appl/EMfld.f90
#gfortran -c -O3 Appl/TWS.f90
#gfortran -c -O3 Appl/BeamLineElem.f90
#gfortran -c -O3 Appl/CompDom.f90
#gfortran -c -O3 Appl/Field.f90
#gfortran -c -O3 Appl/BeamBunch.f90
#gfortran -c -O3 Appl/Distribution.f90

add_library(impact_appl BPM.f90 CCL.f90 CCDTL.f90 DTL.f90 SC.f90 DriftTube.f90 Quadrupole.f90 ConstFoc.f90 SolRF.f90 Sol.f90 Dipole.f90 Multipole.f90 EMfld.f90 TWS.f90 BeamLineElem.f90 CompDom.f90 Field.f90 BeamBunch.f90 Distribution.f90)

target_link_libraries(impact_appl impact_func)
