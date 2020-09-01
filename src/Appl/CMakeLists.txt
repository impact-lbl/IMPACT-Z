add_library(impact_appl BPM.f90 CCL.f90 CCDTL.f90 DTL.f90 SC.f90 DriftTube.f90 Quadrupole.f90 ConstFoc.f90 SolRF.f90 Sol.f90 Dipole.f90 Multipole.f90 EMfld.f90 TWS.f90 BeamLineElem.f90 CompDom.f90 Field.f90 BeamBunch.f90 Distribution.f90)

target_link_libraries(impact_appl impact_func)
if(USE_MPI)
    target_link_libraries(impact_appl ${MPI_Fortran_LIBRARIES})
endif()
