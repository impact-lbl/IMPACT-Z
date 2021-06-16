!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! PhysConstclass: physical constants class in CONSTANTS module of 
!                 DATA STRUCTURE layer.
!
! MODULE  : ... PhysConstclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang, LBNL
!
! DESCRIPTION:
!> This class defines the physical constant parameters used
!> in the simulation.
! Comments:
!----------------------------------------------------------------
      module PhysConstclass
#if USE_MPI != 1
        use mpistub
#endif
        implicit none
      
        !physical parameters and constants ---------------------
        double precision :: Pi
        double precision :: Clight !< speed of light in vaccum
        double precision :: Scxl !< length scale
        double precision :: Rad2deg !< conversion factor from radian to degree
        double precision :: Epsilon0 !< permittivity of vacuum
        double precision :: Scfreq !< time scale
      contains
        subroutine construct_PhysConst(freq)
        double precision, intent(in) :: freq

        Clight = 299792458.0d0 
        Pi = 2*asin(1.0d0)
        Scxl = Clight/(2*pi*freq)
        Rad2deg = 180.0d0/Pi
        Epsilon0 = 8.854187817d-12
        Scfreq = freq

        end subroutine construct_PhysConst
 
      end module PhysConstclass
