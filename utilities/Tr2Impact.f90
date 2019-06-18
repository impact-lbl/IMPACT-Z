!----------------------------------------------------------------
! (c) Copyright, 2011 by the Regents of the University of California.
! Version: 1.0
! Author: Ji Qiang, LBNL
!---------------
! convert from the Trace3d twiss parameter to the IMPACT input distribution parameters.
! Here, the TRACE3d units are, alpha, beta (m/rad), emittance (mm-mrad, rms, unnormalized)
!                              alpha_z,beta_z (deg/keV), emittance (deg-keV,rms)
!
      program Tr2Impact
      implicit none
!-----File Tr2Impact.f90:  Generates IMPACT input distribution from
!     "beam.dat" Trace3D twiss parameter file
 
      double precision alpha, beta, eps  ! C-S parameters from beam.dat.
      double precision sig1, sig2, rho   ! IMPACT parameters for sim.dat.
      double precision mass,kine,bgamma,bbeta,pi,clight,rad2deg,xl,freq
      integer i
 
      open(2,file='sim.dat',status='unknown')
      open(3,file='beam.dat',status='old')
 
      pi=4.d0*atan(1.d0)
      clight = 2.99792458e8
      rad2deg = 180.0/pi
      read(3,*)freq,mass,kine
      bgamma = 1 + kine/mass
      bbeta = sqrt(1.0-1.0/bgamma/bgamma)
      xl=clight/(2.d0*pi*freq)     ! Length scale
      do i=1,3
        read(3,*) alpha, beta, eps
        !-----Change units
        if (i.ne.3) then
           eps=1.0d-6*eps  ! unnormalized emittance
           sig1=sqrt(beta*eps/(1.d0+alpha**2))
           !gammabeta is needed because <x'^2>=emittance*gamma instead of <(px/mc)^2>
           sig2=sqrt(eps/beta)*bgamma*bbeta 
           !sig2=sqrt(eps/beta) !correct for x' or px/p0 as independent variable
           rho=alpha/sqrt(1.d0+alpha**2)
           write(2,'(3e15.7,4f7.3)') sig1/xl, sig2, rho, 1.0, 1.0, 0.0, 0.0
        else
           beta=1000.d0*beta
           eps=1.0d-3*eps
           sig1=sqrt(beta*eps/(1.d0+alpha**2))/rad2deg
           sig2=sqrt(eps/beta)/mass*1.0e6
           !sig2=sqrt(eps/beta)/mass  !//old wrong one
           rho=-alpha/sqrt(1.d0+alpha**2)
           write(2,'(3e15.7,4f7.3)') sig1, sig2, rho, 1.0, 1.0, 0.0, 0.0
        endif
      enddo
 
      close(2)
      close(3)
 
      end program Tr2Impact
!-------------------------------------

