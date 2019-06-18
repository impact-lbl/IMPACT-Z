!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Timerclass: Clock for individual subroutine class in Utility module
!             of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang
! Description: This module is to record time spent in different subroutines.
! Comments:
!----------------------------------------------------------------
      module Timerclass
        use mpistub
        double precision :: t_integ
        double precision :: t_map1
        double precision :: t_map2
        double precision :: t_prnt
        double precision :: t_dist
        double precision :: t_kvdist
        double precision :: t_gaussn
        double precision :: t_normdv
        double precision :: t_correct
        double precision :: t_spch2d
        double precision :: t_boundgeom
        double precision :: t_greenf
        double precision :: t_rhofas
        double precision :: t_ntrslo
        double precision :: t_fft2dhpf
        double precision :: t_mfft_local1
        double precision :: t_field
        double precision :: t_force
        double precision :: t_charge
        double precision :: t_beamln
        double precision :: t_diag
        double precision :: t_ptsmv
        double precision :: t_init
        double precision :: t_transp
        double precision :: t_enlarge
        double precision :: t_shrink
        double precision :: t_guardsum
        double precision :: t_boundint
        double precision :: t_guardexch
      contains
        subroutine construct_Timer(values)
        implicit none
        double precision, intent(in) :: values
  
        t_integ = values
        t_map1 = values
        t_map2 = values
        t_diag = values
        t_dist = values
        t_kvdist = values
        t_gaussn = values
        t_normdv = values
        t_correct = values
        t_spch2d = values
        t_boundgeom = values
        t_greenf = values
        t_rhofas = values
        t_ntrslo = values
        t_fft2dhpf = values
        t_mfft_local1 = values
        t_field = values
        t_force = values
        t_charge = values
        t_beamln = values
        t_ptsmv = values
        t_init = values
        t_transp = values
        t_enlarge = values
        t_shrink = values
        t_guardsum = values
        t_boundint = values
        t_guardexch = values

        end subroutine construct_Timer

        subroutine starttime_Timer( sec0 )
        implicit none
        include 'mpif.h'
        double precision, intent( out ) :: sec0
        sec0 = MPI_WTIME()
        end subroutine starttime_Timer

        function elapsedtime_Timer( sec0 ) result( elapsed )
        implicit none
        include 'mpif.h'
        double precision, intent( in ) :: sec0
        double precision:: elapsed, t1
        t1 = MPI_WTIME()
        elapsed = t1 - sec0
        end function elapsedtime_Timer

        subroutine showtime_Timer()
        implicit none
        include 'mpif.h'
        call showtime_aux( "kvdist",        t_kvdist )
        call showtime_aux( "gaussn",      t_gaussn )
        call showtime_aux( "normdv",      t_normdv )
        call showtime_aux( "correct",     t_correct )
        call showtime_aux( "integ",       t_integ )
        call showtime_aux( "map1",        t_map1 )
        call showtime_aux( "map2",        t_map2 )
        call showtime_aux( "spch2d",      t_spch2d )
        call showtime_aux( "boundgeom",   t_boundgeom )
        call showtime_aux( "greenf",      t_greenf )
        call showtime_aux( "rhofas",      t_rhofas )
        call showtime_aux( "ntrslo",      t_ntrslo )
        call showtime_aux( "fft2dhpf",    t_fft2dhpf )
        call showtime_aux( "mfft_local1", t_mfft_local1 )
        call showtime_aux( "diag",        t_diag )
        call showtime_aux( "ufield",      t_field )
        call showtime_aux( "force",       t_force )
        call showtime_aux( "charge",      t_charge )
        call showtime_aux( "ubeamln",     t_beamln )
        call showtime_aux( "ptsmv",     t_ptsmv )
        call showtime_aux( "init",     t_init )
        call showtime_aux( "transp",     t_transp )
        call showtime_aux( "enlarge",     t_enlarge )
        call showtime_aux( "shrink",     t_shrink )
        call showtime_aux( "guardsum",     t_guardsum )
        call showtime_aux( "boundint",     t_boundint )
        call showtime_aux( "guardexch",     t_guardexch )
        end subroutine showtime_Timer

        subroutine showtime_aux( routine, sec )
        implicit none
        include 'mpif.h'  
        character (len = *), intent( in ) :: routine
        double precision, intent( in ) :: sec
        double precision:: msec
        integer :: myrank,ierr, procs
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        call MPI_REDUCE(sec,msec,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                        MPI_COMM_WORLD,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,procs,ierr)

        msec = msec / procs

        if(myrank.eq.0) then
          write( 6, "(a15, f11.5, a)" ) routine, msec, " seconds"
        endif
        end subroutine showtime_aux

      end module Timerclass
