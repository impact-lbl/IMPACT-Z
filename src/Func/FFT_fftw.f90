!----------------------------------------------------------------
! FFTW3 variant of FFTclass. This file is the drop-in replacement
! for FFT.f90 when the build is configured with USE_FFTW: the CMake
! and Makefile builds compile EITHER FFT.f90 (bundled FFTPACK 5.1)
! OR this file, never both, since both define module FFTclass.
!
! It contains two modules:
!   * FFTbackendclass -- the 1D FFT primitives four1/realft backed
!     by FFTW3 (fftw3.f03), keeping the Numerical-Recipes-style data
!     packing and scaling of the FFTPACK versions in FFT.f90:
!       four1:  isign=+1 -> exp(+i), otherwise exp(-i); unnormalized.
!       realft: isign=+1 -> unnormalized half-complex forward, NR
!               packing (data(1)=DC, data(2)=Nyquist); inverse is
!               n/2 times the true inverse.
!     Plans are cached per transform size and executed on persistent
!     FFTW-allocated (SIMD-aligned) buffers, so planner cost and
!     alignment constraints are confined to the first call for a
!     given size. Grid sizes are fixed during a run, so the cache
!     stays small.
!   * FFTclass -- the high-level transforms, byte-for-byte the same
!     routines as FFT.f90; the only difference from the FFTPACK build
!     is that four1/realft come from the FFTbackendclass above via
!     'use' rather than being defined inline.
!
! The build system selects this file with CMake -DUSE_FFTW=ON or
! make USE_FFTW=1, which also links against libfftw3.
!----------------------------------------------------------------
      module FFTbackendclass
      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f03'
      private
      public :: four1, realft, sinft, cosft1_fftpack

      ! Cached in-place complex-to-complex plan (one per size and
      ! transform direction).
      type :: cplan_t
        integer :: n = 0
        integer(C_INT) :: sign = 0
        type(C_PTR) :: plan = C_NULL_PTR
        complex(C_DOUBLE_COMPLEX), pointer :: buf(:) => null()
      end type cplan_t

      ! Cached real<->half-complex plan (r2c or c2r, one per size).
      type :: rplan_t
        integer :: n = 0
        type(C_PTR) :: plan = C_NULL_PTR
        real(C_DOUBLE), pointer :: rbuf(:) => null()
        complex(C_DOUBLE_COMPLEX), pointer :: cbuf(:) => null()
      end type rplan_t

      integer, parameter :: MAXPLANS = 64
      type(cplan_t), save, target :: cplans(MAXPLANS)
      type(rplan_t), save, target :: r2cplans(MAXPLANS)
      type(rplan_t), save, target :: c2rplans(MAXPLANS)
      integer, save :: ncplans = 0
      integer, save :: nr2cplans = 0
      integer, save :: nc2rplans = 0

      contains

      ! Complex 1D FFT of nn points packed NR-style into a real
      ! array of length 2*nn (alternating real/imaginary parts).
      ! Unnormalized in both directions; isign is the sign of the
      ! exponent, matching four1 in FFT_fftpack.f90.
      subroutine four1(data,nn,isign)
      integer, intent(in) :: nn, isign
      real*8 :: data(2*nn)
      type(cplan_t), pointer :: p
      integer(C_INT) :: sign
      integer :: i

      if (isign.eq.1) then
        sign = FFTW_BACKWARD
      else
        sign = FFTW_FORWARD
      endif
      p => get_cplan(nn,sign)
      do i = 1, nn
        p%buf(i) = dcmplx(data(2*i-1),data(2*i))
      enddo
      call fftw_execute_dft(p%plan,p%buf,p%buf)
      do i = 1, nn
        data(2*i-1) = dble(p%buf(i))
        data(2*i) = dimag(p%buf(i))
      enddo

      end subroutine four1

      ! Real 1D FFT of n points (n even), NR packing:
      ! data(1) = DC term, data(2) = Nyquist term, then
      ! (Re,Im) pairs of the positive frequencies with the
      ! exp(+i) convention of NR/the FFTPACK-based realft --
      ! the conjugate of FFTW's r2c convention, hence the sign
      ! flips on the imaginary parts below.
      ! Forward (isign=1) is unnormalized; the inverse returns
      ! n/2 times the true inverse, so callers multiply by
      ! scale*2 as with the FFTPACK backend.
      subroutine realft(data,n,isign)
      integer, intent(in) :: n, isign
      real*8 :: data(n)
      type(rplan_t), pointer :: p
      integer :: i

      if (isign.eq.1) then
        p => get_rplan(n,.true.)
        do i = 1, n
          p%rbuf(i) = data(i)
        enddo
        call fftw_execute_dft_r2c(p%plan,p%rbuf,p%cbuf)
        data(1) = dble(p%cbuf(1))
        data(2) = dble(p%cbuf(n/2+1))
        do i = 2, n/2
          data(2*i-1) = dble(p%cbuf(i))
          data(2*i) = -dimag(p%cbuf(i))
        enddo
      else
        p => get_rplan(n,.false.)
        p%cbuf(1) = dcmplx(data(1),0.0d0)
        p%cbuf(n/2+1) = dcmplx(data(2),0.0d0)
        do i = 2, n/2
          p%cbuf(i) = dcmplx(data(2*i-1),-data(2*i))
        enddo
        call fftw_execute_dft_c2r(p%plan,p%cbuf,p%rbuf)
        do i = 1, n
          data(i) = 0.5d0*p%rbuf(i)
        enddo
      endif

      end subroutine realft

      ! No callers in IMPACT-T; implemented only in the FFTPACK
      ! backend. Fail loudly rather than return wrong data.
      subroutine sinft(y,ny)
      integer, intent(in) :: ny
      real*8 :: y(ny)

      print*,'sinft is not implemented in the FFTW backend ',&
             '(build without USE_FFTW to use it).'
      stop

      end subroutine sinft

      subroutine cosft1_fftpack(y,n)
      integer, intent(in) :: n
      real*8 :: y(n+1)

      print*,'cosft1_fftpack is not implemented in the FFTW ',&
             'backend (build without USE_FFTW to use it).'
      stop

      end subroutine cosft1_fftpack

      ! Look up (or create on first use) the cached complex
      ! in-place plan for size n and direction sign.
      function get_cplan(n,sign) result(p)
      integer, intent(in) :: n
      integer(C_INT), intent(in) :: sign
      type(cplan_t), pointer :: p
      type(C_PTR) :: pbuf
      integer :: i

      do i = 1, ncplans
        if (cplans(i)%n.eq.n .and. cplans(i)%sign.eq.sign) then
          p => cplans(i)
          return
        endif
      enddo
      if (ncplans.ge.MAXPLANS) then
        print*,'FFT_fftw: plan cache overflow; increase MAXPLANS.'
        stop
      endif
      ncplans = ncplans + 1
      p => cplans(ncplans)
      p%n = n
      p%sign = sign
      pbuf = fftw_alloc_complex(int(n,C_SIZE_T))
      call c_f_pointer(pbuf,p%buf,[n])
      ! The plan is created on, and always executed with, this
      ! buffer, so FFTW's new-array alignment rules never apply.
      p%plan = fftw_plan_dft_1d(int(n,C_INT),p%buf,p%buf,sign,&
                                FFTW_ESTIMATE)

      end function get_cplan

      ! Look up (or create on first use) the cached r2c (forward
      ! .true.) or c2r (forward .false.) plan for size n.
      function get_rplan(n,forward) result(p)
      integer, intent(in) :: n
      logical, intent(in) :: forward
      type(rplan_t), pointer :: p
      type(C_PTR) :: prbuf,pcbuf
      integer :: i

      if (forward) then
        do i = 1, nr2cplans
          if (r2cplans(i)%n.eq.n) then
            p => r2cplans(i)
            return
          endif
        enddo
        if (nr2cplans.ge.MAXPLANS) then
          print*,'FFT_fftw: plan cache overflow; increase MAXPLANS.'
          stop
        endif
        nr2cplans = nr2cplans + 1
        p => r2cplans(nr2cplans)
      else
        do i = 1, nc2rplans
          if (c2rplans(i)%n.eq.n) then
            p => c2rplans(i)
            return
          endif
        enddo
        if (nc2rplans.ge.MAXPLANS) then
          print*,'FFT_fftw: plan cache overflow; increase MAXPLANS.'
          stop
        endif
        nc2rplans = nc2rplans + 1
        p => c2rplans(nc2rplans)
      endif
      p%n = n
      prbuf = fftw_alloc_real(int(n,C_SIZE_T))
      pcbuf = fftw_alloc_complex(int(n/2+1,C_SIZE_T))
      call c_f_pointer(prbuf,p%rbuf,[n])
      call c_f_pointer(pcbuf,p%cbuf,[n/2+1])
      if (forward) then
        p%plan = fftw_plan_dft_r2c_1d(int(n,C_INT),p%rbuf,p%cbuf,&
                                      FFTW_ESTIMATE)
      else
        p%plan = fftw_plan_dft_c2r_1d(int(n,C_INT),p%cbuf,p%rbuf,&
                                      FFTW_ESTIMATE)
      endif

      end function get_rplan

      end module FFTbackendclass

!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! FFTclass: Fourier function class in Math Function module of FUNCTION layer.
! 
! MODULE  : ... FFTclass
! VERSION : ... 2.0
!> @author
!> Ji Qiang 
! 
! DESCRIPTION: 
!> This class defines the 3d FFT transformation subject to
!> open or periodic conditions, Fourier Sine transformation,
!> Complex-Complex, Complex-Real, and Real-Complex FFT.
! Comments:
!----------------------------------------------------------------
      module FFTclass
      use Timerclass
      use Transposeclass
      use FFTbackendclass
      interface fftcrlocal_FFT
        module procedure fftcrlocal1_FFT,fftcrlocal2_FFT
      end interface
      interface fftrclocal_FFT
        module procedure fftrclocal1_FFT,fftrclocal2_FFT
      end interface
      contains
!----------------------------------------------------------------
! FFT for 3D open boundary conditions. 
! The original computational domain is doubled in each dimension
! to apply the FFT for the new domain.
        ! 3_D FFT.
        subroutine fft3d1_FFT(nx,ny,nz,nsizez,nsizey,nsizexy,&
                nsizeyz,ksign,scalex,scaley,scalez,x,xstable,xrtable,&
            ystable,yrtable,nprocrow,commrow,nproccol,commcol,comm2d,&
            myidx,myidy,xout)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer,intent(in) :: nprocrow,commrow,nproccol,commcol,comm2d
        integer,intent(in) :: ksign,myidx,myidy
        double precision, intent(in) :: scalex,scaley,scalez
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
        double precision, dimension(nx/2,nsizey,nsizez), intent(in) &
                            :: x
        double complex, dimension(nz,nsizexy,nsizeyz), intent(out) &
                            :: xout
        double precision, dimension(nx,nsizey) :: tmp1
        double complex, dimension(nx/2+1,nsizey) :: tmp10
        double complex, dimension(ny,nsizexy) :: tmp2
        double complex, dimension(nz,nsizexy) :: tmp3
        integer :: i,j,k
        double precision :: t0
        integer :: ierr,nxx
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        nxx = nx/2 + 1
        !FFTs along x dimensions: could be a lot of cache miss.
        allocate(x0(nx/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, nx/2
              tmp1(i,j) = x(i,j,k)
            enddo
            do i = nx/2+1, nx
              tmp1(i,j) = 0.0
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,nx,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, nx/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(ny/2,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
! yrtable needs to be changed.
! transpose between the x and y dimension.
        call trans3d_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
                     ystable,yrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(ny,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, ny/2
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = ny/2+1,ny
              tmp2(i,j) = (0.0,0.0)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

          do j = 1, nsizexy 
            do i = 1, ny
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(nz/2,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
! xrtable needs to be changed.
! transpose between serial ny (stored in i index) and z.
        call trans3d3_TRANSP(ny,nsizexy,nsizez,nsizeyz,x0,x1,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz/2)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, nz/2
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = nz/2+1,nz
              tmp3(i,j) = (0.0,0.0)
            enddo
          enddo

          !FFT along Z.
          call fftlocal_FFT(ksign,scalez,tmp3,nz,nsizexy)

          do j = 1, nsizexy
            do i = 1, nz
              xout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine fft3d1_FFT

        ! 3_D inverse FFT for open BCs.
        subroutine invfft3d1_FFT(nz,ny,nx,nsizexy,nsizeyz,nsizey,&
                nsizez,ksign,scalex,scaley,scalez,x,xstable,xrtable,&
            ystable,yrtable,nprocrow,commrow,nproccol,commcol,comm2d,&
            myidx,myidy,xout)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer,intent(in) :: nprocrow,commrow,nproccol,commcol,comm2d
        integer,intent(in) :: ksign,myidx,myidy
        double precision, intent(in) :: scalex,scaley,scalez
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
        double complex, dimension(nz,nsizexy,nsizeyz), intent(inout) &
                            :: x
        double precision, dimension(nx/2,nsizey,nsizez), intent(out) &
                            :: xout
        double complex, dimension(nz,nsizexy) :: tmp1
        double complex, dimension(ny,nsizexy) :: tmp2
        double complex, dimension(nx/2+1,nsizey) :: tmp3
        double precision, dimension(nx,nsizey) :: tmp30
        integer :: i,j,k
        double precision :: t0
        integer :: ierr,nxx
        double complex, allocatable, dimension(:,:,:) :: x0
        double complex, allocatable, dimension(:,:,:) :: x1

        call starttime_Timer(t0)

        nxx = nx/2 + 1
        allocate(x0(nz/2,nsizexy,nsizeyz))
        !FFTs along y and z dimensions: could be a lot of cache miss.
        do k = 1, nsizeyz
          do j = 1, nsizexy 
            do i = 1, nz
              tmp1(i,j) = x(i,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftlocal_FFT(ksign,scalez,tmp1,nz,nsizexy)

          do j = 1, nsizexy 
            do i = 1, nz/2
              x0(i,j,k) = tmp1(i,j)
!              x0(i,j,k) = tmp1(i,j)*scalez
            enddo
          enddo
        enddo

        allocate(x1(ny,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
! xrtable needs to be changed.
        call trans3d3_TRANSP(nz/2,nsizexy,nsizeyz,nsizez,x0,x1,nprocrow,&
                      xstable,xrtable,commrow,myidx,ny)

        deallocate(x0)
        allocate(x0(ny/2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, ny
              tmp2(i,j) = x1(i,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

          do j = 1, nsizexy 
            do i = 1, ny/2
              x0(i,j,k) = tmp2(i,j) 
!              x0(i,j,k) = tmp2(i,j)*scaley 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(nxx,nsizey,nsizez))
        call trans3d_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
                     ystable,yrtable,commcol,nsizez)
!        call MPI_BARRIER(commcol,ierr)
        deallocate(x0)

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nxx
              tmp3(i,j) = x1(i,j,k) 
            enddo
          enddo

          call fftcrlocal_FFT(ksign,scalex,tmp3,nx,nsizey,tmp30)

          do j = 1, nsizey
            do i = 1, nx/2
              xout(i,j,k) = tmp30(i,j)
!              xout(i,j,k) = tmp30(i,j)*scalex*2
            enddo
          enddo
        enddo

        deallocate(x1)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine invfft3d1_FFT

        ! 3_D inverse FFT for open BCs.
        subroutine invfft3d1Img_FFT(nz,ny,nx,nsizexy,nsizeyz,nsizey,&
                nsizez,ksign,scalex,scaley,scalez,x,xstable,xrtable,&
            ystable,yrtable,nprocrow,commrow,nproccol,commcol,comm2d,&
            myidx,myidy,xout)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer,intent(in) :: nprocrow,commrow,nproccol,commcol,comm2d
        integer,intent(in) :: ksign,myidx,myidy
        double precision, intent(in) :: scalex,scaley,scalez
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
        double complex, dimension(nz,nsizexy,nsizeyz), intent(inout) &
                            :: x
        double precision, dimension(nx/2,nsizey,nsizez), intent(out) &
                            :: xout
        double complex, dimension(nz,nsizexy) :: tmp1
        double complex, dimension(ny,nsizexy) :: tmp2
        double complex, dimension(nx/2+1,nsizey) :: tmp3
        double precision, dimension(nx,nsizey) :: tmp30
        integer :: i,j,k
        double precision :: t0
        integer :: ierr,nxx
        double complex, allocatable, dimension(:,:,:) :: x0
        double complex, allocatable, dimension(:,:,:) :: x1

        call starttime_Timer(t0)

        nxx = nx/2 + 1
        allocate(x0(nz/2,nsizexy,nsizeyz))
        !FFTs along y and z dimensions: could be a lot of cache miss.
        do k = 1, nsizeyz
          do j = 1, nsizexy 
            do i = 1, nz
              tmp1(i,j) = x(i,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftlocal_FFT(ksign,scalez,tmp1,nz,nsizexy)

          ! reverse the distribution along z.
          do j = 1, nsizexy 
            do i = 1, nz/2
!              x0(i,j,k) = tmp1(i,j)
              x0(nz/2-i+1,j,k) = tmp1(i,j)
!              x0(i,j,k) = tmp1(i,j)*scalez
            enddo
          enddo
        enddo

        allocate(x1(ny,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
! xrtable needs to be changed.
        call trans3d3_TRANSP(nz/2,nsizexy,nsizeyz,nsizez,x0,x1,nprocrow,&
                      xstable,xrtable,commrow,myidx,ny)

        deallocate(x0)
        allocate(x0(ny/2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, ny
              tmp2(i,j) = x1(i,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

          do j = 1, nsizexy 
            do i = 1, ny/2
              x0(i,j,k) = tmp2(i,j) 
!              x0(i,j,k) = tmp2(i,j)*scaley 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(nxx,nsizey,nsizez))
        call trans3d_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
                     ystable,yrtable,commcol,nsizez)
!        call MPI_BARRIER(commcol,ierr)
        deallocate(x0)

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nxx
              tmp3(i,j) = x1(i,j,k) 
            enddo
          enddo

          call fftcrlocal_FFT(ksign,scalex,tmp3,nx,nsizey,tmp30)

          do j = 1, nsizey
            do i = 1, nx/2
              xout(i,j,k) = tmp30(i,j)
!              xout(i,j,k) = tmp30(i,j)*scalex*2
            enddo
          enddo
        enddo

        deallocate(x1)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine invfft3d1Img_FFT

      !used to find the first derivative after FFT.
      ! Subroutine to perform 1D FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftlocal0_FFT(ksign,scale,x,ny,nsizex)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(inout) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           do i = 1, ny
             tempi(i) = x(i,j)
           enddo
           call four1(tempi,ny/2,ksign)
           do i = 1, ny
             x(i,j) = tempi(i)*scale
           enddo
      end do

      end subroutine fftlocal0_FFT

      ! Subroutine to perform 1D FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftlocal_FFT(ksign,scale,x,ny,nsizex)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double complex, dimension(ny,nsizex), intent(inout) :: x
      real*8, dimension(2*ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           do i = 1, ny
             tempi(2*i-1) = real(x(i,j))
             tempi(2*i) = aimag(x(i,j))
           enddo
           call four1(tempi,ny,ksign)
           do i = 1, ny
             x(i,j) = dcmplx(tempi(2*i-1),tempi(2*i))*scale
             !x(i,j) = dcmplx(tempi(2*i-1),tempi(2*i))
           enddo
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftlocal_FFT

      ! Subroutine to perform 1D real to complex
      ! FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftrclocal1_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(in) :: x
      double complex, dimension(ny/2+1,nsizex), intent(out) :: y
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
         tempi = x(:,j)
         call realft(tempi,ny,ksign)
         y(1,j) = dcmplx(tempi(1),0.0d0)*scale
         y(ny/2+1,j) = dcmplx(tempi(2),0.0d0)*scale
         do i = 2,ny/2
           y(i,j) = dcmplx(tempi(2*i-1),tempi(2*i))*scale
         enddo
         !y(1,j) = cmplx(tempi(1),0.0)
         !y(ny/2+1,j) = cmplx(tempi(2),0.0)
         !do i = 2,ny/2
         !  y(i,j) = cmplx(tempi(2*i-1),tempi(2*i))
         !enddo
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftrclocal1_FFT

      subroutine fftrclocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(in) :: x
      double precision, dimension(ny,nsizex), intent(out) :: y
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
         tempi = x(:,j)
         call realft(tempi,ny,ksign)
         y(:,j) = tempi*scale
         !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftrclocal2_FFT

      ! Subroutine to perform 1D complex to real 
      ! FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftcrlocal1_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(out) :: y
      double complex, dimension(ny/2+1,nsizex), intent(in) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           tempi(1) = real(x(1,j))
           tempi(2) = real(x(ny/2+1,j))
           do i = 2, ny/2
             tempi(2*i-1) = real(x(i,j))
             tempi(2*i) = aimag(x(i,j))
           enddo
           call realft(tempi,ny,ksign)
           y(:,j) = tempi*scale*2
           !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftcrlocal1_FFT

      subroutine fftcrlocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(out) :: y
      double precision, dimension(ny,nsizex), intent(in) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           tempi = x(:,j)
           call realft(tempi,ny,ksign)
           y(:,j) = tempi*scale*2
           !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftcrlocal2_FFT

      end module FFTclass
