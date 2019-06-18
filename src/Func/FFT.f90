!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! FFTclass: Fourier function class in Math Function module of FUNCTION layer.
! Version: 2.0
! Author: Ji Qiang 
! Description: This class defines the 3d FFT transformation subject to
!              open or periodic conditions, Fourier Sine transformation,
!              Complex-Complex, Complex-Real, and Real-Complex FFT.
! Comments:
!----------------------------------------------------------------
      module FFTclass
      use Timerclass
      use Transposeclass
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

!File name: fftnr_Fromfftpack.f90
!Implement FFT1D, REAL FFT1D, SINE transform and cosine transform with 
!FFTPACK 5.1 and keep the same setting style as NR
! Create Date: 01/27/2015

! Modify Record:
! <author>       <date>         <version>       <desc>
!----------------------------------------------------------------------------------------------------
! Dawei Zheng    01/27/2015        1.0             create the file!


!Fast fourier transform
subroutine four1(data,nn,isign)
  integer isign, nn
  real*8 data(2*nn)
!  real ( kind = 8 ) c(2*nn)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenhc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  lenwrk = 2 * nn
  lensav = 2 * nn + int ( log ( real ( nn, kind = 8 ) ) / log ( 2.0E+00 ) ) + 4

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )
!Initialization
  call cfft1i ( nn, wsave, lensav, ier )
   inc = 1
  lenhc = nn

if (isign.eq.1) then
!fft backward with fftpack	
       call cfft1b ( nn, inc, data, lenhc, wsave, lensav, work, lenwrk, ier )
else
!fft forward with fftpack
	call cfft1f ( nn, inc, data, lenhc, wsave, lensav, work, lenwrk, ier )
do i =1,2*nn
data(i)= data(i)*dble(nn)
enddo
endif

  deallocate ( work )
  deallocate ( wsave )

return

end subroutine four1

!real fast Fourier tranform 
subroutine realft(data,n,isign)
  integer isign, n
  real*8 data(n)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  lensav = n + int ( log ( real ( n, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4
  lenwrk = n

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )
!Initialization
  call rfft1i ( n, wsave, lensav, ier )
   inc = 1
  lenr = n
      
if (isign.eq.1) then
!transfer the data to the fftpack style
	call Rvec_nr2pack_rfft1f(r,n,data)
!real forward fast Fourier transform with fftpack
       call rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
!transfer the data back to the nr style 
	call Rvec_pack2nr_rfft1f (r, n, data)

else
       call Rvec_nr2pack_rfft1b(r,n,data)
!real backward fast Fourier transform with fftpack
	call rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
!transfer the data back to the nr style 
	call Rvec_pack2nr_rfft1b (r, n, data)
endif

  deallocate ( work )
  deallocate ( wsave )

return

end subroutine realft

!fast discrete sine transform
subroutine sinft(y,ny)
integer ny, n
real*8 y(ny)
 real ( kind = 8 ) r(ny-1)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk 
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

!
!  Set work arrays.
!
       n=ny-1
  lenwrk = 2 * ( n + 1 )
  lensav = n / 2 + n + int ( log ( real ( n, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4

  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )
!Initialization
  call sint1i ( n, wsave, lensav, ier )

  inc = 1
  lenr = n
!transfer the data to the fftpack style
  call Rvec_nr2pack_sint1f(r, n, y)
!fast discrete sine transform with fftpack
  call sint1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
!transfer the data back to the nr style 
  call Rvec_pack2nr_sint1f (r, n, y)

  deallocate ( work )
  deallocate ( wsave )

return
end subroutine sinft

!fast discrete cosine transform
subroutine cosft1_fftpack(y,n)
integer n,ny
real*8 y(n+1)

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n+1)

  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave
  ny=n+1

  lensav = 2 * ny + int ( log ( real ( ny, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4
  lenwrk = ny

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )
!Initialization
  call cost1i ( ny, wsave, lensav, ier )

   inc = 1
  lenr = ny
!transfer the data to the fftpack style
  call Rvec_nr2pack_gnl(r,ny,y)
!fast discrete cosine transform with fftpack
  call cost1f ( ny, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
!transfer the data back to the nr style 
  call Rvec_pack2nr_cosft1(r,ny,y)
  deallocate ( work )
  deallocate ( wsave )

  return

end subroutine cosft1_fftpack

!transfer real vector from fftpack style to nr style for (forward) cost1D
subroutine Rvec_pack2nr_cosft1 (r,ny,y)
integer ny
integer i
real ( kind = 8 ) r(ny)
real*8 y(ny)
y(1)= r(1)*dble(ny-1)
do i= 2,ny-1
  y(i)= r(i)*((ny-1)/2)
enddo
y(ny)=r(ny)*dble(ny-1)
return
end subroutine Rvec_pack2nr_cosft1

!transfer real vector from nr style to fftpack style for (forward) sint1D
subroutine Rvec_nr2pack_sint1f(r,n,y)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 y(n+1)
!r(1)=0.0!
do i= 1, n 
  r(i)= y(i+1)
enddo
return
end subroutine Rvec_nr2pack_sint1f

!transfer real vector from fftpack style to nr style for (forward) sint 1D
subroutine Rvec_pack2nr_sint1f (r,n,y)
integer n, ny
integer i
real ( kind = 8 ) r(n)
real*8 y(n+1)
ny=n+1
y(1)= dble(0.0)
do i= 2,ny
y(i)= r(i-1)*dble(ny/2)
enddo
!y(1)=0.0!
return
end subroutine Rvec_pack2nr_sint1f

!transfer real vector from nr style to fftpack style for forward real fft 1D
subroutine Rvec_nr2pack_rfft1f(r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)
do i= 1, n 
r(i)= data(i)
enddo
return
end subroutine Rvec_nr2pack_rfft1f

!transfer real vector from pack style to nr style for forward real fft 1D
subroutine Rvec_pack2nr_rfft1f (r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)
data(1)= r(1)*dble(n)
data(2)= r(n)*dble(n)
do i= 3,n
data(i)= r(i-1)*dble(n)/dble(2)
enddo
return
end subroutine Rvec_pack2nr_rfft1f

!transfer real vector from nr style to fftpack style for backward real fft 1D
subroutine Rvec_nr2pack_rfft1b(r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)
r(1)=data(1)/dble(n)
r(n)=data(2)/dble(n)
do i= 2, n-1 
r(i)= data(i+1)*dble(2)/dble(n)
enddo
return
end subroutine Rvec_nr2pack_rfft1b

!transfer real vector from pack style to nr style for backward real fft 1D
subroutine Rvec_pack2nr_rfft1b (r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)

do i= 1,n
data(i)= r(i)*dble(n)/dble(2)
enddo
return
end subroutine Rvec_pack2nr_rfft1b

!transfer real vector from nr style to fftpack style generally 
subroutine Rvec_nr2pack_gnl(r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)
do i= 1, n 
r(i)= data(i)
enddo
return
end subroutine Rvec_nr2pack_gnl

!tansfer real vector from fftpack style to nr style generally 
subroutine Rvec_pack2nr_gnl (r,n,data)
integer n
integer i
real ( kind = 8 ) r(n)
real*8 data(n)
do i= 1,n
data(i)= r(i)
enddo
return
end subroutine Rvec_pack2nr_gnl

subroutine c1f2kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F2KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    01 August 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)



  if ( ido <= 1 .and. na /= 1 ) then

    do k=1,l1
      chold1 = cc(1,k,1,1)+cc(1,k,1,2)
      cc(1,k,1,2) = cc(1,k,1,1)-cc(1,k,1,2)
      cc(1,k,1,1) = chold1
      chold2 = cc(2,k,1,1)+cc(2,k,1,2)
      cc(2,k,1,2) = cc(2,k,1,1)-cc(2,k,1,2)
      cc(2,k,1,1) = chold2
    end do

    return

  end if

  do k=1,l1
    ch(1,k,1,1) = cc(1,k,1,1)+cc(1,k,1,2)
    ch(1,k,2,1) = cc(1,k,1,1)-cc(1,k,1,2)
    ch(2,k,1,1) = cc(2,k,1,1)+cc(2,k,1,2)
    ch(2,k,2,1) = cc(2,k,1,1)-cc(2,k,1,2)
  end do

  do i=2,ido
    do k=1,l1
      ch(1,k,1,i) = cc(1,k,i,1)+cc(1,k,i,2)
      tr2 = cc(1,k,i,1)-cc(1,k,i,2)
      ch(2,k,1,i) = cc(2,k,i,1)+cc(2,k,i,2)
      ti2 = cc(2,k,i,1)-cc(2,k,i,2)
      ch(2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
      ch(1,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
    end do
  end do

  return
end subroutine c1f2kb

subroutine c1f2kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F2KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

      if (  1 < ido ) go to 102

      sn = 1.0E+00 / real ( 2 * l1, kind = 8 )

      if (na == 1) go to 106

      do k=1,l1
         chold1 = sn*(cc(1,k,1,1)+cc(1,k,1,2))
         cc(1,k,1,2) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
         cc(1,k,1,1) = chold1
         chold2 = sn*(cc(2,k,1,1)+cc(2,k,1,2))
         cc(2,k,1,2) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
         cc(2,k,1,1) = chold2
      end do

      return

  106 do k=1,l1
         ch(1,k,1,1) = sn*(cc(1,k,1,1)+cc(1,k,1,2))
         ch(1,k,2,1) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
         ch(2,k,1,1) = sn*(cc(2,k,1,1)+cc(2,k,1,2))
         ch(2,k,2,1) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
      end do

      return

  102 do k=1,l1
         ch(1,k,1,1) = cc(1,k,1,1)+cc(1,k,1,2)
         ch(1,k,2,1) = cc(1,k,1,1)-cc(1,k,1,2)
         ch(2,k,1,1) = cc(2,k,1,1)+cc(2,k,1,2)
         ch(2,k,2,1) = cc(2,k,1,1)-cc(2,k,1,2)
      end do

      do i=2,ido
         do k=1,l1
            ch(1,k,1,i) = cc(1,k,i,1)+cc(1,k,i,2)
            tr2 = cc(1,k,i,1)-cc(1,k,i,2)
            ch(2,k,1,i) = cc(2,k,i,1)+cc(2,k,i,2)
            ti2 = cc(2,k,i,1)-cc(2,k,i,2)
            ch(2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
            ch(1,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
         end do
      end do

  return
end subroutine c1f2kf
subroutine c1f3kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F3KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,3)
  real ( kind = 8 ) ch(in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ), parameter :: taui =  0.866025403784439E+00
  real ( kind = 8 ), parameter :: taur = -0.5E+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

      if ( 1 < ido .or. na == 1) go to 102

      do k=1,l1
         tr2 = cc(1,k,1,2)+cc(1,k,1,3)
         cr2 = cc(1,k,1,1)+taur*tr2
         cc(1,k,1,1) = cc(1,k,1,1)+tr2
         ti2 = cc(2,k,1,2)+cc(2,k,1,3)
         ci2 = cc(2,k,1,1)+taur*ti2
         cc(2,k,1,1) = cc(2,k,1,1)+ti2
         cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
         ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
         cc(1,k,1,2) = cr2-ci3
         cc(1,k,1,3) = cr2+ci3
         cc(2,k,1,2) = ci2+cr3
         cc(2,k,1,3) = ci2-cr3
      end do

      return

  102 continue

      do k=1,l1
         tr2 = cc(1,k,1,2)+cc(1,k,1,3)
         cr2 = cc(1,k,1,1)+taur*tr2
         ch(1,k,1,1) = cc(1,k,1,1)+tr2
         ti2 = cc(2,k,1,2)+cc(2,k,1,3)
         ci2 = cc(2,k,1,1)+taur*ti2
         ch(2,k,1,1) = cc(2,k,1,1)+ti2
         cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
         ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
         ch(1,k,2,1) = cr2-ci3
         ch(1,k,3,1) = cr2+ci3
         ch(2,k,2,1) = ci2+cr3
         ch(2,k,3,1) = ci2-cr3
      end do

      do i=2,ido
        do k=1,l1
            tr2 = cc(1,k,i,2)+cc(1,k,i,3)
            cr2 = cc(1,k,i,1)+taur*tr2
            ch(1,k,1,i) = cc(1,k,i,1)+tr2
            ti2 = cc(2,k,i,2)+cc(2,k,i,3)
            ci2 = cc(2,k,i,1)+taur*ti2
            ch(2,k,1,i) = cc(2,k,i,1)+ti2
            cr3 = taui*(cc(1,k,i,2)-cc(1,k,i,3))
            ci3 = taui*(cc(2,k,i,2)-cc(2,k,i,3))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
            ch(1,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
            ch(2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
            ch(1,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
         end do
       end do

  return
end subroutine c1f3kb
subroutine c1f3kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F3KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,3)
  real ( kind = 8 ) ch(in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ), parameter :: taui = -0.866025403784439E+00
  real ( kind = 8 ), parameter :: taur = -0.5E+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 3 * l1, kind = 8 )
      if (na == 1) go to 106

      do k=1,l1
         tr2 = cc(1,k,1,2)+cc(1,k,1,3)
         cr2 = cc(1,k,1,1)+taur*tr2
         cc(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
         ti2 = cc(2,k,1,2)+cc(2,k,1,3)
         ci2 = cc(2,k,1,1)+taur*ti2
         cc(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
         cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
         ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
         cc(1,k,1,2) = sn*(cr2-ci3)
         cc(1,k,1,3) = sn*(cr2+ci3)
         cc(2,k,1,2) = sn*(ci2+cr3)
         cc(2,k,1,3) = sn*(ci2-cr3)
      end do

      return

  106 do k=1,l1
         tr2 = cc(1,k,1,2)+cc(1,k,1,3)
         cr2 = cc(1,k,1,1)+taur*tr2
         ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
         ti2 = cc(2,k,1,2)+cc(2,k,1,3)
         ci2 = cc(2,k,1,1)+taur*ti2
         ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
         cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
         ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
         ch(1,k,2,1) = sn*(cr2-ci3)
         ch(1,k,3,1) = sn*(cr2+ci3)
         ch(2,k,2,1) = sn*(ci2+cr3)
         ch(2,k,3,1) = sn*(ci2-cr3)
      end do

      return

  102 do 103 k=1,l1
         tr2 = cc(1,k,1,2)+cc(1,k,1,3)
         cr2 = cc(1,k,1,1)+taur*tr2
         ch(1,k,1,1) = cc(1,k,1,1)+tr2
         ti2 = cc(2,k,1,2)+cc(2,k,1,3)
         ci2 = cc(2,k,1,1)+taur*ti2
         ch(2,k,1,1) = cc(2,k,1,1)+ti2
         cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
         ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
         ch(1,k,2,1) = cr2-ci3
         ch(1,k,3,1) = cr2+ci3
         ch(2,k,2,1) = ci2+cr3
         ch(2,k,3,1) = ci2-cr3
  103 continue

      do 105 i=2,ido
        do 104 k=1,l1
            tr2 = cc(1,k,i,2)+cc(1,k,i,3)
            cr2 = cc(1,k,i,1)+taur*tr2
            ch(1,k,1,i) = cc(1,k,i,1)+tr2
            ti2 = cc(2,k,i,2)+cc(2,k,i,3)
            ci2 = cc(2,k,i,1)+taur*ti2
            ch(2,k,1,i) = cc(2,k,i,1)+ti2
            cr3 = taui*(cc(1,k,i,2)-cc(1,k,i,3))
            ci3 = taui*(cc(2,k,i,2)-cc(2,k,i,3))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
  104    continue
  105 continue

  return
end subroutine c1f3kf
subroutine c1f4kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F4KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,4)
  real ( kind = 8 ) ch(in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

      if ( 1 < ido .or. na == 1) go to 102

      do k=1,l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,4)-cc(2,k,1,2)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,2)-cc(1,k,1,4)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         cc(1,k,1,1) = tr2+tr3
         cc(1,k,1,3) = tr2-tr3
         cc(2,k,1,1) = ti2+ti3
         cc(2,k,1,3) = ti2-ti3
         cc(1,k,1,2) = tr1+tr4
         cc(1,k,1,4) = tr1-tr4
         cc(2,k,1,2) = ti1+ti4
         cc(2,k,1,4) = ti1-ti4
      end do

      return

  102 do 103 k=1,l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,4)-cc(2,k,1,2)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,2)-cc(1,k,1,4)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = tr2+tr3
         ch(1,k,3,1) = tr2-tr3
         ch(2,k,1,1) = ti2+ti3
         ch(2,k,3,1) = ti2-ti3
         ch(1,k,2,1) = tr1+tr4
         ch(1,k,4,1) = tr1-tr4
         ch(2,k,2,1) = ti1+ti4
         ch(2,k,4,1) = ti1-ti4
  103 continue

      do 105 i=2,ido
         do 104 k=1,l1
            ti1 = cc(2,k,i,1)-cc(2,k,i,3)
            ti2 = cc(2,k,i,1)+cc(2,k,i,3)
            ti3 = cc(2,k,i,2)+cc(2,k,i,4)
            tr4 = cc(2,k,i,4)-cc(2,k,i,2)
            tr1 = cc(1,k,i,1)-cc(1,k,i,3)
            tr2 = cc(1,k,i,1)+cc(1,k,i,3)
            ti4 = cc(1,k,i,2)-cc(1,k,i,4)
            tr3 = cc(1,k,i,2)+cc(1,k,i,4)
            ch(1,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
            ch(2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
            ch(1,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
            ch(2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
            ch(1,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
            ch(2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end subroutine  c1f4kb
subroutine c1f4kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F4KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,4)
  real ( kind = 8 ) ch(in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 4 * l1, kind = 8 )
      if (na == 1) go to 106

      do k=1,l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         cc(1,k,1,1) = sn*(tr2+tr3)
         cc(1,k,1,3) = sn*(tr2-tr3)
         cc(2,k,1,1) = sn*(ti2+ti3)
         cc(2,k,1,3) = sn*(ti2-ti3)
         cc(1,k,1,2) = sn*(tr1+tr4)
         cc(1,k,1,4) = sn*(tr1-tr4)
         cc(2,k,1,2) = sn*(ti1+ti4)
         cc(2,k,1,4) = sn*(ti1-ti4)
      end do

      return

  106 do 107 k=1,l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(tr2+tr3)
         ch(1,k,3,1) = sn*(tr2-tr3)
         ch(2,k,1,1) = sn*(ti2+ti3)
         ch(2,k,3,1) = sn*(ti2-ti3)
         ch(1,k,2,1) = sn*(tr1+tr4)
         ch(1,k,4,1) = sn*(tr1-tr4)
         ch(2,k,2,1) = sn*(ti1+ti4)
         ch(2,k,4,1) = sn*(ti1-ti4)
  107 continue

      return

  102 do 103 k=1,l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = tr2+tr3
         ch(1,k,3,1) = tr2-tr3
         ch(2,k,1,1) = ti2+ti3
         ch(2,k,3,1) = ti2-ti3
         ch(1,k,2,1) = tr1+tr4
         ch(1,k,4,1) = tr1-tr4
         ch(2,k,2,1) = ti1+ti4
         ch(2,k,4,1) = ti1-ti4
  103 continue
      do 105 i=2,ido
         do 104 k=1,l1
            ti1 = cc(2,k,i,1)-cc(2,k,i,3)
            ti2 = cc(2,k,i,1)+cc(2,k,i,3)
            ti3 = cc(2,k,i,2)+cc(2,k,i,4)
            tr4 = cc(2,k,i,2)-cc(2,k,i,4)
            tr1 = cc(1,k,i,1)-cc(1,k,i,3)
            tr2 = cc(1,k,i,1)+cc(1,k,i,3)
            ti4 = cc(1,k,i,4)-cc(1,k,i,2)
            tr3 = cc(1,k,i,2)+cc(1,k,i,4)
            ch(1,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
            ch(2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
            ch(1,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
            ch(2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
            ch(1,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
            ch(2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end subroutine c1f4kf
subroutine c1f5kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F5KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,5)
  real ( kind = 8 ) ch(in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 =  0.9510565162951536E+00
  real ( kind = 8 ), parameter :: ti12 =  0.5877852522924731E+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474E+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474E+00
  real ( kind = 8 ) wa(ido,4,2)

      if ( 1 < ido .or. na == 1) go to 102

      do k=1,l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         chold1 = cc(1,k,1,1)+tr2+tr3
         chold2 = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,k,1,1) = chold1
         cc(2,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,k,1,2) = cr2-ci5
         cc(1,k,1,5) = cr2+ci5
         cc(2,k,1,2) = ci2+cr5
         cc(2,k,1,3) = ci3+cr4
         cc(1,k,1,3) = cr3-ci4
         cc(1,k,1,4) = cr3+ci4
         cc(2,k,1,4) = ci3-cr4
         cc(2,k,1,5) = ci2-cr5
      end do

      return

  102 do 103 k=1,l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
         ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = cr2-ci5
         ch(1,k,5,1) = cr2+ci5
         ch(2,k,2,1) = ci2+cr5
         ch(2,k,3,1) = ci3+cr4
         ch(1,k,3,1) = cr3-ci4
         ch(1,k,4,1) = cr3+ci4
         ch(2,k,4,1) = ci3-cr4
         ch(2,k,5,1) = ci2-cr5
  103 continue

      do 105 i=2,ido
         do 104 k=1,l1
            ti5 = cc(2,k,i,2)-cc(2,k,i,5)
            ti2 = cc(2,k,i,2)+cc(2,k,i,5)
            ti4 = cc(2,k,i,3)-cc(2,k,i,4)
            ti3 = cc(2,k,i,3)+cc(2,k,i,4)
            tr5 = cc(1,k,i,2)-cc(1,k,i,5)
            tr2 = cc(1,k,i,2)+cc(1,k,i,5)
            tr4 = cc(1,k,i,3)-cc(1,k,i,4)
            tr3 = cc(1,k,i,3)+cc(1,k,i,4)
            ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
            ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
            cr2 = cc(1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
            ch(2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
            ch(1,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
            ch(2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
            ch(1,k,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
            ch(2,k,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
            ch(1,k,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
            ch(2,k,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end subroutine c1f5kb
subroutine c1f5kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F5KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,5)
  real ( kind = 8 ) ch(in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 = -0.9510565162951536E+00
  real ( kind = 8 ), parameter :: ti12 = -0.5877852522924731E+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474E+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474E+00
  real ( kind = 8 ) wa(ido,4,2)

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 5 * l1, kind = 8 )
      if (na == 1) go to 106

      do k=1,l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
         chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,k,1,1) = chold1
         cc(2,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,k,1,2) = sn*(cr2-ci5)
         cc(1,k,1,5) = sn*(cr2+ci5)
         cc(2,k,1,2) = sn*(ci2+cr5)
         cc(2,k,1,3) = sn*(ci3+cr4)
         cc(1,k,1,3) = sn*(cr3-ci4)
         cc(1,k,1,4) = sn*(cr3+ci4)
         cc(2,k,1,4) = sn*(ci3-cr4)
         cc(2,k,1,5) = sn*(ci2-cr5)
      end do

      return

  106 do 107 k=1,l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
         ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = sn*(cr2-ci5)
         ch(1,k,5,1) = sn*(cr2+ci5)
         ch(2,k,2,1) = sn*(ci2+cr5)
         ch(2,k,3,1) = sn*(ci3+cr4)
         ch(1,k,3,1) = sn*(cr3-ci4)
         ch(1,k,4,1) = sn*(cr3+ci4)
         ch(2,k,4,1) = sn*(ci3-cr4)
         ch(2,k,5,1) = sn*(ci2-cr5)
  107 continue

      return

  102 do 103 k=1,l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
         ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = cr2-ci5
         ch(1,k,5,1) = cr2+ci5
         ch(2,k,2,1) = ci2+cr5
         ch(2,k,3,1) = ci3+cr4
         ch(1,k,3,1) = cr3-ci4
         ch(1,k,4,1) = cr3+ci4
         ch(2,k,4,1) = ci3-cr4
         ch(2,k,5,1) = ci2-cr5
  103 continue

      do 105 i=2,ido
         do 104 k=1,l1
            ti5 = cc(2,k,i,2)-cc(2,k,i,5)
            ti2 = cc(2,k,i,2)+cc(2,k,i,5)
            ti4 = cc(2,k,i,3)-cc(2,k,i,4)
            ti3 = cc(2,k,i,3)+cc(2,k,i,4)
            tr5 = cc(1,k,i,2)-cc(1,k,i,5)
            tr2 = cc(1,k,i,2)+cc(1,k,i,5)
            tr4 = cc(1,k,i,3)-cc(1,k,i,4)
            tr3 = cc(1,k,i,3)+cc(1,k,i,4)
            ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
            ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
            cr2 = cc(1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
            ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
            ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
            ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
            ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end subroutine c1f5kf
subroutine c1fgkb ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )

!*****************************************************************************80
!
!! C1FGKB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) na
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      ipp2 = ip+2
      ipph = (ip+1)/2

      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
  112    continue
  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = wa(1,l-1,2)*ch1(1,ki,ip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = wa(1,l-1,2)*ch1(2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = wa(1,idlj,2)
            do 114 ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
  114       continue
  115    continue
  116 continue

      if( 1 < ido .or. na == 1) go to 136

      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
            chold1 = cc1(1,ki,j)-cc1(2,ki,jc)
            chold2 = cc1(1,ki,j)+cc1(2,ki,jc)
            cc1(1,ki,j) = chold1
            cc1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
            cc1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
            cc1(1,ki,jc) = chold2
  119    continue
  120 continue
      return

  136 do 137 ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
  137 continue

      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
            ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
  134    continue
  135 continue

      if (ido == 1) then
        return
      end if

      do 131 i=1,ido
         do 130 k=1,l1
            cc(1,k,1,i) = ch(1,k,i,1)
            cc(2,k,1,i) = ch(2,k,i,1)
  130    continue
  131 continue

      do 123 j=2,ip
         do 122 k=1,l1
            cc(1,k,j,1) = ch(1,k,1,j)
            cc(2,k,j,1) = ch(2,k,1,j)
  122    continue
  123 continue

      do 126 j=2,ip
         do 125 i=2,ido
            do 124 k=1,l1
               cc(1,k,j,i) = wa(i,j-1,1)*ch(1,k,i,j) &
                            -wa(i,j-1,2)*ch(2,k,i,j)
               cc(2,k,j,i) = wa(i,j-1,1)*ch(2,k,i,j) &
                            +wa(i,j-1,2)*ch(1,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end subroutine c1fgkb 
subroutine c1fgkf ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )

!*****************************************************************************80
!
!! C1FGKF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      ipp2 = ip+2
      ipph = (ip+1)/2
      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
  112    continue
  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,ip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do 114 ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
  114       continue
  115    continue
  116 continue

      if ( 1 < ido ) go to 136
      sn = 1.0E+00 / real ( ip * l1, kind = 8 )
      if (na == 1) go to 146
      do 149 ki=1,lid
         cc1(1,ki,1) = sn*cc1(1,ki,1)
         cc1(2,ki,1) = sn*cc1(2,ki,1)
  149 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
            chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            cc1(1,ki,j) = chold1
            cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
            cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            cc1(1,ki,jc) = chold2
  119    continue
  120 continue
      return

  146 do 147 ki=1,lid
         ch1(1,ki,1) = sn*cc1(1,ki,1)
         ch1(2,ki,1) = sn*cc1(2,ki,1)
  147 continue
      do 145 j=2,ipph
         jc = ipp2-j
         do 144 ki=1,lid
            ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
  144    continue
  145 continue
      return

  136 do 137 ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
  137 continue
      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
            ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
            ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
  134    continue
  135 continue
      do 131 i=1,ido
         do 130 k=1,l1
            cc(1,k,1,i) = ch(1,k,i,1)
            cc(2,k,1,i) = ch(2,k,i,1)
  130    continue
  131 continue
      do 123 j=2,ip
         do 122 k=1,l1
            cc(1,k,j,1) = ch(1,k,1,j)
            cc(2,k,j,1) = ch(2,k,1,j)
  122    continue
  123 continue
      do 126 j=2,ip
         do 125 i=2,ido
            do 124 k=1,l1
               cc(1,k,j,i) = wa(i,j-1,1)*ch(1,k,i,j) &
                            +wa(i,j-1,2)*ch(2,k,i,j)
               cc(2,k,j,i) = wa(i,j-1,1)*ch(2,k,i,j) &
                            -wa(i,j-1,2)*ch(1,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end subroutine c1fgkf
subroutine c1fm1b ( n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! C1FM1B is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

!  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      inc2 = inc+inc
      nf = fnf
      na = 0
      l1 = 1
      iw = 1

      do k1=1,nf
         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4) 
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call c1f2kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   62    call c1f2kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   53    call c1f3kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   63    call c1f3kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   54    call c1f4kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   64    call c1f4kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   55    call c1f5kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   65    call c1f5kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   56    call c1fgkb (ido,ip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
         go to 120
   66    call c1fgkb (ido,ip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
      end do

  return
end subroutine c1fm1b
subroutine c1fm1f ( n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! C1FM1F is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

!  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      inc2 = inc+inc
      nf = fnf
      na = 0
      l1 = 1
      iw = 1

      do k1=1,nf

         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4)
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call c1f2kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   62    call c1f2kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   53    call c1f3kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   63    call c1f3kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   54    call c1f4kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   64    call c1f4kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   55    call c1f5kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   65    call c1f5kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   56    call c1fgkf (ido,ip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
         go to 120
   66    call c1fgkf (ido,ip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
      end do

  return
end subroutine c1fm1f
subroutine cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! CFFT1B: complex single precision backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1B computes the one-dimensional Fourier transform of a single 
!    periodic sequence within a complex array.  This transform is referred 
!    to as the backward transform or Fourier synthesis, transforming the
!    sequence from spectral to physical space.
!
!    This transform is normalized since a call to CFFT1B followed
!    by a call to CFFT1F (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.  
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to CFFT1I before the first call to routine CFFT1F 
!    or CFFT1B for a given transform length N.  WSAVE's contents may be 
!    re-used for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  !complex ( kind = 8 ) c(lenc)
   real    (kind = 8 ) c(2*lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenc < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'cfft1b ', 4 )
  else if ( lensav < 2 * n + int ( log ( real ( n, kind = 8 ) ) &
    / log ( 2.0E+00 ) ) + 4 ) then
    ier = 2
    call xerfft ( 'cfft1b ', 6 )
  else if ( lenwrk < 2 * n ) then
    ier = 3
    call xerfft ( 'cfft1b ', 8 )
  end if

  if ( n == 1 ) then
    return
  end if

  iw1 = n + n + 1

  call c1fm1b ( n, inc, c, work, wsave, wsave(iw1), wsave(iw1+1) )

  return
end subroutine cfft1b 
subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! CFFT1F: complex single precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1F computes the one-dimensional Fourier transform of a single 
!    periodic sequence within a complex array.  This transform is referred 
!    to as the forward transform or Fourier analysis, transforming the 
!    sequence from physical to spectral space.
!
!    This transform is normalized since a call to CFFT1F followed
!    by a call to CFFT1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.  
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to CFFT1I before the first call to routine CFFT1F 
!    or CFFT1B for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

 ! complex ( kind = 8 ) c(lenc)
   real ( kind = 8 ) c(2*lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenc < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cfft1f ', 4)
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) + 4) then
    ier = 2
    call xerfft ('cfft1f ', 6)
  else if (lenwrk < 2*n) then
    ier = 3
    call xerfft ('cfft1f ', 8)
  end if

  if (n == 1) then
    return
  end if

  iw1 = n+n+1
  call c1fm1f (n,inc,c,work,wsave,wsave(iw1),wsave(iw1+1))

  return
end subroutine cfft1f
subroutine cfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! CFFT1I: initialization for CFFT1B and CFFT1F.
!
!  Discussion:
!
!    CFFT1I initializes array WSAVE for use in its companion routines 
!    CFFT1B and CFFT1F.  Routine CFFT1I must be called before the first 
!    call to CFFT1B or CFFT1F, and after whenever the value of integer 
!    N changes.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product 
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors 
!    of N and  also containing certain trigonometric values which will be used 
!    in routines CFFT1B or CFFT1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.

  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) + 4) then
    ier = 2
    call xerfft ('cfftmi ', 3)
  end if

  if ( n == 1 ) then
    return
  end if

  iw1 = n+n+1

  call r4_mcfti1 (n,wsave,wsave(iw1),wsave(iw1+1))

  return
end subroutine cfft1i



subroutine cmf2kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF2KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,2)
  real ( kind = 8 ) ch(2,in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

  if ( 1 < ido .or. na == 1 ) go to 102

    do k=1,l1
      do m1=1,m1d,im1
        chold1 = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
        cc(1,m1,k,1,2) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
        cc(1,m1,k,1,1) = chold1
        chold2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
        cc(2,m1,k,1,2) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
        cc(2,m1,k,1,1) = chold2
      end do
    end do

    return

  102 continue

    do k=1,l1
      m2 = m2s
      do m1=1,m1d,im1
        m2 = m2+im2
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
        ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
        ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
      end do
    end do

      do i=2,ido
         do k=1,l1
         m2 = m2s
         do m1=1,m1d,im1
         m2 = m2+im2
            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
            tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
            ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
            ch(2,m2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
            ch(1,m2,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
         end do
       end do
  end do

  return
end subroutine cmf2kb
subroutine cmf2kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF2KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,2)
  real ( kind = 8 ) ch(2,in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

      m1d = (lot-1)*im1+1
      m2s = 1-im2

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 2 * l1, kind = 8 )
      if (na == 1) go to 106

      do k=1,l1
        do m1=1,m1d,im1
         chold1 = sn*(cc(1,m1,k,1,1)+cc(1,m1,k,1,2))
         cc(1,m1,k,1,2) = sn*(cc(1,m1,k,1,1)-cc(1,m1,k,1,2))
         cc(1,m1,k,1,1) = chold1
         chold2 = sn*(cc(2,m1,k,1,1)+cc(2,m1,k,1,2))
         cc(2,m1,k,1,2) = sn*(cc(2,m1,k,1,1)-cc(2,m1,k,1,2))
         cc(2,m1,k,1,1) = chold2
        end do
      end do

      return

  106 do 107 k=1,l1
         m2 = m2s
         do 107 m1=1,m1d,im1
         m2 = m2+im2
         ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+cc(1,m1,k,1,2))
         ch(1,m2,k,2,1) = sn*(cc(1,m1,k,1,1)-cc(1,m1,k,1,2))
         ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+cc(2,m1,k,1,2))
         ch(2,m2,k,2,1) = sn*(cc(2,m1,k,1,1)-cc(2,m1,k,1,2))
  107 continue

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
         ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
         ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
         ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
  103 continue

  do i=2,ido
    do k=1,l1
      m2 = m2s
      do m1=1,m1d,im1
        m2 = m2+im2
        ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
        tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
        ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
        ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
        ch(2,m2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
        ch(1,m2,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
      end do
    end do
  end do

  return
end subroutine cmf2kf
subroutine cmf3kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF3KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,3)
  real ( kind = 8 ) ch(2,in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ), parameter :: taui =  0.866025403784439E+00
  real ( kind = 8 ), parameter :: taur = -0.5E+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

  if ( 1 < ido .or. na == 1) go to 102

      do k=1,l1
         do m1=1,m1d,im1
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
         cr2 = cc(1,m1,k,1,1)+taur*tr2
         cc(1,m1,k,1,1) = cc(1,m1,k,1,1)+tr2
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
         ci2 = cc(2,m1,k,1,1)+taur*ti2
         cc(2,m1,k,1,1) = cc(2,m1,k,1,1)+ti2
         cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
         ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
         cc(1,m1,k,1,2) = cr2-ci3
         cc(1,m1,k,1,3) = cr2+ci3
         cc(2,m1,k,1,2) = ci2+cr3
         cc(2,m1,k,1,3) = ci2-cr3
        end do
      end do

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
         cr2 = cc(1,m1,k,1,1)+taur*tr2
         ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
         ci2 = cc(2,m1,k,1,1)+taur*ti2
         ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
         cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
         ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
         ch(1,m2,k,2,1) = cr2-ci3
         ch(1,m2,k,3,1) = cr2+ci3
         ch(2,m2,k,2,1) = ci2+cr3
         ch(2,m2,k,3,1) = ci2-cr3
  103 continue

      do 105 i=2,ido
        do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
            cr2 = cc(1,m1,k,i,1)+taur*tr2
            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
            ci2 = cc(2,m1,k,i,1)+taur*ti2
            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
            cr3 = taui*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
            ci3 = taui*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
            ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
            ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
            ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
  104    continue
  105 continue

  return
end subroutine cmf3kb
subroutine cmf3kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF3KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,3)
  real ( kind = 8 ) ch(2,in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ), parameter :: taui = -0.866025403784439E+00
  real ( kind = 8 ), parameter :: taur = -0.5E+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 3 * l1, kind = 8 )
      if (na == 1) go to 106
      do 101 k=1,l1
         do 101 m1=1,m1d,im1
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
         cr2 = cc(1,m1,k,1,1)+taur*tr2
         cc(1,m1,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
         ci2 = cc(2,m1,k,1,1)+taur*ti2
         cc(2,m1,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
         cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
         ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
         cc(1,m1,k,1,2) = sn*(cr2-ci3)
         cc(1,m1,k,1,3) = sn*(cr2+ci3)
         cc(2,m1,k,1,2) = sn*(ci2+cr3)
         cc(2,m1,k,1,3) = sn*(ci2-cr3)
  101 continue

      return

  106 do 107 k=1,l1
         m2 = m2s
         do 107 m1=1,m1d,im1
         m2 = m2+im2
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
         cr2 = cc(1,m1,k,1,1)+taur*tr2
         ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
         ci2 = cc(2,m1,k,1,1)+taur*ti2
         ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
         cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
         ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
         ch(1,m2,k,2,1) = sn*(cr2-ci3)
         ch(1,m2,k,3,1) = sn*(cr2+ci3)
         ch(2,m2,k,2,1) = sn*(ci2+cr3)
         ch(2,m2,k,3,1) = sn*(ci2-cr3)
  107 continue

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
         cr2 = cc(1,m1,k,1,1)+taur*tr2
         ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
         ci2 = cc(2,m1,k,1,1)+taur*ti2
         ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
         cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
         ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
         ch(1,m2,k,2,1) = cr2-ci3
         ch(1,m2,k,3,1) = cr2+ci3
         ch(2,m2,k,2,1) = ci2+cr3
         ch(2,m2,k,3,1) = ci2-cr3
  103 continue
      do 105 i=2,ido
        do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
            cr2 = cc(1,m1,k,i,1)+taur*tr2
            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
            ci2 = cc(2,m1,k,i,1)+taur*ti2
            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
            cr3 = taui*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
            ci3 = taui*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
  104    continue
  105 continue

  return
end subroutine cmf3kf
subroutine cmf4kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF4KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,4)
  real ( kind = 8 ) ch(2,in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

      if ( 1 < ido .or. na == 1) go to 102

      do k=1,l1
        do m1=1,m1d,im1
          ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
          ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
          tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
          ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
          tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
          tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
          ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
          tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
          cc(1,m1,k,1,1) = tr2+tr3
          cc(1,m1,k,1,3) = tr2-tr3
          cc(2,m1,k,1,1) = ti2+ti3
          cc(2,m1,k,1,3) = ti2-ti3
          cc(1,m1,k,1,2) = tr1+tr4
          cc(1,m1,k,1,4) = tr1-tr4
          cc(2,m1,k,1,2) = ti1+ti4
          cc(2,m1,k,1,4) = ti1-ti4
        end do
      end do

      return

  102 continue

      do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
         ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
         tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
         ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
         tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
         tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
         ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = tr2+tr3
         ch(1,m2,k,3,1) = tr2-tr3
         ch(2,m2,k,1,1) = ti2+ti3
         ch(2,m2,k,3,1) = ti2-ti3
         ch(1,m2,k,2,1) = tr1+tr4
         ch(1,m2,k,4,1) = tr1-tr4
         ch(2,m2,k,2,1) = ti1+ti4
         ch(2,m2,k,4,1) = ti1-ti4
  103 continue

      do 105 i=2,ido
         do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
            ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
            ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
            tr4 = cc(2,m1,k,i,4)-cc(2,m1,k,i,2)
            tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
            tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
            ti4 = cc(1,m1,k,i,2)-cc(1,m1,k,i,4)
            tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
            ch(1,m2,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,m2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,m2,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
            ch(2,m2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
            ch(1,m2,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
            ch(2,m2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
            ch(1,m2,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
            ch(2,m2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end subroutine cmf4kb 
subroutine cmf4kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF4KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,4)
  real ( kind = 8 ) ch(2,in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 4 * l1, kind = 8 )
      if (na == 1) go to 106
      do 101 k=1,l1
         do 101 m1=1,m1d,im1
         ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
         ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
         tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
         tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
         tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
         ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
         tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
         cc(1,m1,k,1,1) = sn*(tr2+tr3)
         cc(1,m1,k,1,3) = sn*(tr2-tr3)
         cc(2,m1,k,1,1) = sn*(ti2+ti3)
         cc(2,m1,k,1,3) = sn*(ti2-ti3)
         cc(1,m1,k,1,2) = sn*(tr1+tr4)
         cc(1,m1,k,1,4) = sn*(tr1-tr4)
         cc(2,m1,k,1,2) = sn*(ti1+ti4)
         cc(2,m1,k,1,4) = sn*(ti1-ti4)
  101 continue

      return

  106 do 107 k=1,l1
         m2 = m2s
         do 107 m1=1,m1d,im1
         m2 = m2+im2
         ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
         ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
         tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
         tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
         tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
         ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
         tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = sn*(tr2+tr3)
         ch(1,m2,k,3,1) = sn*(tr2-tr3)
         ch(2,m2,k,1,1) = sn*(ti2+ti3)
         ch(2,m2,k,3,1) = sn*(ti2-ti3)
         ch(1,m2,k,2,1) = sn*(tr1+tr4)
         ch(1,m2,k,4,1) = sn*(tr1-tr4)
         ch(2,m2,k,2,1) = sn*(ti1+ti4)
         ch(2,m2,k,4,1) = sn*(ti1-ti4)
  107 continue

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
         ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
         tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
         tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
         tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
         ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
         tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = tr2+tr3
         ch(1,m2,k,3,1) = tr2-tr3
         ch(2,m2,k,1,1) = ti2+ti3
         ch(2,m2,k,3,1) = ti2-ti3
         ch(1,m2,k,2,1) = tr1+tr4
         ch(1,m2,k,4,1) = tr1-tr4
         ch(2,m2,k,2,1) = ti1+ti4
         ch(2,m2,k,4,1) = ti1-ti4
  103 continue
      do 105 i=2,ido
         do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
            ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
            ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
            tr4 = cc(2,m1,k,i,2)-cc(2,m1,k,i,4)
            tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
            tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
            ti4 = cc(1,m1,k,i,4)-cc(1,m1,k,i,2)
            tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
            ch(1,m2,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,m2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,m2,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
            ch(2,m2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
            ch(1,m2,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
            ch(2,m2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
            ch(1,m2,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
            ch(2,m2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end subroutine cmf4kf
subroutine cmf5kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF5KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,5)
  real ( kind = 8 ) ch(2,in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 =  0.9510565162951536E+00
  real ( kind = 8 ), parameter :: ti12 =  0.5877852522924731E+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474E+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474E+00
  real ( kind = 8 ) wa(ido,4,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

      if ( 1 < ido .or. na == 1) go to 102
      do 101 k=1,l1
         do 101 m1=1,m1d,im1
         ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
         ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
         tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
         tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
         chold1 = cc(1,m1,k,1,1)+tr2+tr3
         chold2 = cc(2,m1,k,1,1)+ti2+ti3
         cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,m1,k,1,1) = chold1
         cc(2,m1,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,m1,k,1,2) = cr2-ci5
         cc(1,m1,k,1,5) = cr2+ci5
         cc(2,m1,k,1,2) = ci2+cr5
         cc(2,m1,k,1,3) = ci3+cr4
         cc(1,m1,k,1,3) = cr3-ci4
         cc(1,m1,k,1,4) = cr3+ci4
         cc(2,m1,k,1,4) = ci3-cr4
         cc(2,m1,k,1,5) = ci2-cr5
  101 continue

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
         ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
         tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
         tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
         ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
         cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,m2,k,2,1) = cr2-ci5
         ch(1,m2,k,5,1) = cr2+ci5
         ch(2,m2,k,2,1) = ci2+cr5
         ch(2,m2,k,3,1) = ci3+cr4
         ch(1,m2,k,3,1) = cr3-ci4
         ch(1,m2,k,4,1) = cr3+ci4
         ch(2,m2,k,4,1) = ci3-cr4
         ch(2,m2,k,5,1) = ci2-cr5
  103 continue

      do 105 i=2,ido
         do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
            ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
            ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
            tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
            tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
            tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
            cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
            ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
            ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
            ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
            ch(1,m2,k,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
            ch(2,m2,k,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
            ch(1,m2,k,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
            ch(2,m2,k,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end subroutine cmf5kb
subroutine cmf5kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*****************************************************************************80
!
!! CMF5KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(2,in1,l1,ido,5)
  real ( kind = 8 ) ch(2,in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 = -0.9510565162951536E+00
  real ( kind = 8 ), parameter :: ti12 = -0.5877852522924731E+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474E+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474E+00
  real ( kind = 8 ) wa(ido,4,2)

  m1d = (lot-1)*im1+1
  m2s = 1-im2

      if ( 1 < ido ) go to 102
      sn = 1.0E+00 / real ( 5 * l1, kind = 4 )
      if (na == 1) go to 106
      do 101 k=1,l1
         do 101 m1=1,m1d,im1
         ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
         ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
         tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
         tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
         chold1 = sn*(cc(1,m1,k,1,1)+tr2+tr3)
         chold2 = sn*(cc(2,m1,k,1,1)+ti2+ti3)
         cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,m1,k,1,1) = chold1
         cc(2,m1,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,m1,k,1,2) = sn*(cr2-ci5)
         cc(1,m1,k,1,5) = sn*(cr2+ci5)
         cc(2,m1,k,1,2) = sn*(ci2+cr5)
         cc(2,m1,k,1,3) = sn*(ci3+cr4)
         cc(1,m1,k,1,3) = sn*(cr3-ci4)
         cc(1,m1,k,1,4) = sn*(cr3+ci4)
         cc(2,m1,k,1,4) = sn*(ci3-cr4)
         cc(2,m1,k,1,5) = sn*(ci2-cr5)
  101 continue

      return

  106 do 107 k=1,l1
         m2 = m2s
         do 107 m1=1,m1d,im1
         m2 = m2+im2
         ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
         ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
         tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
         tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2+tr3)
         ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2+ti3)
         cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,m2,k,2,1) = sn*(cr2-ci5)
         ch(1,m2,k,5,1) = sn*(cr2+ci5)
         ch(2,m2,k,2,1) = sn*(ci2+cr5)
         ch(2,m2,k,3,1) = sn*(ci3+cr4)
         ch(1,m2,k,3,1) = sn*(cr3-ci4)
         ch(1,m2,k,4,1) = sn*(cr3+ci4)
         ch(2,m2,k,4,1) = sn*(ci3-cr4)
         ch(2,m2,k,5,1) = sn*(ci2-cr5)
  107 continue

      return

  102 do 103 k=1,l1
         m2 = m2s
         do 103 m1=1,m1d,im1
         m2 = m2+im2
         ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
         ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
         ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
         ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
         tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
         tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
         tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
         tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
         ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
         ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
         cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,m2,k,2,1) = cr2-ci5
         ch(1,m2,k,5,1) = cr2+ci5
         ch(2,m2,k,2,1) = ci2+cr5
         ch(2,m2,k,3,1) = ci3+cr4
         ch(1,m2,k,3,1) = cr3-ci4
         ch(1,m2,k,4,1) = cr3+ci4
         ch(2,m2,k,4,1) = ci3-cr4
         ch(2,m2,k,5,1) = ci2-cr5
  103 continue
      do 105 i=2,ido
         do 104 k=1,l1
         m2 = m2s
         do 104 m1=1,m1d,im1
         m2 = m2+im2
            ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
            ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
            ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
            tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
            tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
            tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
            cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
            ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,m2,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
            ch(2,m2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
            ch(1,m2,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
            ch(2,m2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end subroutine cmf5kf
subroutine cmfgkb ( lot, ido, ip, l1, lid, na, cc, cc1, im1, in1, &
  ch, ch1, im2, in2, wa )

!*****************************************************************************80
!
!! CMFGKB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(2,in1,l1,ip,ido)
  real ( kind = 8 ) cc1(2,in1,lid,ip)
  real ( kind = 8 ) ch(2,in2,l1,ido,ip)
  real ( kind = 8 ) ch1(2,in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

  m1d = (lot-1)*im1+1
  m2s = 1-im2
  ipp2 = ip+2
  ipph = (ip+1)/2

  do ki=1,lid
    m2 = m2s
    do m1=1,m1d,im1
      m2 = m2+im2
      ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
      ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
    end do
  end do

  do j=2,ipph
    jc = ipp2-j
    do ki=1,lid
      m2 = m2s
      do m1=1,m1d,im1
        m2 = m2+im2
        ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
        ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
      end do
    end do
  end do

  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
         m2 = m2s
         do 117 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
            cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
         m2 = m2s
         do 113 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
            cc1(1,m1,ki,lc) = wa(1,l-1,2)*ch1(1,m2,ki,ip)
            cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
            cc1(2,m1,ki,lc) = wa(1,l-1,2)*ch1(2,m2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = wa(1,idlj,2)
            do 114 ki=1,lid
               m2 = m2s
               do 114 m1=1,m1d,im1
               m2 = m2+im2
               cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
               cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
               cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
               cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
  114       continue
  115    continue
  116 continue

      if( 1 < ido .or. na == 1) go to 136

      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
         do 119 m1=1,m1d,im1
            chold1 = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
            chold2 = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
            cc1(1,m1,ki,j) = chold1
            cc1(2,m1,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
            cc1(2,m1,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
            cc1(1,m1,ki,jc) = chold2
  119    continue
  120 continue

      return

  136 do 137 ki=1,lid
         m2 = m2s
         do 137 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
  137 continue

      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
         m2 = m2s
         do 134 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
            ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
  134    continue
  135 continue

      if (ido == 1) then
        return
      end if

      do 131 i=1,ido
         do 130 k=1,l1
         m2 = m2s
         do 130 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
            cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
  130    continue
  131 continue

      do 123 j=2,ip
         do 122 k=1,l1
         m2 = m2s
         do 122 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
            cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
  122    continue
  123 continue

      do 126 j=2,ip
         do 125 i=2,ido
            do 124 k=1,l1
               m2 = m2s
               do 124 m1=1,m1d,im1
               m2 = m2+im2
               cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                            -wa(i,j-1,2)*ch(2,m2,k,i,j)
               cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                            +wa(i,j-1,2)*ch(1,m2,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end subroutine cmfgkb
subroutine cmfgkf ( lot, ido, ip, l1, lid, na, cc, cc1, im1, in1, &
  ch, ch1, im2, in2, wa )

!*****************************************************************************80
!
!! CMFGKF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(2,in1,l1,ip,ido)
  real ( kind = 8 ) cc1(2,in1,lid,ip)
  real ( kind = 8 ) ch(2,in2,l1,ido,ip)
  real ( kind = 8 ) ch1(2,in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j 
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      m1d = (lot-1)*im1+1
      m2s = 1-im2
      ipp2 = ip+2
      ipph = (ip+1)/2
      do 110 ki=1,lid
         m2 = m2s
         do 110 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
  110 continue

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
         m2 = m2s
         do 112 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
            ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
  112    continue
  111 continue
      do 118 j=2,ipph
         do 117 ki=1,lid
         m2 = m2s
         do 117 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
            cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
  117    continue
  118 continue
      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
         m2 = m2s
         do 113 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
            cc1(1,m1,ki,lc) = -wa(1,l-1,2)*ch1(1,m2,ki,ip)
            cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
            cc1(2,m1,ki,lc) = -wa(1,l-1,2)*ch1(2,m2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do 114 ki=1,lid
               m2 = m2s
               do 114 m1=1,m1d,im1
               m2 = m2+im2
               cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
               cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
               cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
               cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
  114       continue
  115    continue
  116 continue
      if ( 1 < ido ) go to 136
      sn = 1.0E+00 / real ( ip * l1, kind = 8 )
      if (na == 1) go to 146
      do 149 ki=1,lid
         m2 = m2s
         do 149 m1=1,m1d,im1
         m2 = m2+im2
         cc1(1,m1,ki,1) = sn*cc1(1,m1,ki,1)
         cc1(2,m1,ki,1) = sn*cc1(2,m1,ki,1)
  149 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
         do 119 m1=1,m1d,im1
            chold1 = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
            chold2 = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
            cc1(1,m1,ki,j) = chold1
            cc1(2,m1,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
            cc1(2,m1,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
            cc1(1,m1,ki,jc) = chold2
  119    continue
  120 continue

      return

  146 do 147 ki=1,lid
         m2 = m2s
         do 147 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = sn*cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = sn*cc1(2,m1,ki,1)
  147 continue
      do 145 j=2,ipph
         jc = ipp2-j
         do 144 ki=1,lid
         m2 = m2s
         do 144 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
            ch1(2,m2,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
            ch1(1,m2,ki,jc) = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
            ch1(2,m2,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
  144    continue
  145 continue

      return

  136 do 137 ki=1,lid
         m2 = m2s
         do 137 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
  137 continue
      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
         m2 = m2s
         do 134 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
            ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
  134    continue
  135 continue
      do 131 i=1,ido
         do 130 k=1,l1
         m2 = m2s
         do 130 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
            cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
  130    continue
  131 continue
      do 123 j=2,ip
         do 122 k=1,l1
         m2 = m2s
         do 122 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
            cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
  122    continue
  123 continue
      do 126 j=2,ip
         do 125 i=2,ido
            do 124 k=1,l1
               m2 = m2s
               do 124 m1=1,m1d,im1
               m2 = m2+im2
               cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                            +wa(i,j-1,2)*ch(2,m2,k,i,j)
               cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                            -wa(i,j-1,2)*ch(1,m2,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end subroutine cmfgkf

subroutine cosq1b ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! COSQ1B: real single precision backward cosine quarter wave transform, 1D.
!
!  Discussion:
!
!    COSQ1B computes the one-dimensional Fourier transform of a sequence 
!    which is a cosine series with odd wave numbers.  This transform is 
!    referred to as the backward transform or Fourier synthesis, transforming
!    the sequence from spectral to physical space.
!
!    This transform is normalized since a call to COSQ1B followed
!    by a call to COSQ1F (or vice-versa) reproduces the original
!    array  within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements to be transformed 
!    in the sequence.  The transform is most efficient when N is a 
!    product of small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR); on input, containing the sequence 
!    to be transformed, and on output, containing the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSQ1I before the first call to routine COSQ1F 
!    or COSQ1B for a given transform length N.  WSAVE's contents may be 
!    re-used for subsequent calls to COSQ1F and COSQ1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) ssqrt2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cosq1b', 6)
    return
  else if (lensav < &
    2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosq1b', 8)
    return
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('cosq1b', 10)
    return
  end if

      if (n-2) 300,102,103

 102  ssqrt2 = 1.0E+00 / sqrt ( 2.0E+00 )
      x1 = x(1,1)+x(1,2)
      x(1,2) = ssqrt2*(x(1,1)-x(1,2))
      x(1,1) = x1
      return

  103 call cosqb1 (n,inc,x,wsave,work,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('cosq1b',-5)
      end if

  300 continue

  return
end subroutine cosq1b
subroutine cosq1f ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! COSQ1F: real single precision forward cosine quarter wave transform, 1D.
!
!  Discussion:
!
!    COSQ1F computes the one-dimensional Fourier transform of a sequence 
!    which is a cosine series with odd wave numbers.  This transform is 
!    referred to as the forward transform or Fourier analysis, transforming 
!    the sequence from physical to spectral space.
!
!    This transform is normalized since a call to COSQ1F followed
!    by a call to COSQ1B (or vice-versa) reproduces the original
!    array  within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements to be transformed 
!    in the sequence.  The transform is most efficient when N is a 
!    product of small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR); on input, containing the sequence 
!    to be transformed, and on output, containing the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSQ1I before the first call to routine COSQ1F 
!    or COSQ1B for a given transform length N.  WSAVE's contents may be 
!    re-used for subsequent calls to COSQ1F and COSQ1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) lenx
  real ( kind = 8 ) ssqrt2
  real ( kind = 8 ) tsqx
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cosq1f', 6)
    return
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosq1f', 8)
    return
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('cosq1f', 10)
    return
  end if

      if (n-2) 102,101,103
  101 ssqrt2 = 1.0E+00 / sqrt ( 2.0E+00 )
      tsqx = ssqrt2*x(1,2)
      x(1,2) = 0.5E+00 *x(1,1)-tsqx
      x(1,1) = 0.5E+00 *x(1,1)+tsqx
  102 return
  103 call cosqf1 (n,inc,x,wsave,work,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cosq1f',-5)
  end if

  return
end subroutine cosq1f
subroutine cosq1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! COSQ1I: initialization for COSQ1B and COSQ1F.
!
!  Discussion:
!
!    COSQ1I initializes array WSAVE for use in its companion routines 
!    COSQ1F and COSQ1B.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product 
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of N
!    and also containing certain trigonometric values which will be used 
!    in routines COSQ1B or COSQ1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  real ( kind = 8 ) fk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  real ( kind = 8 ) pih
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosq1i', 3)
    return
  end if

!  pih = 2.0E+00 * atan ( 1.0E+00 )
   pih = 3.141592653589793d0/dble(2)
  dt = pih / real  ( n, kind = 8 )
  fk = 0.0E+00

  do k=1,n
    fk = fk + 1.0E+00
    wsave(k) = cos(fk*dt)
  end do

  lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4
  call rfft1i (n, wsave(n+1), lnsv, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cosq1i',-5)
  end if

  return
end subroutine cosq1i
subroutine cosqb1 ( n, inc, x, wsave, work, ier )

!*****************************************************************************80
!
!! COSQB1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xim1

  ier = 0
  ns2 = (n+1)/2
  np2 = n+2

  do i=3,n,2
    xim1 = x(1,i-1)+x(1,i)
    x(1,i) = 0.5E+00 * (x(1,i-1)-x(1,i))
    x(1,i-1) = 0.5E+00 * xim1
  end do

  x(1,1) = 0.5E+00 * x(1,1)
  modn = mod(n,2)

  if (modn == 0 ) then
    x(1,n) = 0.5E+00 * x(1,n)
  end if

  lenx = inc*(n-1)  + 1
  lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) + 4
  lnwk = n

  call rfft1b(n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cosqb1',-5)
    return
  end if

  do k=2,ns2
    kc = np2-k
    work(k) = wsave(k-1)*x(1,kc)+wsave(kc-1)*x(1,k)
    work(kc) = wsave(k-1)*x(1,k)-wsave(kc-1)*x(1,kc)
  end do

  if (modn == 0) then
    x(1,ns2+1) = wsave(ns2)*(x(1,ns2+1)+x(1,ns2+1))
  end if

  do k=2,ns2
    kc = np2-k
    x(1,k) = work(k)+work(kc)
    x(1,kc) = work(k)-work(kc)
  end do

  x(1,1) = x(1,1)+x(1,1)

  return
end subroutine cosqb1
subroutine cosqf1 ( n, inc, x, wsave, work, ier )

!*****************************************************************************80
!
!! COSQF1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xim1

  ier = 0
  ns2 = (n+1)/2
  np2 = n+2

  do k=2,ns2
    kc = np2-k
    work(k)  = x(1,k)+x(1,kc)
    work(kc) = x(1,k)-x(1,kc)
  end do

  modn = mod(n,2)

  if (modn == 0) then
    work(ns2+1) = x(1,ns2+1)+x(1,ns2+1)
  end if

  do k=2,ns2
    kc = np2-k
    x(1,k)  = wsave(k-1)*work(kc)+wsave(kc-1)*work(k)
    x(1,kc) = wsave(k-1)*work(k) -wsave(kc-1)*work(kc)
  end do

  if (modn == 0) then
    x(1,ns2+1) = wsave(ns2)*work(ns2+1)
  end if

  lenx = inc*(n-1)  + 1
  lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) + 4
  lnwk = n

  call rfft1f(n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cosqf1',-5)
    return
  end if

  do i=3,n,2
    xim1 = 0.5E+00 * (x(1,i-1)+x(1,i))
    x(1,i) = 0.5E+00 * (x(1,i-1)-x(1,i))
    x(1,i-1) = xim1
  end do

  return
end subroutine cosqf1
subroutine cosqmb ( lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, &
  ier )

!*****************************************************************************80
!
!! COSQMB: real single precision backward cosine quarter wave, multiple vectors.
!
!  Discussion:
!
!    COSQMB computes the one-dimensional Fourier transform of multiple
!    sequences, each of which is a cosine series with odd wave numbers.
!    This transform is referred to as the backward transform or Fourier
!    synthesis, transforming the sequences from spectral to physical space.
!
!    This transform is normalized since a call to COSQMB followed
!    by a call to COSQMF (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, 
!    in array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), array containing LOT sequences, 
!    each having length N.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSQMI before the first call to routine COSQMF 
!    or COSQMB for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to COSQMF and COSQMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) ssqrt2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1
!  logical xercon

  ier = 0

  if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cosqmb', 6)
    return
  else if (lensav < &
    2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosqmb', 8)
    return
  else if (lenwrk < lot*n) then
    ier = 3
    call xerfft ('cosqmb', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('cosqmb', -1)
    return
  end if

      lj = (lot-1)*jump+1
      if (n-2) 101,102,103
 101  do m=1,lj,jump
        x(m,1) = x(m,1)
      end do
      return
 102  ssqrt2 = 1.0E+00 / sqrt ( 2.0E+00 )
      do m=1,lj,jump
        x1 = x(m,1)+x(m,2)
        x(m,2) = ssqrt2*(x(m,1)-x(m,2))
        x(m,1) = x1
      end do
      return

  103 call mcsqb1 (lot,jump,n,inc,x,wsave,work,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('cosqmb',-5)
      end if

  return
end subroutine cosqmb 
subroutine cosqmf ( lot, jump, n, inc, x, lenx, wsave, lensav, work, &
  lenwrk, ier )

!*****************************************************************************80
!
!! COSQMF: real single precision forward cosine quarter wave, multiple vectors.
!
!  Discussion:
!
!    COSQMF computes the one-dimensional Fourier transform of multiple 
!    sequences within a real array, where each of the sequences is a 
!    cosine series with odd wave numbers.  This transform is referred to 
!    as the forward transform or Fourier synthesis, transforming the 
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to COSQMF followed
!    by a call to COSQMB (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in 
!    array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), array containing LOT sequences, 
!    each having length N.  R can have any number of dimensions, but the total
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSQMI before the first call to routine COSQMF 
!    or COSQMB for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to COSQMF and COSQMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) ssqrt2
  real ( kind = 8 ) tsqx
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon

  ier = 0

  if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cosqmf', 6)
    return
  else if (lensav < &
    2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosqmf', 8)
    return
  else if (lenwrk < lot*n) then
    ier = 3
    call xerfft ('cosqmf', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('cosqmf', -1)
    return
  end if

  lj = (lot-1)*jump+1

  if (n-2) 102,101,103
  101 ssqrt2 = 1.0E+00 / sqrt ( 2.0E+00 )

      do m=1,lj,jump
        tsqx = ssqrt2*x(m,2)
        x(m,2) = 0.5E+00 * x(m,1)-tsqx
        x(m,1) = 0.5E+00 * x(m,1)+tsqx
      end do

  102 return

  103 call mcsqf1 (lot,jump,n,inc,x,wsave,work,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('cosqmf',-5)
      end if

  return
end subroutine cosqmf
subroutine cosqmi ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! COSQMI: initialization for COSQMB and COSQMF.
!
!  Discussion:
!
!    COSQMI initializes array WSAVE for use in its companion routines 
!    COSQMF and COSQMB.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of 
!    N and also containing certain trigonometric values which will be used 
!    in routines COSQMB or COSQMF.
!
!    Input, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  real ( kind = 8 ) fk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  real ( kind = 8 ) pih
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cosqmi', 3)
    return
  end if

 ! pih = 2.0E+00 * atan ( 1.0E+00 )
   pih = 3.141592653589793d0/dble(2)
  dt = pih/real ( n, kind = 8 )
  fk = 0.0E+00

  do k=1,n
    fk = fk + 1.0E+00
    wsave(k) = cos(fk*dt)
  end do

  lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4

  call rfftmi (n, wsave(n+1), lnsv, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cosqmi',-5)
  end if

  return
end subroutine cosqmi
subroutine cost1b ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! COST1B: real single precision backward cosine transform, 1D.
!
!  Discussion:
!
!    COST1B computes the one-dimensional Fourier transform of an even 
!    sequence within a real array.  This transform is referred to as 
!    the backward transform or Fourier synthesis, transforming the sequence 
!    from spectral to physical space.
!
!    This transform is normalized since a call to COST1B followed
!    by a call to COST1F (or vice-versa) reproduces the original array 
!    within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N-1 is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing the sequence to 
!     be transformed.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COST1I before the first call to routine COST1F 
!    or COST1B for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to COST1F and COST1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N-1.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cost1b', 6)
    return
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cost1b', 8)
    return
  else if (lenwrk < n-1) then
    ier = 3
    call xerfft ('cost1b', 10)
    return
  end if

  if (n == 1) then
    return
  end if

  call costb1 (n,inc,x,wsave,work,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cost1b',-5)
  end if

  return
end subroutine cost1b
subroutine cost1f ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! COST1F: real single precision forward cosine transform, 1D.
!
!  Discussion:
!
!    COST1F computes the one-dimensional Fourier transform of an even 
!    sequence within a real array.  This transform is referred to as the 
!    forward transform or Fourier analysis, transforming the sequence 
!    from  physical to spectral space.
!
!    This transform is normalized since a call to COST1F followed by a call 
!    to COST1B (or vice-versa) reproduces the original array within 
!    roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N-1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing the sequence to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COST1I before the first call to routine COST1F
!    or COST1B for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to COST1F and COST1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N-1.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cost1f', 6)
    return
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cost1f', 8)
    return
  else if (lenwrk < n-1) then
    ier = 3
    call xerfft ('cost1f', 10)
    return
  end if

  if (n == 1) then
    return
  end if

  call costf1(n,inc,x,wsave,work,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cost1f',-5)
  end if

  return
end subroutine cost1f
subroutine cost1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! COST1I: initialization for COST1B and COST1F.
!
!  Discussion:
!
!    COST1I initializes array WSAVE for use in its companion routines 
!    COST1F and COST1B.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N-1 is a product 
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, dimension of WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of
!    N and also containing certain trigonometric values which will be used in
!    routines COST1B or COST1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  real ( kind = 8 ) fk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) pi
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('cost1i', 3)
    return
  end if

  if ( n <= 3 ) then
    return
  end if

  nm1 = n-1
  np1 = n+1
  ns2 = n/2
!  pi = 4.0E+00 * atan ( 1.0E+00 )
   pi= 3.141592653589793d0 
  dt = pi/ real ( nm1, kind = 8 )
  fk = 0.0E+00

  do k=2,ns2
    kc = np1-k
    fk = fk + 1.0E+00
    wsave(k) = 2.0E+00 * sin(fk*dt)
    wsave(kc) = 2.0E+00 * cos(fk*dt)
  end do

  lnsv = nm1 + int(log( real ( nm1, kind = 8 ) )/log( 2.0E+00 )) +4

  call rfft1i (nm1, wsave(n+1), lnsv, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cost1i',-5)
  end if

  return
end subroutine cost1i 
subroutine costb1 ( n, inc, x, wsave, work, ier )

!*****************************************************************************80
!
!! COSTB1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum
  real ( kind = 8 ) fnm1s2
  real ( kind = 8 ) fnm1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) x2
  real ( kind = 8 ) xi

  ier = 0

      nm1 = n-1
      np1 = n+1
      ns2 = n/2

      if (n-2) 106,101,102
  101 x1h = x(1,1)+x(1,2)
      x(1,2) = x(1,1)-x(1,2)
      x(1,1) = x1h
      return
  102 if ( 3 < n ) go to 103
      x1p3 = x(1,1)+x(1,3)
      x2 = x(1,2)
      x(1,2) = x(1,1)-x(1,3)
      x(1,1) = x1p3+x2
      x(1,3) = x1p3-x2
      return
  103 x(1,1) = x(1,1)+x(1,1)
      x(1,n) = x(1,n)+x(1,n)
      dsum = x(1,1)-x(1,n)
      x(1,1) = x(1,1)+x(1,n)

      do k=2,ns2
         kc = np1-k
         t1 = x(1,k)+x(1,kc)
         t2 = x(1,k)-x(1,kc)
         dsum = dsum+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(1,k) = t1-t2
         x(1,kc) = t1+t2
      end do

      modn = mod(n,2)
      if (modn == 0) go to 124
      x(1,ns2+1) = x(1,ns2+1)+x(1,ns2+1)
  124 lenx = inc*(nm1-1)  + 1
      lnsv = nm1 + int(log( real ( nm1, kind = 8 ) )/log( 2.0E+00 )) + 4
      lnwk = nm1

      call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('costb1',-5)
        return
      end if

      fnm1s2 = real ( nm1, kind = 8 ) / 2.0E+00 
      dsum = 0.5D+00 * dsum
      x(1,1) = fnm1s2*x(1,1)
      if(mod(nm1,2) /= 0) go to 30
      x(1,nm1) = x(1,nm1)+x(1,nm1)
   30 fnm1s4 = real ( nm1, kind = 8 ) / 4.0E+00

      do i=3,n,2
         xi = fnm1s4*x(1,i)
         x(1,i) = fnm1s4*x(1,i-1)
         x(1,i-1) = dsum
         dsum = dsum+xi
      end do

      if (modn /= 0) return
         x(1,n) = dsum
  106 continue

  return
end subroutine costb1
subroutine costf1 ( n, inc, x, wsave, work, ier )

!*****************************************************************************80
!
!! COSTF1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) snm1
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tx2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) xi

  ier = 0

  nm1 = n-1
  np1 = n+1
  ns2 = n/2

      if (n-2) 200,101,102
  101 x1h = x(1,1)+x(1,2)
      x(1,2) = 0.5E+00 * (x(1,1)-x(1,2))
      x(1,1) = 0.5E+00 * x1h
      go to 200
  102 if ( 3 < n ) go to 103
      x1p3 = x(1,1)+x(1,3)
      tx2 = x(1,2)+x(1,2)
      x(1,2) = 0.5E+00 * (x(1,1)-x(1,3))
      x(1,1) = 0.25E+00 *(x1p3+tx2)
      x(1,3) = 0.25E+00 *(x1p3-tx2)
      go to 200
  103 dsum = x(1,1)-x(1,n)
      x(1,1) = x(1,1)+x(1,n)
      do k=2,ns2
         kc = np1-k
         t1 = x(1,k)+x(1,kc)
         t2 = x(1,k)-x(1,kc)
         dsum = dsum+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(1,k) = t1-t2
         x(1,kc) = t1+t2
      end do
      modn = mod(n,2)
      if (modn == 0) go to 124
      x(1,ns2+1) = x(1,ns2+1)+x(1,ns2+1)
  124 lenx = inc*(nm1-1)  + 1
      lnsv = nm1 + int(log( real ( nm1, kind = 8 ) )/log( 2.0E+00 )) + 4
      lnwk = nm1

      call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('costf1',-5)
        go to 200
      end if

      snm1 = 1.0E+00 / real ( nm1, kind = 8 )
      dsum = snm1*dsum
      if(mod(nm1,2) /= 0) go to 30
      x(1,nm1) = x(1,nm1)+x(1,nm1)
   30 do i=3,n,2
         xi = 0.5E+00 * x(1,i)
         x(1,i) = 0.5E+00 * x(1,i-1)
         x(1,i-1) = dsum
         dsum = dsum+xi
      end do
      if (modn /= 0) go to 117
      x(1,n) = dsum
  117 x(1,1) = 0.5E+00 * x(1,1)
      x(1,n) = 0.5E+00 * x(1,n)
  200 continue

  return
end subroutine costf1
subroutine costmb ( lot, jump, n, inc, x, lenx, wsave, lensav, work, &
  lenwrk, ier )

!*****************************************************************************80
!
!! COSTMB: real single precision backward cosine transform, multiple vectors.
!
!  Discussion:
!
!    COSTMB computes the one-dimensional Fourier transform of multiple 
!    even sequences within a real array.  This transform is referred to 
!    as the backward transform or Fourier synthesis, transforming the 
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to COSTMB followed
!    by a call to COSTMF (or vice-versa) reproduces the original
!    array  within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in 
!    array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N-1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), array containing LOT sequences, 
!    each having length N.  On input, the data to be transformed; on output,
!    the transormed data.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSTMI before the first call to routine COSTMF 
!    or COSTMB for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to COSTMF and COSTMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*(N+1).
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon

  ier = 0

  if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('costmb', 6)
    return
  else if (lensav < &
    2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('costmb', 8)
    return
  else if (lenwrk < lot*(n+1)) then
    ier = 3
    call xerfft ('costmb', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('costmb', -1)
    return
  end if

  iw1 = lot+lot+1
  call mcstb1(lot,jump,n,inc,x,wsave,work,work(iw1),ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('costmb',-5)
  end if

  return
end subroutine costmb
subroutine costmf ( lot, jump, n, inc, x, lenx, wsave, lensav, work, &
  lenwrk, ier )

!*****************************************************************************80
!
!! COSTMF: real single precision forward cosine transform, multiple vectors.
!
!  Discussion:
!
!    COSTMF computes the one-dimensional Fourier transform of multiple 
!    even sequences within a real array.  This transform is referred to 
!    as the forward transform or Fourier analysis, transforming the 
!    sequences from physical to spectral space.
!
!    This transform is normalized since a call to COSTMF followed
!    by a call to COSTMB (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, 
!    in array R, of the first elements of two consecutive sequences to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N-1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), array containing LOT sequences, 
!    each having length N.  On input, the data to be transformed; on output,
!    the transormed data.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the  R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to COSTMI before the first call to routine COSTMF
!    or COSTMB for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to COSTMF and COSTMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*(N+1).
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon

  ier = 0

  if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('costmf', 6)
    return
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('costmf', 8)
    return
  else if (lenwrk < lot*(n+1)) then
    ier = 3
    call xerfft ('costmf', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('costmf', -1)
    return
  end if

  iw1 = lot+lot+1

  call mcstf1(lot,jump,n,inc,x,wsave,work,work(iw1),ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('costmf',-5)
  end if

  return
end subroutine costmf
subroutine costmi ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! COSTMI: initialization for COSTMB and COSTMF.
!
!  Discussion:
!
!    COSTMI initializes array WSAVE for use in its companion routines 
!    COSTMF and COSTMB.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of N 
!    and also containing certain trigonometric values which will be used 
!    in routines COSTMB or COSTMF.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  real ( kind = 8 ) fk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) pi
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('costmi', 3)
    return
  end if

  if (n <= 3) then
    return
  end if

      nm1 = n-1
      np1 = n+1
      ns2 = n/2
   !   pi = 4.0E+00 * atan ( 1.0E+00 )
       pi= 3.141592653589793d0
      dt = pi/ real ( nm1, kind = 8 )
      fk = 0.0E+00

      do k=2,ns2
         kc = np1-k
         fk = fk + 1.0E+00
         wsave(k) = 2.0E+00 * sin(fk*dt)
         wsave(kc) = 2.0E+00 * cos(fk*dt)
      end do

      lnsv = nm1 + int(log( real ( nm1, kind = 8 ) )/log( 2.0E+00 )) +4

      call rfftmi (nm1, wsave(n+1), lnsv, ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('costmi',-5)
      end if

  return
end subroutine costmi
subroutine mcsqb1 (lot,jump,n,inc,x,wsave,work,ier)

!*****************************************************************************80
!
!! MCSQB1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lot

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lot,*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xim1

  ier = 0
  lj = (lot-1)*jump+1
  ns2 = (n+1)/2
  np2 = n+2

  do i=3,n,2
    do m=1,lj,jump
      xim1 = x(m,i-1)+x(m,i)
      x(m,i) = 0.5E+00 * (x(m,i-1)-x(m,i))
      x(m,i-1) = 0.5E+00 * xim1
    end do
  end do

  do m=1,lj,jump
    x(m,1) = 0.5E+00 * x(m,1)
  end do

  modn = mod(n,2)
  if (modn == 0) then
    do m=1,lj,jump
      x(m,n) = 0.5E+00 * x(m,n)
    end do
  end if

  lenx = (lot-1)*jump + inc*(n-1)  + 1
  lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) + 4
  lnwk = lot*n

  call rfftmb(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('mcsqb1',-5)
    return
  end if

  do k=2,ns2
    kc = np2-k
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      work(m1,k) = wsave(k-1)*x(m,kc)+wsave(kc-1)*x(m,k)
      work(m1,kc) = wsave(k-1)*x(m,k)-wsave(kc-1)*x(m,kc)
    end do
  end do

  if (modn == 0) then
    do m=1,lj,jump
      x(m,ns2+1) = wsave(ns2)*(x(m,ns2+1)+x(m,ns2+1))
    end do
  end if

  do k=2,ns2
    kc = np2-k
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      x(m,k) = work(m1,k)+work(m1,kc)
      x(m,kc) = work(m1,k)-work(m1,kc)
    end do
  end do

  do m=1,lj,jump
    x(m,1) = x(m,1)+x(m,1)
  end do

  return
end subroutine mcsqb1 
subroutine mcsqf1 (lot,jump,n,inc,x,wsave,work,ier)

!*****************************************************************************80
!
!! MCSQF1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lot

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lot,*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xim1

  ier = 0

  lj = (lot-1)*jump+1
  ns2 = (n+1)/2
  np2 = n+2

  do k=2,ns2
    kc = np2-k
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      work(m1,k)  = x(m,k)+x(m,kc)
      work(m1,kc) = x(m,k)-x(m,kc)
    end do
  end do

  modn = mod(n,2)

  if (modn == 0) then
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      work(m1,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
    end do
  end if

       do 102 k=2,ns2
         kc = np2-k
         m1 = 0
         do 302 m=1,lj,jump
         m1 = m1 + 1
         x(m,k)  = wsave(k-1)*work(m1,kc)+wsave(kc-1)*work(m1,k)
         x(m,kc) = wsave(k-1)*work(m1,k) -wsave(kc-1)*work(m1,kc)
 302     continue
  102 continue

      if (modn /= 0) go to 303
      m1 = 0
      do 304 m=1,lj,jump
         m1 = m1 + 1
         x(m,ns2+1) = wsave(ns2)*work(m1,ns2+1)
 304  continue
 303  continue

      lenx = (lot-1)*jump + inc*(n-1)  + 1
      lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) + 4
      lnwk = lot*n

      call rfftmf(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('mcsqf1',-5)
        go to 400
      end if

      do 103 i=3,n,2
         do 203 m=1,lj,jump
            xim1 = 0.5E+00 * (x(m,i-1)+x(m,i))
            x(m,i) = 0.5E+00 * (x(m,i-1)-x(m,i))
            x(m,i-1) = xim1
 203     continue
  103 continue

  400 continue

  return
end subroutine mcsqf1 
subroutine mcstb1(lot,jump,n,inc,x,wsave,dsum,work,ier)

!*****************************************************************************80
!
!! MCSTB1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum(*)
  real ( kind = 8 ) fnm1s2
  real ( kind = 8 ) fnm1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) x2
  real ( kind = 8 ) xi

  ier = 0

      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      lj = (lot-1)*jump+1
      if (n-2) 106,101,102

  101 do 111 m=1,lj,jump
         x1h = x(m,1)+x(m,2)
         x(m,2) = x(m,1)-x(m,2)
         x(m,1) = x1h
  111 continue

      return

  102 if ( 3 < n ) go to 103
      do 112 m=1,lj,jump
         x1p3 = x(m,1)+x(m,3)
         x2 = x(m,2)
         x(m,2) = x(m,1)-x(m,3)
         x(m,1) = x1p3+x2
         x(m,3) = x1p3-x2
  112 continue

      return

 103  do m=1,lj,jump
        x(m,1) = x(m,1)+x(m,1)
        x(m,n) = x(m,n)+x(m,n)
      end do

      m1 = 0

      do m=1,lj,jump
         m1 = m1+1
         dsum(m1) = x(m,1)-x(m,n)
         x(m,1) = x(m,1)+x(m,n)
      end do

      do 104 k=2,ns2
         m1 = 0
         do 114 m=1,lj,jump
           m1 = m1+1
           kc = np1-k
           t1 = x(m,k)+x(m,kc)
           t2 = x(m,k)-x(m,kc)
           dsum(m1) = dsum(m1)+wsave(kc)*t2
           t2 = wsave(k)*t2
           x(m,k) = t1-t2
           x(m,kc) = t1+t2
  114    continue
  104 continue

      modn = mod(n,2)
      if (modn == 0) go to 124
         do 123 m=1,lj,jump
         x(m,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
  123    continue
 124  continue
      lenx = (lot-1)*jump + inc*(nm1-1)  + 1
      lnsv = nm1 + int(log( real ( nm1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = lot*nm1

      call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('mcstb1',-5)
        go to 106
      end if

      fnm1s2 = real ( nm1, kind = 8 ) /  2.0E+00
      m1 = 0
      do 10 m=1,lj,jump
      m1 = m1+1
      dsum(m1) = 0.5D+00 * dsum(m1)
      x(m,1) = fnm1s2 * x(m,1)
   10 continue
      if(mod(nm1,2) /= 0) go to 30
      do 20 m=1,lj,jump
      x(m,nm1) = x(m,nm1)+x(m,nm1)
   20 continue
 30   fnm1s4 = real ( nm1, kind = 8 ) / 4.0E+00
      do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            xi = fnm1s4*x(m,i)
            x(m,i) = fnm1s4*x(m,i-1)
            x(m,i-1) = dsum(m1)
            dsum(m1) = dsum(m1)+xi
  115 continue
  105 continue
      if (modn /= 0) return
      m1 = 0
      do m=1,lj,jump
         m1 = m1+1
         x(m,n) = dsum(m1)
      end do
  106 continue

  return
end subroutine mcstb1
subroutine mcstf1(lot,jump,n,inc,x,wsave,dsum,work,ier)

!*****************************************************************************80
!
!! MCSTF1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) snm1
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tx2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) xi

  ier = 0

  nm1 = n-1
  np1 = n+1
  ns2 = n/2
  lj = (lot-1)*jump+1

      if (n-2) 200,101,102
  101 do 111 m=1,lj,jump
         x1h = x(m,1)+x(m,2)
         x(m,2) = 0.5E+00 * (x(m,1)-x(m,2))
         x(m,1) = 0.5E+00 * x1h
  111 continue
      go to 200
  102 if ( 3 < n ) go to 103
      do 112 m=1,lj,jump
         x1p3 = x(m,1)+x(m,3)
         tx2 = x(m,2)+x(m,2)
         x(m,2) = 0.5E+00 * (x(m,1)-x(m,3))
         x(m,1) = 0.25E+00 * (x1p3+tx2)
         x(m,3) = 0.25E+00 * (x1p3-tx2)
  112 continue
      go to 200
  103 m1 = 0
      do 113 m=1,lj,jump
         m1 = m1+1
         dsum(m1) = x(m,1)-x(m,n)
         x(m,1) = x(m,1)+x(m,n)
  113 continue
      do 104 k=2,ns2
         m1 = 0
         do 114 m=1,lj,jump
         m1 = m1+1
         kc = np1-k
         t1 = x(m,k)+x(m,kc)
         t2 = x(m,k)-x(m,kc)
         dsum(m1) = dsum(m1)+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(m,k) = t1-t2
         x(m,kc) = t1+t2
  114    continue
  104 continue
      modn = mod(n,2)
      if (modn == 0) go to 124
         do 123 m=1,lj,jump
         x(m,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
  123    continue
 124  continue
      lenx = (lot-1)*jump + inc*(nm1-1)  + 1
      lnsv = nm1 + int(log( real ( nm1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = lot*nm1

      call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('mcstf1',-5)
        return
      end if

      snm1 = 1.0E+00 / real ( nm1, kind = 8 )
      do 10 m=1,lot
      dsum(m) = snm1*dsum(m)
   10 continue
      if(mod(nm1,2) /= 0) go to 30
      do 20 m=1,lj,jump
      x(m,nm1) = x(m,nm1)+x(m,nm1)
   20 continue
 30   do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            xi = 0.5E+00 * x(m,i)
            x(m,i) = 0.5E+00 * x(m,i-1)
            x(m,i-1) = dsum(m1)
            dsum(m1) = dsum(m1)+xi
  115 continue
  105 continue
      if (modn /= 0) go to 117

      m1 = 0
      do m=1,lj,jump
         m1 = m1+1
         x(m,n) = dsum(m1)
      end do

 117  continue

      do m=1,lj,jump
        x(m,1) = 0.5E+00 * x(m,1)
        x(m,n) = 0.5E+00 * x(m,n)
      end do

 200  continue

  return
end subroutine mcstf1
subroutine mradb2 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

!*****************************************************************************80
!
!! MRADB2 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,2,l1)
  real ( kind = 8 ) ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) wa1(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2

  do k=1,l1
    m2 = m2s
    do m1=1,m1d,im1
      m2 = m2+im2
      ch(m2,1,k,1) = cc(m1,1,1,k)+cc(m1,ido,2,k)
      ch(m2,1,k,2) = cc(m1,1,1,k)-cc(m1,ido,2,k)
    end do
  end do

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
               m2 = m2s
               do 1002 m1=1,m1d,im1
               m2 = m2+im2
        ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+cc(m1,ic-1,2,k)
        ch(m2,i,k,1) = cc(m1,i,1,k)-cc(m1,ic,2,k)
        ch(m2,i-1,k,2) = wa1(i-2)*(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k)) &
        -wa1(i-1)*(cc(m1,i,1,k)+cc(m1,ic,2,k))

        ch(m2,i,k,2) = wa1(i-2)*(cc(m1,i,1,k)+cc(m1,ic,2,k))+wa1(i-1) &
        *(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k))

 1002          continue
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 do 106 k=1,l1
          m2 = m2s
          do 1003 m1=1,m1d,im1
          m2 = m2+im2
         ch(m2,ido,k,1) = cc(m1,ido,1,k)+cc(m1,ido,1,k)
         ch(m2,ido,k,2) = -(cc(m1,1,2,k)+cc(m1,1,2,k))
 1003     continue
  106 continue
  107 continue

  return
end subroutine mradb2
subroutine mradb3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

!*****************************************************************************80
!
!! MRADB3 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,3,l1)
  real ( kind = 8 ) ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
 ! arg= 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  arg= 6.28318530717959d0 / dble(3.0)
  taur=cos(arg)
  taui=sin(arg)

  do k=1,l1
    m2 = m2s
    do m1=1,m1d,im1
      m2 = m2+im2
      ch(m2,1,k,1) = cc(m1,1,1,k)+ 2.0E+00 *cc(m1,ido,2,k)
      ch(m2,1,k,2) = cc(m1,1,1,k)+( 2.0E+00 *taur)*cc(m1,ido,2,k) &
        -( 2.0E+00 *taui)*cc(m1,1,3,k)
      ch(m2,1,k,3) = cc(m1,1,1,k)+( 2.0E+00 *taur)*cc(m1,ido,2,k) &
        + 2.0E+00 *taui*cc(m1,1,3,k)
    end do
  end do

  if (ido == 1) then
    return
  end if

  idp2 = ido+2

      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
               m2 = m2s
               do 1002 m1=1,m1d,im1
               m2 = m2+im2
        ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
        ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k))

        ch(m2,i-1,k,2) = wa1(i-2)* &
       ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
       (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
                         -wa1(i-1)* &
       ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
       (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

            ch(m2,i,k,2) = wa1(i-2)* &
       ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
       (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                        +wa1(i-1)* &
       ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
       (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

              ch(m2,i-1,k,3) = wa2(i-2)* &
       ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
       (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
         -wa2(i-1)* &
       ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
       (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

            ch(m2,i,k,3) = wa2(i-2)* &
       ((cc(m1,i,1,k)+taur*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
       (taui*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                       +wa2(i-1)* &
       ((cc(m1,i-1,1,k)+taur*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
       (taui*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

 1002          continue
  102    continue
  103 continue

  return
end subroutine mradb3
subroutine mradb4 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! MRADB4 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,4,l1)
  real ( kind = 8 ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) sqrt2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

      m1d = (m-1)*im1+1
      m2s = 1-im2
      sqrt2=sqrt( 2.0E+00 )
      do 101 k=1,l1
          m2 = m2s
          do m1=1,m1d,im1
          m2 = m2+im2
         ch(m2,1,k,3) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
         -(cc(m1,ido,2,k)+cc(m1,ido,2,k))
         ch(m2,1,k,1) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
         +(cc(m1,ido,2,k)+cc(m1,ido,2,k))
         ch(m2,1,k,4) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
         +(cc(m1,1,3,k)+cc(m1,1,3,k))
         ch(m2,1,k,2) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
         -(cc(m1,1,3,k)+cc(m1,1,3,k))
        end do
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
               m2 = m2s
               do 1002 m1=1,m1d,im1
               m2 = m2+im2
        ch(m2,i-1,k,1) = (cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
        +(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
        ch(m2,i,k,1) = (cc(m1,i,1,k)-cc(m1,ic,4,k)) &
        +(cc(m1,i,3,k)-cc(m1,ic,2,k))
        ch(m2,i-1,k,2)=wa1(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
        -(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa1(i-1) &
        *((cc(m1,i,1,k)+cc(m1,ic,4,k))+(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
        ch(m2,i,k,2)=wa1(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
        +(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa1(i-1) &
        *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))-(cc(m1,i,3,k)+cc(m1,ic,2,k)))
        ch(m2,i-1,k,3)=wa2(i-2)*((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
        -(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))-wa2(i-1) &
        *((cc(m1,i,1,k)-cc(m1,ic,4,k))-(cc(m1,i,3,k)-cc(m1,ic,2,k)))
        ch(m2,i,k,3)=wa2(i-2)*((cc(m1,i,1,k)-cc(m1,ic,4,k)) &
        -(cc(m1,i,3,k)-cc(m1,ic,2,k)))+wa2(i-1) &
        *((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k))-(cc(m1,i-1,3,k) &
        +cc(m1,ic-1,2,k)))
        ch(m2,i-1,k,4)=wa3(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
        +(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa3(i-1) &
       *((cc(m1,i,1,k)+cc(m1,ic,4,k))-(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
        ch(m2,i,k,4)=wa3(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
        -(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa3(i-1) &
        *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))+(cc(m1,i,3,k)+cc(m1,ic,2,k)))
 1002          continue
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
               m2 = m2s
               do 1003 m1=1,m1d,im1
               m2 = m2+im2
         ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
         +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
         ch(m2,ido,k,2) = sqrt2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
         -(cc(m1,1,2,k)+cc(m1,1,4,k)))
         ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
         +(cc(m1,1,4,k)-cc(m1,1,2,k))
         ch(m2,ido,k,4) = -sqrt2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
         +(cc(m1,1,2,k)+cc(m1,1,4,k)))
 1003          continue
  106 continue
  107 continue

  return
end subroutine mradb4 
subroutine mradb5 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! MRADB5 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,5,l1)
  real ( kind = 8 ) ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
!  arg= 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00 
   arg= 6.28318530717959d0 / dble(5.0)
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg)
  ti12=sin( 2.0E+00 *arg)

      do 101 k=1,l1
      m2 = m2s
      do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,k,1) = cc(m1,1,1,k)+ 2.0E+00 *cc(m1,ido,2,k)&
           + 2.0E+00 *cc(m1,ido,4,k)
         ch(m2,1,k,2) = (cc(m1,1,1,k)+tr11* 2.0E+00 *cc(m1,ido,2,k) &
         +tr12* 2.0E+00 *cc(m1,ido,4,k))-(ti11* 2.0E+00 *cc(m1,1,3,k) &
         +ti12* 2.0E+00 *cc(m1,1,5,k))
         ch(m2,1,k,3) = (cc(m1,1,1,k)+tr12* 2.0E+00 *cc(m1,ido,2,k) &
         +tr11* 2.0E+00 *cc(m1,ido,4,k))-(ti12* 2.0E+00 *cc(m1,1,3,k) &
         -ti11* 2.0E+00 *cc(m1,1,5,k))
         ch(m2,1,k,4) = (cc(m1,1,1,k)+tr12* 2.0E+00 *cc(m1,ido,2,k) &
         +tr11* 2.0E+00 *cc(m1,ido,4,k))+(ti12* 2.0E+00 *cc(m1,1,3,k) &
         -ti11* 2.0E+00 *cc(m1,1,5,k))
         ch(m2,1,k,5) = (cc(m1,1,1,k)+tr11* 2.0E+00 *cc(m1,ido,2,k) &
         +tr12* 2.0E+00 *cc(m1,ido,4,k))+(ti11* 2.0E+00 *cc(m1,1,3,k) &
         +ti12* 2.0E+00 *cc(m1,1,5,k))
 1001          continue
  101 continue

      if (ido == 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            m2 = m2s
      do 1002 m1=1,m1d,im1
        m2 = m2+im2
        ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))
        ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
        +(cc(m1,i,5,k)-cc(m1,ic,4,k))
        ch(m2,i-1,k,2) = wa1(i-2)*((cc(m1,i-1,1,k)+tr11* &
       (cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))+tr12 &
        *(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti11*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
        -wa1(i-1)*((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
        +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))+(ti11*(cc(m1,i-1,3,k) &
        -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) 

        ch(m2,i,k,2) = wa1(i-2)*((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k) &
        -cc(m1,ic,2,k))+tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
        +(ti11*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))+ti12 &
        *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))+wa1(i-1) &
        *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k) &
        +cc(m1,ic-1,2,k))+tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))) &
        -(ti11*(cc(m1,i,3,k)+cc(m1,ic,2,k))+ti12 &
        *(cc(m1,i,5,k)+cc(m1,ic,4,k))))
        ch(m2,i-1,k,3) = wa2(i-2) &
        *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti12*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
       -wa2(i-1) &
       *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
        cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
        +(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
        *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

        ch(m2,i,k,3) = wa2(i-2) &
       *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
        cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
        +(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
        *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
        +wa2(i-1) &
        *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(ti12*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

        ch(m2,i-1,k,4) = wa3(i-2) &
        *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti12*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
        -wa3(i-1) &
       *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
        cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
        -(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
        *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

        ch(m2,i,k,4) = wa3(i-2) &
       *((cc(m1,i,1,k)+tr12*(cc(m1,i,3,k)- &
        cc(m1,ic,2,k))+tr11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
        -(ti12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-ti11 &
        *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
        +wa3(i-1) &
        *((cc(m1,i-1,1,k)+tr12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti12*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))-ti11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

        ch(m2,i-1,k,5) = wa4(i-2) &
        *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti11*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
        -wa4(i-1) &
        *((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
        +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(ti11*(cc(m1,i-1,3,k) &
        -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

        ch(m2,i,k,5) = wa4(i-2) &
        *((cc(m1,i,1,k)+tr11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
        +tr12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(ti11*(cc(m1,i-1,3,k) &
        -cc(m1,ic-1,2,k))+ti12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
        +wa4(i-1) &
        *((cc(m1,i-1,1,k)+tr11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
        +tr12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(ti11*(cc(m1,i,3,k) &
        +cc(m1,ic,2,k))+ti12*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

 1002      continue
  102    continue
  103 continue

  return
end subroutine mradb5
subroutine mradbg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

!*****************************************************************************80
!
!! MRADBG is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  !tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
   tpi = 6.28318530717959d0
  arg = tpi / real ( ip, kind = 8 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  idp2 = ido + 2
  nbd = (ido-1)/2
  ipp2 = ip+2
  ipph = (ip+1)/2

      if (ido < l1) go to 103

      do k=1,l1
         do i=1,ido
            m2 = m2s
            do m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i,k,1) = cc(m1,i,1,k)
            end do
      end do
  end do

      go to 106

  103 do 105 i=1,ido
         do 104 k=1,l1
            m2 = m2s
            do 1004 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i,k,1) = cc(m1,i,1,k)
 1004       continue
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            m2 = m2s
            do 1007 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = cc(m1,ido,j2-2,k)+cc(m1,ido,j2-2,k)
            ch(m2,1,k,jc) = cc(m1,1,j2-1,k)+cc(m1,1,j2-1,k)
 1007       continue
  107    continue
  108 continue
      if (ido == 1) go to 116
      if (nbd < l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               m2 = m2s
               do 1009 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
               ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
               ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
               ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
 1009          continue
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               m2 = m2s
               do 1013 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
               ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
               ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
               ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
 1013          continue
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            m2 = m2s
            do 1017 m1=1,m1d,im1
            m2 = m2+im2
            c2(m1,ik,l) = ch2(m2,ik,1)+ar1*ch2(m2,ik,2)
            c2(m1,ik,lc) = ai1*ch2(m2,ik,ip)
 1017       continue
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               m2 = m2s
               do 1018 m1=1,m1d,im1
               m2 = m2+im2
               c2(m1,ik,l) = c2(m1,ik,l)+ar2*ch2(m2,ik,j)
               c2(m1,ik,lc) = c2(m1,ik,lc)+ai2*ch2(m2,ik,jc)
 1018          continue
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            m2 = m2s
            do 1021 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,1) = ch2(m2,ik,1)+ch2(m2,ik,j)
 1021       continue
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            m2 = m2s
            do 1023 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = c1(m1,1,k,j)-c1(m1,1,k,jc)
            ch(m2,1,k,jc) = c1(m1,1,k,j)+c1(m1,1,k,jc)
 1023       continue
  123    continue
  124 continue
      if (ido == 1) go to 132
      if (nbd < l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               m2 = m2s
               do 1025 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
               ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
               ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
               ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
 1025          continue
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               m2 = m2s
               do 1029 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
               ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
               ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
               ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
 1029          continue
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido == 1) return
      do 133 ik=1,idl1
         m2 = m2s
         do 1033 m1=1,m1d,im1
         m2 = m2+im2
         c2(m1,ik,1) = ch2(m2,ik,1)
 1033    continue
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            m2 = m2s
            do 1034 m1=1,m1d,im1
            m2 = m2+im2
            c1(m1,1,k,j) = ch(m2,1,k,j)
 1034       continue
  134    continue
  135 continue
      if (l1 < nbd ) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               m2 = m2s
               do 1036 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = wa(idij-1)*ch(m2,i-1,k,j)-wa(idij)* ch(m2,i,k,j)
               c1(m1,i,k,j) = wa(idij-1)*ch(m2,i,k,j)+wa(idij)* ch(m2,i-1,k,j)
 1036          continue
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               m2 = m2s
               do 1040 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = wa(idij-1)*ch(m2,i-1,k,j)-wa(idij)*ch(m2,i,k,j)
               c1(m1,i,k,j) = wa(idij-1)*ch(m2,i,k,j)+wa(idij)*ch(m2,i-1,k,j)
 1040          continue
  140       continue
  141    continue
  142 continue
  143 continue

  return
end subroutine mradbg
subroutine mradf2 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

!*****************************************************************************80
!
!! MRADF2 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,l1,2)
  real ( kind = 8 ) ch(in2,ido,2,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) wa1(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2

      do 101 k=1,l1
         m2 = m2s
         do m1=1,m1d,im1
           m2 = m2+im2
           ch(m2,1,1,k) = cc(m1,1,k,1)+cc(m1,1,k,2)
           ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,2)
         end do
  101 continue

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            m2 = m2s
            do 1003 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i,1,k) = cc(m1,i,k,1)+(wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))
            ch(m2,ic,2,k) = (wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
             cc(m1,i-1,k,2))-cc(m1,i,k,1)
            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+(wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))
            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)-(wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))
 1003       continue
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 do 106 k=1,l1
         m2 = m2s
         do 1006 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,2,k) = -cc(m1,ido,k,2)
         ch(m2,ido,1,k) = cc(m1,ido,k,1)
 1006    continue
  106 continue
  107 continue

  return
end subroutine mradf2
subroutine mradf3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

!*****************************************************************************80
!
!! MRADF3 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,3)
  real ( kind = 8 ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
 ! arg= 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  arg= 6.28318530717959d0 / dble(3.0)
  taur=cos(arg)
  taui=sin(arg)

      do 101 k=1,l1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,2)+cc(m1,1,k,3))
         ch(m2,1,3,k) = taui*(cc(m1,1,k,3)-cc(m1,1,k,2))
         ch(m2,ido,2,k) = cc(m1,1,k,1)+taur*(cc(m1,1,k,2)+cc(m1,1,k,3))
 1001    continue
  101 continue

  if (ido == 1) then
    return
  end if

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            m2 = m2s
            do 1002 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3)))

            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3)))

            ch(m2,i-1,3,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))+(taui*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))-(taui*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

            ch(m2,i,3,k) = (cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))))+(taui*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2))))

            ch(m2,ic,2,k) = (taui*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2))))-(cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))))
 1002       continue
  102    continue
  103 continue

  return
end subroutine mradf3
subroutine mradf4 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! MRADF4 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,l1,4)
  real ( kind = 8 ) ch(in2,ido,4,l1)
  real ( kind = 8 ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  hsqt2=sqrt( dble(2.0) ) / dble(2.0) 
  m1d = (m-1)*im1+1
  m2s = 1-im2

      do 101 k=1,l1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,1,k) = (cc(m1,1,k,2)+cc(m1,1,k,4)) &
            +(cc(m1,1,k,1)+cc(m1,1,k,3))
         ch(m2,ido,4,k) = (cc(m1,1,k,1)+cc(m1,1,k,3)) &
            -(cc(m1,1,k,2)+cc(m1,1,k,4))
         ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,3)
         ch(m2,1,3,k) = cc(m1,1,k,4)-cc(m1,1,k,2)
 1001    continue
  101 continue

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            m2 = m2s
            do 1003 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i-1,1,k) = ((wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4)))+(cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
             wa2(i-1)*cc(m1,i,k,3)))

            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
             wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
             wa3(i-1)*cc(m1,i,k,4)))

            ch(m2,i,1,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
             cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4)))+(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
             wa2(i-1)*cc(m1,i-1,k,3)))

            ch(m2,ic,4,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
             cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4)))-(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
             wa2(i-1)*cc(m1,i-1,k,3)))

            ch(m2,i-1,3,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
             cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4)))+(cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
             wa2(i-1)*cc(m1,i,k,3)))

            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
             wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
             cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4)))

            ch(m2,i,3,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))+(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
             wa2(i-1)*cc(m1,i-1,k,3)))

            ch(m2,ic,2,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))-(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
             wa2(i-1)*cc(m1,i-1,k,3)))

 1003       continue
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
         m2 = m2s
         do 1006 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,ido,1,k) = (hsqt2*(cc(m1,ido,k,2)-cc(m1,ido,k,4)))+ &
             cc(m1,ido,k,1)
            ch(m2,ido,3,k) = cc(m1,ido,k,1)-(hsqt2*(cc(m1,ido,k,2)- &
             cc(m1,ido,k,4)))
            ch(m2,1,2,k) = (-hsqt2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))- &
             cc(m1,ido,k,3)
            ch(m2,1,4,k) = (-hsqt2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))+ &
             cc(m1,ido,k,3)
 1006    continue
  106 continue
  107 continue

  return
end subroutine mradf4
subroutine mradf5 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! MRADF5 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,5)
  real ( kind = 8 ) ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
 ! arg= 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  arg= 6.28318530717959d0 / dble(5.0)
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg)
  ti12=sin( 2.0E+00 *arg)

      do 101 k=1,l1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
          (cc(m1,1,k,4)+cc(m1,1,k,3))
         ch(m2,ido,2,k) = cc(m1,1,k,1)+tr11*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
          tr12*(cc(m1,1,k,4)+cc(m1,1,k,3))
         ch(m2,1,3,k) = ti11*(cc(m1,1,k,5)-cc(m1,1,k,2))+ti12* &
          (cc(m1,1,k,4)-cc(m1,1,k,3))
         ch(m2,ido,4,k) = cc(m1,1,k,1)+tr12*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
          tr11*(cc(m1,1,k,4)+cc(m1,1,k,3))
         ch(m2,1,5,k) = ti12*(cc(m1,1,k,5)-cc(m1,1,k,2))-ti11* &
          (cc(m1,1,k,4)-cc(m1,1,k,3))
 1001    continue
  101 continue

      if (ido == 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            m2 = m2s
            do 1002 m1=1,m1d,im1
            m2 = m2+im2

            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
             cc(m1,i,k,5)))+((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
             wa3(i-1)*cc(m1,i,k,4)))

            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
             cc(m1,i-1,k,5)))+((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4)))

            ch(m2,i-1,3,k) = cc(m1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
             +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+tr12* &
            ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
             +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))+ti11* &
            ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
             -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
             -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4)))

            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
             +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+tr12* &
           ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
            +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))-(ti11* &
            ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
             -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
             -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

            ch(m2,i,3,k) = (cc(m1,i,k,1)+tr11*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
             cc(m1,i-1,k,5)))+tr12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(m1,i-1,k,5)+ &
             wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))+ti12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))))

            ch(m2,ic,2,k) = (ti11*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
             cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))+ti12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))))-(cc(m1,i,k,1)+tr11*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
             cc(m1,i-1,k,5)))+tr12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4))))

            ch(m2,i-1,5,k) = (cc(m1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
             cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
             cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))+(ti12*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
             cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
             cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
             cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
             cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))-(ti12*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
             cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
             cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

            ch(m2,i,5,k) = (cc(m1,i,k,1)+tr12*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
             cc(m1,i-1,k,5)))+tr11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(m1,i-1,k,5)+ &
             wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))-ti11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))))

            ch(m2,ic,4,k) = (ti12*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
             cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2)))-ti11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
             cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))))-(cc(m1,i,k,1)+tr12*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
             cc(m1,i-1,k,5)))+tr11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
             cc(m1,i-1,k,4))))
 1002       continue
  102    continue
  103 continue

  return
end subroutine mradf5
subroutine mradfg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

!*****************************************************************************80
!
!! MRADFG is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
  !tpi= 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  tpi= 6.28318530717959d0
  arg = tpi / real ( ip, kind = 8 )
  dcp = cos(arg)
  dsp = sin(arg)
  ipph = (ip+1)/2
  ipp2 = ip+2
  idp2 = ido+2
  nbd = (ido-1)/2

      if (ido == 1) go to 119

      do 101 ik=1,idl1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch2(m2,ik,1) = c2(m1,ik,1)
 1001    continue
  101 continue

      do 103 j=2,ip
         do 102 k=1,l1
            m2 = m2s
            do 1002 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = c1(m1,1,k,j)
 1002       continue
  102    continue
  103 continue
      if ( l1 < nbd ) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               m2 = m2s
               do 1004 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
               ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
 1004          continue
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               m2 = m2s
               do 1008 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
               ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
 1008          continue
  108       continue
  109    continue
  110 continue
  111 if (nbd < l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               m2 = m2s
               do 1012 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
               c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
 1012          continue
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               m2 = m2s
               do 1016 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
               c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
 1016          continue
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         m2 = m2s
         do 1020 m1=1,m1d,im1
         m2 = m2+im2
         c2(m1,ik,1) = ch2(m2,ik,1)
 1020    continue
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            m2 = m2s
            do 1022 m1=1,m1d,im1
            m2 = m2+im2
            c1(m1,1,k,j) = ch(m2,1,k,j)+ch(m2,1,k,jc)
            c1(m1,1,k,jc) = ch(m2,1,k,jc)-ch(m2,1,k,j)
 1022       continue
  122    continue
  123 continue

      ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            m2 = m2s
            do 1024 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,l) = c2(m1,ik,1)+ar1*c2(m1,ik,2)
            ch2(m2,ik,lc) = ai1*c2(m1,ik,ip)
 1024       continue
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               m2 = m2s
               do 1025 m1=1,m1d,im1
               m2 = m2+im2
               ch2(m2,ik,l) = ch2(m2,ik,l)+ar2*c2(m1,ik,j)
               ch2(m2,ik,lc) = ch2(m2,ik,lc)+ai2*c2(m1,ik,jc)
 1025          continue
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            m2 = m2s
            do 1028 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,1) = ch2(m2,ik,1)+c2(m1,ik,j)
 1028       continue
  128    continue
  129 continue

      if (ido < l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            m2 = m2s
            do 1030 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,i,1,k) = ch(m2,i,k,1)
 1030       continue
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            m2 = m2s
            do 1033 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,i,1,k) = ch(m2,i,k,1)
 1033       continue
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            m2 = m2s
            do 1036 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,ido,j2-2,k) = ch(m2,1,k,j)
            cc(m1,1,j2-1,k) = ch(m2,1,k,jc)
 1036       continue
  136    continue
  137 continue
      if (ido == 1) return
      if (nbd < l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               m2 = m2s
               do 1038 m1=1,m1d,im1
               m2 = m2+im2
               cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
               cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
 1038          continue
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               m2 = m2s
               do 1042 m1=1,m1d,im1
               m2 = m2+im2
               cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
               cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
 1042          continue
  142       continue
  143    continue
  144 continue

  return
end subroutine mradfg 
subroutine mrftb1 (m,im,n,in,c,ch,wa,fac)

!*****************************************************************************80
!
!! MRFTB1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(m,*)
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) half
  real ( kind = 8 ) halfm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) wa(n)

  nf = fac(2)
  na = 0

  do k1=1,nf
      ip = fac(k1+2)
      na = 1-na      
      if(ip <= 5) go to 10
      if(k1 == nf) go to 10
      na = 1-na
   10 continue 
  end do

      half = 0.5E+00
      halfm = -0.5E+00
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na == 0) go to 120
      m2 = 1-im
      do 117 i=1,m
      m2 = m2+im
      ch(i,1) = c(m2,1)
      ch(i,n) = c(m2,n)
  117 continue
      do 118 j=2,nl,2
      m2 = 1-im
      do 118 i=1,m
         m2 = m2+im
	 ch(i,j) = half*c(m2,j)
	 ch(i,j+1) = halfm*c(m2,j+1)
  118 continue
      go to 124
  120 continue
      do 122 j=2,nl,2
      m2 = 1-im
      do 122 i=1,m
         m2 = m2+im
	 c(m2,j) = half*c(m2,j)
	 c(m2,j+1) = halfm*c(m2,j+1)
  122 continue
  124 l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip /= 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
         call mradb4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call mradb4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip /= 2) go to 106
         if (na /= 0) go to 104
	 call mradb2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
         go to 105
  104    call mradb2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip /= 3) go to 109
         ix2 = iw+ido
         if (na /= 0) go to 107
	 call mradb3 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
         go to 108
  107    call mradb3 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip /= 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 110
         call mradb5 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call mradb5 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na /= 0) go to 113
	 call mradbg (m,ido,ip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
         go to 114
  113    call mradbg (m,ido,ip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
  114    if (ido == 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue

  return
end subroutine mrftb1
subroutine mrftf1 (m,im,n,in,c,ch,wa,fac)

!*****************************************************************************80
!
!! MRFTF1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(m,*)
  real ( kind = 8 ) fac(15)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) sn
  real ( kind = 8 ) tsn
  real ( kind = 8 ) tsnm
  real ( kind = 8 ) wa(n)

  nf = fac(2)
  na = 1
  l2 = n
  iw = n

      do 111 k1=1,nf
         kh = nf-k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip /= 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
	     call mradf4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call mradf4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip /= 2) go to 104
         if (na /= 0) go to 103
	     call mradf2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
         go to 110
  103    call mradf2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
         go to 110
  104    if (ip /= 3) go to 106
         ix2 = iw+ido
         if (na /= 0) go to 105
	     call mradf3 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
         go to 110
  105    call mradf3 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
         go to 110
  106    if (ip /= 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 107
         call mradf5(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call mradf5(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido == 1) na = 1-na
         if (na /= 0) go to 109
         call mradfg (m,ido,ip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
         na = 1
         go to 110
  109    call mradfg (m,ido,ip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
         na = 0
  110    l2 = l1
  111 continue

      sn = 1.0E+00 / real ( n, kind = 8 )
      tsn =  2.0E+00 / real ( n, kind = 8 )
      tsnm = -tsn
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na /= 0) go to 120
      m2 = 1-im

      do i=1,m
        m2 = m2+im
        c(m2,1) = sn*ch(i,1)
      end do

      do j=2,nl,2
        m2 = 1-im
        do i=1,m
          m2 = m2+im
	      c(m2,j) = tsn*ch(i,j)
	      c(m2,j+1) = tsnm*ch(i,j+1)
        end do
      end do

      if(modn /= 0) return
      m2 = 1-im
      do 119 i=1,m
         m2 = m2+im
         c(m2,n) = sn*ch(i,n)
  119 continue
      return
  120 m2 = 1-im
      do 121 i=1,m
         m2 = m2+im
         c(m2,1) = sn*c(m2,1)
  121 continue
      do 122 j=2,nl,2
      m2 = 1-im
      do 122 i=1,m
         m2 = m2+im
	 c(m2,j) = tsn*c(m2,j)
	 c(m2,j+1) = tsnm*c(m2,j+1)
  122 continue
      if(modn /= 0) return
      m2 = 1-im
      do i=1,m
         m2 = m2+im
         c(m2,n) = sn*c(m2,n)
      end do

  return
end subroutine mrftf1
subroutine mrfti1 (n,wa,fac)

!*****************************************************************************80
!
!! MRFTI1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization and 
!    other information is needed.
!
!    Output, real ( kind = 8 ) WA(N), trigonometric information.
!
!    Output, real ( kind = 8 ) FAC(15), factorization information.  FAC(1) is 
!    N, FAC(2) is NF, the number of factors, and FAC(3:NF+2) are the factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipm
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry
  integer ( kind = 4 ) ntryh(4)
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(n)

  save ntryh

  data ntryh / 4, 2, 3, 5 /

  nl = n
  nf = 0
  j = 0

  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      fac(nf+2) = ntry
      nl = nq
      if (ntry /= 2) go to 107
      do i=2,nf
         ib = nf-i+2
         fac(ib+2) = fac(ib+1)
      end do
      fac(3) = 2
  107 if (nl /= 1) go to 104
      fac(1) = n
      fac(2) = nf
     ! tpi = 8.0D+00 * atan ( 1.0D+00 )
       tpi = 6.28318530717959d0
      argh = tpi / real ( n, kind = 8 )
      is = 0
      nfm1 = nf-1
      l1 = 1

      do k1=1,nfm1
         ip = fac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do j=1,ipm
            ld = ld+l1
            i = is
            argld = real ( ld, kind = 8 ) * argh
            fi = 0.0E+00
            do ii=3,ido,2
              i = i+2
              fi = fi + 1.0E+00
              arg = fi*argld
	          wa(i-1) = cos ( arg )
	          wa(i) = sin ( arg )
            end do
            is = is+ido
         end do
         l1 = l2
      end do

  return
end subroutine mrfti1
subroutine msntb1(lot,jump,n,inc,x,wsave,dsum,xh,work,ier)

!*****************************************************************************80
!
!! MSNTB1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lot

  real ( kind = 8 ) dsum(*)
  real ( kind = 8 ) fnp1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lnxh
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) srt3s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xh(lot,*)
  real ( kind = 8 ) xhold

  ier = 0
  lj = (lot-1)*jump+1

      if (n-2) 200,102,103
 102  srt3s2 = sqrt( dble(3.0) )/ dble(2.0) 

      do m=1,lj,jump
         xhold = srt3s2*(x(m,1)+x(m,2))
         x(m,2) = srt3s2*(x(m,1)-x(m,2))
         x(m,1) = xhold
      end do

      go to 200
  103 np1 = n+1
      ns2 = n/2
      do 104 k=1,ns2
         kc = np1-k
         m1 = 0
         do 114 m=1,lj,jump
         m1 = m1+1
         t1 = x(m,k)-x(m,kc)
         t2 = wsave(k)*(x(m,k)+x(m,kc))
         xh(m1,k+1) = t1+t2
         xh(m1,kc+1) = t2-t1
  114    continue
  104 continue
      modn = mod(n,2)
      if (modn == 0) go to 124
      m1 = 0
      do 123 m=1,lj,jump
         m1 = m1+1
         xh(m1,ns2+2) =  4.0E+00 * x(m,ns2+1)
  123 continue
  124 do m=1,lot
         xh(m,1) = 0.0E+00
      end do

      lnxh = lot-1 + lot*(np1-1) + 1
      lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = lot*np1

      call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,ier1)     

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('msntb1',-5)
        go to 200
      end if

      if(mod(np1,2) /= 0) go to 30
      do m=1,lot
        xh(m,np1) = xh(m,np1)+xh(m,np1)
      end do
 30   fnp1s4 = real ( np1 ) / 4.0E+00
      m1 = 0
      do 125 m=1,lj,jump
         m1 = m1+1
         x(m,1) = fnp1s4*xh(m1,1)
         dsum(m1) = x(m,1)
  125 continue
      do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            x(m,i-1) = fnp1s4*xh(m1,i)
            dsum(m1) = dsum(m1)+fnp1s4*xh(m1,i-1)
            x(m,i) = dsum(m1)
  115    continue
  105 continue
      if (modn /= 0) go to 200
      m1 = 0
      do 116 m=1,lj,jump
         m1 = m1+1
         x(m,n) = fnp1s4*xh(m1,n+1)
  116 continue

  200 continue

  return
end subroutine msntb1
subroutine msntf1(lot,jump,n,inc,x,wsave,dsum,xh,work,ier)

!*****************************************************************************80
!
!! MSNTF1 is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lot

  real ( kind = 8 ) dsum(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lnxh
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) sfnp1
  real ( kind = 8 ) ssqrt3
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xh(lot,*)
  real ( kind = 8 ) xhold

  ier = 0
  lj = (lot-1)*jump+1

      if (n-2) 101,102,103
 102  ssqrt3 = dble(1.0) / sqrt (dble(3.0) )

      do m=1,lj,jump
         xhold = ssqrt3*(x(m,1)+x(m,2))
         x(m,2) = ssqrt3*(x(m,1)-x(m,2))
         x(m,1) = xhold
      end do

  101  go to 200
  103 np1 = n+1
      ns2 = n/2
      do 104 k=1,ns2
         kc = np1-k
         m1 = 0
         do 114 m=1,lj,jump
         m1 = m1 + 1
         t1 = x(m,k)-x(m,kc)
         t2 = wsave(k)*(x(m,k)+x(m,kc))
         xh(m1,k+1) = t1+t2
         xh(m1,kc+1) = t2-t1
  114    continue
  104 continue
      modn = mod(n,2)
      if (modn == 0) go to 124
      m1 = 0
      do 123 m=1,lj,jump
         m1 = m1 + 1
         xh(m1,ns2+2) =  4.0E+00  * x(m,ns2+1)
  123 continue
  124 do 127 m=1,lot
         xh(m,1) = 0.0E+00
  127 continue 
      lnxh = lot-1 + lot*(np1-1) + 1
      lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = lot*np1

      call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,ier1)     
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('msntf1',-5)
        go to 200
      end if

      if(mod(np1,2) /= 0) go to 30
      do 20 m=1,lot
      xh(m,np1) = xh(m,np1)+xh(m,np1)
   20 continue
   30 sfnp1 = 1.0E+00 / real ( np1, kind = 8 )
      m1 = 0
      do 125 m=1,lj,jump
         m1 = m1+1
         x(m,1) = 0.5E+00 * xh(m1,1)
         dsum(m1) = x(m,1)
  125 continue
      do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            x(m,i-1) = 0.5E+00 * xh(m1,i)
            dsum(m1) = dsum(m1)+ 0.5E+00 * xh(m1,i-1)
            x(m,i) = dsum(m1)
  115    continue
  105 continue
      if (modn /= 0) go to 200

      m1 = 0
      do m=1,lj,jump
         m1 = m1+1
         x(m,n) = 0.5E+00 * xh(m1,n+1)
      end do

  200 continue

  return
end subroutine msntf1
subroutine r1f2kb (ido,l1,cc,in1,ch,in2,wa1)

!*****************************************************************************80
!
!! R1F2KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,2,l1)
  real ( kind = 8 ) ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)

  do k=1,l1
    ch(1,1,k,1) = cc(1,1,1,k)+cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k)-cc(1,ido,2,k)
  end do

      if (ido-2) 107,105,102
 102  idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            
            ch(1,i-1,k,1) = cc(1,i-1,1,k)+cc(1,ic-1,2,k)
            ch(1,i,k,1) = cc(1,i,1,k)-cc(1,ic,2,k)
            
            ch(1,i-1,k,2) = wa1(i-2)*(cc(1,i-1,1,k)-cc(1,ic-1,2,k)) &
                 -wa1(i-1)*(cc(1,i,1,k)+cc(1,ic,2,k))
            ch(1,i,k,2) = wa1(i-2)*(cc(1,i,1,k)+cc(1,ic,2,k))+wa1(i-1) &
                 *(cc(1,i-1,1,k)-cc(1,ic-1,2,k))

 103     continue
 104  continue
      if (mod(ido,2) == 1) return
 105  do 106 k=1,l1
         ch(1,ido,k,1) = cc(1,ido,1,k)+cc(1,ido,1,k)
         ch(1,ido,k,2) = -(cc(1,1,2,k)+cc(1,1,2,k))
 106  continue
 107  continue

  return
end subroutine r1f2kb
subroutine r1f2kf (ido,l1,cc,in1,ch,in2,wa1)

!*****************************************************************************80
!
!! R1F1KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ch(in2,ido,2,l1)
  real ( kind = 8 ) cc(in1,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+cc(1,1,k,2)
    ch(1,ido,2,k) = cc(1,1,k,1)-cc(1,1,k,2)
  end do

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ch(1,i,1,k) = cc(1,i,k,1)+(wa1(i-2)*cc(1,i,k,2) &
              -wa1(i-1)*cc(1,i-1,k,2))
            ch(1,ic,2,k) = (wa1(i-2)*cc(1,i,k,2) &
              -wa1(i-1)*cc(1,i-1,k,2))-cc(1,i,k,1)
            ch(1,i-1,1,k) = cc(1,i-1,k,1)+(wa1(i-2)*cc(1,i-1,k,2) &
              +wa1(i-1)*cc(1,i,k,2))
            ch(1,ic-1,2,k) = cc(1,i-1,k,1)-(wa1(i-2)*cc(1,i-1,k,2) &
              +wa1(i-1)*cc(1,i,k,2))
          end do
  104 continue
      if (mod(ido,2) == 1) return
  105 do 106 k=1,l1
         ch(1,1,2,k) = -cc(1,ido,k,2)
         ch(1,ido,1,k) = cc(1,ido,k,1)
  106 continue
  107 continue

  return
end subroutine r1f2kf
subroutine r1f3kb (ido,l1,cc,in1,ch,in2,wa1,wa2)

!*****************************************************************************80
!
!! R1F3KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,3,l1)
  real ( kind = 8 ) ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

 ! arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00 
   arg = 6.28318530717959d0 / dble(3)
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) + ( 2.0E+00 * taur ) * cc(1,ido,2,k) &
      - ( 2.0E+00 *taui)*cc(1,1,3,k)
    ch(1,1,k,3) = cc(1,1,1,k) + ( 2.0E+00 *taur)*cc(1,ido,2,k) &
      + 2.0E+00 *taui*cc(1,1,3,k)
  end do

  if (ido == 1) then
    return
  end if

  idp2 = ido+2

      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k))

        ch(1,i-1,k,2) = wa1(i-2)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) &
                         -wa1(i-1)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) 

            ch(1,i,k,2) = wa1(i-2)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) &
                        +wa1(i-1)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))

              ch(1,i-1,k,3) = wa2(i-2)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) &
         -wa2(i-1)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))

            ch(1,i,k,3) = wa2(i-2)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) &
                       +wa2(i-1)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))

  102    continue
  103 continue

  return
end subroutine r1f3kb
subroutine r1f3kf (ido,l1,cc,in1,ch,in2,wa1,wa2)

!*****************************************************************************80
!
!! R1F3KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,3)
  real ( kind = 8 ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

!  arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )/ 3.0E+00 
   arg= 6.28318530717959d0 / dble(3.0)
  taur=cos(arg)
  taui=sin(arg)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+(cc(1,1,k,2)+cc(1,1,k,3))
    ch(1,1,3,k) = taui*(cc(1,1,k,3)-cc(1,1,k,2))
    ch(1,ido,2,k) = cc(1,1,k,1)+taur*(cc(1,1,k,2)+cc(1,1,k,3))
  end do

  if (ido == 1) then
    return
  end if

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i

            ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
             wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3)))

            ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3)))

            ch(1,i-1,3,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))+(taui*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))

            ch(1,ic-1,2,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))-(taui*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))

            ch(1,i,3,k) = (cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))))+(taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))))

            ch(1,ic,2,k) = (taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))))-(cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))))
  102    continue
  103 continue

  return
end subroutine r1f3kf
subroutine r1f4kb (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! R1F4KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,4,l1)
  real ( kind = 8 ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) sqrt2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  sqrt2=sqrt( dble(2.0E+00) )

  do k=1,l1
    ch(1,1,k,3) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      -(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,1) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      +(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,4) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      +(cc(1,1,3,k)+cc(1,1,3,k))
    ch(1,1,k,2) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      -(cc(1,1,3,k)+cc(1,1,3,k))
  end do

      if (ido-2) 107,105,102

  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = (cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        +(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = (cc(1,i,1,k)-cc(1,ic,4,k)) &
        +(cc(1,i,3,k)-cc(1,ic,2,k))
        ch(1,i-1,k,2)=wa1(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        -(cc(1,i,3,k)+cc(1,ic,2,k)))-wa1(i-1) &
        *((cc(1,i,1,k)+cc(1,ic,4,k))+(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,2)=wa1(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        +(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa1(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))-(cc(1,i,3,k)+cc(1,ic,2,k)))
        ch(1,i-1,k,3)=wa2(i-2)*((cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        -(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))-wa2(i-1) &
        *((cc(1,i,1,k)-cc(1,ic,4,k))-(cc(1,i,3,k)-cc(1,ic,2,k)))
        ch(1,i,k,3)=wa2(i-2)*((cc(1,i,1,k)-cc(1,ic,4,k)) &
        -(cc(1,i,3,k)-cc(1,ic,2,k)))+wa2(i-1) &
        *((cc(1,i-1,1,k)+cc(1,ic-1,4,k))-(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k)))
        ch(1,i-1,k,4)=wa3(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        +(cc(1,i,3,k)+cc(1,ic,2,k)))-wa3(i-1) &
       *((cc(1,i,1,k)+cc(1,ic,4,k))-(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,4)=wa3(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        -(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa3(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))+(cc(1,i,3,k)+cc(1,ic,2,k)))
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
         ch(1,ido,k,1) = (cc(1,ido,1,k)+cc(1,ido,3,k)) &
         +(cc(1,ido,1,k)+cc(1,ido,3,k))
         ch(1,ido,k,2) = sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         -(cc(1,1,2,k)+cc(1,1,4,k)))
         ch(1,ido,k,3) = (cc(1,1,4,k)-cc(1,1,2,k)) &
         +(cc(1,1,4,k)-cc(1,1,2,k))
         ch(1,ido,k,4) = -sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         +(cc(1,1,2,k)+cc(1,1,4,k)))
  106 continue
  107 continue

  return
end subroutine r1f4kb
subroutine r1f4kf (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! R1F4KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,l1,4)
  real ( kind = 8 ) ch(in2,ido,4,l1)
  real ( kind = 8 ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  hsqt2=sqrt( dble(2.0) )/ dble(2.0) 

  do k=1,l1
    ch(1,1,1,k) = (cc(1,1,k,2)+cc(1,1,k,4))+(cc(1,1,k,1)+cc(1,1,k,3))
    ch(1,ido,4,k) = (cc(1,1,k,1)+cc(1,1,k,3))-(cc(1,1,k,2)+cc(1,1,k,4))
    ch(1,ido,2,k) = cc(1,1,k,1)-cc(1,1,k,3)
    ch(1,1,3,k) = cc(1,1,k,4)-cc(1,1,k,2)
  end do

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(1,i-1,1,k) = ((wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4)))+(cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))
            ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i-1,k,2)+ &
             wa1(i-1)*cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+ &
             wa3(i-1)*cc(1,i,k,4)))
            ch(1,i,1,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))+(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,ic,4,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))-(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,i-1,3,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))+(cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))
            ch(1,ic-1,2,k) = (cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))
            ch(1,i,3,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,ic,2,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
            ch(1,ido,1,k) = (hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))+cc(1,ido,k,1)
            ch(1,ido,3,k) = cc(1,ido,k,1)-(hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))
            ch(1,1,2,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))-cc(1,ido,k,3)
            ch(1,1,4,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))+cc(1,ido,k,3)
  106 continue
  107 continue

  return
end subroutine r1f4kf
subroutine r1f5kb (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! R1F5KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,5,l1)
  real ( kind = 8 ) ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

 ! arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 ) / 5.0E+00
   arg= 6.28318530717959d0 / dble(5.0)
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg )
  ti12=sin( 2.0E+00 *arg )

  do k=1,l1
    ch(1,1,k,1) = cc(1,1,1,k)+ 2.0E+00 *cc(1,ido,2,k)+ 2.0E+00 *cc(1,ido,4,k)
    ch(1,1,k,2) = (cc(1,1,1,k)+tr11* 2.0E+00 *cc(1,ido,2,k) &
      +tr12* 2.0E+00 *cc(1,ido,4,k))-(ti11* 2.0E+00 *cc(1,1,3,k) &
      +ti12* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,3) = (cc(1,1,1,k)+tr12* 2.0E+00 *cc(1,ido,2,k) &
      +tr11* 2.0E+00 *cc(1,ido,4,k))-(ti12* 2.0E+00 *cc(1,1,3,k) &
      -ti11* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,4) = (cc(1,1,1,k)+tr12* 2.0E+00 *cc(1,ido,2,k) &
      +tr11* 2.0E+00 *cc(1,ido,4,k))+(ti12* 2.0E+00 *cc(1,1,3,k) &
      -ti11* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,5) = (cc(1,1,1,k)+tr11* 2.0E+00 *cc(1,ido,2,k) &
      +tr12* 2.0E+00 *cc(1,ido,4,k))+(ti11* 2.0E+00 *cc(1,1,3,k) &
      +ti12* 2.0E+00 *cc(1,1,5,k))
  end do

  if (ido == 1) return

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +(cc(1,i-1,5,k)+cc(1,ic-1,4,k))
        ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +(cc(1,i,5,k)-cc(1,ic,4,k))
        ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)+tr11* &
        (cc(1,i-1,3,k)+cc(1,ic-1,2,k))+tr12 &
        *(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa1(i-1)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))+(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k) &
        -cc(1,ic,2,k))+tr12*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti11*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))+ti12 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))+wa1(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k))+tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k))) &
        -(ti11*(cc(1,i,3,k)+cc(1,ic,2,k))+ti12 &
        *(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,3) = wa2(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
       -wa2(i-1) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,3) = wa2(i-2) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa2(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,4) = wa3(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa3(i-1) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,4) = wa3(i-2) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa3(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,5) = wa4(i-2) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa4(i-1) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,5) = wa4(i-2) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa4(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k))))

  102    continue
  103 continue

  return
end subroutine r1f5kb
subroutine r1f5kf (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! R1F5KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,5)
  real ( kind = 8 ) ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

 ! arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 ) / 5.0E+00
   arg= 6.28318530717959d0 / dble(5.0)
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg)
  ti12=sin( 2.0E+00 *arg)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+(cc(1,1,k,5)+cc(1,1,k,2))+ &
      (cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,ido,2,k) = cc(1,1,k,1)+tr11*(cc(1,1,k,5)+cc(1,1,k,2))+ &
      tr12*(cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,1,3,k) = ti11*(cc(1,1,k,5)-cc(1,1,k,2))+ti12* &
      (cc(1,1,k,4)-cc(1,1,k,3))
    ch(1,ido,4,k) = cc(1,1,k,1)+tr12*(cc(1,1,k,5)+cc(1,1,k,2))+ &
      tr11*(cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,1,5,k) = ti12*(cc(1,1,k,5)-cc(1,1,k,2))-ti11* &
      (cc(1,1,k,4)-cc(1,1,k,3))
  end do

  if (ido == 1) return

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i

            ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
            wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
            cc(1,i,k,5)))+((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
            cc(1,i,k,3))+(wa3(i-2)*cc(1,i-1,k,4)+ &
            wa3(i-1)*cc(1,i,k,4))) 

            ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))

            ch(1,i-1,3,k) = cc(1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
             +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
            ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
             +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))+ti11* &
            ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
             -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
             -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4)))

            ch(1,ic-1,2,k) = cc(1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
             +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
           ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
            +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))-(ti11* &
            ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
             -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
             -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,i,3,k) = (cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(1,i-1,k,5)+ &
             wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))

            ch(1,ic,2,k) = (ti11*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
             cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))-(cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))

            ch(1,i-1,5,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
             cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
             cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))+(ti12*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
             cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
             cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
             cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
             cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))-(ti12*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
             cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
             cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,i,5,k) = (cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(1,i-1,k,5)+ &
             wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))

            ch(1,ic,4,k) = (ti12*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
             cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))-(cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))

  102    continue
  103 continue

  return
end subroutine r1f5kf
subroutine r1fgkb (ido,ip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)

!*****************************************************************************80
!
!! R1FGKB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

 ! tpi= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )
   tpi= 6.28318530717959d0
  arg = tpi / real ( ip, kind = 8 )
  dcp = cos(arg)
  dsp = sin(arg)
  idp2 = ido+2
  nbd = (ido-1)/2
  ipp2 = ip+2
  ipph = (ip+1)/2

  if (ido < l1) go to 103

      do k=1,l1
         do i=1,ido
            ch(1,i,k,1) = cc(1,i,1,k)
         end do
      end do

      go to 106

  103 continue

  do i=1,ido
    do k=1,l1
      ch(1,i,k,1) = cc(1,i,1,k)
    end do
  end do

  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,1,k,j) = cc(1,ido,j2-2,k)+cc(1,ido,j2-2,k)
            ch(1,1,k,jc) = cc(1,1,j2-1,k)+cc(1,1,j2-1,k)
 1007       continue
  107    continue
  108 continue
      if (ido == 1) go to 116
      if (nbd < l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
               ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
               ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
               ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
               ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
               ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
               ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(1,ik,l) = ch2(1,ik,1)+ar1*ch2(1,ik,2)
            c2(1,ik,lc) = ai1*ch2(1,ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(1,ik,l) = c2(1,ik,l)+ar2*ch2(1,ik,j)
               c2(1,ik,lc) = c2(1,ik,lc)+ai2*ch2(1,ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(1,ik,1) = ch2(1,ik,1)+ch2(1,ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,1,k,j) = c1(1,1,k,j)-c1(1,1,k,jc)
            ch(1,1,k,jc) = c1(1,1,k,j)+c1(1,1,k,jc)
  123    continue
  124 continue
      if (ido == 1) go to 132
      if (nbd < l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
               ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
               ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
               ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
               ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
               ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
               ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido == 1) return
      do 133 ik=1,idl1
         c2(1,ik,1) = ch2(1,ik,1)
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,1,k,j) = ch(1,1,k,j)
  134    continue
  135 continue
      if ( l1 < nbd ) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)*ch(1,i,k,j)
               c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)*ch(1,i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)*ch(1,i,k,j)
               c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)*ch(1,i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 continue

  return
end subroutine r1fgkb
subroutine r1fgkf (ido,ip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)

!*****************************************************************************80
!
!! R1FGKF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

 ! tpi= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )
   tpi= 6.28318530717959d0
  arg = tpi/real ( ip, kind = 8 )
  dcp = cos(arg)
  dsp = sin(arg)
  ipph = (ip+1)/2
  ipp2 = ip+2
  idp2 = ido+2
  nbd = (ido-1)/2

      if (ido == 1) go to 119

      do ik=1,idl1
         ch2(1,ik,1) = c2(1,ik,1)
      end do

      do j=2,ip
         do k=1,l1
            ch(1,1,k,j) = c1(1,1,k,j)
         end do
      end do

      if ( l1 < nbd ) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij)*c1(1,i,k,j)
               ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij)*c1(1,i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij)*c1(1,i,k,j)
               ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij)*c1(1,i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd < l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
               c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
               c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
               c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
               c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(1,ik,1) = ch2(1,ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,1,k,j) = ch(1,1,k,j)+ch(1,1,k,jc)
            c1(1,1,k,jc) = ch(1,1,k,jc)-ch(1,1,k,j)
  122    continue
  123 continue

      ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(1,ik,l) = c2(1,ik,1)+ar1*c2(1,ik,2)
            ch2(1,ik,lc) = ai1*c2(1,ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(1,ik,l) = ch2(1,ik,l)+ar2*c2(1,ik,j)
               ch2(1,ik,lc) = ch2(1,ik,lc)+ai2*c2(1,ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(1,ik,1) = ch2(1,ik,1)+c2(1,ik,j)
  128    continue
  129 continue

      if (ido < l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(1,i,1,k) = ch(1,i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(1,i,1,k) = ch(1,i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(1,ido,j2-2,k) = ch(1,1,k,j)
            cc(1,1,j2-1,k) = ch(1,1,k,jc)
  136    continue
  137 continue
      if (ido == 1) return
      if (nbd < l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
               cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
               cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
               cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
               cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
  142       continue
  143    continue
  144 continue

  return
end subroutine r1fgkf
subroutine r2w ( ldr, ldw, l, m, r, w )

!*****************************************************************************80
!
!! R2W copies a 2D array, allowing for different leading dimensions.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldw
  integer ( kind = 4 ) m

  integer ( kind = 4 ) l
  real ( kind = 8 ) r(ldr,m)
  real ( kind = 8 ) w(ldw,m)

  w(1:l,1:m) = r(1:l,1:m)

  return
end subroutine r2w
subroutine r4_factor ( n, nf, fac )

!*****************************************************************************80
!
!! R4_FACTOR factors of an integer for real single precision computations.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization and 
!    other information is needed.
!
!    Output, integer ( kind = 4 ) NF, the number of factors.
!
!    Output, real ( kind = 8 ) FAC(*), a list of factors of N.
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = real ( ntry, kind = 8 )
      nl = nq

    end do

  end do

  return
end subroutine r4_factor
subroutine r4_mcfti1 ( n, wa, fnf, fac )

!*****************************************************************************80
!
!! R4_MCFTI1 sets up factors and tables, real single precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)
!
!  Get the factorization of N.
!
  call r4_factor ( n, nf, fac )
  fnf = real ( nf, kind = 8 )
  iw = 1
  l1 = 1
!
!  Set up the trigonometric tables.
!
  do k1 = 1, nf
    ip = int ( fac(k1) )
    l2 = l1 * ip
    ido = n / l2
    call r4_tables ( ido, ip, wa(iw) )
    iw = iw + ( ip - 1 ) * ( ido + ido )
    l1 = l2
  end do

  return
end subroutine r4_mcfti1
subroutine r4_tables ( ido, ip, wa )

!*****************************************************************************80
!
!! R4_TABLES computes trigonometric tables, real single precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) argz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido,ip-1,2)

 ! tpi = 8.0E+00 * atan ( 1.0E+00 )
   tpi= 6.28318530717959d0 
  argz = tpi / real ( ip, kind = 8 )
  arg1 = tpi / real ( ido * ip, kind = 8 )

  do j = 2, ip

    arg2 = real ( j - 1, kind = 8 ) * arg1

    do i = 1, ido
      arg3 = real ( i - 1, kind = 8 ) * arg2
      wa(i,j-1,1) = cos ( arg3 )
      wa(i,j-1,2) = sin ( arg3 )
    end do

    if ( 5 < ip ) then
      arg4 = real ( j - 1, kind = 8 ) * argz
      wa(1,j-1,1) = cos ( arg4 )
      wa(1,j-1,2) = sin ( arg4 )
    end if

  end do

  return
end subroutine r4_tables
subroutine r8_factor ( n, nf, fac )

!*****************************************************************************80
!
!! R8_FACTOR factors of an integer for real double precision computations.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization and 
!    other information is needed.
!
!    Output, integer ( kind = 4 ) NF, the number of factors.
!
!    Output, real ( kind = 8 ) FAC(*), a list of factors of N.
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = real ( ntry, kind = 8 )
      nl = nq

    end do

  end do

  return
end subroutine r8_factor
subroutine r8_mcfti1 ( n, wa, fnf, fac )

!*****************************************************************************80
!
!! R8_MCFTI1 sets up factors and tables, real double precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)
!
!  Get the factorization of N.
!
  call r8_factor ( n, nf, fac )
  fnf = real ( nf, kind = 8 )
  iw = 1
  l1 = 1
!
!  Set up the trigonometric tables.
!
  do k1 = 1, nf
    ip = int ( fac(k1) )
    l2 = l1 * ip
    ido = n / l2
    call r8_tables ( ido, ip, wa(iw) )
    iw = iw + ( ip - 1 ) * ( ido + ido )
    l1 = l2
  end do

  return
end subroutine r8_mcfti1
subroutine r8_tables ( ido, ip, wa )

!*****************************************************************************80
!
!! R8_TABLES computes trigonometric tables, real double precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) argz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido,ip-1,2)

 ! tpi = 8.0D+00 * atan ( 1.0D+00 )
   tpi= 6.28318530717959d0
  argz = tpi / real ( ip, kind = 8 )
  arg1 = tpi / real ( ido * ip, kind = 8 )

  do j = 2, ip

    arg2 = real ( j - 1, kind = 8 ) * arg1

    do i = 1, ido
      arg3 = real ( i - 1, kind = 8 ) * arg2
      wa(i,j-1,1) = cos ( arg3 )
      wa(i,j-1,2) = sin ( arg3 )
    end do

    if ( 5 < ip ) then
      arg4 = real ( j - 1, kind = 8 ) * argz
      wa(1,j-1,1) = cos ( arg4 )
      wa(1,j-1,2) = sin ( arg4 )
    end if

  end do

  return
end subroutine r8_tables
subroutine rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1B: real single precision backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1B computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the backward
!    transform or Fourier synthesis, transforming the sequence from 
!    spectral to physical space.  This transform is normalized since a 
!    call to RFFT1B followed by a call to RFFT1F (or vice-versa) reproduces
!    the original array within roundoff error. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, the data to be 
!    transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1. 
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine
!    RFFT1F or RFFT1B for a given transform length N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4. 
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(lenr)
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenr < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfft1b ', 6)
  else if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1b ', 8)
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('rfft1b ', 10)
  end if

  if (n == 1) then
    return
  end if

  call rfftb1 (n,inc,r,work,wsave,wsave(n+1))

  return
end subroutine rfft1b
subroutine rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1F: real single precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1F computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the forward 
!    transform or Fourier analysis, transforming the sequence from physical
!    to spectral space.  This transform is normalized since a call to 
!    RFFT1F followed by a call to RFFT1B (or vice-versa) reproduces the 
!    original array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, contains the sequence 
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine RFFT1F
!    or RFFT1B for a given transform length N. 
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK). 
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array. 
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough:
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) r(lenr)

  ier = 0

  if (lenr < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfft1f ', 6)
  else if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1f ', 8)
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('rfft1f ', 10)
  end if

  if (n == 1) then
    return
  end if

  call rfftf1 (n,inc,r,work,wsave,wsave(n+1))

  return
end subroutine rfft1f
subroutine rfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! RFFT1I: initialization for RFFT1B and RFFT1F.
!
!  Discussion:
!
!    RFFT1I initializes array WSAVE for use in its companion routines 
!    RFFT1B and RFFT1F.  The prime factorization of N together with a
!    tabulation of the trigonometric functions are computed and stored
!    in array WSAVE.  Separate WSAVE arrays are required for different
!    values of N. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of
!    N and also containing certain trigonometric values which will be used in
!    routines RFFT1B or RFFT1F.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1i ', 3)
  end if

  if (n == 1) then
    return
  end if

  call rffti1 (n,wsave(1),wsave(n+1))

  return
end subroutine rfft1i

subroutine rfftb1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTB1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) half
  real ( kind = 8 ) halfm
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) wa(n)

      nf = fac(2)
      na = 0
      do 10 k1=1,nf
      ip = fac(k1+2)
      na = 1-na      
      if(ip <= 5) go to 10
      if(k1 == nf) go to 10
      na = 1-na
   10 continue 
      half = dble(0.5)
      halfm = dble(-0.5)
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na == 0) go to 120

      ch(1) = c(1,1)
      ch(n) = c(1,n)
      do j=2,nl,2
	    ch(j) = half*c(1,j)
	    ch(j+1) = halfm*c(1,j+1)
      end do

      go to 124
  120 do 122 j=2,nl,2
	 c(1,j) = half*c(1,j)
	 c(1,j+1) = halfm*c(1,j+1)
  122 continue
  124 l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip /= 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
         call r1f4kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call r1f4kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip /= 2) go to 106
         if (na /= 0) go to 104
	 call r1f2kb (ido,l1,c,in,ch,1,wa(iw))
         go to 105
  104    call r1f2kb (ido,l1,ch,1,c,in,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip /= 3) go to 109
         ix2 = iw+ido
         if (na /= 0) go to 107
	 call r1f3kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2))
         go to 108
  107    call r1f3kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip /= 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 110
         call r1f5kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call r1f5kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na /= 0) go to 113
	 call r1fgkb (ido,ip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw))
         go to 114
  113    call r1fgkb (ido,ip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw))
  114    if (ido == 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue

  return
end subroutine rfftb1
subroutine rfftf1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTF1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(15)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) sn
  real ( kind = 8 ) tsn
  real ( kind = 8 ) tsnm
  real ( kind = 8 ) wa(n)

      nf = int ( fac(2) )
      na = 1
      l2 = n
      iw = n

      do 111 k1=1,nf
         kh = nf-k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip /= 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
	     call r1f4kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call r1f4kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip /= 2) go to 104
         if (na /= 0) go to 103
	 call r1f2kf (ido,l1,c,in,ch,1,wa(iw))
         go to 110
  103    call r1f2kf (ido,l1,ch,1,c,in,wa(iw))
         go to 110
  104    if (ip /= 3) go to 106
         ix2 = iw+ido
         if (na /= 0) go to 105
	 call r1f3kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2))
         go to 110
  105    call r1f3kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2))
         go to 110
  106    if (ip /= 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 107
         call r1f5kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call r1f5kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido == 1) na = 1-na
         if (na /= 0) go to 109
	 call r1fgkf (ido,ip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw))
         na = 1
         go to 110
  109    call r1fgkf (ido,ip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw))
         na = 0
  110    l2 = l1
  111 continue

      sn = 1.0E+00 / real ( n, kind = 8 )
      tsn =  2.0E+00  / real ( n, kind = 8 )
      tsnm = -tsn
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na /= 0) go to 120
      c(1,1) = sn*ch(1)
      do 118 j=2,nl,2
	 c(1,j) = tsn*ch(j)
	 c(1,j+1) = tsnm*ch(j+1)
  118 continue
      if(modn /= 0) return
      c(1,n) = sn*ch(n)
      return
  120 c(1,1) = sn*c(1,1)
      do 122 j=2,nl,2
	 c(1,j) = tsn*c(1,j)
	 c(1,j+1) = tsnm*c(1,j+1)
  122 continue
      if(modn /= 0) return
      c(1,n) = sn*c(1,n)

  return
end subroutine rfftf1
subroutine rffti1 ( n, wa, fac )

!*****************************************************************************80
!
!! RFFTI1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization 
!    and other information is needed.
!
!    Output, real ( kind = 8 ) WA(N), trigonometric information.
!
!    Output, real ( kind = 8 ) FAC(15), factorization information.  
!    FAC(1) is N, FAC(2) is NF, the number of factors, and FAC(3:NF+2) are the 
!    factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipm
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry
  integer ( kind = 4 ) ntryh(4)
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(n)

  save ntryh

  data ntryh / 4, 2, 3, 5 /

      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      fac(nf+2) = ntry
      nl = nq
      if (ntry /= 2) go to 107
      if (nf == 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         fac(ib+2) = fac(ib+1)
  106 continue
      fac(3) = 2
  107 if (nl /= 1) go to 104
      fac(1) = n
      fac(2) = nf
    !  tpi = 8.0D+00 * atan ( 1.0D+00 )
       tpi= 6.28318530717959d0
      argh = tpi / real ( n, kind = 8 )
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 == 0) return
      do 110 k1=1,nfm1
         ip = fac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = real ( ld, kind = 8 ) * argh
            fi = 0.0E+00
            do 108 ii=3,ido,2
               i = i+2
               fi = fi + 1.0E+00
               arg = fi*argld
	       wa(i-1) = cos ( arg )
	       wa(i) = sin ( arg )
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue

  return
end subroutine rffti1
subroutine rfftmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, &
  ier )

!*****************************************************************************80
!
!! RFFTMB: real single precision backward FFT, 1D, multiple vectors.
!
!  Discussion:
!
!    RFFTMB computes the one-dimensional Fourier transform of multiple 
!    periodic sequences within a real array.  This transform is referred 
!    to as the backward transform or Fourier synthesis, transforming the
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to RFFTMB followed
!    by a call to RFFTMF (or vice-versa) reproduces the original
!    array  within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in
!    array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), real array containing LOT 
!    sequences, each having length N.  R can have any number of dimensions, 
!    but the total number of locations must be at least LENR.  On input, the
!    spectral data to be transformed, on output the physical data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array. 
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFTMI before the first call to routine RFFTMF 
!    or RFFTMB for a given transform length N.  
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array. 
!    LENSAV must  be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC, JUMP, N, LOT are not consistent.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(lenr)
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
!  logical xercon

  ier = 0

  if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfftmb ', 6)
    return
  else if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfftmb ', 8)
    return
  else if (lenwrk < lot*n) then
    ier = 3
    call xerfft ('rfftmb ', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('rfftmb ', -1)
    return
  end if

  if (n == 1) then
    return
  end if

  call mrftb1 (lot,jump,n,inc,r,work,wsave,wsave(n+1))

  return
end subroutine rfftmb
subroutine rfftmf ( lot, jump, n, inc, r, lenr, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFTMF: real single precision forward FFT, 1D, multiple vectors.
!
!  Discussion:
!
!    RFFTMF computes the one-dimensional Fourier transform of multiple 
!    periodic sequences within a real array.  This transform is referred 
!    to as the forward transform or Fourier analysis, transforming the 
!    sequences from physical to spectral space.
!
!    This transform is normalized since a call to RFFTMF followed
!    by a call to RFFTMB (or vice-versa) reproduces the original array
!    within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed 
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in 
!    array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), real array containing LOT 
!    sequences, each having length N.  R can have any number of dimensions, but 
!    the total number of locations must be at least LENR.  On input, the
!    physical data to be transformed, on output the spectral data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFTMI before the first call to routine RFFTMF 
!    or RFFTMB for a given transform length N.  
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC, JUMP, N, LOT are not consistent.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(lenr)
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
!  logical xercon

  ier = 0

  if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfftmf ', 6)
    return
  else if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfftmf ', 8)
    return
  else if (lenwrk < lot*n) then
    ier = 3
    call xerfft ('rfftmf ', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('rfftmf ', -1)
    return
  end if

  if (n == 1) then
    return
  end if

  call mrftf1 (lot,jump,n,inc,r,work,wsave,wsave(n+1))

  return
end subroutine rfftmf
subroutine rfftmi ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! RFFTMI: initialization for RFFTMB and RFFTMF.
!
!  Discussion:
!
!    RFFTMI initializes array WSAVE for use in its companion routines 
!    RFFTMB and RFFTMF.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), work array containing the prime 
!    factors of N and also containing certain trigonometric 
!    values which will be used in routines RFFTMB or RFFTMF.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfftmi ', 3)
    return
  end if

  if (n == 1) then
    return
  end if

  call mrfti1 (n,wsave(1),wsave(n+1))

  return
end subroutine rfftmi
subroutine sinq1b ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! SINQ1B: real single precision backward sine quarter wave transform, 1D.
!
!  Discussion:
!
!    SINQ1B computes the one-dimensional Fourier transform of a sequence 
!    which is a sine series with odd wave numbers.  This transform is 
!    referred to as the backward transform or Fourier synthesis, 
!    transforming the sequence from spectral to physical space.
!
!    This transform is normalized since a call to SINQ1B followed
!    by a call to SINQ1F (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, the sequence to be 
!    transformed.  On output, the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINQ1I before the first call to routine SINQ1F 
!    or SINQ1B for a given transform length N.  WSAVE's contents may be 
!    re-used for subsequent calls to SINQ1F and SINQ1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, the error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xhold

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('sinq1b', 6)
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('sinq1b', 8)
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('sinq1b', 10)
  end if

  if ( 1 < n ) go to 101
!
!     x(1,1) = 4.*x(1,1) line disabled by dick valent 08/26/2010
!
      return
  101 ns2 = n/2
      do 102 k=2,n,2
         x(1,k) = -x(1,k)
  102 continue
      call cosq1b (n,inc,x,lenx,wsave,lensav,work,lenwrk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinq1b',-5)
        return
      end if

      do 103 k=1,ns2
         kc = n-k
         xhold = x(1,k)
         x(1,k) = x(1,kc+1)
         x(1,kc+1) = xhold
  103 continue

  return
end subroutine sinq1b
subroutine sinq1f ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! SINQ1F: real single precision forward sine quarter wave transform, 1D.
!
!  Discussion:
! 
!    SINQ1F computes the one-dimensional Fourier transform of a sequence 
!    which is a sine series of odd wave numbers.  This transform is 
!    referred to as the forward transform or Fourier analysis, transforming 
!    the sequence from physical to spectral space.
!
!    This transform is normalized since a call to SINQ1F followed
!    by a call to SINQ1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, the sequence to be 
!    transformed.  On output, the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINQ1I before the first call to routine SINQ1F 
!    or SINQ1B for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINQ1F and SINQ1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, the error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xhold

      ier = 0

      if (lenx < inc*(n-1) + 1) then
        ier = 1
        call xerfft ('sinq1f', 6)
        go to 300
      else if (lensav < 2*n + int(log( real ( n, kind = 8 ) )&
        /log( 2.0E+00 )) +4) then
        ier = 2
        call xerfft ('sinq1f', 8)
        go to 300
      else if (lenwrk < n) then
        ier = 3
        call xerfft ('sinq1f', 10)
        go to 300
      end if

      if (n == 1) return
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(1,k)
         x(1,k) = x(1,kc+1)
         x(1,kc+1) = xhold
  101 continue
      call cosq1f (n,inc,x,lenx,wsave,lensav,work,lenwrk,ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinq1f',-5)
        go to 300
      end if
      do 102 k=2,n,2
         x(1,k) = -x(1,k)
  102 continue
  300 continue

  return
end subroutine sinq1f
subroutine sinq1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! SINQ1I: initialization for SINQ1B and SINQ1F.
!
!  Discussion:
!
!    SINQ1I initializes array WSAVE for use in its companion routines 
!    SINQ1F and SINQ1B.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors 
!    of N and also containing certain trigonometric values which will be used 
!   in routines SINQ1B or SINQ1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

      ier = 0

      if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
        /log( 2.0E+00 )) +4) then
        ier = 2
        call xerfft ('sinq1i', 3)
        go to 300
      end if

      call cosq1i (n, wsave, lensav, ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinq1i',-5)
      end if
  300 continue

  return
end subroutine sinq1i
subroutine sinqmb ( lot, jump, n, inc, x, lenx, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! SINQMB: real single precision backward sine quarter wave, multiple vectors.
!
!  Discussion:
!
!    SINQMB computes the one-dimensional Fourier transform of multiple 
!    sequences within a real array, where each of the sequences is a 
!    sine series with odd wave numbers.  This transform is referred to as 
!    the backward transform or Fourier synthesis, transforming the 
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to SINQMB followed
!    by a call to SINQMF (or vice-versa) reproduces the original
!    array  within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in
!    array R, of the first elements of two consecutive sequences to be 
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing LOT sequences, each 
!    having length N.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINQMI before the first call to routine SINQMF
!    or SINQMB for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINQMF and SINQMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon
  real ( kind = 8 ) xhold

      ier = 0

      if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
        ier = 1
        call xerfft ('sinqmb', 6)
      else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
        /log( 2.0E+00 )) +4) then
        ier = 2
        call xerfft ('sinqmb', 8)
      else if (lenwrk < lot*n) then
        ier = 3
        call xerfft ('sinqmb', 10)
      else if (.not. xercon(inc,jump,n,lot)) then
        ier = 4
        call xerfft ('sinqmb', -1)
      end if

      lj = (lot-1)*jump+1
      if (1 < n ) go to 101

      do m=1,lj,jump
         x(m,1) =  4.0E+00 * x(m,1)
      end do

      return
  101 ns2 = n/2

    do k=2,n,2
       do m=1,lj,jump
         x(m,k) = -x(m,k)
      end do
    end do

      call cosqmb (lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinqmb',-5)
        go to 300
      end if
      do 103 k=1,ns2
         kc = n-k
         do 203 m=1,lj,jump
         xhold = x(m,k)
         x(m,k) = x(m,kc+1)
         x(m,kc+1) = xhold
 203     continue
  103 continue
  300 continue

  return
end subroutine sinqmb
subroutine sinqmf ( lot, jump, n, inc, x, lenx, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! SINQMF: real single precision forward sine quarter wave, multiple vectors.
!
!  Discussion:
!
!    SINQMF computes the one-dimensional Fourier transform of multiple 
!    sequences within a real array, where each sequence is a sine series 
!    with odd wave numbers.  This transform is referred to as the forward
!    transform or Fourier synthesis, transforming the sequences from 
!    spectral to physical space.
!
!    This transform is normalized since a call to SINQMF followed
!    by a call to SINQMB (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations,
!    in array R, of the first elements of two consecutive sequences to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing LOT sequences, each 
!    having length N.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINQMI before the first call to routine SINQMF 
!    or SINQMB for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINQMF and SINQMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon
  real ( kind = 8 ) xhold

      ier = 0

      if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
        ier = 1
        call xerfft ('sinqmf', 6)
        return
      else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
        /log( 2.0E+00 )) +4) then
        ier = 2
        call xerfft ('sinqmf', 8)
        return
      else if (lenwrk < lot*n) then
        ier = 3
        call xerfft ('sinqmf', 10)
        return
      else if (.not. xercon(inc,jump,n,lot)) then
        ier = 4
        call xerfft ('sinqmf', -1)
        return
      end if

      if (n == 1) then
        return
      end if

      ns2 = n/2
      lj = (lot-1)*jump+1
      do 101 k=1,ns2
         kc = n-k
         do m=1,lj,jump
           xhold = x(m,k)
           x(m,k) = x(m,kc+1)
           x(m,kc+1) = xhold
         end do
  101 continue
      call cosqmf (lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinqmf',-5)
        return
      end if

  do k=2,n,2
    do m=1,lj,jump
      x(m,k) = -x(m,k)
    end do
  end do

  300 continue

  return
end subroutine sinqmf
subroutine sinqmi ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! SINQMI: initialization for SINQMB and SINQMF.
!
!  Discussion:
!
!    SINQMI initializes array WSAVE for use in its companion routines 
!    SINQMF and SINQMB.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors 
!    of N and also containing certain trigonometric values which will be used 
!    in routines SINQMB or SINQMF.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('sinqmi', 3)
    return
  end if

  call cosqmi (n, wsave, lensav, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('sinqmi',-5)
  end if

  return
end subroutine sinqmi
subroutine sint1b ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! SINT1B: real single precision backward sine transform, 1D.
!
!  Discussion:
!
!    SINT1B computes the one-dimensional Fourier transform of an odd 
!    sequence within a real array.  This transform is referred to as 
!    the backward transform or Fourier synthesis, transforming the 
!    sequence from spectral to physical space.
!
!    This transform is normalized since a call to SINT1B followed
!    by a call to SINT1F (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N+1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, contains the sequence 
!    to be transformed, and on output, the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINT1I before the first call to routine SINT1F 
!    or SINT1B for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINT1F and SINT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N+2.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('sint1b', 6)
    return
  else if ( lensav < n / 2 + n + int ( log ( real ( n, kind = 8 ) ) &
    / log ( 2.0E+00 ) ) + 4 ) then
    ier = 2
    call xerfft ('sint1b', 8)
    return
  else if (lenwrk < (2*n+2)) then
    ier = 3
    call xerfft ('sint1b', 10)
    return
  end if

  call sintb1(n,inc,x,wsave,work,work(n+2),ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('sint1b',-5)
  end if

  return
end subroutine sint1b
subroutine sint1f ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! SINT1F: real single precision forward sine transform, 1D.
!
!  Discussion:
!
!    SINT1F computes the one-dimensional Fourier transform of an odd 
!    sequence within a real array.  This transform is referred to as the 
!    forward transform or Fourier analysis, transforming the sequence
!    from physical to spectral space.
!
!    This transform is normalized since a call to SINT1F followed
!    by a call to SINT1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N+1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, contains the sequence 
!    to be transformed, and on output, the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINT1I before the first call to routine SINT1F 
!    or SINT1B for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to SINT1F and SINT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N+2.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('sint1f', 6)
    return
  else if (lensav < n/2 + n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('sint1f', 8)
    return
  else if (lenwrk < (2*n+2)) then
    ier = 3
    call xerfft ('sint1f', 10)
    return
  end if

  call sintf1(n,inc,x,wsave,work,work(n+2),ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('sint1f',-5)
  end if

  return
end subroutine sint1f 
subroutine sint1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! SINT1I: initialization for SINT1B and SINT1F.
!
!  Discussion:
!
!    SINT1I initializes array WSAVE for use in its companion routines 
!    SINT1F and SINT1B.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N+1 is a product 
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors 
!    of N and also containing certain trigonometric values which will be used 
!    in routines SINT1B or SINT1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) pi
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < n/2 + n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('sint1i', 3)
    return
  end if

 ! pi =  4.0E+00 * atan ( 1.0E+00 )
  pi= 3.141592653589793d0
  if (n <= 1) then
    return
  end if

  ns2 = n/2
  np1 = n+1
  dt = pi / real ( np1, kind = 8 )

  do k=1,ns2
    wsave(k) =  2.0E+00 *sin(k*dt)
  end do

  lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0E+00 )) +4

  call rfft1i (np1, wsave(ns2+1), lnsv, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('sint1i',-5)
  end if

  return
end subroutine sint1i
subroutine sintb1 ( n, inc, x, wsave, xh, work, ier )

!*****************************************************************************80
!
!! SINTB1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum
  real ( kind = 8 ) fnp1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lnxh
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) srt3s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xh(*)
  real ( kind = 8 ) xhold

      ier = 0
      if (n-2) 200,102,103
  102 srt3s2 = sqrt( dble(3.0) ) / dble(2.0)
      xhold = srt3s2*(x(1,1)+x(1,2))
      x(1,2) = srt3s2*(x(1,1)-x(1,2))
      x(1,1) = xhold
      return

  103 np1 = n+1
      ns2 = n/2
      do 104 k=1,ns2
         kc = np1-k
         t1 = x(1,k)-x(1,kc)
         t2 = wsave(k)*(x(1,k)+x(1,kc))
         xh(k+1) = t1+t2
         xh(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn == 0) go to 124
      xh(ns2+2) =  4.0E+00 * x(1,ns2+1)
  124 xh(1) = 0.0E+00
      lnxh = np1
      lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = np1

      call rfft1f(np1,1,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,ier1)    
 
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sintb1',-5)
        return
      end if

      if(mod(np1,2) /= 0) go to 30
      xh(np1) = xh(np1)+xh(np1)
 30   fnp1s4 = real ( np1, kind = 8 ) / 4.0E+00
         x(1,1) = fnp1s4*xh(1)
         dsum = x(1,1)
      do i=3,n,2
            x(1,i-1) = fnp1s4*xh(i)
            dsum = dsum+fnp1s4*xh(i-1)
            x(1,i) = dsum
      end do

      if ( modn == 0 ) then
         x(1,n) = fnp1s4*xh(n+1)
      end if

  200 continue

  return
end subroutine sintb1
subroutine sintf1 ( n, inc, x, wsave, xh, work, ier )

!*****************************************************************************80
!
!! SINTF1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) inc

  real ( kind = 8 ) dsum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lnxh
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) sfnp1
  real ( kind = 8 ) ssqrt3
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xh(*)
  real ( kind = 8 ) xhold

      ier = 0
      if (n-2) 200,102,103
  102 ssqrt3 = 1.0E+00 / sqrt ( dble(3.0) )
      xhold = ssqrt3*(x(1,1)+x(1,2))
      x(1,2) = ssqrt3*(x(1,1)-x(1,2))
      x(1,1) = xhold
      go to 200
  103 np1 = n+1
      ns2 = n/2
      do k=1,ns2
         kc = np1-k
         t1 = x(1,k)-x(1,kc)
         t2 = wsave(k)*(x(1,k)+x(1,kc))
         xh(k+1) = t1+t2
         xh(kc+1) = t2-t1
      end do

      modn = mod(n,2)
      if (modn == 0) go to 124
      xh(ns2+2) =  4.0E+00 * x(1,ns2+1)
  124 xh(1) = 0.0E+00
      lnxh = np1
      lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0E+00 )) + 4
      lnwk = np1

      call rfft1f(np1,1,xh,lnxh,wsave(ns2+1),lnsv,work, lnwk,ier1)     
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sintf1',-5)
        go to 200
      end if

      if(mod(np1,2) /= 0) go to 30
      xh(np1) = xh(np1)+xh(np1)
   30 sfnp1 = 1.0E+00 / real ( np1, kind = 8 )
         x(1,1) = 0.5E+00 * xh(1)
         dsum = x(1,1)

      do i=3,n,2
            x(1,i-1) = 0.5E+00 * xh(i)
            dsum = dsum + 0.5E+00 * xh(i-1)
            x(1,i) = dsum
      end do

      if (modn /= 0) go to 200
      x(1,n) = 0.5E+00 * xh(n+1)
  200 continue

  return
end subroutine sintf1
subroutine sintmb ( lot, jump, n, inc, x, lenx, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! SINTMB: real single precision backward sine transform, multiple vectors.
!
!  Discussion:
!
!    SINTMB computes the one-dimensional Fourier transform of multiple 
!    odd sequences within a real array.  This transform is referred to as 
!    the backward transform or Fourier synthesis, transforming the 
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to SINTMB followed
!    by a call to SINTMF (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed
!    within the array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in 
!    array R, of the first elements of two consecutive sequences.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N+1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing LOT sequences, each 
!    having length N.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINTMI before the first call to routine SINTMF 
!    or SINTMB for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINTMF and SINTMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*(2*N+4).
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) iw2
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon

  ier = 0

  if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
    ier = 1
    call xerfft ('sintmb', 6)
    return
  else if ( lensav < n / 2 + n + int ( log ( real ( n, kind = 8 ) ) &
    /log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('sintmb', 8)
    return
  else if (lenwrk < lot*(2*n+4)) then
    ier = 3
    call xerfft ('sintmb', 10)
    return
  else if (.not. xercon(inc,jump,n,lot)) then
    ier = 4
    call xerfft ('sintmb', -1)
    return
  end if

  iw1 = lot+lot+1
  iw2 = iw1+lot*(n+1)

  call msntb1(lot,jump,n,inc,x,wsave,work,work(iw1),work(iw2),ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('sintmb',-5)
    return
  end if

  return
end subroutine sintmb
subroutine sintmf ( lot, jump, n, inc, x, lenx, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! SINTMF: real single precision forward sine transform, multiple vectors.
!
!  Discussion:
!
!    SINTMF computes the one-dimensional Fourier transform of multiple 
!    odd sequences within a real array.  This transform is referred to as 
!    the forward transform or Fourier analysis, transforming the sequences 
!    from physical to spectral space.
!
!    This transform is normalized since a call to SINTMF followed
!    by a call to SINTMB (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be 
!    transformed within.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, 
!    in array R, of the first elements of two consecutive sequences.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N+1 is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing LOT sequences, each 
!    having length N.  R can have any number of dimensions, but the total 
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to SINTMI before the first call to routine SINTMF
!    or SINTMB for a given transform length N.  WSAVE's contents may be re-used 
!    for subsequent calls to SINTMF and SINTMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least LOT*(2*N+4).
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) iw2
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
!  logical xercon

  ier = 0

  if ( lenx < ( lot - 1) * jump + inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'sintmf', 6 )
    return
  else if ( lensav < n / 2 + n + int ( log ( real ( n, kind = 8 ) ) &
    / log ( 2.0E+00 ) ) + 4 ) then
    ier = 2
    call xerfft ( 'sintmf', 8 )
    return
  else if ( lenwrk < lot * ( 2 * n + 4 ) ) then
    ier = 3
    call xerfft ( 'sintmf', 10 )
    return
  else if ( .not. xercon ( inc, jump, n, lot ) ) then
    ier = 4
    call xerfft ( 'sintmf', -1 )
    return
  end if

  iw1 = lot + lot + 1
  iw2 = iw1 + lot * ( n + 1 )
  call msntf1 ( lot, jump, n, inc, x, wsave, work, work(iw1), work(iw2), ier1 )

  if ( ier1 /= 0 ) then
    ier = 20
    call xerfft ( 'sintmf', -5 )
  end if

  return
end subroutine sintmf
subroutine sintmi ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! SINTMI: initialization for SINTMB and SINTMF.
!
!  Discussion:
!
!    SINTMI initializes array WSAVE for use in its companion routines 
!    SINTMF and SINTMB.  The prime factorization of N together with a 
!    tabulation of the trigonometric functions are computed and stored 
!    in array WSAVE.  Separate WSAVE arrays are required for different 
!    values of N.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors
!    of N and also containing certain trigonometric values which will be used
!    in routines SINTMB or SINTMF.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  real ( kind = 8 ) dt
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) pi
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if ( lensav < n / 2 + n + int ( log ( real ( n, kind = 8 ) ) &
    / log ( 2.0E+00 ) ) + 4 ) then
    ier = 2
    call xerfft ( 'sintmi', 3 )
    return
  end if

 ! pi = 4.0E+00 * atan ( 1.0E+00 )
   pi= 3.141592653589793d0
  if ( n <= 1 ) then
    return
  end if

  ns2 = n / 2
  np1 = n + 1
  dt = pi / real ( np1, kind = 8 )

  do k = 1, ns2
    wsave(k) = 2.0E+00 * sin ( k * dt )
  end do

  lnsv = np1 + int ( log ( real ( np1, kind = 8 ) ) / log ( 2.0E+00 ) ) + 4
  call rfftmi ( np1, wsave(ns2+1), lnsv, ier1 )

  if ( ier1 /= 0 ) then
    ier = 20
    call xerfft ( 'sintmi', -5 )
  end if

  return
end subroutine sintmi
subroutine w2r ( ldr, ldw, l, m, r, w )

!*****************************************************************************80
!
!! W2R copies a 2D array, allowing for different leading dimensions.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldw
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real ( kind = 8 ) r(ldr,m)
  real ( kind = 8 ) w(ldw,m)

  do j = 1, m
    do i = 1, l
      r(i,j) = w(i,j)
    end do
  end do

  return
end subroutine w2r

function xercon( inc, jump, n, lot )

!*****************************************************************************80
!
!! XERCON checks INC, JUMP, N and LOT for consistency.
!
!  Discussion:
!
!    Positive integers INC, JUMP, N and LOT are "consistent" if,
!    for any values I1 and I2 < N, and J1 and J2 < LOT,
!
!      I1 * INC + J1 * JUMP = I2 * INC + J2 * JUMP
!
!    can only occur if I1 = I2 and J1 = J2.
!
!    For multiple FFT's to execute correctly, INC, JUMP, N and LOT must
!    be consistent, or else at least one array element will be
!    transformed more than once.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INC, JUMP, N, LOT, the parameters to check.
!
!    Output, logical XERCON, is TRUE if the parameters are consistent.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnew
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lcm
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  logical xercon

  i = inc
  j = jump

  do while ( j /= 0 )
    jnew = mod ( i, j )
    i = j
    j = jnew
  end do
!
!  LCM = least common multiple of INC and JUMP.
!
  lcm = ( inc * jump ) / i

  if ( lcm <= ( n - 1 ) * inc .and. lcm <= ( lot - 1 ) * jump ) then
    xercon = .false.
  else
    xercon = .true.
  end if

  return
end function xercon

subroutine xerfft ( srname, info )

!*****************************************************************************80
!
!! XERFFT is an error handler for the FFTPACK routines.
!
!  Discussion:
!
!    XERFFT is an error handler for FFTPACK version 5.1 routines.
!    It is called by an FFTPACK 5.1 routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!    Installers may consider modifying the stop statement in order to
!    call system-specific exception-handling facilities.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, character ( len = * ) SRNAME, the name of the calling routine.
!
!    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid 
!    parameter in the parameter list of the calling routine has been detected, 
!    INFO is the position of that parameter.  In the case when an illegal 
!    combination of LOT, JUMP, N, and INC has been detected, the calling 
!    subprogram calls XERFFT with INFO = -1.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)') '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write( *, '(a,a,a,a)') '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end subroutine xerfft


    end module FFTclass
