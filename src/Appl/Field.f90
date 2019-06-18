!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! FieldQuantclass: 3D field quantity class in Field module of APPLICATION
!                 layer.
! Version: 2.0
! Author: Ji Qiang, LBNL
! Description: This class defines a 3-D field quantity in the accelerator.
!              The field quantity can be updated at each step. 
! Comments:
!----------------------------------------------------------------
      module FieldQuantclass
        use Timerclass
        use CompDomclass
        use Pgrid2dclass
        use FFTclass
        use Transposeclass
        use PhysConstclass
        use Dataclass
        type FieldQuant
!          private
          !# of mesh points in x and y directions.
          integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal
          !# field quantity array.
          double precision, pointer, dimension(:,:,:) :: FieldQ
        end type FieldQuant
      contains
        !Initialize field class.
        subroutine construct_FieldQuant(this,innx,inny,innz,geom,grid) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz 
        integer :: myid, myidx, myidy, nptot,nproccol,nprocrow
        integer, allocatable, dimension(:,:,:) :: LocalTable

        call getsize_Pgrid2d(grid,nptot,nproccol,nprocrow) 
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        allocate(LocalTable(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(nproccol.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprocrow.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
  
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0

        deallocate(LocalTable)

        end subroutine construct_FieldQuant
   
        ! set field quantity.
        subroutine set_FieldQuant(this,innx,inny,innz,geom,grid,&
                                  nprx,npry) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz, nprx, npry 
        integer :: myid, myidx, myidy
        integer, dimension(2,0:nprx-1,0:npry-1)::LocalTable

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(npry.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprx.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
        this%Nxlocal = innx
  
        deallocate(this%FieldQ) 
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal)) 
        this%FieldQ = 0.0

        end subroutine set_FieldQuant

!----------------------------------------------------------------------
! update potential (solving Possion's equation) with 3D isolated 
! boundary conditions.
!
        subroutine update3O_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (FieldQuant), intent(inout) :: this
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0
! Sherry
        double precision :: t1, t_openbc3d
        integer :: iam


        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

! Sherry
!        call starttime_Timer(t1)

        ! Open boundary conditions!
        call openBC3D(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)
! Sherry
!        call MPI_COMM_RANK(comm2d,iam,ierr)
!        t_openbc3d = elapsedtime_Timer(t1)
!        if ( iam == 0 ) then
!           print*, 'openBC3D time = ', t_openbc3d
!        endif

        end subroutine update3O_FieldQuant


        ! Solving Poisson's equation with open BCs.
        ! Sherry: now, innx, inny, innz are local dimensions
        subroutine openBC3D(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr
! Sherry
        double precision :: t1, t_fft
        integer :: iam

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

! Sherry
!        call starttime_Timer(t1)

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

! Sherry
!        call MPI_COMM_RANK(comm2d,iam,ierr)
!        t_fft = elapsedtime_Timer(t1)
!        if ( iam == 0 ) then
!           print*, 'fft3d1_FFT time = ', t_fft
!        endif

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        call greenf1tIntnew(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green 
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0d0/dble(n1)
        scaley = 1.0d0/dble(n2)
        scalez = 1.0d0/dble(n3)
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        deallocate(rho2out)

        end subroutine openBC3D

!----------------------------------------------------------------------
        ! green function for extended array.
        subroutine greenf1tIntnew(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision, dimension (nx+2,nsizey+1,nsizez+1) :: grntmp
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: rr,aa,bb,cc,dd,ee,ff,ss
        double complex :: gg,gg2
        double complex :: ggrr
        double precision, dimension(2) :: xx,yy,zz
        double precision, dimension(3) :: vv
        integer :: n,i0,j0,k0
        double precision :: recfourpi

        recfourpi = 1.0d0/(8.0d0*asin(1.0d0))

        call starttime_Timer(t0)

!        if(myidx.eq.0 .and. myidy.eq.0) then
!          print*,"into integrated Green function....."
!        endif
        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k0 = 1, nsizez+1
          do j0 = 1, nsizey+1
            do i0 = 1, nx+2
              jj = j0 + nblocky
              kk = k0 + nblockx
              kkk = kk - 1
              jjj = jj - 1
              iii = i0 - 1

              vv(1) = iii*hx-hx/2
              vv(2) = jjj*hy-hy/2
              vv(3) = kkk*hz-hz/2
            
                rr = sqrt(vv(1)**2+vv(2)**2+vv(3)**2)
                aa = vv(1)**2*vv(3)+(vv(2)**2+vv(3)**2)*vv(3) + &
                            vv(1)*vv(3)*rr
                bb = vv(1)**2*vv(2) + vv(1)*vv(2)*rr
                cc = vv(1)*(vv(1)**2+vv(2)**2)+vv(3)*vv(1)*(vv(3)+rr)
                dd = vv(3)*vv(2)*(vv(3) + rr)
                ee = vv(1)**2*vv(2)+vv(2)*(vv(2)**2+vv(3)*(vv(3)+rr))
                ff = vv(1)*vv(3)*(vv(3)+rr)
                ss = 4*vv(2)*vv(3)*log(vv(1)+rr) + &
                        4*vv(1)*vv(3)*log(vv(2)+rr) + &
                        4*vv(1)*vv(2)*log(vv(3)+rr)

                gg2 = cmplx(0.0,vv(3)**2)*log(cmplx(aa**2-bb**2,2*aa*bb)/ &
                        (aa**2+bb**2) ) + cmplx(0.0,vv(1)**2)*&
                        log(cmplx(cc**2-dd**2,2*cc*dd )/(cc**2+dd**2))+&
                        cmplx(0.0,vv(2)**2)*log(cmplx(ee**2-ff**2,2*ee*ff)/ &
                        (ee**2+ff**2) ) + ss
               
              !grntmp(i0,j0,k0) = recfourpi*real(gg2)/(4*hx*hy*hz) !wrong in Z code
              grntmp(i0,j0,k0) = real(gg2)/(4*hx*hy*hz)

            enddo
          enddo
        enddo

        do k0 = 1, nsizez
          do j0 = 1, nsizey
            do i0 = 1, nx+1
              grn(i0,j0,k0) = grntmp(i0+1,j0+1,k0+1)-grntmp(i0,j0+1,k0+1)-grntmp(i0+1,j0,k0+1)+&
                              grntmp(i0,j0,k0+1)-grntmp(i0+1,j0+1,k0)+grntmp(i0,j0+1,k0)+&
                              grntmp(i0+1,j0,k0)-grntmp(i0,j0,k0)
            enddo
          enddo
        enddo

!              xx(1) = iii*hx-hx/2
!              xx(2) = iii*hx+hx/2
!              yy(1) = jjj*hy-hy/2
!              yy(2) = jjj*hy+hy/2
!              zz(1) = kkk*hz-hz/2
!              zz(2) = kkk*hz+hz/2
!       
!              !find the integrated Green function.
!              n = 0
!              do k = 1, 2
!                do j = 1, 2
!                  do i = 1, 2
!                    n = n+1
!                    rr(n) = sqrt(xx(i)**2+yy(j)**2+zz(k)**2)
!                    aa(n) = xx(i)**2*zz(k)+(yy(j)**2+zz(k)**2)*zz(k) + &
!                            xx(i)*zz(k)*rr(n)
!                    bb(n) = xx(i)**2*yy(j) + xx(i)*yy(j)*rr(n)
!                    cc(n) = xx(i)*(xx(i)**2+yy(j)**2)+zz(k)*xx(i)*(zz(k)+rr(n))
!                    dd(n) = zz(k)*yy(j)*(zz(k) + rr(n))
!                    ee(n) = xx(i)**2*yy(j)+yy(j)*(yy(j)**2+zz(k)*(zz(k)+rr(n)))
!                    ff(n) = xx(i)*zz(k)*(zz(k)+rr(n))
!                    ss(n) = 4*yy(j)*zz(k)*log(xx(i)+rr(n)) + &
!                            4*xx(i)*zz(k)*log(yy(j)+rr(n)) + &
!                            4*xx(i)*yy(j)*log(zz(k)+rr(n))
!                    gg(n) = cmplx(0.0,zz(k)**2)*log(cmplx(aa(n)**2-bb(n)**2,2*aa(n)*bb(n))/ &
!                            (aa(n)**2+bb(n)**2) ) + cmplx(0.0,xx(i)**2)*&
!                            log(cmplx(cc(n)**2-dd(n)**2,2*cc(n)*dd(n) )/(cc(n)**2+dd(n)**2))+&
!                            cmplx(0.0,yy(j)**2)*log(cmplx(ee(n)**2-ff(n)**2,2*ee(n)*ff(n))/ &
!                            (ee(n)**2+ff(n)**2) )
!                    gg2(n) = ss(n) +  gg(n)
!                  enddo
!                enddo
!              enddo
!              ggrr = (-gg2(1)+gg2(2)+gg2(3)-gg2(4)+gg2(5d0)-gg2(6)-gg2(7)+gg2(8))/4
!              grn(i0,j0,k0) = real(ggrr)/(hx*hy*hz)
!            enddo
!          enddo
!        enddo


        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0d0
        scaley = 1.0d0
        scalez = 1.0d0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1tIntnew
!--------------------------------------------------------------------

!-------------------------------------------------------------------------
        subroutine setval_FieldQuant(this,i,j,k,value)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        integer, intent(in) :: i, j, k
        double precision, intent(in) :: value

        this%FieldQ(i,j,k) = value

        end subroutine setval_FieldQuant

        double precision function get_FieldQuant(this,i,j,k)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(in) :: i, j, k

        get_FieldQuant = this%FieldQ(i,j,k)

        end function get_FieldQuant

        subroutine getglb_FieldQuant(this,temp)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        type (FieldQuant), intent(out) :: temp
        integer :: i, j, k, lcnz,lcny,lcnx
        double precision :: value

        lcnx = this%Nxlocal
        lcny = this%Nylocal
        lcnz = this%Nzlocal
    
        do k = 1, lcnz
          do j = 1, lcny
            do i = 1, lcnx
              value = get_FieldQuant(this,i,j,k)
              call setval_FieldQuant(temp,i,j,k,value)
            enddo
          enddo 
        enddo

        end subroutine getglb_FieldQuant

        subroutine getlcgrid_FieldQuant(this,nxlc,nylc,nzlc)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(out) :: nxlc,nylc,nzlc

        nxlc = this%Nxlocal
        nylc = this%Nylocal
        nzlc = this%Nzlocal

        end subroutine

        subroutine destruct_FieldQuant(this)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this

        deallocate(this%FieldQ) 

        end subroutine destruct_FieldQuant

!-------------------------------------------------------------------------
       !longitudinal and transverse wakefield
       subroutine wakefield_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                                       hz,aa,gg,leng)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz
       double precision, intent(in) :: hz, aa, gg, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i
       double precision :: scale,zz00,pilc,Z0,alpha1,alpha,zz
       double precision :: coef1,coef2,densconst,offset,zzmax
       real*8 :: tmptmp

       pilc = 2*asin(1.0)
  
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       !longitudinal wakefield function
       !the following longitudinal wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003.
       pilc = 2*asin(1.0d0)
       alpha1 = 0.4648
       alpha = 1.-alpha1*sqrt(gg/leng)-(1.-2*alpha1)*(gg/leng)
       !slac wake
       !zz00 = gg*(aa/(alpha*leng))**2/8 
       !fermi wake
       zz00 = 0.41*(aa/leng)**0.8*(gg/leng)**1.6*aa
       Z0 = 120*pilc
       !greenwake(1,1) = 0.0
       zz = 0.0d0
       !The 0.5 is from S+
       greenwake(1,1) = 0.5d0*Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
       !do kz = 1, Nz+1
       do kz = 2, Nz+1
         zz = (kz-1)*hz
         greenwake(kz,1) = Z0*Clight*exp(-sqrt(zz/zz00))/(pilc*aa*aa)
       enddo
       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       !do kz = 1, twonz
       !  greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       !enddo
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz
       enddo
       
!------------------------------------------------------------------------------
       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       !tranverse wakefield
       !the following transverse wakefield function is from
       !"Short-range dipole wakefields in accelerating structures for the NLC",
       !K. L. F. Bane, SLAC-PUB-9663, 2003. here, 1.1 is an average as suggested       !on page 11 for cell 45.
       !for LCLS slac
       !zz00 = 1.1*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       !for Fermi Elettra
       zz00 = 1.0*0.169*aa*(aa/leng)**1.17*(gg/aa)**0.38
       coef1 = 4*Z0*Clight*zz00/(pilc*aa*aa*aa*aa)
       !densconst = 0.5e-6
       !offset = 0.001
       !zzmax = 0.002
       !coef2 = coef1/zz00*densconst*offset
       !print*,"aa wake: ",aa,leng,gg,zz00,Z0,coef1,coef2,offset,densconst
       greenwake(1,1) = 0.0
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         greenwake(kz,1) = coef1*(1.0-(1.+sqrt(zz/zz00))*exp(-sqrt(zz/zz00)))
       enddo

       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)

       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz
       enddo

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)
 
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz
       enddo

       end subroutine wakefield_FieldQuant

!-------------------------------------------------------------------------
       !longitudinal and transverse wakefield using read in discrete data
       subroutine wakefieldread_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                                       hz,aa,gg,leng,flagbtw)
       implicit none
       include 'mpif.h'
       integer, intent(in) :: Nz,flagbtw
       double precision, intent(in) :: hz, aa, gg, leng
       double precision, dimension(Nz,2), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(in) :: xwakez,ywakez
       double precision, dimension(Nz), intent(out) :: exwake,ezwake,eywake
       double precision, dimension(2*Nz,1) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i,iwk,iwk1
       double precision :: scale,zz00,pilc,Z0,alpha1,alpha,zz
       double precision :: coef1,coef2,densconst,offset,zzmax
       real*8 :: tmptmp,pbtw1,pbtw2,pbtw3,pbtw4
       real*8 :: leng1,leng2,hwk,zwk,ab

       pilc = 2*dasin(1.0d0)
  

       if(flagbtw.eq.4) then !skip longitudinal wake
         ezwake = 0.0d0
         goto 50
       endif

       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid 
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(edatwk(iwk)*ab + edatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = edatwk(iwk)*ab + edatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       
       do kz = Nz+2, twonz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 2, Nz
         greenwake(kz,1) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       !do kz = 1, twonz
       !  greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       !enddo
       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz,1)*hz
       enddo

50     continue
       
       if((flagbtw.eq.2) .or. (flagbtw.eq.3)) then !no transverse wake for standing wave cavity so far
         exwake = 0.0d0
         eywake = 0.0d0
         goto 100
       endif

!------------------------------------------------------------------------------
       !calculate the transverse wakefield effects
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*xwakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo
 
       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(epdatwk(iwk)*ab + epdatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = epdatwk(iwk)*ab + epdatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   densz2nout)

       do kz = 1, Nz
         exwake(kz) = densz2nout(kz,1)*hz
       enddo

       !vertical "Y" wakefields
       do kz = 1, Nz
          densz2n(kz,1) = recvdensz(kz,1)*ywakez(kz)
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz,1) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       hwk = zdatwk(2)-zdatwk(1) !assumed uniform grid
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         zwk = zz - zdatwk(1)
         iwk = zwk/hwk + 1
         iwk1 = iwk+1
         ab = iwk*hwk-zwk
         if(kz.eq.1) then
           greenwake(kz,1) = 0.5d0*(eppdatwk(iwk)*ab + eppdatwk(iwk1)*(1.0d0-ab))
         else
           greenwake(kz,1) = eppdatwk(iwk)*ab + eppdatwk(iwk1)*(1.0d0-ab)
         endif
       enddo
       do kz = Nz+2, 2*Nz
         greenwake(kz,1) = greenwake(twonz-kz+2,1)
       enddo
       do kz = 1, Nz
         greenwake(kz,1) = 0.0
       enddo

       ksign = 1
       !scale = 1./(twonz)
       scale = 1.0d0
       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
                greenwakeout)

       do kz = 1, 2
          greenwake(kz,1) = densz2nout(kz,1)*greenwakeout(kz,1)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz-1,1)-&
                            densz2nout(2*kz,1)*greenwakeout(2*kz,1)
          greenwake(2*kz,1) = densz2nout(2*kz-1,1)*greenwakeout(2*kz,1)+&
                            densz2nout(2*kz,1)*greenwakeout(2*kz-1,1)
       enddo

       !scale = 1.0
       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
                   greenwakeout)

       do kz = 1, Nz
         eywake(kz) = greenwakeout(kz,1)*hz
       enddo

100    continue

       end subroutine wakefieldread_FieldQuant

        !This subroutine calculates the 1d csr wakefield including
        !entrance, stead-state, and transitions effects.
        !Here, hx, rho,...are in real units. The return ezwake is in V/m.
        !The current version uses IGF corresponding to the four cases of Saldin et al.
        subroutine csrwakeTrIGF_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,rhonewp,&
                             rhonewpp,gam,ezwake)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength,gam
        real*8, dimension(Nx) :: rhonew,rhonewp,rhonewpp,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp,xslpN
        integer :: Ni,il,i,Nmax,j,islpN,islpN1
        integer :: myidlc,ierr,islp,islp1,islp0,jstart
        real*8 :: aa,uuh,uuh24,xslp,csrss,xxp,ssh,xk2,tcoef
        real*8 :: bb,cc,yh,phih,phpy2,csrtr,xblg,xxbarh,phimh,psimax
        real*8 :: psip2x2,psipx2,csrdrm1,csrdr1,phpxy2,phpy,csrss1,csrss2,&
                  csrtr1,csrtr2,csrdrm2,csrdr2
!        real*8 :: IcsrCaseA,IcsrCaseB,IcsrCaseC,IcsrCaseD

!        call MPI_COMM_RANK(MPI_COMM_WORLD, myidlc, ierr)

!        epstol = 1.0d-10 !tolerance for root finding
        epstol = 2.0d-9
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817d-12)
        xk2 = gam/r0

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24

           Ni = i
           ezwake(i) = 0.0

           xblg = (i-1)*hx

           !if(myidlc.eq.0) print*,"ii:",i,xxl,xx

           if(xx.le.0) then
           else if(xx.le.blength) then

             xslp = (xx/r0)*r0/2/gam/gam + xxl
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1

             !Case A 
             phih = xx/r0*gam
             ! IGF integral over the interior sample points for Case A:
             do j = 2, islp1-1
                 !write(*,*) 'Inside Case A!'
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrtr1-csrtr2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case A
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             !Add the lower end IGF integral for case A
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrtr1-csrtr2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             endif

             ! Case B (steady-state regime)
             jstart = max(islp1,1)
             ! IGF integral over the interior sample points for Case B:
             do j = jstart+2,Ni-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             enddo
             !add end integrals
             if(jstart.le.Ni-2) then
             !Add the upper end IGF integral for case B
                 xxp = ptmin + (Ni-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (Ni-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(Ni)*(csrss1-csrss2)
             !Add the lower end IGF integral for case B
                 if(islp1.lt.1) then
                   j = 1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)  
                   xxp = ptmin + (j-1)*hx 
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 else
                   j = islp1+1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)
                   xxp = xx-xslp
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 endif
             else if(jstart.eq.Ni-1) then
                 j = Ni
                 xxp = ptmin + (j-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = xx-xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             endif

           !  ezwake(i) = ezwake(i)*xconst

             !if(myidlc.eq.0) print*,"xxl: ",xxl,xx,ezwake(i),r0,ptmin,phim
             !if(myidlc.eq.0) print*,"xxl: ",xx,ezwake(i),xxl,xx-ptmin,xslp,hx

           else

             !if(myidlc.eq.0) print*,"CD:",i,islp,xslp,hx,xx,ezwake(i)

             xxbar = (xx - blength)
             xslp= (r0*phim+xxbar)/2/gam/gam + r0*phim**3/24*&
                   (r0*phim+4*xxbar)/(r0*phim+xxbar)
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1
             xslpN = xxbar/2.d0/gam/gam
             islpN = (xx-xslpN-ptmin)/hx
             islpN1 = islpN + 1

             ! Case C
             xxbarh = xxbar*gam/r0
             phimh = phim*gam
             ! IGF integral over the interior sample points for Case C:
             do j = 2, islp1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdrm1-csrdrm2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case C
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             !Add the lower end IGF integral for case C
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrdrm1-csrdrm2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             endif


             !Case D 
             psimax = phim*gam
             xxbarh = xxbar*gam/r0
             jstart = max(islp1,0)
             !write(13,*) 'xslp and xslpN',xslp,xslpN
             ! IGF integral over the interior sample points for Case D:
             do j = jstart+2,islpN1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdr1-csrdr2)
             enddo
             if(islpN1.ge.(jstart+2)) then
             !Add the upper end IGF integral for case D
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (islpN1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(islpN1)*(csrdr1-csrdr2)
             !Add the lower end IGF integral for case D
                 xxp = ptmin + (jstart)*hx + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             else if(islpN1.eq.(jstart+1)) then
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             endif


           endif
           ezwake(i) = ezwake(i)*xconst
        enddo

        end subroutine csrwakeTrIGF_FieldQuant

           function IcsrCaseA(phih,ssh,xk2)
           implicit none
           double precision:: phih,ssh,xk2,IcsrCaseA
           double precision:: bb,cc,yh,phpy
           bb = 2*phih+phih**3/3-2*ssh
           cc = phih**2+phih**4/12-2*ssh*phih
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpy = (phih+yh)
           IcsrCaseA = xk2*(-(2*phpy+phih**3)/&
                    (phpy**2+phih**4/4)+1.0d0/ssh)
           end function IcsrCaseA


           function IcsrCaseB(ssh,xk2)
           implicit none
           double precision:: ssh,xk2,IcsrCaseB
           double precision:: aa,uuh
           aa = sqrt(64.0+144.0d0*ssh**2)
           uuh = (aa+12*ssh)**(1.0d0/3.0d0)-(aa-12*ssh)**(1.0d0/3.0d0)
           IcsrCaseB = xk2*(-4*uuh*(uuh**2+8)/ &
                    ((uuh**2+4.0d0)*(uuh**2+12.0d0)))
           end function IcsrCaseB


           function IcsrCaseC(phimh,xxbarh,ssh,xk2)
           implicit none
           double precision:: phimh,xxbarh,ssh,xk2,IcsrCaseC
           double precision:: bb,cc,yh,phpxya,phpxyb
           bb = 2*(phimh+xxbarh)-2*ssh+phimh**3/3+phimh**2*xxbarh
           cc = (phimh+xxbarh)**2+phimh**2*(phimh**2+4*phimh*xxbarh)/12-&
                      2*ssh*(phimh+xxbarh)
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpxya = (phimh+xxbarh+yh)
           phpxyb = phimh*xxbarh+phimh*phimh/2.d0
           IcsrCaseC = xk2*(-2.d0*(phpxya+phimh*phpxyb)/ &
                          (phpxya**2+phpxyb**2)+1.0d0/ssh)
           end function IcsrCaseC

           function IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
           implicit none
           double precision:: xxbarh,ssh,xk2,epstol,IcsrCaseD
           double precision:: psi,psipxa,psipxb,psimax,x1,x2
           integer:: Nmax
           x1 = -epstol
           x2 = psimax*(1.d0+epstol)
           call root(x1,x2,epstol,Nmax,ssh,xxbarh,psi)
           psipxa = (psi+xxbarh)
           psipxb = psi*(xxbarh+psi/2.d0)
           IcsrCaseD = xk2*(-2.d0*(psipxa+psi*psipxb)/&
                    (psipxa**2+psipxb**2)+1.0d0/ssh)
           end function IcsrCaseD

     subroutine root(x1,x2,xacc,maxit,zeta,xxh,rtsafe)
!************************************************************************
!  This routine computes the root of the function that is evaluated in 
!  the subroutine 'funcd'. It is based on the subroutine 'root' of 
!  Numerical Recipes 9.4, which makes use of a Newton-Raphson method 
!  with root bracketing.  It has been modified to handle the two bracket 
!  endpoints carefully. The routine searches for a root in the interval
!  [x1,x2] with a tolerance given by 'xacc', and returns this value
!  as 'rtsafe'.  The maximum number of iterations allowed is 'maxit'.
!  C.E.M.
!***********************************************************************
     implicit none
     double precision:: rtsafe,x1,x2,xacc
     double precision:: xxh,zeta
     integer:: j,maxit
     double precision:: df,dx,dxold,f,fh,fl,temp,xh,xl
     call funcd(x1,xxh,zeta,fl,df)
     call funcd(x2,xxh,zeta,fh,df)
     if((fl>0.d0.and.fh>0.d0).or.(fl<0.d0.and.fh<0.d0)) then
           pause 'root must be bracketed in rtsafe'
           write(*,*) 'psimax,fl,fh = ',x2,fl,fh
     endif
     if(dabs(fl)< xacc) then
       rtsafe=x1
       return
     else if(dabs(fh)< xacc) then
       rtsafe=x2
       return
     else if(fl<0.d0) then
       xl=x1
       xh=x2
     else
       xh=x1
       xl=x2
     endif
     rtsafe=0.5d0*(x1+x2)
     dxold=dabs(x2-x1)
     dx=dxold
     call funcd(rtsafe,xxh,zeta,f,df)
     do j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0.d0.or. &
&         dabs(2.d0*f)>dabs(dxold*df)) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl==rtsafe) return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp==rtsafe) return
        endif
        if(abs(dx)<xacc) return
        call funcd(rtsafe,xxh,zeta,f,df)
        if(f<0.d0) then
           xl=rtsafe
        else
           xh=rtsafe
        endif
     enddo
     pause 'root finding exceeding maximum iterations'
     return
     end subroutine

     subroutine funcd(psi,xxh,deltas,f,derivf)
!**********************************************************
!  This routine evaluates the function whose root produces
!  the retarded angle psi that is required for evaluating
!  the CSR kernel in Case D of Saldin et al.  The value
!  of the function is output as 'f', and its derivative
!  is output as 'derivf'.  C.E.M.
!*********************************************************
     implicit none
     double precision:: deltas,psi,term1,xxh,gamma,f,derivf
     double precision:: alpha,kappa,tau,theta
     f = psi**4/12+xxh*psi**3/3+psi**2+(2*xxh-2*deltas)*psi-&
               2*deltas*xxh+xxh*xxh
     derivf = psi**3/3+xxh*psi**2+2*psi+2*xxh-2*deltas
     end subroutine

      end module FieldQuantclass
