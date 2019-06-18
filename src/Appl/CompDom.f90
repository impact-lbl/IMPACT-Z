!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! CompDomclass: 3D global and local computational domain class in 
!               Computational Domain module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines 3-D global and local computational domain
!              in the parallel simulation.
! Comments:
!----------------------------------------------------------------
      module CompDomclass
        use Timerclass
        use Pgrid2dclass
        integer, private, parameter :: Ndim = 3   !spatial dimension
        interface setlctab_CompDom
          module procedure setlctab1_CompDom, setlctab2_CompDom
        end interface
        type CompDom
!          private
          !spatial range in each dimenion.
          double precision, dimension(2*Ndim) :: SpatRange
          double precision, dimension(2*Ndim) :: Sptrnglocal
          !mesh size in each dimension.
          double precision, dimension(Ndim) :: Meshsize
          !# of mesh points in each dimension.
          integer, dimension(Ndim) :: Meshnum
          integer, dimension(Ndim) :: Mshlocal
          !ymin, ymax, # of y grids on local processor.
          double precision, pointer, dimension(:,:,:) :: LcTabrg
          integer, pointer, dimension(:,:,:) :: LcTabnm
        end type CompDom
        interface construct_CompDom
          module procedure init_CompDom
        end interface
      contains
        ! calculate the initial computational geometry parameters
        ! and double precision parameters. Here, computational domain is mapped
        ! onto a one dimension processor array in y direction. 
        subroutine init_CompDom(this,distparam,nparam, &
               flg,nx,ny,nz,grid2d,nprocrow,nproccol,Flagbc,xrad,yrad,perd)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid2d
        type (CompDom), intent(out) :: this
        integer, intent(in) :: flg,nx,ny,nz,nprocrow,nproccol,nparam,Flagbc
        double precision, dimension(nparam), intent(in) :: distparam
        double precision, intent(in) :: xrad,yrad,perd
        integer :: myid,myidx,myidy,comm2d, &
                   commcol, commrow, ierr, totnp,npxx,npyy
        integer :: i,j,ih,il,nsum
        double precision :: sig1x,sig2x,ux,scalex,sig1y,sig2y,uy,&
        scaley,sig1z,sig2z,uz,scalez,scalepx,scalepy,scalepz
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: pi,cons1,cons2,rxsq2,rysq2,rzsq2,xmin,&
                 xmax,ymin,ymax, &
                zmin,zmax,sig11x,sig22x,sig11y,sig22y,sig11z,sig22z
        double precision :: deltaz,deltay,hz,hy
        double precision :: deltaint,sumint,interval1,interval0,tmp
        double precision, dimension(2,0:nprocrow-1) :: lctabrgz
        double precision, dimension(2,0:nproccol-1) :: lctabrgy
        integer :: ip,numint
        integer :: i0,j0,iflag
     
        sig1x = distparam(1)
        sig2x = distparam(2)
        ux = distparam(3)
        scalex = distparam(4)
        scalepx = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sig1y = distparam(8)
        sig2y = distparam(9)
        uy = distparam(10)
        scaley = distparam(11)
        scalepy = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sig1z = distparam(15)
        sig2z = distparam(16)
        uz = distparam(17)
        scalez = distparam(18)
        scalepz = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)
        
        sig11x = sig1x*scalex
        sig22x = sig2x*scalepx
        sig11y = sig1y*scaley
        sig22y = sig2y*scalepy
        sig11z = sig1z*scalez
        sig22z = sig2z*scalepz

        pi = 2.0*asin(1.0)
        if(flg.eq.1) then  !6D uniform
          xmin = xmu1 -sqrt(3.0)*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + sqrt(3.0)*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -sqrt(3.0)*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + sqrt(3.0)*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
        elseif(flg.eq.2) then  !6D Gaussian
          xmin = xmu1 -4.0*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + 4.0*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -4.0*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + 4.0*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -4.0*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + 4.0*sig11z/sqrt(1.0-uz*uz)
        elseif(flg.eq.3) then  !Waterbag
          xmin = xmu1 -sqrt(8.0)*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + sqrt(8.0)*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -sqrt(8.0)*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + sqrt(8.0)*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -sqrt(8.0)*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + sqrt(8.0)*sig11z/sqrt(1.0-uz*uz)
        elseif(flg.eq.4) then  !Semi-Gaussian
          xmin = xmu1 -sqrt(5.0)*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + sqrt(5.0)*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -sqrt(5.0)*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + sqrt(5.0)*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -sqrt(5.0)*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + sqrt(5.0)*sig11z/sqrt(1.0-uz*uz)
        elseif(flg.eq.5) then !3D KV
          xmin = xmu1 -2*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + 2*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -2*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + 2*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
        elseif(flg.eq.6) then !read in
          xmin = xmu1 -2*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + 2*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -2*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + 2*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -2*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + 2*sig11z/sqrt(1.0-uz*uz)
        else
          xmin = xmu1 -sqrt(3.0)*sig11x/sqrt(1.0-ux*ux)
          xmax = xmu1 + sqrt(3.0)*sig11x/sqrt(1.0-ux*ux)
          ymin = xmu3 -sqrt(3.0)*sig11y/sqrt(1.0-uy*uy)
          ymax = xmu3 + sqrt(3.0)*sig11y/sqrt(1.0-uy*uy)
          zmin = xmu5 -sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
          zmax = xmu5 + sqrt(3.0)*sig11z/sqrt(1.0-uz*uz)
          !print*,"wrong option of initial distributions!"
        endif

        if(Flagbc.eq.3) then
          xmin = 0.0
          xmax = xrad
          ymin = 0.0
          ymax = 4*asin(1.0)
        else if(Flagbc.eq.4) then
          xmin = 0.0
          xmax = xrad
          ymin = 0.0
          ymax = 4*asin(1.0)
          zmin = -perd/2
          zmax = perd/2
        endif
        
        this%SpatRange(1) = xmin
        this%SpatRange(2) = xmax
        this%SpatRange(3) = ymin
        this%SpatRange(4) = ymax
        this%SpatRange(5) = zmin
        this%SpatRange(6) = zmax

        this%Meshnum(1) = nx
        this%Meshnum(2) = ny
        this%Meshnum(3) = nz

        this%Meshsize(1) = (xmax-xmin)/float(nx-1)
        this%Meshsize(2) = (ymax-ymin)/float(ny-1)
        this%Meshsize(3) = (zmax-zmin)/float(nz-1)

        call getsize_Pgrid2d(grid2d,totnp,npyy,npxx)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)

        allocate(this%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
        allocate(this%LcTabnm(2,0:nprocrow-1,0:nproccol-1))

        deltaz = (zmax-zmin)/nprocrow
        deltay = (ymax-ymin)/nproccol
        hz = this%Meshsize(3)
        hy = this%Meshsize(2)
   
        if(nprocrow.gt.1) then

        if(flg.eq.2) then
          numint = 50001

          deltaint = (zmax-zmin)/(numint-1)
          cons1 = sqrt(1.0-uz*uz)/sig11z/sqrt(2.0*pi)
          cons2 = -0.5*(1.0-uz*uz)/sig11z/sig11z

          sumint = 0.0
          do i=1, numint-1
            tmp = zmin + (i-1)*deltaint
            sumint = sumint + deltaint*cons1*exp(cons2*tmp**2)  
          enddo

          interval0 = sumint/float(nprocrow)
          interval1 = interval0
          sumint = 0.0
          do i=1, numint-1
            tmp = zmin + (i-1)*deltaint
            sumint = sumint + deltaint*cons1*exp(cons2*tmp**2)  
            if(sumint.gt.interval1) then
              interval1 = interval1 + interval0
              ip = int(sumint/interval0)
              lctabrgz(2,ip-1) = tmp + 0.5*deltaint
            else
            endif
          enddo
          lctabrgz(2,nprocrow-1) = zmax

          lctabrgz(1,0) = zmin
          do i = 1, nprocrow-1
            lctabrgz(1,i) = lctabrgz(2,i-1)
          enddo
        endif

        do j = 0, nproccol-1
          do i = 0, nprocrow-1
            if(flg.eq.2) then
              this%LcTabrg(1,i,j) = lctabrgz(1,i) 
              this%LcTabrg(2,i,j) = lctabrgz(2,i)
            else 
              this%LcTabrg(1,i,j) = zmin + i*deltaz
              this%LcTabrg(2,i,j) = zmin + (i+1)*deltaz
            endif
            ih=int((this%Lctabrg(2,i,j)-zmin)/hz) + 1
            il=int((this%Lctabrg(1,i,j)-zmin)/hz) + 1
            if(i.eq.0) then
              this%LcTabnm(1,i,j)=ih - il + 1
            else
              this%LcTabnm(1,i,j)=ih - il
            endif
          enddo
        enddo

        do j = 0, nproccol - 1
          nsum = 0
          do i = 0, nprocrow -2
            nsum = nsum + this%LcTabnm(1,i,j)
          enddo
          this%LcTabnm(1,nprocrow-1,j) = nz - nsum
        enddo

        endif

        if(nproccol.gt.1) then
          if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
            do i = 0, nprocrow-1
              do j = 0, nproccol-1
                ! uniform partition.
                this%LcTabrg(3,i,j) = ymin + j*deltay
                this%LcTabrg(4,i,j) = ymin + (j+1)*deltay
                ih=int((this%Lctabrg(4,i,j)-ymin)/hy) + 1
                il=int((this%Lctabrg(3,i,j)-ymin)/hy) + 1
                !print*,"ih: ",ih,il
                if(j.eq.0) then
                  this%LcTabnm(2,i,j)=ih - il + 1
                else
                  this%LcTabnm(2,i,j)=ih - il
                endif
              enddo
            enddo
          else
            if(flg.eq.2) then
              numint = 50001
              deltaint = (ymax-ymin)/(numint-1)
              cons1 = sqrt(1.0-uy*uy)/sig11y/sqrt(2.0*pi)
              cons2 = -0.5*(1.0-uy*uy)/sig11y/sig11y

              sumint = 0.0
              do i=1, numint-1
                tmp = ymin + (i-1)*deltaint
                sumint = sumint + deltaint*cons1*exp(cons2*tmp**2)
              enddo

              interval0 = sumint/float(nproccol)
              interval1 = interval0
              sumint = 0.0

              do i=1, numint-1
                tmp = ymin + (i-1)*deltaint
                sumint = sumint + deltaint*cons1*exp(cons2*tmp**2)
                if(sumint.gt.interval1) then
                  interval1 = interval1 + interval0
                  ip = int(sumint/interval0)
                  lctabrgy(2,ip-1) = tmp + 0.5*deltaint
                else
                endif
              enddo
              lctabrgy(2,nproccol-1) = ymax

              lctabrgy(1,0) = ymin
              do i = 1, nproccol-1
                lctabrgy(1,i) = lctabrgy(2,i-1)
              enddo
            endif

            do i = 0, nprocrow-1
              do j = 0, nproccol-1
                if(flg.eq.2) then
                  this%LcTabrg(3,i,j) = lctabrgy(1,j)
                  this%LcTabrg(4,i,j) = lctabrgy(2,j)
                else 
                  this%LcTabrg(3,i,j) = ymin + j*deltay
                  this%LcTabrg(4,i,j) = ymin + (j+1)*deltay
                endif
                ih=int((this%Lctabrg(4,i,j)-ymin)/hy) + 1
                il=int((this%Lctabrg(3,i,j)-ymin)/hy) + 1
                if(j.eq.0) then
                  this%LcTabnm(2,i,j)=ih - il + 1
                else
                  this%LcTabnm(2,i,j)=ih - il
                endif
              enddo
            enddo
          endif

          do i = 0, nprocrow - 1
            nsum = 0
            do j = 0, nproccol-2
              nsum = nsum + this%LcTabnm(2,i,j)
            enddo
            !print*,"nsum: ",nsum
            this%LcTabnm(2,i,nproccol-1) = ny - nsum
          enddo
        endif

        if(nprocrow.gt.1) then

        do i = 0, nprocrow - 1
          if(this%LcTabnm(1,i,0).eq.0) then
            this%LcTabnm(1,i,0) = 1
            this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) - hz
            iflag = 0
            do i0 = i-1,0,-1
              if(this%LcTabnm(1,i0,0).gt.1) then
                this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
                iflag = 1
                exit
              else
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)-hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
              endif
            enddo
            if(iflag.eq.0) then
              this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) + hz
              do i0 = 0, i-1
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
              enddo
              this%LcTabrg(2,i,0) = this%LcTabrg(2,i,0) + hz
              do i0 = i+1, nprocrow-1
                if(this%LcTabnm(1,i0,0).gt.1) then
                  this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  iflag = 1
                  exit
                else
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in z!"
              stop
            endif
          endif
        enddo

        do j = 1, nproccol -1
          do i = 0, nprocrow - 1
            this%LcTabnm(1,i,j) = this%LcTabnm(1,i,0)
            this%LcTabrg(1,i,j) = this%LcTabrg(1,i,0)
            this%LcTabrg(2,i,j) = this%LcTabrg(2,i,0)
          enddo
        enddo

        else
          do j = 0, nproccol-1
          this%LcTabrg(1,0,j) = zmin 
          this%LcTabrg(2,0,j) = zmin + deltaz
          this%LcTabnm(1,0,j) = nz 
          enddo
        endif

        if(nproccol.gt.1) then

        do j = 0, nproccol - 1
          if(this%LcTabnm(2,0,j).eq.0) then
            this%LcTabnm(2,0,j) = 1
            this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) - hy
            iflag = 0
            do j0 = j-1,0,-1
              if(this%LcTabnm(2,0,j0).gt.1) then
                this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
                iflag = 1
                exit
              else
                this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)-hy
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
              endif
            enddo
            if(iflag.eq.0) then
              this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) + hy
              do j0 = 0, j-1
                this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
              enddo
              this%LcTabrg(4,0,j) = this%LcTabrg(4,0,j) + hy
              do j0 = j+1, nproccol-1
                if(this%LcTabnm(2,0,j0).gt.1) then
                  this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                  this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                  iflag = 1
                  exit
                else
                  this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                  this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in y!"
              stop
            endif
          endif
        enddo

        do i = 1, nprocrow - 1
          do j = 0, nproccol - 1
            this%LcTabnm(2,i,j) = this%LcTabnm(2,0,j)
            this%LcTabrg(3,i,j) = this%LcTabrg(3,0,j)
            this%LcTabrg(4,i,j) = this%LcTabrg(4,0,j)
            !print*,"LcTabY: ",this%LcTabnm(2,i,j)
          enddo
        enddo

        else
          do i = 0, nprocrow-1
          this%LcTabrg(3,i,0) = ymin 
          this%LcTabrg(4,i,0) = ymin + deltay
          this%LcTabnm(2,i,0) = ny 
          enddo
        endif

        this%Sptrnglocal(5) = this%LcTabrg(1,myidx,myidy)
        this%Sptrnglocal(6) = this%LcTabrg(2,myidx,myidy)
        this%Sptrnglocal(3) = this%LcTabrg(3,myidx,myidy)
        this%Sptrnglocal(4) = this%LcTabrg(4,myidx,myidy)
        this%Sptrnglocal(1) = this%SpatRange(1)
        this%Sptrnglocal(2) = this%SpatRange(2)

        this%Mshlocal(3) = this%LcTabnm(1,myidx,myidy)
        this%Mshlocal(2) = this%LcTabnm(2,myidx,myidy)
        this%Mshlocal(1) = this%Meshnum(1)
        if(this%Mshlocal(3).lt.0 .or. this%Mshlocal(3).gt.nz) then
          print*,"wrong in the z direction local mesh number!", &
          this%Mshlocal(3),myidx,myidy,nz
        else if(this%Mshlocal(2).lt.0 .or. this%Mshlocal(2).gt.ny) then
          print*,"wrong in the y direction local mesh number!", &
          this%Mshlocal(2),myidx,myidy,ny
        else
        endif

!        call MPI_BARRIER(commrow,ierr)
        
        end subroutine init_CompDom

        !update geometry parameters using new particle positions.
        subroutine update_CompDom(this, ptrange, grid2d, Flagbc)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid2d
        integer, intent(in) :: Flagbc
        double precision, dimension(6), intent(in) :: ptrange
        double precision, dimension(3) :: localrange1, &
                              localrange2,temp1,temp2
        integer :: i,j,comm2d,commrow,commcol,nproc,nproccol,nprocrow,&
                   myid,myidx,myidy
        integer :: ierr,nsum
        double precision :: xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,&
                          epsx,epsy,epsz,scale,length
        integer :: nz,ny
        integer :: i0,j0,iflag
        double precision :: tp1,tp2,tp3,t0

        call starttime_Timer(t0)

        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,nproc,nproccol,nprocrow)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)

        do i = 1, 3
          localrange1(i) = ptrange(2*i-1)
          localrange2(i) = ptrange(2*i)
        enddo

        call MPI_ALLREDUCE(localrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,comm2d,ierr)
        call MPI_ALLREDUCE(localrange2,temp2,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,comm2d,ierr) 

        !print*,"temp: ",temp1(1),temp2(1),temp1(2),temp2(2),temp1(3),temp2(3)
        !epsx = 0.001
        !epsy = 0.001
        !epsz = 0.001
        epsx = 1.0/(this%Meshnum(1)-3.0) 
        epsy = 1.0/(this%Meshnum(2)-3.0) 
        epsz = 1.0/(this%Meshnum(3)-3.0) 

        if(Flagbc.eq.1) then
          xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
          xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
          ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
          ymax = temp2(2)+epsy*(temp2(2)-temp1(2))
          zmin = temp1(3)-epsz*(temp2(3)-temp1(3))
          zmax = temp2(3)+epsz*(temp2(3)-temp1(3))
        else if(Flagbc.eq.2) then
          xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
          xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
          ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
          ymax = temp2(2)+epsy*(temp2(2)-temp1(2))
          zmin = temp1(3)
          zmax = temp2(3)
        else if((Flagbc.eq.3).or.(Flagbc.eq.5)) then
          xmin = temp1(1)
          xmax = temp2(1)
          ymin = temp1(2)
          ymax = temp2(2)
          zmin = temp1(3)-epsz*(temp2(3)-temp1(3))
          zmax = temp2(3)+epsz*(temp2(3)-temp1(3))
        else if((Flagbc.eq.4).or.(Flagbc.eq.6)) then
          xmin = temp1(1)
          xmax = temp2(1)
          ymin = temp1(2)
          ymax = temp2(2)
          zmin = temp1(3)
          zmax = temp2(3)
        else
          print*,"no such boundary condition!!!"
          stop
        endif

        !print*,"xmin ",xmin,xmax,ymin,ymax,zmin,zmax

        scale = (zmax-zmin)/(this%SpatRange(6)-this%SpatRange(5))
        do j = 0, nproccol-1
          length = this%LcTabrg(2,0,j)-this%LcTabrg(1,0,j)
          this%LcTabrg(1,0,j) = zmin
          this%LcTabrg(2,0,j) = zmin + scale*length
          do i = 1, nprocrow-1
            length = this%LcTabrg(2,i,j)-this%LcTabrg(1,i,j)
            this%LcTabrg(1,i,j) = this%LcTabrg(2,i-1,j)
            this%LcTabrg(2,i,j) = this%LcTabrg(1,i,j) + &
                                      scale*length
          enddo
        enddo
            
        if(Flagbc.eq.3) then
        else

        scale = (ymax-ymin)/(this%SpatRange(4)-this%SpatRange(3))
        do i = 0, nprocrow-1
          length = this%LcTabrg(4,i,0)-this%LcTabrg(3,i,0)
          this%LcTabrg(3,i,0) = ymin
          this%LcTabrg(4,i,0) = ymin + scale*length
          do j = 1, nproccol-1
            length = this%LcTabrg(4,i,j)-this%LcTabrg(3,i,j)
            this%LcTabrg(3,i,j) = this%LcTabrg(4,i,j-1)
            this%LcTabrg(4,i,j) = this%LcTabrg(3,i,j) + &
                                      scale*length
          enddo
        enddo
        hx = (xmax-xmin)/(this%Meshnum(1)-1) 
        hy = (ymax-ymin)/(this%Meshnum(2)-1) 
        this%Meshsize(1) = hx
        this%Meshsize(2) = hy

        endif

        hz = (zmax-zmin)/(this%Meshnum(3)-1) 
        this%Meshsize(3) = hz

        if(nprocrow.gt.1) then

        do i = 0, nprocrow-1
          do j = 0, nproccol-1
            if(i.eq.0) then
              this%LcTabnm(1,i,j)=int((this%Lctabrg(2,i,j)-zmin)/hz)- &
                                  int((this%Lctabrg(1,i,j)-zmin)/hz)+1
            else
              this%LcTabnm(1,i,j)=int((this%Lctabrg(2,i,j)-zmin)/hz)- &
                                  int((this%Lctabrg(1,i,j)-zmin)/hz)
            endif
          enddo
        enddo

        do j = 0, nproccol - 1
          nsum = 0
          do i = 0, nprocrow -2
            nsum = nsum + this%LcTabnm(1,i,j)
          enddo
          this%LcTabnm(1,nprocrow-1,j) = this%Meshnum(3) - nsum
        enddo

        endif

        if(Flagbc.eq.3) then
        else
          if(nproccol.gt.1) then

            do i = 0, nprocrow-1
              do j = 0, nproccol-1
                if(j.eq.0) then
                  this%LcTabnm(2,i,j)=int((this%Lctabrg(4,i,j)-ymin)/hy)- &
                                      int((this%Lctabrg(3,i,j)-ymin)/hy)+1
                else
                  this%LcTabnm(2,i,j)=int((this%Lctabrg(4,i,j)-ymin)/hy)- &
                                      int((this%Lctabrg(3,i,j)-ymin)/hy)
                endif
              enddo
            enddo

            do i = 0, nprocrow - 1
              nsum = 0
              do j = 0, nproccol-2
                nsum = nsum + this%LcTabnm(2,i,j)
              enddo
              this%LcTabnm(2,i,nproccol-1) = this%Meshnum(2) - nsum
            enddo
          endif
        endif

        if(nprocrow.gt.1) then

        do i = 0, nprocrow - 1
          if(this%LcTabnm(1,i,0).eq.0) then
            this%LcTabnm(1,i,0) = 1
            this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) - hz
            iflag = 0
            do i0 = i-1,0,-1
              if(this%LcTabnm(1,i0,0).gt.1) then
                this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
                iflag = 1
                exit
              else
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)-hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
              endif
            enddo
            if(iflag.eq.0) then
              this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) + hz
              do i0 = 0, i-1
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
              enddo
              this%LcTabrg(2,i,0) = this%LcTabrg(2,i,0) + hz
              do i0 = i+1, nprocrow-1
                if(this%LcTabnm(1,i0,0).gt.1) then
                  this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  iflag = 1
                  exit
                else
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in z!"
              stop
            endif
          endif
        enddo

        do j = 1, nproccol -1
          do i = 0, nprocrow - 1
            this%LcTabnm(1,i,j) = this%LcTabnm(1,i,0)
            this%LcTabrg(1,i,j) = this%LcTabrg(1,i,0)
            this%LcTabrg(2,i,j) = this%LcTabrg(2,i,0)
          enddo
        enddo

        else
          do j = 0, nproccol-1
          this%LcTabnm(1,0,j) = this%Meshnum(3) 
          enddo
        endif
  
        if(Flagbc.eq.3) then
        else
          if(nproccol.gt.1) then
            do j = 0, nproccol - 1
              if(this%LcTabnm(2,0,j).eq.0) then
                this%LcTabnm(2,0,j) = 1
                this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) - hy
                iflag = 0
                do j0 = j-1,0,-1
                  if(this%LcTabnm(2,0,j0).gt.1) then
                    this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                    this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
                    iflag = 1
                    exit
                  else
                    this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)-hy
                    this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
                  endif
                enddo
                if(iflag.eq.0) then
                  this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) + hy
                  do j0 = 0, j-1
                    this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                    this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
                  enddo
                  this%LcTabrg(4,0,j) = this%LcTabrg(4,0,j) + hy
                  do j0 = j+1, nproccol-1
                    if(this%LcTabnm(2,0,j0).gt.1) then
                      this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                      this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                      iflag = 1
                      exit
                    else
                      this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                      this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
                    endif
                  enddo
                endif
                if(iflag.eq.0) then
                  print*,"Number of grids less than number of processors", &
                          "in y!"
                  stop
                endif
              endif
            enddo

            do i = 1, nprocrow - 1
              do j = 0, nproccol - 1
                this%LcTabnm(2,i,j) = this%LcTabnm(2,0,j)
                this%LcTabrg(3,i,j) = this%LcTabrg(3,0,j)
                this%LcTabrg(4,i,j) = this%LcTabrg(4,0,j)
              enddo
            enddo
          else
            do i = 0, nprocrow-1
              this%LcTabnm(2,i,0) = this%Meshnum(2) 
            enddo
          endif
          this%Sptrnglocal(3) = this%LcTabrg(3,myidx,myidy)
          this%Sptrnglocal(4) = this%LcTabrg(4,myidx,myidy)
          this%Sptrnglocal(1) = xmin
          this%Sptrnglocal(2) = xmax
          this%Mshlocal(3) = this%LcTabnm(1,myidx,myidy)
          this%Mshlocal(2) = this%LcTabnm(2,myidx,myidy)
          this%SpatRange(1) = xmin
          this%SpatRange(2) = xmax
          this%SpatRange(3) = ymin
          this%SpatRange(4) = ymax
        endif

        this%Sptrnglocal(5) = this%LcTabrg(1,myidx,myidy)
        this%Sptrnglocal(6) = this%LcTabrg(2,myidx,myidy)

        this%Mshlocal(1) = this%Meshnum(1)
        nz = this%Meshnum(3)
        ny = this%Meshnum(2)
        if(this%Mshlocal(3).lt.0 .or. this%Mshlocal(3).gt.nz) then
          print*,"wrong in the z direction local mesh number!",&
          this%Mshlocal(3),nz,myidx,myidy
        else if(this%Mshlocal(2).lt.0 .or. this%Mshlocal(2).gt.ny) then
          print*,"wrong in the y direction local mesh number!",&
          this%Mshlocal(2),ny,myidx,myidy
        else
        endif

        this%SpatRange(5) = zmin
        this%SpatRange(6) = zmax

        t_boundgeom = t_boundgeom + elapsedtime_Timer(t0)

        end subroutine update_CompDom

        !update geometry parameters using new particle positions.
        subroutine updateold_CompDom(this, inrange, grid2d, nplc)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid2d
        double precision, dimension(:,:), intent(in) :: inrange
        integer, intent(in) :: nplc
        double precision, dimension(3) :: localrange1, &
                              localrange2,temp1,temp2
        integer :: i,j,comm2d,commrow,commcol,nproc,nproccol,nprocrow,&
                   myid,myidx,myidy
        integer :: ierr,nsum
        double precision :: xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,&
                          epsx,epsy,epsz,scale,length
        integer :: nz,ny
        integer :: i0,j0,iflag
        double precision :: tp1,tp2,tp3,t0

        call starttime_Timer(t0)

!        do i = 1, Ndim
!          localrange1(i) = minval(inrange(2*i-1,1:nplc),1)
!          localrange2(i) = maxval(inrange(2*i-1,1:nplc),1)
!        enddo
        ! for cache optimization
        if(nplc.ne.0) then
          do i = 1, Ndim
            localrange1(i) = inrange(2*i-1,1)
            localrange2(i) = inrange(2*i-1,1)
          enddo 
        else
          do i = 1, Ndim
            localrange1(i) = 10.0
            localrange2(i) = -10.0
          enddo 
        endif
        do j = 1, nplc
          do i = 1, Ndim
            if(localrange1(i).gt.inrange(2*i-1,j)) then
              localrange1(i) = inrange(2*i-1,j)
            endif
            if(localrange2(i).lt.inrange(2*i-1,j)) then
              localrange2(i) = inrange(2*i-1,j)
            endif
          enddo
        enddo

        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,nproc,nproccol,nprocrow)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)

        call MPI_ALLREDUCE(localrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,comm2d,ierr)
        call MPI_ALLREDUCE(localrange2,temp2,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,comm2d,ierr) 

        !print*,"temp: ",temp1(1),temp2(1),temp1(2),temp2(2),temp1(3),temp2(3)
        epsx = 0.001
        epsy = 0.001
        epsz = 0.001
        !epsx = 0.5*((this%Meshnum(1)-1.0)/(this%Meshnum(1)-3.0)-1) &
        !       + 0.000001
        !epsy = 0.5*((this%Meshnum(2)-1.0)/(this%Meshnum(2)-3.0)-1) &
        !       + 0.000001
        !epsz = 0.5*((this%Meshnum(3)-1.0)/(this%Meshnum(3)-3.0)-1) &
        !       + 0.000001

        xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
        xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
        ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
        ymax = temp2(2)+epsy*(temp2(2)-temp1(2))
        zmin = temp1(3)-epsz*(temp2(3)-temp1(3))
        zmax = temp2(3)+epsz*(temp2(3)-temp1(3))
        !print*,"xmin ",xmin,xmax,ymin,ymax,zmin,zmax

        scale = (zmax-zmin)/(this%SpatRange(6)-this%SpatRange(5))
        do j = 0, nproccol-1
          length = this%LcTabrg(2,0,j)-this%LcTabrg(1,0,j)
          this%LcTabrg(1,0,j) = zmin
          this%LcTabrg(2,0,j) = zmin + scale*length
          do i = 1, nprocrow-1
            length = this%LcTabrg(2,i,j)-this%LcTabrg(1,i,j)
            this%LcTabrg(1,i,j) = this%LcTabrg(2,i-1,j)
            this%LcTabrg(2,i,j) = this%LcTabrg(1,i,j) + &
                                      scale*length
          enddo
        enddo
            
        scale = (ymax-ymin)/(this%SpatRange(4)-this%SpatRange(3))
        do i = 0, nprocrow-1
          length = this%LcTabrg(4,i,0)-this%LcTabrg(3,i,0)
          this%LcTabrg(3,i,0) = ymin
          this%LcTabrg(4,i,0) = ymin + scale*length
          do j = 1, nproccol-1
            length = this%LcTabrg(4,i,j)-this%LcTabrg(3,i,j)
            this%LcTabrg(3,i,j) = this%LcTabrg(4,i,j-1)
            this%LcTabrg(4,i,j) = this%LcTabrg(3,i,j) + &
                                      scale*length
          enddo
        enddo

        hx = (xmax-xmin)/(this%Meshnum(1)-1) 
        hy = (ymax-ymin)/(this%Meshnum(2)-1) 
        hz = (zmax-zmin)/(this%Meshnum(3)-1) 
        this%Meshsize(1) = hx
        this%Meshsize(2) = hy
        this%Meshsize(3) = hz

        if(nprocrow.gt.1) then

        do i = 0, nprocrow-1
          do j = 0, nproccol-1
            if(i.eq.0) then
              this%LcTabnm(1,i,j)=int((this%Lctabrg(2,i,j)-zmin)/hz)- &
                                  int((this%Lctabrg(1,i,j)-zmin)/hz)+1
            else
              this%LcTabnm(1,i,j)=int((this%Lctabrg(2,i,j)-zmin)/hz)- &
                                  int((this%Lctabrg(1,i,j)-zmin)/hz)
            endif
          enddo
        enddo

        do j = 0, nproccol - 1
          nsum = 0
          do i = 0, nprocrow -2
            nsum = nsum + this%LcTabnm(1,i,j)
          enddo
          this%LcTabnm(1,nprocrow-1,j) = this%Meshnum(3) - nsum
        enddo

        endif

        if(nproccol.gt.1) then

        do i = 0, nprocrow-1
          do j = 0, nproccol-1
            if(j.eq.0) then
              this%LcTabnm(2,i,j)=int((this%Lctabrg(4,i,j)-ymin)/hy)- &
                                  int((this%Lctabrg(3,i,j)-ymin)/hy)+1
            else
              this%LcTabnm(2,i,j)=int((this%Lctabrg(4,i,j)-ymin)/hy)- &
                                  int((this%Lctabrg(3,i,j)-ymin)/hy)
            endif
          enddo
        enddo

        do i = 0, nprocrow - 1
          nsum = 0
          do j = 0, nproccol-2
            nsum = nsum + this%LcTabnm(2,i,j)
          enddo
          this%LcTabnm(2,i,nproccol-1) = this%Meshnum(2) - nsum
        enddo

        endif

        if(nprocrow.gt.1) then

        do i = 0, nprocrow - 1
          if(this%LcTabnm(1,i,0).eq.0) then
            this%LcTabnm(1,i,0) = 1
            this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) - hz
            iflag = 0
            do i0 = i-1,0,-1
              if(this%LcTabnm(1,i0,0).gt.1) then
                this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
                iflag = 1
                exit
              else
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)-hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)-hz
              endif
            enddo
            if(iflag.eq.0) then
              this%LcTabrg(1,i,0) = this%LcTabrg(1,i,0) + hz
              do i0 = 0, i-1
                this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
              enddo
              this%LcTabrg(2,i,0) = this%LcTabrg(2,i,0) + hz
              do i0 = i+1, nprocrow-1
                if(this%LcTabnm(1,i0,0).gt.1) then
                  this%LcTabnm(1,i0,0) = this%LcTabnm(1,i0,0)-1
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  iflag = 1
                  exit
                else
                  this%LcTabrg(1,i0,0) = this%LcTabrg(1,i0,0)+hz
                  this%LcTabrg(2,i0,0) = this%LcTabrg(2,i0,0)+hz
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in z!"
              stop
            endif
          endif
        enddo

        do j = 1, nproccol -1
          do i = 0, nprocrow - 1
            this%LcTabnm(1,i,j) = this%LcTabnm(1,i,0)
            this%LcTabrg(1,i,j) = this%LcTabrg(1,i,0)
            this%LcTabrg(2,i,j) = this%LcTabrg(2,i,0)
          enddo
        enddo

        else
          do j = 0, nproccol-1
          this%LcTabnm(1,0,j) = this%Meshnum(3) 
          enddo
        endif
  
        if(nproccol.gt.1) then

        do j = 0, nproccol - 1
          if(this%LcTabnm(2,0,j).eq.0) then
            this%LcTabnm(2,0,j) = 1
            this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) - hy
            iflag = 0
            do j0 = j-1,0,-1
              if(this%LcTabnm(2,0,j0).gt.1) then
                this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
                iflag = 1
                exit
              else
                this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)-hy
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)-hy
              endif
            enddo
            if(iflag.eq.0) then
              this%LcTabrg(3,0,j) = this%LcTabrg(3,0,j) + hy
              do j0 = 0, j-1
                this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
              enddo
              this%LcTabrg(4,0,j) = this%LcTabrg(4,0,j) + hy
              do j0 = j+1, nproccol-1
                if(this%LcTabnm(2,0,j0).gt.1) then
                  this%LcTabnm(2,0,j0) = this%LcTabnm(2,0,j0)-1
                  this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                  iflag = 1
                  exit
                else
                  this%LcTabrg(3,0,j0) = this%LcTabrg(3,0,j0)+hy
                  this%LcTabrg(4,0,j0) = this%LcTabrg(4,0,j0)+hy
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in y!"
              stop
            endif
          endif
        enddo

        do i = 1, nprocrow - 1
          do j = 0, nproccol - 1
            this%LcTabnm(2,i,j) = this%LcTabnm(2,0,j)
            this%LcTabrg(3,i,j) = this%LcTabrg(3,0,j)
            this%LcTabrg(4,i,j) = this%LcTabrg(4,0,j)
          enddo
        enddo

        else
          do i = 0, nprocrow-1
            this%LcTabnm(2,i,0) = this%Meshnum(2) 
          enddo
        endif

        this%Sptrnglocal(5) = this%LcTabrg(1,myidx,myidy)
        this%Sptrnglocal(6) = this%LcTabrg(2,myidx,myidy)
        this%Sptrnglocal(3) = this%LcTabrg(3,myidx,myidy)
        this%Sptrnglocal(4) = this%LcTabrg(4,myidx,myidy)
        this%Sptrnglocal(1) = xmin
        this%Sptrnglocal(2) = xmax

        this%Mshlocal(3) = this%LcTabnm(1,myidx,myidy)
        this%Mshlocal(2) = this%LcTabnm(2,myidx,myidy)
        this%Mshlocal(1) = this%Meshnum(1)
        nz = this%Meshnum(3)
        ny = this%Meshnum(2)
        if(this%Mshlocal(3).lt.0 .or. this%Mshlocal(3).gt.nz) then
          print*,"wrong in the z direction local mesh number!",&
          this%Mshlocal(3),nz,myidx,myidy
        else if(this%Mshlocal(2).lt.0 .or. this%Mshlocal(2).gt.ny) then
          print*,"wrong in the y direction local mesh number!",&
          this%Mshlocal(2),ny,myidx,myidy
        else
        endif

        this%SpatRange(1) = xmin
        this%SpatRange(2) = xmax
        this%SpatRange(3) = ymin
        this%SpatRange(4) = ymax
        this%SpatRange(5) = zmin
        this%SpatRange(6) = zmax

        t_boundgeom = t_boundgeom + elapsedtime_Timer(t0)

        end subroutine updateold_CompDom

        ! find the balanced local domain geometry so that the number 
        ! of particles on this domain about equal.
        subroutine balance_CompDom(source,lctabnmz,lctabnmy,&
        lctabrgz,lctabrgy,npz,npy,commrow,commcol,innx,inny,innz,inyglb,&
        inzglb,hy,hz,ymin,zmin)
        implicit none
        include 'mpif.h'
        integer, parameter :: Ndiv = 10
        double precision, dimension(:,:,:) :: source
        double precision, intent(in) :: hy,hz,ymin,zmin
        integer, intent(in) :: npz,npy,commrow,commcol
        integer, dimension(0:npz-1), intent(inout) :: lctabnmz
        integer, dimension(0:npy-1), intent(inout) :: lctabnmy
        double precision, dimension(2,0:npz-1), intent(out) :: lctabrgz
        double precision, dimension(2,0:npy-1), intent(out) :: lctabrgy
        integer, dimension(0:npz-1) :: ztable,zdisp
        integer, dimension(0:npy-1) :: ytable,ydisp
        double precision:: temp
        double precision, dimension(innz) :: rhozlocal, rhoz
        double precision, dimension(inzglb) :: recvrhoz
        double precision, dimension(Ndiv*inzglb) :: rhoznew
        double precision, dimension(inny) :: rhoylocal, rhoy
        double precision, dimension(inyglb) :: recvrhoy
        double precision, dimension(Ndiv*inyglb) :: rhoynew
        integer :: i,j,k,inyglb,inzglb,innx,inny,innz,ierr
        integer :: icount,ip,sumcount
        double precision :: interval0,interval1,sumint
        integer, dimension(npz) :: nozeronmz
        integer, dimension(npy) :: nozeronmy
        integer :: i0, nozero,iflag,jadd,kadd
        double precision :: upper,hnew1,hnew2

        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npz.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        do i = 1, innz
          rhozlocal(i) = 0.0
        enddo
        do i = 1, inny
          rhoylocal(i) = 0.0
        enddo
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
               rhoylocal(j) = rhoylocal(j) + source(i,j+jadd,k+kadd)
               rhozlocal(k) = rhozlocal(k) + source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        !start the z direction.
        if(npz.gt.1) then

        call MPI_ALLREDUCE(rhozlocal,rhoz,innz,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,commcol,ierr)

        do i = 0, npz-1
          ztable(i) = lctabnmz(i)
        enddo

        zdisp(0) = 0
        do i = 1, npz-1
          zdisp(i) = zdisp(i-1)+ztable(i-1)
        enddo

        call MPI_ALLGATHERV(rhoz,innz,MPI_DOUBLE_PRECISION,recvrhoz,&
                            ztable,&
                            zdisp,MPI_DOUBLE_PRECISION,commrow,ierr)

        do i = 1, Ndiv*inzglb
          rhoznew(i) = recvrhoz(int((i-1)/Ndiv)+1)/Ndiv
        enddo
        hnew1 = 0.5*hz/(Ndiv-1)
        hnew2 = hz/(Ndiv)

        sumint = 0.0
        temp = sum(recvrhoz)
        interval0 = temp/dble(npz)
        interval1 = interval0 
        icount = 0

!        temp = sum(rhoznew)
!        print*, "temp: ",temp

        lctabnmz = 0

        do i = 1, Ndiv*inzglb
          sumint = sumint + rhoznew(i)
          ip = int(sumint/interval0)
          if(ip.ne.npz) then
            if(i.le.Ndiv) then
              lctabrgz(2,ip) = zmin + (i-1)*hnew1 + 0.5*hnew1
            else if((i.gt.Ndiv).and.(i.le.Ndiv*inzglb-Ndiv)) then
              lctabrgz(2,ip) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
            else
              lctabrgz(2,ip) = zmin + (inzglb-1)*hz - 0.5*hz + &
                               (i+Ndiv-Ndiv*inzglb)*hz/dble(Ndiv)
            endif
          endif
        enddo
        lctabrgz(2,npz-1) = zmin + (inzglb-1)*hz
        lctabrgz(1,0) = zmin
        do i = 1,npz-1
          lctabrgz(1,i) = lctabrgz(2,i-1)
        enddo

        do i = 0, npz-1
          if(i.eq.0) then
            lctabnmz(i)=int((lctabrgz(2,i)-zmin)/hz)- &
                                  int((lctabrgz(1,i)-zmin)/hz)+1
          else
            lctabnmz(i)=int((lctabrgz(2,i)-zmin)/hz)- &
                                  int((lctabrgz(1,i)-zmin)/hz)
          endif
        enddo

        nozero = 0
        do i = 0, npz-1
          if(lctabnmz(i).gt.0) then
            nozero = nozero + 1
            nozeronmz(nozero) = i
          endif
        enddo

        if(nozero.lt.npz) then
        do i = 1, nozero-1
          if((nozeronmz(i+1)-nozeronmz(i)).gt.1) then
            upper=zmin+(int((lctabrgz(2,nozeronmz(i))-zmin)/hz)+1)*hz
            do i0 = 1, nozeronmz(i+1)-nozeronmz(i)-1
              lctabrgz(2,nozeronmz(i)+i0) = lctabrgz(2,nozeronmz(i))+ &
              i0*(upper-lctabrgz(2,nozeronmz(i)))/ &
              dble((nozeronmz(i+1)-nozeronmz(i)))
            enddo
          endif
        enddo
        endif

        sumcount = 0
        do i = 0, npz-2
          sumcount = sumcount + lctabnmz(i)
        enddo

        lctabnmz(npz-1) = inzglb - sumcount 
        lctabrgz(2,npz-1) = zmin + (inzglb-1)*hz

        lctabrgz(1,0) = zmin
        do i = 1, npz-1
          lctabrgz(1,i) = lctabrgz(2,i-1)
        enddo

        do i = 0, npz - 1
          if(lctabnmz(i).eq.0) then
            lctabnmz(i) = 1
            lctabrgz(1,i) = lctabrgz(1,i) - hz
            iflag = 0
            do i0 = i-1,0,-1
              if(lctabnmz(i0).gt.1) then
                lctabnmz(i0) = lctabnmz(i0) - 1
                lctabrgz(2,i0) = lctabrgz(2,i0) - hz
                iflag = 1
                exit
              else
                lctabrgz(1,i0) = lctabrgz(1,i0) - hz
                lctabrgz(2,i0) = lctabrgz(2,i0) - hz
              endif
            enddo
            if(iflag.eq.0) then
              lctabrgz(1,i) = lctabrgz(1,i) + hz
              do i0 = 0, i-1
                lctabrgz(1,i0) = lctabrgz(1,i0) + hz
                lctabrgz(2,i0) = lctabrgz(2,i0) + hz
              enddo
              lctabrgz(2,i) = lctabrgz(2,i) + hz
              do i0 = i+1, npz - 1
                if(lctabnmz(i0).gt.1) then
                  lctabnmz(i0) = lctabnmz(i0) - 1
                  lctabrgz(1,i0) = lctabrgz(1,i0) + hz
                  iflag = 1
                  exit
                else
                  lctabrgz(1,i0) = lctabrgz(1,i0) + hz
                  lctabrgz(2,i0) = lctabrgz(2,i0) + hz
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in z!"
              stop
            endif
          endif
        enddo

        endif

        !start the y direction.
!        print*, "start y"
        
        if(npy.gt.1) then

        call MPI_ALLREDUCE(rhoylocal,rhoy,inny,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,&
                           commrow,ierr)

        do i = 0, npy-1
          ytable(i) = lctabnmy(i)
        enddo

        ydisp(0) = 0
        do i = 1, npy-1
          ydisp(i) = ydisp(i-1)+ytable(i-1)
        enddo

        call MPI_ALLGATHERV(rhoy,inny,MPI_DOUBLE_PRECISION,recvrhoy,&
                            ytable,&
                            ydisp,MPI_DOUBLE_PRECISION,commcol,ierr)

        do i = 1, Ndiv*inyglb
          rhoynew(i) = recvrhoy(int((i-1)/Ndiv)+1)/Ndiv
        enddo
        hnew1 = 0.5*hy/(Ndiv-1)
        hnew2 = hy/(Ndiv)

        sumint = 0.0
        temp = sum(recvrhoy)
        interval0 = temp/dble(npy)
        interval1 = interval0 
        icount = 0

!        temp = sum(recvrhoy)

        lctabnmy = 0
        do i = 1, Ndiv*inyglb
          sumint = sumint + rhoynew(i)
          ip = int(sumint/interval0)
          if(ip.ne.npy) then
            if(i.le.Ndiv) then
              lctabrgy(2,ip) = ymin + (i-1)*hnew1 + 0.5*hnew1
            else if((i.gt.Ndiv).and.(i.le.Ndiv*inyglb-Ndiv)) then
              lctabrgy(2,ip) = ymin +0.5*hy+(i-Ndiv)*hnew2 + 0.5*hnew2
            else
              lctabrgy(2,ip) = ymin + (inyglb-1)*hy - 0.5*hy + &
                               (i+Ndiv-Ndiv*inyglb)*hy/dble(Ndiv)
            endif
          endif
        enddo
        lctabrgy(2,npy-1) = ymin + (inyglb-1)*hy
        lctabrgy(1,0) = ymin
        do i = 1,npy-1
          lctabrgy(1,i) = lctabrgy(2,i-1)
        enddo

        do i = 0, npy-1
          if(i.eq.0) then
            lctabnmy(i)=int((lctabrgy(2,i)-ymin)/hy)- &
                                  int((lctabrgy(1,i)-ymin)/hy)+1
          else
            lctabnmy(i)=int((lctabrgy(2,i)-ymin)/hy)- &
                                  int((lctabrgy(1,i)-ymin)/hy)
          endif
        enddo

        nozero = 0
        do i = 0, npy-1
          if(lctabnmy(i).gt.0) then
            nozero = nozero + 1
            nozeronmy(nozero) = i
          endif
        enddo

        if(nozero.lt.npy) then
        do i = 1, nozero-1
          if((nozeronmy(i+1)-nozeronmy(i)).gt.1) then
            upper=ymin+(int((lctabrgy(2,nozeronmy(i))-ymin)/hy)+1)*hy
            do i0 = 1, nozeronmy(i+1)-nozeronmy(i)-1
              lctabrgy(2,nozeronmy(i)+i0) = lctabrgy(2,nozeronmy(i))+ &
              i0*(upper-lctabrgy(2,nozeronmy(i)))/ &
              dble((nozeronmy(i+1)-nozeronmy(i)))
            enddo
          endif
        enddo
        endif

        sumcount = 0
        do i = 0, npy-2
          sumcount = sumcount + lctabnmy(i)
        enddo

        lctabnmy(npy-1) = inyglb - sumcount
        lctabrgy(2,npy-1) = ymin + (inyglb-1)*hy

        lctabrgy(1,0) = ymin
        do i = 1, npy-1
          lctabrgy(1,i) = lctabrgy(2,i-1)
        enddo

        do i = 0, npy - 1
          if(lctabnmy(i).eq.0) then
            lctabnmy(i) = 1
            lctabrgy(1,i) = lctabrgy(1,i) - hy
            iflag = 0
            do i0 = i-1,0,-1
              if(lctabnmy(i0).gt.1) then
                lctabnmy(i0) = lctabnmy(i0) - 1
                lctabrgy(2,i0) = lctabrgy(2,i0) - hy
                iflag = 1
                exit
              else
                lctabrgy(1,i0) = lctabrgy(1,i0) - hy
                lctabrgy(2,i0) = lctabrgy(2,i0) - hy
              endif
            enddo
            if(iflag.eq.0) then
              lctabrgy(1,i) = lctabrgy(1,i) + hy
              do i0 = 0, i-1
                lctabrgy(1,i0) = lctabrgy(1,i0) + hy
                lctabrgy(2,i0) = lctabrgy(2,i0) + hy
              enddo
              lctabrgy(2,i) = lctabrgy(2,i) + hy
              do i0 = i+1, npy - 1
                if(lctabnmy(i0).gt.1) then
                  lctabnmy(i0) = lctabnmy(i0) - 1
                  lctabrgy(1,i0) = lctabrgy(1,i0) + hy
                  iflag = 1
                  exit
                else
                  lctabrgy(1,i0) = lctabrgy(1,i0) + hy
                  lctabrgy(2,i0) = lctabrgy(2,i0) + hy
                endif
              enddo
            endif
            if(iflag.eq.0) then
              print*,"Number of grids less than number of processors", &
                      "in y!"
              stop
            endif
          endif
        enddo

        endif

        if(npz.le.1) then
          lctabnmz(0) = innz
          lctabrgz(1,0) = zmin
          lctabrgz(2,0) = zmin + (innz-1)*hz
        endif

        if(npy.le.1) then
          lctabnmy(0) = inny
          lctabrgy(1,0) = ymin
          lctabrgy(2,0) = ymin + (inny-1)*hy
        endif

        if(minval(lctabnmz).lt.0 .or. maxval(lctabnmz).gt.inzglb) then
          print*,"wrong in the z number of grid in load balance!",&
          minval(lctabnmz),maxval(lctabnmz)
        else if(minval(lctabnmy).lt.0 .or. maxval(lctabnmy).gt.inyglb) &
        then
          print*,"wrong in the y number of grid in load balance!",&
          minval(lctabnmy),maxval(lctabnmy)
        else
        endif

        end subroutine balance_CompDom

        subroutine getmsize_CompDom(this,msize)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        double precision, dimension(:), intent(out) :: msize

        msize = this%Meshsize

        end subroutine getmsize_CompDom

        subroutine getmnum_CompDom(this,mnum)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        integer, dimension(:), intent(out) :: mnum

        mnum = this%Meshnum

        end subroutine getmnum_CompDom

        subroutine getlcmnum_CompDom(this,mnum)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        integer, dimension(:), intent(out) :: mnum

        mnum = this%Mshlocal

        end subroutine getlcmnum_CompDom

        subroutine getrange_CompDom(this,range)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        double precision, dimension(:), intent(out) :: range

        range = this%SpatRange

        end subroutine getrange_CompDom
        
        subroutine getlcrange_CompDom(this,range)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        double precision, dimension(:), intent(out) :: range

        range = this%Sptrnglocal

        end subroutine getlcrange_CompDom

        subroutine getlctabrg_CompDom(this,lctable)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        double precision, dimension(:,:,:), intent(out) :: lctable

        lctable = this%LcTabrg

        end subroutine getlctabrg_CompDom

        subroutine getlctabnm_CompDom(this,lctable)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(in) :: this
        integer, dimension(:,:,:), intent(out) :: lctable

        lctable = this%LcTabnm

        end subroutine getlctabnm_CompDom

        subroutine setlctab1_CompDom(this,lctabnmz,lctabnmy,&
                      lctabrgz,lctabrgy,npx,npy,myidx,myidy)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(inout) :: this
        integer, intent(in) :: npx,npy,myidx,myidy
        integer, dimension(0:npx-1), intent(in) :: lctabnmz
        integer, dimension(0:npy-1), intent(in) :: lctabnmy
        double precision, dimension(2,0:npx-1), intent(in) :: lctabrgz
        double precision, dimension(2,0:npy-1), intent(in) :: lctabrgy
        integer :: i,j

        do j = 0, npy-1
          do i = 0, npx-1
            this%LcTabrg(1,i,j) = lctabrgz(1,i)
            this%LcTabrg(2,i,j) = lctabrgz(2,i)
            this%LcTabrg(3,i,j) = lctabrgy(1,j)
            this%LcTabrg(4,i,j) = lctabrgy(2,j)
            this%LcTabnm(1,i,j) = lctabnmz(i)
            this%LcTabnm(2,i,j) = lctabnmy(j)
          enddo
        enddo
        
        this%Sptrnglocal(5) = this%LcTabrg(1,myidx,myidy)
        this%Sptrnglocal(6) = this%LcTabrg(2,myidx,myidy)
        this%Sptrnglocal(3) = this%LcTabrg(3,myidx,myidy)
        this%Sptrnglocal(4) = this%LcTabrg(4,myidx,myidy)

        this%Mshlocal(3) = this%LcTabnm(1,myidx,myidy)
        this%Mshlocal(2) = this%LcTabnm(2,myidx,myidy)

        end subroutine setlctab1_CompDom

        subroutine setlctab2_CompDom(this,lctabnmz,lctabrgz,npx,npy,&
                                     myidx,myidy)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(inout) :: this
        integer, intent(in) :: npx,npy,myidx,myidy
        integer, dimension(0:npx-1), intent(in) :: lctabnmz
        double precision, dimension(2,0:npx-1), intent(in) :: lctabrgz
        integer :: i,j

        do j = 0, npy-1
          do i = 0, npx-1
            this%LcTabrg(1,i,j) = lctabrgz(1,i)
            this%LcTabrg(2,i,j) = lctabrgz(2,i)
            this%LcTabnm(1,i,j) = lctabnmz(i)
          enddo
        enddo
        
        this%Sptrnglocal(5) = this%LcTabrg(1,myidx,myidy)
        this%Sptrnglocal(6) = this%LcTabrg(2,myidx,myidy)

        this%Mshlocal(3) = this%LcTabnm(1,myidx,myidy)

        end subroutine setlctab2_CompDom

        subroutine destruct_CompDom(this)
        implicit none
        include 'mpif.h'
        type (CompDom), intent(out) :: this

        deallocate(this%LcTabrg)
        deallocate(this%LcTabnm)

        end subroutine destruct_CompDom

      end module CompDomclass
