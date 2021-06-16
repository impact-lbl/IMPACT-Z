!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Dataclass: Field data class in DATA STRUCTURE layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class stores the rf cavity data Ez, Ez', Ez'' on the 
!              axis; Fourier coefficients of Ez on the axis; Ez(r,z),
!              Er(r,z), Htheta(r,z) on the r-z grid plane; and Ex(x,y,z),
!              Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z) on
!              uniform x, y, z grid.
! Comments: 
!----------------------------------------------------------------
      module Dataclass
#if USE_MPI != 1
        use mpistub
#endif
        save
!-----------------------------------------------------------------------
! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) directly.
        !number of grid points along x, y, and z direction.
        integer :: NxIntvRfg = 1
        integer :: NyIntvRfg = 1
        integer :: NzIntvRfg = 1
        !range in x, y, and zdirections.
        double precision :: XmaxRfg,XminRfg,YmaxRfg,YminRfg,ZmaxRfg,ZminRfg
        ! discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z) and Bx(x,y,z), By(x,y,z), and 
        ! Bz(x,y,z) rf data. Here, the grid is uniform in x, y and z.
        double precision,allocatable,dimension(:,:,:) :: &
               Exgrid,Eygrid,Ezgrid,Bxgrid,Bygrid,Bzgrid
!-----------------------------------------------------------------------
! using the r-z field data (Er,Ez,Htheta) directly.
        !number of grid points along r direction.
        integer :: NrIntvRf = 1
        !number of grid points along z direction.
        integer :: NzIntvRf = 1
        !range in r and z directions.
        double precision :: RmaxRf,RminRf,ZmaxRf,ZminRf
        ! discrete Ez(r,z), Er(r,z) and Htheta(r,z) rf data. Here, the grid
        ! is uniform in r and z.
        double precision,allocatable,dimension(:,:) :: &
               ezdata,erdata,btdata
!-----------------------------------------------------------------------
! using only on-axis field data and its derivities.
        !initial number of grid points on the axis.
        integer, parameter :: Ndataini = 5000
        ! discrete Ez(0,z), Ez'(0,z), Ez''(0,z) rf data.
        double precision,dimension(Ndataini) :: zdat,edat,epdat,eppdat
        ! discrete wake function, z, longitudinal, x and y wake function
        double precision,dimension(Ndataini) :: zdatwk,edatwk,epdatwk,eppdatwk
!----------------------------------------------------------------------
! using the Fourier coefficients
        !number of Fourier expansion coefficients.
        integer, parameter :: NcoefF = 401
        double precision,dimension(NcoefF) :: Fcoef
        !Fcoef(1): constant
        !Fcoef(2k): cosine term
        !Fcoef(2k+1): sine term
!----------------------------------------------------------------------
        ! practical number of grid data on the axis or Fourier coefficients.
        integer :: Ndata,Ndatawk
      contains
        !Initialize the data storage arrays.
        subroutine init_Data()
        implicit none
        include 'mpif.h' 
        integer :: i,j

        NzIntvRf = 1
        NrIntvRf = 1
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        do j = 1, NrIntvRf+1
          do i = 1, NzIntvRf+1
            ezdata(i,j) = 0.0
            erdata(i,j) = 0.0
            btdata(i,j) = 0.0
          enddo
        enddo

        do i = 1, Ndataini
          zdat(i) = 0.0
          edat(i) = 0.0
          epdat(i) = 0.0
          eppdat(i) = 0.0
        enddo

        Ndata = 1
        RminRf = 0.0
        RmaxRf = 1.0
        ZminRf = 0.0
        ZmaxRf = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        NxIntvRfg = 1
        NyIntvRfg = 1
        NzIntvRfg = 1
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0
        XminRfg = 0.0
        XmaxRfg = 1.0
        YminRfg = 0.0
        YmaxRfg = 1.0
        ZminRfg = 0.0
        ZmaxRfg = 1.0

        end subroutine init_Data

        subroutine destruct_Data()
        implicit none
        include 'mpif.h' 

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        deallocate(Exgrid)
        deallocate(Eygrid)
        deallocate(Ezgrid)
        deallocate(Bxgrid)
        deallocate(Bygrid)
        deallocate(Bzgrid)
         
        end subroutine destruct_Data

        !read in the discrete field data Ez(0,z), Ez'(0,z), Ez''(0,z) 
        !distribution along axis zdat from files "rfdatax or rfdataxx 
        !or rfdataxxx".
        subroutine read1_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1
        character*10 name1
        character*11 name2
        character*12 name3

        name1 = 'rfdatax.in'
        name2 = 'rfdataxx.in'
        name3 = 'rfdataxxx.in'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
!            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
            n = n + 1
            zdat(n) = tmp1
            edat(n) = tmp2
            epdat(n) = tmp3
            eppdat(n) = tmp4
            !write(15,100)zdat(n),edat(n),epdat(n),eppdat(n)
            !divided by 100 is due to the unit in rf data is cm.
          goto 50
77        continue
          close(14)
!          close(15)
          Ndata = n
          zdat1 = zdat(1)
          do i = 1, Ndata
            zdat(i) = zdat(i) - zdat1
          enddo
        endif
100     format(6(1x,1pe15.8))

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(epdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(eppdat,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read1_Data

        ! read in discrete Ez(r,z), Er(r,z) and Btheta(r,z) rf data from
        ! files "rfdatax.in or rfdataxx.in or rfdataxxx.in. Here, the grid
        ! is uniform in r and z.
        subroutine read2_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0
        character*10 name1
        character*11 name2
        character*12 name3

        name1 = 'rfdatax.in'
        name2 = 'rfdataxx.in'
        name3 = 'rfdataxxx.in'
   
        mu0 = 4*2*asin(1.0d0)*1.0d-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are cm
          read(14,*,end=33)tmp1,tmp2,tmpint
          ZminRf = tmp1/100.0
          ZmaxRf = tmp2/100.0
          NzIntvRf = tmpint
          if(tmpint.ne.NzIntvRf) then
            print*,"input data wrong in Z: ",NzIntvRf,tmpint
            stop
          endif
          ! the input range units are cm
          read(14,*,end=33)tmp1
          read(14,*,end=33)tmp1,tmp2,tmpint
          RminRf = tmp1/100.0
          RmaxRf = tmp2/100.0
          NrIntvRf = tmpint
          if(tmpint.ne.NrIntvRf) then
            print*,"input data wrong in R: ",NrIntvRf,tmpint
            stop
          endif
          !print*,"Nz: ",NzIntvRf,ZminRf,ZmaxRf
          !print*,"Nr: ",NrIntvRf,RminRf,RmaxRf
        endif
33      continue
        call MPI_BCAST(NrIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        if(myrank.eq.0) then
          n = 0
50        continue
            if(mod(n,2).eq.0) then
              read(14,*,end=77)tmp1,tmp2,tmp3
              nn = n/2+1
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              ezdata(i,j) = tmp1
              erdata(i,j) = tmp2
              n = n + 1
              write(15,100)float(i-1),ezdata(i,j)
            else
              read(14,*,end=77)tmp1
              nn = (n+1)/2
              j  = (nn-1)/(NzIntvRf+1) + 1
              i = mod((nn-1),NzIntvRf+1) + 1
              btdata(i,j) = tmp1
              n = n + 1
            endif
          goto 50
77        continue
          close(14)
          close(15)
          Ndatalc = n/2
          !print*,"Ndata in 0: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RminRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(RmaxRf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ezdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(erdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(btdata,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndatalc
        end subroutine read2_Data

        !readin the Fourier coefficients of the RF field from files 
        !"rfdatax or rfdataxx or rfdataxxx".
        subroutine read3_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1
        character*10 name1
        character*11 name2
        character*12 name3

        name1 = 'rfdatax.in'
        name2 = 'rfdataxx.in'
        name3 = 'rfdataxxx.in'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1
            n = n + 1
            Fcoef(n) = tmp1
          goto 50
77        continue
          close(14)
          Ndata = n
        endif

        call MPI_BCAST(Ndata,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Fcoef,Ndata,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        !enforce the total # of coefficients to be an odd number.
        if(mod(Ndata,2).eq.0) then
          Ndata = Ndata + 1
          Fcoef(Ndata) = 0.0
        endif

        !print*,"Ndata: ",Ndata
        end subroutine read3_Data

        ! read in discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z),
        ! Bz(x,y,z) rf data from
        ! files "rfdatax.in or rfdataxx.in or rfdataxxx.in. Here, the grid
        ! is uniform in x, y and z.
        subroutine read4_Data(ifile)
        implicit none
        include 'mpif.h' 
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n,nn,Ndatalc,j,tmpint,k
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1,mu0,tmp5,tmp6
        character*10 name1
        character*11 name2
        character*12 name3

        name1 = 'rfdatax.in'
        name2 = 'rfdataxx.in'
        name3 = 'rfdataxxx.in'
   
        mu0 = 4*2*asin(1.0)*1.0e-7

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr) 

        !print*,"into read: "
        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
!            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          ! the input range units are m
          read(14,*,end=33)tmp1,tmp2,tmpint
          XminRfg = tmp1
          XmaxRfg = tmp2
          NxIntvRfg = tmpint
          read(14,*,end=33)tmp1,tmp2,tmpint
          YminRfg = tmp1
          YmaxRfg = tmp2
          NyIntvRfg = tmpint
          read(14,*,end=33)tmp1,tmp2,tmpint
          ZminRfg = tmp1
          ZmaxRfg = tmp2
          NzIntvRfg = tmpint
          !print*,"Nz: ",NzIntvRf,ZminRf,ZmaxRf
          !print*,"Nr: ",NrIntvRf,RminRf,RmaxRf
        endif
33      continue
        call MPI_BCAST(NxIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NyIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(NzIntvRfg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(YmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZminRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ZmaxRfg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        deallocate(Exgrid)
        deallocate(Eygrid)
        deallocate(Ezgrid)
        deallocate(Bxgrid)
        deallocate(Bygrid)
        deallocate(Bzgrid)
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0

        if(myrank.eq.0) then
          n = 0
50        continue
              read(14,*,end=77)tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
              n = n+1
              k = (n-1)/((NxIntvRfg+1)*(NyIntvRfg+1))+1
              j = (n-1-(k-1)*(NxIntvRfg+1)*(NyIntvRfg+1))/(NxIntvRfg+1) + 1
              i = n - (k-1)*(NxIntvRfg+1)*(NyIntvRfg+1) - (j-1)*(NxIntvRfg+1)
              Exgrid(i,j,k) = tmp1
              Eygrid(i,j,k) = tmp2
              Ezgrid(i,j,k) = tmp3
              Bxgrid(i,j,k) = tmp4
              Bygrid(i,j,k) = tmp5
              Bzgrid(i,j,k) = tmp6
          goto 50
77        continue
          close(14)
          Ndatalc = n
          print*,"Ndata in 0 PE: ",Ndatalc
        endif
100     format(2(1x,1pe15.8))

        call MPI_BCAST(Ndatalc,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Exgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Eygrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Ezgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bxgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bygrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Bzgrid,Ndatalc,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

!        if(myrank.eq.0) then
!          do i = 1,NzIntvRfg+1
!            write(14,1000)ZminRfg+(i-1)*(ZmaxRfg-ZminRfg)/NzIntvRfg,&
!            Ezgrid(NxIntvRfg/2+1,NyIntvRfg/2+1,i)
!          enddo
!        endif
!1000     format(2(1x,1pe15.8))
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read4_Data

        !read in the discrete wake field data in z, x, y.
        !distribution along axis zdat from files "rfdatax or rfdataxx 
        !or rfdataxxx".
        subroutine read1wk_Data(ifile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifile
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1
        character*10 name1
        character*11 name2
        character*12 name3

        name1 = 'rfdatax.in'
        name2 = 'rfdataxx.in'
        name3 = 'rfdataxxx.in'

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        if(myrank.eq.0) then
          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
!            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
            n = n + 1
            zdatwk(n) = tmp1
            edatwk(n) = tmp2
            epdatwk(n) = tmp3
            eppdatwk(n) = tmp4
          goto 50
77        continue
          close(14)
!          close(15)
          Ndatawk = n
          zdat1 = zdatwk(1)
          do i = 1, Ndatawk
            zdatwk(i) = zdatwk(i) - zdat1
          enddo
        endif
100     format(6(1x,1pe15.8))

        call MPI_BCAST(Ndatawk,1,MPI_INTEGER,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zdatwk,Ndatawk,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edatwk,Ndatawk,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(epdatwk,Ndatawk,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)
        call MPI_BCAST(eppdatwk,Ndatawk,MPI_DOUBLE_PRECISION,0,&
             MPI_COMM_WORLD,ierr)

        !print*,"Ndata: ",Ndata
        end subroutine read1wk_Data

      end module Dataclass
