!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in 
!                    Beam module of APPLICATION layer.
! Version: 2.1
! Author: Ji Qiang, LBNL
! Description: This class defines initial distributions for the charged 
!              particle beam bunch information in the accelerator.
! Comments: we have added three attributes to each particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Distributionclass
        use Pgrid2dclass
        use CompDomclass
        use BeamBunchclass      
        use Timerclass
        use NumConstclass
        use PhysConstclass
      contains
        ! sample the particles with intial distribution.
        subroutine sample_Dist(this,distparam0,nparam,flagdist,geom,grid,Flagbc,&
                               nchrg,nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nparam,Flagbc,nchrg
        double precision, dimension(nparam) :: distparam,distparam0
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist
        type (BeamBunch), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: flagdist
        integer :: myid, myidx, myidy,seedsize,i,isize
        !integer seedarray(1)
        integer, allocatable, dimension(:) :: seedarray
        real rancheck
        real*8 :: gam,gambet,alpha,beta,eps

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(100001+myid)*(myid+7)
!        call random_seed(put=seedarray(1:1))
!        write(6,*)'seedarray=',seedarray

        call random_seed(SIZE=seedsize)
        allocate(seedarray(seedsize))
        do i = 1, seedsize
!          seedarray(i) = (1000+5*myid)*(myid+7)+i-1 !//original one.
          !//2nd group
!          seedarray(i) = (2000+5*myid)*(myid+7)+i-1
          !//3rd group
          !seedarray(i) = (3000+5*myid)*(myid+7)+i-1
          seedarray(i) = (3001+5001*myid)*(myid+7)+i-1
        enddo
        call random_seed(PUT=seedarray)
        call random_number(rancheck)
        !the following is added new for 2nd random group ....
        do i = 1, 3000
          call random_number(rancheck)
        enddo
!        write(6,*)'myid,rancheck=',seedarray,myid,rancheck

        distparam = distparam0

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0d0)
!X-Px
        alpha = distparam0(1)     
        beta = distparam0(2)
        eps = distparam0(3)/gambet !unnormalized emittance m-rad
        distparam(1) = sqrt(beta*eps/(1.d0+alpha**2))/Scxl
        distparam(2) = sqrt(eps/beta)*gambet 
        distparam(3) = alpha/sqrt(1.d0+alpha**2)
!Y-Py
        alpha = distparam0(8)     
        beta = distparam0(9)
        eps = distparam0(10)/gambet
        distparam(8) = sqrt(beta*eps/(1.d0+alpha**2))/Scxl
        distparam(9) = sqrt(eps/beta)*gambet 
        distparam(10) = alpha/sqrt(1.d0+alpha**2)
!T-Pt
        alpha = distparam0(15)     
        beta = distparam0(16) !degree/MeV
        eps = distparam0(17) !degree-MeV
        distparam(15) = sqrt(beta*eps/(1.d0+alpha**2))/Rad2deg
        distparam(16) = sqrt(eps/beta)/(this%Mass/1.0d6)
        distparam(17) = -alpha/sqrt(1.d0+alpha**2)

        if(flagdist.eq.1) then
          call Uniform_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.2) then
          call Gauss3_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.3) then
          call Waterbag_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.4) then
          call Semigauss_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.5) then
          call KV3d_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.16) then
          call WaterbagMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.17) then
          call GaussMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.19) then
          call readin_Dist(this)
        else if(flagdist.eq.22) then
          call readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.35) then
          call readimpt_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else
          print*,"Initial distribution not available!!"
          stop
        endif

        deallocate(seedarray)

        end subroutine sample_Dist
       
        subroutine Gauss3_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,1).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
!        intvsamp = 10
        intvsamp = avgpts
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3_Dist

        subroutine normdv(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(2), intent(out) :: y
        double precision :: twopi,x1,x2,epsilon

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
10      call random_number(x1)
!        x1 = 0.5
!10      x2 = 0.6
        if(x1.eq.0.0) goto 10
!        if(x1.eq.0.0) x1 = epsilon
        y(1) = sqrt(-2.0*log(x1))*cos(twopi*x2)
        y(2) = sqrt(-2.0*log(x1))*sin(twopi*x2)

        end subroutine normdv

        subroutine normVec(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
        call random_number(x1)
        do i = 1, num
          if(x1(i).eq.0.0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
          y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
        enddo

        end subroutine normVec

        ! sample the particles with intial distribution 
        ! using rejection method. 
        subroutine Waterbag_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call random_number(x11)
        !print*,"x11: ",myid,x11

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
        !print*,"avgpts: ",avgpts
       
!        print*,avgpts,isamz,isamy
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Waterbag_Dist

        subroutine KV3d_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,j,pid
!        integer seedarray(1)
        double precision :: t0,x11,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        twopi = 4*asin(1.0)

        do numpts = 1, avgpts
          call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          x1 = 2*r4*cos(r2)
          x2 = 2*r4*sin(r2)
          x3 = 2*r5*cos(r3)
          x4 = 2*r5*sin(r3)
!x-px:
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation.
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          call random_number(r5)
          r5 = 2*r5 - 1.0
          call random_number(r6)
          r6 = 2*r6 - 1.0
          x5 = r5*sqrt(3.0)
          x6 = r6*sqrt(3.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo
          
        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KV3d_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Uniform_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(6,1) :: a
        double precision, dimension(2) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        this%Pts1 = 0.0
    
        do ii = 1, this%Npt
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(1,1) = xmu1 + sig1*r1/sq12
          a(2,1) = xmu2 + sig2*(-muxpx*r2/sq12+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(3,1) = xmu3 + sig3*r1/sq34
          a(4,1) = xmu4 + sig4*(-muypy*r2/sq34+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(5,1) = xmu5 + sig5*r1/sq56
          a(6,1) = xmu6 + sig6*(-muzpz*r2/sq56+r2)

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
          endif

          enddo
        enddo
          
!        call MPI_BARRIER(comm2d,ierr)

        this%Nptlocal = numpts0
       
!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo


        end subroutine Uniform_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Semigauss_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(6,2) :: a
        double precision, dimension(3) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2, r, r3
        double precision,allocatable,dimension(:,:) :: ptstmp
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
! The performance inside this loop might be improved due to
! a lot of subroutine call in this loop.    
        do ii = 1, this%Npt
          ! rejection sample.
10        call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r1 = 2.0*r1-1.0
          r2 = 2.0*r2-1.0
          r3 = 2.0*r3-1.0
          if(r1**2+r2**2+r3**2.gt.1.0) goto 10
          x2(1) = r1
          x2(2) = r2
          x2(3) = r3
          call normdv2(x1)

          !x-px:
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x2(1)/sq12*sqrt(5.0)
          a(2,1) = xmu2 + sig2*(-muxpx*x2(1)/sq12+x1(1))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x2(1)*sqrt(5.0)
          !a(2,1) = xmu2 + sig2*(muxpx*x2(1)+sq12*x1(1))
          !y-py
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x2(2)/sq34*sqrt(5.0)
          a(4,1) = xmu4 + sig4*(-muypy*x2(2)/sq34+x1(2))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x2(2)*sqrt(5.0)
          !a(4,1) = xmu4 + sig4*(muypy*x2(2)+sq34*x1(2))
          !z-pz
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x2(3)/sq56*sqrt(5.0)
          a(6,1) = xmu6 + sig6*(-muzpz*x2(3)/sq56+x1(3))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x2(3)*sqrt(5.0)
          !a(6,1) = xmu6 + sig6*(muzpz*x2(3)+sq56*x1(3))

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
            this%Pts1(7,numpts0) = this%Charge/this%mass
            this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
            this%Pts1(9,numpts0) = ii

          endif

          enddo
        enddo
          
!        call MPI_BARRIER(comm2d,ierr)
        allocate(ptstmp(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            ptstmp(j,i0) = this%Pts1(j,i0)
          enddo
        enddo
        deallocate(this%Pts1)
        allocate(this%Pts1(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            this%Pts1(j,i0) = ptstmp(j,i0)
          enddo
        enddo
        deallocate(ptstmp)

        this%Nptlocal = numpts0
!        print*,"numpts0: ",numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Semigauss_Dist

        subroutine normdv2(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(3), intent(out) :: y
        double precision :: sumtmp,x
        integer :: i

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(1) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(2) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(3) = sumtmp - 6.0

        end subroutine normdv2

        ! sample the particles with intial distribution 
        ! using rejection method for multi-charge state. 
        subroutine WaterbagMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j,i,ii,kk
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call random_number(x11)
        !print*,"x11: ",myid,x11

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
        !print*,"avgpts: ",avgpts
       
!        print*,avgpts,isamz,isamy
        
        ii = 0
        avgpts = sum(nptlist)/totnp
        do i = 1, nchrg
          kk = nptlist(i)/totnp
          do j = 1, kk
            ii = ii + 1
            this%Pts1(7,ii) = qmcclist(i)
            this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
            pid = ii + myid*totnp
            this%Pts1(9,ii) = pid
          enddo
        enddo

!        do j = 1, avgpts
!          pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
!          this%Pts1(9,j) = pid
!        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine WaterbagMC_Dist

        subroutine GaussMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer, dimension(nchrg) :: nptlist
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid,ii,kk
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

!        if(mod(avgpts,10).ne.0) then
!          print*,"The number of particles has to be an integer multiple of 10Nprocs"
!          stop
!        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
!        intvsamp = 10
        intvsamp = avgpts
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        ii = 0
        avgpts = sum(nptlist)/totnp
        do i = 1, nchrg
          kk = nptlist(i)/totnp
          do j = 1, kk
            ii = ii + 1
            this%Pts1(7,ii) = qmcclist(i)
            this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
            pid = ii + myid*totnp
            this%Pts1(9,ii) = pid
          enddo
        enddo
 
!        do j = 1, avgpts
!          pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
!          this%Pts1(9,j) = pid
!        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussMC_Dist

        subroutine readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        real*8 :: beta0,sumeng
        integer :: jlow,jhigh,nleft,avgpts

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='particle.in',status='old')
        read(12,*)inipts
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(7,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        sumeng = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(2,i)
          sumx2 = sumx2 + Ptcl(2,i)*Ptcl(2,i)
          sumy = sumy + Ptcl(4,i)
          sumy2 = sumy2 + Ptcl(4,i)*Ptcl(4,i)
          sumeng = sumeng + Ptcl(7,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        sumeng = sumeng/inipts
        !print*,"sumx2: ",sumx2,sumy2,sumeng
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+sumeng*1.0e6/mccq  
        gamma0 = -this%refptcl(6) 
        !print*,"gamma0: ",gamma0
        beta0 = sqrt(1.0-1./gamma0**2)
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0

        avgpts = inipts/totnp
        nleft = nptot - avgpts*totnp
        if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
        else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
        endif

        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        do j = 1, inipts
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
            i = j - jlow + 1
            this%Pts1(1,i) = Ptcl(2,j)*xscale 
            this%Pts1(3,i) = Ptcl(4,j)*yscale 
            this%Pts1(5,i) = Ptcl(6,j)*zscale 
            gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
            this%Pts1(2,i) = Ptcl(3,j)*gammabet*pxscale + xmu2
            this%Pts1(4,i) = Ptcl(5,j)*gammabet*pyscale + xmu4
            !this%Pts1(6,i) = gammabet*pzscale + xmu6
            this%Pts1(6,i) = gamma0 - sqrt(1.0 + Ptcl(7,j)**2) 
            this%Pts1(9,i) = j
          endif
        enddo

        this%Nptlocal = avgpts

        do i = 1, this%Nptlocal
          this%Pts1(1,i) = this%Pts1(1,i)/Scxl + xmu1
          !this%Pts1(2,i) = 0.0
          this%Pts1(3,i) = this%Pts1(3,i)/Scxl + xmu3
          !this%Pts1(4,i) = 0.0
          this%Pts1(5,i) = this%Pts1(5,i)*Scfreq*2*Pi + xmu5
          !!this%Pts1(6,i) = xmu6
        enddo

        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
        enddo

        end subroutine readElegant_Dist

        subroutine readin_Dist(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nptot,nleft
        double precision, dimension(9) :: tmptcl
        double precision :: sum1,sum2,x1,x2
 
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
 
        open(unit=12,file='particle.in',status='old')
 
        sum1 = 0.0
        sum2 = 0.0
 
          read(12,*)nptot
          avgpts = nptot/nproc
          nleft = nptot - avgpts*nproc
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
          !jlow = myid*avgpts + 1
          !jhigh = (myid+1)*avgpts
          !print*,"avgpts, jlow, and jhigh: ",avgpts,jlow,jhigh
          do j = 1, nptot
            read(12,*)tmptcl(1:9)
            sum1 = sum1 + tmptcl(1)
            sum2 = sum2 + tmptcl(3)
            if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              this%Pts1(1:9,i) = tmptcl(1:9)
            endif
!            if(myid.eq.0) print*,i,sum1,sum2
          enddo
          !print*,"sumx1,sumy1: ",sum1/nptot,sum2/nptot
 
          close(12)
 
          this%Nptlocal = avgpts
 
        end subroutine readin_Dist

        !read particle distribution from ImpactT output
        subroutine readimpt_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sum1,sum2
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(6) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz,phi0lc
        real*8 :: xmax,ymax,zmax,pxmax,pymax,pzmax

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        phi0lc = this%refptcl(5)
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='particle.in',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(6,avgpts))

        xmax = 0.0
        ymax = 0.0
        zmax = 0.0
        pxmax = 0.0
        pymax = 0.0
        pzmax = 0.0

        sum1 = 0.0d0
        sum2 = 0.0d0
        do j = 1, inipts
          read(12,*)tmptcl(1:6)
          sum1 = sum1 + tmptcl(5)
          gamma = sqrt(1.0+tmptcl(2)**2+tmptcl(4)**2+tmptcl(6)**2)
          sum2 = sum2 + gamma
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:6,i) = tmptcl(1:6)
          endif
        enddo
        sum1 = sum1/inipts
        sum2 = sum2/inipts
        print*,"mean z location and energy: ",sum1, sum2
        close(12)
        call MPI_BARRIER(comm2d,ierr)
        this%refptcl(6) = -sum2 

        do j = 1, avgpts
          Ptcl(1,j) = (Ptcl(1,j)/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/xl)*yscale + xmu3
          Ptcl(5,j) = -((Ptcl(5,j)-sum1)/(xl*sqrt(1.d0-1.d0/sum2**2)))*zscale + xmu5
          gamma = sqrt(1.0+Ptcl(2,j)**2+Ptcl(4,j)**2+Ptcl(6,j)**2)
          Ptcl(6,j) = (sum2-gamma)*pzscale + xmu6
          Ptcl(2,j) = Ptcl(2,j)*pxscale + xmu2
          Ptcl(4,j) = Ptcl(4,j)*pyscale + xmu4
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                r = (2.0*r-1.0)*0.03*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                r = (2.0*r-1.0)*0.03*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.005*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                r = (2.0*r-1.0)*0.01*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                !r = (2.0*r-1.0)*0.004*pzmax
                r = (2.0*r-1.0)*0.01*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j + jlow
        enddo

        end subroutine readimpt_Dist

      end module Distributionclass
