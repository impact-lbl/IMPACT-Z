!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer. 
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments: We have added 3 more attributes to the particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Outputclass
!        use Timer_class
        use BeamBunchclass
        use Pgrid2dclass
        use Fldmgerclass
        use PhysConstclass
      contains
        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnostic1_Output(z,this,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: nchrg
        integer, dimension(nchrg), intent(in) :: nptlist
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: xpx2local,xpx3local,x2pxlocal,x2px2local,x3pxlocal
        double precision:: ypy2local,ypy3local,y2pylocal,y2py2local,y3pylocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: cubx2px,cubxpx2,fthx3px,fthx2px2,fthxpx3
        double precision:: x3px,xpx3,x2px2
        double precision:: I2x,I4x,Hx
        double precision:: cuby2py,cubypy2,fthy3py,fthy2py2,fthypy3
        double precision:: y3py,ypy3,y2py2
        double precision:: I2y,I4y,Hy
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(37) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
!        double precision :: alphax,alphay,alphaz
!       for the calculation of bunched beam current.
!        double precision, dimension(3) :: localaa, aa, localbb, bb 
!        double precision :: pi,tmp5,tmp55
!        integer :: ntmp5,nf
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0d6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

        den1 = 1.0d0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        xpx2local  = 0.0
        xpx3local  = 0.0
        x2pxlocal  = 0.0
        x2px2local = 0.0
        x3pxlocal  = 0.0
        ypy2local  = 0.0
        ypy3local  = 0.0
        y2pylocal  = 0.0
        y2py2local = 0.0
        y3pylocal  = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

!        localaa = 0.0
!        localbb = 0.0
!        aa = 0.0
!        bb = 0.0
!        pi = 2*asin(1.0)

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          xpx2local  = xpx2local + this%Pts1(1,i)*this%Pts1(2,i)*&
                       this%Pts1(2,i)
          xpx3local  = xpx3local + this%Pts1(1,i)*this%Pts1(2,i)*&
                       this%Pts1(2,i)*this%Pts1(2,i)
          x2pxlocal  = x2pxlocal + this%Pts1(1,i)*this%Pts1(1,i)*&
                       this%Pts1(2,i)
          x2px2local = x2px2local + this%Pts1(1,i)*this%Pts1(1,i)*&
                       this%Pts1(2,i)*this%Pts1(2,i)
          x3pxlocal  = x3pxlocal + this%Pts1(1,i)*this%Pts1(1,i)*&
                       this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          ypy2local  = ypy2local + this%Pts1(3,i)*this%Pts1(4,i)*&
                       this%Pts1(4,i)
          ypy3local  = ypy3local + this%Pts1(3,i)*this%Pts1(4,i)*&
                       this%Pts1(4,i)*this%Pts1(4,i)
          y2pylocal  = y2pylocal + this%Pts1(3,i)*this%Pts1(3,i)*&
                       this%Pts1(4,i)
          y2py2local = y2py2local + this%Pts1(3,i)*this%Pts1(3,i)*&
                       this%Pts1(4,i)*this%Pts1(4,i)
          y3pylocal  = y3pylocal + this%Pts1(3,i)*this%Pts1(3,i)*&
                       this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = xpx2local
        tmplc(29) = xpx3local
        tmplc(30) = x2pxlocal
        tmplc(31) = x2px2local
        tmplc(32) = x3pxlocal
        tmplc(33) = ypy2local
        tmplc(34) = ypy3local
        tmplc(35) = y2pylocal
        tmplc(36) = y2py2local
        tmplc(37) = y3pylocal
        
        call MPI_REDUCE(tmplc,tmpgl,37,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          cubxpx2 = tmpgl(28)*den1
          cubx2px = tmpgl(30)*den1
          fthx = tmpgl(17)*den1
          fthx2px2 = tmpgl(31)*den1
          fthx3px  = tmpgl(32)*den1
          fthxpx3  = tmpgl(29)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          cubypy2 = tmpgl(33)*den1
          cuby2py = tmpgl(35)*den1
          fthy2py2 = tmpgl(36)*den1
          fthy3py  = tmpgl(37)*den1
          fthypy3  = tmpgl(34)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          x2px2 =(fthx2px2+px0*px0*sqx+4*x0*px0*(xpx+x0*px0)+x0*x0*sqpx-&
            2*px0*cubx2px-2*x0*cubxpx2-3*x0*x0*px0*px0)
          xpx3  = (fthxpx3+3*px0*px0*(xpx+x0*px0)+3*x0*px0*sqpx-&
                  x0*cubpx-3*px0*cubxpx2-3*x0*px0*px0*px0)
          x3px  = (fthx3px+3*x0*px0*sqx+3*x0*x0*(xpx+x0*px0)-&
                   px0*cubx-3*x0*cubx2px-3*x0*x0*x0*px0) 
          y2py2 = fthy2py2-2*py0*cuby2py+py0*py0*sqy-2*y0*cubypy2+&
                  4*y0*py0*(ypy+y0*py0)+y0*y0*sqpy-3*y0*y0*py0*py0
          ypy3  = fthypy3-3*py0*cubypy2+3*py0*py0*(ypy+y0*py0)-&
                  y0*cubpy+3*y0*py0*sqpy-3*y0*py0*py0*py0
          y3py  = fthy3py-py0*cuby-3*y0*cuby2py+3*y0*py0*sqy+&
                  3*y0*y0*(ypy+y0*py0)-3*y0*y0*y0*py0
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0/gambet,pxrms/gambet,-xpx/epx,epx*xl
          write(25,100)z,y0*xl,yrms*xl,py0/gambet,pyrms/gambet,-ypy/epy,epy*xl
          write(26,100)z,z0*xt,zrms*xt,pz0*qmc,pzrms*qmc,-zpz/epz,epz*qmc*xt
          write(27,100)z,glmax(1)*xl,glmax(2)/gambet,glmax(3)*xl,&
                       glmax(4)/gambet,glmax(5)*xt,glmax(6)*qmc
          write(28,101)z,npctmin,npctmax,nptot
          write(29,100)z,x03*xl,px03/gambet,y03*xl,py03/gambet,z03*xt,&
                       pz03*qmc
          write(30,100)z,x04*xl,px04/gambet,y04*xl,py04/gambet,z04*xt,&
                       pz04*qmc
          write(32,*)z,nptlist(1:nchrg)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
          call flush(32)
        endif

99      format(8(1x,e13.6))
100      format(7(1x,e13.6))
101     format(1x,e13.6,3I10)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1_Output

        subroutine diagnostic2_Output(this,z,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        integer, intent(in) :: nchrg
        integer, dimension(nchrg), intent(in) :: nptlist
        integer, parameter :: nbin=200
        integer, dimension(nbin) :: tmpbin,glbin
        integer :: innp,i,my_rank,ierr,j,npctmin,npctmax,nfreq
        integer :: nii,iitmp,ip1,ip2,ip3
        integer :: nptot
        double precision :: qmc,qmass,xl,xt,tg,pi,t0,gb,bgend,blg
        double precision :: ez1,ezp1,ezpp1,gam,energy,bet,gambet
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                           epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                           xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local,x0lc,x0,px0lc,px0
        double precision:: y0lc,y0,py0lc,py0,z0lc,z0,pz0lc,pz0
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(15) :: tmplc,tmpgl
        double precision :: lcrmax,glrmax,tmp1,tmp2
        double precision :: f90,f95,f99,ex90,ex95,ex99,ey90,ey95,ey99,&
        ez90,ez95,ez99,r90,r95,r99,rrms,ravg,rrmax,rtmp2,rrmaxlc,ravg2,&
        ravglc,ravg2lc
        double precision, dimension(3) :: Ealpha,Ebeta,Egamma
        double precision :: epsmxlc,epsmylc,epsmzlc,xtmp,pxtmp,ytmp,pytmp,&
                   ztmp,pztmp,epsmx,epsmy,epsmz,eps,hxeps,hyeps,hzeps,hreps
        double precision,allocatable,dimension(:) :: epsiontmp
        real*8 :: sqsum7,sqsum7local,sqsum8,sqsum8local,sqsum9,sqsum9local

        call starttime_Timer(t0)

        qmc = this%Mass/1.0d6
        qmass =this%Mass
        xl = Scxl
        xt = Rad2deg
        tg = this%refptcl(5)
        pi = 2*asin(1.0d0)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt

        den1 = 1.0d0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          tmp1 = this%Pts1(2,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*tmp1
          px0lc = px0lc + tmp1
          sqsum2local = sqsum2local + tmp1*tmp1
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          tmp2 = this%Pts1(4,i)
          ypylocal = ypylocal + this%Pts1(3,i)*tmp2
          py0lc = py0lc + tmp2
          sqsum4local = sqsum4local + tmp2*tmp2
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        
        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        bet = sqrt(gam**2-1.0)/gam
        gambet = sqrt(gam**2-1.0)
        energy = (gam-1.)*qmc

        x0 = tmpgl(1)
        px0 = tmpgl(2)
        y0 = tmpgl(3)
        py0 = tmpgl(4)
        z0 = tmpgl(5)
        pz0 = tmpgl(6)
        sqsum1 = tmpgl(7) - x0*x0*den1
        sqsum2 = tmpgl(8) - px0*px0*den1
        sqsum3 = tmpgl(9) - y0*y0*den1
        sqsum4 = tmpgl(10) - py0*py0*den1
        sqsum5 = tmpgl(11) - z0*z0*den1
        sqsum6 = tmpgl(12) - pz0*pz0*den1
        xpx = tmpgl(13) - x0*px0*den1
        ypy = tmpgl(14) - y0*py0*den1
        zpz = tmpgl(15) - z0*pz0*den1
        epsx2 = (sqsum1*sqsum2-xpx*xpx)*den2
        epsy2 = (sqsum3*sqsum4-ypy*ypy)*den2
        epsz2 = (sqsum5*sqsum6-zpz*zpz)*den2
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1*den1)
        pxrms = sqrt(sqsum2*den1)
        yrms = sqrt(sqsum3*den1)
        pyrms = sqrt(sqsum4*den1)
        zrms = sqrt(sqsum5*den1)
        pzrms = sqrt(sqsum6*den1)
        xpx = xpx*den1
        ypy = ypy*den1
        zpz = zpz*den1
        xpxfac = 0.0
        ypyfac = 0.0
        zpzfac = 0.0
        if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
        if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
        if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
        x0 = x0*den1
        px0 = px0*den1
        y0 = y0*den1
        py0 = py0*den1
        z0 = z0*den1
        pz0 = pz0*den1

        Ealpha(1) = -xpx/epx
        Ealpha(2) = -ypy/epy
        Ealpha(3) = -zpz/epz
        Ebeta(1) = xrms*xrms*gambet/epx
        Ebeta(2) = yrms*yrms*gambet/epy
        Ebeta(3) = zrms*zrms*bet*bet/(epz/gam**3/bet)
        Egamma(:) = (1.0+Ealpha(:)*Ealpha(:))/Ebeta(:)

        allocate(epsiontmp(innp))
        epsmxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
          tmp1 = this%Pts1(2,i)
          pxtmp = (tmp1 - px0)/gambet
          epsiontmp(i)=Egamma(1)*xtmp*xtmp+2*Ealpha(1)*xtmp*pxtmp+&
                     Ebeta(1)*pxtmp*pxtmp
          if(epsmxlc.le.epsiontmp(i)) epsmxlc = epsiontmp(i)
        enddo

        call MPI_ALLREDUCE(epsmxlc,epsmx,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        eps = 1.0e-8
        hxeps = (epsmx+eps)/nbin

        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          if (.not. isnan(epsiontmp(i))) then
            iitmp = epsiontmp(i)/hxeps
            nii = iitmp+1
            tmpbin(nii) = tmpbin(nii) + 1 
          endif
        enddo

        glbin = 0
        call MPI_ALLREDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        f90 = 0.9d0*nptot
        f95 = 0.95d0*nptot
        f99 = 0.99d0*nptot
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ex90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i1: ",i,nbin,glbin(i-1),glbin(i),f90,f95,f99
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ex95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i2: ",i,nbin,glbin(i-1),glbin(i)
          do i =1, nbin
            if(glbin(i).gt.f99) then
              ex99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i3: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hxeps: ",glbin(nbin),nptot,epsmx,ex90,ex95,ex99
          !compute the 90%,95%,99% rms emittance
          sqsum1local = 0.0d0
          sqsum2local = 0.0d0
          sqsum3local = 0.0d0
          sqsum4local = 0.0d0
          sqsum5local = 0.0d0
          sqsum6local = 0.0d0
          sqsum7local = 0.0d0
          sqsum8local = 0.0d0
          sqsum9local = 0.0d0
          ip1 = 0
          ip2 = 0
          ip3 = 0
          do i = 1, innp
            if(epsiontmp(i).le.ex90) then
              sqsum1local = sqsum1local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(1,i)-x0)
              sqsum2local = sqsum2local + (this%Pts1(2,i)-px0)*&
                                          (this%Pts1(2,i)-px0)
              sqsum3local = sqsum3local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(2,i)-px0)
              ip1 = ip1 + 1
            else if(epsiontmp(i).le.ex95) then
              sqsum4local = sqsum4local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(1,i)-x0)
              sqsum5local = sqsum5local + (this%Pts1(2,i)-px0)*&
                                          (this%Pts1(2,i)-px0)
              sqsum6local = sqsum6local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(2,i)-px0)
              ip2 = ip2 + 1
            else if(epsiontmp(i).le.ex99) then
              sqsum7local = sqsum7local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(1,i)-x0)
              sqsum8local = sqsum8local + (this%Pts1(2,i)-px0)*&
                                          (this%Pts1(2,i)-px0)
              sqsum9local = sqsum9local + (this%Pts1(1,i)-x0)*&
                                          (this%Pts1(2,i)-px0)
              ip3 = ip3 + 1
            endif 
          enddo
          sqsum4local = sqsum4local + sqsum1local
          sqsum5local = sqsum5local + sqsum2local
          sqsum6local = sqsum6local + sqsum3local
          ip2 = ip2 + ip1
          sqsum7local = sqsum4local + sqsum7local
          sqsum8local = sqsum5local + sqsum8local
          sqsum9local = sqsum6local + sqsum9local
          ip3 = ip3 + ip2
          tmplc(1) = sqsum1local
          tmplc(2) = sqsum2local
          tmplc(3) = sqsum3local
          tmplc(4) = sqsum4local
          tmplc(5) = sqsum5local
          tmplc(6) = sqsum6local
          tmplc(7) = sqsum7local
          tmplc(8) = sqsum8local
          tmplc(9) = sqsum9local
          tmplc(10) = ip1
          tmplc(11) = ip2
          tmplc(12) = ip3
          tmpgl = 0.0d0
          call MPI_REDUCE(tmplc,tmpgl,12,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          sqsum1 = tmpgl(1)/tmpgl(10)
          sqsum2 = tmpgl(2)/tmpgl(10)
          sqsum3 = tmpgl(3)/tmpgl(10)
          sqsum4 = tmpgl(4)/tmpgl(11)
          sqsum5 = tmpgl(5)/tmpgl(11)
          sqsum6 = tmpgl(6)/tmpgl(11)
          sqsum7 = tmpgl(7)/tmpgl(12)
          sqsum8 = tmpgl(8)/tmpgl(12)
          sqsum9 = tmpgl(9)/tmpgl(12)
          ex90 = sqrt(sqsum1*sqsum2-sqsum3**2)
          ex95 = sqrt(sqsum4*sqsum5-sqsum6**2)
          ex99 = sqrt(sqsum7*sqsum8-sqsum9**2)
        endif

        epsmylc = -1.0e10
        do i = 1, innp
          ytmp = this%Pts1(3,i) - y0
!          tmp2 = this%Pts1(4,i)-(xl*xl*ezp1*this%Pts1(3,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp2 = this%Pts1(4,i)
          pytmp = (tmp2 - py0)/gambet
          epsiontmp(i)=Egamma(2)*ytmp*ytmp+2*Ealpha(2)*ytmp*pytmp+&
                     Ebeta(2)*pytmp*pytmp
          if(epsmylc.le.epsiontmp(i)) epsmylc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmylc,epsmy,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hyeps = (epsmy+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          if (.not. isnan(epsiontmp(i))) then
            iitmp = epsiontmp(i)/hyeps
            nii = iitmp+1
            tmpbin(nii) = tmpbin(nii) + 1
          endif
        enddo
        glbin = 0
        call MPI_ALLREDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ey90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i4: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ey95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i5: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ey99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !compute the 90%,95%,99% rms emittance
          sqsum1local = 0.0d0
          sqsum2local = 0.0d0
          sqsum3local = 0.0d0
          sqsum4local = 0.0d0
          sqsum5local = 0.0d0
          sqsum6local = 0.0d0
          sqsum7local = 0.0d0
          sqsum8local = 0.0d0
          sqsum9local = 0.0d0
          ip1 = 0
          ip2 = 0
          ip3 = 0
          do i = 1, innp
            if(epsiontmp(i).le.ey90) then
              sqsum1local = sqsum1local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(3,i)-y0)
              sqsum2local = sqsum2local + (this%Pts1(4,i)-py0)*&
                                          (this%Pts1(4,i)-py0)
              sqsum3local = sqsum3local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(4,i)-py0)
              ip1 = ip1 + 1
            else if(epsiontmp(i).le.ey95) then
              sqsum4local = sqsum4local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(3,i)-y0)
              sqsum5local = sqsum5local + (this%Pts1(4,i)-py0)*&
                                          (this%Pts1(4,i)-py0)
              sqsum6local = sqsum6local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(4,i)-py0)
              ip2 = ip2 + 1
            else if(epsiontmp(i).le.ey99) then
              sqsum7local = sqsum7local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(3,i)-y0)
              sqsum8local = sqsum8local + (this%Pts1(4,i)-py0)*&
                                          (this%Pts1(4,i)-py0)
              sqsum9local = sqsum9local + (this%Pts1(3,i)-y0)*&
                                          (this%Pts1(4,i)-py0)
              ip3 = ip3 + 1
            endif 
          enddo
          sqsum4local = sqsum4local + sqsum1local
          sqsum5local = sqsum5local + sqsum2local
          sqsum6local = sqsum6local + sqsum3local
          ip2 = ip2 + ip1
          sqsum7local = sqsum4local + sqsum7local
          sqsum8local = sqsum5local + sqsum8local
          sqsum9local = sqsum6local + sqsum9local
          ip3 = ip3 + ip2
          tmplc(1) = sqsum1local
          tmplc(2) = sqsum2local
          tmplc(3) = sqsum3local
          tmplc(4) = sqsum4local
          tmplc(5) = sqsum5local
          tmplc(6) = sqsum6local
          tmplc(7) = sqsum7local
          tmplc(8) = sqsum8local
          tmplc(9) = sqsum9local
          tmplc(10) = ip1
          tmplc(11) = ip2
          tmplc(12) = ip3
          tmpgl = 0.0d0
          call MPI_REDUCE(tmplc,tmpgl,12,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          sqsum1 = tmpgl(1)/tmpgl(10)
          sqsum2 = tmpgl(2)/tmpgl(10)
          sqsum3 = tmpgl(3)/tmpgl(10)
          sqsum4 = tmpgl(4)/tmpgl(11)
          sqsum5 = tmpgl(5)/tmpgl(11)
          sqsum6 = tmpgl(6)/tmpgl(11)
          sqsum7 = tmpgl(7)/tmpgl(12)
          sqsum8 = tmpgl(8)/tmpgl(12)
          sqsum9 = tmpgl(9)/tmpgl(12)
          ey90 = sqrt(sqsum1*sqsum2-sqsum3**2)
          ey95 = sqrt(sqsum4*sqsum5-sqsum6**2)
          ey99 = sqrt(sqsum7*sqsum8-sqsum9**2)
        endif


        epsmzlc = -1.0e10
        do i = 1, innp
          ztmp = (this%Pts1(5,i) - z0)*bet
          pztmp = (this%Pts1(6,i) - pz0)/gam**3/bet**2
          epsiontmp(i)=Egamma(3)*ztmp*ztmp+2*Ealpha(3)*ztmp*pztmp+&
                     Ebeta(3)*pztmp*pztmp
          if(epsmzlc.le.epsiontmp(i)) epsmzlc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmzlc,epsmz,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hzeps = (epsmz+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          if (.not. isnan(epsiontmp(i))) then
            iitmp = epsiontmp(i)/hzeps
            nii = iitmp+1
            tmpbin(nii) = tmpbin(nii) + 1
          endif
        enddo
        glbin = 0
        call MPI_ALLREDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)

          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ez90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i7: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ez95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i8: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ez99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i9: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hzeps: ",glbin(nbin),nptot,epsmz,ez90,ez95,ez99
          !ez90 = ez90*gam**3*bet
          !ez95 = ez95*gam**3*bet
          !ez99 = ez99*gam**3*bet

          !compute the 90%,95%,99% rms emittance
          sqsum1local = 0.0d0
          sqsum2local = 0.0d0
          sqsum3local = 0.0d0
          sqsum4local = 0.0d0
          sqsum5local = 0.0d0
          sqsum6local = 0.0d0
          sqsum7local = 0.0d0
          sqsum8local = 0.0d0
          sqsum9local = 0.0d0
          ip1 = 0
          ip2 = 0
          ip3 = 0
          do i = 1, innp
            if(epsiontmp(i).le.ez90) then
              sqsum1local = sqsum1local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(5,i)-z0)
              sqsum2local = sqsum2local + (this%Pts1(6,i)-pz0)*&
                                          (this%Pts1(6,i)-pz0)
              sqsum3local = sqsum3local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(6,i)-pz0)
              ip1 = ip1 + 1
            else if(epsiontmp(i).le.ez95) then
              sqsum4local = sqsum4local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(5,i)-z0)
              sqsum5local = sqsum5local + (this%Pts1(6,i)-pz0)*&
                                          (this%Pts1(6,i)-pz0)
              sqsum6local = sqsum6local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(6,i)-pz0)
              ip2 = ip2 + 1
            else if(epsiontmp(i).le.ez99) then
              sqsum7local = sqsum7local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(5,i)-z0)
              sqsum8local = sqsum8local + (this%Pts1(6,i)-pz0)*&
                                          (this%Pts1(6,i)-pz0)
              sqsum9local = sqsum9local + (this%Pts1(5,i)-z0)*&
                                          (this%Pts1(6,i)-pz0)
              ip3 = ip3 + 1
            endif 
          enddo
          sqsum4local = sqsum4local + sqsum1local
          sqsum5local = sqsum5local + sqsum2local
          sqsum6local = sqsum6local + sqsum3local
          ip2 = ip2 + ip1
          sqsum7local = sqsum4local + sqsum7local
          sqsum8local = sqsum5local + sqsum8local
          sqsum9local = sqsum6local + sqsum9local
          ip3 = ip3 + ip2
          tmplc(1) = sqsum1local
          tmplc(2) = sqsum2local
          tmplc(3) = sqsum3local
          tmplc(4) = sqsum4local
          tmplc(5) = sqsum5local
          tmplc(6) = sqsum6local
          tmplc(7) = sqsum7local
          tmplc(8) = sqsum8local
          tmplc(9) = sqsum9local
          tmplc(10) = ip1
          tmplc(11) = ip2
          tmplc(12) = ip3
          tmpgl = 0.0d0
          call MPI_REDUCE(tmplc,tmpgl,12,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          sqsum1 = tmpgl(1)/tmpgl(10)
          sqsum2 = tmpgl(2)/tmpgl(10)
          sqsum3 = tmpgl(3)/tmpgl(10)
          sqsum4 = tmpgl(4)/tmpgl(11)
          sqsum5 = tmpgl(5)/tmpgl(11)
          sqsum6 = tmpgl(6)/tmpgl(11)
          sqsum7 = tmpgl(7)/tmpgl(12)
          sqsum8 = tmpgl(8)/tmpgl(12)
          sqsum9 = tmpgl(9)/tmpgl(12)
          ez90 = sqrt(sqsum1*sqsum2-sqsum3**2)
          ez95 = sqrt(sqsum4*sqsum5-sqsum6**2)
          ez99 = sqrt(sqsum7*sqsum8-sqsum9**2)
        endif

        deallocate(epsiontmp)

        ravglc = 0.0
        ravg2lc = 0.0
        rrmaxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
          ytmp = this%Pts1(3,i) - y0
          rtmp2 = xtmp*xtmp + ytmp*ytmp
          ravglc = ravglc + sqrt(rtmp2)
          ravg2lc = ravg2lc + rtmp2
          if(rrmaxlc.le.rtmp2) rrmaxlc = rtmp2
        enddo
        call MPI_REDUCE(ravglc,ravg,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ravg2lc,ravg2,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(rrmaxlc,rrmax,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hreps = (sqrt(rrmax)+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = sqrt((this%Pts1(1,i)-x0)**2+(this%Pts1(3,i)-y0)**2)/hreps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              r90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i10: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              r95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i11: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              r99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i12: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hreps: ",glbin(nbin),nptot,sqrt(rrmax),r90,r95,r99
        endif

        if(my_rank.eq.0) then
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0/gambet,pxrms/gambet,-xpx/epx,&
                       epx*xl,ex90*xl,ex95*xl,ex99*xl
          write(25,100)z,y0*xl,yrms*xl,py0/gambet,pyrms/gambet,-ypy/epy,&
                       epy*xl,ey90*xl,ey95*xl,ey99*xl
          write(26,100)z,z0*xt,zrms*xt,pz0*qmc,pzrms*qmc,-zpz/epz,epz*qmc*xt,&
                       ez90*qmc*xt,ez95*qmc*xt,ez99*qmc*xt
          write(27,102)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
                     glmax(5)*xt,glmax(6)
          write(28,101)z,npctmin,npctmax,nptot
          ravg = ravg/nptot
          rrms = sqrt(ravg2/nptot - ravg*ravg)
          write(29,102)z,ravg*xl,rrms*xl,r90*xl,r95*xl,r99*xl,sqrt(rrmax)*xl
          write(32,*)z,nptlist(1:nchrg)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(32)
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,xpx*xpxfac,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,ypy*ypyfac,epy*xl
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,zpz*zpzfac,epz*qmc*xt
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl/gambet
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl/gambet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,&
!                       epz*xl/gam**3/bet
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*xl/(gam**3*bet)
!        write(24,99)z,xrms,pxrms,xpx*xpxfac,epx
!        write(25,99)z,yrms,pyrms,ypy*ypyfac,epy
!        write(26,99)z,zrms,pzrms,zpz*zpzfac,epz
!          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
!                     glmax(5)*bet*xl,glmax(6)
!        write(27,100)z,glmax(1),glmax(2),glmax(3),glmax(4), &
!                     glmax(5),glmax(6)
        endif

99      format(6(1x,e13.6))
100      format(10(1x,e13.6))
101     format(1x,e13.6,3I13)
102      format(7(1x,e13.6))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic2_Output

        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnosticT_Output(z,this)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

!        do i = 1, innp
!          write(70,999)this%Pts1(1:9,i)
!        enddo
!999     format(9(1x,e14.6))

        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        
        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          write(48,100)z,this%refptcl(1)*xl,this%refptcl(2),this%refptcl(3)*xl,&
                      this%refptcl(4),this%refptcl(5)*xl,this%refptcl(6)
          write(34,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
          write(35,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
          write(36,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl

          write(37,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(38,101)z,npctmin,npctmax,nptot
          write(39,100)z,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(40,100)z,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(48)
          call flush(34)
          call flush(35)
          call flush(36)
          call flush(37)
          call flush(38)
          call flush(39)
          call flush(40)
        endif

99      format(6(1x,e13.6))
100      format(7(1x,e13.6))
101     format(1x,e13.6,3I10)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnosticT_Output

        subroutine inpoint_Output(nfile,this,z,inb,jstp,ibal,nstep,nprocrow,nproccol,&
                                  geom,nx,ny,nz,myidx,myidy,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(out) :: this
        double precision, intent(out) :: z
        type(CompDom), intent(out) :: geom
        integer, intent(out) :: inb,jstp,ibal,nstep
        integer, intent(inout) :: nptot
        integer, intent(in) :: nprocrow,nproccol,nx,ny,nz,myidx,myidy
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: nplc,ierr

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))

        open(nfile,status="old",form="unformatted")

        read(nfile)z
        read(nfile)inb,jstp,ibal,nstep
        read(nfile)msize(1:3)
        read(nfile)range(1:6)
        read(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        read(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        read(nfile)this%Nptlocal
        read(nfile)this%refptcl(1:6)
        allocate(this%Pts1(9,this%Nptlocal))
        read(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        geom%Meshsize = msize
        geom%SpatRange = range
        allocate(geom%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
        allocate(geom%LcTabnm(2,0:nprocrow-1,0:nproccol-1))
        geom%lcTabnm = Localnum
        geom%LcTabrg = Localrange

        deallocate(Localnum)
        deallocate(Localrange)

        geom%Meshnum(1) = nx
        geom%Meshnum(2) = ny
        geom%Meshnum(3) = nz

        geom%Mshlocal(3) = geom%LcTabnm(1,myidx,myidy)
        geom%Mshlocal(2) = geom%LcTabnm(2,myidx,myidy)
        geom%Mshlocal(1) = geom%Meshnum(1)

        geom%Sptrnglocal(5) = geom%LcTabrg(1,myidx,myidy)
        geom%Sptrnglocal(6) = geom%LcTabrg(2,myidx,myidy)
        geom%Sptrnglocal(3) = geom%LcTabrg(3,myidx,myidy)
        geom%Sptrnglocal(4) = geom%LcTabrg(4,myidx,myidy)
        geom%Sptrnglocal(1) = geom%SpatRange(1)
        geom%Sptrnglocal(2) = geom%SpatRange(2)

        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine inpoint_Output

        subroutine phaseout_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: my_rank,ierr
        integer :: i,j,k,l,m,n,ioerr
        character*8 name1
        character*9 name2
        character*10 name3
        character*11 name4
        character*3 name
        integer gethostname
        real*4, dimension(:,:), allocatable :: ptout

        name1 = 'ph0000sx'
        name2 = 'ph0000sxx'
        name3 = 'ph0000sxxx'
        name4 = 'ph0000sxxxx'
!        ioerr = gethostname(name,3)
!        name = 'n01'

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
 
        if(nfile < 10) then
          if(my_rank.lt.10) then
            name1(6:6) = char(nfile+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            name2(6:6) = char(nfile+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            name3(6:6) = char(nfile+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            name4(6:6) = char(nfile+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.10).and.(nfile.lt.100)) then
          if(my_rank.lt.10) then
            i = nfile/10
            j = nfile - 10*i
            name1(5:5) = char(i+48)
            name1(6:6) = char(j+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(5:5) = char(i+48)
            name2(6:6) = char(j+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/10
            j = nfile - 10*i
            name3(5:5) = char(i+48)
            name3(6:6) = char(j+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/10
            j = nfile - 10*i
            name4(5:5) = char(i+48)
            name4(6:6) = char(j+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
          if(my_rank.lt.10) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name1(4:4) = char(i+48)
            name1(5:5) = char(k+48)
            name1(6:6) = char(l+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name2(4:4) = char(i+48)
            name2(5:5) = char(k+48)
            name2(6:6) = char(l+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(l+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else
          if(my_rank.lt.10) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name1(3:3) = char(i+48)
            name1(4:4) = char(k+48)
            name1(5:5) = char(m+48)
            name1(6:6) = char(n+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name2(3:3) = char(i+48)
            name2(4:4) = char(k+48)
            name2(5:5) = char(m+48)
            name2(6:6) = char(n+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name3(3:3) = char(i+48)
            name3(4:4) = char(k+48)
            name3(5:5) = char(m+48)
            name3(6:6) = char(n+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name4(3:3) = char(i+48)
            name4(4:4) = char(k+48)
            name4(5:5) = char(m+48)
            name4(6:6) = char(n+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        endif

!        do i = 1, this%Nptlocal
!          write(9,100)this%Pts1(1:9,i)
!          write(9)this%Pts1(1:6,i)
!        enddo
!         write(9)this%refptcl(1:6)

        allocate(ptout(9,this%Nptlocal))
        do i = 1, this%Nptlocal
           ptout(:,i) = this%Pts1(:,i)
        enddo
        write(9)this%Nptlocal
        write(9)ptout(1:9,1:this%Nptlocal)

        close(9)
100     format(9(1x,e14.7))

        deallocate(ptout)

        end subroutine phaseout_Output

        subroutine outpoint_Output(nfile,this,z,inb,jstp,ibal,nstep,&
                  nprocrow,nproccol,geom)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: inb,jstp,ibal,nstep,nprocrow,nproccol
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)
  
        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

        open(nfile,status="unknown",form="unformatted")

        write(nfile)z
        write(nfile)inb,jstp,ibal,nstep
        write(nfile)msize(1:3)
        write(nfile)range(1:6)
        write(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        write(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        deallocate(Localnum)
        deallocate(Localrange)

        write(nfile)this%Nptlocal
        write(nfile)this%refptcl(1:6)
        write(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        end subroutine outpoint_Output

        subroutine phase_Output(nfile,this,samplePeriod)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          if(samplePeriod.ge.0) then
           do i = 1, this%Nptlocal,samplePeriod
            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i),&
                            this%Pts1(7,i),this%Pts1(8,i),this%Pts1(9,i)
           enddo
           do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/9,samplePeriod
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j),&
                              recvbuf(7,j),recvbuf(8,j),recvbuf(9,j)
            enddo
           enddo
          else
           write (*,*) 'dumping phase space Impact-T format'
           do i = 1, this%Nptlocal,abs(samplePeriod)
            write(nfile,101)this%Pts1(1,i)*Scxl,this%Pts1(2,i), &
                  this%Pts1(3,i)*Scxl,this%Pts1(4,i), &
                 -this%Pts1(5,i)*Scxl,-this%Pts1(6,i)-this%refptcl(6)
           enddo
           do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/9,abs(samplePeriod)
              write(nfile,101)recvbuf(1,j)*Scxl,recvbuf(2,j),&
                    recvbuf(3,j)*Scxl,recvbuf(4,j),&
                   -recvbuf(5,j)*Scxl,-recvbuf(6,j)-this%refptcl(6)
            enddo
           enddo
          endif
          close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(9(g0,1x))
101     format(6(g0,1x))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phase_Output

        !//output the 3D particle number density
        subroutine dens3d_Output(nstep,nfile,this,totnptcls,xmni,xmxi,ymni,&
          ymxi,zmni,zmxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: totnptcls
        double precision, intent(in) :: xmni,xmxi,ymni,ymxi,zmni,zmxi
        integer, parameter :: nxcell = 100
        integer, parameter :: nycell = 100
        integer, parameter :: nzcell = 100
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,jj,kk,icc
        integer :: ix,iy,iz,ic,nsend,req
        integer, dimension(nxcell,nycell,nzcell) :: count
        integer, dimension(4,nxcell*nycell*nzcell) :: sendbuf,recvbuf
        integer, allocatable,dimension(:) :: nptlist
        double precision :: hxx,hyy,hzz,invvol,tmp1,tmp2,tmp3,tmp4,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: xl,xt
        double precision  :: xmnin,xmxin,ymnin,ymxin,zmnin,zmxin,xmin,xmax,&
        ymin,ymax,zmin,zmax

        xl = Scxl
        xt = 90.0/asin(1.0)

        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl
        zmnin = zmni/xt
        zmxin = zmxi/xt

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        eps = 1.0e-4
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + eps*temp1(i)
        enddo 
        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + eps*temp1(i)
        enddo 

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hzz = (zmax-zmin)/nzcell
        count = 0 
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          count(ix,iy,iz) = count(ix,iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
                sendbuf(4,ic) = count(i,j,k)
              enddo
            enddo
          enddo
        else
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
!                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
!                sendbuf(4,ic) = 0
                if(count(i,j,k).gt.0) then
                  icc = icc + 1
                  sendbuf(1,icc) = i
                  sendbuf(2,icc) = j
                  sendbuf(3,icc) = k
                  sendbuf(4,icc) = count(i,j,k)
                endif
              enddo
            enddo
          enddo
        endif
        nsend = 4*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/4
              ii = recvbuf(1,j)
              jj = recvbuf(2,j)
              kk = recvbuf(3,j)
              icc = (kk-1)*nycell*nxcell+(jj-1)*nxcell+ii
              sendbuf(4,icc) = sendbuf(4,icc) + recvbuf(4,j)
            enddo
          enddo

          open(nfile,status='unknown')
          tmp4 = invvol/dble(totnptcls)
          do k = 1, nzcell
            do j = 1, nycell
              do i = 1, nxcell
                tmp1 = (range(1) + i*hxx - 0.5*hxx)*xl
                tmp2 = (range(3) + j*hyy - 0.5*hyy)*xl
                tmp3 = (range(5) + k*hzz - 0.5*hzz)*xt
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell+i
                write(nfile,100)tmp1,tmp2,tmp3,sendbuf(4,ic)*tmp4
              enddo
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(4(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens3d_Output

        !//output 2D particle number density.
        subroutine dens2d_Output(nstep,nfile,this,totnptcls,xmnin,xmxin,&
        pxmnin,pxmxin,ymnin,ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: totnptcls
        double precision, intent(in) :: xmnin,xmxin,pxmnin,pxmxin,ymnin,&
        ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin
        integer, parameter :: nxcell = 128
        integer, parameter :: npxcell = 128
        integer, parameter :: nycell = 128
        integer, parameter :: npycell = 128
        integer, parameter :: nzcell = 128
        integer, parameter :: npzcell = 128
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,req,jj,kk,icc
        integer :: ix,ipx,iy,ipy,iz,ipz,ic,nsend,nfile1,nfile2
        integer :: nfile3,nfile4,nfile5
        integer, dimension(nxcell,npxcell) :: count1
        integer, dimension(nxcell,npxcell) :: count2
        integer, dimension(nxcell,npxcell) :: count3
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(nxcell,nzcell) :: count5
        integer, dimension(nycell,nzcell) :: count6
        integer, dimension(3,nxcell*npxcell) :: sendbuf1,recvbuf1
        integer, dimension(3,nycell*npycell) :: sendbuf2,recvbuf2
        integer, dimension(3,nzcell*npzcell) :: sendbuf3,recvbuf3
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, dimension(3,nxcell*nzcell) :: sendbuf5,recvbuf5
        integer, dimension(3,nycell*nzcell) :: sendbuf6,recvbuf6
        integer, allocatable, dimension(:) :: nptlist
        double precision, dimension(6) :: range,prange
        double precision, dimension(6) :: lcrange1,lcrange2,temp1
        double precision :: hxx,hyy,hzz,hpxx,hpyy,hpzz,invvol
        double precision ::eps,tmp1,tmp2,tmp3,xmax,xmin,pxmax,pxmin,&
        ymax,ymin,pymax,pymin,zmin,zmax,pzmin,pzmax
        character*3 name
        integer :: ioerr
        integer gethostname

!        ioerr = gethostname(name,3)
        eps = 1.0e-4
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange1(i+3) = this%Pts1(2*i,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
            lcrange2(i+3) = this%Pts1(2*i,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange1(i+3) = 5000.0
            lcrange2(i) = -5000.0
            lcrange2(i+3) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
            ii = i + 3
            if(lcrange1(ii).gt.this%Pts1(2*i,j)) then
              lcrange1(ii) = this%Pts1(2*i,j)
            endif
            if(lcrange2(ii).lt.this%Pts1(2*i,j)) then
              lcrange2(ii) = this%Pts1(2*i,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)

        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i-1) = temp1(i+3) + temp1(i+3)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i) = temp1(i+3) + temp1(i+3)*eps
        enddo 

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))
        pxmin = min(pxmnin,prange(1))
        pxmax = max(pxmxin,prange(2))
        pymin = min(pymnin,prange(3))
        pymax = max(pymxin,prange(4))
        pzmin = min(pzmnin,prange(5))
        pzmax = max(pzmxin,prange(6))

        hxx = (xmax-xmin)/nxcell
        hpxx = (pxmax-pxmin)/npxcell
        hyy = (ymax-ymin)/nycell
        hpyy = (pymax-pymin)/npycell
        hzz = (zmax-zmin)/nzcell
        hpzz = (pzmax-pzmin)/npzcell

        count1 = 0 
        count2 = 0 
        count3 = 0 
        count4 = 0
        count5 = 0
        count6 = 0
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          ipx = int((this%Pts1(2,i)-prange(1))/hpxx) + 1
          count1(ix,ipx) = count1(ix,ipx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          ipy = int((this%Pts1(4,i)-prange(3))/hpyy) + 1
          count2(iy,ipy) = count2(iy,ipy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          ipz = int((this%Pts1(6,i)-prange(5))/hpzz) + 1
          count3(iz,ipz) = count3(iz,ipz) + 1
          count4(ix,iy) = count4(ix,iy) + 1
          count5(ix,iz) = count5(ix,iz) + 1
          count6(iy,iz) = count6(iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hpxx)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf1(3,ic) = count1(i,j)
            enddo
          enddo
        else
          do j = 1, npxcell
            do i = 1, nxcell
              if(count1(i,j).gt.0) then
                icc = icc + 1
                sendbuf1(1,icc) = i
                sendbuf1(2,icc) = j
                sendbuf1(3,icc) = count1(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf1(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf1(1,j)
              jj = recvbuf1(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf1(3,icc) = sendbuf1(3,icc) + recvbuf1(3,j)
            enddo
          enddo

!          open(unit=nfile,file="/n/"//name//&
!               "/scratch/jiqiang/XPx.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile,file='XPx.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,prange(1)+0.5*hpxx,prange(1)+npxcell*hpxx-&
           0.5*hpxx
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = prange(1) + j*hpxx - 0.5*hpxx
!              write(nfile)sendbuf1(3,ic)
!              write(nfile)tmp1,tmp2,sendbuf1(3,ic)*tmp3
              write(nfile,100)tmp1,tmp2,sendbuf1(3,ic)*tmp3
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf1(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hpyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf2(3,ic) = count2(i,j) 
            enddo
          enddo
        else
          do j = 1, npycell
            do i = 1, nycell
              if(count2(i,j).gt.0) then
                icc = icc + 1
                sendbuf2(1,icc) = i
                sendbuf2(2,icc) = j
                sendbuf2(3,icc) = count2(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf2(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf2(1,j)
              jj = recvbuf2(2,j)
              icc = (jj-1)*nycell+ii
              sendbuf2(3,icc) = sendbuf2(3,icc) + recvbuf2(3,j)
            enddo
          enddo

          nfile1 = nfile + 1
!          open(unit=nfile1,file="/n/"//name//&
!               "/scratch/jiqiang/YPy.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile1,file='YPy.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile1,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,prange(3)+0.5*hpyy,prange(3)+npycell*hpyy-&
           0.5*hpyy
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = prange(3) + j*hpyy - 0.5*hpyy
!              write(nfile1)sendbuf2(3,ic)
!              write(nfile1)tmp1,tmp2,sendbuf2(3,ic)*tmp3
              write(nfile1,100)tmp1,tmp2,sendbuf2(3,ic)*tmp3
            enddo
          enddo
          close(nfile1)
        else
          call MPI_ISEND(sendbuf2(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hzz*hpzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              sendbuf3(3,ic) = count3(i,j) 
            enddo
          enddo
        else
          do j = 1, npzcell
            do i = 1, nzcell
              if(count3(i,j).gt.0) then
                icc = icc + 1
                sendbuf3(1,icc) = i
                sendbuf3(2,icc) = j
                sendbuf3(3,icc) = count3(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf3(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf3(1,j)
              jj = recvbuf3(2,j)
              icc = (jj-1)*nzcell + ii
              sendbuf3(3,icc) = sendbuf3(3,icc) + recvbuf3(3,j)
            enddo
          enddo

          nfile2 = nfile + 2
!          open(unit=nfile2,file="/n/"//name//&
!               "/scratch/jiqiang/ZPz.data",status="unknown",&
!               position="append",form="unformatted")
         open(unit=nfile2,file='ZPz.data',status='unknown',&
              position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile2,*)nstep,range(5)+0.5*hzz,range(5)+nzcell*hzz &
           -0.5*hzz,prange(5)+0.5*hpzz,prange(5)+npzcell*hpzz-&
           0.5*hpzz
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              tmp1 = range(5) + i*hzz - 0.5*hzz
              tmp2 = prange(5) + j*hpzz - 0.5*hpzz
!              write(nfile2)sendbuf3(3,ic)
!              write(nfile2)tmp1,tmp2,sendbuf3(3,ic)*tmp3
              write(nfile2,100)tmp1,tmp2,sendbuf3(3,ic)*tmp3
            enddo
          enddo
          close(nfile2)
        else
          call MPI_ISEND(sendbuf3(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
!          open(unit=nfile3,file="/n/"//name//&
!                "/scratch/jiqiang/XY.data",status="unknown",&
!                position="append",form="unformatted")
          open(unit=nfile3,file='XY.data',status='unknown',&
                position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile3,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(3)+0.5*hyy,range(3)+nycell*hyy-&
           0.5*hyy
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(3) + j*hyy - 0.5*hyy
!              write(nfile3)sendbuf4(3,ic)
!              write(nfile3)tmp1,tmp2,sendbuf4(3,ic)*tmp3
              write(nfile3,100)tmp1,tmp2,sendbuf4(3,ic)*tmp3
            enddo
          enddo
          close(nfile3)
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf5(3,ic) = count5(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nxcell
              if(count5(i,j).gt.0) then
                icc = icc + 1
                sendbuf5(1,icc) = i
                sendbuf5(2,icc) = j
                sendbuf5(3,icc) = count5(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf5(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf5(1,j)
              jj = recvbuf5(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf5(3,icc) = sendbuf5(3,icc) + recvbuf5(3,j)
            enddo
          enddo

          nfile4 = nfile + 4
!          open(unit=nfile4,file="/n/"//name//&
!               "/scratch/jiqiang/XZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile4,file='XZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile4,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile4)sendbuf5(3,ic)
!              write(nfile4)tmp1,tmp2,sendbuf5(3,ic)*tmp3
              write(nfile4,100)tmp1,tmp2,sendbuf5(3,ic)*tmp3
            enddo
          enddo
          close(nfile4)
        else
          call MPI_ISEND(sendbuf5(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf6(3,ic) = count6(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nycell
              if(count6(i,j).gt.0) then
                icc = icc + 1
                sendbuf6(1,icc) = i
                sendbuf6(2,icc) = j
                sendbuf6(3,icc) = count6(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf6(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf6(1,j)
              jj = recvbuf6(2,j)
              icc = (jj-1)*nycell + ii
              sendbuf6(3,icc) = sendbuf6(3,icc) + recvbuf6(3,j)
            enddo
          enddo

          nfile5 = nfile + 5
!          open(unit=nfile5,file="/n/"//name//&
!               "/scratch/jiqiang/YZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile5,file='YZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile5,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile5)sendbuf6(3,ic)
!              write(nfile5)tmp1,tmp2,sendbuf6(3,ic)*tmp3
              write(nfile5,100)tmp1,tmp2,sendbuf6(3,ic)*tmp3
            enddo
          enddo
          close(nfile5)
        else
          call MPI_ISEND(sendbuf6(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        deallocate(nptlist)

100     format(3(1x,e14.7))

        end subroutine dens2d_Output

        !//output the 1D accumulated particle density
        subroutine accdens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                 ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nptot,nstep
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 128
        integer, parameter :: nxcell = 128
        integer, parameter :: nycell = 128
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer, dimension(nxcell) :: densx,acdensx
        integer, dimension(nycell) :: densy,acdensy
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin,eps

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl

        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))
        call MPI_BARRIER(MPI_COMM_WORLD,ierr);
        !print*,"xmin: ",xmin,xmax,ymin,ymax,rmax

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0 
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nrcell
            glcount(i) = glcount(i) + glcount(i-1)
          enddo
          !print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens.data',status='unknown',&
               position='append')

          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)dble(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,dble(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)dble(sumy)/nptot
            !testsum = testsum + sumy
            densy(j) = sumy
          enddo
          acdensy = 0
          do j = nycell/2,1,-1
            acdensy(j) = densy(j) + densy(nycell-j+1) + acdensy(j+1)
          enddo
          do j = 1,nycell/2
            acdensy(j+nycell/2) = acdensy(nycell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nycell
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,dble(acdensy(j))/nptot
          enddo

          close(nfile3)
          !print*,"testsumy: ",testsum
        
          open(unit=nfile3,file='Xprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)dble(sumx)/nptot
            !testsum = testsum + sumx
            densx(i) = sumx
          enddo
          acdensx = 0
          do j = nxcell/2,1,-1
            acdensx(j) = densx(j) + densx(nxcell-j+1) + acdensx(j+1)
          enddo
          do j = 1,nxcell/2
            acdensx(j+nxcell/2) = acdensx(nxcell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nxcell
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,dble(acdensx(j))/nptot
          enddo
          close(nfile3)
          !print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine accdens1d_Output

        !//output 1d particle density
        subroutine dens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                  ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nptot,nstep
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 64
        integer, parameter :: nxcell = 64
        integer, parameter :: nycell = 64
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl
        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.0001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0 
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
!          do i = 2, nrcell
!            glcount(i) = glcount(i) + glcount(i-1)
!          enddo
!          print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens2.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens2.data',status='unknown',&
               position='append')
          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)dble(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,dble(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,dble(sumy)/nptot
            testsum = testsum + sumy
          enddo
          close(nfile3)
!          print*,"testsumy: ",testsum
        
          open(unit=nfile3,file='Xprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,dble(sumx)/nptot
            testsum = testsum + sumx
          enddo
          close(nfile3)
!          print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens1d_Output

        ! Terminate MPI
        subroutine end_Output(time)
        implicit none
        include 'mpif.h'
        double precision, intent(inout) :: time
        double precision :: endtime, mtime
        integer :: my_rank,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          close(24)
          close(25)
          print*,"time: ",mtime
        endif

        !for measurement of memory
        !call system_stats()

        call MPI_Finalize(ierr)

        end subroutine end_Output

!This program calculate the current profile, slice emittance
! and uncorrelated energy spread using linear deposition.
!J. Q. 04/26/12
      subroutine sliceprocdep_Output(pts,innp,npt,nslice,qchg,pmass,nfile,&
                                     alphax0,betax0,alphay0,betay0,gam0)
      implicit none
      include 'mpif.h'
      integer :: nfile,innp,nslice
      integer :: npt
      real*8, pointer, dimension(:,:) :: pts
      real*8, dimension(nslice) :: epx,epy,gam,xx,xx2,&
      px,px2,yy,yy2,py,py2,count,xpx,ypy,gam2uncor,gam2uncor2,&
      bmatchx,bmatchy
      real*8, dimension(nslice,12) :: tmparry,tmpglb
      real*8 :: qchg,pmass,alphax0,betax0,alphay0,betay0,gam0
      integer :: i,iz,iz1,nsend,ierr,my_rank
      real*8 :: zmin,zmax,hz,zz,zavg,detaabb,&
                deltagamma,gammaz,ab,clite,sclcur,zmingl,zmaxgl
      real*8 :: alphax,betax,gamx,alphay,betay,gamy,gamx0,gamy0,gambet

      clite = 2.99792d8
      sclcur = clite*qchg/npt

      gamx0 = (1.0d0+alphax0**2)/betax0
      gamy0 = (1.0d0+alphay0**2)/betay0

      zmin = 1.e12
      zmax = -1.e12
      do i = 1, innp
        if(zmin.ge.pts(5,i)) then
          zmin = pts(5,i)
        endif
        if(zmax.le.pts(5,i)) then
          zmax = pts(5,i)
        endif
      enddo
      call MPI_ALLREDUCE(zmin,zmingl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zmax,zmaxgl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MAX,MPI_COMM_WORLD,ierr)

      !avoid index overflow
      zmin = zmingl - (zmaxgl-zmingl)*0.5e-5
      zmax = zmaxgl + (zmaxgl-zmingl)*0.5e-5
      hz = (zmax-zmin)/(nslice-1)
      !print*,"zmin, zmax and hz: ",zmin*scxlt,zmax*scxlt,hz*scxlt
      count = 0.0
      xx = 0.0
      xx2 = 0.0
      px = 0.0
      px2 = 0.0
      xpx = 0.0
      yy = 0.0
      yy2 = 0.0
      py = 0.0 
      py2 = 0.0 
      ypy = 0.0
      gam = 0.0
      zavg = 0.0
      zavg = (zmin+zmax)/2
      !calculate the mean of each slice by linear deposition
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        count(iz) = count(iz) + ab
        count(iz1) = count(iz1) + 1.0d0 - ab
        xx(iz) = xx(iz) + pts(1,i)*ab
        xx(iz1) = xx(iz1) + pts(1,i)*(1.0d0-ab) 
        xx2(iz) = xx2(iz) + pts(1,i)**2*ab
        xx2(iz1) = xx2(iz1) + pts(1,i)**2*(1.0d0-ab)
        px(iz) = px(iz) + pts(2,i)*ab
        px(iz1) = px(iz1) + pts(2,i)*(1.0d0-ab)
        px2(iz) = px2(iz) + pts(2,i)**2*ab
        px2(iz1) = px2(iz1) + pts(2,i)**2*(1.0d0-ab)
        xpx(iz) = xpx(iz) + pts(1,i)*pts(2,i)*ab
        xpx(iz1) = xpx(iz1) + pts(1,i)*pts(2,i)*(1.0d0-ab)
        yy(iz) = yy(iz) + pts(3,i)*ab
        yy(iz1) = yy(iz1) + pts(3,i)*(1.0d0-ab)
        yy2(iz) = yy2(iz) + pts(3,i)**2*ab
        yy2(iz1) = yy2(iz1) + pts(3,i)**2*(1.0d0-ab)
        py(iz) = py(iz) + pts(4,i)*ab
        py(iz1) = py(iz1) + pts(4,i)*(1.0d0-ab)
        py2(iz) = py2(iz) + pts(4,i)**2*ab
        py2(iz1) = py2(iz1) + pts(4,i)**2*(1.0d0-ab)
        ypy(iz) = ypy(iz) + pts(3,i)*pts(4,i)*ab
        ypy(iz1) = ypy(iz1) + pts(3,i)*pts(4,i)*(1.0d0-ab)
        !gam(iz) = gam(iz) + sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1)*ab
        !gam(iz1) = gam(iz1) + &
        !           sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1)*(1.0-ab)
        gam(iz) = gam(iz) + pts(6,i)*ab
        gam(iz1) = gam(iz1) + pts(6,i)*(1.0-ab)
      enddo 
      tmparry(:,1) = count
      tmparry(:,2) = xx
      tmparry(:,3) = xx2
      tmparry(:,4) = px
      tmparry(:,5) = px2
      tmparry(:,6) = xpx
      tmparry(:,7) = yy
      tmparry(:,8) = yy2
      tmparry(:,9) = py
      tmparry(:,10) = py2
      tmparry(:,11) = ypy
      tmparry(:,12) = gam

      nsend = nslice*12

      tmpglb = 0.0d0
      call MPI_ALLREDUCE(tmparry,tmpglb,nsend,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)
      count = tmpglb(:,1)
      xx = tmpglb(:,2)
      xx2 = tmpglb(:,3)
      px = tmpglb(:,4)
      px2 = tmpglb(:,5)
      xpx = tmpglb(:,6)
      yy = tmpglb(:,7)
      yy2 = tmpglb(:,8)
      py = tmpglb(:,9)
      py2 = tmpglb(:,10)
      ypy = tmpglb(:,11)
      gam = tmpglb(:,12)

      !print*,"sum count: ",sum(count)
      do i = 1, nslice
        if(count(i).gt.0) then
          xx(i) = xx(i)/count(i)
          xx2(i) = xx2(i)/count(i)
          px(i) = px(i)/count(i)
          px2(i) = px2(i)/count(i)
          xpx(i) = xpx(i)/count(i)
          yy(i) = yy(i)/count(i)
          yy2(i) = yy2(i)/count(i)
          py(i) = py(i)/count(i)
          py2(i) = py2(i)/count(i)
          ypy(i) = ypy(i)/count(i)
          gam(i) = gam(i)/count(i)
        else
          xx(i) = 0.0d0
          xx2(i) = 0.0d0
          px(i) = 0.0d0
          px2(i) = 0.0d0
          yy(i) = 0.0d0
          yy2(i) = 0.0d0
          py(i) = 0.0d0
          py2(i) = 0.0d0
          xpx(i) = 0.0d0
          ypy(i) = 0.0d0
          gam(i) = 0.0d0
        endif
      enddo

      do i = 1, nslice
        if(count(i).gt.0) then
          epx(i) = sqrt((xx2(i)-xx(i)**2)*(px2(i)-px(i)**2)-&
                      (xpx(i)-xx(i)*px(i))**2)
          epy(i) = sqrt((yy2(i)-yy(i)**2)*(py2(i)-py(i)**2)-&
                      (ypy(i)-yy(i)*py(i))**2)
          alphax = -xpx(i)/epx(i)
          gambet = sqrt((gam0-gam(i))**2-1.0d0)
          betax = (xx2(i)-xx(i)**2)*gambet/epx(i)
          gamx = (1.0d0+alphax**2)/betax
          bmatchx(i) = (betax*gamx0-2*alphax*alphax0+gamx*betax0)
          alphay = -ypy(i)/epy(i)
          betay = (yy2(i)-yy(i)**2)*gambet/epy(i)
          gamy = (1.0d0+alphay**2)/betay
          bmatchy(i) = (betay*gamy0-2*alphay*alphay0+gamy*betay0)
        else
          epx(i) = 0.0d0
          epy(i) = 0.0d0
          bmatchx(i) = 0.0d0
          bmatchy(i) = 0.0d0
        endif
      enddo

      gam2uncor = 0.0
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        gammaz = gam(iz)*ab+ gam(iz1)*(1.0d0-ab)
        !deltagamma = sqrt(pts(2,i)**2+pts(4,i)**2+pts(6,i)**2+1.0)-gammaz
        deltagamma = -(pts(6,i)-gammaz)
        !write(11,111)pts(5,i),deltagamma,gammaz
        gam2uncor(iz) = gam2uncor(iz) + deltagamma**2*ab
        gam2uncor(iz1) = gam2uncor(iz1) + deltagamma**2*(1.0d0-ab)
      enddo
!111   format(3(1x,e15.7))
      gam2uncor2 = 0.0d0
      call MPI_ALLREDUCE(gam2uncor,gam2uncor2,nslice,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

      do i = 1, nslice
        gam2uncor2(i) = sqrt(gam2uncor2(i)/count(i))
      enddo

      call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

      if(my_rank.eq.0) then
        do i = 1, nslice
          zz = zmin + (i-1)*hz 
          write(nfile,777)-zz*scxl,count(i),count(i)/(hz*scxl)*sclcur,epx(i)*scxl,&
                  epy(i)*scxl,-gam(i)/gam0,gam2uncor2(i)*pmass,xx(i)*scxl,yy(i)*scxl,&
                  bmatchx(i),bmatchy(i)
        enddo
      endif

      call flush(nfile)

777   format(11(1x,e15.7))

      end subroutine sliceprocdep_Output

      end module Outputclass
