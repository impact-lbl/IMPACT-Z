!----------------------------------------------------------------
! (c) Copyright, 2011 by the Regents of the University of California.
! Version: 1.1
! Author: Ji Qiang, LBNL
! Description: find the Fourier coeficients of input discrete RF data.
!----------------------------------------------------------------

      program rfcoef
      integer, parameter :: Ndata = 50000
      integer, parameter :: Ncoef = 100
      double precision, dimension(Ndata) :: zdata,edata
      double precision, dimension(Ncoef) :: Fcoef,Fcoef2
      double precision :: emax

      emax = 0.0d0
      open(3,file="rfdata.in",status="old")
      n = 0
10    continue
        read(3,*,end=100)tmp1,tmp2
        n = n+1
        zdata(n) = tmp1
        edata(n) = tmp2
        if(emax.le.abs(tmp2)) then
          emax = abs(tmp2)
        endif
      goto 10
100   continue
      close(3)
      ndatareal = n

      print*,"Emax is: ",emax

      print*,"How many Fourier coeficients you want?"
      read(*,*)ncoefreal

      zlen = zdata(ndatareal)-zdata(1)
      zhalf = zlen/2.0
      zmid = (zdata(ndatareal)+zdata(1))/2
      h = zlen/(ndatareal-1)
      pi = 2*asin(1.0)
      print*,"The RF data number is: ",ndatareal,zlen,zmid,h

      do j = 1, ncoefreal
        zz = zdata(1) - zmid
        Fcoef(j) = (-0.5*edata(1)*cos((j-1)*2*pi*zz/zlen)*h)/zhalf
        Fcoef2(j) = (-0.5*edata(1)*sin((j-1)*2*pi*zz/zlen)*h)/zhalf
        zz = zdata(ndatareal) - zmid
        Fcoef(j) = Fcoef(j)-(0.5*edata(ndatareal)*cos((j-1)*2*pi*zz/zlen)*h)&
                            /zhalf
        Fcoef2(j) = Fcoef2(j)-(0.5*edata(ndatareal)*sin((j-1)*2*pi*zz/zlen)*h)&
                            /zhalf
      enddo

      do i = 1, ndatareal

        zz = (i-1)*h + zdata(1)
!        zz = zdata(i) - zmid
        klo=1
        khi=ndatareal
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdata(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdata(khi)-zdata(klo)
        slope=(edata(khi)-edata(klo))/hstep
        ez1 =edata(klo)+slope*(zz-zdata(klo))
 
        !zz = (i-1)*h - zmid
        !zz = (i-1)*h - half
        zz = zdata(1) + (i-1)*h - zmid
!        zz = zdata(i) - zmid
        do j = 1, ncoefreal
          !Fcoef(j) = Fcoef(j) + (edata(i)*cos((j-1)*2*pi*zz/zlen)*h)/zhalf
          !Fcoef2(j) = Fcoef2(j) + (edata(i)*sin((j-1)*2*pi*zz/zlen)*h)/zhalf
          Fcoef(j) = Fcoef(j) + (ez1*cos((j-1)*2*pi*zz/zlen)*h)/zhalf
          Fcoef2(j) = Fcoef2(j) + (ez1*sin((j-1)*2*pi*zz/zlen)*h)/zhalf
        enddo
      enddo

      open(7,file="rfcoef.out",status="unknown")
      do j = 1, ncoefreal
        write(7,*)j,Fcoef(j),Fcoef2(j)
      enddo
      close(7)

      emax = 1.0d0
      open(8,file="rfdatax",status="unknown")
      write(8,*)Fcoef(1)/emax
      do j = 2, ncoefreal
        write(8,*)Fcoef(j)/emax
        write(8,*)Fcoef2(j)/emax
      enddo
      close(8)

      open(8,file="rfdata.out",status="unknown")
      do i = 1, ndatareal
        zz = zdata(i) - zmid
        tmpsum = 0.5*Fcoef(1)
        tmpsump = 0.0
        tmpsumpp = 0.0
        do j = 2,ncoefreal
         tmpsum = tmpsum + Fcoef(j)*cos((j-1)*2*pi*zz/zlen) + &
                  Fcoef2(j)*sin((j-1)*2*pi*zz/zlen)
         tmpsump = tmpsump-(j-1)*2*pi*Fcoef(j)*sin((j-1)*2*pi*zz/zlen)/zlen +&
                  (j-1)*2*pi*Fcoef2(j)*cos((j-1)*2*pi*zz/zlen)/zlen
         tmpsumpp = tmpsumpp-((j-1)*2*pi*Fcoef(j)/zlen)**2*cos((j-1)*2*pi*zz/zlen) -&
                  ((j-1)*2*pi*Fcoef2(j)/zlen)**2*sin((j-1)*2*pi*zz/zlen)
        enddo
        write(8,*)zdata(i),tmpsum,tmpsump,tmpsumpp
      enddo
      close(8)

      stop
      end program rfcoef

