!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! DTLclass: Drift-tube-linac beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and RF field
!              for the DTL beam line elment.
! Comments:
!----------------------------------------------------------------
      module DTLclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 25
        type DTL
          !Itype = 101
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : quad 1 length
          !      (8) : quad 1 gradient
          !      (9) : quad 2 length
          !      (10) : quad 2 gradient
          !      (11) : x misalignment error for Quad 1
          !      (12) : y misalignment error for Quad 1
          !      (13) : rotation error x for Quad 1
          !      (14) : rotation error y for Quad 1
          !      (15) : rotation error z for Quad 1
          !      (16) : x misalignment error for Quad 2
          !      (17) : x misalignment error for Quad 2
          !      (18) : rotation error x for Quad 2
          !      (19) : rotation error y for Quad 2
          !      (20) : rotation error z for Quad 2
          !      (21) : x misalignment error for RF cavity
          !      (22) : y misalignment error for RF cavity
          !      (23) : rotation error x for RF cavity
          !      (24) : rotation error y for RF cavity
          !      (25) : rotation error z for RF cavity
        end type DTL
        interface getparam_DTL
          module procedure getparam1_DTL,  &
                          getparam2_DTL,   &
                          getparam3_DTL
        end interface
        interface setparam_DTL
          module procedure setparam1_DTL,  &
                          setparam2_DTL, setparam3_DTL
        end interface
      contains
        subroutine construct_DTL(this,numseg,nmpstp,type,blength)
        implicit none
        type (DTL), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_DTL
   
        subroutine setparam1_DTL(this,i,value)
        implicit none
        type (DTL), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_DTL

        subroutine setparam2_DTL(this,values)
        implicit none
        type (DTL), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_DTL

        subroutine setparam3_DTL(this,numseg,nmpstp,type,blength)
        implicit none
        type (DTL), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_DTL
   
        subroutine getparam1_DTL(this,i,blparam) 
        implicit none 
        type (DTL), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_DTL
  
        subroutine getparam2_DTL(this,blparams)
        implicit none
        type (DTL), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_DTL

        subroutine getparam3_DTL(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (DTL), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_DTL
       
        subroutine maplinear_DTL(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (DTL), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,ww,theta0,qmcc
        double precision :: ein,e1in,e2in,ef,e1f,e2f
        double precision :: sinthi,costhi,uprimi,qpwi,tfin,sinthf,costhf
        double precision :: uprimf,qpwf
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90
        qmcc = Bchg/Bmass
        mpstp = this%Mapstp

        y(1)=0.0
        y(2)=0.0
        y(3)=0.0
        y(4)=0.0
        y(5)=refpt(5)
        y(6)=refpt(6)
        y(7)=1.0
        y(8)=0.0
        y(9)=0.0
        y(10)=1.0
        y(11)=1.0
        y(12)=0.0
        y(13)=0.0
        y(14)=1.0
        y(15)=1.0
        y(16)=0.0
        y(17)=0.0
        y(18)=1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3

        call getaxfldE_DTL(t,this,ein,e1in,e2in)
        sinthi=sin(ww*refpt(5)+theta0)
        costhi=cos(ww*refpt(5)+theta0)
        uprimi=qmcc*ein/betai*costhi
        qpwi=0.5*qmcc*Scxl/(ui*ww)
        dlti=Scxl*(0.5*uprimi/ui-qpwi*e1in*sinthi)
        thli=1.5*Scxl*(uprimi/ui)

        call rk6i_DTL(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin = t + tau
        call getaxfldE_DTL(tfin,this,ef,e1f,e2f)
        sinthf=sin(ww*refpt(5)+theta0)
        costhf=cos(ww*refpt(5)+theta0)
        uprimf=qmcc*ef/betaf*costhf
        qpwf=0.5*qmcc*Scxl/(uf*ww)
        dltf=Scxl*(0.5*uprimf/uf-qpwf*e1f*sinthf)
        thlf=1.5*Scxl*(uprimf/uf)

        xm(1,1)= y(7)*squi/squf+y(9)*squi/squf*dlti
        xm(2,1)=(y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti
        xm(1,2)= y(9)/(squi*squf)
        xm(2,2)=(y(10)-y(9)*dltf)*squf/squi
        xm(3,3)= y(11)*squi/squf+y(13)*squi/squf*dlti
        xm(4,3)= &
        (y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti
        xm(3,4)= y(13)/(squi*squf)
        xm(4,4)=(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3

        !print*,"xm: ",xm(1,1),xm(2,1),xm(1,2),xm(2,2),xm(3,3),xm(4,3),&
        !        xm(4,4),xm(5,5),xm(5,6),xm(6,6)

        end subroutine maplinear_DTL

        subroutine rk6i_DTL(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (DTL), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_DTL(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_DTL(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_DTL

        subroutine intfunc1_DTL(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (DTL), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55,s11tmp
        double precision :: zedge,escale,ez1,ezp1,ezpp1,ww,theta0,qmcc,&
                            sinphi,cosphi,rfdsgn,brho,bgrad
        integer :: my_rank, ierr

        zedge = this%Param(1)
        call getaxfldE_DTL(t,this,ez1,ezp1,ezpp1)
        call getBgradfld_DTL(t,this,bgrad)
        ww = this%Param(3)/Scfreq 
        theta0 = this%Param(4)*asin(1.0)/90 
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))
        sinphi=sin(ww*y(5)+theta0)
        cosphi=cos(ww*y(5)+theta0)
        rfdsgn=ez1*cosphi
        f(5)=1.0/(beta0*Scxl)
        f(6)=-qmcc*rfdsgn

        ! matrix elements
        brho=gbet/Clight/qmcc
        s11tmp=0.5e0*(1.0+0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2 + &
          qmcc*0.5/beta0**3/gamma0**3*ez1*sinphi*ww/Scxl
        s11=(s11tmp+bgrad/brho)*Scxl
        s33=(s11tmp-bgrad/brho)*Scxl
        s55=  &
         -1.5e0*qmcc/beta0**2/gamma0*ezp1*cosphi &
         +(beta0**2+0.5)/beta0**3/gamma0*qmcc*ez1*sinphi*ww/Scxl &
         +1.5e0*(1.0-0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2
        s55=Scxl*s55

        f(7)=y(8)/Scxl
        f(8)=-s11*y(7)
        f(9)=y(10)/Scxl
        f(10)=-s11*y(9)
        f(11)=y(12)/Scxl
        f(12)=-s33*y(11)
        f(13)=y(14)/Scxl
        f(14)=-s33*y(13)
        f(15)=y(16)/Scxl
        f(16)=-s55*y(15)
        f(17)=y(18)/Scxl
        f(18)=-s55*y(17)

        end subroutine intfunc1_DTL

        !interpolate the field from the DTL rf cavity onto bunch location.
        subroutine getaxfldE_DTL(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (DTL), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale
        integer :: klo,khi,k
        integer :: my_rank,ierr

        zedge = this%Param(1)
        escale = this%Param(2)
        zz=z-zedge
        klo=1
        khi=Ndata
1       if(khi-klo.gt.1) then
          k=(khi+klo)/2
          if(zdat(k).gt.zz)then
             khi=k
          else
             klo=k
          endif
          goto 1
        endif
        hstep=zdat(khi)-zdat(klo)
        slope=(edat(khi)-edat(klo))/hstep
        ez1 =edat(klo)+slope*(zz-zdat(klo))
        slope=(epdat(khi)-epdat(klo))/hstep
        ezp1=epdat(klo)+slope*(zz-zdat(klo))
        ezpp1 = 0.0
        ez1 = ez1*escale
        ezp1 = ezp1*escale

        !print*,"ez1: ",ez1,escale,zedge,zz
!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_DTL

        subroutine getBgradfld_DTL(z,this,bgrad)
        implicit none
        include 'mpif.h'
        type (DTL), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: bgrad
        double precision:: zz,bgrad1,bgrad2,len1,len2,len,zedge

        len = this%Length
        zedge = this%Param(1)
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)
        zz=z-zedge
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        end subroutine getBgradfld_DTL

        !get external field with displacement and rotation errors.
        subroutine  getflderrold_DTL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision :: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        double precision :: dx1,dy1,anglex1,angley1,anglez1
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        !zz=pos(3)-zedge

        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)

!        dx = this%Param(11)
!        dy = this%Param(12)
!        anglex = this%Param(13)
!        angley = this%Param(14)
!        anglez = this%Param(15)
        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = pos(3) - zedge
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*tmp(2)
        extfld(5) = bgrad*tmp(1)
        extfld(6) = 0.0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        dx1 = this%Param(16)
        dy1 = this%Param(17)
        anglex1 = this%Param(18)
        angley1 = this%Param(19)
        anglez1 = this%Param(20)

        zz=pos(3)-zmid
        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = zz
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = tmp(1)**2+tmp(2)**2
        tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0

        !for E field
        tmp(1) = tmpex
        tmp(2) = tmpey*cos(anglex)-tmpez*sin(anglex)
        tmp(3) = tmpey*sin(anglex)+tmpez*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = extfld(1)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = extfld(2)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = extfld(3)+temp(3)

        !for B field
        tmp(1) = tmpbx
        tmp(2) = tmpby*cos(anglex)-tmpbz*sin(anglex)
        tmp(3) = tmpby*sin(anglex)+tmpbz*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = extfld(4)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = extfld(5)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = extfld(6)+temp(3)

        end subroutine getflderrold_DTL
        
        !get external field with different field and Quad displacement and rotation errors.
        !Here, there 2 displacement errors and rotation errors for Quad 
        subroutine  getflderr_DTL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision :: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        double precision :: dx1,dy1,anglex1,angley1,anglez1
        double precision :: dx0,dy0,anglex0,angley0,anglez0,zz0
        double precision, dimension(3) :: temp,tmp
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        zz0=pos(3)-zedge

        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)

        if(zz0.lt.len1) then
          bgrad = bgrad1
        else if(zz0.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        if(zz0.lt.len1) then
          bgrad = bgrad1
          temp(1) = pos(1) - dx
          temp(2) = pos(2) - dy
          tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
          tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
          tmp(3) = pos(3) - zedge
          temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
          tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
          zz = tmp(3)

          extfld(1) = 0.0
          extfld(2) = 0.0
          extfld(3) = 0.0
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0

          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
          tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(6) = temp(3)
        else if(zz0.lt.len2) then
          bgrad = 0.0
          extfld = 0.0
        else
          dx0 = this%Param(16)
          dy0 = this%Param(17)
          anglex0 = this%Param(18)
          angley0 = this%Param(19)
          anglez0 = this%Param(20)
          temp(1) = pos(1) - dx0
          temp(2) = pos(2) - dy0
          tmp(1) = temp(1)*cos(anglez0) + temp(2)*sin(anglez0)
          tmp(2) = -temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          tmp(3) = pos(3) - zedge
          temp(1) = tmp(1)*cos(angley0)+tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex0)+temp(3)*sin(anglex0)
          tmp(3) = -temp(2)*sin(anglex0)+temp(3)*cos(anglex0)
          zz = tmp(3)

          extfld(1) = 0.0
          extfld(2) = 0.0
          extfld(3) = 0.0
          extfld(4) = bgrad*tmp(2)
          extfld(5) = bgrad*tmp(1)
          extfld(6) = 0.0

          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex0)-extfld(6)*sin(anglex0)
          tmp(3) = extfld(5)*sin(anglex0)+extfld(6)*cos(anglex0)
          temp(1) = tmp(1)*cos(angley0)-tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          extfld(4) = temp(1)*cos(anglez0) - temp(2)*sin(anglez0)
          extfld(5) = temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          extfld(6) = temp(3)
        endif
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        dx1 = this%Param(21)
        dy1 = this%Param(22)
        anglex1 = this%Param(23)
        angley1 = this%Param(24)
        anglez1 = this%Param(25)

        zz=pos(3)-zmid
        temp(1) = pos(1) - dx
        temp(2) = pos(2) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = zz
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = tmp(1)**2+tmp(2)**2
        tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0

        !for E field
        tmp(1) = tmpex
        tmp(2) = tmpey*cos(anglex)-tmpez*sin(anglex)
        tmp(3) = tmpey*sin(anglex)+tmpez*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1) = extfld(1)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2) = extfld(2)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3) = extfld(3)+temp(3)

        !for B field
        tmp(1) = tmpbx
        tmp(2) = tmpby*cos(anglex)-tmpbz*sin(anglex)
        tmp(3) = tmpby*sin(anglex)+tmpbz*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = extfld(4)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = extfld(5)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = extfld(6)+temp(3)

        end subroutine getflderr_DTL

        !get external field without displacement and rotation errors.
        subroutine  getfld_DTL(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        integer :: i

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = bgrad*pos(2)
        extfld(5) = bgrad*pos(1)
        extfld(6) = 0.0
        !This search is based on the assumption that the data is uniformly
        !distributed along z.
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        zz=pos(3)-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
!        f1 = 0.0
!        f1p = 0.0
        r2 = pos(1)**2+pos(2)**2
        tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0
        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = extfld(4) + tmpbx
        extfld(5) = extfld(5) + tmpby
        extfld(6) = extfld(6) + tmpbz

        end subroutine getfld_DTL

        !get external field without displacement and rotation errors.
        subroutine  getfldpts_DTL(pos,extfld,this,nplc)
        implicit none
        include 'mpif.h'
        integer,intent(in) ::nplc
        double precision, dimension(4,nplc), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, dimension(6,nplc), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        integer :: i,ipt

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)


        !move into the tank local coordinate.
        zz=pos(3,1)-zedge
        if(zz.lt.len1) then
          bgrad = bgrad1
        else if(zz.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        do ipt =1,nplc

        extfld(1,ipt) = 0.0
        extfld(2,ipt) = 0.0
        extfld(3,ipt) = 0.0
        extfld(4,ipt) = bgrad*pos(2,ipt)
        extfld(5,ipt) = bgrad*pos(1,ipt)
        extfld(6,ipt) = 0.0

        end do

        !This search is based on the assumption that the data is uniformly
        !distributed along z.
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        zz=pos(3,1)-zmid
        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
!        f1 = 0.0
!        f1p = 0.0

        do ipt=1,nplc

        tt = pos(4,ipt)
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        r2 = pos(1,ipt)**2+pos(2,ipt)**2
        tmpex = -pos(1,ipt)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -pos(2,ipt)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = pos(2,ipt)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -pos(1,ipt)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0
        extfld(1,ipt) = extfld(1,ipt) + tmpex
        extfld(2,ipt) = extfld(2,ipt) + tmpey
        extfld(3,ipt) = extfld(3,ipt) + tmpez
        extfld(4,ipt) = extfld(4,ipt) + tmpbx
        extfld(5,ipt) = extfld(5,ipt) + tmpby
        extfld(6,ipt) = extfld(6,ipt) + tmpbz

        enddo

        end subroutine getfldpts_DTL

        !get external field with different field and Quad displacement and rotation errors.
        !Here, there 2 displacement errors and rotation errors for Quad 
        subroutine  getflderrpts_DTL(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez,nplc)
        implicit none
        include 'mpif.h'
        integer,intent(in)::nplc
        double precision, dimension(4,nplc), intent(in) :: pos
        type (DTL), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6,nplc), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,zz,xl
        double precision :: ez1,ezp1,ezpp1,ezppp,f1,f1p,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision :: bgrad,bgrad1,bgrad2,len1,len2,r2,zmid
        double precision :: dx1,dy1,anglex1,angley1,anglez1
        double precision :: dx0,dy0,anglex0,angley0,anglez0,zz0
        double precision, dimension(3) :: temp,tmp
        integer :: i,ipt

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        zmid = zedge + len/2
        ww = this%Param(3)/Scfreq 
        xl = Scxl/ww  !real frequency has to be used here
        theta0 = this%Param(4)*asin(1.0)/90
        len1 = this%Param(7)
        bgrad1 = this%Param(8)
        len2 = len - this%Param(9)
        bgrad2 = this%Param(10)

        !move into the tank local coordinate.
        zz0=pos(3,1)-zedge

        if(zz0.lt.len1) then
          bgrad = bgrad1
        else if(zz0.lt.len2) then
          bgrad = 0.0
        else
          bgrad = bgrad2
        endif

        do ipt=1,nplc

        if(zz0.lt.len1) then
          bgrad = bgrad1
          temp(1) = pos(1,ipt) - dx
          temp(2) = pos(2,ipt) - dy
          tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
          tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
          tmp(3) = pos(3,ipt) - zedge
          temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
          tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
          zz = tmp(3)

          extfld(1,ipt) = 0.0
          extfld(2,ipt) = 0.0
          extfld(3,ipt) = 0.0
          extfld(4,ipt) = bgrad*tmp(2)
          extfld(5,ipt) = bgrad*tmp(1)
          extfld(6,ipt) = 0.0

          tmp(1) = extfld(4,ipt)
          tmp(2) = extfld(5,ipt)*cos(anglex)-extfld(6,ipt)*sin(anglex)
          tmp(3) = extfld(5,ipt)*sin(anglex)+extfld(6,ipt)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(4,ipt) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(5,ipt) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(6,ipt) = temp(3)
        else if(zz0.lt.len2) then
          bgrad = 0.0
          extfld(:,ipt) = 0.0
        else
          dx0 = this%Param(16)
          dy0 = this%Param(17)
          anglex0 = this%Param(18)
          angley0 = this%Param(19)
          anglez0 = this%Param(20)
          temp(1) = pos(1,ipt) - dx0
          temp(2) = pos(2,ipt) - dy0
          tmp(1) = temp(1)*cos(anglez0) + temp(2)*sin(anglez0)
          tmp(2) = -temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          tmp(3) = pos(3,ipt) - zedge
          temp(1) = tmp(1)*cos(angley0)+tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = -tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          tmp(1) = temp(1)
          tmp(2) = temp(2)*cos(anglex0)+temp(3)*sin(anglex0)
          tmp(3) = -temp(2)*sin(anglex0)+temp(3)*cos(anglex0)
          zz = tmp(3)

          extfld(1,ipt) = 0.0
          extfld(2,ipt) = 0.0
          extfld(3,ipt) = 0.0
          extfld(4,ipt) = bgrad*tmp(2)
          extfld(5,ipt) = bgrad*tmp(1)
          extfld(6,ipt) = 0.0

          tmp(1) = extfld(4,ipt)
          tmp(2) = extfld(5,ipt)*cos(anglex0)-extfld(6,ipt)*sin(anglex0)
          tmp(3) = extfld(5,ipt)*sin(anglex0)+extfld(6,ipt)*cos(anglex0)
          temp(1) = tmp(1)*cos(angley0)-tmp(3)*sin(angley0)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley0)+tmp(3)*cos(angley0)
          extfld(4,ipt) = temp(1)*cos(anglez0) - temp(2)*sin(anglez0)
          extfld(5,ipt) = temp(1)*sin(anglez0) + temp(2)*cos(anglez0)
          extfld(6,ipt) = temp(3)
        endif
! get external RF field on axis from analytical function
!-----------------------------------------------------------------
        dx1 = this%Param(21)
        dy1 = this%Param(22)
        anglex1 = this%Param(23)
        angley1 = this%Param(24)
        anglez1 = this%Param(25)

        zz=pos(3,ipt)-zmid
        temp(1) = pos(1,ipt) - dx
        temp(2) = pos(2,ipt) - dy
        tmp(1) = temp(1)*cos(anglez) + temp(2)*sin(anglez)
        tmp(2) = -temp(1)*sin(anglez) + temp(2)*cos(anglez)
        tmp(3) = zz
        temp(1) = tmp(1)*cos(angley)+tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = -tmp(1)*sin(angley)+tmp(3)*cos(angley)
        tmp(1) = temp(1)
        tmp(2) = temp(2)*cos(anglex)+temp(3)*sin(anglex)
        tmp(3) = -temp(2)*sin(anglex)+temp(3)*cos(anglex)
        zz = tmp(3)

        ez1 = Fcoef(1)/2
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        do i = 2, (Ndata-1)/2 + 1
          ez1 = ez1 + Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len) + &
                Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len)
          ezp1 = ezp1 + (i-1)*2*pi/len*(-Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)+&
                Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
          ezpp1 = ezpp1+((i-1)*2*pi/len)**2*(-Fcoef(2*i-2)*cos((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*sin((i-1)*2*pi*zz/len))
          ezppp = ezppp+((i-1)*2*pi/len)**3*(Fcoef(2*i-2)*sin((i-1)*2*pi*zz/len)&
                  - Fcoef(2*i-1)*cos((i-1)*2*pi*zz/len))
        enddo
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=ezpp1*escale
        ezppp = ezppp*escale
!-----------------------------------------------------------------
        tt = pos(4,ipt) 
        tmpcos = cos(ww*tt+theta0)
        tmpsin = sin(ww*tt+theta0)
        f1 = -(ezpp1+ez1/xl/xl)/4
        f1p = -(ezppp+ezp1/xl/xl)/4
        r2 = tmp(1)**2+tmp(2)**2
        tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
        tmpez = (ez1+f1*r2)*tmpcos
        tmpbx = tmp(2)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpby = -tmp(1)/(xl*clite)*(ez1/2+f1*r2/4)*tmpsin
        tmpbz = 0.0

        !for E field
        tmp(1) = tmpex
        tmp(2) = tmpey*cos(anglex)-tmpez*sin(anglex)
        tmp(3) = tmpey*sin(anglex)+tmpez*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(1,ipt) = extfld(1,ipt)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(2,ipt) = extfld(2,ipt)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(3,ipt) = extfld(3,ipt)+temp(3)

        !for B field
        tmp(1) = tmpbx
        tmp(2) = tmpby*cos(anglex)-tmpbz*sin(anglex)
        tmp(3) = tmpby*sin(anglex)+tmpbz*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4,ipt) = extfld(4,ipt)+temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5,ipt) = extfld(5,ipt)+temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6,ipt) = extfld(6,ipt)+temp(3)

        enddo

        end subroutine getflderrpts_DTL
      end module DTLclass
