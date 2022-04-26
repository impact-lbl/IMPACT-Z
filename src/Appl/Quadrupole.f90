!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Quadrupoleclass: Quadrupole beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL, 2016
! Description: This class defines the linear transfer map and field
!              for the quadrupole beam line elment.
! Comments:
!----------------------------------------------------------------
      module Quadrupoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 9
        type Quadrupole
          !Itype = 1
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : quad gradient
          !      (3) : file ID
          !      (4) : radius
          !      (5) : x misalignment error
          !      (6) : y misalignment error
          !      (7) : rotation error x
          !      (8) : rotation error y
          !      (9) : rotation error z
        end type Quadrupole
        interface getparam_Quadrupole
          module procedure getparam1_Quadrupole,  &
                          getparam2_Quadrupole,   &
                          getparam3_Quadrupole
        end interface
        interface setparam_Quadrupole
          module procedure setparam1_Quadrupole,  &
                          setparam2_Quadrupole, setparam3_Quadrupole
        end interface
      contains
        subroutine construct_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Quadrupole
   
        subroutine setparam1_Quadrupole(this,i,value)
        implicit none
        type (Quadrupole), intent(out) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Quadrupole

        subroutine setparam2_Quadrupole(this,values)
        implicit none
        type (Quadrupole), intent(out) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param(1:Nparam) = values(1:Nparam)

        end subroutine setparam2_Quadrupole

        subroutine setparam3_Quadrupole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Quadrupole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Quadrupole
   
        subroutine getparam1_Quadrupole(this,i,blparam) 
        implicit none 
        type (Quadrupole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Quadrupole
  
        subroutine getparam2_Quadrupole(this,blparams)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams(1:Nparam) = this%Param(1:Nparam)

        end subroutine getparam2_Quadrupole

        subroutine getparam3_Quadrupole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Quadrupole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Quadrupole
       
        subroutine maplinear_Quadrupole(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Quadrupole), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,tfin
        integer  :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
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

        dlti=0.0
        thli=0.0

        call rk6i_Quadrupole(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        tfin = t + tau

        dltf = 0.0
        thlf = 0.0

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
        !print*,"xm: ",xm(1,1),xm(2,1),xm(1,2),xm(2,2),xm(3,3),xm(4,3),xm(4,4)

        end subroutine maplinear_Quadrupole

        subroutine rk6i_Quadrupole(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Quadrupole), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Quadrupole(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Quadrupole(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Quadrupole

        subroutine intfunc1_Quadrupole(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Quadrupole), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,zz

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = t - zedge
          call getfldfrg_Quadrupole(zz,this,bgrad)
        else
          bgrad = this%Param(2)
        endif
        qmcc = Bchg/Bmass

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0
        gbet=sqrt((gamma0-1.0)*(gamma0+1.0))

        f(5)=1.0/(beta0*Scxl)
        f(6)=0.0

        ! matrix elements
        brho=gbet/Clight/qmcc
        s11=bgrad/brho*Scxl
        s33=-bgrad/brho*Scxl
        s55=  0.0

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

        end subroutine intfunc1_Quadrupole

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Quadrupole(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp
        real*8 :: bgradp,bgradpp

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = pos(3)-zedge
          !call getfldfrg_Quadrupole(zz,this,bgrad)
          !call getfldfrgAna_Quadrupole(zz,this,bgrad)
          call getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
        else
          bgrad = this%Param(2)
        endif

!        dx = this%Param(5)
!        dy = this%Param(6)
!        anglex = this%Param(7)
!        angley = this%Param(8)
!        anglez = this%Param(9)

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

        end subroutine getflderr_Quadrupole
        
        !get external field without displacement and rotation errors.
        subroutine getfld_Quadrupole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Quadrupole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge,bgradp,bgradpp

        zedge = this%Param(1)
        zz=pos(3)-zedge
        if(this%Param(3).gt.1.0e-5) then
          !call getfldfrg_Quadrupole(zz,this,bgrad)
          !call getfldfrgAna_Quadrupole(zz,this,bgrad)
          call getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
          !bgradp = 0.0
          !bgradpp = 0.0
        else
          bgrad = this%Param(2)
        endif

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        if(this%Param(3).gt.0.0) then
          extfld(4) = bgrad*pos(2) - bgradpp*(pos(2)**3+3*pos(1)**2*pos(2))/12
          extfld(5) = bgrad*pos(1) - bgradpp*(pos(1)**3+3*pos(1)*pos(2)**2)/12
          extfld(6) = bgradp*pos(1)*pos(2)
        else
          extfld(4) = bgrad*pos(2)
          extfld(5) = bgrad*pos(1)
          extfld(6) = 0.0
        endif

        end subroutine getfld_Quadrupole

        !interpolate the field from the SC rf cavity onto bunch location.
        subroutine getfldfrg_Quadrupole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision:: hstep,slope
        integer :: klo,khi,k
        integer :: my_rank,ierr

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
        bgrad =edat(klo)+slope*(zz-zdat(klo))

        end subroutine getfldfrg_Quadrupole

        subroutine getfldfrgAna_Quadrupole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision :: bb,dd,z1,z2

        bb = this%Param(2)
        z1 = 2*(this%Length - this%Param(3))/2
        z2 = this%Length - z1
        dd = 2*this%Param(4) 
        if(zz.lt.z1) then
          bgrad = bb/(1.0+exp(-0.00004+4.518219*(-(zz-z1)/dd-1.5)))
        else if(zz.gt.z2) then
          bgrad = bb/(1.0+exp(-0.00004+4.518219*((zz-z2)/dd-1.5)))
        else 
          bgrad = bb
        endif

        !write(1,*)zz,bgrad

        end subroutine getfldfrgAna_Quadrupole

        subroutine getfldfrgAna2_Quadrupole(zz,this,bgrad,bgradp,bgradpp)
        implicit none
        include 'mpif.h'
        type (Quadrupole), intent(in) :: this
        double precision, intent(in) :: zz
        double precision, intent(out) :: bgrad
        double precision :: bgradp,bgradpp
        double precision :: bb,dd,z10,z20,z3,z4,tmpz
        double precision :: c1,c2,s1,s2,s3,s4
        double precision :: tmp1
 
        c1 = -0.00004d0
        c2 = 4.518219d0
        bb = this%Param(2)
        dd = 2*this%Param(4)
        !coordinate 0 points for Enge function.
        z10 = (this%Length - this%Param(3))/2 !for entrance
        z20 = this%Length - z10               !for exit
        !fringe field range inside coordinate 0 point.
        !this could be different for different Enge coefficients.
        z3 = z10 + 1.5*dd  !for entrance
        z4 = z20 - 1.5*dd  !for exit
        if(zz.lt.z3) then
          tmpz = -(zz - z10)/dd
          s1 = c1 + c2*tmpz
          bgrad = bb/(1.0+exp(s1))
          s2 = -c2/dd !first derivative
          tmp1 = exp(s1)*s2
          bgradp = -bb*tmp1/(1.0+exp(s1))**2
          s3 = 0.0 !second derivative
          s4 = exp(s1)*(s3+s2*s2)
          bgradpp = bb*(-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)
        else if(zz.gt.z4) then
          tmpz = (zz - z20)/dd
          s1 = c1 + c2*tmpz
          bgrad = bb/(1.0+exp(s1))
          s2 = c2/dd
          tmp1 = exp(s1)*s2
          bgradp = -bb*tmp1/(1.0+exp(s1))**2
          s3 = 0.0
          s4 = exp(s1)*(s3 + s2*s2)
          bgradpp = bb*(-s4/(1.0+exp(s1))**2+2*tmp1*tmp1/(1.0+exp(s1))**3)
        else
          bgrad = bb
          bgradp = 0.0
          bgradpp = 0.0
        endif

        !write(1,*)zz,bgrad,bgradp,bgradpp
 
        end subroutine getfldfrgAna2_Quadrupole

        !quad strength using K
        subroutine transfmapK_Quadrupole(tt,tau,this,refpt,Nplc,pts,qmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nplc
        double precision, intent(inout) :: tt
        double precision, intent(in) :: tau,qmass
        double precision, dimension(6), intent(inout) :: refpt
        double precision, dimension(6) :: tmp
        type (Quadrupole), intent(in) :: this
        double precision, pointer, dimension(:,:) :: pts
        real*8 :: xm11,xm12,xm21,xm22,xm33,xm34,xm43,xm44,gam,gambetz,&
                  betaz,beta0,rtkstrzz,rtkstr,kstr,gambet0
        integer  :: i
        real*8 :: t,cs,ss

        gambet0 = sqrt(refpt(6)**2-1.0d0)
        beta0 = sqrt(1.0d0-1.0d0/(refpt(6)**2))
        do i = 1, Nplc
          gam = -refpt(6) - pts(6,i)
          gambetz = sqrt(gam**2-1.0d0-pts(2,i)**2-pts(4,i)**2)
          betaz = gambetz/gam
          !each reference particle momentum
          gambetz = sqrt(gam**2-1.0d0) 
!          kstr = pts(7,i)*2.997928e8/gambetz*this%Param(2)
!Param(2) is the K defined in MAD, i.e. G/Brho
          kstr = pts(7,i)/qmass*gambet0/gambetz*this%Param(2)
          rtkstr = sqrt(abs(kstr))
          rtkstrzz = rtkstr*tau
          if(kstr.gt.0.0) then
            cs = cos(rtkstrzz)
            ss = sin(rtkstrzz)
            xm11 = cs
            xm12 = ss/rtkstr
            xm21 = -rtkstr*ss
            xm22 = cs
            !xm11 = cos(rtkstrzz)
            !xm12 = sin(rtkstrzz)/rtkstr
            !xm21 = -rtkstr*sin(rtkstrzz)
            !xm22 = cos(rtkstrzz)
            t = exp(rtkstrzz)
            cs= (t*t +1.0d0) / (t + t)
            ss= (t*t -1.0d0) / (t + t) 
            xm33 = cs
            xm34 = ss/rtkstr
            xm43 = ss*rtkstr
            xm44 = cs
            !xm33 = cosh(rtkstrzz)
            !xm34 = sinh(rtkstrzz)/rtkstr
            !xm43 = sinh(rtkstrzz)*rtkstr
            !xm44 = cosh(rtkstrzz)
          else if(kstr.lt.0.0) then
            t = exp(rtkstrzz)
            cs= (t*t +1.0d0) / (t + t)
            ss= (t*t -1.0d0) / (t + t) 
            xm11 = cs
            xm12 = ss/rtkstr
            xm21 = rtkstr*ss
            xm22 = cs
            !xm11 = cosh(rtkstrzz)
            !xm12 = sinh(rtkstrzz)/rtkstr
            !xm21 = rtkstr*sinh(rtkstrzz)
            !xm22 = cosh(rtkstrzz)
            cs = cos(rtkstrzz)
            ss = sin(rtkstrzz)
            xm33 = cs
            xm34 = ss/rtkstr
            xm43 = -ss*rtkstr
            xm44 = cs
            !xm33 = cos(rtkstrzz)
            !xm34 = sin(rtkstrzz)/rtkstr
            !xm43 = -sin(rtkstrzz)*rtkstr
            !xm44 = cos(rtkstrzz)
          else
            xm11 = 1.0d0
            xm12 = tau
            xm21 = 0.0d0
            xm22 = 1.0d0
            xm33 = 1.0d0
            xm34 = tau
            xm43 = 0.0
            xm44 = 1.0d0
          endif
          tmp(1) = xm11*pts(1,i)+xm12*pts(2,i)/gambetz/Scxl
          tmp(2) = gambetz*Scxl*xm21*pts(1,i)+xm22*pts(2,i)
          tmp(3) = xm33*pts(3,i)+xm34*pts(4,i)/gambetz/Scxl
          tmp(4) = gambetz*Scxl*xm43*pts(3,i)+xm44*pts(4,i)
          tmp(5) = pts(5,i) + (1.0/betaz-1.0/beta0)*tau/Scxl 
          tmp(6) = pts(6,i)
          pts(1,i) = tmp(1)
          pts(2,i) = tmp(2)
          pts(3,i) = tmp(3)
          pts(4,i) = tmp(4)
          pts(5,i) = tmp(5)
          pts(6,i) = tmp(6)
        enddo
        refpt(5) = refpt(5) + tau/beta0/Scxl
!        tt = tt + tau

        end subroutine transfmapK_Quadrupole

      end module Quadrupoleclass
