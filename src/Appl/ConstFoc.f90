!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! ConstFocclass: 3D constant focusing beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and field
!              for the 3d constant focusing beam line elment.
! Comments:
!----------------------------------------------------------------
      module ConstFocclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 5
        type ConstFoc
          !Itype = 2
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : x focusing gradient: kx0^2
          !      (3) : y focusing gradient: ky0^2 
          !      (4) : z focusing gradient: kz0^2 
          !      (5) : radius 
        end type ConstFoc
        interface getparam_ConstFoc
          module procedure getparam1_ConstFoc,  &
                          getparam2_ConstFoc,   &
                          getparam3_ConstFoc
        end interface
        interface setparam_ConstFoc
          module procedure setparam1_ConstFoc,  &
                          setparam2_ConstFoc, setparam3_ConstFoc
        end interface
      contains
        subroutine construct_ConstFoc(this,numseg,nmpstp,type,blength)
        implicit none
        type (ConstFoc), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_ConstFoc
   
        subroutine setparam1_ConstFoc(this,i,value)
        implicit none
        type (ConstFoc), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_ConstFoc

        subroutine setparam2_ConstFoc(this,values)
        implicit none
        type (ConstFoc), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_ConstFoc

        subroutine setparam3_ConstFoc(this,numseg,nmpstp,type,blength)
        implicit none
        type (ConstFoc), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_ConstFoc
   
        subroutine getparam1_ConstFoc(this,i,blparam) 
        implicit none 
        type (ConstFoc), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_ConstFoc
  
        subroutine getparam2_ConstFoc(this,blparams)
        implicit none
        type (ConstFoc), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_ConstFoc

        subroutine getparam3_ConstFoc(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (ConstFoc), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_ConstFoc
       
        subroutine maplinear_ConstFoc(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (ConstFoc), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,tfin
        integer  :: mpstp

        xm = 0.0

        zedge = this%Param(1)
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

        call rk6i_ConstFoc(h,mpstp,t0,y,18,this,Bchg,Bmass)

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

        end subroutine maplinear_ConstFoc

        subroutine rk6i_ConstFoc(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (ConstFoc), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_ConstFoc(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_ConstFoc(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_ConstFoc

        subroutine intfunc1_ConstFoc(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (ConstFoc), intent(in) :: this
        double precision :: gamma0,beta0,s11,s33,s55
        double precision :: zedge,qmcc,xkperp021,xkperp022,xklong02

        xkperp021 = this%Param(2)
        xkperp022 = this%Param(3)
        xklong02 = this%Param(4)
!        xkperp021=(4.*asin(1.0))**2
!        xkperp022=(4.*asin(1.0))**2
!        xklong02=(4.*asin(1.0)*0.565)**2

        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        f(4) = 0.0

        ! synchronous particle:
        gamma0=-y(6)
        beta0=sqrt(gamma0**2-1.)/gamma0

        f(5)=1.0/(beta0*Scxl)
        f(6)=0.0

        ! matrix elements
        s11=xkperp021*Scxl
        s33=xkperp022*Scxl
        s55=xklong02*Scxl

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

        end subroutine intfunc1_ConstFoc

        subroutine  getfld_ConstFoc(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (ConstFoc), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zedge
        double precision :: xkperp02x,xkperp02y,xklong02,qmcc,gambetz,pz

        zedge = this%Param(1)
        xkperp02x = this%Param(2)
        xkperp02y = this%Param(3)
        xklong02 = this%Param(4)
        !needs to be set
        qmcc = 1.0
        gambetz = 1.0
        pz = gambetz 

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = xklong02*pos(4)*pz*pz*pz*Scxl/qmcc
        extfld(4) = -xkperp02x*pos(2)*pz/qmcc/Clight
        extfld(5) = xkperp02y*pos(1)*pz/qmcc/Clight
        extfld(6) = 0.0

        end subroutine getfld_ConstFoc
        
      end module ConstFocclass
