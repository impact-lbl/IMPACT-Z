!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Solclass: Solenoid beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Authors: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and field
!              for the Solenoid beam line elment.
! Comments: The linear map does NOT work.
!----------------------------------------------------------------
      module Solclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 9
        type Sol
          !Itype = 3
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : Bz0
          !      (3) : file ID
          !      (4) : radius
          !      (5) : x misalignment error
          !      (6) : y misalignment error
          !      (7) : rotation error x
          !      (8) : rotation error y
          !      (9) : rotation error z
        end type Sol
        interface getparam_Sol
          module procedure getparam1_Sol,  &
                          getparam2_Sol,   &
                          getparam3_Sol
        end interface
        interface setparam_Sol
          module procedure setparam1_Sol,  &
                          setparam2_Sol, setparam3_Sol
        end interface
      contains
        subroutine construct_Sol(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Sol
   
        subroutine setparam1_Sol(this,i,value)
        implicit none
        type (Sol), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Sol

        subroutine setparam2_Sol(this,values)
        implicit none
        type (Sol), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Sol

        subroutine setparam3_Sol(this,numseg,nmpstp,type,blength)
        implicit none
        type (Sol), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Sol
   
        subroutine getparam1_Sol(this,i,blparam) 
        implicit none 
        type (Sol), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Sol
  
        subroutine getparam2_Sol(this,blparams)
        implicit none
        type (Sol), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Sol

        subroutine getparam3_Sol(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Sol), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Sol
       
        subroutine maplinear_Sol(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Sol), intent(in) :: this
        double precision, dimension(22) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,qmcc
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
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
        y(19) = 1.0
        y(20) = 0.0
        y(21) = 0.0
        y(22) = 1.0
        h=tau/mpstp
        t0=t
        gi=-refpt(6)
        betai=sqrt(gi**2-1.)/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3
        dlti = 0.0
        thli = 0.0

        call rk6i_Sol(h,mpstp,t0,y,22,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3

        dltf = 0.0
        thlf = 0.0

        xm(1,1)= y(19)*(y(7)*squi/squf+y(9)*squi/squf*dlti)
        xm(2,1)=y(19)*((y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti)
        xm(1,2)= y(19)*y(9)/(squi*squf)
        xm(2,2)=y(19)*(y(10)-y(9)*dltf)*squf/squi
        xm(1,3)= y(20)*(y(7)*squi/squf+y(9)*squi/squf*dlti)
        xm(2,3)=y(20)*((y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*&
                 dlti)
        xm(1,4)= y(20)*y(9)/(squi*squf)
        xm(2,4)=y(20)*(y(10)-y(9)*dltf)*squf/squi
        xm(3,1)= y(21)*(y(11)*squi/squf+y(13)*squi/squf*dlti)
        xm(4,1)= &
        y(21)*((y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti)
        xm(3,2)= y(21)*y(13)/(squi*squf)
        xm(4,2)=y(21)*(y(14)-y(13)*dltf)*squf/squi
        xm(3,3)= y(22)*(y(11)*squi/squf+y(13)*squi/squf*dlti)
        xm(4,3)= y(22)* &
        ((y(12)-y(11)*dltf)*squi*squf+(y(14)-y(13)*dltf)*squi*squf*dlti)
        xm(3,4)= y(22)*y(13)/(squi*squf)
        xm(4,4)=y(22)*(y(14)-y(13)*dltf)*squf/squi
        xm(5,5)= y(15)*sqi3/sqf3+y(17)*sqi3/sqf3*thli
        xm(6,5)= &
        (y(16)-y(15)*thlf)*sqi3*sqf3+(y(18)-y(17)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(17)/(sqi3*sqf3)
        xm(6,6)=(y(18)-y(17)*thlf)*sqf3/sqi3

        end subroutine maplinear_Sol

        subroutine rk6i_Sol(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Sol), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Sol(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Sol(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Sol

        subroutine intfunc1_Sol(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Sol), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,C2alpha,b0
        integer :: my_rank, ierr

        zedge = this%Param(1)
        b0 = this%Param(2)
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

        s11=((0.5*b0/brho)**2)*Scxl
        s33=((0.5*b0/brho)**2)*Scxl
        s55= 0.0

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
        C2alpha = -b0/brho/2 !rotation matrix of Solenoid
        f(19)=-C2alpha*y(21)
        f(20)=-C2alpha*y(22)
        f(21)=C2alpha*y(19)
        f(22)=C2alpha*y(20)
 
        end subroutine intfunc1_Sol

        subroutine getBgradfld_Sol(z,this,b0,bgrad)
        implicit none
        include 'mpif.h'
        type (Sol), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: b0,bgrad

        !uniform bz field.
        b0 = this%Param(2)
        bgrad = 0.0

        end subroutine getBgradfld_Sol

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Sol(pos,extfld,this,dx,dy,anglex,angley,&
                                  anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,b0
        double precision, dimension(3) :: temp,tmp

        zedge = this%Param(1)
        b0 = this%Param(2)
        !move into the tank local coordinate.
!        zz=pos(3)-zedge
        !uniform bz field.

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

        bgrad = 0.0
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = -bgrad/2*tmp(1)
        extfld(5) = -bgrad/2*tmp(2)
        extfld(6) = b0

        tmp(1) = extfld(4)
        tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
        tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
        temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
        temp(2) = tmp(2)
        temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
        extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
        extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
        extfld(6) = temp(3)

        end subroutine getflderr_Sol
        
        !get external field without displacement and rotation errors and
        !with fringe field of Solenoid. (f(z) = b0 + bb*z)
        !here, the length of the solenoid has to the effective hard-edge length + 2*aperature 
        subroutine  getfld_Sol(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,b0,bb,cc,zshift,x0,x02

        zedge = this%Param(1)
        x0 = this%Param(4) 
        x02 = 2*x0 !range of fringe field (x0 outside the length)
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        !uniform bz field.
        b0 = this%Param(2)
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0

        if((zz.ge.0.0).and.(zz.lt.x02)) then
          bb = b0/x02
          bgrad = bb
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = bb*zz
        else if((zz.ge.x02).and.(zz.le.this%Length-x02)) then
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = b0
        else if((zz.gt.this%Length-x02).and.(zz.le.this%Length)) then
          bb = -b0/x02
          bgrad = bb
          zshift = this%Length - x02
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = b0 + bb*(zz-zshift)
        else
          !stop
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = 0.0
        endif

        end subroutine getfld_Sol

        !get external field without displacement and rotation errors and
        !with fringe field of Solenoid. (f(z) = b0 + bb*z^2 + cc*z^4,f(x0)=0,f'(x0)=0)
        subroutine  getfldold_Sol(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Sol), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge
        double precision:: zz,bgrad,b0,bb,cc,zshift,x0

        zedge = this%Param(1)
        x0 = 2*this%Param(4)
        !move into the tank local coordinate.
        zz=pos(3)-zedge
        !uniform bz field.
        b0 = this%Param(2)
        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0

        if((zz.ge.0.0).and.(zz.lt.x0)) then
          bb = -2*b0/(x0*x0)
          cc = b0/(x0*x0*x0*x0)
          bgrad = 2*bb*(zz-x0)+4*cc*(zz-x0)*(zz-x0)*(zz-x0)
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = b0 + bb*(zz-x0)*(zz-x0) + cc*(zz-x0)* &
                              (zz-x0)*(zz-x0)*(zz-x0)
        else if((zz.ge.x0).and.(zz.lt.this%Length-x0)) then
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = b0
        else if((zz.ge.this%Length-x0).and.(zz.le.this%Length)) then
          bb = -2*b0/(x0*x0)
          cc = b0/(x0*x0*x0*x0)
          zshift = this%Length - x0
          bgrad = 2*bb*(zz-zshift)+4*cc*(zz-zshift)*(zz-zshift)*(zz-zshift)
          extfld(4) = -bgrad/2*pos(1)
          extfld(5) = -bgrad/2*pos(2)
          extfld(6) = b0 + bb*(zz-zshift)*(zz-zshift) + cc*(zz-zshift)* &
                             (zz-zshift)*(zz-zshift)*(zz-zshift)
        else
          !stop
          extfld(4) = 0.0
          extfld(5) = 0.0
          extfld(6) = 0.0
        endif

        end subroutine getfldold_Sol
      end module Solclass
