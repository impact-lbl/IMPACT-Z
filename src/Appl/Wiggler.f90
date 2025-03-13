!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Wigglerclass: Wiggler beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and field
!              for the wiggler/undulator beam line elment.
! Comments:
!----------------------------------------------------------------
      module Wigglerclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 12
        type Wiggler
          !Itype = 6
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : id for planar(1), helical(2)
          !      (3) : max. field strength
          !      (4) : file ID
          !      (5) : radius
          !      (6) : kx
          !      (7) : wiggler period
          !      (8) : x misalignment error
          !      (9) : y misalignment error
          !      (10) : rotation error x
          !      (11) : rotation error y
          !      (12) : rotation error z
        end type Wiggler
        interface getparam_Wiggler
          module procedure getparam1_Wiggler,  &
                          getparam2_Wiggler,   &
                          getparam3_Wiggler
        end interface
        interface setparam_Wiggler
          module procedure setparam1_Wiggler,  &
                          setparam2_Wiggler, setparam3_Wiggler
        end interface
      contains
        subroutine construct_Wiggler(this,numseg,nmpstp,type,blength)
        implicit none
        type (Wiggler), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Wiggler
   
        subroutine setparam1_Wiggler(this,i,value)
        implicit none
        type (Wiggler), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Wiggler

        subroutine setparam2_Wiggler(this,values)
        implicit none
        type (Wiggler), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param(1:Nparam) = values(1:Nparam)

        end subroutine setparam2_Wiggler

        subroutine setparam3_Wiggler(this,numseg,nmpstp,type,blength)
        implicit none
        type (Wiggler), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Wiggler
   
        subroutine getparam1_Wiggler(this,i,blparam) 
        implicit none 
        type (Wiggler), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Wiggler
  
        subroutine getparam2_Wiggler(this,blparams)
        implicit none
        type (Wiggler), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams(1:Nparam) = this%Param(1:Nparam)

        end subroutine getparam2_Wiggler

        subroutine getparam3_Wiggler(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Wiggler), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Wiggler
       
!Warning: the linear map does not work for this case.
!You have to use Lorentz integrator 
!----------------------------------------------------------------
        subroutine maplinear_Wiggler(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Wiggler), intent(in) :: this
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

        call rk6i_Wiggler(h,mpstp,t0,y,18,this,Bchg,Bmass)

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

        end subroutine maplinear_Wiggler

        subroutine rk6i_Wiggler(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Wiggler), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Wiggler(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Wiggler(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Wiggler

        subroutine intfunc1_Wiggler(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Wiggler), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,zz

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = t - zedge
          call getfldfrg_Wiggler(zz,this,bgrad)
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

        end subroutine intfunc1_Wiggler
!------------------------------------------------------------------------

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Wiggler(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Wiggler), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp
        integer :: idpole

        print*,"wrong subroutine in Wiggler getflderr:"
        return

        end subroutine getflderr_Wiggler
        
        !get external field without displacement and rotation errors.
        subroutine  getfld_Wiggler(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Wiggler), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,b0,zedge
        integer :: idpole
        real*8 :: xkx,xky,xku,leng,phiz

        leng = this%Length
        zedge = this%Param(1)
        !flag for planary or helical undulator
        idpole = this%Param(2)+0.1 
        zz=pos(3)-zedge
        if(this%Param(4).gt.0.0) then !should not be used now
          call getfldfrg_Wiggler(zz,this,b0)
        else
          b0 = this%Param(3)
        endif

        extfld = 0.0d0

        pi = 2*asin(1.0d0)
        xkx = this%Param(6) !kx
        xku = 2*pi/this%Param(7) !kz
        phiz = -xku*leng/2
        xky = sqrt(xkx**2+xku**2)
       
        if(idpole.eq.1) then !planar
          extfld(4) = -b0*xkx/xky*sin(xkx*pos(1))*sinh(xky*pos(2))*cos(xku*zz+phiz)
          extfld(5) = b0*cos(xkx*pos(1))*cosh(xky*pos(2))*cos(xku*zz+phiz)
          extfld(6) = -b0*xku/xky*cos(xkx*pos(1))*sinh(xky*pos(2))*sin(xku*zz+phiz)
        else if(idpole.eq.2) then !helical
          extfld(4) = -b0*cosh(xku*pos(1))*sin(xku*zz+phiz)
          extfld(5) = b0*cosh(xku*pos(2))*cos(xku*zz+phiz)
          extfld(6) = -b0*(sinh(xku*pos(1))*cos(xku*zz+phiz)+&
                           sinh(xku*pos(2))*sin(xku*zz+phiz))
        else
          print*,"wrong input in wiggler element!"
          return
        endif

        end subroutine getfld_Wiggler

        !interpolate the field from the SC rf cavity onto bunch location.
        subroutine getfldfrg_Wiggler(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Wiggler), intent(in) :: this
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

        end subroutine getfldfrg_Wiggler

      end module Wigglerclass
