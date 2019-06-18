!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Multipoleclass: Multipole beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and field
!              for the multipole (sextupole, octupole, decapole)
!              beam line elment.
! Comments:
!----------------------------------------------------------------
      module Multipoleclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 10
        type Multipole
          !Itype = 5
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : id for sextupole(2), octupole(3), decapole(4)
          !      (3) : field strength
          !      (4) : file ID
          !      (5) : radius
          !      (6) : x misalignment error
          !      (7) : y misalignment error
          !      (8) : rotation error x
          !      (9) : rotation error y
          !      (10) : rotation error z
        end type Multipole
        interface getparam_Multipole
          module procedure getparam1_Multipole,  &
                          getparam2_Multipole,   &
                          getparam3_Multipole
        end interface
        interface setparam_Multipole
          module procedure setparam1_Multipole,  &
                          setparam2_Multipole, setparam3_Multipole
        end interface
      contains
        subroutine construct_Multipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Multipole), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_Multipole
   
        subroutine setparam1_Multipole(this,i,value)
        implicit none
        type (Multipole), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_Multipole

        subroutine setparam2_Multipole(this,values)
        implicit none
        type (Multipole), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_Multipole

        subroutine setparam3_Multipole(this,numseg,nmpstp,type,blength)
        implicit none
        type (Multipole), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_Multipole
   
        subroutine getparam1_Multipole(this,i,blparam) 
        implicit none 
        type (Multipole), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Multipole
  
        subroutine getparam2_Multipole(this,blparams)
        implicit none
        type (Multipole), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Multipole

        subroutine getparam3_Multipole(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (Multipole), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Multipole
       
!Warning: the linear map does not work for this case.
!You have to use Lorentz integrator 
!----------------------------------------------------------------
        subroutine maplinear_Multipole(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (Multipole), intent(in) :: this
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

        call rk6i_Multipole(h,mpstp,t0,y,18,this,Bchg,Bmass)

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

        end subroutine maplinear_Multipole

        subroutine rk6i_Multipole(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (Multipole), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_Multipole(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_Multipole(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_Multipole

        subroutine intfunc1_Multipole(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (Multipole), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55
        double precision :: zedge,qmcc,brho,bgrad,zz

        zedge = this%Param(1)
        if(this%Param(3).gt.1.0e-5) then
          zz = t - zedge
          call getfldfrg_Multipole(zz,this,bgrad)
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

        end subroutine intfunc1_Multipole
!------------------------------------------------------------------------

        !get external field with displacement and rotation errors.
        subroutine  getflderr_Multipole(pos,extfld,this,dx,dy,anglex,&
                                         angley,anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Multipole), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        double precision, dimension(3) :: temp,tmp
        integer :: idpole

        zedge = this%Param(1)
        idpole = this%Param(2)+0.1
        if(this%Param(4).gt.1.0e-5) then
          zz = pos(3)-zedge
          call getfldfrg_Multipole(zz,this,bgrad)
        else
          bgrad = this%Param(3)
        endif

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
        if(idpole.eq.2) then !sextupole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*pos(2)*pos(1)
          extfld(5) = bgrad*(pos(1)**2-pos(2)**2)/2
! MaryLie and Transport notation of bgrad
!          extfld(4) = 2*bgrad*pos(2)*pos(1)
!          extfld(5) = bgrad*(pos(1)**2-pos(2)**2)
        else if(idpole.eq.3) then !octupole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*(3*pos(1)**2*pos(2)-pos(2)**3)/6
          extfld(5) = bgrad*(pos(1)**3-3*pos(1)*pos(2)**2)/6
! MaryLie and Transport notation of bgrad
!          extfld(4) = bgrad*(3*pos(1)**2*pos(2)-pos(2)**3)
!          extfld(5) = bgrad*(pos(1)**3-3*pos(1)*pos(2)**2)
        else if(idpole.eq.4) then !decapole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*(pos(1)**3*pos(2)-pos(1)*pos(2)**3)/6
          extfld(5) = bgrad*(pos(1)**4-6*pos(1)**2*pos(2)**2+pos(2)**4)/24
! MaryLie and Transport notation of bgrad
!          extfld(4) = 4*bgrad*(pos(1)**3*pos(2)-pos(1)*pos(2)**3)
!          extfld(5) = bgrad*(pos(1)**4-6*pos(1)**2*pos(2)**2+pos(2)**4)
        else
          print*,"type of multipole not implemented...."
          stop
        endif
!        if(idpole.eq.2) then !sextupole
!          extfld(4) = bgrad*tmp(2)*tmp(1)
!          extfld(5) = bgrad*(tmp(1)**2-tmp(2)**2)/2
!        else if(idpole.eq.3) then !octupole
!          extfld(4) = bgrad*(3*tmp(1)**2*tmp(2)-tmp(2)**3)/6
!          extfld(5) = bgrad*(tmp(1)**3-3*tmp(1)*tmp(2)**2)/6
!        else if(idpole.eq.4) then !decapole
!          extfld(4) = bgrad*(tmp(1)**3*tmp(2)-tmp(1)*tmp(2)**3)/6
!          extfld(5) = bgrad*(tmp(1)**4-6*tmp(1)**2*tmp(2)**2+tmp(2)**4)/24
!        else
!          print*,"type of multipole not implemented...."
!          stop
!        endif
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

        end subroutine getflderr_Multipole
        
        !get external field without displacement and rotation errors.
        subroutine  getfld_Multipole(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (Multipole), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,bgrad,zedge
        integer :: idpole

        zedge = this%Param(1)
        idpole = this%Param(2)+0.1
        zz=pos(3)-zedge
        if(this%Param(4).gt.0.0) then
          call getfldfrg_Multipole(zz,this,bgrad)
        else
          bgrad = this%Param(3)
        endif

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        if(idpole.eq.2) then !sextupole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*pos(2)*pos(1)
          extfld(5) = bgrad*(pos(1)**2-pos(2)**2)/2
! MaryLie and Transport notation of bgrad
!          extfld(4) = 2*bgrad*pos(2)*pos(1)
!          extfld(5) = bgrad*(pos(1)**2-pos(2)**2)
        else if(idpole.eq.3) then !octupole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*(3*pos(1)**2*pos(2)-pos(2)**3)/6
          extfld(5) = bgrad*(pos(1)**3-3*pos(1)*pos(2)**2)/6
! MaryLie and Transport notation of bgrad
!          extfld(4) = bgrad*(3*pos(1)**2*pos(2)-pos(2)**3)
!          extfld(5) = bgrad*(pos(1)**3-3*pos(1)*pos(2)**2)
        else if(idpole.eq.4) then !decapole
! MAD and Wiedemann notation of bgrad
          extfld(4) = bgrad*(pos(1)**3*pos(2)-pos(1)*pos(2)**3)/6
          extfld(5) = bgrad*(pos(1)**4-6*pos(1)**2*pos(2)**2+pos(2)**4)/24
! MaryLie and Transport notation of bgrad
!          extfld(4) = 4*bgrad*(pos(1)**3*pos(2)-pos(1)*pos(2)**3)
!          extfld(5) = bgrad*(pos(1)**4-6*pos(1)**2*pos(2)**2+pos(2)**4)
        else 
          print*,"type of multipole not implemented...."
          stop
        endif
        extfld(6) = 0.0

        end subroutine getfld_Multipole

        !interpolate the field from the SC rf cavity onto bunch location.
        subroutine getfldfrg_Multipole(zz,this,bgrad)
        implicit none
        include 'mpif.h'
        type (Multipole), intent(in) :: this
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

        end subroutine getfldfrg_Multipole

      end module Multipoleclass
