!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! TWSclass: Traveling wave structure beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and RF field
!              for the TWS beam line elment. Here, the TWS does not include
!              input and output cells. The TWS is simulated using the
!              superposition of two standing wave structure.
!              (G. A. Loew, et al. SLAC-PUB-2295, 1979.)
! Comments:
!----------------------------------------------------------------
      module TWSclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 15
        type TWS
          !Itype = 106
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : RF frequency
          !      (4) : theta0
          !      (5) : file ID
          !      (6) : radius
          !      (7) : x misalignment error
          !      (8) : y misalignment error
          !      (9) : rotation error x
          !      (10) : rotation error y
          !      (11) : rotation error z
          !      (12) : theta1 (pi - beta * d) phase difference B and A
          !      (13) : aawk: aperture size for the wakefield 
          !      (14) : ggwk: gap size for the wakefield
          !      (15) : lengwk: length for the wakefield
        end type TWS
        interface getparam_TWS
          module procedure getparam1_TWS,  &
                          getparam2_TWS,   &
                          getparam3_TWS
        end interface
        interface setparam_TWS
          module procedure setparam1_TWS,  &
                          setparam2_TWS, setparam3_TWS
        end interface
      contains
        subroutine construct_TWS(this,numseg,nmpstp,type,blength)
        implicit none
        type (TWS), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_TWS
   
        subroutine setparam1_TWS(this,i,value)
        implicit none
        type (TWS), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_TWS

        subroutine setparam2_TWS(this,values)
        implicit none
        type (TWS), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_TWS

        subroutine setparam3_TWS(this,numseg,nmpstp,type,blength)
        implicit none
        type (TWS), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_TWS
   
        subroutine getparam1_TWS(this,i,blparam) 
        implicit none 
        type (TWS), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_TWS
  
        subroutine getparam2_TWS(this,blparams)
        implicit none
        type (TWS), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_TWS

        subroutine getparam3_TWS(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (TWS), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_TWS
       
        subroutine maplinear_TWS(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (TWS), intent(in) :: this
        double precision, dimension(18) :: y
        double precision :: h,t0,gi,betai,ui,squi,sqi3
        double precision :: dlti,thli,gf,betaf,uf,squf,sqf3
        double precision :: dltf,thlf,zedge,escale,ww,theta0,qmcc
        double precision :: ein,e1in,e2in,ef,e1f,e2f
        double precision :: ein2,e1in2,e2in2,ef2,e1f2,e2f2
        double precision :: sinthi,costhi,uprimi,qpwi,tfin,sinthf,costhf
        double precision :: sinthi2,costhi2,sinthf2,costhf2,theta2
        double precision :: uprimf,qpwf
        integer :: mpstp

        xm = 0.0

        zedge = this%Param(1)
        escale = this%Param(2)
        ww = this%Param(3)/Scfreq
        theta0 = this%Param(4)*asin(1.0)/90
        theta2 = this%Param(12)*asin(1.0)/90
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

        call getaxfldE_TWS(t,this,ein,e1in,e2in,ein2,e1in2,e2in2)
        sinthi=sin(ww*refpt(5)+theta0)
        costhi=cos(ww*refpt(5)+theta0)
        sinthi2=sin(ww*refpt(5)+theta0+theta2)
        costhi2=cos(ww*refpt(5)+theta0+theta2)
        uprimi=qmcc*(ein*costhi+ein2*costhi2)/betai
        qpwi=0.5*qmcc*Scxl/(ui*ww)
        dlti=Scxl*(0.5*uprimi/ui-qpwi*(e1in*sinthi+e1in2*sinthi2))
        thli=1.5*Scxl*(uprimi/ui)

        call rk6i_TWS(h,mpstp,t0,y,18,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin = t + tau
        call getaxfldE_TWS(tfin,this,ef,e1f,e2f,ef2,e1f2,e2f2)
        sinthf=sin(ww*refpt(5)+theta0)
        costhf=cos(ww*refpt(5)+theta0)
        sinthf2=sin(ww*refpt(5)+theta0+theta2)
        costhf2=cos(ww*refpt(5)+theta0+theta2)
        uprimf=qmcc*(ef*costhf+ef2*costhf2)/betaf
        qpwf=0.5*qmcc*Scxl/(uf*ww)
        dltf=Scxl*(0.5*uprimf/uf-qpwf*(e1f*sinthf+e1f2*sinthf2))
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

        end subroutine maplinear_TWS

        subroutine rk6i_TWS(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (TWS), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_TWS(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_TWS(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_TWS

        subroutine intfunc1_TWS(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (TWS), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55,s11tmp
        double precision :: zedge,escale,ez1,ezp1,ezpp1,ww,theta0,qmcc,&
                            sinphi,cosphi,rfdsgn
        double precision :: ez12,ezp12,ezpp12,theta2,cosphi2,sinphi2
        integer :: my_rank, ierr

        zedge = this%Param(1)
        call getaxfldE_TWS(t,this,ez1,ezp1,ezpp1,ez12,ezp12,ezpp12)
        ww = this%Param(3)/Scfreq 
        theta0 = this%Param(4)*asin(1.0)/90
        theta2 = this%Param(12)*asin(1.0)/90
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
        sinphi2=sin(ww*y(5)+theta0+theta2)
        cosphi2=cos(ww*y(5)+theta0+theta2)
        rfdsgn=ez1*cosphi + ez12*cosphi2
        f(5)=1.0/(beta0*Scxl)
        f(6)=-qmcc*rfdsgn

        ! matrix elements
        s11tmp=0.5e0*(1.0+0.5*gamma0**2)* &
         (qmcc*rfdsgn/beta0**2/gamma0**2)**2 + &
          qmcc*0.5/beta0**3/gamma0**3*(ez1*sinphi+ez12*sinphi2)*ww/Scxl
        s11=s11tmp*Scxl
        s33=s11tmp*Scxl
        s55=  &
         -1.5e0*qmcc/beta0**2/gamma0*(ezp1*cosphi+ezp12*cosphi2) &
         +(beta0**2+0.5)/beta0**3/gamma0*qmcc*&
          (ez1*sinphi+ez12*sinphi2)*ww/Scxl &
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

        end subroutine intfunc1_TWS

        !get the E field on the axis
        subroutine getaxfldE_TWS(z,this,ez1,ezp1,ezpp1,ez12,ezp12,ezpp12)
        implicit none
        include 'mpif.h'
        type (TWS), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1,&
                                         ez12,ezp12,ezpp12
        double precision :: zedge,escale,len,tt,xl
        double precision :: f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc
        integer :: i,ntmp,numpar1,numpar2
        double precision :: epstol

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
        epstol = 1.0e-10
    
        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ez12 = 0.0
        ezp12 = 0.0
        ezpp12 = 0.0
        !if((z.gt.zedge).and.(z.lt.(zedge+len))) then
        !if((z.ge.zedge).and.(z.le.(zedge+len))) then
        if(((z.ge.zedge).and.(z.le.(zedge+len))) .or. &
           (abs(z-zedge).le.epstol) .or. &
           (abs(z-zedge-len).le.epstol) ) then
          !move into the tank local coordinate.
          zlc=z-zedge
          numpar1 = Fcoef(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = Fcoef(2)
          zend1 = Fcoef(3)
          zlength1 = Fcoef(4)
          !if( (zlc.ge.zstart1).and.(zlc.le.zend1)) then
          if( ((zlc.ge.zstart1).and.(zlc.le.zend1)) .or. &
              (abs(zlc-zstart1).le.epstol) .or. &
              (abs(zlc-zend1).le.epstol) ) then
            zmid = zlength1/2
            zz = z-zedge-zstart1-zmid
            zlen = zlength1
            !first 4 parameters in Fcoeft are not Forier coefficients.
            ntmp = 4
            !//find the field on the axis and its 1rst, 2nd, and 3rd
            !//derivatives which will be used to calculate the off-axis
            !//field distribution. 
            ez1 = Fcoef(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            do i = 2,(numpar1-1)/2+1
              ez1 = ez1 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
            enddo
            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
          else
          endif

          !//# of parameters for solenoid B fields.
          numpar2 = Fcoef(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = Fcoef(6+numpar1) 
          zend2 = Fcoef(7+numpar1) 
          zlength2 = Fcoef(8+numpar1) 
          !if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
          if( ((zlc.ge.zstart2).and.(zlc.le.zend2)) .or. &
              (abs(zlc-zstart2).le.epstol) .or. &
              (abs(zlc-zend2).le.epstol) ) then
            zmid = zlength2/2
            zz = z - zedge - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez12 = Fcoef(ntmp+1)/2
            ezp12 = 0.0
            ezpp12 = 0.0
            do i = 2,(numpar2-1)/2+1
              ez12 = ez12 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp12 = ezp12 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp12 = ezpp12+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
            enddo

            ez12=ez12*escale
            ezp12=ezp12*escale
            ezpp12=ezpp12*escale
          else
          endif
        else
        endif

        !write(14,*)z,ez1,ezp1,ezpp1,ez12,ezp12,ezpp12 

        end subroutine getaxfldE_TWS

        !get external field without displacement and rotation errors.
        subroutine  getfld_TWS(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (TWS), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,theta2
        integer :: i,ntmp,numpar1,numpar2

        clite = 299792458.e0
        pi = 2*asin(1.0)

        zedge = this%Param(1)
        escale = this%Param(2)
        len = this%Length
!        print*,"zedge: ",zedge,len,pos(3)
    
        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        if((pos(3).gt.zedge).and.(pos(3).lt.(zedge+len))) then
          ww = this%Param(3)*2*pi
!          xl = clite/ww  !real frequency has to be used here
          xlrep = ww/clite  !real frequency has to be used here
          theta0 = this%Param(4)*asin(1.0)/90
          !move into the tank local coordinate.
          zlc=pos(3)-zedge
          numpar1 = Fcoef(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = Fcoef(2)
          zend1 = Fcoef(3)
          zlength1 = Fcoef(4)
          extfld = 0.0
!          print*,"zstart1: ",zstart1,zend1,zlc,numpar1
          if( (zlc.ge.zstart1).and.(zlc.le.zend1)) then
            zmid = zlength1/2
            zz = pos(3)-zedge-zstart1-zmid
            zlen = zlength1
            !first 4 parameters in Fcoeft are not Forier coefficients.
            ntmp = 4
            !//find the field on the axis and its 1rst, 2nd, and 3rd
            !//derivatives which will be used to calculate the off-axis
            !//field distribution. 
            ez1 = Fcoef(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar1-1)/2+1
              ez1 = ez1 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
!               print*,Fcoef(2*(i-ntmp)-2)
!               print*,Fcoef(2*(i-ntmp)-1)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)/(2*pi*Scfreq)
            tmpcos = cos(ww*tt+theta0)
            tmpsin = sin(ww*tt+theta0)
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
            r2 = pos(1)**2+pos(2)**2
            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = pos(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -pos(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
            extfld = 0.0
          endif

!          print*,"strange",zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,ww,tt,theta0
!          write(10,101)pos(3),zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,tt,theta0
!101       format(10(1x,e14.5))
!          call flush_(10)

          !//# of parameters for solenoid B fields.
          numpar2 = Fcoef(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = Fcoef(6+numpar1) 
          zend2 = Fcoef(7+numpar1) 
          zlength2 = Fcoef(8+numpar1) 
          theta2 = this%Param(12)*pi/180
          if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
            zmid = zlength2/2
            zz = pos(3) - zedge - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez1 = Fcoef(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar2-1)/2+1
              ez1 = ez1 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo

            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)/(2*pi*Scfreq)
            tmpcos = cos(ww*tt+theta0+theta2)
            tmpsin = sin(ww*tt+theta0+theta2)
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
            r2 = pos(1)**2+pos(2)**2
            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = pos(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -pos(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
          endif
        else
          extfld = 0.0
        endif

        end subroutine getfld_TWS

        !get external field with displacement and rotation errors.
        subroutine  getflderr_TWS(pos,extfld,this,dx,dy,anglex,angley,&
                                    anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (TWS), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,theta2
        integer :: i,ntmp,numpar1,numpar2
        double precision :: dx,dy,anglex,angley,anglez
        double precision, dimension(3) :: temp,tmp

        clite = 299792458.e0
        pi = 2*asin(1.0)

        len = this%Length
        zedge = this%Param(1)
        escale = this%Param(2)
        !dx = this%Param(7)
        !dy = this%Param(8)
        !anglex = this%Param(9)
        !angley = this%Param(10)
        !anglez = this%Param(11)
!        print*,"zedge: ",zedge,len,pos(3)

        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
        ezppp = 0.0
        if((pos(3).gt.zedge).and.(pos(3).lt.(zedge+len))) then

          ww = this%Param(3)*2*pi
!          xl = clite/ww  !real frequency has to be used here
          xlrep = ww/clite  !real frequency has to be used here
          theta0 = this%Param(4)*asin(1.0)/90
          !move into the tank local coordinate.
          !zlc=pos(3)-zedge
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
          zlc = tmp(3)

          numpar1 = Fcoef(1)+0.1
          !//Here, zstart1 is the starting RF field location with 
          !//respect to the "zedge " as origin.
          !zstart1 can be negative
          zstart1 = Fcoef(2)
          zend1 = Fcoef(3)
          zlength1 = Fcoef(4)
          extfld = 0.0
!          print*,"zstart1: ",zstart1,zend1,zlc,numpar1
          if( (zlc.ge.zstart1).and.(zlc.le.zend1)) then
            zmid = zlength1/2
            !zz = pos(3)-zedge-zstart1-zmid
            zz = zlc-zstart1-zmid
            zlen = zlength1
            !first 4 parameters in Fcoeft are not Forier coefficients.
            ntmp = 4
            !//find the field on the axis and its 1rst, 2nd, and 3rd
            !//derivatives which will be used to calculate the off-axis
            !//field distribution. 
            ez1 = Fcoef(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar1-1)/2+1
              ez1 = ez1 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
!               print*,Fcoef(2*(i-ntmp)-2)
!               print*,Fcoef(2*(i-ntmp)-1)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo
            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)/(2*pi*Scfreq)
            tmpcos = cos(ww*tt+theta0)
            tmpsin = sin(ww*tt+theta0)
!            f1 = -(ezpp1+ez1/xl/xl)/4
!            f1p = -(ezppp+ezp1/xl/xl)/4
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
!            f1 = 0.0
!            f1p = 0.0
            !r2 = pos(1)**2+pos(2)**2
            r2 = tmp(1)**2+tmp(2)**2
            tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = tmp(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -tmp(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
!            tmpex = -pos(1)*(ezp1/2+f1p*r2/4)
!            tmpey = -pos(2)*(ezp1/2+f1p*r2/4)
!            tmpez = (ez1+f1*r2)
!            tmpbx = pos(2)/(xl*clite)*(ez1/2+f1*r2/4)
!            tmpby = -pos(1)/(xl*clite)*(ez1/2+f1*r2/4)
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
            extfld = 0.0
          endif

!          print*, zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,ww,tt,theta0
!          write(10,101)pos(3),zlc,ez1,ezp1,ezpp1,tmpez,tmpcos,tmpsin,tt,theta0
!101       format(10(1x,e14.5))
!          call flush(10)

          !//# of parameters for solenoid B fields.
          numpar2 = Fcoef(5+numpar1) + 0.1
          !zstart2 can be negative
          zstart2 = Fcoef(6+numpar1) 
          zend2 = Fcoef(7+numpar1) 
          zlength2 = Fcoef(8+numpar1) 
          theta2 = this%Param(12)
          if( (zlc.ge.zstart2) .and. (zlc.le.zend2)) then
            zmid = zlength2/2
            !zz = pos(3) - zedge - zstart2 - zmid
            zz = zlc - zstart2 - zmid
            zlen = zlength2
            ntmp = numpar1+8
            !//find the B field on axis and its 1st,2nd,3rd derivatives
            !//which will be used to calculate the offaxis B field
            ez1 = Fcoef(ntmp+1)/2
            ezp1 = 0.0
            ezpp1 = 0.0
            ezppp = 0.0
            do i = 2,(numpar2-1)/2+1
              ez1 = ez1 + Fcoef(2*i-2+ntmp)*&
                cos((i-1)*2*pi*zz/zlen) + &
                Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen)
              ezp1 = ezp1 + (i-1)*2*pi/zlen*&
                (-Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)+&
                Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
              ezpp1 = ezpp1+((i-1)*2*pi/zlen)**2*&
                (-Fcoef(2*i-2+ntmp)*cos((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*sin((i-1)*2*pi*zz/zlen))
              ezppp = ezppp+((i-1)*2*pi/zlen)**3*&
                (Fcoef(2*i-2+ntmp)*sin((i-1)*2*pi*zz/zlen)&
                - Fcoef(2*i-1+ntmp)*cos((i-1)*2*pi*zz/zlen))
            enddo

            ez1=ez1*escale
            ezp1=ezp1*escale
            ezpp1=ezpp1*escale
            ezppp = ezppp*escale
            tt = pos(4)/(2*pi*Scfreq)
            tmpcos = cos(ww*tt+theta0+theta2)
            tmpsin = sin(ww*tt+theta0+theta2)
            f1 = -(ezpp1+ez1*xlrep*xlrep)/4
            f1p = -(ezppp+ezp1*xlrep*xlrep)/4
            !r2 = pos(1)**2+pos(2)**2
            r2 = tmp(1)**2+tmp(2)**2
            tmpex = -tmp(1)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpey = -tmp(2)*(ezp1/2+f1p*r2/4)*tmpcos
            tmpez = (ez1+f1*r2)*tmpcos
            tmpbx = tmp(2)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpby = -tmp(1)*xlrep/clite*(ez1/2+f1*r2/4)*tmpsin
            tmpbz = 0.0
            extfld(1) = extfld(1) + tmpex
            extfld(2) = extfld(2) + tmpey
            extfld(3) = extfld(3) + tmpez
            extfld(4) = extfld(4) + tmpbx
            extfld(5) = extfld(5) + tmpby
            extfld(6) = extfld(6) + tmpbz
          else
!            extfld = 0.0
          endif
!          print*,zlc,extfld
          !for E field
          tmp(1) = extfld(1)
          tmp(2) = extfld(2)*cos(anglex)-extfld(3)*sin(anglex)
          tmp(3) = extfld(2)*sin(anglex)+extfld(3)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(1) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(2) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(3) = temp(3)
 
          !for B field
          tmp(1) = extfld(4)
          tmp(2) = extfld(5)*cos(anglex)-extfld(6)*sin(anglex)
          tmp(3) = extfld(5)*sin(anglex)+extfld(6)*cos(anglex)
          temp(1) = tmp(1)*cos(angley)-tmp(3)*sin(angley)
          temp(2) = tmp(2)
          temp(3) = tmp(1)*sin(angley)+tmp(3)*cos(angley)
          extfld(4) = temp(1)*cos(anglez) - temp(2)*sin(anglez)
          extfld(5) = temp(1)*sin(anglez) + temp(2)*cos(anglez)
          extfld(6) = temp(3)
        else
          extfld = 0.0
        endif

        end subroutine getflderr_TWS
      end module TWSclass
