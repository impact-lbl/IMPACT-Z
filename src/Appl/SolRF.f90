!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! SolRFclass: Solenoid with imbeded RF field beam line element class
!             in Lattice module of APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the linear transfer map and RF field
!              for the Sol-RF beam line elment.
! Comments:
!----------------------------------------------------------------
      module SolRFclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 15
        type SolRF
          !Itype = 105
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
          !      (12) : Bz0
          !      (13) : aawk: aperture size for the wakefield
          !      (14) : ggwk: gap size for the wakefield
          !      (15) : lengwk: length for the wakefield
        end type SolRF
        interface getparam_SolRF
          module procedure getparam1_SolRF,  &
                          getparam2_SolRF,   &
                          getparam3_SolRF
        end interface
        interface setparam_SolRF
          module procedure setparam1_SolRF,  &
                          setparam2_SolRF, setparam3_SolRF
        end interface
      contains
        subroutine construct_SolRF(this,numseg,nmpstp,type,blength)
        implicit none
        type (SolRF), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_SolRF
   
        subroutine setparam1_SolRF(this,i,value)
        implicit none
        type (SolRF), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_SolRF

        subroutine setparam2_SolRF(this,values)
        implicit none
        type (SolRF), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_SolRF

        subroutine setparam3_SolRF(this,numseg,nmpstp,type,blength)
        implicit none
        type (SolRF), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_SolRF
   
        subroutine getparam1_SolRF(this,i,blparam) 
        implicit none 
        type (SolRF), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_SolRF
  
        subroutine getparam2_SolRF(this,blparams)
        implicit none
        type (SolRF), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_SolRF

        subroutine getparam3_SolRF(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (SolRF), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_SolRF
       
        subroutine maplinear_SolRF(t,tau,xm,this,refpt,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,tau,Bchg,Bmass
        double precision, dimension(6,6), intent(out) :: xm
        double precision, dimension(6), intent(inout) :: refpt
        type (SolRF), intent(in) :: this
        double precision, dimension(22) :: y
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

        call getaxfldE_SolRF(t,this,ein,e1in,e2in)
        sinthi=sin(ww*refpt(5)+theta0)
        costhi=cos(ww*refpt(5)+theta0)
        uprimi=qmcc*ein/betai*costhi
        qpwi=0.5*qmcc*Scxl/(ui*ww)
        dlti=Scxl*(0.5*uprimi/ui-qpwi*e1in*sinthi)
        thli=1.5*Scxl*(uprimi/ui)

        call rk6i_SolRF(h,mpstp,t0,y,22,this,Bchg,Bmass)

        refpt(5)=y(5)
        refpt(6)=y(6)
        gf=-refpt(6)
        betaf=sqrt(gf**2-1.)/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin = t + tau
        call getaxfldE_SolRF(tfin,this,ef,e1f,e2f)
        sinthf=sin(ww*refpt(5)+theta0)
        costhf=cos(ww*refpt(5)+theta0)
        uprimf=qmcc*ef/betaf*costhf
        qpwf=0.5*qmcc*Scxl/(uf*ww)
        dltf=Scxl*(0.5*uprimf/uf-qpwf*e1f*sinthf)
        thlf=1.5*Scxl*(uprimf/uf)

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

        end subroutine maplinear_SolRF

        subroutine rk6i_SolRF(h,ns,t,y,nvar,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ns,nvar
        double precision, intent(inout) :: t
        double precision, intent(in) :: h,Bchg,Bmass
        double precision, dimension(nvar), intent(inout) :: y
        type (SolRF), intent(in) :: this
        double precision, dimension(nvar) :: yt,a,b,c,d,e,f,g,o,p 
        double precision:: tint,tt
        integer :: i
        tint=t
        do i=1,ns
          call intfunc1_SolRF(t,y,f,this,Bchg,Bmass)
          a=h*f
          yt=y+a/9.d0
          tt=t+h/9.d0
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          b=h*f
          yt=y + (a + 3.d0*b)/24.d0
          tt=t+h/6.d0
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          c=h*f
          yt=y+(a-3.d0*b+4.d0*c)/6.d0
          tt=t+h/3.d0
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          d=h*f
          yt=y + (278.d0*a - 945.d0*b + 840.d0*c + 99.d0*d)/544.d0
          tt=t+.5d0*h
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          e=h*f
          yt=y + (-106.d0*a+273.d0*b-104.d0*c-107.d0*d+48.d0*e)/6.d0
          tt = t+2.d0*h/3.d0
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          g=h*f
          yt = y+(110974.d0*a-236799.d0*b+68376.d0*c+ &
                  103803.d0*d-10240.d0*e + 1926.d0*g)/45648.d0
          tt = t + 5.d0*h/6.d0
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          o=h*f
          yt = y+(-101195.d0*a+222534.d0*b-71988.d0*c-  &
                  26109.d0*d-2.d4*e-72.d0*g+22824.d0*o)/25994.d0
          tt = t + h
          call intfunc1_SolRF(tt,yt,f,this,Bchg,Bmass)
          p=h*f
          y = y+(41.d0*a+216.d0*c+27.d0*d+ &
                 272.d0*e+27.d0*g+216.d0*o+41.d0*p)/840.d0
          t=tint+i*h
        enddo

        end subroutine rk6i_SolRF

        subroutine intfunc1_SolRF(t,y,f,this,Bchg,Bmass)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: t,Bchg,Bmass
        double precision, dimension(:), intent(in) :: y
        double precision, dimension(:), intent(out) :: f
        type (SolRF), intent(in) :: this
        double precision :: gamma0,beta0,gbet,s11,s33,s55,s11tmp
        double precision :: zedge,escale,ez1,ezp1,ezpp1,ww,theta0,qmcc,&
                            sinphi,cosphi,rfdsgn,brho,bgrad,C2alpha,b0
        integer :: my_rank, ierr

        zedge = this%Param(1)
        call getaxfldE_SolRF(t,this,ez1,ezp1,ezpp1)
        call getBgradfld_SolRF(t,this,b0,bgrad)
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
        s11=(s11tmp+(0.5*b0/brho)**2/2)*Scxl
        s33=(s11tmp+(0.5*b0/brho)**2/2)*Scxl
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
        C2alpha = -b0/brho/2 !rotation matrix of Solenoid
        f(19)=-C2alpha*y(21)
        f(20)=-C2alpha*y(22)
        f(21)=C2alpha*y(19)
        f(22)=C2alpha*y(20)
 
        end subroutine intfunc1_SolRF

        !interpolate the field from the SolRF rf cavity onto bunch location.
        subroutine getaxfldE_SolRF(z,this,ez1,ezp1,ezpp1)
        implicit none
        include 'mpif.h'
        type (SolRF), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1
        double precision:: zz,hstep,slope,zedge,escale,zlen
        integer :: klo,khi,k
        integer :: my_rank,ierr
        integer :: i, ntmp,numpar1
        double precision :: len,zstart1,zend1,zlength1,zlc,zmid
        double precision :: epstol

        zedge = this%Param(1)
        escale = this%Param(2)
!        zz=z-zedge
        len = this%Length
        epstol = 1.0e-10
   
        ez1 = 0.0
        ezp1 = 0.0
        ezpp1 = 0.0
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
        endif

!-------------------------------------------------------------------
! on axis field from interpolation
!        klo=1
!        khi=Ndata
!1       if(khi-klo.gt.1) then
!          k=(khi+klo)/2
!          if(zdat(k).gt.zz)then
!             khi=k
!          else
!             klo=k
!          endif
!          goto 1
!        endif
!        hstep=zdat(khi)-zdat(klo)
!        slope=(edat(khi)-edat(klo))/hstep
!        ez1 =edat(klo)+slope*(zz-zdat(klo))
!        slope=(epdat(khi)-epdat(klo))/hstep
!        ezp1=epdat(klo)+slope*(zz-zdat(klo))
!        ezpp1 = 0.0
!-------------------------------------------------------------------
! on axis field from analytical function
!        zlen = this%Length
!        pi = 2*asin(1.0)
!        ez1 = 30.0*sin(26.72*2*pi*zz/zlen)
!        ezp1 = 30.0*26.72*2*pi/zlen*cos(26.72*2*pi*zz/zlen)
!        ezpp1 = 0.0
!
!        ez1 = ez1*escale
!        ezp1 = ezp1*escale
!        ezpp1 = ezpp1*escale
!
        !print*,"ez1: ",ez1,escale,zedge,zz
!        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!        if(my_rank.eq.1) then
!         write(19,101)z,ez1/1.0e6,ezp1/1.0e6,zz,float(m),float(init)
!        endif
!        if(my_rank.eq.1) then
!         write(19,101)z,gt
!        endif
  101   format(6(1x,1pe13.6))

        end subroutine getaxfldE_SolRF

        subroutine getBgradfld_SolRF(z,this,b0,bgrad)
        implicit none
        include 'mpif.h'
        type (SolRF), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: b0,bgrad

        !uniform bz field.
        b0 = this%Param(12)
        bgrad = 0.0

        end subroutine getBgradfld_SolRF

        !get external field without displacement and rotation errors.
        subroutine  getfld_SolRF(pos,extfld,this)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,bscale
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
          bscale = this%Param(12)
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
            r2 = pos(1)*pos(1)+pos(2)*pos(2)
            extfld(4) = extfld(4)-bscale*ezp1/2*pos(1)+bscale*ezppp*pos(1)*r2/16  
            extfld(5) = extfld(5)-bscale*ezp1/2*pos(2)+bscale*ezppp*pos(2)*r2/16  
            extfld(6) = extfld(6)+bscale*ez1-bscale*ezpp1*r2/4
!            write(11,102)pos(3),zlc,ez1*bscale,ezp1*bscale,ezpp1,ezppp
!102         format(6(1x,e14.5))
!            call flush_(11)
          else
!            extfld = 0.0
          endif
!          print*,zlc,extfld
        else
          extfld = 0.0
        endif

        end subroutine getfld_SolRF

        !get external field with displacement and rotation errors.
        subroutine  getflderr_SolRF(pos,extfld,this,dx,dy,anglex,angley,&
                                    anglez)
        implicit none
        include 'mpif.h'
        double precision, dimension(4), intent(in) :: pos
        type (SolRF), intent(in) :: this
        double precision, dimension(6), intent(out) :: extfld
        double precision :: zedge,escale,theta0,ww,len,tt,xl,xlrep
        double precision :: ez1,ezp1,ezpp1,f1,clite,tmpsin,tmpcos,pi
        double precision :: tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: zz,bgrad,b0,zmid,rr,zstart1,zend1,zlength1,&
                   zstart2,zend2,zlength2,zlen,f1p,r2,zlc,ezppp,bscale
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
          bscale = this%Param(12)
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
            r2 = tmp(1)*tmp(1)+tmp(2)*tmp(2)
            extfld(4) = extfld(4)-bscale*ezp1/2*tmp(1)+bscale*ezppp*tmp(1)*r2/16  
            extfld(5) = extfld(5)-bscale*ezp1/2*tmp(2)+bscale*ezppp*tmp(2)*r2/16  
            extfld(6) = extfld(6)+bscale*ez1-bscale*ezpp1*r2/4
!            write(11,102)pos(3),zlc,ez1*bscale,ezp1*bscale,ezpp1,ezppp
!102         format(6(1x,e14.5))
!            call flush_(11)
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

        end subroutine getflderr_SolRF

      end module SolRFclass
