!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BPMclass: Beam position monitor class in Lattice module of APPLICATION 
!           layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the different beam diagnostics at given
!              beam position.
! Comments:
!  1) Itype = -1, shift the transverse centroid position to 0.
!  2) Itype = -2, shift the transverse centroid position and angle to 0.
!                 (this one not work yet due to conflict of definition)
!  3) Itype = -10, mismatch the beam distribution by the amount given in
!                  Param(3) - Param(10).  
!  4) Itype = -13, collimator slit
!  5) Itype = -21, shift the beam centroid in 6D phase space by the amount
!                  given in Param(3) - Param(10).
!----------------------------------------------------------------
      module BPMclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 10
        type BPM
          !Itype < 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
          !      (3) : xmax
          !      (4) : pxmax
          !      (5) : ymax
          !      (6) : pymax
          !      (7) : zmax
          !      (8) : pzmax
        end type BPM
        interface getparam_BPM
          module procedure getparam1_BPM,  &
                          getparam2_BPM,   &
                          getparam3_BPM
        end interface
        interface setparam_BPM
          module procedure setparam1_BPM,  &
                           setparam2_BPM, setparam3_BPM
        end interface
      contains
        subroutine construct_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_BPM
   
        subroutine setparam1_BPM(this,i,value)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_BPM

        subroutine setparam2_BPM(this,values)
        implicit none
        type (BPM), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param(1:Nparam) = values(1:Nparam)

        end subroutine setparam2_BPM

        subroutine setparam3_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_BPM
   
        subroutine getparam1_BPM(this,i,blparam) 
        implicit none 
        type (BPM), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_BPM
  
        subroutine getparam2_BPM(this,blparams)
        implicit none
        type (BPM), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams(1:Nparam) = this%Param(1:Nparam)

        end subroutine getparam2_BPM

        subroutine getparam3_BPM(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (BPM), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_BPM

        subroutine shift_BPM(Pts1,itype,innp,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itype,innp,nptot
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: x0lc,px0lc,y0lc,py0lc
        double precision, dimension(4) :: tmplc,tmpgl
        integer :: i,j,ierr

        tmplc = 0.0
        tmpgl = 0.0
        if(itype.eq.(-2)) then
          x0lc = 0.0
          px0lc = 0.0
          y0lc = 0.0
          py0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            px0lc = px0lc + Pts1(2,i)
            y0lc = y0lc + Pts1(3,i)
            py0lc = py0lc + Pts1(4,i)
          enddo

          tmplc(1) = x0lc
          tmplc(2) = px0lc
          tmplc(3) = y0lc
          tmplc(4) = py0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(2) = tmpgl(2)/nptot
          tmpgl(3) = tmpgl(3)/nptot
          tmpgl(4) = tmpgl(4)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(2,i) = Pts1(2,i) - tmpgl(2)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
            Pts1(4,i) = Pts1(4,i) - tmpgl(4)
          enddo
        else if(itype.eq.(-1)) then
          x0lc = 0.0
          y0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            y0lc = y0lc + Pts1(3,i)
          enddo

          tmplc(1) = x0lc
          tmplc(3) = y0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(3) = tmpgl(3)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
          enddo
        else
        endif

        end subroutine shift_BPM

        !mismatch the beam at given location.
        !Here, the storage Param(3:8) is used to store the mismatch factors
        subroutine scale_BPM(Pts1,innp,xmis,pxmis,ymis,pymis,zmis,pzmis)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xmis,pxmis,ymis,pymis,zmis,pzmis
        integer :: i
 
        do i = 1, innp
            Pts1(1,i) = xmis*Pts1(1,i)
            Pts1(2,i) = pxmis*Pts1(2,i)
            Pts1(3,i) = ymis*Pts1(3,i)
            Pts1(4,i) = pymis*Pts1(4,i)
            Pts1(5,i) = zmis*Pts1(5,i)
            Pts1(6,i) = pzmis*Pts1(6,i)
        enddo
 
        end subroutine scale_BPM

        !shift the beam centroid in the 6D phase space.
        !This element can be used to model steering magnet etc.
        !Here, the storage Param(3:8) is used to store the amount of shift.
        !drange(3); shift in x (m)
        !drange(4); shift in Px (rad)
        !drange(5); shift in y (m)
        !drange(6); shift in Py (rad)
        !drange(7); shift in z (deg)
        !drange(8); shift in Pz (MeV)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_BPM(Pts1,innp,xshift,pxshift,yshift,pyshift,zshift,&
                   pzshift,gam,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xshift,pxshift,yshift,pyshift,zshift,pzshift
        double precision, intent(in) :: gam,mass
        integer :: i
        double precision :: gambetz
 
        do i = 1, innp
            gambetz = sqrt((gam-Pts1(6,i))**2-Pts1(2,i)**2-Pts1(4,i)**2-1.0)
            Pts1(1,i) = Pts1(1,i)+xshift/Scxl
            Pts1(2,i) = Pts1(2,i)+pxshift*gambetz
            Pts1(3,i) = Pts1(3,i)+yshift/Scxl
            Pts1(4,i) = Pts1(4,i)+pyshift*gambetz
            Pts1(5,i) = Pts1(5,i)+zshift/Rad2deg
            Pts1(6,i) = Pts1(6,i)+pzshift*1.0e6/mass
        enddo
 
        end subroutine kick_BPM

        !kick the beam longitudinally by the rf nonlinearity (the linear
        !part has been included in the map integrator and substracted.)
        !drange(3); vmax (V)
        !drange(4); phi0 (degree)
        !drange(5); horm number of rf
        subroutine kickRF_BPM(Pts1,innp,vmax,phi0,horm,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: vmax,phi0,mass,horm
        integer :: i
        real*8 :: vtmp,phi0lc,sinphi0,cosphi0  
 
        vtmp = vmax/mass 
        phi0lc = phi0*asin(1.0)/90
        sinphi0 = sin(phi0lc)
        cosphi0 = cos(phi0lc)
        do i = 1, innp
            Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)*horm+phi0lc)+ &
                                         cosphi0*Pts1(5,i)*horm)
            !Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)/horm+phi0lc)+ &
            !                             cosphi0*Pts1(5,i)/horm)
            !Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)/horm+phi0lc))
        enddo
 
        end subroutine kickRF_BPM

        !rotate the beam along vertical x-axis at given location
        !btype = -16
        subroutine xrot_BPM(Pts1,innp,phi,gam0)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        real*8, intent(in) :: gam0
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: phi
        integer :: i
        real*8 :: tmpx,tmpy,tmppx,tmppy,gam,gambetz

        do i = 1, innp
          gam = gam0 - Pts1(6,i)
          gambetz = sqrt(gam**2-1.0-Pts1(2,i)**2-Pts1(4,i)**2)

          tmpx = Pts1(3,i)*cos(phi)
          tmpy = -Pts1(3,i)*sin(phi)
          tmppx = Pts1(4,i)*cos(phi)+gambetz*sin(phi)
          tmppy = -Pts1(4,i)*sin(phi)+gambetz*cos(phi)

          Pts1(1,i) = Pts1(1,i) -tmpy*Pts1(2,i)/tmppy
          Pts1(3,i) = tmpx - tmpy*tmppx/tmppy
          Pts1(4,i) = tmppx
          Pts1(5,i) = Pts1(5,i) + tmpy*gam/tmppy
        enddo

        end subroutine xrot_BPM


        !rotate the beam along vertical y-axis at given location
        !btype = -17
        subroutine yrot_BPM(Pts1,innp,phi,gam0)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        real*8, intent(in) :: gam0
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: phi
        integer :: i
        real*8 :: tmpx,tmpy,tmppx,tmppy,gam,gambetz

        do i = 1, innp
          gam = gam0 - Pts1(6,i)
          gambetz = sqrt(gam**2-1.0-Pts1(2,i)**2-Pts1(4,i)**2)

          tmpx = Pts1(1,i)*cos(phi)
          tmpy = -Pts1(1,i)*sin(phi)
          tmppx = Pts1(2,i)*cos(phi)+gambetz*sin(phi)
          tmppy = -Pts1(2,i)*sin(phi)+gambetz*cos(phi)
          Pts1(1,i) = tmpx - tmpy*tmppx/tmppy
          Pts1(2,i) = tmppx
          Pts1(3,i) = Pts1(3,i) -tmpy*Pts1(4,i)/tmppy
          Pts1(5,i) = Pts1(5,i) + tmpy*gam/tmppy
        enddo

        end subroutine yrot_BPM

        !rotate the beam along longitudinal s-axis at given location
        !using the MAD8 notation.
        !btype = -18
        subroutine srot_BPM(Pts1,innp,phi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: phi
        integer :: i
        real*8 :: tmpx,tmpy,tmppx,tmppy
 
        do i = 1, innp
            tmpx = Pts1(1,i)*cos(phi)+Pts1(3,i)*sin(phi)
            tmpy = -Pts1(1,i)*sin(phi)+Pts1(3,i)*cos(phi)
            tmppx = Pts1(2,i)*cos(phi)+Pts1(4,i)*sin(phi)
            tmppy = -Pts1(2,i)*sin(phi)+Pts1(4,i)*cos(phi)
            Pts1(1,i) = tmpx
            Pts1(2,i) = tmppx
            Pts1(3,i) = tmpy
            Pts1(4,i) = tmppy
        enddo
 
        end subroutine srot_BPM

        !shift the centroid of the longitudinal phase space
        subroutine shiftlong_BPM(Pts1,itype,innp,nptot,tc,ptc)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itype,innp
        integer, intent(in) :: nptot
        real*8, intent(inout) :: tc,ptc
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: x0lc,px0lc,y0lc,py0lc
        double precision, dimension(2) :: tmplc,tmpgl
        integer :: i,j,ierr

        tmplc = 0.0
        tmpgl = 0.0
          x0lc = 0.0
          y0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(5,i)
            y0lc = y0lc + Pts1(6,i)
          enddo

          tmplc(1) = x0lc
          tmplc(2) = y0lc

!          print*,"tmplc: ",tmplc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,2,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(2) = tmpgl(2)/nptot

          do i = 1, innp
            Pts1(5,i) = Pts1(5,i) - tmpgl(1)
            Pts1(6,i) = Pts1(6,i) - tmpgl(2)
          enddo

!          print*,"tmpgl: ",tmpgl

          tc = tc + tmpgl(1)
          ptc = ptc + tmpgl(2)

        end subroutine shiftlong_BPM

        !The following subroutine increases the uncorrelated energy spread by heating
        subroutine engheater_BPM(Pts1,innp,b0,qmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        real*8 :: b0,qmass
        real*8, dimension(innp) :: rd
        integer :: i
        real*8 :: sigeng

        sigeng = b0*qmass
        call normVec1(rd,innp)

        do i = 1, innp
          Pts1(6,i) = Pts1(6,i) - sigeng*rd(i)
        enddo

        end subroutine engheater_BPM

        subroutine normVec1(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i
 
        epsilon = 1.0d-18
 
        twopi = 4.0d0*asin(1.0d0)
        call random_number(x2)
        call random_number(x1)
        do i = 1, num
          if(x1(i).eq.0.0d0) x1(i) = epsilon
          y(i) = sqrt(-2.0d0*log(x1(i)))*cos(twopi*x2(i))
        enddo
 
        end subroutine normVec1

!        k0 = this%Param(3) !dipole
!        k1 = this%Param(4) !quad
!        k2 = this%Param(5) !sext
!        k3 = this%Param(6) !oct
!        k4 = this%Param(7) !dec
!        k5 = this%Param(8) !dodec
        subroutine thinMultipole_BPM(Pts1in,innp,k0,k1,k2,k3,k4,k5,gam0,&
                                     qmass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1in
        double precision, intent(in) :: k0,k1,k2,k3,k4,k5,gam0,qmass
        integer :: i
        real*8 :: gam,gambet,dpp,bet,xx,yy,gambet0

        gambet0 = sqrt(gam0**2-1.0d0)

        do i = 1, innp
          gam = gam0-Pts1in(6,i)
          gambet = sqrt(gam**2-1.0d0)
          dpp = (gambet-gambet0)/gambet0
          bet = gambet/gam
          xx = Pts1in(1,i)*Scxl
          yy = Pts1in(3,i)*Scxl
          !for single charge state
          !Pts1in(2,i) = Pts1in(2,i) + gambet0* &
          Pts1in(2,i) = Pts1in(2,i) + Pts1in(7,i)/qmass*gambet0* &
          (-k0+k0*dpp-k1*xx-k2/2*(xx**2-yy**2)-&
          k3/6*(xx**3-3*xx*yy**2)-k4/24*((xx**2-yy**2)**2-4*xx**2*yy**2)&
           - k5/120.0d0*(((xx**2-yy**2)**2-4*xx**2*yy**2)*xx - &
                      4*xx*yy**2*(xx**2-yy**2)))
          !for single charge state
          !Pts1in(4,i) = Pts1in(4,i) + gambet0* &
          Pts1in(4,i) = Pts1in(4,i) + Pts1in(7,i)/qmass*gambet0* &
                        (k1*yy+k2*xx*yy+k3/6*(3*yy*xx**2-yy**3)+&
                        k4/6*xx*yy*(xx**2-yy**2) + &
                        k5/120.0d0*(4*xx**2*yy*(xx**2-yy**2)+yy*(xx**2-yy**2)**2-&
                        4*xx**2*yy**3) )
!Lg frozen
          Pts1in(5,i) = Pts1in(5,i) + k0*xx/bet/Scxl
        enddo

        end subroutine thinMultipole_BPM

        !RF deflecting cavity thin lens map
        !drange(2); sk
        !drange(3); flaghv (flag for horizontal or vertical deflecting)
        !thin lens RF deflection kick
        subroutine kickthindef_BPM(Pts1,innp,sk,flaghv,gam0)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,flaghv
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: sk,gam0
        integer :: i
        real*8 :: gam,gambetz,x,z

!        print*,"sk:",sk,flaghv,innp

        gam = gam0
        gambetz = sqrt(gam**2-1.0d0)
        if(flaghv.ge.0) then !H-Z deflecting
          do i = 1, innp
!            gam = gam0 - Pts1(6,i)
!            gambetz = sqrt(gam**2-1.0d0-Pts1(2,i)**2-Pts1(4,i)**2)
            x = Pts1(1,i)*Scxl
            z = -Pts1(5,i)*Scxl
            Pts1(2,i) = Pts1(2,i)+sk*z*gambetz
            Pts1(6,i) = Pts1(6,i)-sk*x*gam
          enddo
        else !V-Z deflecting
          do i = 1, innp
!            gam = gam0 - Pts1(6,i)
!            gambetz = sqrt(gam**2-1.0d0-Pts1(2,i)**2-Pts1(4,i)**2)
            x = Pts1(3,i)*Scxl
            z = -Pts1(5,i)*Scxl
            Pts1(4,i) = Pts1(4,i)+sk*z*gambetz
            Pts1(6,i) = Pts1(6,i)-sk*x*gam
          enddo
        endif

        end subroutine kickthindef_BPM

        !instant kick by external linear map
        subroutine kickextmap_BPM(pts1,innp,gam0,nfile)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,nfile
        double precision, pointer, dimension(:,:) :: pts1
        real*8, dimension(6,6) :: extmap
        double precision, intent(in) :: gam0
        integer :: i,j,k
        real*8 :: gam,gambet,bet
        real*8, dimension(6) :: tmp,temp

!        open(17,file="linearmap.in",status="old")
        do i = 1, 6
          read(nfile,*)extmap(i,1),extmap(i,2),extmap(i,3),&
                    extmap(i,4),extmap(i,5),extmap(i,6)
        enddo
        close(nfile)
        gam = gam0
        gambet = sqrt(gam**2-1.0d0)
        bet = gambet/gam
!        print*,"gam0: ",gam0,extmap(1,1),extmap(6,6),scxl,gambet,bet
        do i = 1, innp
          !convert to x,x',y,y',z,dp/p coordinates
          tmp(1) = pts1(1,i)*scxl
          tmp(2) = pts1(2,i)/gambet
          tmp(3) = pts1(3,i)*scxl
          tmp(4) = pts1(4,i)/gambet
          tmp(5) = -pts1(5,i)*scxl*bet
          tmp(6) = -pts1(6,i)/(gam*bet**2)
          do k = 1, 6
            temp(k) = 0.0d0
            do j = 1, 6
              temp(k) = temp(k)+extmap(k,j)*tmp(j)
            enddo
          enddo
          pts1(1,i) = temp(1)/scxl
          pts1(2,i) = temp(2)*gambet 
          pts1(3,i) = temp(3)/scxl
          pts1(4,i) = temp(4)*gambet 
          pts1(5,i) = -temp(5)/(scxl*bet)
          pts1(6,i) = -temp(6)*gam*bet**2
        enddo

        end subroutine kickextmap_BPM
      end module BPMclass
