!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Ptclmgerclass: Particles moving manager class in Communication module 
!                of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines functions to transport particles to 
!              their local compuatation processor domain through an
!              iterative neighboring processor communcation process.
! Comments: We have added 3 more attributes to the particle array:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
        module Ptclmgerclass
          use Timerclass
          use Pgrid2dclass
        contains
        ! move particles from one processor to 4 neighboring processors.
        subroutine ptsmv2_ptclmger(Ptsl,Nptlocal,grid,pdim,npmax,lcrange)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(in) :: grid
        integer, intent(inout) :: Nptlocal
        integer, intent(in) :: pdim,npmax
        double precision, pointer, dimension(:,:) :: Ptsl
!        double precision, dimension(pdim,npmax) :: Ptsl
        double precision, dimension(:),intent(in) :: lcrange
        !integer, parameter :: nptmv = 100000
        integer, parameter :: nptmv = 3000000
!        integer, parameter :: nptmv = 300000
        double precision, dimension(9,nptmv) :: left,right,up,down
        !integer, parameter :: ntmpstr = 400000
        integer, parameter :: ntmpstr = 10000000
!        integer, parameter :: ntmpstr = 600000
        double precision, dimension(9,ntmpstr) :: temp1
        double precision, allocatable, dimension(:,:) :: recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer :: nsmall,i,j,numpts,ic
        integer, dimension(2) :: tmpcoord
        logical, dimension(Nptlocal) :: msk
        logical, allocatable, dimension(:) :: mmsk
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv
        integer :: msid,ierr,nrecv
        integer status(MPI_STATUS_SIZE) 
        integer statarry(MPI_STATUS_SIZE,4), req(4)
        integer :: flag,Nptlocal0,nst,nout0,iileft,iiright,iiup,iidown
        integer :: totflag
        double precision :: t0
        real*8, dimension(6) :: tmpctlc,tmpct,tmpct2

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        if(myidx.ne.(npx-1)) then
          myright = myidx + 1
        else
          myright = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          myleft = myidx - 1
        else
          myleft = MPI_PROC_NULL
        endif 

        if(myidy.ne.npy-1) then
          myup = myidy + 1
        else
          myup = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          mydown = myidy -1
        else
          mydown = MPI_PROC_NULL
        endif

!        call MPI_BARRIER(comm2d,ierr)

        tmpctlc = 0.0
!        do i = 1, Nptlocal
!          tmpctlc(1) = tmpctlc(1) + abs(Ptsl(1,i))
!          tmpctlc(2) = tmpctlc(2) + abs(Ptsl(2,i))
!          tmpctlc(3) = tmpctlc(3) + abs(Ptsl(3,i))
!          tmpctlc(4) = tmpctlc(4) + abs(Ptsl(4,i))
!          tmpctlc(5) = tmpctlc(5) + abs(Ptsl(5,i))
!          tmpctlc(6) = tmpctlc(6) + abs(Ptsl(6,i))
!        enddo
!        call MPI_ALLREDUCE(tmpctlc,tmpct,6,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                           comm2d,ierr)
!        if(myid.eq.0) then
!          print*,"before particle move: ",tmpct
!        endif

        flag = 0
        Nptlocal0 = Nptlocal
        nout0 = 0

        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        iileft = 0
        iiright = 0
        iiup = 0
        iidown = 0
        do i = 1, Nptlocal0 - nout0
          msk(i) = .true.
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              iileft = iileft + 1
              if(iileft.le.nptmv) then
                ileft = ileft + 1
                left(:,ileft) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iiright = iiright + 1
              if(iiright.le.nptmv) then
                iright = iright + 1
                right(:,iright) = Ptsl(:,i)
                msk(i) = .false.
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                  iiup = iiup + 1
                  if(iiup.le.nptmv) then
                    iup = iup + 1
                    up(:,iup) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                  iidown = iidown + 1
                  if(iidown.le.nptmv) then
                    idown = idown + 1
                    down(:,idown) = Ptsl(:,i)
                    msk(i) = .false.
                  endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iiup = iiup + 1
              if(iiup.le.nptmv) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              iidown = iidown + 1
              if(iidown.le.nptmv) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
              endif
            endif
          else
          endif

        enddo

        if((iileft.gt.nptmv).or.(iiright.gt.nptmv).or.(iiup.gt.nptmv) &
            .or.(iidown.gt.nptmv)) then
           print*,"too small buffer size, nptmv, for particle move!"
           flag = 1
           stop
        else
           flag = 0
        endif

        nout = ileft+iright+iup+idown

        allocate(recv(9,1))
        allocate(mmsk(1))
        recv = 0.0d0

        call MPI_BARRIER(comm2d,ierr)
!        temp1 = 0.0d0

        ic = 0

        do

        jleft = 0
        jright = 0

        call MPI_IRECV(jleft,1,MPI_INTEGER,myright,0,commrow,req(1),&
                      ierr)
        call MPI_IRECV(jright,1,MPI_INTEGER,myleft,0,commrow,req(2),&
                        ierr)
        call MPI_ISEND(ileft,1,MPI_INTEGER,myleft,0,commrow,req(3),&
                       ierr)
        call MPI_ISEND(iright,1,MPI_INTEGER,myright,0,commrow,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 1:"
!        endif

        jup = 0
        jdown = 0

        call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
                      ierr)
        call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
                      ierr)
        call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
                       ierr)
        call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
                       ierr)
        call MPI_WAITALL(4,req,statarry,ierr) 

        call MPI_BARRIER(commcol,ierr)
!        if(myid.eq.0) then
!          print*,"pass 2:"
!        endif

        numbuf = jleft+jright+jup+jdown 

        deallocate(recv)
        allocate(recv(9,numbuf))

        nmv0 = 0
        nst = nmv0 + 1
        !send outgoing particles to left neibhoring processor.
        jleft = 9*jleft
        ileft = 9*ileft
        call MPI_IRECV(recv(1,nst),jleft,MPI_DOUBLE_PRECISION,myright,&
                       0,commrow,msid,ierr)
        call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        ileft = ileft/9
        jleft = jleft/9
        nmv0 = nmv0+jleft
        
        nst = nmv0 + 1
        !send outgoing particles to right neibhoring processor.
        jright = 9*jright
        iright = 9*iright
        call MPI_IRECV(recv(1,nst),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,msid,ierr)
        call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                      0,commrow,ierr)
        call MPI_WAIT(msid,status,ierr) 
        iright = iright/9
        jright = jright/9
        nmv0 = nmv0 + jright

        call MPI_BARRIER(commrow,ierr)
!        if(myid.eq.0) then
!          print*,"pass 3:"
!        endif

        nst = nmv0 + 1
        !send outgoing particles to down neibhoring processor.
        jdown = 9*jdown
        idown = 9*idown
        call MPI_IRECV(recv(1,nst),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,msid,ierr)
        call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        idown = idown/9
        jdown = jdown/9
        nmv0 = nmv0 + jdown

        nst = nmv0 + 1
        !send outgoing particles to up neibhoring processor.
        jup = 9*jup
        iup = 9*iup
        call MPI_IRECV(recv(1,nst),jup,MPI_DOUBLE_PRECISION,mydown,&
                      0,commcol,msid,ierr)
        call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                      0,commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        iup = iup/9
        jup = jup/9

        call MPI_BARRIER(commcol,ierr)

        deallocate(mmsk)
        allocate(mmsk(numbuf))

        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = 0

        do i = 1, numbuf
          mmsk(i) = .true.
          ii = i+nmv0
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
            endif
          else
          endif
        enddo

        if((ileft.gt.nptmv).or.(iright.gt.nptmv).or.(iup.gt.nptmv) &
            .or.(idown.gt.nptmv)) then
           print*,"too small buffer size, nptmv, for particle move!"
           flag = 1
           stop
        endif

        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        do i = 1, numbuf
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic) = recv(:,i)
          endif
        enddo

        if(ic.gt.ntmpstr) then
          print*,"temporary storage size, ntmpstr, is too small!"
          stop
        endif

        if(totnmv.eq.0) then
          exit
        endif

        call MPI_BARRIER(comm2d,ierr)

        enddo

        nrecv = ic
        !print*,"nrecv: ",nrecv,myid,nout,Nptlocal

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        deallocate(recv)
        deallocate(mmsk)
        allocate(recv(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            do j = 1, 9
              recv(j,ic) = Ptsl(j,i)
            enddo
          endif
        enddo

        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nrecv

        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          do j = 1, 9
            Ptsl(j,i) = recv(j,i)
          enddo
        enddo
        do i = 1, nrecv
          ii = i + numpts
          do j = 1, 9
            Ptsl(j,ii) = temp1(j,i)
          enddo
        enddo

        deallocate(recv)

        tmpctlc = 0.0
!        do i = 1, Nptlocal
!          tmpctlc(1) = tmpctlc(1) + abs(Ptsl(1,i))
!          tmpctlc(2) = tmpctlc(2) + abs(Ptsl(2,i))
!          tmpctlc(3) = tmpctlc(3) + abs(Ptsl(3,i))
!          tmpctlc(4) = tmpctlc(4) + abs(Ptsl(4,i))
!          tmpctlc(5) = tmpctlc(5) + abs(Ptsl(5,i))
!          tmpctlc(6) = tmpctlc(6) + abs(Ptsl(6,i))
!        enddo
!        call MPI_ALLREDUCE(tmpctlc,tmpct2,6,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                           comm2d,ierr)
!        if(myid.eq.0) then
!          print*,"after particle move: ",tmpct2
!          print*,"difference of move: ",sum(tmpct2-tmpct)
!        endif

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv2_ptclmger
      end module Ptclmgerclass
