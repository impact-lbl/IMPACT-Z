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
        double precision, allocatable, dimension(:,:) :: left,right,up,down
        double precision, allocatable, dimension(:,:) :: temp1,recv
        integer :: myid,myidx,myidy,totnp,npy,npx, &
                   comm2d,commcol,commrow
        integer :: ileft,iright,iup,idown,iupright,iupleft,&
                   idownleft,idownright
        integer :: jleft,jright,jup,jdown,jupright,jupleft,&
                   jdownleft,jdownright
        integer :: myleft,myright,myup,mydown,myupright,myupleft,&
                   mydwnleft,mydwnright
        integer status(MPI_STATUS_SIZE)
        integer :: nsmall,i,j,numpts,ic,ierr
        integer, dimension(2) :: tmpcoord
        logical, allocatable, dimension(:) :: msk,mmsk
        double precision :: t0,t1
        integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij

        call starttime_Timer(t0)
!        call starttime_Timer(t1)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        if(Nptlocal.gt.480) then
          nsmall = Nptlocal/16
        else
          nsmall = Nptlocal
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))

        left = 0.0
        right = 0.0
        up = 0.0
        down = 0.0
        ileft = 0
        iright = 0
        iup = 0
        idown = 0

        allocate(msk(Nptlocal))
        msk = .true.

        do i = 1, Nptlocal
          if(Ptsl(5,i).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(5,i).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(Ptsl(3,i).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(Ptsl(3,i).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = Ptsl(:,i)
                msk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(Ptsl(3,i).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(Ptsl(3,i).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = Ptsl(:,i)
              msk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif

        enddo

!        t_ptsmv1 = t_ptsmv1 + elapsedtime_Timer(t1)
        nmv0 = 0
        nout = ileft+iright+iup+idown
!        print*,"in comm: ",Nptlocal,nout
        allocate(recv(9,nmv0+1))
        allocate(temp1(9,nmv0+1))
        ij = 0

        do

        ij = ij + 1
        myleft = myidx - 1
        myright = myidx + 1
        jleft = 0
        jright = 0

        if(myidx.ne.0) then
          call MPI_SEND(ileft,1,MPI_INTEGER,myleft,0,commrow,ierr)
        endif
        if(myidx.ne.(npx-1)) then
          call MPI_RECV(jleft,1,MPI_INTEGER,myright,0,commrow,status,&
                        ierr)
        endif

        if(myidx.ne.(npx-1)) then
          call MPI_SEND(iright,1,MPI_INTEGER,myright,0,commrow,ierr)
        endif
        if(myidx.ne.0) then
          call MPI_RECV(jright,1,MPI_INTEGER,myleft,0,commrow,status,&
                        ierr)
        endif

        myup = myidy + 1
        mydown = myidy - 1
        jup = 0
        jdown = 0

        if(myidy.ne.0) then
          call MPI_SEND(idown,1,MPI_INTEGER,mydown,0,commcol,ierr)
        endif
        if(myidy.ne.(npy-1)) then
          call MPI_RECV(jdown,1,MPI_INTEGER,myup,0,commcol,status,&
                        ierr)
        endif
        
        if(myidy.ne.(npy-1)) then
          call MPI_SEND(iup,1,MPI_INTEGER,myup,0,commcol,ierr)
        endif
        if(myidy.ne.0) then
          call MPI_RECV(jup,1,MPI_INTEGER,mydown,0,commcol,status,&
                        ierr)
        endif

        numbuf = jleft+jright+jup+jdown 
        
        deallocate(recv)
        allocate(recv(9,numbuf+nmv0+1))
        do i = 1, nmv0
          recv(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)

        !send outgoing particles to left neibhoring processor.
        allocate(temp1(9,jleft))
        temp1 = 0.0
        if(myidx.ne.0) then
          ileft = 9*ileft
          call MPI_SEND(left(1,1),ileft,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,ierr)
          ileft = ileft/9
        endif
        if(myidx.ne.(npx-1)) then
          jleft = 9*jleft
          call MPI_RECV(temp1(1,1),jleft,MPI_DOUBLE_PRECISION,myright,&
                        0,commrow,&
                        status,ierr)
          jleft = jleft/9
        endif
        do i = 1, jleft
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0+jleft
        deallocate(temp1)
        
!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to right neibhoring processor.
        allocate(temp1(9,jright))
        temp1 = 0.0
        if(myidx.ne.(npx-1)) then
          iright = 9*iright
          call MPI_SEND(right(1,1),iright,MPI_DOUBLE_PRECISION,myright,&
                        0,commrow,&
                        ierr)
          iright = iright/9
        endif
        if(myidx.ne.0) then
          jright = 9*jright
          call MPI_RECV(temp1(1,1),jright,MPI_DOUBLE_PRECISION,myleft,&
                        0,commrow,&
                        status,ierr)
          jright = jright/9
        endif
        do i = 1, jright
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jright
        deallocate(temp1)

!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to down neibhoring processor.
        allocate(temp1(9,jdown))
        temp1 = 0.0
        if(myidy.ne.0) then
          idown = 9*idown
          call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,&
                        ierr)
          idown = idown/9
        endif
        if(myidy.ne.(npy-1)) then
          jdown = 9*jdown
          call MPI_RECV(temp1(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,&
                        status,ierr)
          jdown = jdown/9
        endif
        do i = 1, jdown
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jdown
        deallocate(temp1)

!        call MPI_BARRIER(comm2d,ierr)
        !send outgoing particles to up neibhoring processor.
        allocate(temp1(9,jup))
        temp1 = 0.0
        if(myidy.ne.(npy-1)) then
          iup = 9*iup
          call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
                        0,commcol,&
                        ierr)
          iup = iup/9
        endif
        if(myidy.ne.0) then
          jup = 9*jup
          call MPI_RECV(temp1(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
                        0,commcol,&
                        status,ierr)
          jup = jup/9
        endif
        do i = 1, jup
          recv(:,i+nmv0) = temp1(:,i)
        enddo
        nmv0 = nmv0 + jup
        deallocate(temp1)
 
        deallocate(up)
        deallocate(down)
        deallocate(left)
        deallocate(right)

        if(numbuf.gt.480) then
          nsmall = numbuf/16
        else
          nsmall = numbuf
        endif
        allocate(left(9,nsmall))
        allocate(right(9,nsmall))
        allocate(up(9,nsmall))
        allocate(down(9,nsmall))
        allocate(mmsk(numbuf))
        mmsk = .true.
        ileft = 0
        iright = 0
        iup = 0
        idown = 0
        nmv0 = nmv0 - numbuf
        do i = 1, numbuf
          ii = i+nmv0
          if(recv(5,ii).le.lcrange(5)) then
            if(myidx.ne.0) then
              ileft = ileft + 1
              left(:,ileft) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(ileft,nsmall).eq.0) then
                allocate(temp1(9,ileft))
                do j = 1,ileft
                  temp1(:,j) = left(:,j)
                enddo
                deallocate(left)
                allocate(left(9,ileft+nsmall))
                do j = 1,ileft
                  left(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(5,ii).gt.lcrange(6)) then
            if(myidx.ne.(npx-1)) then
              iright = iright + 1
              right(:,iright) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iright,nsmall).eq.0) then
                allocate(temp1(9,iright))
                do j =1, iright
                  temp1(:,j) = right(:,j)
                enddo
                deallocate(right)
                allocate(right(9,iright+nsmall))
                do j =1, iright
                  right(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            else
              if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(iup,nsmall).eq.0) then
                  allocate(temp1(9,iup))
                  do j = 1, iup
                    temp1(:,j) = up(:,j)
                  enddo
                  deallocate(up)
                  allocate(up(9,iup+nsmall))
                  do j = 1, iup
                    up(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recv(:,ii)
                mmsk(i) = .false.
                if(mod(idown,nsmall).eq.0) then
                  allocate(temp1(9,idown))
                  do j = 1, idown
                    temp1(:,j) = down(:,j)
                  enddo
                  deallocate(down)
                  allocate(down(9,idown+nsmall))
                  do j = 1, idown
                    down(:,j) = temp1(:,j) 
                  enddo
                  deallocate(temp1)
                endif
                endif
              else
              endif
            endif
          else if(recv(3,ii).gt.lcrange(4)) then
            if(myidy.ne.(npy-1)) then
              iup = iup + 1
              up(:,iup) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(iup,nsmall).eq.0) then
                allocate(temp1(9,iup))
                do j = 1, iup
                  temp1(:,j) = up(:,j)
                enddo
                deallocate(up)
                allocate(up(9,iup+nsmall))
                do j = 1, iup
                  up(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else if(recv(3,ii).le.lcrange(3)) then
            if(myidy.ne.0) then
              idown = idown + 1
              down(:,idown) = recv(:,ii)
              mmsk(i) = .false.
              if(mod(idown,nsmall).eq.0) then
                allocate(temp1(9,idown))
                do j = 1, idown
                  temp1(:,j) = down(:,j)
                enddo
                deallocate(down)
                allocate(down(9,idown+nsmall))
                do j = 1, idown
                  down(:,j) = temp1(:,j) 
                enddo
                deallocate(temp1)
              endif
            endif
          else
          endif
        enddo
        nmv = ileft+iright+idown+iup
        call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
                           comm2d,ierr)
        if(totnmv.eq.0) then
          deallocate(mmsk)
          deallocate(up)
          deallocate(down)
          deallocate(left)
          deallocate(right)
          nmv0 = nmv0 + numbuf - nmv
          exit
        endif

        ic = 0
        allocate(temp1(9,nmv0+numbuf-nmv))
        do i = 1, nmv0
          temp1(:,i) = recv(:,i)
        enddo
        do i = 1, numbuf
          ii = i + nmv0
          if(mmsk(i)) then
            ic = ic + 1
            temp1(:,ic+nmv0) = recv(:,ii)
          endif
        enddo
        deallocate(mmsk)
        nmv0 = nmv0 + numbuf - nmv

!        print*," loop ", ij,numbuf,nmv,nmv0

        enddo

        !copy the remaining local particles into a temporary array.
        numpts = Nptlocal-nout
        allocate(temp1(9,numpts))
        ic = 0
        do i = 1, Nptlocal
          if(msk(i)) then
            ic = ic + 1
            temp1(:,ic) = Ptsl(:,i)
          endif
        enddo

        deallocate(msk)

!        call MPI_BARRIER(comm2d,ierr)
        !recopy the remaining local particles back to Ptsl which has 
        !a new size now.
        Nptlocal = numpts+nmv0 
        deallocate(Ptsl)
        allocate(Ptsl(9,Nptlocal))
        do i = 1, numpts
          Ptsl(:,i) = temp1(:,i)
        enddo
        deallocate(temp1)
        do i = 1, nmv0
          ii = i + numpts
          Ptsl(:,ii) = recv(:,i)
        enddo

        deallocate(recv)

        t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

        end subroutine ptsmv2_ptclmger
      end module Ptclmgerclass
