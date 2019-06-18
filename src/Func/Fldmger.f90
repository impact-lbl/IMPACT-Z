!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Fldmgerclass: Exchange guard grid data between neighboring processors 
! class in Communication module of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the functions to sum up the particle
!              contribution from neighboring processor domain, exchange
!              the potential, exchange the field for interpolation between
!              neighboring processors.
! Comments:
!----------------------------------------------------------------
        module Fldmgerclass
          use Timerclass
          use Pgrid2dclass
        contains
!----------------------------------------------------------------
! neighboring grid communication for the 3D open boundary conditions
        ! sum up the contributions from neighboring guard cells.
        subroutine guardsum1_Fldmger(x,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        double precision, dimension(innz)  :: tmp1,sendtmp
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npy,&
                   npx
        integer, dimension(2) :: tmpcoord
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        if(myidx.ne.(npx-1)) then
          right = myidx + 1
        else
          right = MPI_PROC_NULL
        endif
        if(myidx.ne.0) then
          left = myidx - 1
        else
          left = MPI_PROC_NULL
        endif 
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.ne.npy-1) then
          top = myidy + 1
        else
          top = MPI_PROC_NULL
        endif
        if(myidy.ne.0) then
          bottom = myidy -1
        else
          bottom = MPI_PROC_NULL
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        nsend = innz*innx
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,inny,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.0) then
          do k = 1, innz
            do i = 1, innx
              x(i,2,k) = x(i,2,k) + recvbuf1(i,k)
            enddo
          enddo
        endif
          
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                         msid,ierr)
          do k = 1, innz
            do i = 1, innx
              sendbuf1(i,k) = x(i,1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
        if(myidy.ne.(npy-1)) then
          do k = 1, innz
            do i = 1, innx
              x(i,inny-1,k) = x(i,inny-1,k) + recvbuf1(i,k)
            enddo
          enddo
        endif

        endif

        if(npx.gt.1) then
          nsend = innx*(inny-2*jadd)
          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,innz)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.0) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,2) = x(i,j,2) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif

          call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              sendbuf2(i,j-jadd) = x(i,j,1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow,&
                        ierr)
          call MPI_WAIT(msid,status,ierr)
          if(myidx.ne.(npx-1)) then
          do j = 1+jadd, inny-jadd
            do i = 1, innx
              x(i,j,innz-1) = x(i,j,innz-1) + recvbuf2(i,j-jadd)
            enddo
          enddo
          endif
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardsum = t_guardsum + elapsedtime_Timer(t0)

        end subroutine guardsum1_Fldmger

        ! exchange grid information between neighboring guard cell
        ! to calculate E from phi. 
        subroutine guardexch1_Fldmger(x,innx,inny,innz,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(innx,innz) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top
        integer :: leftx,rightx,bottomy,topy
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,npx,&
                   npy,msid
        integer, dimension(2) :: tmpcoord
        integer :: i,j,k,ierr
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
        
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        !define left,right,top,and bottom for cyclic shift purpose.
        if(myidx.eq.0) then
          left = npx - 1
          right = myidx + 1
        else if(myidx.eq.(npx-1)) then
          left = myidx - 1
          right = 0
        else
          left = myidx - 1
          right = myidx + 1
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        if(myidy.eq.0) then
          bottom = npy - 1
          top = myidy + 1
        else if(myidy.eq.(npy-1)) then
          bottom = myidy - 1
          top = 0
        else
          bottom = myidy - 1
          top = myidy + 1
        endif
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif

        if(npy.gt.1) then

        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,inny-1,k)
          enddo
        enddo
        nsend = (innz-2*kadd)*innx
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,1,k) = recvbuf1(i,k-kadd)
          enddo
        enddo
          
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            sendbuf1(i,k-kadd) = x(i,2,k)
          enddo
        enddo
        call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,&
                      commcol,msid,ierr)
        call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,&
                      commcol,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+kadd, innz-kadd
          do i = 1, innx
            x(i,inny,k) = recvbuf1(i,k-kadd)
          enddo
        enddo

        endif

        if(npx.gt.1) then

        nsend = innx*(inny-2*jadd)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            sendbuf2(j,k-jadd) = x(j,k,innz-1)
          enddo
        enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,3,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,3,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,1) = recvbuf2(j,k-jadd)
          enddo
        enddo
          
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            sendbuf2(j,k-jadd) = x(j,k,2)
          enddo
        enddo
        call MPI_IRECV(recvbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,right,4,&
                      commrow,msid,ierr)
        call MPI_SEND(sendbuf2(1,1),nsend,MPI_DOUBLE_PRECISION,left,4,&
                      commrow,ierr)
        call MPI_WAIT(msid,status,ierr)
        do k = 1+jadd, inny-jadd
          do j = 1, innx
            x(j,k,innz) = recvbuf2(j,k-jadd)
          enddo
        enddo
 
        endif

        call MPI_BARRIER(comm2d,ierr)
        t_guardexch = t_guardexch + elapsedtime_Timer(t0)

        end subroutine guardexch1_Fldmger

!----------------------------------------------------------------
        ! exchange the information for interpolation 
        ! from neighboring guard cells.
        subroutine boundint4_Fldmger(x1,x2,x3,innx,inny,innz,grid) 
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz
        double precision, dimension(innx,inny,innz), intent(inout) :: x1,x2,x3
        type (Pgrid2d), intent(in) :: grid
!        double precision, dimension(innx,innz-2,3) :: sendbuf1, recvbuf1
        double precision, allocatable, dimension(:,:,:) :: sendbuf1, recvbuf1
        double precision, dimension(innx,inny,3) :: sendbuf2, recvbuf2
        integer :: left,right,bottom,top,msid,&
                   ierr,i,j,k,nx,ny
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,totnp,&
                   npx,npy
        integer status(MPI_STATUS_SIZE)
        double precision :: t0
        integer :: nsend,jadd,kadd
    
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)
        if(npy.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        if(npx.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif

        top = myidy + 1
        bottom = myidy - 1
        left = myidx - 1
        right = myidx + 1

        ny = inny - 2*jadd

! This could be modified to improve speed
        if(npy.gt.1) then

        nsend = innx*(innz-2*kadd)*3
        allocate(recvbuf1(innx,innz-2*kadd,3))
        allocate(sendbuf1(innx,innz-2*kadd,3))

        if(myidy.ne.0) then
          call MPI_IRECV(recvbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.(npy-1)) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,1) = x1(i,inny-1,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,2) = x2(i,inny-1,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,3) = x3(i,inny-1,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
                        commcol,ierr)
        endif
        if(myidy.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x1(i,1,k) = recvbuf1(i,k-kadd,1)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x2(i,1,k) = recvbuf1(i,k-kadd,2)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x3(i,1,k) = recvbuf1(i,k-kadd,3)
            enddo
          enddo
        endif

        if(myidy.ne.(npy-1)) then
          call MPI_IRECV(recvbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
                        msid,ierr)
        endif
        if(myidy.ne.0) then
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,1) = x1(i,2,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,2) = x2(i,2,k)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              sendbuf1(i,k-kadd,3) = x3(i,2,k)
            enddo
          enddo
          call MPI_SEND(sendbuf1(1,1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
                        ierr)
        endif
        if(myidy.ne.(npy-1)) then
          call MPI_WAIT(msid,status,ierr)
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x1(i,inny,k) = recvbuf1(i,k-kadd,1)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x2(i,inny,k) = recvbuf1(i,k-kadd,2)
            enddo
          enddo
          do k = 1+kadd, innz-kadd
            do i = 1, innx
              x3(i,inny,k) = recvbuf1(i,k-kadd,3)
            enddo
          enddo
        endif

        deallocate(recvbuf1)
        deallocate(sendbuf1)

        endif

        if(npx.gt.1) then

        nsend = innx*inny*3
        if(myidx.ne.0) then
          call MPI_IRECV(recvbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,left,3,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.(npx-1)) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,1) = x1(k,j,innz-1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,2) = x2(k,j,innz-1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,3) = x3(k,j,innz-1)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,right,3,commrow,&
                        ierr)
        endif
        if(myidx.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x1(k,j,1) = recvbuf2(k,j,1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x2(k,j,1) = recvbuf2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x3(k,j,1) = recvbuf2(k,j,3)
            enddo
          enddo
        endif

        if(myidx.ne.(npx-1)) then
          call MPI_IRECV(recvbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,right,4,commrow,&
                         msid,ierr)
        endif
        if(myidx.ne.0) then
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,1) = x1(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,2) = x2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              sendbuf2(k,j,3) = x3(k,j,2)
            enddo
          enddo
          call MPI_SEND(sendbuf2(1,1,1),nsend,MPI_DOUBLE_PRECISION,left,4,commrow, &
                        ierr)
        endif
        if(myidx.ne.(npx-1)) then
          call MPI_WAIT(msid,status,ierr)
          do j = 1, inny
            do k = 1, innx
              x1(k,j,innz) = recvbuf2(k,j,1)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x2(k,j,innz) = recvbuf2(k,j,2)
            enddo
          enddo
          do j = 1, inny
            do k = 1, innx
              x3(k,j,innz) = recvbuf2(k,j,3)
            enddo
          enddo
        endif

        endif

        call MPI_BARRIER(comm2d,ierr)
        t_boundint = t_boundint + elapsedtime_Timer(t0)

        end subroutine boundint4_Fldmger

      end module Fldmgerclass
