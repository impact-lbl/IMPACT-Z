!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Pgrid2dclass: logical 2D Cartesian processor grid class in DATA 
!               STRUCTURE layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class construct a logical 2-D Cartesian processor grid.
! Comments:
!----------------------------------------------------------------
      module Pgrid2dclass
        use mpistub
        type Pgrid2d
          private
          integer comm_2d      ! comunicator for entire grid
          integer col_comm     ! comunicator for my col
          integer row_comm     ! comunicator for my row
          integer np           ! total number of processors
          integer npcol        ! number of processors in my col
          integer nprow        ! number of processors in my row
          integer myrank       ! my rank in the grid comm
          integer my_col       ! my column number
          integer my_row       ! my row number
        end type Pgrid2d
      contains
        ! set up 2-D Cartisian grid.
        subroutine construct_Pgrid2d(this,comm,nrow,ncol)
        implicit none
        include 'mpif.h'
        type (Pgrid2d), intent(out) :: this
        integer, intent(in) :: comm,nrow,ncol
        integer, dimension(2) :: dims,local
        integer :: ierr
        logical, dimension(2) :: period,remaindims
       
        ! set up global grid information.
        call MPI_COMM_SIZE(comm,this%np,ierr)

        ! set up the size of 2-D logical processor grid.
        this%npcol = ncol
        this%nprow = nrow

        ! create 2-D grid communicator.
        dims(1) = ncol
        dims(2) = nrow
        period(1) = .false.
        period(2) = .false.

        call MPI_CART_CREATE(comm,2,dims,period,.true., &
                             this%comm_2d,ierr)

        ! get my rank, my column number and my row number in grid.
        call MPI_COMM_RANK(this%comm_2d, this%myrank, ierr)
        call MPI_CART_COORDS(this%comm_2d,this%myrank,2,local,ierr)
        this%my_col = local(1)
        this%my_row = local(2)

        ! set up column communicators.
        remaindims(1) = .true.
        remaindims(2) = .false.
        call MPI_CART_SUB(this%comm_2d,remaindims,this%col_comm,ierr)

        ! set up row communciators.
        remaindims(1) = .false.
        remaindims(2) = .true.
        call MPI_CART_SUB(this%comm_2d,remaindims,this%row_comm,ierr)

        end subroutine construct_Pgrid2d

        ! get 2D, column, and row communicators.
        subroutine getcomm_Pgrid2d(this,comm2d,commcol,commrow)
        implicit none
        type (Pgrid2d),intent(in) :: this
        integer,intent(out) :: comm2d,commcol,commrow

        comm2d = this%comm_2d
        commcol = this%col_comm
        commrow = this%row_comm

        end subroutine getcomm_Pgrid2d

        !get my rank position in 2d, column, and row communicators.
        subroutine getpost_Pgrid2d(this,mypos,mycolpos,myrowpos)
        implicit none
        type (Pgrid2d),intent(in) :: this
        integer,intent(out) :: mypos,mycolpos,myrowpos

        mypos = this%myrank
        mycolpos = this%my_col
        myrowpos = this%my_row

        end subroutine getpost_Pgrid2d

        !get the total number of PEs, PEs in column and row communicators.
        subroutine getsize_Pgrid2d(this,totsize,colsize,rowsize)
        implicit none
        type (Pgrid2d),intent(in) :: this
        integer,intent(out) :: totsize,colsize,rowsize

        totsize = this%np
        colsize = this%npcol
        rowsize = this%nprow
        
        end subroutine getsize_Pgrid2d

      end module Pgrid2dclass
