      module mpistub
        interface MPI_ALLGATHERV
          module procedure MPI_ALLGATHERV1,MPI_ALLGATHERV2
        end interface
        interface MPI_ALLTOALLV
          module procedure MPI_ALLTOALLV1,MPI_ALLTOALLV2,MPI_ALLTOALLV3
        end interface
        interface MPI_ALLGATHER
          module procedure MPI_ALLGATHER1,MPI_ALLGATHER2,MPI_ALLGATHER3,&
                           MPI_ALLGATHER4,MPI_ALLGATHER5,MPI_ALLGATHER6
        end interface
        interface MPI_GATHER
          module procedure MPI_GATHER1,MPI_GATHER2,MPI_GATHER3,MPI_GATHER4
        end interface
        interface MPI_ISEND
          module procedure MPI_ISEND1,MPI_ISEND2,MPI_ISEND3,MPI_ISEND4
        end interface
        interface MPI_SEND
          module procedure MPI_SEND1,MPI_SEND2,MPI_SEND3,MPI_SEND4,&
                           MPI_SEND5,MPI_SEND6,MPI_SEND7
        end interface
        interface MPI_IRECV
          module procedure MPI_IRECV1,MPI_IRECV2,MPI_IRECV3,MPI_IRECV4,&
                           MPI_IRECV5,MPI_IRECV6
        end interface
        interface MPI_RECV
          module procedure MPI_RECV1,MPI_RECV2,MPI_RECV3,MPI_RECV4
        end interface
        interface MPI_REDUCE
          module procedure MPI_REDUCE1,MPI_REDUCE2,MPI_REDUCE3,MPI_REDUCE4
        end interface
        interface MPI_ALLREDUCE
          module procedure MPI_ALLREDUCE1,MPI_ALLREDUCE2,MPI_ALLREDUCE3,&
                           MPI_ALLREDUCE4,MPI_ALLREDUCE5,MPI_ALLREDUCE6
        end interface
        interface MPI_BCAST
          module procedure MPI_BCAST1,MPI_BCAST2,MPI_BCAST3,MPI_BCAST4,&
          MPI_BCAST5,MPI_BCAST6    
        end interface
      contains
        !mpi initialization
        subroutine MPI_INIT(ierr)
        end subroutine MPI_INIT

        !mpi initialization
        subroutine MPI_INITIALIZED(ierr,ierr2)
        end subroutine MPI_INITIALIZED

        !mpi end
        subroutine MPI_Finalize(ierr)
        end subroutine MPI_Finalize

        !global sum
        subroutine  MPI_ALLREDUCE1(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          double precision, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE1

        subroutine  MPI_ALLREDUCE2(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          double precision :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE2

        subroutine  MPI_ALLREDUCE3(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          integer, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE3

        subroutine  MPI_ALLREDUCE4(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          integer :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE4

        subroutine  MPI_ALLREDUCE5(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          double precision, dimension(:,:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE5

        subroutine  MPI_ALLREDUCE6(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,MPI_COMM_WORLD,ierr)
          double precision, dimension(:,:,:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLREDUCE6

        !synchronize communication 
        subroutine MPI_BARRIER(comm2d,ierr)
          integer :: comm2d
        end subroutine MPI_BARRIER

        !processor ID
        subroutine MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
          my_rank = 0
        end subroutine MPI_COMM_RANK

        !mpi timing
        double precision function MPI_WTIME()
          MPI_WTIME = 0.0
        end function MPI_WTIME

        !mpi broadcast
        subroutine MPI_BCAST1(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          double precision, dimension(:) :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST1

        subroutine MPI_BCAST2(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          double precision :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST2

        subroutine MPI_BCAST3(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          integer, dimension(:) :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST3

        subroutine MPI_BCAST4(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          integer :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST4

        subroutine MPI_BCAST5(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          double precision, dimension(:,:) :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST5

        subroutine MPI_BCAST6(rffile,num1,MPI_INTEGER,num2,comm2d,ierr)
          double precision, dimension(:,:,:) :: rffile
          integer :: comm2d
        end subroutine MPI_BCAST6

        !total number of processors
        subroutine MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
          np = 1
        end subroutine MPI_COMM_SIZE

        !sum to local processor
        subroutine  MPI_REDUCE1(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,num2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_REDUCE1
 
        subroutine  MPI_REDUCE2(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,num2,MPI_COMM_WORLD,ierr)
          double precision :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_REDUCE2
 
        subroutine  MPI_REDUCE3(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,num2,MPI_COMM_WORLD,ierr)
          integer, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_REDUCE3
 
        subroutine  MPI_REDUCE4(tmplc,tmpgl,num,MPI_DOUBLE_PRECISION,&
                                  MPI_SUM,num2,MPI_COMM_WORLD,ierr)
          integer :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_REDUCE4

        !mpi send command
        subroutine  MPI_SEND1(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:) :: tmplc
        end subroutine MPI_SEND1
 
        subroutine  MPI_SEND2(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          double precision :: tmplc
        end subroutine MPI_SEND2

        subroutine  MPI_SEND3(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          integer, dimension(:) :: tmplc
        end subroutine MPI_SEND3

        subroutine  MPI_SEND4(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          integer :: tmplc
        end subroutine MPI_SEND4

        subroutine  MPI_SEND5(tmplc,num,MPI_DOUBLE_COMPLEX,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          double complex, dimension(:,:) :: tmplc
        end subroutine MPI_SEND5

        subroutine  MPI_SEND6(tmplc,num,MPI_DOUBLE_COMPLEX,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          double complex, dimension(:,:,:) :: tmplc
        end subroutine MPI_SEND6

        subroutine  MPI_SEND7(tmplc,num,MPI_DOUBLE_COMPLEX,&
                                  num1,num2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:,:,:) :: tmplc
        end subroutine MPI_SEND7

        !mpi isend command
        subroutine  MPI_ISEND1(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,num3,ierr)
          double precision, dimension(:) :: tmplc
        end subroutine MPI_ISEND1
 
        subroutine  MPI_ISEND2(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,num3,ierr)
          double precision :: tmplc
        end subroutine MPI_ISEND2
 
        subroutine  MPI_ISEND3(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,num3,ierr)
          integer, dimension(:) :: tmplc
        end subroutine MPI_ISEND3
 
        subroutine  MPI_ISEND4(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,num3,ierr)
          integer :: tmplc
        end subroutine MPI_ISEND4

        !mpi wait command
        subroutine MPI_WAIT(num3,status,ierr)
          integer, dimension(:) :: status
        end subroutine MPI_WAIT

        !mpi wait all command
        subroutine MPI_WAITALL(num3,req,status,ierr)
          integer, dimension(:) :: req
          integer, dimension(:,:) :: status
        end subroutine MPI_WAITALL

        !mpi recv command
        subroutine  MPI_RECV1(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,status,ierr)
          double precision :: tmplc
          integer, dimension(:) :: status
        end subroutine MPI_RECV1
 
        subroutine  MPI_RECV2(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,status,ierr)
          integer :: tmplc
          integer, dimension(:) :: status
        end subroutine MPI_RECV2

        subroutine  MPI_RECV3(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,status,ierr)
          double precision, dimension(:) :: tmplc
          integer, dimension(:) :: status
        end subroutine MPI_RECV3
 
        subroutine  MPI_RECV4(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,status,ierr)
          integer, dimension(:) :: tmplc
          integer, dimension(:) :: status
        end subroutine MPI_RECV4
 
        !mpi irecv command
        subroutine  MPI_IRECV1(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          double precision :: tmplc
        end subroutine MPI_IRECV1
 
        subroutine  MPI_IRECV2(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          integer :: tmplc
        end subroutine MPI_IRECV2

        subroutine  MPI_IRECV3(tmplc,num,MPI_DOUBLE_COMPLEX,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          double complex, dimension(:,:) :: tmplc
        end subroutine MPI_IRECV3
 
        subroutine  MPI_IRECV4(tmplc,num,MPI_DOUBLE_COMPLEX,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          double complex, dimension(:,:,:) :: tmplc
        end subroutine MPI_IRECV4

        subroutine  MPI_IRECV5(tmplc,num,MPI_DOUBLE_PRECISION,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          double precision, dimension(:,:,:) :: tmplc
        end subroutine MPI_IRECV5
 
        subroutine  MPI_IRECV6(tmplc,num,MPI_INTEGER,&
                                  num1,num2,MPI_COMM_WORLD,msid,ierr)
          integer, dimension(:) :: tmplc
        end subroutine MPI_IRECV6

        !mpi gather command
        subroutine  MPI_GATHER1(tmplc,num,MPI_DOUBLE_PRECISION,&
                               tmpgl,num1,MPI_DOUBLE_PRECISION2,num2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_GATHER1

        subroutine  MPI_GATHER2(tmplc,num,MPI_DOUBLE_PRECISION,&
                               tmpgl,num1,MPI_DOUBLE_PRECISION2,num2,MPI_COMM_WORLD,ierr)
          double precision :: tmplc
          double precision, dimension(:) :: tmpgl
          tmpgl = tmplc
        end subroutine MPI_GATHER2
 
        subroutine  MPI_GATHER3(tmplc,num,MPI_INTEGER,&
                               tmpgl,num1,MPI_INTEGER2,num2,MPI_COMM_WORLD,ierr)
          integer, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_GATHER3

        subroutine  MPI_GATHER4(tmplc,num,MPI_INTEGER,&
                               tmpgl,num1,MPI_INTEGER2,num2,MPI_COMM_WORLD,ierr)
          integer :: tmplc
          integer, dimension(:) :: tmpgl
          tmpgl = tmplc
        end subroutine MPI_GATHER4

        !mpi allgather command
        subroutine  MPI_ALLGATHER1(tmplc,num,MPI_DOUBLE_PRECISION,&
                               tmpgl,num1,MPI_DOUBLE_PRECISION2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER1

        subroutine  MPI_ALLGATHER2(tmplc,num,MPI_DOUBLE_PRECISION,&
                               tmpgl,num1,MPI_DOUBLE_PRECISION2,MPI_COMM_WORLD,ierr)
          double precision :: tmplc
          double precision, dimension(:) :: tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER2
 
        subroutine  MPI_ALLGATHER3(tmplc,num,MPI_INTEGER,&
                               tmpgl,num1,MPI_INTEGER2,MPI_COMM_WORLD,ierr)
          integer, dimension(:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER3

        subroutine  MPI_ALLGATHER4(tmplc,num,MPI_INTEGER,&
                               tmpgl,num1,MPI_INTEGER2,MPI_COMM_WORLD,ierr)
          integer :: tmplc
          integer, dimension(:) :: tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER4

        subroutine  MPI_ALLGATHER5(tmplc,num,MPI_DOUBLE_PRECISION,&
                               tmpgl,num1,MPI_DOUBLE_PRECISION2,MPI_COMM_WORLD,ierr)
          double precision, dimension(:,:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER5

        subroutine  MPI_ALLGATHER6(tmplc,num,MPI_INTEGER,&
                               tmpgl,num1,MPI_INTEGER2,MPI_COMM_WORLD,ierr)
          integer, dimension(:,:) :: tmplc,tmpgl
          tmpgl = tmplc
        end subroutine MPI_ALLGATHER6

        subroutine  MPI_CART_CREATE(comm,num1,dims,period,tt, &
                             comm_2d,ierr)
        integer :: comm,comm_2d
        integer, dimension(:) :: dims
        logical, dimension(:) :: period
        logical :: tt

        comm_2d = 1

        end subroutine MPI_CART_CREATE

        subroutine MPI_CART_COORDS(comm_2d,myrank,num,local,ierr)
          integer :: comm_2d
          integer, dimension(:) :: local

          local = 0
        end subroutine MPI_CART_COORDS

        subroutine MPI_CART_SUB(comm_2d,remaindims,col_comm,ierr)
          integer :: comm_2d,col_comm
          logical, dimension(:) :: remaindims

          col_comm = 1
        end subroutine MPI_CART_SUB

        !mpi alltoallv1 command
        subroutine MPI_ALLTOALLV1(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX2,comm,ierr)
          double complex, dimension(:) :: sendbuf,recvbuf 
          integer, dimension(:) :: sendcount,recvcount,senddisp,recvdisp
          integer :: comm
        end subroutine MPI_ALLTOALLV1

        subroutine MPI_ALLTOALLV2(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX2,comm,ierr)
          double precision, dimension(:) :: sendbuf,recvbuf 
          integer, dimension(:) :: sendcount,recvcount,senddisp,recvdisp
          integer :: comm
        end subroutine MPI_ALLTOALLV2

        subroutine MPI_ALLTOALLV3(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
           recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX2,comm,ierr)
          integer, dimension(:) :: sendbuf,recvbuf 
          integer, dimension(:) :: sendcount,recvcount,senddisp,recvdisp
          integer :: comm
        end subroutine MPI_ALLTOALLV3

        !mpi allgatherv command
        subroutine MPI_ALLGATHERV1(rhoz,innz,MPI_DOUBLE_PRECISION,recvrhoz,&
                            ztable,zdisp,MPI_DOUBLE_PRECISION2,commrow,ierr)
          double precision, dimension(:) :: rhoz,recvrhoz
          integer, dimension(:) :: ztable,zdisp
          integer :: commrow

          recvrhoz = rhoz

        end subroutine MPI_ALLGATHERV1

        subroutine MPI_ALLGATHERV2(rhoz,innz,MPI_INTEGER,recvrhoz,&
                            ztable,zdisp,MPI_INTEGER2,commrow,ierr)
          integer, dimension(:) :: rhoz,recvrhoz
          integer, dimension(:) :: ztable,zdisp
          integer :: commrow

          recvrhoz = rhoz

        end subroutine MPI_ALLGATHERV2

      end module mpistub
