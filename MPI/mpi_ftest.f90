        program ftest
        use ISO_C_BINDING
        implicit none
        include 'mpif.h'
        integer error,noprocs,nid
        integer buf(10),status(10),buf2(10)
        real *8 :: buf8(4)

        call MPI_Init(error)
        call MPI_Comm_size(MPI_COMM_WORLD, noprocs, error)
        call MPI_Comm_rank(MPI_COMM_WORLD, nid, error)
        buf8=-1.0
        buf8(3)=654321.0
        buf8(4)=123456.0
        print *,'I am PE',nid+1,' of',noprocs
        if(noprocs /= 2) goto 777
        call mpi_bcast(buf,9,MPI_INTEGER,0,MPI_COMM_WORLD,error)
        call mpi_reduce(buf,status,5,MPI_INTEGER,MPI_BOR,0,MPI_COMM_WORLD,error)
        call mpi_allreduce(buf,status,7,MPI_INTEGER,MPI_BOR,MPI_COMM_WORLD,error)
        call mpi_alltoall(buf,3,MPI_INTEGER,buf2,3,MPI_INTEGER,MPI_COMM_WORLD,error)
        call mpi_alltoall(buf,5,MPI_INTEGER,buf2,5,MPI_INTEGER,MPI_COMM_WORLD,error)
        if(nid==0)then
          call mpi_send(buf,10,MPI_INTEGER,1,0,MPI_COMM_WORLD,error)
          call mpi_send(buf,8,MPI_INTEGER,1,0,MPI_COMM_WORLD,error)
          call mpi_recv(buf,6,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          call mpi_recv(buf,4,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          call mpi_sendrecv(buf,9,MPI_INTEGER,1,0,buf2,9,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          call mpi_sendrecv(buf,7,MPI_INTEGER,1,0,buf2,7,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          print *,'PE 0 sending to PE 1'
        else
          call mpi_recv(buf,10,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          call mpi_recv(buf,8,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          call mpi_send(buf,6,MPI_INTEGER,0,0,MPI_COMM_WORLD,error)
          call mpi_send(buf,4,MPI_INTEGER,0,0,MPI_COMM_WORLD,error)
          call mpi_sendrecv(buf,9,MPI_INTEGER,0,0,buf2,9,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          call mpi_sendrecv(buf,7,MPI_INTEGER,0,0,buf2,7,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          print *,'PE 1 receiving from PE 0'
        endif
        call mpi_bcast(buf,9,MPI_INTEGER,1,MPI_COMM_WORLD,error)
        call mpi_reduce(buf,status,7,MPI_INTEGER,MPI_BOR,1,MPI_COMM_WORLD,error)
        call mpi_allreduce(buf,status,9,MPI_INTEGER,MPI_BOR,MPI_COMM_WORLD,error)
        call mpi_allreduce(buf8(3),buf8(1),2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,error)
        print *,'MPI routines time, error=',buf8,error
777     call MPI_Finalize(error)
        stop
        end
 
