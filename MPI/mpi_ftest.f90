        program ftest
        implicit none
        include 'mpif.h'
        integer error,noprocs,nid
        integer buf(10),status(10)

        call MPI_Init(error)
        call MPI_Comm_size(MPI_COMM_WORLD, noprocs, error)
        call MPI_Comm_rank(MPI_COMM_WORLD, nid, error)
        print *,'I am PE',nid+1,' of',noprocs
        call mpi_bcast(buf,9,MPI_INTEGER,0,MPI_COMM_WORLD,error)
        if(nid==0)then
          call mpi_send(buf,10,MPI_INTEGER,1,0,MPI_COMM_WORLD,error)
          call mpi_send(buf,8,MPI_INTEGER,1,0,MPI_COMM_WORLD,error)
          call mpi_recv(buf,6,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          call mpi_recv(buf,4,MPI_INTEGER,1,0,MPI_COMM_WORLD,status,error)
          print *,'PE 0 sending to PE 1'
        else
          call mpi_recv(buf,10,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          call mpi_recv(buf,8,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,error)
          call mpi_send(buf,6,MPI_INTEGER,0,0,MPI_COMM_WORLD,error)
          call mpi_send(buf,4,MPI_INTEGER,0,0,MPI_COMM_WORLD,error)
          print *,'PE 1 receiving from PE 0'
        endif
        call mpi_bcast(buf,9,MPI_INTEGER,0,MPI_COMM_WORLD,error)
!        call dump_mpi_stats
        call MPI_Finalize(error)
        stop
        end
 
