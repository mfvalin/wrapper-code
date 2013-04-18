#include <unistd.h>                                                                                                                                                                                                                                                     
#include <stdlib.h>                                                                                                                                                                                                                                                     
#include <stdio.h>                                                                                                                                                                                                                                                      
#include <mpi.h>                                                                                                                                                                                                                                                        
#include <mpi.h>
                                                                                                                                                                                                                                                                        
void main(int argc, char **argv)                                                                                                                                                                                                                                        
{                                                                                                                                                                                                                                                                       
 int my_rank=-1;                                                                                                     
 int my_size=-1;                                                                                                     
 char hostname[1204];                                                                                                                                                                                                                                                   
 char *MP_CHILD=getenv("MP_CHILD");                                                                                                                                                                                                                                 
 int buf[10];
 int buf2[10];
 MPI_Status status;
                                                                                                                                                                                                                                                                        
 if(MP_CHILD==NULL) MP_CHILD="-1";                                                                                                                                                                                                                                      
 gethostname(hostname, 1023);                                                                                                                                                                                                                                           
 MPI_Init(&argc,&argv);                                                                                                                                                                                                                                                 
 MPI_Comm_rank(MPI_COMM_WORLD , &my_rank);                                                                                                                                                                                                                              
 MPI_Comm_size(MPI_COMM_WORLD , &my_size);
 if(my_size != 2) goto the_end ;
MPI_Bcast(&buf,9,MPI_INTEGER,0,MPI_COMM_WORLD);
MPI_Reduce(&buf,&status,5,MPI_INTEGER,MPI_BOR,0,MPI_COMM_WORLD);
MPI_Allreduce(&buf,&buf2,7,MPI_INTEGER,MPI_BOR,MPI_COMM_WORLD);
MPI_Alltoall(&buf,3,MPI_INTEGER,&buf2,3,MPI_INTEGER,MPI_COMM_WORLD);
MPI_Alltoall(&buf,5,MPI_INTEGER,&buf2,5,MPI_INTEGER,MPI_COMM_WORLD);
 if(my_rank==0) {
   MPI_Send(&buf,10,MPI_INTEGER,1,0,MPI_COMM_WORLD);
   MPI_Send(&buf,8,MPI_INTEGER,1,0,MPI_COMM_WORLD);
   MPI_Recv(&buf,6,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
   MPI_Recv(&buf,4,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(&buf,9,MPI_INTEGER,1,0,&buf2,9,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(&buf,7,MPI_INTEGER,1,0,&buf2,7,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
 }else{
   MPI_Recv(&buf,10,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
   MPI_Recv(&buf,8,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
   MPI_Send(&buf,6,MPI_INTEGER,0,0,MPI_COMM_WORLD);
   MPI_Send(&buf,4,MPI_INTEGER,0,0,MPI_COMM_WORLD);
   MPI_Sendrecv(&buf,9,MPI_INTEGER,0,0,&buf2,9,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(&buf,7,MPI_INTEGER,0,0,&buf2,7,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
 }
MPI_Bcast(&buf,7,MPI_INTEGER,1,MPI_COMM_WORLD);
MPI_Reduce(&buf,&status,7,MPI_INTEGER,MPI_BOR,0,MPI_COMM_WORLD);
MPI_Allreduce(&buf,&buf2,9,MPI_INTEGER,MPI_BOR,MPI_COMM_WORLD);
the_end:
 printf("host = %s, rank = %d, MP_CHILD=%s\n",hostname,my_rank,MP_CHILD);
 MPI_Finalize();
}
