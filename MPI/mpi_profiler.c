#include <unistd.h>                                                                                                                                                                                                                                                     
#include <stdlib.h>                                                                                                                                                                                                                                                     
#include <stdio.h>                                                                                                                                                                                                                                                      
#include <string.h>                                                                                                                                                                                                                                                        
#include <mpi.h>

static int MPI_Irecv_ctr=0;
static int MPI_Recv_ctr=0;
static int MPI_Isend_ctr=0;
static int MPI_Send_ctr=0;
static int MPI_Irecv_elem=0;
static int MPI_Recv_elem=0;
static int MPI_Isend_elem=0;
static int MPI_Send_elem=0;

struct {
  double time;
  int calls;
} pmpi_r_statistics ={0.0,0};

typedef struct {
  char *name;          /* name of tracked function */
  long long sbytes;          /* sum of bytes moved */
  double sbytes2;      /* sum of squares for bytes moved */
  long long sbytes_coll;  /* sum of bytes * communicator_size for collectives */
  double sbytes_coll2; /* sum of squares for collectives */
  double time;         /* time spent in routine */
  int calls;           /* number of calls to function*/
} stat_table_entry;

static stat_table_entry stat_table[]={
  {"MPI_Recv"     ,0,0.0,0,0.0,0.0,0},
  {"MPI_Irecv"    ,0,0.0,0,0.0,0.0,0},
  {"MPI_Send"     ,0,0.0,0,0.0,0.0,0},
  {"MPI_Sendrecv" ,0,0.0,0,0.0,0.0,0},
  {"MPI_Isend"    ,0,0.0,0,0.0,0.0,0},
  {"MPI_Waitall"  ,0,0.0,0,0.0,0.0,0},
  {"MPI_Testall"  ,0,0.0,0,0.0,0.0,0},
  {"MPI_Alltoall" ,0,0.0,0,0.0,0.0,0},
  {"MPI_Alltoallv",0,0.0,0,0.0,0.0,0},
  {"MPI_Reduce"   ,0,0.0,0,0.0,0.0,0},
  {"MPI_Allreduce",0,0.0,0,0.0,0.0,0},
  {"MPI_Gather"   ,0,0.0,0,0.0,0.0,0},
  {"MPI_Gatherv"  ,0,0.0,0,0.0,0.0,0},
  {"MPI_Scatter"  ,0,0.0,0,0.0,0.0,0},
  {"MPI_Scatterv" ,0,0.0,0,0.0,0.0,0},
  {"MPI_Bcast"    ,0,0.0,0,0.0,0.0,0},
  {"MPI_Barrier"  ,0,0.0,0,0.0,0.0,0},
  {NULL           ,0,0.0,0,0.0,0.0,0}
};

static FILE *listfile=NULL;

static int buffer_sizes[33];
static int root_buffer_sizes[33];

static void buf_size_stat(int bsize)
{
  int nbits=0;
  if(bsize <= 0) return;
  while(bsize & ~0xFF) { nbits+=8 ; bsize >>=8; }
  while(bsize & ~0x07) { nbits+=3 ; bsize >>=3; }
  while(bsize & ~0x01) { nbits+=1 ; bsize >>=1; }
  if(nbits>32) return;
  buffer_sizes[nbits]++;
}

static void print_buf_size_stat(char *mesg,int *buffer_sizes)
{
  int i,start;
  start=1;
  for (i=0 ; i<30 ; i++) {
    if(buffer_sizes[i] > 0){
      fprintf(listfile,"%s%6d buffers of size %9d -> %9d\n",mesg,buffer_sizes[i],start,start*2-1);
    }
    start *=2;
  }
}

static int find_table_entry(const char *name){
  int i=0;
  while(stat_table[i].name != NULL){
    if(strcmp(name,stat_table[i].name)==0) return i;
    i++;
  }
  return -1;
}

static void add_to_entry(int me,int bytes,int commsize,double time){
  float rbytes;
  if(me<0)return;
  pmpi_r_statistics.time+=time;
  pmpi_r_statistics.calls++;
  stat_table[me].calls++;
  stat_table[me].sbytes+=bytes;
  stat_table[me].time+=time;
  rbytes=bytes;
  stat_table[me].sbytes2+=(rbytes*rbytes);
  if(commsize>1){   /* for collectives only */
    rbytes=rbytes*commsize;
    stat_table[me].sbytes_coll+=rbytes;
    stat_table[me].sbytes_coll2+=(rbytes*rbytes);
  }
}

void reset_mpi_stats()
{
 int i=0;
 while(stat_table[i].name != NULL){
   stat_table[i].sbytes=0;
   stat_table[i].calls=0;
   stat_table[i].time=0;
   stat_table[i].sbytes2=0;
   stat_table[i].sbytes_coll=0;
   stat_table[i].sbytes_coll2=0;
   i++;
 }
}
void reset_mpi_stats_(){ reset_mpi_stats();}
void reset_mpi_stats__(){ reset_mpi_stats();}

void dump_mpi_stats()
{
 int my_rank=-1;
 int size=1;
 int i=0;
 double AVG=0.0;
 double tmax,tmin,tmean;
 int tcalls;
 long long tbytes;
 
 MPI_Comm_rank(MPI_COMM_WORLD , &my_rank);
 MPI_Comm_size(MPI_COMM_WORLD , &size);
 while(stat_table[i].name != NULL){
   if(stat_table[i].calls>0) {
     AVG=stat_table[i].sbytes;
     AVG=AVG/stat_table[i].calls;
     fprintf(listfile,"%5.5d: %-20s %6d messages %9Ld [%9Ld] bytes (avg=%f), %12.6f seconds\n",
	    my_rank,stat_table[i].name,stat_table[i].calls,stat_table[i].sbytes,
	    stat_table[i].sbytes_coll,AVG,stat_table[i].time);
   }
   i++;
 }
 i=0;
 while(stat_table[i].name != NULL){
     PMPI_Allreduce(&stat_table[i].calls,&tcalls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
     PMPI_Allreduce(&stat_table[i].sbytes,&tbytes,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD);
     PMPI_Allreduce(&stat_table[i].time,&tmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD);
     PMPI_Allreduce(&stat_table[i].time,&tmin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD);
     PMPI_Allreduce(&stat_table[i].time,&tmean,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
     if(my_rank==0 && tcalls>0) {
       fprintf(listfile,"TOTAL: %-20s %6d messages %9Ld bytes, min/max/avg= %f/%f/%f seconds\n",
	      stat_table[i].name,tcalls,tbytes,tmin,tmax,tmean/size);
     }
   i++;
 }
 print_buf_size_stat("local: ",buffer_sizes);
 PMPI_Allreduce(buffer_sizes,root_buffer_sizes,33,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
 if(my_rank==0 )print_buf_size_stat("TOTAL: ",root_buffer_sizes);
 return;
}
void dump_mpi_stats_(){ dump_mpi_stats();}
void dump_mpi_stats__(){ dump_mpi_stats();}

static int my_rank=-1;

int MPI_Init(int *argc, char ***argv) {
  int rc;
  char fname[4096];
  char *Fname = &fname[0] ;
  char *envfile;
  char *mode="w";

  rc = PMPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD , &my_rank);
  envfile=getenv("PMPI_OUT_FILE");
  if(envfile != NULL) {
    sprintf(Fname,"%s_%5.5d",envfile,my_rank);
    if( *Fname == '+' ) { Fname++ ; mode="a"; }
    listfile=fopen(Fname,mode);
  }
  if(listfile == NULL ) listfile=stdout;
  if(my_rank==0) printf("INFO: entering profiling layer MPI_Init...\n");
  return(rc);
}

int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype,int root, MPI_Comm comm){
  int dsize, size;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Bcast");
  PMPI_Type_size( datatype, &dsize );
  MPI_Comm_size(comm,&size);
  t0=PMPI_Wtime();
  status=PMPI_Bcast(buffer,count,datatype,root,comm);
  add_to_entry(me,count*dsize,size,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Reduce(void *sendbuf, void *recvbuf, int count,
            MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
  int dsize, size;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Reduce");
  PMPI_Type_size( datatype, &dsize );
  MPI_Comm_size(comm,&size);
  t0=PMPI_Wtime();
  status=PMPI_Reduce(sendbuf,recvbuf,count,datatype,op,root,comm);
  add_to_entry(me,count*dsize,size,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int dsize, size;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Allreduce");
  PMPI_Type_size( datatype, &dsize );
  MPI_Comm_size(comm,&size);
  t0=PMPI_Wtime();
  status=PMPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  add_to_entry(me,count*dsize,size,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Barrier( MPI_Comm comm )
{
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Barrier");
  t0=PMPI_Wtime();
  status=PMPI_Barrier( comm );
  add_to_entry(me,0,1,PMPI_Wtime()-t0);
  return status;
}

int MPI_Testall(int count,MPI_Request *array_of_requests,int *flag,MPI_Status *array_of_statuses){
  int dsize=1;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Testall");
  t0=PMPI_Wtime();
  status=PMPI_Testall(count,array_of_requests,flag,array_of_statuses);
  add_to_entry(me,0,1,PMPI_Wtime()-t0);
  return status;
}

int MPI_Waitall(int count,MPI_Request *array_of_requests,MPI_Status *array_of_statuses){
  int dsize=1;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Waitall");
  t0=PMPI_Wtime();
  status=PMPI_Waitall(count,array_of_requests,array_of_statuses);
  add_to_entry(me,0,1,PMPI_Wtime()-t0);
  return status;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *Status)
{
  int dsize;
  int status;
  double t0;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Recv");
  PMPI_Type_size( datatype, &dsize );
  t0=PMPI_Wtime();
  status = PMPI_Recv(buf,count,datatype,source,tag,comm,Status);
  add_to_entry(me,count*dsize,1,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
{
  int dsize;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Irecv");
  PMPI_Type_size( datatype, &dsize );
  t0=PMPI_Wtime();
  status = PMPI_Irecv(buf,count,datatype,source,tag,comm,request);
  add_to_entry(me,count*dsize,1,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
{
  int dsize;
  int status;
  double t0;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Isend");
  PMPI_Type_size( datatype, &dsize );
  t0=PMPI_Wtime();
  status = PMPI_Isend(buf,count,datatype,dest,tag,comm,request);
  add_to_entry(me,count*dsize,1,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)
{
  int dsize;
  int status;
  double t0;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Send");
  PMPI_Type_size( datatype, &dsize );
  t0=PMPI_Wtime();
  status = PMPI_Send(buf,count,datatype,dest,tag,comm);
  add_to_entry(me,count*dsize,1,PMPI_Wtime()-t0);
  buf_size_stat(count*dsize);
  return status;
}

int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
            int dest, int sendtag, void *recvbuf, int recvcount,
            MPI_Datatype recvtype, int source, int recvtag,
            MPI_Comm comm, MPI_Status *Status)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Sendrecv");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Sendrecv(sendbuf,sendcount,sendtype,dest,sendtag,recvbuf,recvcount,recvtype,source,recvtag,comm,Status);
  add_to_entry(me,sendcount*ssize+recvcount*rsize,1,PMPI_Wtime()-t0);
  buf_size_stat(sendcount*ssize);
  buf_size_stat(recvcount*rsize);
  return status;
}

int MPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  int rank;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if(me==-1) me=find_table_entry("MPI_Scatter");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Scatter(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  if(rank!=root){
    add_to_entry(me,recvcnt*rsize,1,PMPI_Wtime()-t0);
    buf_size_stat(recvcnt*rsize);
  }else{
    add_to_entry(me,sendcnt*ssize*size,1,PMPI_Wtime()-t0);
    while(size--) buf_size_stat(sendcnt*ssize);
  }
  return status;
}

int MPI_Scatterv( void *sendbuf, int *sendcnts, int *displs, 
                 MPI_Datatype sendtype, void *recvbuf, int recvcnt,
                 MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  int rank;
  int sendcnt, i;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if(me==-1) me=find_table_entry("MPI_Scatterv");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Scatterv(sendbuf,sendcnts,displs,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  if(rank==root) {
    sendcnt=0;
    for (i=0;i<size;i++) sendcnt+=sendcnts[i];
    add_to_entry(me,sendcnt*ssize,1,PMPI_Wtime()-t0);
    while(size--) buf_size_stat(sendcnts[i]*rsize);
  }else{
    add_to_entry(me,recvcnt*rsize,1,PMPI_Wtime()-t0);
    buf_size_stat(recvcnt*rsize);
  }
  return status;
}

int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  int rank;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if(me==-1) me=find_table_entry("MPI_Gather");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);
  if(rank==root) {
    add_to_entry(me,recvcnt*rsize*size,1,PMPI_Wtime()-t0);
    while(size--) buf_size_stat(recvcnt*rsize);
  }else{
    add_to_entry(me,sendcnt*ssize,1,PMPI_Wtime()-t0);
    buf_size_stat(sendcnt*ssize);
  }
  return status;
}

int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  int rank;
  int recvcnt, i;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  if(me==-1) me=find_table_entry("MPI_Gatherv");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Gatherv(sendbuf,sendcnt,sendtype,recvbuf,recvcnts,displs,recvtype,root,comm);
  if(rank==root) {
    recvcnt=0;
    for (i=0;i<size;i++) recvcnt+=recvcnts[i];
    add_to_entry(me,recvcnt*rsize,1,PMPI_Wtime()-t0);
    while(size--) buf_size_stat(recvcnts[i]*rsize);
  }else{
    add_to_entry(me,sendcnt*ssize,1,PMPI_Wtime()-t0);
    buf_size_stat(sendcnt*ssize);
  }
  return status;
}

int MPI_Alltoall(void *sendbuf, int sendcount,
            MPI_Datatype sendtype, void *recvbuf, int recvcount,
            MPI_Datatype recvtype, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  
  MPI_Comm_size(comm,&size);
  if(me==-1) me=find_table_entry("MPI_Alltoall");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  t0=PMPI_Wtime();
  status = PMPI_Alltoall(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm);
  add_to_entry(me,sendcount*ssize+recvcount*rsize,size,PMPI_Wtime()-t0);
  while(size--) { buf_size_stat(sendcount*ssize); buf_size_stat(recvcount*rsize); }
  return status;
}

int MPI_Alltoallv(void *sendbuf, int *sendcnts, int *sdispls, 
                  MPI_Datatype sendtype, void *recvbuf, int *recvcnts, 
                  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
{
  int rsize,ssize;
  int status;
  double t0;
  static int me=-1;
  int size;
  int sendcnt,recvcnt,i;
  
  MPI_Comm_size(comm,&size);
  if(me==-1) me=find_table_entry("MPI_Alltoallv");
  PMPI_Type_size( sendtype, &ssize );
  PMPI_Type_size( recvtype, &rsize );
  recvcnt=0;
  for (i=0;i<size;i++) { recvcnt+=recvcnts[i]; buf_size_stat(recvcnts[i]*rsize) ; }
  sendcnt=0;
  for (i=0;i<size;i++) { sendcnt+=sendcnts[i]; buf_size_stat(sendcnts[i]*ssize) ; }
  t0=PMPI_Wtime();
  status = PMPI_Alltoallv(sendbuf,sendcnts,sdispls,sendtype,recvbuf,recvcnts,rdispls,recvtype,comm);
  add_to_entry(me,sendcnt*ssize+recvcnt*rsize,1,PMPI_Wtime()-t0);
  return status;
}

int MPI_Finalize( void )
{
  dump_mpi_stats();
  fprintf(listfile,"INFO: exiting from profiling layer MPI_Init...\n====================================================================\n");
  if(listfile!=stdout) fclose(listfile);
  return(PMPI_Finalize());
}

void mpi_finalize(int *error){ *error=MPI_Finalize() ;}
void mpi_finalize_(int *error){ *error=MPI_Finalize() ;}
void mpi_finalize__(int *error){ *error=MPI_Finalize() ;}

#ifdef TEST
                                                                                                                                                                                                                                                                        
void main(int argc, char **argv)                                                                                                                                                                                                                                        
{                                                                                                                                                                                                                                                                       
 char hostname[1204];                                                                                                                                                                                                                                                   
 char *MP_CHILD=getenv("MP_CHILD");                                                                                                                                                                                                                                     
 int buf[10];
 MPI_Status status;
 int position=-2;
 
 position=find_table_entry("MPI_Irecv");
 printf("entry %s found at position %d\n","MPI_Irecv",position);
                                                                                                                                                                                                                                                                        
 if(MP_CHILD==NULL) MP_CHILD="-1";                                                                                                                                                                                                                                      
 gethostname(hostname, 1023);                                                                                                                                                                                                                                           
 MPI_Init(&argc,&argv);    
 MPI_Bcast(&buf,9,MPI_INTEGER,0,MPI_COMM_WORLD);
 if(my_rank==0) {
   MPI_Send(&buf,10,MPI_INTEGER,1,0,MPI_COMM_WORLD);
   MPI_Send(&buf,8,MPI_INTEGER,1,0,MPI_COMM_WORLD);
   MPI_Recv(&buf,6,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
   MPI_Recv(&buf,4,MPI_INTEGER,1,0,MPI_COMM_WORLD,&status);
 }else{
   MPI_Recv(&buf,10,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
   MPI_Recv(&buf,8,MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
   MPI_Send(&buf,6,MPI_INTEGER,0,0,MPI_COMM_WORLD);
   MPI_Send(&buf,4,MPI_INTEGER,0,0,MPI_COMM_WORLD);
 }
 MPI_Bcast(&buf,7,MPI_INTEGER,0,MPI_COMM_WORLD);
 printf("host = %s, rank = %d, MP_CHILD=%s\n",hostname,my_rank,MP_CHILD);
 MPI_Finalize();
}

#endif
