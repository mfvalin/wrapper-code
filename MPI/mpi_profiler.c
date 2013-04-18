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
  {"MPI_Reduce"   ,0,0.0,0,0.0,0.0,0},
  {"MPI_Allreduce",0,0.0,0,0.0,0.0,0},
  {"MPI_Bcast"    ,0,0.0,0,0.0,0.0,0},
  {NULL           ,0,0.0,0,0.0,0.0,0}
};

static FILE *listfile=NULL;

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
  t0=MPI_Wtime();
  status=PMPI_Bcast(buffer,count,datatype,root,comm);
  add_to_entry(me,count*dsize,size,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status=PMPI_Reduce(sendbuf,recvbuf,count,datatype,op,root,comm);
  add_to_entry(me,count*dsize,size,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status=PMPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
  add_to_entry(me,count*dsize,size,MPI_Wtime()-t0);
  return status;
}

int MPI_Testall(int count,MPI_Request *array_of_requests,int *flag,MPI_Status *array_of_statuses){
  int dsize=1;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Testall");
  t0=MPI_Wtime();
  status=PMPI_Testall(count,array_of_requests,flag,array_of_statuses);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
  return status;
}

int MPI_Waitall(int count,MPI_Request *array_of_requests,MPI_Status *array_of_statuses){
  int dsize=1;
  int status;
  double t0;
  static int me=-1;
  
  if(me==-1) me=find_table_entry("MPI_Waitall");
  t0=MPI_Wtime();
  status=PMPI_Waitall(count,array_of_requests,array_of_statuses);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Recv(buf,count,datatype,source,tag,comm,Status);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Irecv(buf,count,datatype,source,tag,comm,request);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Isend(buf,count,datatype,dest,tag,comm,request);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Send(buf,count,datatype,dest,tag,comm);
  add_to_entry(me,count*dsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Sendrecv(sendbuf,sendcount,sendtype,dest,sendtag,recvbuf,recvcount,recvtype,source,recvtag,comm,Status);
  add_to_entry(me,sendcount*ssize+recvcount*rsize,1,MPI_Wtime()-t0);
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
  t0=MPI_Wtime();
  status = PMPI_Alltoall(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm);
  add_to_entry(me,sendcount*ssize+recvcount*rsize,size,MPI_Wtime()-t0);
  return status;
}

int MPI_Finalize( void )
{
  dump_mpi_stats();
  fprintf(listfile,"INFO: exiting from profiling layer MPI_Init...\n====================================================================\n");
  if(listfile!=stdout) fclose(listfile);
  return(PMPI_Finalize());
}

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
