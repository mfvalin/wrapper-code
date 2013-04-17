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
  long long bytes;          /* sum of bytes moved */
  double sbytes2;      /* sum of squares for bytes moved */
  double sbytes_coll;  /* sum of bytes * communicator_size for collectives */
  double sbytes_coll2; /* sum of squares for collectives */
  int calls;           /* number of calls to function*/
} stat_table_entry;

static stat_table_entry stat_table[]={
  {"MPI_Init" ,0,0.0,0.0,0.0,0},
  {"MPI_Recv" ,0,0.0,0.0,0.0,0},
  {"MPI_Irecv",0,0.0,0.0,0.0,0},
  {"MPI_Send" ,0,0.0,0.0,0.0,0},
  {"MPI_Isend",0,0.0,0.0,0.0,0},
  {NULL       ,0,0.0,0.0,0.0,0}
};

static int find_table_entry(const char *name){
  int i=0;
  while(stat_table[i].name != NULL){
    if(strcmp(name,stat_table[i].name)==0) return i;
    i++;
  }
  return -1;
}

static void add_to_entry(int me,int bytes,int commsize){
  float rbytes;
  if(me<0)return;
  stat_table[me].calls++;
  stat_table[me].bytes+=bytes;
  rbytes=bytes;
  stat_table[me].sbytes2+=(rbytes*rbytes);
  if(commsize>1){   /* for collectives only */
    rbytes=rbytes*commsize;
    stat_table[me].sbytes_coll+=rbytes;
    stat_table[me].sbytes_coll2+=(rbytes*rbytes);
  }
}

void dump_mpi_stats()
{
 int my_rank=-1;
 int i=0;
 double AVG=0.0;
 MPI_Comm_rank(MPI_COMM_WORLD , &my_rank);
 while(stat_table[i].name != NULL){
   if(stat_table[i].calls>0) {
     AVG=stat_table[i].bytes;
     AVG=AVG/stat_table[i].calls;
     printf("%4.4d: %-20s %6d messages %9Ld bytes (avg=%f)\n",my_rank,stat_table[i].name,stat_table[i].calls,stat_table[i].bytes,AVG);
   }
   i++;
 }
 return;
}
void dump_mpi_stats_(){ dump_mpi_stats();}
void dump_mpi_stats__(){ dump_mpi_stats();}

int MPI_Init(int *argc, char ***argv) {
  int rc;

  printf("INFO: entering profiling layer MPI_Init...\n");
  rc = PMPI_Init(argc, argv);
  return(rc);
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)
{
  int dsize;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Recv");
  PMPI_Type_size( datatype, &dsize );
  add_to_entry(me,count*dsize,1);
  return( PMPI_Recv(buf,count,datatype,source,tag,comm,status));
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
{
  int dsize;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Irecv");
  PMPI_Type_size( datatype, &dsize );
  add_to_entry(me,count*dsize,1);
  return( PMPI_Irecv(buf,count,datatype,source,tag,comm,request));
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
{
  int dsize;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Isend");
  PMPI_Type_size( datatype, &dsize );
  add_to_entry(me,count*dsize,1);
  return( PMPI_Isend(buf,count,datatype,dest,tag,comm,request));
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)
{
  int dsize;
  static int me=-1;
  if(me==-1) me=find_table_entry("MPI_Send");
  PMPI_Type_size( datatype, &dsize );
  add_to_entry(me,count*dsize,1);
  return( PMPI_Send(buf,count,datatype,dest,tag,comm));
}

int MPI_Finalize( void )
{
  dump_mpi_stats();
  printf("INFO: exiting from profiling layer MPI_Init...\n");
  return(PMPI_Finalize());
}

#ifdef TEST
                                                                                                                                                                                                                                                                        
void main(int argc, char **argv)                                                                                                                                                                                                                                        
{                                                                                                                                                                                                                                                                       
 int my_rank=-1;                                                                                                                                                                                                                                                        
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
 MPI_Comm_rank(MPI_COMM_WORLD , &my_rank);                                                                                                                                                                                                                              
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
 printf("host = %s, rank = %d, MP_CHILD=%s\n",hostname,my_rank,MP_CHILD);
 MPI_Finalize();
}

#endif
