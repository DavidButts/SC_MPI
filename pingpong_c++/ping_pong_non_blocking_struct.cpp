#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

struct object{
  int a;
  double b[20];
  char c;
};

int main(void){

  //This code will write times out to a file of this
  // name.
  string fileName="Cpp_NoBlock_Struct.txt";

  // number of times to average over
  int commNum = 5000;

  int totalSize = pow(2,24);
  int rank;
  int size;
  real* bufferSend = NULL;
  real* bufferRecv = NULL;
  double start, end, totalTime;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  MPI_Request recvreq,sendreq;
  MPI_Status status, status2;

  //init datatype
  MPI_Datatype struc;
  int structlen = 3;
  int blocklengths[structlen];
  MPI_Datatype types[structlen];
  MPI_Aint displacements[structlen];
  //int a = 5;
  blocklengths[0] = 1;
  types[0] = MPI_INT;
  displacements[0] = offsetof(object,a);
  //double b[20];
  blocklengths[1] = 20;
  types[1] = MPI_DOUBLE;
  displacements[1] = offsetof(object,b);
  //char c = 'b';
  blocklengths[2] = 1;
  types[2] = MPI_CHAR;
  displacements[2] = offsetof(object,c);
  //create and commit
  MPI_Type_create_struct(structlen,blocklengths,displacements,types,&struc);
  MPI_Type_commit(&struc);


  //Allocate arrays only on two ranks to cut down on
  // memory size.
  // also, I allocate two arrays of max datasize (2^24 bytes)
  // and only send first n of those bytes from n =1 to 2^24.
  ofstream myFile;
  if(rank==0){
    myFile.open(fileName);
  }

  //Perform non-blocking and and recv calls
  for( int n = 1; n <= totalSize; n*=2){

    //ranks 1 and 0 allocate memory
    if(rank==0){
      bufferSend = (real* ) malloc(n*sizeof(object));
      bufferRecv = (real* ) malloc(n*sizeof(object));
    }
    if(rank==size-1){
      bufferSend = (real* ) malloc(n*sizeof(object));
      bufferRecv = (real* ) malloc(n*sizeof(object));
    }

    //loop that will perform non-block send/recv
    //only actual send/recv is timed.
    totalTime=0;
    for (int i =0; i < commNum; i++){
      start = MPI_Wtime();
      if(rank ==0){
        MPI_Isend(bufferSend,n,struc,size-1,1,MPI_COMM_WORLD,&sendreq);
        MPI_Irecv(bufferRecv,n,struc,size-1,0,MPI_COMM_WORLD,&recvreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
      if(rank == size-1){
        MPI_Irecv(bufferRecv,n,struc,0,1,MPI_COMM_WORLD,&recvreq);
        MPI_Isend(bufferSend,n,struc,0,0,MPI_COMM_WORLD,&sendreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
      end = MPI_Wtime();
      totalTime += end - start;
    }

    //have the two ranks that
    //initally allocated memory
    //free thier memory.
    if(rank==0){
      free(bufferSend);
      free(bufferRecv);
    }
    if(rank==size-1){
      free(bufferSend);
      free(bufferRecv);
    }

    //compute average for rank0/rank1
    double avgTime = (totalTime/((double)commNum));

    //rank1 send rank0 its average
    double rank1Avg = 0;
    if(rank == 0){
      MPI_Recv(&rank1Avg,1,MPI_DOUBLE,size-1,2,MPI_COMM_WORLD,&status2);
    }
    if (rank == size-1){
      MPI_Ssend(&avgTime,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    }

    //rank0 write time average of rank0/rank1 out to txt file.
    if(rank == 0){
      cout << "dataSize =  " << n << endl;
      myFile << n << " " << (totalTime/((double)commNum)) << " " << rank1Avg << endl;
    }
  }

  MPI_Type_free(&struc);

  //rank0 close the file it initially opened
  if(rank == 0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
