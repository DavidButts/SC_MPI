#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

int main(void){

  //This code will write times out to a file of this
  // name.
  string fileName="Cpp_NoBlock_Double.txt";

  // number of times to average over
  int commNum = 5000;

  int totalSize = pow(2,4);
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

    //Allocate memory and first/last ranks
    if(rank==0){
      bufferSend = (real* ) malloc(n*sizeof(real));
      bufferRecv = (real* ) malloc(n*sizeof(real));
    }
    if(rank==size-1){
      bufferSend = (real* ) malloc(n*sizeof(real));
      bufferRecv = (real* ) malloc(n*sizeof(real));
    }

    //loop that will perform non-block send/recv
    //only actual send/recv is timed.
    totalTime=0;
    for (int i =0; i < commNum; i++){
      start = MPI_Wtime();
      if(rank ==0){
        MPI_Isend(bufferSend,n,MPI_DOUBLE,size-1,0,MPI_COMM_WORLD,&sendreq);
        MPI_Irecv(bufferRecv,n,MPI_DOUBLE,size-1,0,MPI_COMM_WORLD,&recvreq);
        MPI_Wait(&recvreq,&status);
      }
      if(rank == size-1){
        MPI_Isend(bufferSend,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&sendreq);
        MPI_Irecv(bufferRecv,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&recvreq);
        MPI_Wait(&recvreq,&status);
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
      MPI_Recv(&rank1Avg,1,MPI_DOUBLE,size-1,1,MPI_COMM_WORLD,&status2);
    }
    if (rank == size-1){
      MPI_Ssend(&avgTime,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
    }

    //rank0 write time average of rank0/rank1 out to txt file.
    if(rank == 0){
      cout << "dataSize =  " << n << endl;
      myFile << n << " " << (totalTime/((double)commNum)) << " " << rank1Avg << endl;
    }
  }

  //rank0 close the file it initially opened
  if(rank == 0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
