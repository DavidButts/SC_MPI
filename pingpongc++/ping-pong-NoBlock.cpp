#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

int main(void){

  //This code will write times out to a file of this
  // name.
  string fileName="singleNodeNoBlock.txt";

  // number of times to average over
  int commNum = 10000;

  int totalSize = pow(2,24);
  int rank;
  int size;
  char* bufferSend = NULL;
  char* bufferRecv = NULL;
  double start, end, totalTime;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  MPI_Request recvreq,sendreq;
  MPI_Status status;

  //Allocate arrays only on two ranks to cut down on
  // memory size.
  // also, I allocate two arrays of max datasize (2^24 bytes)
  // and only send first n of those bytes from n =1 to 2^24.
  ofstream myFile;
  if(rank==0){
    myFile.open(fileName);
    bufferSend = (char* ) malloc(totalSize*sizeof(char));
    bufferRecv = (char* ) malloc(totalSize*sizeof(char));
  }
  if(rank==size-1){
    bufferSend = (char* ) malloc(totalSize*sizeof(char));
    bufferRecv = (char* ) malloc(totalSize*sizeof(char));
  }

  //Perform non-blocking and and recv calls
  for( int n = 1; n <= totalSize; n*=2){
    int dataSize = n; // only send first n elements of bufferSend/Recv
    start = MPI_Wtime();
    for (int i =0; i < commNum; i++){
      if(rank ==0){
        MPI_Isend(bufferSend,dataSize,MPI_CHAR,size-1,0,MPI_COMM_WORLD,&sendreq);
        MPI_Irecv(bufferRecv,dataSize,MPI_CHAR,size-1,0,MPI_COMM_WORLD,&recvreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
      if(rank == size-1){
        MPI_Irecv(bufferRecv,dataSize,MPI_CHAR,0,0,MPI_COMM_WORLD,&recvreq);
        MPI_Isend(bufferSend,dataSize,MPI_CHAR,0,0,MPI_COMM_WORLD,&sendreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
    }

    //measure time on rank 0
    end = MPI_Wtime();
    totalTime = end - start;
    if(rank == 0){
      myFile << (double)dataSize << " " << (totalTime/((double)commNum)) << endl;
    }
  }
  // free on only ranks that have memory allocated
  if(rank == 0){
    myFile.close();
    free(bufferSend);
    free(bufferRecv);
  }
  if(rank==size-1){
    free(bufferSend);
    free(bufferRecv);
  }

  MPI_Finalize();
  return 0;
}
