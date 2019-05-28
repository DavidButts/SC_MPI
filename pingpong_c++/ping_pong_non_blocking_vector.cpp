#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

int main(void){

  //This code will write times out to a file of this
  // name.
  string fileName="Cpp_NoBlock_Vector.txt";

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
    int stride = 50;
    if(rank==0){
      bufferSend = (real* ) malloc(n*stride*sizeof(real));
      bufferRecv = (real* ) malloc(n*sizeof(real));
    }
    if(rank==size-1){
      bufferSend = (real* ) malloc(n*stride*sizeof(real));
      bufferRecv = (real* ) malloc(n*sizeof(real));
    }

    MPI_Datatype vec;
    MPI_Type_vector(n,1,stride,MPI_DOUBLE,&vec);
    MPI_Type_commit(&vec);


    //loop that will perform non-block send/recv
    //only actual send/recv is timed.
    totalTime=0;
    for (int i =0; i < commNum; i++){
      start = MPI_Wtime();
      if(rank ==0){
        MPI_Isend(bufferSend,1,vec,size-1,1,MPI_COMM_WORLD,&sendreq);
        MPI_Irecv(bufferRecv,n,MPI_DOUBLE,size-1,0,MPI_COMM_WORLD,&recvreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
      if(rank == size-1){
        MPI_Irecv(bufferRecv,n,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&recvreq);
        MPI_Isend(bufferSend,1,vec,0,0,MPI_COMM_WORLD,&sendreq);
        MPI_Wait(&recvreq,&status);
        MPI_Wait(&sendreq,&status);
      }
      end = MPI_Wtime();
      totalTime += end - start;
    }

    MPI_Type_free(&vec);

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

  //rank0 close the file it initially opened
  if(rank == 0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
