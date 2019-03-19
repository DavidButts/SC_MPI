#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

int main(void){


  // write results out to a text file with this name
  string fileName="singleNode.txt";

  // number of runs to perform to average out the noise.
  int commNum = 10000;

  int totalSize = pow(2,24);
  int rank;

  // here I allocate a single array of max size
  // and send only parts of it.
  char* buffer = (char* ) calloc(totalSize*sizeof(char));

  double start, end, totalTime;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Status status;

  //only rank 0 writes to txt file
  ofstream myFile;
  if(rank==0){
    myFile.open(fileName);
  }

  //main loop to comm and average
  for( int n = 1; n <= totalSize; n*=2){
    int dataSize = n; //only send the first n elements of buffer
    start = MPI_Wtime();
    for (int i =0; i < commNum; i++){
      if(rank ==0){
        MPI_Ssend(buffer,dataSize,MPI_CHAR,1,0,MPI_COMM_WORLD);
        MPI_Recv(buffer,dataSize,MPI_CHAR,1,0,MPI_COMM_WORLD,&status);
      }
      else{
        MPI_Recv(buffer,dataSize,MPI_CHAR,0,0,MPI_COMM_WORLD,&status);
        MPI_Ssend(buffer,dataSize,MPI_CHAR,0,0,MPI_COMM_WORLD);
      }
    }

    //measure only on rank 0
    end = MPI_Wtime();
    totalTime = end - start;
    if(rank == 0){
      myFile << (double)dataSize << " " << (totalTime/((double)commNum)) << endl;
    }
  }
  if(rank == 0){
    myFile.close();
  }
  free(buffer);
  MPI_Finalize();
  return 0;
}
