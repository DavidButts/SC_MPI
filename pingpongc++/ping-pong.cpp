#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

int main(void){

  // write results out to a text file with this name
  string fileName="doubleDataO3.txt";

  // number of runs to perform to average out the noise.
  int commNum = 5000;

  int totalSize = pow(2,24);
  int rank;
  int numRanks;
  double start, end, totalTime;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numRanks);

  int rank0 = 0;
  int rank1 = (numRanks-1);
  cout << rank1 << endl;
  cout << "myrank = " << rank << endl;

  double rank1Avg;

  MPI_Status status1;
  MPI_Status status2;

  //only rank 0 writes to txt file
  ofstream myFile;
  if(rank == rank0){
    myFile.open(fileName);
  }

  //main loop to comm and average
  for( int n = 1; n <= totalSize; n*=2){
    totalTime = 0;
    real* buffer = (real* ) malloc(n*sizeof(real));
    for (int i =0; i < commNum; i++){
      start = MPI_Wtime();
      if(rank == rank0){
        MPI_Ssend(buffer,n,MPI_DOUBLE,rank1,0,MPI_COMM_WORLD);
        MPI_Recv(buffer,n,MPI_DOUBLE,rank1,0,MPI_COMM_WORLD,&status1);
      }
      else if(rank == rank1){
        MPI_Recv(buffer,n,MPI_DOUBLE,rank0,0,MPI_COMM_WORLD,&status1);
        MPI_Ssend(buffer,n,MPI_DOUBLE,rank0,0,MPI_COMM_WORLD);
      }
      //measure only on rank 0
      end = MPI_Wtime();
      totalTime += end - start;
    }
    free(buffer);
    double avgTime = (totalTime/((double)commNum));
    if(rank == rank0){
      MPI_Recv(&rank1Avg,1,MPI_DOUBLE,rank1,1,MPI_COMM_WORLD,&status2);
    }
    else if (rank == rank1){
      MPI_Ssend(&avgTime,1,MPI_DOUBLE,rank0,1,MPI_COMM_WORLD);
    }
    if(rank == rank0){
      cout << "dataSize =  " << n << endl;
      myFile << n << " " << (totalTime/((double)commNum)) << " " << rank1Avg << endl;
    }
  }

  if(rank == rank0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
