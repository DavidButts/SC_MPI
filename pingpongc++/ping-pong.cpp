#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

int main(void){

  // write results out to a text file with this name
  string fileName="Cpp_Block_Double.txt";

  // number of runs to perform to average out the noise.
  int commNum = 5000;

  int totalSize = pow(2,24);
  int rank;
  int numRanks;
  double start, end, totalTime;
  real* buffer=NULL;

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

  //loop over size of data
  for( int n = 1; n <= totalSize; n*=2){

     //ranks 1 and 0 allocate memory
     if(rank==rank0){
       buffer = (real* ) malloc(n*sizeof(real));
     }
     if(rank==rank1){
       buffer = (real* ) malloc(n*sizeof(real));
     }

    //perform average of commNum send and recvs
    totalTime = 0;
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
      end = MPI_Wtime();
      totalTime += end - start;
    }

    //ranks 1 and 0 free memory
     if(rank==rank0){
       free(buffer);
     }
     if(rank==rank1){
       free(buffer);
     }

    //all ranks compute thier average
    double avgTime = (totalTime/((double)commNum));

    //rank1 sends rank0 its average time
    if(rank == rank0){
      MPI_Recv(&rank1Avg,1,MPI_DOUBLE,rank1,1,MPI_COMM_WORLD,&status2);
    }
    else if (rank == rank1){
      MPI_Ssend(&avgTime,1,MPI_DOUBLE,rank0,1,MPI_COMM_WORLD);
    }

    //rank0 writes its time average and rank1's time avg out to file
    if(rank == rank0){
      cout << "dataSize =  " << n << endl;
      myFile << n << " " << (totalTime/((double)commNum)) << " " << rank1Avg << endl;
    }
  }

  //rank0 close file since it was
  //only rank to open it.
  if(rank == rank0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
