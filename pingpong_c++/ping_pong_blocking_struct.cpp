#include <iostream>
#include <fstream>
#include <stddef.h>
#include <math.h>
#include <mpi.h>

using namespace std;

typedef double real;

struct object{
  double a[20];
  int b;
  char c;
};

int main(void){

  // write results out to a text file with this name
  string fileName="Cpp_Block_Struct.txt";

  // number of runs to perform to average out the noise.
  int commNum = 5000;

  int totalSize = pow(2,24);
  int rank;
  int numRanks;
  double start, end, totalTime;
  object* bufferSend=NULL;
  object* bufferRecv=NULL;

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

  //init datatype
  MPI_Datatype struc;
  int structlen = 3;
  int blocklengths[structlen];
  MPI_Datatype types[structlen];
  MPI_Aint displacements[structlen];
  //double b[20];
  blocklengths[0] = 20;
  types[0] = MPI_DOUBLE;
  displacements[0] = offsetof(object,a);
  //int a = 5;
  blocklengths[1] = 1;
  types[1] = MPI_INT;
  displacements[1] = offsetof(object,b);
  //char c = 'b';
  blocklengths[2] = 1;
  types[2] = MPI_CHAR;
  displacements[2] = offsetof(object,c);
  //create and commit
  MPI_Type_create_struct(structlen,blocklengths,displacements,types,&struc);
  MPI_Type_commit(&struc);

  //only rank 0 writes to txt file
  ofstream myFile;
  if(rank == rank0){
    myFile.open(fileName);
  }

  //loop over size of data
  for( int n = 1; n <= totalSize; n*=2){

    //ranks 1 and 0 allocate memory
    if(rank==rank0){
      bufferSend = (object* ) malloc(n*sizeof(object));
      bufferRecv = (object* ) malloc(n*sizeof(object));
    }
    if(rank==rank1){
      bufferSend = (object* ) malloc(n*sizeof(object));
      bufferRecv = (object* ) malloc(n*sizeof(object));
    }

    //perform average of commNum send and recvs
    totalTime = 0;
    for (int i =0; i < commNum; i++){
      start = MPI_Wtime();
      if(rank == rank0){
        MPI_Ssend(bufferSend,n,struc,rank1,0,MPI_COMM_WORLD);
        MPI_Recv(bufferRecv,n,struc,rank1,0,MPI_COMM_WORLD,&status1);
      }
      else if(rank == rank1){
        MPI_Recv(bufferRecv,n,struc,rank0,0,MPI_COMM_WORLD,&status1);
        MPI_Ssend(bufferSend,n,struc,rank0,0,MPI_COMM_WORLD);
      }
      end = MPI_Wtime();
      totalTime += end - start;
    }

    //ranks 1 and 0 free memory
    if(rank==rank0){
       free(bufferSend);
       free(bufferRecv);
    }
    if(rank==rank1){
       free(bufferSend);
       free(bufferRecv);
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

  MPI_Type_free(&struc);

  //rank0 close file since it was
  //only rank to open it.
  if(rank == rank0){
    myFile.close();
  }

  MPI_Finalize();

  return 0;
}
