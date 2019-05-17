#include <iostream>
#include <stdlib.h>
#include <cstring>
//#include "hdf5.h"
#include <mpi.h>
#include <math.h>
#include <vector>

//timing library
#include "get_walltime.c"
#include "poisson_utils.cpp"

//some physical constants
#define inv_eo 1.2941e11
#define e 1.60217622e-19

// Number Of Boundary nodes in
// a dimension.
// dont change!!
#define NG 2

// Define Global Grid Size (number of cells)
#define NUMCELL_X  4000 //2^9-3
#define NUMCELL_Y  4000 //2^9-3

// number of iterations to perform
// poisson update
#define ITERATIONS 500

// defining real type so I can easily
// switch between the two.
typedef double real;

// only namespace using in program
using namespace std;

// Indexing scheme
inline int index(int i, int j, int numNodeX){
  return i + j*(numNodeX);
}


int main(int argc, char** argv){

  int p;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&p);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //arbitrary grid size
  //assumeing square
  real gridLength = 1.0;

  //convert from cell to node indexing
  // since phi is defined on the nodes.
  // Note: not including boundary nodes.
  int numNodeX = (NUMCELL_X + 1);
  int numNodeY = (NUMCELL_Y + 1);
  int totalNode = numNodeX*numNodeY;

  vector<int> decomp = gridDecomp(size);

  int localNodeX, localNodeY, localTotalNode;

  //distributed remaining grid nodes amongst ranks (load balance)
  localNodeX = numNodeX/decomp[0];
  localNodeX = (rank%decomp[0] < numNodeX%decomp[0]) ? localNodeX : (localNodeX + 1);

  if( decomp[1] == 0){
    localNodeY = numNodeY;
  }
  else{
    localNodeY = numNodeY/decomp[1];
    localNodeY = (int (rank/decomp[0]) < numNodeY%decomp[1]) ? localNodeY : (localNodeY + 1);
  }

  localTotalNode = localNodeX*localNodeY;

  cout << "rank:" << rank << " localNodeY = " << localNodeY << endl;
  cout << "rank:" << rank << " localNodeX = " << localNodeX << endl;

  // grid spacing
  // dx = dy
  real dx = gridLength/NUMCELL_X;
  real dy = gridLength/NUMCELL_Y;

  // allocate grid on the heap and init values to 0
  // this correct since I am only considering
  // conductive boundaries.
  real* phi_0 = (real*) calloc(localTotalNode,sizeof(real));
  real* phi_1 = (real*) calloc(localTotalNode,sizeof(real));
  real* rho   = (real*) calloc(localTotalNode,sizeof(real));

  //allocate MPI read Buffers
  real* phi_N = (real*) calloc(localNodeX,sizeof(real));
  real* phi_S = (real*) calloc(localNodeX,sizeof(real));
  real* phi_E = (real*) calloc(localNodeY,sizeof(real));
  real* phi_W = (real*) calloc(localNodeY,sizeof(real));


  //Intialize a point charge on grid with
  //conductive boundary (V=0).
  //
  // Point Charge distribution located
  // in the middle of each rank's grid.
  if( rank == 0){
    rho[(localNodeX/2) + (localNodeY/2)*localNodeX] = -e;
  }

  //Requests for boundary ranks
  MPI_Request recvEastRank, recvWestRank,recvNorthRank,recvSouthRank;
  MPI_Request sendEastRank, sendWestRank,sendNorthRank,sendSouthRank;

  //Define Neighbors
  int westRank=rank-1;
  int eastRank=rank+1;
  int southRank=rank+decomp[0];
  int northRank=rank-decomp[0];

  // check if neihboring ranks are valid
  int hasNorthRank=rankIsValid(rank,northRank,decomp,0);
  int hasSouthRank=rankIsValid(rank,southRank,decomp,0);
  int hasEastRank=rankIsValid(rank,eastRank,decomp,1);
  int hasWestRank=rankIsValid(rank,westRank,decomp,1);

  //Column Vector data type
  MPI_Datatype col;
  MPI_Type_vector(localNodeY,1,localNodeX,MPI_DOUBLE,&col);
  MPI_Type_commit(&col);

  //timing stamps
  double startTime;
  double endTime;
  get_walltime(&startTime);

  for(int n=0; n<ITERATIONS; ++n){

    //send ghost zones to neighbors
    //
    //North
    if(hasNorthRank){
      MPI_Isend(phi_0,localNodeX,MPI_DOUBLE,northRank,0,MPI_COMM_WORLD,&sendNorthRank);
      MPI_Irecv(phi_N,localNodeX,MPI_DOUBLE,northRank,1,MPI_COMM_WORLD,&recvNorthRank);
    }
    //South
    if(hasSouthRank){
      MPI_Isend(&phi_0[index(0,(localNodeY-1),localNodeX)],localNodeX,MPI_DOUBLE
                ,southRank,1,MPI_COMM_WORLD,&sendSouthRank);
      MPI_Irecv(phi_S,localNodeX,MPI_DOUBLE,southRank,0,MPI_COMM_WORLD,&recvSouthRank);
    }
    //East
    if(hasEastRank){
      MPI_Isend(&phi_0[index((localNodeX-1),0,localNodeX)],1,col,eastRank,2,MPI_COMM_WORLD,&recvEastRank);
      MPI_Irecv(phi_E,localNodeY,MPI_DOUBLE,eastRank,3,MPI_COMM_WORLD,&sendEastRank);
    }
    //West
    if(hasWestRank){
      MPI_Isend(phi_0,1,col,westRank,3,MPI_COMM_WORLD,&recvWestRank);
      MPI_Irecv(phi_W,localNodeY,MPI_DOUBLE,westRank,2,MPI_COMM_WORLD,&sendWestRank);
    }

    //interior jacobi poisson solve
    for(int j=1;j<localNodeY-1;++j){
      for(int i=1;i<localNodeX-1;++i){
        phi_1[index(i,j,localNodeX)] = 0.25*(
                               phi_0[i + 1 + j*localNodeX]
                              +phi_0[i - 1 + j*localNodeX]
                              +phi_0[i + j*localNodeX + localNodeX]
                              +phi_0[i + j*localNodeX - localNodeX]
                              +inv_eo*(dx*dx)*rho[i + j*localNodeX]);
      }
    }

    //interior jacobi poisson solve
    for( int k = 1; k<20; ++k){
      for(int j=1;j<localNodeY-1;++j){
        for(int i=1;i<localNodeX-1;++i){
          phi_1[index(i,j,localNodeX)] = 0.25*(
                               phi_1[i + 1 + j*localNodeX]
                              +phi_1[i - 1 + j*localNodeX]
                              +phi_1[i + j*localNodeX + localNodeX]
                              +phi_1[i + j*localNodeX - localNodeX]
                              +inv_eo*(dx*dx)*rho[i + j*localNodeX]);
        }
      }
    }


    //recieve ghost zones from neighbors
    //North
    if(hasNorthRank){
      MPI_Wait(&recvNorthRank,MPI_STATUS_IGNORE);
      MPI_Wait(&sendNorthRank,MPI_STATUS_IGNORE);
    }
    //South
    if(hasSouthRank){
      MPI_Wait(&sendSouthRank,MPI_STATUS_IGNORE);
      MPI_Wait(&recvSouthRank,MPI_STATUS_IGNORE);
    }
    //East
    if(hasEastRank){
      MPI_Wait(&sendEastRank,MPI_STATUS_IGNORE);
      MPI_Wait(&recvEastRank,MPI_STATUS_IGNORE);
    }
    //West
    if(hasWestRank){
      MPI_Wait(&sendWestRank,MPI_STATUS_IGNORE);
      MPI_Wait(&recvWestRank,MPI_STATUS_IGNORE);
    }

    //update left and right exterior nodes
    for(int j=1;j<localNodeY-1;++j){

      phi_1[j*localNodeX] = 0.25*(
                       phi_E[j]//phi_0[index(1,j,localNodeX)]
                      +phi_0[1+j*localNodeX]
                      +phi_0[(j+1)*localNodeX]
                      +phi_0[(j-1)*localNodeX] //seg_fault
                      +inv_eo*(dx*dx)*rho[j*localNodeX]);

      phi_1[(localNodeX-1)+j*localNodeX] = 0.25*(
                       phi_W[j]//phi_0[index((localNodeX),j,localNodeX)]
                      +phi_0[(localNodeX-2)+j*localNodeX]
                      +phi_0[(localNodeX-1)+(j+1)*localNodeX]
                      +phi_0[(localNodeX-1)+(j-1)*localNodeX]//seg_fault
                      +inv_eo*(dx*dx)*rho[(localNodeX-1)+j*localNodeX]);

    }


    //update top and bottom exterior ndoes
    for(int i=1;i<localNodeX-1;++i){

      phi_1[index(i,0,localNodeX)] = 0.25*(
                             phi_0[i+1]
                            +phi_0[i-1]//seg-faults
                            +phi_0[i+localNodeX]
                            +phi_N[i]//+phi_0[index(i,-1,localNodeX)]
                            +inv_eo*(dx*dx)*rho[i]);

      phi_1[index(i,(localNodeY-1),localNodeX)] = 0.25*(
                             phi_0[i+1 + (localNodeY-1)*localNodeX]
                            +phi_0[i-1 + (localNodeY-1)*localNodeX]//seg-faults
                            +phi_S[i] //+phi_0[index(i,(localNodeY),localNodeX)]
                            +phi_0[i + (localNodeY-2)*localNodeX]
                            +inv_eo*(dx*dx)*rho[i + (localNodeY-1)*localNodeX]);


    }

    //corners
    //NorthWest
    phi_1[index(0,0,localNodeX)] = 0.25*(
                             phi_0[index(1,0,localNodeX)]
                            +phi_W[0]
                            +phi_0[index(0,1,localNodeX)]
                            +phi_N[0]
                            +inv_eo*(dx*dx)*rho[index(0,0,localNodeX)]);
    //NorhtEast
    phi_1[index((localNodeX-1),0,localNodeX)] = 0.25*(
                             phi_E[(localNodeY-1)]
                            +phi_0[index((localNodeX-2),0,localNodeX)]
                            +phi_0[index((localNodeX-1),1,localNodeX)]
                            +phi_N[(localNodeX-1)]
                            +inv_eo*(dx*dx)*rho[index((localNodeX-1),0,localNodeX)]);
    //SouthWest
    phi_1[index(0,(localNodeY-1),localNodeX)] = 0.25*(
                             phi_0[index(1,(localNodeY-1),localNodeX)]
                            +phi_W[(localNodeY-1)]
                            +phi_0[index(0,(localNodeY-2),localNodeX)]
                            +phi_S[0]
                            +inv_eo*(dx*dx)*rho[index(0,(localNodeY-1),localNodeX)]);

    //SouthEast
    phi_1[index((localNodeX-1),(localNodeY-1),localNodeX)] = 0.25*(
                             phi_0[index((localNodeX-2),(localNodeY-1),localNodeX)]
                            +phi_E[(localNodeY-1)]
                            +phi_0[index((localNodeX-2),(localNodeY-1),localNodeX)]
                            +phi_S[(localNodeX-1)]
                            +inv_eo*(dx*dx)*rho[index((localNodeX-1),(localNodeY-1),localNodeX)]);

    //swap pointers
    real* tmp = phi_0;
    phi_0 = phi_1;
    phi_1 = tmp;

  }

  //calculate runtime/(number of iterations)
  get_walltime(&endTime);
  double runTimePerIter = (endTime - startTime)/((double)ITERATIONS);
  cout << "time per iteration = " <<  runTimePerIter << " on rank:" << rank << endl;


  //print final phi field to HDF5 file
  if(rank==0){
    printH5(phi_0,&runTimePerIter,&localNodeX,&localNodeY);
  }
  // a few quick checks
  if(rank==0){
    cout << "particle at x = " << localNodeX/2 << " and y = " << localNodeY/2 << endl;
    cout << "on Paricle = " << phi_0[index(localNodeX/2,localNodeY/2,localNodeX)] << endl;
    cout << "further  away = " << phi_0[index(localNodeX/4,localNodeY/4,localNodeX)] << endl;
    cout << "furthest away = " << phi_0[index(localNodeX/8,localNodeY/8,localNodeX)] << endl;
    cout << "furthest away = " << phi_0[index(2,2,localNodeX)] << endl;
  }
  if(rank==1){
   cout << "particle at x = " << localNodeX/2 << " and y = " << localNodeY/2 << endl;
  }

  //free memory
  free(phi_0);
  //free(phi_1);
  free(rho);

  //frew MPI type
  MPI_Type_free(&col);

  //MPI finalize
  MPI_Finalize();

  return 0;
}
