export HDF5_CXX=mpic++
export HDF5_CLINKER=mpic++

#NOTE:
# to build and run you need the following commands:
#
#
# module purge
# module load GCC/7.3.0-2.30
# module load OpenMPI/3.1.1-CUDA
# module load HDF5/1.10.2
# module load CUDA/9.2.88



CXX=h5c++

NCC=nvcc

poissonSolvePar.o:poissonSolvePar.cpp
	$(CXX) -O3 -ffast-math  poissonSolvePar.cpp
all:poissonSolvePar.o
	$(CXX) -O3 -ffast-math -o poissonSolvePar poissonSolvePar.cpp
clean:
	rm *.h5
	rm *.o
	rm poissonSolvePar
	rm poissonSolveSer
