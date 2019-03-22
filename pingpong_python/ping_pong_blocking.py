from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD

size = comm.size
rank = comm.rank

mess = np.ones((1000),dtype='i')
if rank == 1:
        mess[0] = 74

steps = 5000
tot_time = 0

for t in range(steps):
	start = MPI.Wtime()
	if rank == 0:
        	comm.Send([mess[0:10], MPI.INT], dest = 1, tag = 0)
        	comm.Recv([mess[0:10], MPI.INT], source = 1, tag = 1)
	if rank == 1:
        	comm.Recv([mess[0:10], MPI.INT], source = 0, tag = 0)
        	comm.Send([mess[0:10], MPI.INT], dest = 0, tag = 1)
	end = MPI.Wtime()
	tot_time += end-start

print('Time per send-recv:',tot_time/steps,'s')
