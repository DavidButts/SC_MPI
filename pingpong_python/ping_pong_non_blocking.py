from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD

size = comm.size
rank = comm.rank
j = 1 # other rank
#mess = np.array([0],dtype='i')
mess = np.ones((1000),dtype='i')
if rank == 1:
	mess[0] = 74

steps = 100

tot_time = 0
for t in range(steps):
	start = time.time()
	if rank == 0:
		#print('Rank:',rank,'has data',mess)
		comm.Isend([mess, MPI.INT], dest = 1, tag = 0)
		comm.Irecv([mess, MPI.INT], source = 1, tag = 1)
	if rank == j:
		#print('Rank:',rank,'has data',mess)
		comm.Irecv([mess, MPI.INT], source = 0, tag = 0)
		comm.Isend([mess, MPI.INT], dest = 0, tag = 1)
	end = time.time()
	tot_time += end - start
print('Time per send-recv:',tot_time/steps,'s')
