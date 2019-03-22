from mpi4py import MPI
import numpy as np
import time
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD

size = comm.size
rank = comm.rank

def trial(steps,size):
    mess = np.ones((size),dtype='d')
    print(mess.shape)
    tot_time = 0

    for t in range(steps):
        start = MPI.Wtime()
		
        if rank == 0:
        	comm.Isend(mess, dest = 1, tag = 0)
        	req = comm.Irecv(mess, source = 1, tag = 1)
        if rank == 1:
            req = comm.Irecv(mess, source = 0, tag = 0)
            comm.isend(mess, dest = 0, tag = 1)

        end = MPI.Wtime()
        tot_time += end-start

        return tot_time/steps

data = []
for n in range(0,25,2):
    t = trial(50000,2**n)
    if rank == 0:
        data.append(t)
        # print('n =',2**n)
        # print('Time per send-recv:',t)
if rank == 0:
    plt.loglog(range(0,25,2),data)
    plt.show()
