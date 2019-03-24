from mpi4py import MPI
import numpy as np
import time
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD

size = comm.size
rank = comm.rank


def trial(steps,size):
    mess_in = np.empty((size),dtype=np.float64)
    mess_out = np.ones((size),dtype=np.float64)
    for t in range(steps):
        start = MPI.Wtime()
        if rank == 0:
            req = comm.Isend([mess_out, MPI.DOUBLE], dest = 1, tag = 0)
            comm.Recv([mess_in, MPI.DOUBLE], source = 1, tag = 1)
        if rank == 1:
            comm.Recv([mess_in, MPI.DOUBLE], source = 0, tag = 0)
            req = comm.Isend([mess_out, MPI.DOUBLE], dest = 0, tag = 1)
        req.Wait()
        end = MPI.Wtime()
        tot_time = end - start
    return tot_time/steps

data = []
for n in range(0,25,2):
    t = trial(10,2**n)
    if rank == 0:
        data.append(t)
        print('n =',2**n)
        # print('Time per send-recv:',t)
if rank == 0:
    plt.loglog(range(0,25,2),data)
    plt.show()
