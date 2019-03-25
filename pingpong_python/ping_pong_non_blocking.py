from mpi4py import MPI
import numpy as np
import time
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD

size = comm.size
rank = comm.rank

f=open("python_nonblocking.txt","w")

def trial(steps,size):
    mess_in = np.empty((size),dtype=np.float64)
    mess_out = np.ones((size),dtype=np.float64)
    tot_time=0
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
        tot_time += end - start

    return tot_time/steps

for n in range(0,25,2):
    t = trial(5000,2**n)

    if rank == 1:
        comm.send(t, dest = 0, tag = 2)

    if rank == 0:
        rank1Time = comm.recv(source = 1,tag=2)
        print('n = %d' %(2**n))
        strOut=("%d %f %f\n" %(2**n,t,rank1Time) )
        f.write(strOut)

f.close()
