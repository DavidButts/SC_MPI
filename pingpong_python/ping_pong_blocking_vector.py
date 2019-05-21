from mpi4py import MPI
import numpy as np
import time


comm = MPI.COMM_WORLD

global rank
global size
size = comm.size
rank = comm.rank

#only rank0 needs to open file.
if rank==0:
    f=open("Python_Blocking_Vector.txt","w")

#function to ping pong a message steps number of times
# size of each message is "size".
def trial(steps,n):

    #only rank0 and rank1 allocate memory
    stride=n
    if rank == 0 or rank == (size-1):
        bufferSend = np.ones(n*stride,dtype=np.float64)
        bufferRecv = np.ones(n,dtype=np.float64)

    vec = MPI.DOUBLE.Create_vector(n,1,stride)
    vec.Commit()

    #perform averaging
    tot_time = 0
    for t in range(steps):
        start = MPI.Wtime()
        if rank == 0:
            comm.Send([bufferSend, 1, vec], dest = (size-1), tag = 0)
            comm.Recv([bufferRecv, MPI.DOUBLE], source = (size-1), tag = 1)
        if rank == (size-1):
            comm.Recv([bufferRecv, MPI.DOUBLE], source = 0, tag = 0)
            comm.Send([bufferSend, 1, vec], dest = 0, tag = 1)
        end = MPI.Wtime()
        tot_time += end-start

    vec.Free()

    return np.float64(tot_time/steps)

#loop over message sizes
for n in range(0,25,1):
    t = trial(5000,2**n)

    #rank1 sends its average time to rank0
    if rank == (size-1):
        comm.send(t, dest = 0, tag = 2)

    #rank0 writes results out to text file.
    if rank == 0:
        rank1Time = comm.recv(source = (size-1),tag=2)
        print('n = %d' %(2**n))
        strOut=("%d %f %f\n" %(2**n,t,rank1Time) )
        f.write(strOut)

#only rank0 needs to close file.
if rank==0:
    f.close()

