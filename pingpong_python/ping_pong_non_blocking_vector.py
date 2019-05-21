from mpi4py import MPI
import numpy as np
import time


comm = MPI.COMM_WORLD

global size
global rank
size = comm.size
rank = comm.rank

#only rank0 opens files
if rank==0:
  f=open("Python_NoBlocking_Vector.txt","w")

#function
def trial(steps,n):

    #only rank 1 or rank 2 allocates memory
    stride=n
    if rank==0 or rank==(size-1):
        bufferSend = np.ones(n*stride,dtype=np.float64)
        bufferRecv = np.ones(n,dtype=np.float64)

    vec = MPI.DOUBLE.Create_vector(n,1,stride)
    vec.Commit()

    #loop that performs time averaging of ping pong
    tot_time=0
    for t in range(steps):
        start = MPI.Wtime()
        if rank == 0:
            reqS = comm.Isend([bufferSend, 1, vec], dest = (size-1), tag = 0)
            reqR = comm.Irecv([bufferRecv, MPI.DOUBLE], source = (size-1), tag = 1)
            reqR.Wait()
            reqS.Wait()

        if rank == size-1:
            reqR = comm.Irecv([bufferRecv, MPI.DOUBLE], source = 0, tag = 0)
            reqS = comm.Isend([bufferSend, 1, vec], dest = 0, tag = 1)
            reqR.Wait()
            reqS.Wait()

        end = MPI.Wtime()
        tot_time += end - start

    vec.Free()

    return tot_time/steps

#loop over data size
for n in range(0,25,1):
    t = trial(5000,2**n)

    #rank1 sends rank0 its time average
    if rank == (size-1):
        comm.send(t, dest = 0, tag = 2)

    #rank0 writes all results out to txt file.
    if rank == 0:
        rank1Time = comm.recv(source = (size-1),tag=2)
        print('n = %d' %(2**n))
        strOut=("%d %f %f\n" %(2**n,t,rank1Time) )
        f.write(strOut)

#only rank 0 closes file
if rank==0:
  f.close()
