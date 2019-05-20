from mpi4py import MPI
import numpy as np
import time
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD

global size
global rank
size = comm.size
rank = comm.rank

#only rank0 opens files
if rank==0:
  f=open("Python_NoBlocking_Struct.txt","w")

# This code defines a C struct type in python
#
# Struct Object{
#   double a[20];
#   int b;
#   char c;
# };
#
# This is exactly 168 Bytes which can be Byte aligned like
# C-struct is!
Object = np.dtype([('a','(20,)f8'),('b',np.int32),('c','U1')])

#init MPI STRUCT TYPE
offsets = [Object.fields['a'][1], Object.fields['b'][1],Object.fields['c'][1]]
types = [MPI.DOUBLE,MPI.INT,MPI.BYTE]
struc = MPI.Datatype.Create_struct([1]*3, offsets, types)
struc.Commit()

#function
def trial(steps,n):

    #only rank0 and rank1 allocate memory
    if rank == 0 or rank == (size-1):
        bufferSend = np.zeros(n,dtype=Object)
        bufferRecv = np.zeros(n,dtype=Object)

    #loop that performs time averaging of ping pong
    tot_time=0
    for t in range(steps):
        start = MPI.Wtime()
        if rank == 0:
            reqS = comm.Isend([bufferSend, struc], dest = (size-1), tag = 0)
            reqR = comm.Irecv([bufferRecv, struc], source = (size-1), tag = 1)
            reqR.Wait()
            reqS.Wait()

        if rank == size-1:
            reqR = comm.Irecv([bufferRecv, struc], source = 0, tag = 0)
            reqS = comm.Isend([bufferSend, struc], dest = 0, tag = 1)
            reqR.Wait()
            reqS.Wait()

        end = MPI.Wtime()
        tot_time += end - start

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

struc.Free()

