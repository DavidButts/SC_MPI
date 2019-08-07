import time
import numpy as np
from mpi4py import MPI
import numba as nb

# Create a communicator
comm = MPI.COMM_WORLD

# Get rank and size of comm
rank = comm.rank
size = comm.size

#Set the dimension of square grid
dim = size*100

# apply z**2 + c to each pixel in z
@nb.jit(cache=True, nogil=True, nopython=True, parallel=True)
def julia(c,z):
	it = 0
	max_iter = 100
	while(it < max_iter):
		for y in range(z.shape[0]):
			for x in range(z.shape[1]):
				if abs(z[y][x]) < 10:
					z[y][x] = z[y][x]**2 + c
		it += 1
	return z

# Make the choce that this is always an int (no remainder)
rows_per_core = int(dim/size)


x_min = -1.8
x_max = 1.8
y_min = -1.8j
y_max = 1.8


iters = 10
tot_time = 0.0
for t in range(iters):
	# Create your local peice of the board
	my_z = np.zeros((rows_per_core,dim), dtype='complex128')
	for row in range(rows_per_core):
       		my_z[row] = np.linspace(x_min,x_max,dim) - np.linspace(y_min,y_max,dim)[(rank*rows_per_core)+row]
	        # Each row has all x values and then the corresponsing y values are selected out

	# Starting a timer
	start = time.perf_counter()
	# Each rank start solving
	julia((-.4+.6j),my_z)
	# Stop the timer
	end = time.perf_counter()
	tot_time += end-start


# Calculate the total time taken on each rank
tot_time = np.array(tot_time/iters)

# Buffer for max time
global_time = np.zeros(1)

# Reduce maximum time taken to rank 0
comm.Reduce(tot_time, global_time, op=MPI.MAX, root=0)


if rank == 0:
	# Print out the max time taken
	print(global_time[0])
	
	# Create the solution board on rank 0
	#global_board = np.zeros((dim,dim),dtype='complex128')
	
	# Rank 0 puts in its data
	#global_board[0:rows_per_core] = my_z

# Flag to bring data to rank 0
construct_sol = False

if construct_sol == True:
	# Loop through ranks to send data to rank 0
	for r in range(1,size):
		if rank == 0:
			# Create a buffer for recieving the other ranks data
			buff = np.zeros((rows_per_core,dim),dtype='complex128')

			# Rank 0 recive rank r's data
			comm.Recv([buff, MPI.FLOAT], source=r, tag = r)
		
			# Place the data into the global board
			for row in range(buff.shape[0]):
				global_board[(r*rows_per_core)+row] = buff[row]
		else:
			# Rank r send data to rank 0
			comm.Send([my_z, MPI.FLOAT], dest=0, tag=r)

# Rank 0 print result
#if rank == 0:
#	print(global_board)
