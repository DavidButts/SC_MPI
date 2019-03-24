from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD

rank = comm.rank
size = comm.size

# h = .1
# m = int(1/h - 1)
m = 3
u = np.zeros(shape=(m+2,m+2))
my_rows = int(m/size)
remaining_rows = m%size
if rank < remaining_rows:
    my_rows += 1
print('Rank:',rank,'has',my_rows)

# # Solving u''(x,y) = -2 pi^2 sin(pi x)sin(pi y)
# def jacobi(h,iterations,u,f):
#     temp_u = np.zeros_like(u) # Create a temporary array to store new locations
#     for k in range(iterations):
#         for i in range(1,m+1):
#             for j in range(1,m+1):
#                 temp_u[i,j] = (1/4) * ( u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] ) - (h**2/4)*(-2*np.pi**2*np.sin(np.pi*(i*h))*np.sin(np.pi*(j*h)))
#         u = np.copy(temp_u)
#     return u
#
# if rank == 0:
#     plt.imshow(jacobi(h,100,u,np.ones_like(u)))
#     plt.colorbar()
#     plt.show()
