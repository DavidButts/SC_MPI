import sys
import decimal as dec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

inFile = open("Cpp_Block_Vector.txt",'r')
inFile2 = open("Cpp_NoBlock_Vector.txt",'r')
inFile3 = open("Python_Blocking_Vector.txt",'r')
inFileNB = open("Python_NoBlocking_Vector.txt",'r')

x, y0, y1 = np.genfromtxt(inFile, unpack=True)
x_03, y0_03, y1_03 = np.genfromtxt(inFile2, unpack=True)
x_py, y0_py, y1_py = np.genfromtxt(inFile3, unpack=True)
x_NB, y0_NB, y1_NB = np.genfromtxt(inFileNB, unpack=True)

yAvg = (y1+y0)/2.0
yAvg_03 = (y1_03+y0_03)/2.0
yAvg_py = (y1_py+y0_py)/2.0
yAvg_NB = (y1_NB+y0_NB)/2.0

font = {'size'   : 22}
lw = {'linewidth':4}

plt.figure(figsize=(7,6))
plt.plot(x*8,yAvg,label="C++ Blocking",**lw)
plt.plot(x_03*8,yAvg_03,label="C++ Non-Blocking",**lw)
plt.plot(x_py*8,yAvg_py,'--',label="Python Blocking",**lw)
plt.plot(x_NB*8,yAvg_NB,'--',label="Python Non-Blocking",**lw)
plt.legend(loc=2,prop={'size':16})
plt.xlabel('Bytes',**font)
plt.xscale('log')
plt.yscale('log')
plt.grid(alpha=0.8)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.ylabel('Time (Seconds)',**font)
plt.suptitle('Time vs. Bytes Ping-Pong Type: MPI Vector',**font)
plt.savefig('Time_Bytes_log.png',bbox_inches='tight')
plt.show()


plt.plot(x*8,yAvg,label="C++ Blocking")
plt.plot(x_03*8,yAvg_03,label="C++ Non-Blocking")
plt.plot(x_py*8,yAvg_py,label="Python Blocking")
plt.plot(x_NB*8,yAvg_NB,label="Python Non-Blocking")
plt.legend(loc=2)
plt.xlabel('Bytes')
#plt.xscale('log')
#plt.yscale('log')
plt.ylabel('Time (Seconds)')
plt.suptitle('Time V.s. Bytes Ping-Pong Type: MPI Vector')
plt.savefig('Time_Bytes_linear.png')
plt.show()

