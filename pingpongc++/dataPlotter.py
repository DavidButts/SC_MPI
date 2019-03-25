import sys
import decimal as dec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

inFile = open("doubleDataNoFlag.txt",'r')
inFile2 = open("doubleDataO3.txt",'r')
inFile3 = open("python_blocking.txt",'r')
inFileNB = open("cpp_noBlock_Double.txt",'r')

x, y0, y1 = np.genfromtxt(inFile, unpack=True)
x_03, y0_03, y1_03 = np.genfromtxt(inFile2, unpack=True)
x_py, y0_py, y1_py = np.genfromtxt(inFile3, unpack=True)
x_NB, y0_NB, y1_NB = np.genfromtxt(inFileNB, unpack=True)

yAvg = (y1+y0)/2.0
yAvg_03 = (y1_03+y0_03)/2.0
yAvg_py = (y1_py+y0_py)/2.0
yAvg_NB = (y1_NB+y0_NB)/2.0

plt.loglog(x*8,yAvg,label="C++ No Flag")
plt.loglog(x_03*8,yAvg_03,label="C++ -03")
plt.loglog(x_py*8,yAvg_py,label="Python")
plt.loglog(x_NB*8,yAvg_NB,label="C++ Non-Blocking")
plt.legend(loc=2)
plt.xlabel('Bytes')
plt.ylabel('Time (Seconds)')
plt.suptitle('Time V.s. Bytes Ping-Pong')
plt.show()

