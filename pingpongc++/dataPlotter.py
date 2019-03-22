import sys
import decimal as dec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

inFile = open("doubleDataNoFlag.txt",'r')
inFile2 = open("doubleDataO3.txt",'r')
x, y0, y1 = np.genfromtxt(inFile, unpack=True)
x_03, y0_03, y1_03 = np.genfromtxt(inFile2, unpack=True)

yAvg = (y1+y0)/2.0
yAvg_03 = (y1_03+y0_03)/2.0

plt.loglog(x*8,yAvg,label="No Flag")
plt.loglog(x_03*8,yAvg_03,label="03")
plt.legend(loc=2)
plt.xlabel('Bytes')
plt.ylabel('Time (Seconds)')
plt.suptitle('Time V.s. Bytes Ping-Pong')
plt.show()

