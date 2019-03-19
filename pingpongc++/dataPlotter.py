import sys
import decimal as dec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

#inFile = open("twoNodeNoBlock.txt",'r')
inFile2 = open('singleNodeNoBlock.txt','r')

x, y = np.genfromtxt(inFile2, unpack=True)
#x2,y2= np.genfromtxt(inFile2, unpack=True)

ub=15
lb=len(x)

par =  np.polyfit(x[ub:lb], y[ub:lb], 1)
latency = np.average(y[1:5])
print par
#bestFit=[par[1] + par[0]*i for i in x[ub:lb]]
L=[latency for i in x]
#bestFit=[  + m*i for i in x[ub:lb]]
m=(y[-1] - y[-2])/(x[-1] - x[-2])
print m
print str("{:.2E}".format(dec.Decimal(1.0/m)))
bestFit=[ 4.3e-6+m*i for i in x[ub:lb]]


plt.loglog(x, y,'.b',label="Comm single node", ms = 15)
plt.loglog(x[ub:lb],bestFit,'-r',label="BandWidth = "+str("{:.2E}".format(dec.Decimal(1.0/m)))+"Bytes/Sec.",ms = 15,linewidth=2.5)
plt.loglog(x,L,'-g',label="Latency = "+str(latency)+"Sec.",ms = 15,linewidth=2.5)
plt.legend(loc=2)
plt.xlabel('Bytes')
plt.ylabel('Tpp')
plt.suptitle('Tpp Vs. Bytes Non-Blocking Comm. (Intra Node)')
plt.show()

