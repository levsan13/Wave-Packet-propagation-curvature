# coding: utf-8
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
import sys

def plotgraf(file,labelx,labely,labelbar):
    data1=np.loadtxt(file)
    x=data1[:,0]
    y=data1[:,1]    
    row = int(np.sqrt(len(x)))
    z =np.zeros(row**2)
    x_=np.zeros(row**2)
    y_=np.zeros(row**2)
    for r in range(0,row**2):
         x_[r]=(data1[:,0][r])
         y_[r]=(data1[:,1][r])
         z[r]=(data1[:,2][r])
    x_=x_.reshape((row,row))
    y_=y_.reshape((row,row))
    z=z.reshape((row,row))
    cs1 = plt.contourf(x_, y_, z, 25)
    #plt.title(title)
    plt.xlabel(labelx,fontsize=15)
    plt.ylabel(labely,fontsize=15)
    #plt.axis([-100,100,-100,100]) #Caso espcifico
    plt.colorbar().set_label(labelbar)
    plt.show()
    return

metod=sys.argv[1]

ar1=sys.argv[2]
ar2=sys.argv[3]
ar3=sys.argv[4]
ar4=sys.argv[5]
plotgraf(ar1,ar2,ar3,ar4)
