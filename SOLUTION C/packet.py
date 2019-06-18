from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np   
from math import sqrt

data=np.loadtxt('pc.dat')

fig=plt.figure()

ax= fig.add_subplot(111) # 2d Contorno
x=data[:,0] # numero de linhas
y=data[:,1]
z=data[:,2] # fator multiplicativo para meV

xi=np.linspace(-800,800,500) # selecionando alguns pontos para o grid
yi=np.linspace(-800,800,500)

X,Y = np.meshgrid(xi,yi)
Z = griddata(x, y, z, xi, yi, interp='linear') #criando a matriz do eixo Z

#---- contorno-----------------------------------
surf=ax.contourf(X, Y, Z,600,alpha=0.80,cmap=cm.jet)
#linha=ax.contour(X, Y, Z, 100, colors='black', linewidths=0.6)
fig.colorbar(surf)
#------------------------------------------------


#----------ajustando eixos-----------------------
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)

ax.set_xlabel(u"$X$")
ax.set_ylabel(u"$Y$")
plt.show()

