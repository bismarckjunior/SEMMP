# -*- coding: iso-8859-1 -*-
"""
Created on Thu Jul 05 19:10:11 2012

@author: Bismarck Gomes Souza Júnior
"""
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from pylab import *
from numpy import *
import os

a = 'SPE_A1_1_000010.vtk'
b = 'solucao_real.vtk'
c = 'geometrybc.dat'
simu = [[float(a) for a in li.split() ]for li in open(a, 'r').readlines()]
real = [[float(a) for a in li.split() ]for li in open(b, 'r').readlines()]
#geo  = [[float(a) for a in li.split() ]for li in open(c, 'r').readlines()]
nx = len(simu)
ny = len(simu[-1])
X, Y = meshgrid(range(nx),range(ny))

################################################################################
erro = abs(array(simu)-array(real))/array(real)

#ax = plt.figure().add_subplot(111, projection='3d')
#ax.plot_surface(X,Y,erro)
#ax.plot_wireframe(X,Y,erro)
#plt.imshow(-array(geo), cmap=get_cmap('hot'), interpolation='nearest')
plt.imshow(erro*100)#, interpolation = 'nearest')
plt.colorbar()
plt.xticks(range(ny))
plt.yticks(range(nx))
plt.grid()
plt.title(u'Diferença relativa percentual (%)')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(9,3,'ko', markersize=8)
plt.plot(4,4, 'k^', markersize=8)
plt.legend(['Produtor', 'Injetor'], 'lower left', numpoints=1)
plt.savefig('SPE_A1-1_diferenca_relativa.png')
plt.show()

