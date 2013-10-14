# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Junior
@date:   Mon Oct 14 10:06:18 2013
@email:  bismarckjunior@outlook.com
@brief:  Reservatório Fraturado - Construção do espaçamento horizontal
"""
import numpy as np
from math import factorial as fat
import matplotlib.pyplot as plt
        
#Fratura (direção x)
xf = 100 
kx_fra = 100
nx_fra = 50
px_fra = 1E2

#Formação (direção x)
L = 2*xf
kx_for = 25
nx_for = 50
px_for = 4E2

#Fratura (direção y)
bf_2 = 0.05
ky_fra = kx_fra
ny_fra = 10

#Formação (direção y)
H = 50
ky_for = kx_for
ny_for = 15
py_for = 4E1

###############################################################################

#Fratura na direção x
dx_fra = np.exp(-(np.arange(1.,nx_fra+1)-nx_fra/2)**2/px_fra)
dx_fra = xf*dx_fra/sum(dx_fra)

#Formação direção x
dx_for = np.exp(-(np.arange(1.,nx_for+1)-nx_for)**2/px_for)
dx_for = L*dx_for/sum(dx_for)

#Fratura direção y
dy_fra = bf_2*np.ones(ny_fra)/ny_fra

#Formação direção y
dy_for = np.exp(-(np.arange(1.,ny_for+1)-ny_for)**2/py_for)
dy_for = H*dy_for/sum(dy_for)


#Direção x
dx = np.append(dx_fra, dx_for)
acum_x = [0]
for a in dx:
    acum_x.append(a+acum_x[-1])

#Direção y
dy = np.append(dy_fra, dy_for)[::-1]
acum_y = [0]
for a in dy:
    acum_y.append(a+acum_y[-1])

#Direção x
plt.figure()
plt.title(u'Direção $x$')
plt.plot(dx)
plt.xlabel('Blocos')
plt.ylabel('$\Delta x$', rotation=0, fontsize=16)
plt.grid()

#Direção y
plt.figure()
plt.title(u'Direção $y$')
plt.plot(dy)
plt.xlabel('Blocos')
plt.ylabel('$\Delta y$', rotation=0, fontsize=16)
plt.grid()

#Reservatorio 
plt.figure() 
plt.vlines(acum_x, 0, H+bf_2)
plt.hlines(acum_y, 0, L+xf)
plt.xlabel('$\Delta x$', fontsize=16)
plt.ylabel('$\Delta y$', rotation=0, fontsize=16)
plt.xlim(0, L+xf)
plt.ylim(0, H+bf_2)

plt.show()

for fname, d in zip(['dx', 'dy'], [dx, dy]):
    f = open(fname+'_fraturado.dat', 'w')
    for d_ in d:
        print >> f, "%.5f" % d_
f.close()

for fname, (k_fra, k_for) in zip(['kx', 'ky'], [(kx_fra, kx_for), (ky_fra, ky_for)]):
    f = open(fname+'_fraturado.dat', 'w')
    for row in range(ny_fra):
        for col in range(nx_fra):
            f.write("%7.1f" % k_fra)
        for col in range(nx_for+1):
            f.write("%7.1f" % k_for)
        f.write('\n')
        
    for row in range(ny_for+1):
        for col in range(nx_fra+nx_for+1):
            f.write("%7.1f" % k_for)
        f.write('\n')    
        
        
        
        
        