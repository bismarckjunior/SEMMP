# -*- coding: iso-8859-1 -*-
"""
Created on Mon Dec 17 17:00:24 2012

@author: Bismarck Gomes Souza Júnior
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import os

vtks = filter(lambda x: x.endswith('.vtk'), os.listdir('.'))
vtks.sort()

ims, ims2,ims3 = [], [],[]
min_, max_ = 1E10, -1E10
pressures = []
mxs = []
for vtk in vtks:
  pressures.append([[float(p) for p in line.split()[:]] for line in open(vtk).readlines()])
  mx = max([max(p) for p in pressures[-1]])
  mn = min([min(p) for p in pressures[-1]])
  mxs.append(mx)
  print '%s: %g / %g' % (vtk, mn, mx)
  if mx > max_: max_ = mx
  if mn < min_ and mn != 0: min_ = mn
  
ny = len(pressures[-1])
nx = len(pressures[-1][-1])


#'''
# Reservoir Performance #######################################################
fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.grid(True)
text = ax.text(0,len(pressures[-1])//50+1, '')
plt.xlim(0, nx-1)
plt.ylim(0, ny-1)

def animate(i):
  cax = ax.imshow(pressures[i], norm=plt.Normalize(min_,max_)) #interpolation='nearest'
  #plt.contour(pressures[i])
  fig.savefig('fig_%04d.png' % (i))
  text.set_text(vtks[i])
  return ax, text

ani = animation.FuncAnimation(fig, animate, range(len(pressures)), interval=1000, 
                              blit=True, repeat_delay=2000)

cax = ax.imshow([[min_ for i in range(len(p))] for p in pressures[0]], 
                  norm=plt.Normalize(min_,max_))
text.set_text('')
plt.title(u'Distribuição de Pressão')
plt.colorbar(cax)
###############################################################################

plt.show()