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

'''fig=plt.figure()
#plt.title('Reservoir Performance')
fig2 = plt.figure(2)
ax = fig.add_subplot(111)
ax2 = fig2.add_subplot(121)
ax3 = fig2.add_subplot(122)

ax.set_title('Reservoir Performance')
ax2.set_title('Reservoir Distribuition')
ax3.set_title('Max Pressure Performance')
'''
ims, ims2,ims3 = [], [],[]
min_, max_ = 1E10, -1E10
pressures = []
mxs = []
for vtk in vtks:
  pressures.append([[float(p) for p in line.split()[:]] for line in open(vtk).readlines()])
  mx = max([max(p) for p in pressures[-1]])
  mn = min([min(p) for p in pressures[-1]])
  mxs.append(mx)
  #print '%s: %g / %g' % (vtk, mn, mx)
  if mx > max_: max_ = mx
  if mn < min_ and mn != 0: min_ = mn
  
min_=7300.
max_=7500.

ny = len(pressures[-1])
nx = len(pressures[-1][-1])

'''  
text = ax.text(10, 10, '')#, transform=ax.transAxes)
for vtk, pressure in zip(vtks, pressures):
  cax = ax.imshow(pressure, norm=plt.Normalize(min_,max_))  
  cax2= ax2.imshow(pressure)
  cax3= ax3.imshow(pressure)
  text.set_text(vtk)
  fig.savefig(vtk[:-3]+'png')
  ims.append([cax])
  ims2.append([cax2])
  ims3.append([cax3])
  
#ax3.plot(mxs,'-o')  
#ax3.set_xlim(-0.5, len(mxs)-0.5)
#ax3.set_xticks(range(len(mxs)))
ax3.grid(True)

ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=2000)
ani2 = animation.ArtistAnimation(fig2, ims2, interval=1000, repeat_delay=2000)
ani3 = animation.ArtistAnimation(fig2, ims3, interval=1000, repeat_delay=2000)
fig.colorbar(cax)
plt.show()'''

#'''
# Reservoir Performance #######################################################
fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.grid(True)
text = ax.text(0,len(pressures[-1])//50+1, '')
plt.xlim(0, nx-1)
#plt.ylim(0, ny-1)

def animate(i):
  p = [[p if p!=0.0 else (min_+max_)/2. for p in pr] for pr in pressures[i]]
  cax = ax.imshow(p, norm=plt.Normalize(min_,max_)) #interpolation='nearest'
  #plt.contour(pressures[i])
  text.set_text(vtks[i])
  return ax, text

ani = animation.FuncAnimation(fig, animate, range(len(pressures)), interval=500, 
                              blit=True, repeat_delay=2000)

cax = ax.imshow([[min_ for i in range(len(p))] for p in pressures[0]], 
                  norm=plt.Normalize(min_,max_))
text.set_text('')
plt.title('Reservoir Performance')
plt.colorbar(cax)
###############################################################################
#'''
# Reservoir Distribuition #####################################################
'''fig2 = plt.figure(1)
ax11 = fig2.add_subplot(121)
ax11.set_title('Reservoir Distribuition')
ax12 = fig2.add_subplot(122)
ax12.set_title('Max Pressure Performance')
ax12.set_xlim(-0.5, len(pressures)-0.5)
ax12.set_xticks(range(len(pressures)))
ax12.set_ylim(0.95*min(mxs), 1.05*max(mxs))
ax12.grid(True)
text2 = ax11.text(0, len(pressures[-1])//30+1, '')
'''
def ani_reservoir(i):
  ax11.imshow(pressures[i])
  print i, mxs[i]
  ax12.plot(range(i+1), mxs[:i+1], 'o-')
  text2.set_text(vtks[i])  
  return ax11, ax12, text2
  
#ani2 = animation.FuncAnimation(fig2, ani_reservoir, range(len(pressures)), 
#                               interval=1000, repeat_delay=2000)
                               
##############################################################################

plt.show()