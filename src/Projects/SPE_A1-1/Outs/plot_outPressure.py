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

fname = 'SPE_A1_1_000010.vtk'

min_, max_ = 1E10, 1E-5
data = [[float(p) for p in line.split()[:]] for line in open(fname).readlines()]
 
for dat in data:
  for d in dat:
    if 0<d<min_: min_=d
    elif d>max_: max_=d

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
pressures = np.array(pressures)
pressures[pressures==0] = 0

def animate(i):
  cax = ax.imshow(pressures[i], norm=plt.Normalize(min_,max_)) #interpolation='nearest'
  #plt.contour(pressures[i])
  text.set_text(vtks[i])
  return ax, text

ani = animation.FuncAnimation(fig, animate, range(len(pressures)), interval=1000, 
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
plt.figure()
plt.imshow( (pressures[0]-pressures[1])/pressures[1] )
#print (pressures[0]-pressures[1])/pressures[1] 
plt.title('Erro relativo ( Fluxo Incompressivel)')
plt.colorbar()
plt.savefig('erro_relativo_incompressivel.png')
plt.show()