# -*- coding: iso-8859-1 -*-
"""
Created on Tue Jan 08 09:23:59 2013

@author: Bismarck Gomes Souza Júnior
"""

import matplotlib.pyplot as plt
from numpy import *

file_table = 'dirichlet_report005.txt'
table = [[float(d) for d in line.split()] for line in open(file_table, 'r').readlines()[47:]]
table = array(table)
style = ['o', 's', 'v', '^']

p = []
pn = []
k = 100/1000.   #Darcy
nx = ny = 50
dx = dy = 50
dt = 0.05
phi = 0.3
uB = 0.91#1.0#0.98245
c = (5.61458*1.127)*k/(phi*9E-6)/uB
D = (1./((dx*nx)**2) + 4./((dy*ny)**2))*pi*pi*c

#D = D*0.85
ex = exp(-D*dt)

legend = []
time = table[:,0]
for i in range(10):
  pn.append(table[:,i+1])
  x=y=i*5
  bloco = u'(%2i,%2i): ' % (x,y)
  legend.extend([bloco + u'Analítica', bloco+u'Numérica'])
  p.append([6500+1500*ex**(j)*sin(pi*(0.5+x)/nx)*sin(2*pi*(0.5+y)/ny) for j in range(len(time))])

a,b = 1,5

plt.figure()
for i in range(a,b):
  plt.plot(time, abs(array(p[i])), time, abs(array(pn[i])), style.pop())
  #plt.plot(time, abs(array(p[i])-array(pn[i]))/abs(array(p[i]))*100)
plt.title(u'Dirichlet: Análise da Queda de Pressão ($\Delta$t=%g)'%dt)
plt.xlabel('t [dias]')
plt.ylabel('p [psia]')
plt.legend(legend[2*a:2*b+1], 'best', numpoints=1)
plt.grid()
plt.hlines(6500, 0, time[-1], linestyles='dashed')
plt.xlim(0,1)
plt.savefig('dirichlet_comparacao.png')

print sum([abs(a-b) for a,b in zip(p[-1],pn[-1])])
plt.show()
