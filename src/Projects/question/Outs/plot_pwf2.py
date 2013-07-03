# -*- coding: iso-8859-1 -*-
"""
Created on Wed Dec 26 18:39:58 2012

@author: Bismarck Gomes Souza Júnior
"""
import matplotlib.pyplot as plt
from numpy import *
import re, os
#import ade

ac = 5.614583
bc = 1.127

p_i = 7500.
q = 100.
u = 0.96
B = 0.92
k = 50/1000.
h = 10
phi = 0.3
ct = 9E-6
dx = 400.
  
def td(t):
  return ac*bc*k*t/(phi*u*ct*0.25**2)
  
def pd(p):
  return (2*pi*bc*k*h)/(q*B*u)*(p_i - p)
  
  
def df_dlogt(f,t):
  n = len(f)
  f_ = concatenate( ([f[0]],f,[f[-1]]) )
  t_ = concatenate( ([t[0]],t,[t[-1]]) )
  dfdt = []
  for i in range(n):
    dfdt.append( (f_[i+2]-f_[i])/(t_[i+2]-t_[i])*t_[i+1] )
  return array(dfdt[1:-1])
  
  
txts = filter(lambda x: x.endswith('.txt'), os.listdir('.'))
#txts = filter(lambda x: 'Peaceman' in x, txts)
txts.sort()
style = ['bo-', 'r^-', 'kv-', 'cs-']
names = []#[u'Numérica']
plot=[]

#
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
#plt.figure()

for txt in txts[::-1]:
  print txt,':',
  
  lines = open(txt,'r').readlines()[12:]

  while ('-'*30 not in lines[0]):
    lines.pop(0)

  lines[3] = re.sub(r'\[.+?\]', '', lines[3])

  tab_header = zip(lines[3].split(), ['']+lines[4].split())
  tab_matrix = array([[float(v) for v in l.split()] for l in lines[6:]])
  i = 0
  tab = {tab_header[0][0]: tab_matrix[:,i]}
  for title, pos in tab_header[1:]:
    i+=1
    if not tab.has_key(title): tab[title]={}
    tab[title][pos] = tab_matrix[:,i]
  
  time = tab['Time'][1:]
  if tab.has_key('Pwf'):
    n = 1
    items = tab['Pwf'].items()
    #max_pressure = max([max(item[1]) for item in items])
    #min_pressure = min([min(item[1]) for item in items])
    for pos, pressures in items:
      styl = style.pop(0) 
      plot.append(plt.loglog(td(time), pd(pressures[1:]), 'bo-')[0])
      dflogt = df_dlogt(pd(pressures[1:]),td(time))
      plt.loglog(td(time[1:-1]), dflogt , 'g^-')
      
      x = td(time)[td(time)<0.05]
      y = dflogt[:len(x)]
      #a,b = ade.NLR(x, y, ex ).var
      #print a
      #plt.loglog(x, 6.5E6*x, 'k-')
      #plt.loglog(x, 10**b*x**a)
    #ax.plot((time), (pressures[1:]), styl)
    #
'''
# Solução da linha fonte
from scipy.special import expi

def p(t, r=0.25):
  return  p_i+q*u*B/(2*pi*bc*k*h)*0.5*expi(1./(ac*bc)*phi*u*ct*r**2/(4*k*t))

#
ax.plot((time), (p(time)), 'b')

plot = [plt.loglog(td(time), pd(p(time)), 'b')[0]] + plot
plt.loglog(td(time[1:-1]), df_dlogt(pd(p(time[:])),td(time[:])), 'b')

plt.xlabel('$t_D$')
plt.ylabel(" "*20 + "${\partial Pwf_D}/{\partial t_D}$" + " "*20 + "$Pwf_D$")
plt.grid(True, which='both')
plt.legend(plot, [u'Analítica']+names, 'lower right')
plt.title(u'Modelo de Peaceman X Modelo de Ding et al')
plt.savefig('pseudo_estocagem_peacemanXding.png')


#######################################'''


plt.xlabel('$t_D$')
plt.ylabel(" "*20 + "${\partial Pwf_D}/{\partial t_D}$" + " "*20 + "$Pwf_D$")
plt.grid(True, which='both')
plt.show()







################################################################################
'''
def p2(t, r=0.1):
  p_i = 7500.
  q = 35*19.03
  u = 3.
  k = 100.
  bc = 1.127
  ac = 5.614583
  h = 4.
  phi = 0.2
  ct = 130E-6
  return q*u*1.25/(k*h)*0.5*expi(phi*u*ct*r**2/(4*k*24*0.0003484*t), 1)

print -3.122*expi(1./1.781/2407636./0.1)'''
