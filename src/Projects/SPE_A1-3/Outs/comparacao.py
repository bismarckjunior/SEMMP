# -*- coding: iso-8859-1 -*-
"""
Created on Mon Dec 31 18:59:48 2012

@author: Bismarck Gomes Souza Júnior
"""
import matplotlib.pyplot as plt
from numpy import *

file_table = 'SPE_A1_3_comparacao_mmscf.txt'
table = [[float(d) for d in line.split()] for line in open(file_table, 'r').readlines()[1:]]
table = array(table)
#'''
time = table[:,0]
p44  = table[:,1]
p44n = table[:,2]
p39  = table[:,3]
p39n = table[:,4]
pwf44  = table[:,5]
pwf44n = table[:,6]
pwf39  = table[:,7]
pwf39n = table[:,8]
qw44  = table[:,9]
qw44n = table[:,10]
qw39  = table[:,11]
qw39n = table[:,12]
np  = table[:,13]
npn = table[:,14]
plt.figure()
plt.plot(time,p39,'--',time,p39n, 'o', time, p44, time, p44n, '^')
plt.title('Pressao no bloco (P)')
plt.xlabel('t [dias]')
plt.ylabel('psia')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.grid()
plt.savefig('SPE_A1_3_pressao_bloco.png')

plt.figure()
plt.plot(time, pwf39, '--', time, pwf39n, 'o', time, pwf44, time, pwf44n, '^')
plt.title('Pressao de fundo (Pwf)')
plt.xlabel('t [dias]')
plt.ylabel('psia')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.grid()
plt.savefig('SPE_A1_3_pressao_poco.png')

plt.figure()
plt.plot(time, qw39, '--', time, -qw39n, 'o', time, qw44, time, qw44n,'^')
plt.title(u'Vazão(Qw)')
plt.xlabel('t [dias]')
plt.ylabel('Mscf')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.ylim(2000, 12000)
plt.grid()
plt.savefig('SPE_A1_3_vazao.png')

plt.figure()
plt.plot(time, np/1000., time, npn/1000., 'o')
plt.title('Np')
plt.xlabel('t [dias]')
plt.ylabel('MMscf')
plt.legend(['Ertekin et al.', 'simulado'], 'best')
plt.grid()
plt.savefig('SPE_A1_3_prod_acum.png')

plt.figure()
e_np = abs(np-npn)/np*100
e_qw44 = abs(qw44-qw44n)/qw44*100
e_pwf39 = abs(pwf39-pwf39n)/pwf39*100
e_p39 = abs(p39-p39n)/p39*100
e_p44 = abs(p44-p44n)/p44*100
plt.plot(time, e_p39, time, e_p44,'s-',time, e_pwf39,'.-', time, e_qw44,'k--', time, e_np,'b:')
plt.title(u'Diferença Relativa Percentual')
plt.xlabel('t [dias]')
plt.ylabel('%')
plt.legend(['P: W-2', 'P: W-3', 'Pwf: W-2', 'Qw: W-3', 'Np'], 'best')
plt.grid()
#plt.ylim(-0.5, 6)
plt.savefig('SPE_A1_3_diferencas_relativas.png')

plt.show()#'''


################################################################################
a = 'out_compressible_000005.vtk'
b = 'ertekin_5dias.vtk'
#a = 'out_compressible_000100.vtk'
#b = 'ertekin_100dias.vtk'
simu = [[float(a) for a in li.split() ]for li in open(a, 'r').readlines()]
real = [[float(a) for a in li.split() ]for li in open(b, 'r').readlines()]
nx = len(simu)
ny = len(simu[-1])
X, Y = meshgrid(range(nx),range(ny))

################################################################################
erro = abs(array(simu)-array(real))/array(real)
print simu[0][0], real[0][0]
plt.imshow(erro*100)#, interpolation = 'nearest')
plt.colorbar()
plt.xticks(range(ny))
plt.yticks(range(nx))
plt.grid()
plt.title(u'Diferença percentual (100 dias)')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(9,3,'ko', markersize=8)
plt.plot(4,4, 'ko', markersize=8)
plt.legend(['Produtor'], 'lower left', numpoints=1)
#plt.savefig('SPE_A1_3_diferenca_relativa_100dias.png')
plt.show()