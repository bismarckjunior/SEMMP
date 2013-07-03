# -*- coding: iso-8859-1 -*-
"""
Created on Mon Dec 31 18:59:48 2012

@author: Bismarck Gomes Souza Júnior
"""
import matplotlib.pyplot as plt
from numpy import *

file_table = 'SPE_A1_2_comparacao.txt'
table = [[float(d) for d in line.split()] for line in open(file_table, 'r').readlines()[1:]]
table = array(table)

time = table[:,0]
p39  = table[:,1]
p39n = table[:,2]
p44  = table[:,3]
p44n = table[:,4]
pwf39  = table[:,5]
pwf39n = table[:,6]
qw44  = table[:,7]
qw44n = table[:,8]
np  = table[:,9]
npn = table[:,10]
plt.figure()
plt.plot(time,p39,'--',time,p39n, 'o', time, p44, time, p44n, '^')
plt.title('Pressao no bloco (P)')
plt.xlabel('t [dias]')
plt.ylabel('psia')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.grid()
plt.savefig('SPE_A1_2_pressao_bloco.png')

plt.figure()
pwf = [5600]*len(time)
plt.plot(time, pwf39, '--', time, pwf39n, 'o', time, pwf, time, pwf, '^')
plt.title('Pressao de fundo (Pwf)')
plt.xlabel('t [dias]')
plt.ylabel('psia')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.grid()
plt.savefig('SPE_A1_2_pressao_poco.png')

plt.figure()
qw = [650]*len(time)
plt.plot(time, qw, '--', time, qw, 'o', time, qw44, time, qw44n,'^')
plt.title(u'Vazão(Qw)')
plt.xlabel('t [dias]')
plt.ylabel('STB/D')
plt.legend(['W-2: Ertekin et al.', 'W-2', 'W-3: Ertekin et al.', 'W-3'])
plt.ylim(-1000, 20000)
plt.grid()
plt.savefig('SPE_A1_2_vazao.png')

plt.figure()
plt.plot(time, np, time, npn, 'o')
plt.title('Np')
plt.xlabel('t [dias]')
plt.ylabel('STB')
plt.legend(['Ertekin et al.', 'simulado'], 'upper left')
plt.grid()
plt.savefig('SPE_A1_2_prod_acum.png')

plt.figure()
e_np = abs(np-npn)/np*100
e_qw44 = abs(qw44-qw44n)/qw44*100
e_pwf39 = abs(pwf39-pwf39n)/pwf39*100
e_p39 = abs(p39-p39n)/p39*100
e_p44 = abs(p44-p44n)/p44*100
plt.plot(time[:-1], e_p39[:-1], time[:-1], e_p44[:-1],time[:-1], e_pwf39[:-1],'.-', time[:-1], e_qw44[:-1],'k--', time[:-1], e_np[:-1],'b:')
plt.title(u'Diferença Relativa Percentual')
plt.xlabel('t [dias]')
plt.legend(['P(3,9)', 'P(4,4)', 'Pwf(3,9)', 'Qw(4,4)', 'Np'],'best')
plt.grid()
#plt.ylim(-0.5, 6)
plt.savefig('SPE_A1_2_diferencas_relativas.png')


plt.show()