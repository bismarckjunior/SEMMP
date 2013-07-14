# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Júnior
@email:  bismarckjunior@outlook.com
"""
import re
import os
import matplotlib.pyplot as plt
from numpy import array, sign, concatenate


key = 'Pwf'  # P, Pwf, Qw
logx = True
logy = True
dp_dlogt = True  # if true: logx=logy=True
style_ = ['bo-', 'g^-', 'rv-', 'ks-', 'c<-', 'm>-']*2  # max plots = 5*2 = 10

'''
#Dimensionless
ac = 5.614583
bc = 1.127

p_i = 7500.
q = 500.
rw=0.25
u = 0.96
B = 0.92
k = 10/1000.
h = 10
phi = 0.3
ct = 9E-6
dx = 50.

def td(t):
  return ac*bc*k*t/(phi*u*ct*rw**2)

def pd(p):
  return (2*pi*bc*k*h)/(q*B*u)*(p_i - p)
#'''


def df_dlogt(f, t):
    n = len(f)
    f_ = concatenate(([f[0]], f, [f[-1]]))
    t_ = concatenate(([t[0]], t, [t[-1]]))
    dfdt = []
    for i in range(n):
        dfdt.append((f_[i+2]-f_[i])/(t_[i+2]-t_[i])*t_[i+1])
    return array(dfdt[1:-1])


if os.path.lexists('Outs'):
    folder = 'Outs'
    print 'Reading "Outs" folder...\n'
else:
    folder = '.'
    print 'Reading current folder...\n'

txts = filter(lambda x: x.endswith('.txt'), os.listdir(folder))
txts.sort()

if not txts:
    print 'There is not output: report (txt format)'

for txt in txts:
    fig = plt.figure(txt[:-4])
    ax = fig.add_subplot(111)
    ax.set_title('$%s_{%s}$' % (key[0], key[1:]))
    ax.grid(True, which='both')
    ax.set_xlabel('$t$ [Day]')
    style = style_[:]
    if key == 'Pwf':
        if dp_dlogt:
            ylabel = "${\partial Pwf}/{\partial ln(t)}$"+" "*30+"$Pwf$"
        else:
            ylabel = '$P_{wf}$ [psia]'
    elif key == 'Qw':
        ylabel = '$Q_w$ [bbl/day]'
    elif key == 'P':
        if dp_dlogt:
            ylabel = "${\partial P}/{\partial ln(t)}$" + " "*30 + "$P$"
        else:
            ylabel = '$P$ [psia]'
    else:
        print '"key" does not make sense.'
        break
    ax.set_ylabel(ylabel)

    lines = open(folder+'/'+txt, 'r').readlines()[12:]

    while ('-'*30 not in lines[0]):
        lines.pop(0)

    lines[3] = re.sub(r'\[.+?\]', '', lines[3])

    tab_header = zip(lines[3].split(), ['']+lines[4].split())
    tab_matrix = array([[float(v) for v in l.split()] for l in lines[6:] if l])
    tab = {tab_header[0][0]: tab_matrix[:, 0]}

    for i, header in enumerate(tab_header[1:]):
        title, pos = header
        if not title in tab:
            tab[title] = {}
        tab[title][pos] = tab_matrix[:, i+1]

    time = tab['Time']
    plots = []
    if key in tab:
        legends = []
        for pos, pressures in tab[key].items():
            styl = style.pop(0)
            #Dimensionless
            #time = td(time)
            #pressures = pd(pressures)

            legends.append(key+' '+str(pos))

            if dp_dlogt and key != 'Qw':
                plots.append(plt.loglog(time, pressures, styl)[0])
                sgn = sign(sum(pressures[1:]-pressures[:-1]))
                plt.loglog(time[1:-1], df_dlogt(sgn*pressures, time), styl)
                continue

            if logx:
                plot = plt.loglog if logy else plt.semilogx
            else:
                plot = plt.semilogy if logy else plt.plot

            plots.append(plot(time, pressures, styl)[0])

        plt.legend(plots, legends, 'best')

plt.show()

'''
# Solução da linha fonte ######################################################
from scipy.special import expi

def p(t, r=0.25):
  return  p_i+q*u*B/(2*pi*bc*k*h)*0.5*expi(1./(ac*bc)*phi*u*ct*r**2/(4*k*t))

legends = [u'Analítica']+legends

ax.plot(time, p(time), 'b')
plot = [plt.loglog(td(time), pd(p(time)), 'b')[0]] + plot
plt.loglog(td(time[1:-1]), df_dlogt(pd(p(time)),td(time)), 'b')

plt.xlabel('$t_D$')
plt.ylabel(" "*20 + "${\partial Pwf_D}/{\partial ln(t_D)}$"+" "*20 +"$Pwf_D$")
plt.title(title)
#''

plt.grid(True, which='both')
plt.legend(plot, legends, 'best')
plt.savefig('%s.png' % title.replace(' ', '_'))
plt.show()
#'''
