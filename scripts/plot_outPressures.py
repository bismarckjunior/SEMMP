# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Júnior
@email:  bismarckjunior@outlook.com
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

min_ = None
max_ = None
nullvalue = 0.5  # choose the color for null blocks, value between min and max
interval = 500


if os.path.lexists('Outs'):
    folder = 'Outs'
    print 'Reading "Outs" folder...\n'
else:
    folder = '.'
    print 'Reading current folder...\n'

vtks = filter(lambda x: x.endswith('.vtk'), os.listdir('.'))
vtks.sort()

if not vtks:
    print 'There is not output: outpressure (vtk format)'

pressures = []
mn, mx = [], []
for vtk in vtks:
    pressures.append([[float(p) for p in line.split()[:]] for line in
                      open(vtk).readlines()])
    mx.append(max(map(max, pressures[-1])))
    mn.append(min(map(lambda x: min(filter(lambda y: y != 0., x+[mx[-1]])),
                      pressures[-1])))


if min_ is None:
    min_ = min(mn)
if max_ is None:
    max_ = max(mx)

print 'min/max:', min_, max_

# Reservoir Performance #######################################################
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid(True)
ny = len(pressures[-1])
nx = len(pressures[-1][-1])
text = ax.text(nx//50, ny//40+1, '')
plt.xlim(0, nx)


def animate(i):
    p = [[p if p != 0.0 else min_+(max_-min_)*nullvalue for p in pr] for pr in
         pressures[i]]
    ax.imshow(p, norm=plt.Normalize(min_, max_))  # interpolation='nearest'
    text.set_text(vtks[i][-10:-4])
    return ax, text

ani = animation.FuncAnimation(fig, animate, range(len(pressures)),
                              interval=interval,
                              blit=True, repeat_delay=2000)

cax = ax.imshow([[min_ for i in range(len(p))] for p in pressures[0]],
                  norm=plt.Normalize(min_, max_))

plt.title('Reservoir Performance\n')
plt.colorbar(cax)

plt.show()
