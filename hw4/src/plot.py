#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

datafile = sys.argv[1]
scaletype = sys.argv[2]

data = np.loadtxt(datafile)
x = data[:,0]
t1 = data[:,1]
tp = data[:,1]
teps = data[:, -1]

if scaletype == 'scalecore':  pe = t1 / (x * tp)
if scaletype == 'scalegraph': pe = t1 / (48 * tp)

xlabels = {'scalecore': 'number of threads', 'scalegraph': 'size of torus graph'}

fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (10,5))

ax1 = fig.add_subplot(1,2,1)
ax1.semilogx(x, pe, 'k-', marker = 'o', linewidth = 3, markersize = 6, basex = 2)
ax1.set_xlabel(xlabels[scaletype], fontsize = 20)
ax1.set_ylabel('parallel efficiency', fontsize = 20)
ax1.grid(True)

ax2 = fig.add_subplot(1,2,2)
ax2.loglog(x, teps, 'k-', marker = 'o', linewidth = 3, markersize = 6, basex = 2, basey = 10)
ax2.set_xlabel(xlabels[scaletype], fontsize = 20)
ax2.set_ylabel('TEPS score', fontsize = 20)
ax2.grid(True)

plt.tight_layout()
plt.savefig(scaletype+'.png')
plt.show()