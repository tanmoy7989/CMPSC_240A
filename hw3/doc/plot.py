#!/usr/bin/env python
import numpy as np 
import matplotlib.pyplot as plt

data = np.loadtxt('data', skiprows = 1)

rec = data[:,1] ; loop = data[:,2]; hyper = data[:,3]
x = data[:,0]

fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.semilogx(x, rec, basex = 10, color = 'red', linewidth = 3, marker = 'o', markersize = 7, label = 'rec_cilkified()')
ax.semilogx(x, loop, basex = 10, color = 'blue', linewidth = 3, marker = 'o', markersize = 7, label = 'loop_cilkified()')
ax.semilogx(x, hyper, basex = 10, color = 'green', linewidth = 3, marker = 'o', markersize = 7, label = 'hyperobject_cilkified()')

ax.legend(loc = 'best', prop = {'size': 15})
ax.set_xlabel('Size of array', fontsize = 20)
ax.set_ylabel('Parallel Efficiency', fontsize = 20)
ax.grid(True)

plt.tight_layout()
plt.show()