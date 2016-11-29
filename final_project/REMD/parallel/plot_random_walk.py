#!/usr/bin/env python
import numpy as np
import os, sys
import matplotlib.pyplot as plt

statfile = sys.argv[1]
walkfile = sys.argv[2]

temps = []
with open(statfile, 'r') as of:
        for line in of: temps.append(int ( float(line.split()[0]) ) )

allwalk = np.loadtxt(walkfile)
inds = [0, 2]
clrs = {0: 'red', 2: 'blue'}

fig = plt.figure(figsize = (10,5), facecolor = 'w', edgecolor = 'w')
axs = []
ax_ind = 1
for i in inds:
        ax = fig.add_subplot(1,len(inds), ax_ind) ; axs.append(ax)
        walk = allwalk[:,i]
        time = range(len(walk))
        walk_temps = [temps[int(j)] for j in walk]
        ax.plot(time, walk_temps, 'k-', linewidth = 2, color = clrs[i], label = 'Replica %d' % i)
        ax_ind += 1

for ax in axs:
        ax.legend(loc = 'best', prop = {'size': 15})
        ax.set_xlabel('REX iteration', fontsize = 20)
        ax.set_ylabel('Temperature (K)', fontsize = 20)
        ax.set_yticks(temps)

plt.tight_layout()
plt.show()        
