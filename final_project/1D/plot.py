#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, sys

statefile = sys.argv[1]
enefile = sys.argv[2]
enesurf_filename = sys.argv[3]
makeMovie = bool(int(sys.argv[4]))
movie_filename = sys.argv[5]

pos = np.loadtxt(statefile)
ene = np.loadtxt(enefile)

# test function -- > E(x) = 3 * sin(x) + (0.1*x - 3)**2, 0 <=x <= 60.
# source : http://www.lindonslog.com/programming/stochastic-optimization-r-rmpi-parallel-tempering/
Efunc = lambda x:  0.6 * (3 * np.sin(x) + (0.1*x - 3)**2)
pos_known = np.linspace(0, 60, 100)
ene_known = np.array([Efunc(x) for x in pos_known])

# sampling
fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.plot(pos_known, ene_known, 'k-', linewidth = 3, label = 'Exact')
ax.plot(pos, ene, linestyle = 'none', color = 'red', marker = 'o', markersize = 5, label = 'Monte Carlo')
ax.set_xlabel(r'$x$', fontsize = 20); ax.set_ylabel(r'$E(x)$', fontsize = 20)
ax.legend(loc = 'best', prop = {'size': 15})
plt.tight_layout()
plt.savefig(enesurf_filename)

# movie
if makeMovie:
	fig = plt.figure(figsize = (10,5), facecolor = 'w', edgecolor = 'w')
	framestep = 10
	for i in range(len(pos)):
		if i % framestep: continue
	
		timeslice = range(i)
		posslice = pos[:i]
		eneslice = ene[:i]

		ax1 = fig.add_subplot(1,2,1)
		ax1.plot(np.linspace(0,60,10), np.zeros(10), linestyle = 'solid', color = 'black')
		ax1.plot(posslice, np.zeros(len(posslice)), linestyle = 'none', color = 'red', marker = 'o', markersize = 5)
		if i == 0: ax1.set_xlabel(r'$x$', fontsize = 20)
		ax1.hold(True)

		ax2 = fig.add_subplot(1,2,2)
		ax2.plot(pos_known, ene_known, 'k-', linewidth = 3, label = 'Exact')
		ax2.plot(posslice, eneslice, linestyle = 'none', color = 'red', marker = 'o', markersize = 5, label = 'Monte Carlo')
		if i == 0:
			ax2.set_xlabel(r'$x$', fontsize = 20); ax2.set_ylabel(r'$E(x)$', fontsize = 20)
			ax2.legend(loc = 'best', prop = {'size': 15})
		ax2.hold(True)
	
		plt.tight_layout()
		plt.savefig( 'frame_%d.png' % (int(i/framestep) + 1) )

	cmdstring = 'avconv -i "frame_%d.png" -r 25 -c:v libx264 -crf 10 -pix_fmt yuv420p ' + movie_filename
	os.system(cmdstring)

	for i in range(len(pos)):
		if i % framestep: continue
		os.remove( 'frame_%d.png' % (int(i/framestep) + 1) )

