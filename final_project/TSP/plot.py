#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, sys

coords = np.loadtxt('p01_xy.txt', skiprows = 1) 
dist = np.loadtxt('p01_d.txt', skiprows = 1)
ncities = len(coords)

def plotTour(ax, tour):
	ax.hold(True)
	for ii in range(ncities):
		i = tour[ii]
		p1 = [coords[i,0], coords[i,1]]

		if ii < ncities - 1:
			j = tour[ii+1]
			p2 = [coords[j,0], coords[j,1]]
			ax.quiver(p1[0], p1[1], p2[0]-p1[0], p2[1]-p1[1], scale_units = 'xy', angles = 'xy', scale = 1)

		color = 'red' if ii == 0 or ii == ncities - 1 else 'cyan'
		ax.scatter( [ p1[0] ], [ p1[1] ], s = 60, marker = 'o', color = color)


statefile = sys.argv[1]
enefile = sys.argv[2]
enesurf_filename = sys.argv[3]
makeMovie = bool(int(sys.argv[4]))
movie_filename = sys.argv[5]

state = np.loadtxt(statefile)[:10]
ene = np.loadtxt(enefile)[:10]

# sampling in distance space
fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.plot(ene, 'k-', linewidth = 3)
ax.set_xlabel('iteration', fontsize = 20); ax.set_ylabel('tour length', fontsize = 20)
plt.tight_layout()
plt.savefig(enesurf_filename)

# movie
if makeMovie:
	fig = plt.figure(figsize = (10,5), facecolor = 'w', edgecolor = 'w')
	framestep = 1
	for i in range(len(state)):
		if i % framestep: continue
	
		timeslice = range(i)
		stateslice = state[:i]
		eneslice = ene[:i]

		ax1 = fig.add_subplot(1,2,1)
		for tour in stateslice:
			plotTour(ax1, tour)
			ax1.hold(False)

		ax2 = fig.add_subplot(1,2,2)
		ax2.plot(timeslice, eneslice, linestyle = 'none', color = 'red', marker = 'o', markersize = 5)
		if i == 0:
			ax2.set_xlabel('iteration', fontsize = 20); ax2.set_ylabel('tour length', fontsize = 20)
		ax2.hold(True)
	
		plt.tight_layout()
		plt.savefig( 'frame_%d.png' % (int(i/framestep) + 1) )

	cmdstring = 'avconv -i "frame_%d.png" -r 25 -c:v libx264 -crf 10 -pix_fmt yuv420p ' + movie_filename
	os.system(cmdstring)

	for i in range(len(state)):
		if i % framestep: continue
		#os.remove( 'frame_%d.png' % (int(i/framestep) + 1) )

