#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, sys

coords = np.loadtxt('p01_xy.txt', skiprows = 1) 
dist = np.loadtxt('p01_d.txt', skiprows = 1)
ncities = len(coords)

def plotTour(ax, tour):
	del ax.collections[:]
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
Prefix = sys.argv[3]
makeMovie = bool(int(sys.argv[4]))

state = np.loadtxt(statefile)
ene = np.loadtxt(enefile)

# sampling in distance space
fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.plot(ene, 'k-', linewidth = 3)
ax.set_xlabel('iteration', fontsize = 20); ax.set_ylabel('tour length', fontsize = 20)
plt.tight_layout()
plt.savefig(Prefix+'.png')

# movie
if makeMovie:
	fig = plt.figure(figsize = (10,5), facecolor = 'w', edgecolor = 'w')
	nframes = 200
	framestep = int(len(state) / nframes)
	
	timeslice = []
	eneslice = []
	
	for i in range(nframes):
	
		print 'Rendering frame', i
		
		timeslice.append(i*framestep)
		tour = state[i*framestep]
		eneslice.append(ene[i*framestep])

		ax1 = fig.add_subplot(1,2,1)
	        plotTour(ax1, tour)

		ax2 = fig.add_subplot(1,2,2)
		if i == 0: ax2.set_xlabel('iteration', fontsize = 20); ax2.set_ylabel('tour length', fontsize = 20)
		ax2.plot(timeslice, eneslice, 'k-', linewidth = 3, color = 'red', marker = 'o', markersize = 8)
		ax2.hold(True)
	
		plt.tight_layout()
		plt.savefig( 'frame_%d.png' % i )

	cmdstring = 'avconv -i "frame_%d.png" -r 25 -c:v libx264 -crf 10 -pix_fmt yuv420p ' + Prefix + '.avi'
	os.system(cmdstring)

	os.system('rm frame*.png')

