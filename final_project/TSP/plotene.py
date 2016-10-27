#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, sys

coords = np.loadtxt('p01_xy.txt', skiprows = 1) 
dist = np.loadtxt('p01_d.txt', skiprows = 1)
ncities = len(coords)


def calcEne(tour):
	d = 0.0
	for ii in range(ncities):
		i = tour[ii] ; j = tour[ii+1]
		d += dist[i, j]
	return d


def plotTour(ax, tour):
	for ii in range(ncities + 1):
		i = tour[ii]
		p1 = [coords[i,0], coords[i,1]]
		
		color = 'red' if ii == 0 or ii == ncities else 'cyan'
		ax.scatter( [ p1[0] ], [ p1[1] ], s = 60, marker = 'o', color = color)

		if ii < ncities:
			j = tour[ii+1]
			p2 = [coords[j,0], coords[j,1]]
			ax.quiver(p1[0], p1[1], p2[0]-p1[0], p2[1]-p1[1], scale_units = 'xy', angles = 'xy', scale = 1)

		ax.hold(True)


fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)

s = [0, 12, 1, 14, 8, 4, 6, 2, 11, 13, 9, 7, 5, 3, 10, 0]
plotTour(ax, s)
print calcEne(s)


plt.show()

