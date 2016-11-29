#!/usr/bin/env python

import numpy as np 
import os, sys, pickle

kB = 0.001987
Temp = 300.0
scanFreq = 1

normalize = True

enefile = sys.argv[1]
Prefix = sys.argv[2]
datacol = -1

whamfile = Prefix + '.dat'
binqueryfile = Prefix + '.query'

def wham(nbins = 50):
	data = np.loadtxt(enefile)[:, datacol]
	bin_min = data.min() * 0.98
	bin_max = data.max() * 1.02

	hist = np.zeros((nbins, 3))
	bin_query = []

	delta = (bin_max - bin_min) / float(nbins)
	for i in range(nbins): hist[i,0] = bin_min + (i + 0.5) * delta

	FrameRange = range(0, len(data), scanFreq)
	NFrames = len(FrameRange)
	for frame in FrameRange:
		this_data = data[frame]
		bin_index = int( (this_data - bin_min) / delta )
		hist[bin_index, 1] += 1.0
		bin_query.append(bin_index)

	norm = np.trapz(x = hist[:,0], y = hist[:,1], dx = delta) if normalize else 1.0
	
	for i in range(nbins):
		hist[i,1] /= norm
		hist[i,2] = - kB * Temp * hist[i,1]

	np.savetxt(whamfile, hist, fmt = '%g')
	np.savetxt(binqueryfile, bin_query, fmt = '%d')


if not os.path.isfile(whamfile): wham()
