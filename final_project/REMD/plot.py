#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt

kB = 0.001987
Temp = 300.0

whamfile = sys.argv[1]
binqueryfile = sys.argv[2]

enefile = sys.argv[3]
trjfile = sys.argv[4]
datacol = -1

Prefix = sys.argv[5]
makeMovie = int(sys.argv[6])

# sampling
whamdata = np.loadtxt(whamfile)

fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)
ax.plot(whamdata[:,0], whamdata[:,2], 'k-', linewidth = 3)
ax.set_xlabel(r'$R_g$', fontsize = 20); ax.set_ylabel(r'$\Delta F(R_g)$', fontsize = 20)
plt.tight_layout()
plt.savefig(Prefix + '.png')


# movie

tclscript = '''
set trajfile   [lindex $argv 0]
set start      0
set stop      -1
set freq       [lindex $argv 1]

set prefix     [lindex $argv 2]
	
mol new $trajfile type lammpstrj first $start last $stop step $freq waitfor all

color Display Background white
display projection orthographic
axes location Off
mol modstyle 0 0 VDW 1.0 42.0

set nframes [molinfo top get numframes]
for {set i 0} { $i <= $nframes} {incr i} {
	animate goto $i
	set outname [format {%s%03i} $prefix $i]
	render Tachyon ${outname}.dat
	/usr/local/lib/vmd/tachyon_LINUXAMD64 ${outname}.dat -o ${outname}.tga
}
quit
'''

freeenergy = lambda x : - kB * Temp * x

def renderVMD(framestep):
	print 'Rendering in vmd...'
	file('render.tcl', 'w').write(tclscript)
	vmdstring = 'vmd -dispdev text -e render.tcl -args %s %d %s' % (trjfile, framestep, 'posframe')
	os.system(vmdstring)
	os.remove('render.tcl')
	os.system('rm posframe*.dat')


if makeMovie:
	measure = np.loadtxt(enefile)[:, datacol]
	binquery = np.loadtxt(binqueryfile)

	fig = plt.figure(figsize = (10,5), facecolor = 'w', edgecolor = 'w')
	nframes = 300
	framestep = int(len(measure) / nframes)

	renderVMD(framestep)
	

	for i in range(nframes):
		
		print "Rendering snapshot", i

		ax1 = fig.add_subplot(1,2,1)
		im = plt.imread('posframe%03d.tga' % i)
		ax1.imshow(im)
		ax1.set_xticklabels([]) ; ax1.set_yticklabels([])
		ax1.hold(False)

		ax2 = fig.add_subplot(1,2,2)
		ax2.plot(whamdata[:,0], whamdata[:,2], 'k-', linewidth = 2)
		ax2.plot(whamdata[binquery[i*framestep], 0], whamdata[binquery[i*framestep], 2], marker = 'o', markersize = 8, color = 'red')
		if i == 0:
			ax2.set_xlabel(r'$R_g$', fontsize = 20); ax2.set_ylabel(r'$\Delta F(R_g)$', fontsize = 20)
		ax2.hold(True)
	
		plt.tight_layout()
		plt.savefig( 'frame_%03d.png' % i)

	cmdstring = 'avconv -i "frame_%03d.png" -r 10 -c:v libx264 -crf 10 -pix_fmt yuv420p ' + Prefix + '.avi'
	os.system(cmdstring)

	os.system('rm posframe*.tga frame*.png')


	
