 #!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt


fdict = {0:'matvec', 1:'ddot', 2:'daxpy', 3:'gen_B', 4:'save_vec', 
		 5:'MPI_Send', 6:'MPI_Recv', 7:'MPI_Init', 8:'MPI_Barrier', 9:'MPI_Finalize', 10:'MPI_Allreduce'}
Nfunc = len(fdict.keys())
data = np.loadtxt('tauprof.dat')

figname = 'tauprof.png'
fig = plt.figure(facecolor = 'w', edgecolor = 'w')
data = np.loadtxt('tauprof.dat', skiprows = 2)
pvals = data[:,0]

ind = 1
axlist = []
for key in range(Nfunc):
	val = fdict[key]
	x = pvals ; y = data[:, key+2]
	ax = fig.add_subplot(3,5,ind)
	ax.plot(x, y, linewidth = 3, marker = 'o', markersize = 8)
	ax.set_title(val)

	axlist.append(ax)
	ind += 1

plt.figtext(0.45, 0.02, 'number of processors', fontsize = 25)
plt.figtext(0.02, 0.6, 'total inclusive time (msec)', fontsize = 25, rotation = 90)
plt.tight_layout()
plt.savefig(figname)

plt.show()