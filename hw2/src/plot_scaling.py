 #!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

exptType = sys.argv[1]

data = np.loadtxt('%s.dat' % exptType)
kvals = data[:,0]; pvals = data[:,1]; tvals = data[:,2] ; t1 = tvals[0]
par_eff = t1/(tvals * pvals)

figname = '%s.png' % exptType
fig = plt.figure(figsize = (7,5), facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,2,1); ax2 = fig.add_subplot(1,2,2)
ax1.plot(pvals, tvals, 'k-', linewidth = 3, marker = 'o', markersize = 10)
ax2.plot(pvals, par_eff, linewidth = 3, marker = 'o', markersize = 10)
ax1.set_xlabel('number of processors', fontsize = 15); ax1.set_ylabel('elapsed time (sec)', fontsize = 15)
ax2.set_xlabel('number of processors', fontsize = 15); ax2.set_ylabel('parallel efficiency', fontsize = 15)
if exptType == 'strongscale': ax1.set_ylim([0.,0.3])
plt.tight_layout()
plt.savefig(figname)
plt.show()
