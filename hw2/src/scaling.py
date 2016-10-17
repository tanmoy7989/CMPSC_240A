#!/usr/bin/env python
import math, os, sys, subprocess, time
import numpy as np
import matplotlib.pyplot as plt

job_template = '''#!/bin/bash
#SBATCH --job-name=%(Prefix)s
#SBATCH --output="%(Prefix)s_%(p)d.out"
#SBATCH --partition=compute
#SBATCH --nodes=%(NNodes)d
#SBATCH --ntasks-per-node=%(NCoresPerNode)d
#SBATCH --export=ALL
#SBATCH -t 01:00:00
#SBATCH -A TG-ASC160059
#ibrun in verbose mode will give binding detail
cd /home/$USER/hw2
ibrun -np %(p)d -v ./cgsolve %(k)d 1 %(verifySol)d
touch %(Prefix)s_%(p)d.done
'''

isComputed = False
DelTempFiles = True

MaxCoresPerNode = 24

exptType = sys.argv[1]
if exptType == 'strongscale':
	Prefix = 'strongscale'
	pList = [1,2,4] #[1, 2, 4, 8, 12, 16, 24, 48, 72] # powers of 2 for less than 1 node and multiples of max cores per node for greater
	kList = [ 24*64 ] * len(pList) # elapsed time for this k with p = 1 is ~ 32 sec

if exptType == 'weakscale':
	# make sure that niters = 100 or some low value in main.c or this will take forever
	Prefix = 'weakscale'
	pList = [1, 2, 4, 8, 12, 16, 24, 48, 72]
	factor  = 24 * 64 # this maintains k^2/p = integer (can choose the LCM of p values, but will lead to too small runtimes)
	kList = [24 * 64]; [kList.append( factor * int(math.sqrt(p)) ) for p in pList]

if not isComputed:
	if os.path.isfile('xApprox.txt'): os.remove('xApprox.txt')
	if not os.path.isdir('.trash'): os.mkdir('.trash')
	NCases = len(pList)
	print 'Starting %s experiment....\n' % exptType
	for i in range(NCases):
		p = pList[i] ; k = kList[i]
		NNodes = int(p/MaxCoresPerNode)
		if NNodes < 1: NNodes = 1
		NCoresPerNode = p if p < MaxCoresPerNode else MaxCoresPerNode
		verifySol = 0 # don't verify solution in scaling runs (as this requires global vector assembly)
		paramdict = {'NNodes': NNodes, 'NCoresPerNode': NCoresPerNode, 'k': k, 'p': p, 'verifySol': verifySol, 'Prefix': Prefix}
		jobfilename = '%s_%d.sh' % (Prefix, p)
		file(jobfilename, 'w').write(job_template % paramdict)
	
		# simple job completion monitor
		print 'p = %d job running...' % p
		proc = subprocess.Popen(['sbatch',  jobfilename], stdout = subprocess.PIPE)
		stdoutdata, stderrdata = proc.communicate()
		while not os.path.isfile('%s_%d.done' % (Prefix, p) ): time.sleep(1.0) # monitor
		print 'p = %d job completed\n' % p
		
		# delete temporary files (jobscript, slurm output, completion token)
		if DelTempFiles:
			for filename in ['%s_%d.out', '%s_%d.sh', '%s_%d.done']:			
				if os.path.isfile(filename % (Prefix, p) ): os.system('mv ' + filename % (Prefix, p) + ' ./.trash')


# gather data when all experiments are complete
print 'Gathering time data...'
lines = file('xApprox.txt').readlines()
ind = [i for i, line in enumerate(lines) if line.startswith('#k')]
of = file('%s.dat' % exptType, 'w') ; of.write('k\tp\ttime elapsed\n')
for ii in ind:
	l = lines[ii+1].split()
	k = int(l[0]) ; p = int(l[1]) ; t = float(l[-1])
	of.write('%d\t%d\t%f\n' % (k, p, t))
of.close()

# calculate metrics and plot
figname = '%s.png' % exptType
data = np.loadtxt('%s.dat' % exptType)
kvals = data[:,0]; pvals = data[:,1]; tvals = data[:,2] ; t1 = tvals[0]
par_eff = t1/(pvals * tvals)
fig = plt.figure(figsize = (4,4), facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,2,1); ax2 = fig.add_subplot(1,2,2)
ax1.plot(pvals, tvals, 'k-', linewidth = 3, marker = 'o', markersize = 6)
ax2.plot(pvals, par_eff, linediwth = 3, marker = 'o', markersizer = 6)
ax1.set_xlabel('number of processors', fontsize = 10); ax1.set_ylabel('elapsed time (sec)', fontsize = 10)
ax2.set_xlabel('number of processors', fontsize = 10); ax2.set_ylabel('parallel efficiency', fontsize = 10)
plt.savefig(figname)

plt.show()

