#!/usr/bin/env python
import os, sys, subprocess, time

job_template = '''#!/bin/bash
#SBATCH --job-name=tauprof
#SBATCH --output="tauprof.out"
#SBATCH --partition=compute
#SBATCH --nodes=%(NNodes)d
#SBATCH --ntasks-per-node=%(NCoresPerNode)d
#SBATCH --export=ALL
#SBATCH -t 01:00:00
#SBATCH -A TG-ASC160059
#ibrun in verbose mode will give binding detail
cd /home/$USER/hw2
ibrun -np %(p)d -v ./cgsolve %(k)d 1 0
touch tauprof.done
'''

MaxCoresPerNode = 24
pList = [1, 2, 4, 8, 12, 24, 48, 72] # powers of 2 for less than 1 node and multiples of max cores per node for greater
kList = [ 24*64 ] * len(pList) # elapsed time for this k with p = 1 is ~ 32 sec

if os.path.isfile('xApprox.txt'): os.remove('xApprox.txt')
NCases = len(pList)
print 'Starting tau profiling experiment....\n'
for i in range(NCases):
	p = pList[i] ; k = kList[i]
	NNodes = int(p/MaxCoresPerNode)
	if NNodes < 1: NNodes = 1
	NCoresPerNode = p if p < MaxCoresPerNode else MaxCoresPerNode
	paramdict = {'NNodes': NNodes, 'NCoresPerNode': NCoresPerNode, 'k': k, 'p': p}
	jobfilename = 'tauprof_%d.sh' % p
	file(jobfilename, 'w').write(job_template % paramdict)

	# simple job completion monitor
	print 'p = %d job running...' % p
	proc = subprocess.Popen(['sbatch',  jobfilename], stdout = subprocess.PIPE)
	stdoutdata, stderrdata = proc.communicate()
	while not os.path.isfile('taupof.done'): time.sleep(1.0) # monitor
	print 'p = %d job completed\n' % p

	cmdstring = 'pprof > tauprof_%d.dat ; mv profile.* ./%d; rm tauprof.done; rm tauprof.out; rm xApprox.txt' % (p, p)
	os.system(cmdstring)