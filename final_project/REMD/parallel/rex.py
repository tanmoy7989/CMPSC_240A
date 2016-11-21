#!/usr/bin/env python
import os, sys, subprocess, time, copy, shutil
import numpy as np
from mpi4py import MPI

sys.path.append(os.path.abspath('../..')) ; import rexlib

LammpsExec = os.path.expanduser('~/lammps-30Oct14/src/lmp_mpi')
LammpsTemplate = os.path.abspath('../lammps_template.txt')

STARTTOKEN = 'Step Temp Atoms PotEng rg'
STOPTOKEN = 'Loop time'

MPI_TYPE = 'ompi'
MPI_ENV_TOKENS = {'ompi': ['MPI'],
				  'mvapich2': ['MV2_', 'MPI', 'PAMID', 'PMI']
				  }


def setNewEnv(currenv):
	newenv = {}
	tokens = MPI_ENV_TOKENS[MPI_TYPE]
	for k,v in currenv.iteritems():
		remove = [k.__contains__(x) for x in tokens]
		if any(remove): continue
		else: newenv[k] = v
	return newenv
 

class MD(rexlib.Replica):

	def Run(self, Files, RunSteps, StepFreq):
		self.kB = 0.001987

		Prefix = Files[1].split('.trj')[0]

		InFile = Prefix + '.in'
		LogFile = Prefix + '.log'
		
		d = {'InitDataFile': Files[0], 'TrjFile': Files[1], 'EneFile': Files[2], 'OutputDataFile': Files[3],
			 'RunSteps': int(RunSteps), 'StepFreq': int(StepFreq), 'Temp': self.Temp}
		
		s = file(LammpsTemplate).read()
		file(InFile, 'w').write(s % d)

		# back up the current MPI environment and create a new one (since LAMMPS will call MPI_Init)
		newenv = setNewEnv(copy.copy(os.environ))

		# spawn a system call to LAMMPS
		p = subprocess.Popen('%s -in %s -log %s' % (LammpsExec, InFile, LogFile), shell = True, env = newenv, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		p.communicate()

		# Enefile is same as LogFile
		shutil.copyfile(LogFile, Files[2])

	def getState(self, filename):
		return file(filename).read()

	def getEne(self, filename):
		with open(filename, 'r') as of:
			lines = of.readlines()
			stop = [lines.index(line) for line in lines if line.startswith(STOPTOKEN)][0]
			l = lines[stop - 1].split()
		return float(l[-2])

	def getAllEne(self, filename):
		with open(filename, 'r') as of:
			lines = of.readlines()
			start = [lines.index(line) for line in lines if line.startswith(STARTTOKEN)][0] + 1
			stop = [lines.index(line) for line in lines if line.startswith(STOPTOKEN)][0]
		ret = ''.join(lines[start:stop])
		return ret


### MAIN
Temps = [270, 300, 400, 500, 600]

EquilSteps = 1000 #1000000
ProdSteps = 2000 #20000000
StepFreq = 10 #400
SwapSteps = 100
SwapsPerCycle = 1 #5
DATADIR = os.path.abspath('./test')
INITFILE = os.path.abspath('../init.data')

rex = rexlib.REX(ReplicaClass = MD, Temps = Temps, EquilSteps = EquilSteps, ProdSteps = ProdSteps, StepFreq = StepFreq, SwapSteps = SwapSteps, SwapsPerCycle = SwapsPerCycle, 
				 DATADIR = DATADIR, INITFILE = INITFILE, Verbose = True)
rex.Run()

rex.demux(300)

