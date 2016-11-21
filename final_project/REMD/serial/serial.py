#!/usr/bin/env python
import os, sys, subprocess, shutil
import numpy as np

sys.path.append(os.path.abspath('../..')) ; import rexlib

LammpsExec = os.path.expanduser('~/lammps-30Oct14/src/./lmp_mpi')
LammpsTemplate = os.path.abspath('../lammps_template.txt')

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
		
		p = subprocess.Popen('%s -in %s -log %s' % (LammpsExec, InFile, LogFile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		p.communicate()

		# Enefile is same as LogFile
		shutil.copyfile(LogFile, Files[2])

	def getState(self, filename):
		pass

	def getEne(self, filename):
		pass

	def getAllEne(self, filename):
		pass


### MAIN
TempSet = 300.0
EquilSteps = 1000000
ProdSteps = 20000000
StepFreq = 400
INITFILE = os.path.abspath('../init.data')

r = MD(Temp = TempSet)

# equil runs
Files = [INITFILE, 'equil.trj', 'equil.ene', 'prodinit.dat']
r.Run(Files, EquilSteps, StepFreq)

# prod runs
Files = ['prodinit.dat', 'prod.trj', 'prod.ene', 'prodout.dat']
r.Run(Files, ProdSteps, StepFreq)

