#!/usr/bin/env python
import os, sys, subprocess
import numpy as np

sys.path.append(os.path.abspath('../')) ; import rexlib

LammpsExec = 'lmp'
LammpsTemplate = os.path.abspath('./lammps_template.txt')
INITFILE = os.path.abspath('./init.dat')


class MD(rexlib.Replica):
	def Run(self, Files, RunSteps, StepFreq):
		self.kB = 0.001987

		Prefix = Files[1].split('.trj')[0]

		d = {'InitDataFile': Files[0], 'TrjFile': Files[1], 'EneFile': Files[2], 'OutPutDataFile': Files[3],
			 'RunSteps': int(RunSteps), 'StepFreq': int(StepFreq), 'Temp': self.Temp}

		InputFile = Prefix + '.in'
		LogFile = Prefix + '.log'
		file(InputFile, 'w').write(LammpsTemplate % d)
		cmd = '%s -in %s -log %s' % (LammpsExec, InputFile, LogFile)
		p = subprocess.Popen(cmd.split(), shell = True) # TODO: make sure this works
		retcode, reterr = p.communicate()

	def getState(self, filename):
		pass

	def getEne(self, filename):
		return np.loadtxt(filename)[-1, 1]

	def getAllEne(self, filename):
		pass


### MAIN
Temps = [300, 400, 500, 600, 700]

EquilSteps = 1000000
ProdSteps = 20000000
StepFreq = 400
SwapSteps = 1000
SwapsPerCycle = 5
DATADIR = os.path.abspath('./test')

rex = rexlib.REX(ReplicaClass = MD, Temps = Temps, EquilSteps = EquilSteps, ProdSteps = ProdSteps, StepFreq = StepFreq, SwapSteps = SwapSteps, SwapsPerCycle = SwapsPerCycle, 
				 DATADIR = DATADIR, Verbose = True)
rex.Run()
rex.demux(300)

