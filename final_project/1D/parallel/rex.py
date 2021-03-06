#!/usr/bin/env python
import os, sys, random, copy
import numpy as np

sys.path.append(os.path.abspath('../..')) ; import rexlib

# test function -- > E(x) = 3 * sin(x) + (0.1*x - 3)**2, 0 <=x <= 60.
# source : http://www.lindonslog.com/programming/stochastic-optimization-r-rmpi-parallel-tempering/
Efunc = lambda x:  0.6 * (3 * np.sin(x) + (0.1*x - 3)**2)

class myReplica(rexlib.Replica):
	def Run(self, Files, RunSteps, StepFreq):
		self.kB = 0.001987
		
		RunSteps = int(RunSteps) ; StepFreq = int(StepFreq)
		InitDataFile, StateFile, EneFile, OutputDataFile = Files

		x = float(file(InitDataFile).read())
		Ene = Efunc(x)
		MAXDELTA = 5.0
		Box = [0, 60]; BoxL = Box[1] - Box[0] ; iBoxL = 1./float(BoxL)
		
		# write initial co-ordinate
		file(StateFile, 'a').write('%g\n' % x)
		file(EneFile, 'a').write('%g\n' % Ene)
		
		for n in range(RunSteps):
			x_new = x + random.uniform( -MAXDELTA, MAXDELTA )
			x_new -= BoxL * np.floor(x_new * iBoxL)
			Ene_new = Efunc(x_new)
			Delta = (Ene_new - Ene) / (self.kB * self.Temp)
			pacc = min(1.0, np.exp(-Delta))
			if pacc > random.random():
				x = x_new
				Ene = Ene_new
			
			if not n % StepFreq:
				file(StateFile, 'a').write('%g\n' % x)
				file(EneFile, 'a').write('%g\n' % Ene)

		file(OutputDataFile, 'w').write('%g\n' % x)

		for thisfile in Files: file(thisfile, 'r').close()


### MAIN
INITFILE = os.path.abspath('../init.dat')
DATADIR = os.path.abspath('./test')
Temps = rexlib.getTemps(270, 1000, 24) ; print Temps
rex = rexlib.REX(ReplicaClass = myReplica, Temps = Temps, EquilSteps = 10000, ProdSteps = 20000, StepFreq = 1, SwapFreq = 100, SwapsPerCycle = 2,
				 Verbose = False, Profile = True, DATADIR = DATADIR, INITFILE = INITFILE)

rex.Run()

rex.demux(3)

