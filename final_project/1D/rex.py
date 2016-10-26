#!/usr/bin/env python
import os, sys, random, copy
import numpy as np

sys.path.append(os.path.abspath('../')) ; from rexlib import *

# test function -- > E(x) = 3 * sin(x) + (0.1*x - 3)**2, 0 <=x <= 60.
# source : http://www.lindonslog.com/programming/stochastic-optimization-r-rmpi-parallel-tempering/
Efunc = lambda x:  0.6 * (3 * np.sin(x) + (0.1*x - 3)**2)

class myReplica(Replica):
	def Run(self, RunSteps, StepFreq):
		RunSteps = int(RunSteps)
		x = float(file(self.InitStateFile).read())
		Ene = Efunc(x)
		MAXDELTA = 5.0
		Box = [0, 60]; BoxL = Box[1] - Box[0] ; iBoxL = 1./float(BoxL)
		kB = 0.001987
		
		# write initial co-ordinate
		file(self.CurrStateFile, 'a').write('%g\n' % x)
		file(self.CurrEneFile, 'a').write('%g\n' % Ene)
		
		for n in range(RunSteps):
			x_new = x + random.uniform( -MAXDELTA, MAXDELTA )
			x_new -= BoxL * np.floor(x_new * iBoxL)
			Ene_new = Efunc(x_new)
			Delta = (Ene_new - Ene) / (kB * self.Temp)
			pacc = min(1.0, np.exp(-Delta))
			if pacc > random.random():
				x = x_new
				Ene = Ene_new
			
			if not n % StepFreq:
				file(self.CurrStateFile, 'a').write('%g\n' % x)
				file(self.CurrEneFile, 'a').write('%g\n' % Ene)

		if not self.isTerm: file(self.NextInitStateFile, 'a').write('%g\n' % x)

	def getState(self, filename):
		s = file(filename, 'r').read()
		return s




### MAIN
x0 = 30.0 ; file('init.dat', 'w').write('%g\n' % x0)
Temps = [270, 300, 350, 400, 450, 500, 600, 750]

rex = REX(ReplicaClass = myReplica, Temps = Temps, 
		  EquilSteps = 1e4, ProdSteps = 2e4, StepFreq = 10, SwapSteps = 100)
#rex.Run()

demux(ReplicaClass = myReplica, Temps = Temps, thisTemp = 300, SwapSteps = 100, SwapFile = 'swap.txt')

