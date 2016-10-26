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
		for n in range(RunSteps):
			if not n % StepFreq: print 'Step: %d\tx: %g\tEne: %g' % (n, x, Ene)
			
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
		return np.loadtxt(filename)


### MAIN
TempSet = 300.0
x0 = 30.0 ; file('init.dat', 'w').write('%g\n' % x0)
EquilSteps = 2e4
ProdSteps = 2e4
StepFreq = 10

r = myReplica(RInd = 0, Temp = TempSet)

# equil runs
r.isEquil = True
r.isTerm = False
r.genFiles()
r.Run(EquilSteps, StepFreq)

# prod runs
r.isEquil = False
r.isTerm = True
r.NRexIter = 0
r.genFiles()
r.Run(ProdSteps, StepFreq)