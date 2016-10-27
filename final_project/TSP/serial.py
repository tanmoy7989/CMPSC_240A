#!/usr/bin/env python
import os, sys, random, copy
import numpy as np

sys.path.append(os.path.abspath('../')) ; import rexlib

coords = np.loadtxt('p01_xy.txt', skiprows = 1) 
dist = np.loadtxt('p01_d.txt', skiprows = 1)
ncities = len(coords)

def Efunc(tour):
	return sum( [ dist[tour[i], tour[i+1] ] for i in range(ncities) ] )

def tour2str(tour):
	return ' '.join(str(x) for x in tour) + '\n'

class TSP(rexlib.Replica):
	def Run(self, Files, RunSteps, StepFreq):
		self.kB = 1.0
		
		RunSteps = int(RunSteps) ; StepFreq = int(StepFreq)
		InitDataFile, StateFile, EneFile, OutputDataFile = Files

		tour = [int(x) for x in np.loadtxt(InitDataFile)]
		Ene = Efunc(tour)
		MAXSWAP = 5

		# write initial co-ordinate
		file(StateFile, 'a').write(tour2str(tour))
		file(EneFile, 'a').write('%g\n' % Ene)
		
		for n in range(RunSteps):
			
			tour_new = copy.copy(tour)
			isSwapped = []
			for swap in range(MAXSWAP):
				ii = np.random.choice(range(ncities-1))
				jj = np.random.choice([ii-1, ii + 1])
				if jj > ncities: jj = 0
				elif jj < 0: jj = ncities-1
				if isSwapped and ( isSwapped.__contains__((ii, jj)) or isSwapped.__contains__((jj, ii)) ): continue
				tour_new[ii], tour_new[jj] = tour_new[jj], tour_new[ii]
			
			Ene_new = Efunc(tour_new)
			Delta = (Ene_new - Ene) / (self.kB * self.Temp)
			pacc = min(1.0, np.exp(-Delta))
			if pacc > random.random():
				tour = tour_new
				Ene = Ene_new
			
			if not n % StepFreq:
				str_tour  = ' '.join(str(x) for x in tour) + '\n'
				file(StateFile, 'a').write(tour2str(tour))
				file(EneFile, 'a').write('%g\n' % Ene)

		file(OutputDataFile, 'w').write(tour2str(tour))
		for thisfile in Files: file(thisfile, 'r').close()


### MAIN
TempSet = 10.0
str_tour0 = tour2str(range(ncities) + [0]) ; file('init.dat', 'w').write(str_tour0)

EquilSteps = 1e3
ProdSteps = 2e3
StepFreq = 10

r = TSP(Temp = TempSet)

# equil runs
Files = ['init.dat', 'equil.trj', 'equil.ene', 'prod.dat']
r.Run(Files, EquilSteps, StepFreq)

# prod runs
Files = ['prod.dat', 'prod.trj', 'prod.ene', 'prod1.dat']

