#!/usr/bin/env python
import numpy as np
import random, os
from mpi4py import MPI

class Replica(object):
	def __init__(self, Temp, rID = 0, NRexIter = 0, TempScaleRatio = 1.0, InitStateFile = None, isEquil = False, isTerm = False):
		self.rID = rID
		self.NRexIter = NRexIter

		self.Temp = Temp
		self.TempScaleRatio = TempScaleRatio
		self.Ene = None

		self.InitStateFile = InitStateFile
		self.isEquil = isEquil
		self.isTerm = False
		
		self.DataDir = os.path.join(os.getcwd(), str(self.rID))
		if not os.path.isdir(self.DataDir): os.mkdir(self.DataDir)

		self.StateFiles = []
		self.EneFiles = []

		self.kB = 0.001987 # Boltzmann constant
	
	def genFiles(self):
		if self.isEquil:
			self.NRexIter = -1
			if self.InitStateFile is None: self.InitStateFile = os.path.join(os.getcwd(), 'init.dat')
			self.CurrStateFile = os.path.join(self.DataDir, 'equil.%d.trj' % self. rID)
			self.CurrEneFile = os.path.join(self.DataDir, 'equil.%d.ene' % self.rID)	
		else:
			if self.InitStateFile is None: self.InitStateFile = os.path.join( self.DataDir, 'init.%d.%d.dat' %(self.rID, self.NRexIter) )
			self.CurrStateFile =  os.path.join( self.DataDir, '%d.%d.trj' % (self.rID, self.NRexIter) )
			self.CurrEneFile = os.path.join( self.DataDir, '%d.%d.ene' % (self.rID, self.NRexIter) )
			self.StateFiles.append(self.CurrStateFile) ; self.EneFiles.append(self.CurrEneFile)

		if not self.isTerm: self.NextInitStateFile = os.path.join(self.DataDir, 'init.%d.%d.dat' % (self.rID, self.NRexIter+1))


	def getEne(self):
		self.Ene = np.loadtxt(self.CurrEneFile)[-1]

	def getState(self, filename):
		pass
	
	def Run(self, RunSteps, StepFreq):
		# should write appropriate energy, state and initstate files
		# usually handled well in LAMMPS
		pass


class Ensemble(list):
	def __init__(self, SwapsPerCycle, Verbose = False, LogFile = None):
		list.__init__(self)
		
		self.SortedTempList = sorted([Replica.Temp for Replica in self])
		
		self.SwapFile = os.path.join(os.getcwd(), 'swap.txt')
		
		self.SwapsPerCycle = SwapsPerCycle
		self.allAttempts = np.zeros((len(self), len(self)), np.int32)
		self.accAttempts = np.zeros((len(self), len(self)), np.int32)

		self.Verbose = Verbose
		self.LogFile = LogFile

	def __Update(self):
		self.SortedTempList = sorted([Replica.Temp for Replica in self])
		self.allAttempts = np.zeros((len(self), len(self)), np.int32)
		self.accAttempts = np.zeros((len(self), len(self)), np.int32)

	def Log(self, msg):
		if not self.LogFile is None: file(self.LogFile, 'a').write(msg)
		if self.Verbose: print msg

	def MetHast(self, Replica1, Replica2):
		kB = Replica1.kB
		Delta = ( 1./(kB * Replica2.Temp) - 1./(kB * Replica1.Temp) ) * (Replica1.Ene - Replica2.Ene)
		pacc = min(1, np.exp(- Delta))
		return pacc > random.random()

	def getNeigh(self, i):
		if i == 0: return [1]
		elif i == len(self) - 1: return [len(self) - 2]
		else: return [i-1, i+1]

	def swap(self, NRexIter):
		self.Log('\n\n========== REX ITER %d ============\n\n' % NRexIter)
		
		self.__Update()

		acclist = []
		for  n in range(self.SwapsPerCycle):
				self.Log("Cycle: %d\n" % n)

				i = np.random.choice(len(self)) 
				neighs = self.getNeigh(i) 
				j = neighs[0] if len(neighs) == 1 else np.random.choice(neighs) 
				
				if acclist.__contains__((i,j)) or acclist.__contains__((j,i)): continue

				Replica_i = self[i] ; Replica_j = self[j]
				TInd_i = self.SortedTempList.index(Replica_i.Temp) ; TInd_j = self.SortedTempList.index(Replica_j.Temp)
				self.allAttempts[ min(TInd_i, TInd_j), max(TInd_i, TInd_j) ] += 1

				accept = self.MetHast(Replica_i, Replica_j)
				if accept:
					Replica_i.Temp, Replica_j.Temp = Replica_j.Temp, Replica_i.Temp
					Replica_i.TempScaleRatio = float(Replica_j.Temp) / float(Replica_i.Temp)
					Replica_j.TempScaleRatio = 1.0 / Replica_i.TempScaleRatio
					acclist.append((i, j))
					self.accAttempts[ min(TInd_i, TInd_j), max(TInd_i, TInd_j) ] += 1

	
	def save_permutation(self):
		TempList = np.array([Replica.Temp for Replica in self])
		permutation = np.argsort(TempList)
		print TempList, permutation
		s = '\t'.join([str(x) for x in permutation]) + '\n'
		with open(self.SwapFile, 'a') as of:
			of.write(s)


	def getStats(self):
		# this is written in serial since the data used is relatively small
		self.StatFile = os.path.join(os.getcwd(), 'statistics.txt')
		self.WalkFile = os.path.join(os.getcwd(), 'random_walk.txt')
		
		self.Log('\nCalculating acceptance ratios...')
		of = open(self.StatFile, 'w')
		for i in range(len(self)):
			for j in range(i+1, len(self)-1):
				Ti = self.SortedTempList[i]
				Tj = self.SortedTempList[j]
				acc_ratio = float(self.accAttempts[i,j]) / float(self.allAttempts[i,j]) if self.allAttempts[i,j] else 0.0
				of.write('%d K <--> %d K: %g\n' % (Ti, Tj, acc_ratio))

		of.close()

		self.Log('\nCalculating random walk in temperature space...')
		of = open(self.WalkFile, 'w')
		of.write('\t'.join( [str(x) for x in range(len(self))] ) + '\n')
		Walk = list(list(x) for x in np.loadtxt(self.SwapFile))
		for row in Walk:
			s = ''
			for i in range(len(self)):
				if row.__contains__(i): s += '\t%d' % row.index(i)
			s += '\n'
			of.write(s)

		of.close()




class REX(object):
	def __init__(self, ReplicaClass, Temps = None, EquilSteps = 2e5, ProdSteps = 2e5, StepFreq = 100, SwapSteps = 1000, SwapsPerCycle = 5, Verbose = False):
		self.comm = MPI.COMM_WORLD
		self.nproc = self.comm.size
		self.rank = self.comm.rank
		self.status = MPI.Status()

		self.EquilSteps = EquilSteps
		self.SwapSteps = int(SwapSteps)
		self.ProdSteps = ProdSteps / SwapSteps
		self.StepFreq = StepFreq
		self.SwapsPerCycle = SwapsPerCycle

		self.ReplicaClass = ReplicaClass
		self.Temps = Temps

		self.SIGACTIVE = 1
		self.SIGKILL = 0

		self.Verbose = Verbose
		self.LogFile = 'log.txt'

	def Log(self, msg):
		if not self.LogFile is None: file(self.LogFile, 'a').write(msg)
		if self.Verbose: print msg

	def Run(self):
		if self.rank == 0: self.Master()
		else: self.Slave()

	def Submit(self, sendbuf, destID):
		self.comm.send(sendbuf, dest = destID, tag = self.SIGACTIVE)
		rID =sendbuf[1]
		self.Log('\nStarting Replica %d on proc %d...\n' % (rID, destID))

	def Qsub(self, RunSteps):
		# set up flag arrays
		isRunning = []
		isQueued =range(len(self.server))

		# send off the first batch of jobs in order 
		for proc in range(1, self.nproc):
			rID = proc - 1
			sendbuf = (self.server[rID], rID, RunSteps)
			self.Submit(sendbuf, proc)
			isRunning.append(rID)
			isQueued.remove(rID)

		# set up a listener until all jobs are finished
		while isRunning:
			r, rID, proc = self.comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = self.status)
			self.Log('\nReplica %d on proc %d finished\n' % (rID, proc))
			self.server[rID] = r
			isRunning.remove(rID)
			if isQueued:
					rID = isQueued[0]
					sendbuf = (self.server[rID], rID, RunSteps)
					self.Submit(sendbuf, proc)
					isRunning.append(rID)
					isQueued.remove(rID)


	def Master(self):
		self.server = Ensemble(SwapsPerCycle = self.SwapsPerCycle, Verbose = self.Verbose, LogFile = self.LogFile)
		for i, Temp in enumerate(self.Temps): 
			r = self.ReplicaClass(Temp = Temp, rID = i)
			self.server.append(r)

		# equilbration
		self.Log('========== EQUILIBRATION ============\n\n')
		for r in self.server: r.isEquil = True
		self.Qsub(self.EquilSteps)

		# production with exchange
		self.Log('\n\n========== PRODUCTION ==============\n\n')
		for r in self.server: r.isEquil = False

		for iter in range(0, self.SwapSteps):
			for r in self.server:
				r.NRexIter = iter
				if iter == self.SwapSteps - 1: r.isTerm = True
				r.InitStateFile  = r.NextInitStateFile

			self.Qsub(self.ProdSteps)
			
			for r in self.server: r.getEne()

			self.server.save_permutation()

			self.server.swap(iter)

		self.server.getStats()

		for proc in range(1, self.nproc):
			sendbuf = (None, None, None)
			self.comm.send(sendbuf, dest = proc, tag = self.SIGKILL)


	def Slave(self):
		while True:
			r, rID, RunSteps = self.comm.recv(source = 0, tag = MPI.ANY_TAG, status = self.status)
			if self.status.Get_tag() == self.SIGKILL: break
			r.genFiles()
			r.Run(RunSteps, self.StepFreq)
			self.comm.send((r, rID, self.rank), dest = 0, tag = self.SIGACTIVE)




# This is a demuxer to order the traj by temperature
def demux(ReplicaClass, Temps, thisTemp, SwapSteps, SwapFile = 'swap.txt'):
		# this is written in parallel since the state and ene files for each replica are large

		comm = MPI.COMM_WORLD
		rank = comm.rank

		Ind = sorted(Temps).index(thisTemp)
		r = ReplicaClass(rID = Ind, Temp = thisTemp)

		MasterDir = os.path.join(os.getcwd(), str(thisTemp) + 'K')
		MasterStateFile = os.path.join(MasterDir, str(thisTemp) + '.trj')
		MasterEneFile = os.path.join(MasterDir, str(thisTemp) + '.ene')

		if rank == 0:
			if not os.path.isdir(MasterDir): os.mkdir(MasterDir)

		if rank: return # serial version for now

		Walk = list(np.loadtxt(SwapFile)[:, Ind])
		for NRexIter, i in enumerate(Walk):
			i = int(i) # numpy loads floats
			
			thisStateFile =  os.path.join('%d' % i, '%d.%d.trj' % (i, NRexIter))  
			s_state = r.getState(thisStateFile)
			with open(MasterStateFile, 'a') as of:
				of.write(s_state)

			thisEneFile =  os.path.join('%d' % i, '%d.%d.ene' % (i, NRexIter))
			s_ene = file(thisEneFile, 'r').read()
			with open(MasterEneFile, 'a') as of:
				of.write(s_ene)


# exponential schedule of temps.
def getTemps(Min, Max, NTemps):
	Min = float(Min) ; Max = float(Max)
	rate = (Max/Min) ** (1./NTemps)
	return [Min * (rate)**i for i in range(NTemps)]