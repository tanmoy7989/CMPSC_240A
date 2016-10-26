#!/usr/bin/env python
import numpy as np
import random, os
from mpi4py import MPI

class Replica(object):
	def __init__(self, Temp, RInd = 0, NRexIter = 0, TempScaleRatio = 1.0, InitStateFile = None, isEquil = False, isTerm = False):
		self.RInd = RInd
		self.NRexIter = NRexIter

		self.Temp = Temp
		self.TempScaleRatio = TempScaleRatio
		self.Ene = None

		self.InitStateFile = InitStateFile
		self.isEquil = isEquil
		self.isTerm = False
		
		self.DataDir = os.path.join(os.getcwd(), str(self.RInd))
		if not os.path.isdir(self.DataDir): os.mkdir(self.DataDir)

		self.StateFiles = []
		self.EneFiles = []
	
	def genFiles(self):
		if self.isEquil:
			self.NRexIter = -1
			if self.InitStateFile is None: self.InitStateFile = os.path.join(os.getcwd(), 'init.dat')
			self.CurrStateFile = os.path.join(self.DataDir, 'equil.%d.trj' % self. RInd)
			self.CurrEneFile = os.path.join(self.DataDir, 'equil.%d.ene' % self.RInd)	
		else:
			if self.InitStateFile is None: self.InitStateFile = os.path.join( self.DataDir, 'init.%d.%d.dat' %(self.RInd, self.NRexIter) )
			self.CurrStateFile =  os.path.join( self.DataDir, '%d.%d.trj' % (self.RInd, self.NRexIter) )
			self.CurrEneFile = os.path.join( self.DataDir, '%d.%d.ene' % (self.RInd, self.NRexIter) )
			self.StateFiles.append(self.CurrStateFile) ; self.EneFiles.append(self.CurrEneFile)

		if not self.isTerm: self.NextInitStateFile = os.path.join(self.DataDir, 'init.%d.%d.dat' % (self.RInd, self.NRexIter+1))


	def getEne(self):
		#TODO: check if this is the instantaneous energy or the trj average
		#for this_Enefile in self.EneFiles:
		#	self.Ene += np.loadtxt(this_EneFile, skiprows = 1)[:,0]
		#self.Ene = np.mean(self.Ene)

		# Instantaneous energy
		self.Ene = np.loadtxt(self.CurrEneFile)[-1]

	def getState(self, filename):
		pass
	
	def Run(self, RunSteps, StepFreq):
		# should write appropriate energy, state and initstate files
		# usually handled well in LAMMPS
		pass


class Ensemble(list):
	def __init__(self, SwapsPerCycle = None, Verbose = False):
		list.__init__(self)
		
		self.SortedTempList = sorted([Replica.Temp for Replica in self])
		
		self.SwapFile = os.path.join(os.getcwd(), 'swap.txt')
		
		self.SwapsPerCycle = len(self) if SwapsPerCycle is None else SwapsPerCycle
		self.allAttempts = np.zeros((len(self), len(self)), np.int32)
		self.accAttempts = np.zeros((len(self), len(self)), np.int32)

		self.Verbose = Verbose

	def __Update(self):
		self.SortedTempList = sorted([Replica.Temp for Replica in self])
		self.allAttempts = np.zeros((len(self), len(self)), np.int32)
		self.accAttempts = np.zeros((len(self), len(self)), np.int32)

	def MetHast(self, Replica1, Replica2):
		kB = 1.0 # Boltzmann constant in kcal/mol
		Delta = ( 1./(kB * Replica2.Temp) - 1./(kB * Replica1.Temp) ) * (Replica1.Ene - Replica2.Ene)
		pacc = min(1, np.exp(- Delta))
		return pacc > random.random()

	def getNeigh(self, i):
		if i == 0: return [1]
		elif i == len(self) - 1: return [len(self) - 2]
		else: return [i-1, i+1]

	def swap(self, NRexIter):
		if self.Verbose:
			s = "\nStarting REX ITER %d\n------------------------" % NRexIter
			print s
		
		self.__Update()

		acclist = []
		for  n in range(self.SwapsPerCycle):
				if self.Verbose: print "Cycle: %d" % n
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
		
		if self.Verbose: print "\nCalculating acceptance ratios..."	
		of = open(self.StatFile, 'w')
		for i in range(len(self)):
			for j in range(i+1, len(self)-1):
				Ti = self.SortedTempList[i]
				Tj = self.SortedTempList[j]
				acc_ratio = float(self.accAttempts[i,j]) / float(self.allAttempts[i,j]) if self.allAttempts[i,j] else 0.0
				of.write('%d K <--> %d K: %g\n' % (Ti, Tj, acc_ratio))

		of.close()

		if self.Verbose: print "\nCalculating random walk in temperature space..."
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
	def __init__(self, ReplicaClass, Temps = None, EquilSteps = 2e5, ProdSteps = 2e5, StepFreq = 100, SwapSteps = 1000):
		self.comm = MPI.COMM_WORLD
		self.nproc = self.comm.size
		self.rank = self.comm.rank
		self.status = MPI.Status()

		self.EquilSteps = EquilSteps
		self.SwapSteps = int(SwapSteps)
		self.ProdSteps = ProdSteps / SwapSteps
		self.StepFreq = StepFreq

		self.ReplicaClass = ReplicaClass
		self.Temps = Temps

		self.SIGACTIVE = 1
		self.SIGKILL = 0


	def Run(self):
		if self.rank == 0: self.Master()
		else: self.Slave()


	def Qsub(self, RunSteps):
		# set up flag arrays
		isRunning = []
		isQueued =range(len(self.server))

		# send off the first batch of jobs in order 
		for proc in range(1, self.nproc):
			rID = proc - 1
			sendbuf = (self.server[rID], rID, RunSteps)
			self.comm.send(sendbuf, dest = proc, tag = self.SIGACTIVE)
			isRunning.append(rID)
			isQueued.remove(rID)

		# set up a listener until all jobs are finished
		while isRunning:
			r, rID, proc = self.comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = self.status)
			print 'Replica %d on proc %d finished' % (rID, proc)
			self.server[rID] = r
			isRunning.remove(rID)
			if isQueued:
					rID = isQueued[0]
					sendbuf = (self.server[rID], rID, RunSteps)
					self.comm.send(sendbuf, dest = proc, tag = self.SIGACTIVE)
					isRunning.append(rID)
					isQueued.remove(rID)


	def Master(self):
		self.server = Ensemble(SwapsPerCycle = 5, Verbose = True)
		for i, Temp in enumerate(self.Temps): 
			r = self.ReplicaClass(Temp = Temp, RInd = i)
			self.server.append(r)
		
		# restart capabilities (backs up from 2nd last REX step) (#TODO: debug)
		startiter = -1 
		if os.path.isfile(self.server.SwapFile):
			data = np.loadtxt(self.server.SwapFile)
			startiter = len(data)
			for i, Temp in enumerate(data[-2,:]):
				self.server[i].Temp = Temp

			for iter in range(startiter):
				for i, r in enumerate(self.server):
					self.CurrStateFile =  os.path.join( self.DataDir, '%d.%d.trj' % (i, iter) )
					self.CurrEneFile = os.path.join( self.DataDir, '%d.%d.ene' % (i, iter) )
					self.StateFiles.append(self.CurrStateFile) ; self.EneFiles.append(self.CurrEneFile)

		# equilbration
		if startiter == -1:
			print 'Starting Equilibration....\n'
			for r in self.server: r.isEquil = True
			self.Qsub(self.EquilSteps)

		# production with exchange
		if startiter == -1: startiter = 0
		print '\nStarting Production....\n'
		for r in self.server: r.isEquil = False

		for iter in range(startiter, self.SwapSteps):
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
		r = ReplicaClass(RInd = Ind, Temp = thisTemp)

		MasterDir = os.path.join(os.getcwd(), str(thisTemp) + 'K')
		MasterStateFile = os.path.join(MasterDir, str(thisTemp) + '.trj')
		MasterEneFile = os.path.join(MasterDir, str(thisTemp) + '.ene')

		if rank == 0:
			if not os.path.isdir(MasterDir): os.mkdir(MasterDir)

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