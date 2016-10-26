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

	def genFiles(self):
		if self.isEquil:
			self.NRexIter = -1
			if self.InitStateFile is None: self.InitStateFile = os.path.join(os.getcwd(), 'init.dat')
			self.StateFiles = []
			self.EneFiles = []
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

	def MetHast(self, Replica1, Replica2):
		kB = 1.0 # Boltzmann constant in kcal/mol
		Delta = ( 1./(kB * Replica2.Temp) - 1./(kB * Replica1.Temp) ) * (Replica1.Ene - Replica2.Ene)
		pacc = min(1, np.exp(- Delta))
		return pacc > random.random()

	def getNeigh(self, i):
		if i == 0: return [1]
		elif i == len(self) - 1: return [len(self) - 2]
		else: return [i-1, i+1]

	def swap(self):
		if self.Verbose: print "Starting REX\n"
		acclist = []
		for  n in range(self.SwapsPerCycle):
				if Verbose: print "Swap Cycle: %d" % n 
				i = np.random.choice(len(self)) 
				neighs = getNeigh(i) 
				j = neighs[0] if len(i) == 1 else np.random.choice(neighs) 
				
				if acclist.__contains__((i,j)) or acclist.__contains__((j,i)): continue

				Replica_i = self[i] ; Replica_j = self[j];
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

	def demux(self, TInd = 0):
		# this is written in parallel since the state and ene files for each replica are large
		Temp = self.SortedTempList[TInd]
		MasterDir = os.path.join(os.getcwd(), str(Temp) + 'K')
		if not os.path.isdir(MasterDir): os.mkdir(MasterDir)
		MasterStateFile = os.path.join(MasterDir, str(Temp) + '.trj.txt')
		MasterEneFile = os.path.join(MasterDir, str(Temp) + '.ene.txt')

		comm = MPI.COMM_WORLD
		rank = comm.rank

		if rank == 0:
			Walk = list(np.loadtxt(self.SwapFile)[:, TInd])
			for NRexIter, i in enumerate(Walk):
				Replica_StateFiles = self[i].StateFiles
				Replica_EneFiles = self[i].EneFiles

				if (rank == i + 1):
					#TODO: get state and store it in s_state
					#with of as open(MasterStateFile, 'a'):
					#	of.write(s_state)

					s_ene = file(Replica_EneFiles[NRexIter], 'r').read()
					with open(MasterEneFile, 'a') as of:
						of.write(s_ene)

					comm.Barrier()


	def getStats(self):
		# this is written in serial since the data used is relatively small
		self.StatFile = os.path.join(os.getcwd(), 'statistics.txt')
		self.WalkFile = os.path.join(os.getcwd(), 'random_walk.txt')
		
		if self.Verbose: print "Calculating acceptance ratios..."	
		of = open(self.StatFile, 'w')
		for i in range(len(self)):
			for j in range(i+1, len(self)-1):
				Ti = self.SortedTempList[i]
				Tj = self.SortedTempList[j]
				acc_ratio = float(self.accAttempts[i,j]) / float(self.allAttempts[i,j])
				of.write('%d K <--> %d K: %g\n' % (Ti, Tj, acc_ratio))

		of.close()

		if self.Verbose: print "Calculating random walk in temperature space..."
		of = open(self.WalkFile, 'w')
		of.write('%\t'.join( [x for x in range(len(self))] ) + '\n')
		Walk = list(list(x) for x in np.loadtxt(self.SwapFile))
		for row in walk:
			s = ''
			for i in range(len(self)):
				if row.__contains__(i): s += '\t%d' % row.index(i)
			s += '\n'
			of.write(s)

		of.close()




class REX(object):
	def __init__(self, ReplicaClass, comm, Temps = None, EquilSteps = 2e5, ProdSteps = 2e5, StepFreq = 100, SwapSteps = 1000):
		self.comm = comm
		self.nproc = self.comm.size
		self.rank = self.comm.rank

		self.EquilSteps = EquilSteps
		self.SwapSteps = int(SwapSteps)
		self.ProdSteps = ProdSteps / SwapSteps
		self.StepFreq = StepFreq

		self.ReplicaClass = ReplicaClass
		self.Temps = Temps

		self.SIGINIT = 1
		self.SIGKILL = 0


	def Run(self):
		if self.rank == 0: self.Master()
		else: self.Slave()


	def Qsub(self, RunSteps):
		# set up flag arrays
		isRunning = []
		isQueued =range(len(self.server))

		status = MPI.Status()

		# send off the first batch of jobs in order 
		for proc in range(1, self.nproc):
			rID = proc - 1
			sendbuf = (self.server[rID], rID, RunSteps)
			self.comm.send(sendbuf, dest = proc, tag = self.SIGINIT)
			isRunning.append(rID)
			isQueued.remove(rID)

		# set up a listener until all jobs are finished
		while isRunning:
			rID, proc = self.comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = status)
			print 'Replica %d on proc %d finished' % (rID, proc)
			isRunning.remove(rID)
			if isQueued:
					rID = isQueued[0]
					sendbuf = (self.server[rID], rID, RunSteps)
					self.comm.send(sendbuf, dest = proc, tag = self.SIGINIT)
					isRunning.append(rID)
					isQueued.remove(rID)


		for proc in range(1, self.nproc):
			sendbuf = (None, None, None)
			self.comm.send(sendbuf, dest = proc, tag = self.SIGKILL)


	def Master(self):
		self.server = Ensemble(SwapsPerCycle = 5)
		for i, Temp in enumerate(self.Temps): 
			r = self.ReplicaClass(Temp = Temp, RInd = i)
			self.server.append(r)
		
		# equilbration
		for r in self.server: r.isEquil = True
		self.Qsub(self.EquilSteps)

		# production with exchange



	def Slave(self):
		status = MPI.Status()
		while True:
			r, rID, RunSteps = self.comm.recv(source = 0, tag = MPI.ANY_TAG, status = status)
			if status.Get_tag() == self.SIGINIT:
				print 'Starting replica %d on proc %d...' % (rID, self.rank)
				r.genFiles()
				r.Run(RunSteps, self.StepFreq)
				self.comm.send((rID, self.rank), dest = 0, tag = 1)

			elif status.Get_tag() == self.SIGKILL: break