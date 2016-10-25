#!/usr/bin/env python
import numpy as np
import random, os
from mpi4py import MPI

class Replica(object):
	def __init__(self, RInd, NRexIter, Temp, TempScaleRatio = 1):
		self.RInd = RInd
		self.NRexIter = NRexIter

		self.Temp = Temp
		self.TempScaleRatio = TempScaleRatio
		
		self.DataDir = os.path.join(os.getcwd(), self.Name)
		if not os.path.isdir(self.DataDir): os.mkdir(self.DataDir)
		if self.NRexIter == 0:
			self.StateFiles = []
			self.EneFiles = []
			self.AvgEne = None

		this_StateFile =  os.path.join( self.DataDir, '%d.%d.trj.txt' % (self.RInd, self.NRexIter) )
		this_EneFile = os.path.join( self.DataDir, '%d.%d.ene.txt' % (self.RInd, self.NRexIter) )
		self.StateFiles.append(this_StateFile) ; self.EneFiles.append(this_EneFile)

	def getAvgEne(self):
		for this_Enefile in self.EneFiles:
			self.AvgEne += np.loadtxt(this_EneFile, skiprows = 1)[:,0]
		self.AvgEne = np.mean(self.AvgEne)

	def getState(self):
		pass
	
	def seqWork(self):
		pass


class Ensemble(list):
	def __init__(self, SwapsPerCycle = None):
		list.__init__(self)
		
		self.NReps = len(self)
		self.SortedTempList = sorted([Replica.Temp for Replica in self])
		
		self.SwapFile = os.path.join(os.getcwd(), 'state.txt')
		
		self.SwapsPerCycle = len(self) if SwapsPerCycle is None else SwapsPerCycle
		self.allAttempts = np.zeros((self.NReps, self.NReps), np.int32)
		self.accAttempts = np.zeros((self.NReps, self.NReps), np.int32)  

	def MetHast(self, Replica1, Replica2):
		kB = 1.0 # Boltzmann constant in kcal/mol
		Delta = ( 1./(kB * Replica2.Temp) - 1./(kB * Replica1.Temp) ) * (Replica1.AvgEne - Replica2.AvgEne)
		pacc = min(1, np.exp(- Delta))
		return pacc > random.random()

	def getNeigh(self, i):
		if i == 0: return [1]
		elif i == self.NReps - 1: return [self.NReps - 2]
		else: return [i-1, i+1]

	def swap(self):
		acclist = []
		for  n in range(self.SwapsPerCycle):
				i = np.random.choice(range(NReps)) 
				neighs = getNeigh(i) 
				j = neighs[0] if len(i) == 1 else np.random.choice(neighs) 
				
				if acclist.__contains__((i,j)) or acclist.__contains__((j,i)): continue

				Replica_i = self[i] ; Replica_j = self[j];
				TInd_i = self.SortedTempList.index(Replica_i.Temp) ; TInd_j = self.SortedTempList.index(Replica_j.Temp)
				self.allAttempts[ min(TInd_i, TInd_j), max(TInd_i, TInd_j) ] += 1

				accept = self.MetHast(Replica_i, Replica_j)
				if accept:
					Replica_i.Temp, Replica_j.Temp = Replica_j.Temp, Replica_i.Temp
					Replica_i.TempScaleRatio = float(Replica_j.Temp / Replica_i.Temp)
					Replica_j.TempScaleRatio = 1.0 / Replica_i.TempScaleRatio
					acclist.append((i, j))
					self.accAttempts[ min(TInd_i, TInd_j), max(TInd_i, TInd_j) ] += 1

	
	def save_permutation(self):
		TempList = np.array([Replica.Temp for Replica in self])
		permutation = np.argsort(TempList)
		s = '%t'.join([str(x) for x in permutation]) + '\n'
		with open(self.SwapFile, 'a') as of:
			of.write(s)

	def demux(self, TInd = 0):
		# this is written in parallel I/0 since the state and ene files for each replica are large
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


	def genStatistics(self):
		# this is written in serial since the data used is relatively small
		self.StatFile = os.path.join(os.getcwd(), 'statistics.txt')
		self.WalkFile = os.path.join(os.getcwd(), 'random_walk.txt')
			
		of = open(self.StatFile, 'w')
		for i in range(self.NReps):
			for j in range(i+1, self.NReps-1):
				Ti = self.SortedTempList[i]
				Tj = self.SortedTempList[j]
				acc_ratio = float(self.accAttempts[i,j]) / float(self.allAttempts[i,j])
				of.write('%d K <--> %d K: %g\n' % (Ti, Tj, acc_ratio))

		of.close()

		of = open(self.WalkFile, 'w')
		of.write('%\t'.join( [x for x in range(self.NReps)] ) + '\n')
		Walk = list(list(x) for x in np.loadtxt(self.SwapFile))
		for row in walk:
			s = ''
			for i in range(self.NReps):
				if row.__contains__(i): s += '\t%d' % row.index(i)
			s += '\n'
			of.write(s)

		of.close()

