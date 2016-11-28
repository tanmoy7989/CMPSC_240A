#!/usr/bin/env python
import numpy as np
import random, os, time
from mpi4py import MPI

class Replica(object):
    def __init__(self, Temp):
        self.Temp = Temp
        self.kB = 0.001987 # Boltzmann constant

    def getEne(self, filename):
        '''instantaneous energy from last sample of state space'''
        return np.loadtxt(filename)[-1]

    def getAllEne(self, filename):
        '''all energies for one sampling run'''
        with open(filename, 'r') as of:
            return of.read()

    def getState(self, filename):
        '''complete state space for one sampling run'''
        with open(filename, 'r') as of:
            return of.read()

    def Run(self, Files, RunSteps, StepFreq):
        '''the workhorse that runs the sampling algorithm
           can and will vary according to problem at hand'''
        pass



class REX(object):
    def __init__(self, ReplicaClass, Temps = None, EquilSteps = 2e5, ProdSteps = 2e5, StepFreq = 100, SwapFreq = 1000, SwapsPerCycle = 5,
                 DATADIR = os.getcwd(), INITFILE = 'init.dat', Verbose = False, Profile = False):
        # Initialize MPI
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.size
        self.rank = self.comm.rank
        self.status = MPI.Status()
        self.SIGACTIVE = 1
        self.SIGKILL = 0

        # Initialize serial code data
        self.ReplicaClass = ReplicaClass
        self.EquilSteps = int(EquilSteps)
        self.SwapFreq = int(SwapFreq)
        self.SwapSteps = int(ProdSteps) / self.SwapFreq
        self.ProdSteps = self.SwapFreq
        self.StepFreq = int(StepFreq)
        self.SwapsPerCycle = int(SwapsPerCycle)

        # Initialize rex data
        self.Ensemble = []
        self.Ene = []
        self.Temps = Temps
        self.niter = -1
        self.allAttempts = np.zeros(len(self.Temps))
        self.accAttempts = np.zeros(len(self.Temps))

        # Initialize logging
        self.Verbose = Verbose

        # Profiling
        self.Profile = Profile
        self.ServerTime = 0.0 ; self.RunTime = 0.0 ; self.FetchTime = 0.0

        # File system
        self.DATADIR = DATADIR
        self.INITFILE = INITFILE
        self.SwapFile = os.path.join(self.DATADIR, 'swap.txt')
        self.LogFile = os.path.join(self.DATADIR, 'log.txt')
        self.StatFile = os.path.join(self.DATADIR, 'statistics.txt')
        self.WalkFile = os.path.join(self.DATADIR, 'random_walk.txt')
        if self.Profile: self.ProfileFile = os.path.join(self.DATADIR, 'profile.txt')



    def createFileSys(self):
        '''create all necessary files and directories'''
        if not os.path.isdir(self.DATADIR): os.mkdir(self.DATADIR)

        with open(self.LogFile, 'w') as of: pass

        for rID in range(len(self.Temps)):
            for niter in range(self.SwapSteps):
                Dir = os.path.join(self.DATADIR, str(rID))
                if not os.path.isdir(Dir): os.mkdir(Dir)

    def get_rFiles(self, rID):
        '''get all necessary files needed for a replica to run its sampling code'''
        if self.niter == -1:
            initdatafile = self.INITFILE
            statefile = 'equil.%d.trj' % rID
            enefile = 'equil.%d.ene' % rID
            outputdatafile = '%d.%d.initdata' % (rID, 0)
        else:
            initdatafile = '%d.%d.initdata' % (rID, self.niter)
            statefile = '%d.%d.trj' % (rID, self.niter)
            enefile = '%d.%d.ene' % (rID, self.niter)
            outputdatafile = '%d.%d.initdata' % (rID, self.niter+1)

        files = []
        for thisfile in [initdatafile, statefile, enefile, outputdatafile]:
            Dir = "" if thisfile == self.INITFILE and self.niter == -1 else os.path.join(self.DATADIR, str(rID))
            files.append(os.path.join(Dir, thisfile))

        return files

    def demux(self, thisTemp):
        '''demultiplexer to concatenate all state files at a given temperature.
           This is written in parallel since the state and ene files for each replica may be large'''
        Dir = os.path.join(self.DATADIR, str(thisTemp) + 'K')
        statefile = os.path.join(Dir, str(thisTemp) + '.trj')
        enefile = os.path.join(Dir, str(thisTemp) + '.ene')

        # TODO: write this in parallel
        if self.rank == 0:
            self.Log('\nDemuxing %g trajectory...' % thisTemp)
            if not os.path.isdir(Dir): os.mkdir(Dir)
            for thisfile in [statefile, enefile]:
                with open(thisfile, 'w') as of: pass

            Ind = self.Temps.index(thisTemp)
            Walk = list(np.loadtxt(self.SwapFile)[:,Ind])
            dummy_r = self.ReplicaClass(Temp = None)
            self.niter = 0
            for rID in Walk:
                rID = int(rID) # since numpy.loadtxt returns a float
                #procID = rID - self.nproc * int(rID / self.nproc)
                #if self.rank == procID:
                rFiles = self.get_rFiles(rID)
                s_state = dummy_r.getState(rFiles[1])
                s_ene = dummy_r.getAllEne(rFiles[2])
                with open(statefile, 'a') as of: of.write(s_state)
                with open(enefile, 'a') as of: of.write(s_ene)
                self.niter += 1

    def Log(self, msg):
        '''screen and file logging'''
        if not self.LogFile is None: file(self.LogFile, 'a').write(msg)
        if self.Verbose: print msg

    def writeProfile(self):
        for x in [self.ServerTime, self.RunTime, self.FetchTime]:
            x /= self.SwapSteps
        file(self.ProfileFile, 'w').write('Mean Server Time = %g\nMean Run Time = %g\nMean Fetch Time = %g\n' % (self.ServerTime, self.RunTime, self.FetchTime))



    def Submit(self, rID, destID):
        '''send jobs to client processes'''
        r = self.Ensemble[rID]
        rFiles = self.get_rFiles(rID)
        RunSteps = self.EquilSteps if self.niter == -1 else self.ProdSteps
        sendbuf = (r, rID, rFiles, RunSteps)
        self.comm.send(sendbuf, dest = destID, tag = self.SIGACTIVE)
        self.Log('\nStarting Replica %d on proc %d...\n' % (rID, destID))

    def Qsub(self):
        '''manage a queue of jobs'''
        # set up flag arrays
        isRunning = []
        isQueued =range(len(self.Ensemble))

        # send off the first batch of jobs in order
        for proc in range(1, self.nproc):
            rID = proc - 1
            self.Submit(rID, proc)
            isRunning.append(rID)
            isQueued.remove(rID)

        # set up a listener until all jobs are finished
        while isRunning:
            rID, rEne, proc, clienttime = self.comm.recv(source = MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = self.status)
            self.Log('\nReplica %d on proc %d finished\n' % (rID, proc))
            self.Ene[rID] = rEne
            isRunning.remove(rID)

            if self.Profile:
                self.RunTime += clienttime[0]
                self.FetchTime += clienttime[1]

            if isQueued:
                    rID = isQueued[0]
                    self.Submit(rID, proc)
                    isRunning.append(rID)
                    isQueued.remove(rID)



    def MetHast(self, r1, r2, Ene1, Ene2):
        '''implement the Metropolis-Hastings criteria'''
        kB = r1.kB
        Delta = ( 1./(kB * r2.Temp) - 1./(kB * r1.Temp) ) * (Ene1 - Ene2)
        pacc = min(1, np.exp(- Delta))
        return pacc > random.random()

    def getNeigh(self, i):
        ''' get a random adjancent neighbor'''
        if i == 0: return 1
        elif i == len(self.Temps) - 1: return len(self.Temps) - 2
        else: return np.random.choice([i-1, i+1])

    def Swap(self):
        '''perform exchange between the replicas
           try self.SwapsPerCycle times'''

        self.Log('\n\n========== REX ITER %d ============\n\n' % self.niter)
        acc = []

        t1 = time.time() if self.Profile else None
        for  n in range(self.SwapsPerCycle):
                self.Log("Cycle: %d\n" % n)

                i = np.random.choice(len(self.Ensemble)) ; j = self.getNeigh(i)
                if acc and acc.__contains__((i,j)) or acc.__contains__((j,i)): continue
                acc.append((i, j))

                r1 = self.Ensemble[i] ; r2 = self.Ensemble[j]
                Ene1 = self.Ene[i] ; Ene2 = self.Ene[j]

                Ti = self.Temps.index(r1.Temp) ; Tj = self.Temps.index(r2.Temp)
                self.allAttempts[min(Ti, Tj)] += 1

                accept = self.MetHast(r1, r2, Ene1, Ene2)
                if accept:
                    r1.Temp, r2.Temp = r2.Temp, r1.Temp
                    self.accAttempts[min(Ti, Tj)] += 1

        t2 = time.time() if self.Profile else None
        if self.Profile: self.ServerTime += (t2 - t1)

    def save_permutation(self):
        '''save a permutation of the replicas in temperature space'''
        TempList = np.array([r.Temp for r in self.Ensemble])
        permutation = np.argsort(TempList)
        s = '\t'.join([str(x) for x in permutation]) + '\n'
        with open(self.SwapFile, 'a') as of: of.write(s)

    def getStats(self):
        '''generate statistics about the Metropolis Hastings acceptance
           and how a replica performs a random walk in temperature space.
           this is written in serial since the data used is relatively small'''
        self.Log('\nCalculating acceptance ratios...')
        with open(self.StatFile, 'w') as of:
            for i in range(len(self.Ensemble)-1):
                    Ti = self.Temps[i]
                    Tj = self.Temps[i+1]
                    acc_ratio = float(self.accAttempts[i]) / float(self.allAttempts[i]) if self.allAttempts[i] else 0.0
                    of.write('%g <--> %g : %g\n' % (Ti, Tj, acc_ratio))

        self.Log('\nCalculating random walk in temperature space...')
        with open(self.WalkFile, 'w') as of:
            Walk = list(list(x) for x in np.loadtxt(self.SwapFile))
            for row in Walk:
                s = ''
                for i in range(len(self.Ensemble)):
                    if row.__contains__(i): s += '\t%d' % row.index(i)
                s += '\n'
                of.write(s)



    def Run(self):
        '''wrapper'''
        if self.rank == 0: self.Parent()
        else: self.Child()

    def Parent(self):
        '''server that conducts all the replica exchanges'''

        self.createFileSys()

        for i in range(len(self.Temps)):
            self.Ensemble.append(self.ReplicaClass(self.Temps[i]))
            self.Ene.append(0.0)

        # equilbration
        self.Log('========== EQUILIBRATION ============\n\n')
        self.Qsub()

        # rex iteration
        self.niter = 0
        self.Log('\n\n========== PRODUCTION 0 ==============\n\n')
        while self.niter < self.SwapSteps:
            self.save_permutation()
            self.Qsub()
            self.Swap()
            self.niter += 1

        self.getStats()

        if self.Profile:
            self.RunTime /= self.SwapSteps
            self.FetchTime /= self.SwapSteps
            self.ServerTime /= self.Swapsteps
            self.writeProfile()

        # disconnect all clients
        for proc in range(1, self.nproc):
            sendbuf = (None, None, None, None)
            self.comm.send(sendbuf, dest = proc, tag = self.SIGKILL)

    def Child(self):
        '''client that runs the Run method of the replica'''
        while True:
            r, rID, rFiles, RunSteps = self.comm.recv(source = 0, tag = MPI.ANY_TAG, status = self.status)
            if self.status.Get_tag() == self.SIGKILL: break

            t1 = time.time() if self.Profile else None

            r.Run(rFiles, RunSteps, self.StepFreq)

            t2 = time.time() if self.Profile else None

            rEne = r.getEne(rFiles[2])

            t3 = time.time() if self.Profile else None

            if self.Profile:
                runtime = t2 - t1
                fetchtime = t3 - t2
                clienttime = [runtime, fetchtime]
            else: clienttime = None

            self.comm.send((rID, rEne, self.rank, clienttime), dest = 0, tag = self.SIGACTIVE)



# exponential schedule of temps.
def getTemps(Min, Max, NTemps):
    '''returns set of temps distributed exponentially'''
    Min = float(Min) ; Max = float(Max)
    rate = (Max/Min) ** (1./NTemps)
    return [Min * (rate)**i for i in range(NTemps)]
