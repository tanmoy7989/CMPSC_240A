#!/usr/bin/env python

import os, sys

runJobs = True
collectProfiles = False

def writeSettingsFile(NReplicas, NCores):
    d = {
        'NReplicas': NReplicas,
        'NCores': NCores,
        'EquilSteps': 1000, #100000,
        'ProdSteps': 2000, #200000,
        'StepFreq': 5, #100,
        'SwapFreq': 100, #1000,
        'SwapsPerCycle': 5,
        'DataDir' : os.path.abspath('./profile/r%dp%d' % (NReplicas, NCores)),
        'InitFile': os.path.abspath('../init.data'),
    }

    filename = os.path.abspath('./profile/r%dp%d_settings.txt' % (NReplicas, NCores))
    file(filename, 'w').write(str(d))
    return filename

def writeJobScript(NReplicas, NCores, settingsfile):
    s = '''
#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N remd_%(r)d_%(p)d
# -pe ompi %(p)d
date
cd $HOME/projects/CMPSC_240A/final_project/REMD/parallel
mpirun -np %(p)d python rex.py %(settingsfile)s
'''
    d = {'r': NReplicas, 'p': NCores, 'settingsfile': settingsfile}
    scriptname = os.path.abspath( './remd_%d_%d.sh' % ( d['r'], d['p'] ) )
    file(scriptname, 'w').write(s % d)
    return scriptname

def readProfileFile(filename):
    with open(filename, 'r') as of:
        lines = of.readlines()
        t_s = float(lines[0][-1])
        t_MD = float(lines[1][-1])
        t_fetch = float(lines[2][-1])
    return (t_s, t_MD, t_fetch)



###### MAIN

if not os.path.isdir(os.path.abspath('./profile')): os.mkdir(os.path.abspath('./profile'))

rscale = {'NCores': 8, # use 1 full node
          'NReplicas': [4, 8, 12, 16], #[4, 8, 12, 16, 24, 32, 40, 48, 64, 96, 128],
          'filename': os.path.abspath('./rscale.txt') }

pscale = {'NReplicas': 20, #72,
          'NCores': [4, 8, 12, 16], #[4, 8, 12, 16, 24, 32, 36, 40, 48, 64, 72, 96],
          'filename': os.path.abspath('./pscale.txt') }

if runJobs:
    #1) scaling with number of replicas
    for r in rscale['NReplicas']:
        settingsfile_1 = writeSettingsFile(NReplicas = r, NCores = 1)
        jobscript_1 = writeJobScript(NReplicas = r, NCores = 1, settingsfile = settingsfile_1)

        settingsfile_p = writeSettingsFile(NReplicas = r, NCores = rscale['NCores'])
        jobscript_p = writeJobScript(NReplicas = r, NCores = rscale['NCores'], settingsfile = settingsfile_p)

        os.system('qsub ' + jobscript_1)
        os.system('qsub ' + jobscript_p)

    #2) scaling with number of cores
    settingsfile_1 = writeSettingsFile(NReplicas = pscale['NReplicas'], NCores = 1)
    jobscript_1 = writeJobScript(NReplicas = pscale['NReplicas'], NCores = 1, settingsfile = settingsfile_1)
    os.system('qsub ' + jobscript_1)

    for p in pscale['NCores']:
        settingsfile_p = writeSettingsFile(NReplicas = pscale['NReplicas'], NCores = p)
        jobscript_p = writeJobScript(NReplicas = pscale['NReplicas'], NCores = p, settingsfile = settingsfile_p)
        os.system('qsub ' + jobscript_p)


if collectProfiles:
    with open(rscale['filename'], 'w') as of:
        of.write('#NRep\tNProc\tt1_s\tt1_MD\tt1_fetch\tt_s\tt_MD\t_fetch\n')
        for r in rscale['NReplicas']:
            profilefile_1 = os.path.abspath('./profile/r%dp1/profile.txt' % r)
            t1_s, t1_MD, t1_fetch = readProfileFile(profilefile_1)

            profilefile_p  = os.path.abspath('./profile/r%dp%d/profile.txt' % (r, rscale['NCores']) )
            t_s, t_MD, t_fetch = readProfileFile(profilefile_p)

            of.write('%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n' % (r, rscale['NCores'], t1_s, t1_MD, t1_fetch, t_s, t_MD, t_fetch) )

    with open(pscale['filename'], 'w') as of:
        of.write('#NRep\tNProc\tt1_s\tt1_MD\tt1_fetch\tt_s\tt_MD\t_fetch\n')
        profilefile_1 = os.path.abspath('./profile/r%dp1/profile.txt' % pscale['NReplicas'])
        t1_s, t1_MD, t1_fetch = readProfileFile(profilefile_1)

        for p in rscale['NCores']:
            profilefile_p  = os.path.abspath('./profile/r%dp%d/profile.txt' % (pscale['NReplicas'], p) )
            t_s, t_MD, t_fetch = readProfileFile(profilefile_p)
            of.write('%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n' % (pscale['NReplicas'], p, t1_s, t1_MD, t1_fetch, t_s, t_MD, t_fetch) )
