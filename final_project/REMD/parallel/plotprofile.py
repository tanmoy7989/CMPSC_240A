#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt

def get_time(filename):
    data = np.loadtxt(filename)
    t_server = []
    t_MD = []
    t_fetch = []
    p = []
    r = []
    for i, x in enumerate(data):
        r.append(x[0])
        p.append(x[1])
        t_server.append(x[2])
        t_MD.append(x[3])
        t_fetch.append(x[4])
    return (r, p, t_server, t_MD, t_fetch)

def plot(scaletype = 'rscale'):

    if scaletype == 'rscale':
        datafilename = os.path.abspath('./profile/rscale.txt')
        xlabel = 'number of replicas'
        figname = 'rscale.png'
    elif scaletype == 'pscale':
        datafilename = os.path.abspath('./profile/pscale.txt')
        xlabel = 'number of processors'
        figname = 'pscale.png'

    fig = plt.figure(figsize = (10, 5), facecolor = 'w', edgecolor = 'w')
    
    r, p, t_server, t_MD, t_fetch = get_time(datafilename)

    if scaletype == 'rscale':
        x = r
        lbl = 'N_proc = ' + str(int(p[0]))
    elif scaletype == 'pscale':
        x = p
        lbl = 'N_replica = ' + str(int(r[0]))
    
    speedup = np.zeros(len(x), np.float64)
    runtime = np.zeros(len(x), np.float64)
    for i in range(len(x)):
        runtime[i] = t_server[i] + t_MD[i] / float(p[i]) + t_fetch[i]
        speedup[i] = (t_server[i] + t_MD[i]) / (t_server[i] + t_MD[i] / float(p[i]) + t_fetch[i])


    ax1 = fig.add_subplot(1, 2, 1)
    ax1.semilogx(x, runtime, 'k-', basex = 2, linewidth = 3, marker = 'o', color = 'red', markersize = 8, label = lbl)
    ax1.set_xlabel(xlabel, fontsize = 20)
    ax1.set_ylabel('run time', fontsize = 20)
    ax1.grid(True)
    ax1.legend(loc = 'best', prop = {'size': 15})

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.semilogx(x, speedup, 'k-', basex = 2, linewidth = 3, marker = 'o', color = 'red', markersize = 8)
    ax2.set_xlabel(xlabel, fontsize = 20)
    ax2.set_ylabel('speedup', fontsize = 20)
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(figname)


#### MAIN
plot('rscale')
plot('pscale')
plt.show()
