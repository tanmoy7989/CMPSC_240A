#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def get_observed(filename):
    data = np.loadtxt(filename)
    t1 = np.array([], float64)
    tp = np.array([], float64)
    p = np.array([], int32)
    r = np.array([], int32)
    for i,x in enumerate(data):
        r = np.append(r, x[0])
        p = np.append(p, x[1])
        t1 = np.append(t1, x[2] + x[3] + x[4])
        tp = np.append(tp, x[5] + x[6] + x[7])
        return (r, p, t1, tp)

def get_theoretical(filename):
    data = np.loadtxt(filename)
    t1 = np.array([], float64)
    tp = np.array([], float64)
    p = np.array([], int32)
    r = np.array([], int32)
    for i, x in enumerate(data):
        r = np.append(r, x[0])
        p = np.append(p, x[1])
        t1 = np.append(t1, x[5] + x[6])
        tp = np.append(tp, x[5] + (x[6] + x[7]) / float(x[1]) )
    return (r, p, t1, tp)

def plot(scaletype = 'rscale'):

    if scaletype == 'rscale':
        datafilename = 'rscale.txt'
        xlabel = 'number of replicas'
        figname = 'rscale.png'
    elif scaletype == 'pscale':
        datafilename = 'pscale.txt'
        xlabel = 'number of processors'
        figname = 'pscale.png'

    fig = plt.figure(figsize = (10, 5), facecolor = 'w', edgecolor = 'w')
    ro, po, t1o, tpo = get_observed(datafilename)
    rt, pt, t1t, tpt = get_theoretical(datafilename)

    if scaletype == 'rscale': x = (ro, rt)
    elif scaletype == 'pscale': x = (po, pt)

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.semilogx(x[0], tpo, basex = 2, 'ro', linestyle = 'none', markersize = 6, legend = 'observed')
    ax1.semilogx(x[1], tpt, basex = 2, 'k-', linedwidth = 3, legend = 'theoretical')
    ax1.set_xlabel(xlabel, fontsize = 20)
    ax1.set_ylabel('run time', fontsize = 20)
    ax1.legend(loc = 'best', prop = {'size': 15})
    ax1.grid(True)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.semilogx(x[0], t1o/tpo, basex = 2, 'ro', linestyle = 'none', markersize = 6, legend = 'observed')
    ax2.semilogx(x[1], t1t/tpt, basex = 2, 'k-', linewidth = 3, legend = 'theoretical')
    ax2.set_xlabel(xlabel, fontsize = 20)
    ax2.set_ylabel('speedup', fontsize = 20)
    ax2.legend(loc = 'best', prop = {'size': 15})
    ax2.grid(True)

    plt.savefig(figname)


#### MAIN
plot('rscale')
plot('pscale')
plt.show()
