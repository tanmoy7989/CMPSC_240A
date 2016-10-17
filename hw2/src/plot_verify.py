#!/usr/bin/env python
import numpy as np 
import os, sys
import matplotlib.pyplot as plt

# please just have 1 run data saved in the output file for this to work

xApproxFile = os.path.abspath(sys.argv[1])
l = file(xApproxFile, 'r').readlines()[1]
k = int(l.split()[0]); p = int(l.split()[1]);
x = np.loadtxt(xApproxFile, skiprows = 2).reshape((k,k))

x_expected = np.zeros(k*k)
for i in range(k*k):
	rowi = np.mod(i,k)
	x_expected[i] = float(rowi) / float(k+1)
x_expected = x_expected.reshape(k,k)

r = c = np.linspace(0,1,k)
R,C = np.meshgrid(r,c)
fig = plt.figure(figsize = (8,5), facecolor = 'w', edgecolor = 'w')
ax1 = fig.add_subplot(1,2,1); ax2 = fig.add_subplot(1,2,2)
ax1.contourf(C,R,x_expected) ; ax2.contourf(C,R,x)
ax1.set_title('Expected'); ax2.set_title('Calculated')
for ax in [ax1, ax2]:
	ax.set_xlabel(r'$x$', fontsize = 15)
	ax.set_ylabel(r'$y$', fontsize = 15)
	ax.grid(True)

plt.tight_layout()
plt.savefig('verify.png')
plt.show() 