#!/usr/bin/env python3

""" This plots the velocity profile after 1 timestep of MUSIC evolution
and compares it to a result from Mathematica."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
 

""" Read data from MUSIC results """

results_file = '../../evolution_xyeta.dat'

tvals = np.array([])
xvals = np.array([])
yvals = np.array([])
tempvals = np.array([])
evals = np.array([])
vvals = np.array([])
Tpropvals = np.array([])
pvals = np.array([])
vxvals = np.array([])
vyvals = np.array([])
vzvals = np.array([])

xmax = 15.0
ymax = 15.0
t0 = 0.4
tmax = t0+0.025
dx = 0.1
dy = 0.1
dt = 0.025*1


with open(results_file, 'r') as f:
    for t in np.arange(t0, tmax+dt, dt):
        if t > tmax:
            break
        print("Reading timestep: " + str(t) + ' fm')
        for y in np.arange(-ymax, ymax, dy):
            for x in np.arange(-xmax, xmax, dx):
                line = f.readline()
                tvals = np.append(tvals, t)
                xvals = np.append(xvals, x)
                yvals = np.append(yvals, y)
                tempvals = np.append(tempvals, 1000*float(line.split()[0]))
                vxvals = np.append(vxvals, float(line.split()[2]))
                vyvals = np.append(vyvals, float(line.split()[3]))
                vzvals = np.append(vzvals, float(line.split()[4]))
    f.close()
    
rvals = np.sqrt(xvals**2 + yvals**2)
vvals = np.sqrt(vxvals**2 + vyvals**2)

""" Filter to get one slice of y, t values """

grid_step = 0.15
step_fraction = 0.1
xtarget = 0.0
ytarget = 0.0
ttarget = 0.425

newyvals = yvals[abs(yvals-ytarget)<grid_step]
newxvals = xvals[abs(yvals-ytarget)<grid_step]
newtvals = tvals[abs(yvals-ytarget)<grid_step]
newtempvals = tempvals[abs(yvals-ytarget)<grid_step]
newvvals = vvals[abs(yvals-ytarget)<grid_step]
newvxvals = vxvals[abs(yvals-ytarget)<grid_step]
newvyvals = vyvals[abs(yvals-ytarget)<grid_step]
newvzvals = vzvals[abs(yvals-ytarget)<grid_step]

newtempvals = newtempvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvvals = newvvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvxvals = newvxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvyvals = newvyvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvzvals = newvzvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newxvals = newxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newtvals = newtvals[abs(newtvals-ttarget)<grid_step*step_fraction]

""" Plot against Mathematica data """

plt.rcdefaults()
plt.style.use(['seaborn-darkgrid', 'seaborn-deep', 'seaborn-notebook'])
plt.rcParams.update({
    'lines.linewidth': 1.5,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Lato'],
    'mathtext.fontset': 'custom',
    'mathtext.default': 'it',
    'mathtext.rm': 'sans',
    'mathtext.it': 'sans:italic:medium',
    'mathtext.cal': 'sans',
    'pdf.fonttype': 42,
})
plt.figure(figsize=(10,7))

""" 2D scatter plot of temperature (tempvals) over space and time: """
plt.ylabel(r'$v$')
plt.xlabel(r'$x$ (fm)')
#plt.scatter(rvals, tvals, c=tempvals, cmap=plt.cm.jet, vmin=0, vmax=400)
plt.plot(newxvals, newvvals, c='b', label = 'MUSIC')
mathematica_data = np.fromfile('mathematica_plot.dat').reshape(-1,2)
plt.plot(mathematica_data[:,0],mathematica_data[:,1], c='g', label = 'Mathematica')
plt.plot(-mathematica_data[:,0],mathematica_data[:,1], c='g')
plt.legend()
plt.title(r'$y = 0, \tau = 0.425 fm/c$')
plt.savefig('plots/VelocityComparison.png')
