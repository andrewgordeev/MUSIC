#!/usr/bin/env python3

""" This plots the velocity profile after 1 timestep of MUSIC evolution
and compares it to a result from Mathematica."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
 

""" Read data from MUSIC results """

results_file = '../../evolution_xyeta.dat'

xmax = 15.0
ymax = 15.0
t0 = 0.4
tmax = t0+0.0025
dx = 0.1
dy = 0.1
dt = 0.0025
nx = int(2*xmax/dx)
ny = int(2*ymax/dy)
nt = int((tmax-t0)/dt + 1)

tvals = np.array([])
yvals = np.array([])

for t in np.arange(t0, tmax+dt, dt):
    if t > tmax:
        break
    tvals = np.append(tvals, np.repeat(t,nx*ny))

for y in np.arange(-ymax, ymax, dy):
    yvals = np.append(yvals, np.repeat(y, nx))
yvals = np.tile(yvals, nt)
        
xvals = np.tile(np.arange(-xmax,xmax,dx),ny*nt)

results = np.loadtxt(results_file)
tempvals = 1000 * results[:,0]
vxvals = results[:,2]
vyvals = results[:,3]
vzvals = results[:,4]
#evals = results[:,5]
   
rvals = np.sqrt(xvals**2 + yvals**2)
vvals = np.sqrt(vxvals**2 + vyvals**2)

""" For MUSIC mode 666 """

# results_file = '../../../music_mod/evolution_bjorken.dat'
# tvals = results[:,0]
# xvals = results[:,1]
# evals = results[:,3]
# tempvals = results[:,4]
# uxvals = results[:,5]
# uetavals = results[:,6]
# utvals = np.sqrt(1+uxvals**2)
# vxvals = np.abs(uxvals/utvals)

""" Filter to get one slice of y, t values """

grid_step = 0.000001
step_fraction = 0.00001
xtarget = 0.0
ytarget = 0.0
ttarget = 0.4025

newyvals = yvals[abs(yvals-ytarget)<grid_step]
newxvals = xvals[abs(yvals-ytarget)<grid_step]
newtvals = tvals[abs(yvals-ytarget)<grid_step]
newtempvals = tempvals[abs(yvals-ytarget)<grid_step]
newvvals = vvals[abs(yvals-ytarget)<grid_step]
newvxvals = vxvals[abs(yvals-ytarget)<grid_step]
newvyvals = vyvals[abs(yvals-ytarget)<grid_step]
newvzvals = vzvals[abs(yvals-ytarget)<grid_step]
#newevals = evals[abs(yvals-ytarget)<grid_step]

newtempvals = newtempvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvvals = newvvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvxvals = newvxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvyvals = newvyvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvzvals = newvzvals[abs(newtvals-ttarget)<grid_step*step_fraction]
#newevals = newevals[abs(newtvals-ttarget)<grid_step*step_fraction]
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
#plt.plot(xvals[tvals==0.405], vxvals[tvals==0.405], label='MUSIC')
mathematica_data = np.fromfile('mathematica_plot.dat').reshape(-1,2)
plt.plot(mathematica_data[:,0],mathematica_data[:,1], c='g', label = 'Mathematica')
plt.plot(-mathematica_data[:,0],mathematica_data[:,1], c='g')
plt.legend()
plt.title(r'$y = 0, \tau = 0.4025 fm/c$')
#plt.yscale('log')
#plt.ylim(0,0.01)
#plt.savefig('plots/VelocityComparison.png')