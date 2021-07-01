#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
 
results_file = '../../evolution_xyeta2000.dat'

tvals = np.array([])
xvals = np.array([])
yvals = np.array([])
tempvals = np.array([])
evals = np.array([])
vvals = np.array([])
Tpropvals = np.array([])

xmax = 15.0
ymax = 15.0
t0 = 0.4
tmax = t0+5.0
dx = 0.06*5
dy = 0.06*5
dt = 0.025*5


with open(results_file, 'r') as f:
    for t in np.arange(t0, tmax+dt, dt):
        print(t)
        for y in np.arange(-ymax, ymax, dy):
            for x in np.arange(-xmax, xmax, dx):
                line = f.readline()
                tvals = np.append(tvals, t)
                xvals = np.append(xvals, x)
                yvals = np.append(yvals, y)
                tempvals = np.append(tempvals, 1000*float(line.split()[0]))
                evals = np.append(evals, float(line.split()[1]))
                vvals = np.append(vvals, np.sqrt(float(line.split()[3])**2 + float(line.split()[4])**2))
                Tpropvals = np.append(Tpropvals, float(line.split()[6]))
    f.close()
    
rvals = np.sqrt(xvals**2 + yvals**2)

""" Comparing to OSU-hydro data """

surface = np.fromfile('../../../osu-hydro-pce/test/surface.dat', dtype='f8').reshape(-1, 8)

OSUresults = dict(
    zip(['x','v'], np.hsplit(surface, [3,5])),
    #pi=dict(zip(['xx','xy','yy'],surface.T[11:14])),
    #Pi=surface.T[15],
    Temp = surface.T[5],
    ed = surface.T[6],
    Tprop=surface.T[7])#,
   # sd = surface.T[19],
   # intersect = surface.T[20])

OSUtvals = OSUresults['x'][:,0]
OSUxvals = OSUresults['x'][:,1]
OSUyvals = OSUresults['x'][:,2]
OSUrvals = np.sqrt(OSUxvals**2 + OSUyvals**2)
OSUtempvals = 1000*OSUresults['Temp']
OSUevals = OSUresults['ed']
OSUTpropvals = OSUresults['Tprop']
OSUvvals = np.sqrt((OSUresults['v']**2).sum(axis=1))

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
plt.figure(figsize=(7,5))

""" Initial condition tests """
#ic = np.fromfile('ed.dat').reshape(grid_n + 1, grid_n + 1)
#plt.plot(ic[51,:])

# rg = np.array([])
# for i in np.linspace(-8.9, 8.9, 90):
#     for j in np.linspace(-8.9, 8.9, 90):
#         rg = np.append(rg, np.sqrt(i**2 + j**2))

# glauber = np.zeros([90,90])
# with open('GlauberMC.dat') as f:
#     glaubervals = f.readlines()[1:]
#     f.close()
    
# for vals in glaubervals:
#     x = float(vals[5:15])
#     y = float(vals[15:25])
#     e = float(vals[25:])
#     if (x not in np.linspace(-8.9, 8.9, 90).round(2)):
#         print("x = ", x)
#         break
#     elif not (y in np.linspace(-8.9, 8.9, 90).round(2)):
#         print("y = ", y)
#         break
#     else:
#         xpos = int(round((x+8.9)/0.2))
#         ypos = int(round((y+8.9)/0.2))
#         glauber[xpos, ypos] = e

grid_step = 0.15
step_fraction = 0.1
xtarget = 0.0
ytarget = 0.0
ttarget = 1.1

newyvals = yvals[abs(yvals-ytarget)<grid_step]
newxvals = xvals[abs(yvals-ytarget)<grid_step]
newtvals = tvals[abs(yvals-ytarget)<grid_step]
newtempvals = tempvals[abs(yvals-ytarget)<grid_step]
newevals = evals[abs(yvals-ytarget)<grid_step]
newvvals = vvals[abs(yvals-ytarget)<grid_step]
newTpropvals = Tpropvals[abs(yvals-ytarget)<grid_step]

newtempvals = newtempvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newevals = newevals[abs(newtvals-ttarget)<grid_step*step_fraction]
newTpropvals = newTpropvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvvals = newvvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newxvals = newxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newtvals = newtvals[abs(newtvals-ttarget)<grid_step*step_fraction]

newOSUyvals = OSUyvals[abs(OSUyvals-ytarget)<grid_step]
newOSUxvals = OSUxvals[abs(OSUyvals-ytarget)<grid_step]
newOSUtvals = OSUtvals[abs(OSUyvals-ytarget)<grid_step]
newOSUtempvals = OSUtempvals[abs(OSUyvals-ytarget)<grid_step]
newOSUevals = OSUevals[abs(OSUyvals-ytarget)<grid_step]
newOSUvvals = OSUvvals[abs(OSUyvals-ytarget)<grid_step]
newOSUTpropvals = OSUTpropvals[abs(OSUyvals-ytarget)<grid_step]

newOSUtempvals = newOSUtempvals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]
newOSUevals = newOSUevals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]
newOSUTpropvals = newOSUTpropvals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]
newOSUvvals = newOSUvvals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]
newOSUxvals = newOSUxvals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]
newOSUtvals = newOSUtvals[abs(newOSUtvals-ttarget)<grid_step*step_fraction]


""" 2D scatter plot of temperature (tempvals) over space and time: """
plt.ylabel(r'$\tau$ (fm/c)')
plt.xlabel(r'$x$ (fm)')
plt.scatter(rvals, tvals, c=tempvals, cmap=plt.cm.jet, vmin=0, vmax=400)
#plt.scatter(OSUrvals[::50], OSUtvals[::50], c = OSUtempvals[::50], cmap=plt.cm.jet, vmin = 0, vmax = 400)
plt.colorbar(label=r'$T$ (MeV)', extend='both')
#plt.scatter(newxvals, newTpropvals, c='b', label = 'MUSIC',s=1)
#plt.scatter(np.sort(newOSUxvals), newOSUTpropvals[newOSUxvals.argsort()], c='r', label = 'OSU-hydro',s=1) 
#profile = np.load('../../../osu-hydro-pce/test/profiles/profilex15n65p0grid300central20.npy')
#plt.plot(np.linspace(-15,15,300),profile[150,:],c='g', label='TRENTO average')
# mathematica_data = np.fromfile('../../../ideal_hydro_cylindrical/mathematica_plot.dat').reshape(-1,2)
# plt.plot(mathematica_data[:,0],mathematica_data[:,1], c='g', label = 'Mathematica')
# plt.plot(-mathematica_data[:,0],mathematica_data[:,1], c='g')
plt.legend()
plt.title(r'$T_{eq} = 5 fm/c$')


""" Plotting contours """
r = np.linspace(0, 15, 50)
t = np.linspace(0, 10, 50)


# temp = griddata((rvals,tvals),tempvals,(r[None,:],t[:,None]),method='nearest')
# cs = plt.contour(r, t, temp, levels=[155., 200., 270., 350.], colors='k', linewidths = 0.5, extend='both')
# plt.clabel(cs, inline=0, fontsize=10)


""" 2D scatter plot of fugacity (fugvals) over space and time: """
# plt.xlabel('x (fm)')
# plt.ylabel(r'$\tau$ (fm/c)')
# plt.scatter(rvals, tvals, c = fugvals, cmap=plt.cm.jet, vmin = 0, vmax = 1)
# plt.colorbar(label=r'$\lambda$')


# """ Plotting contours """
# r = np.linspace(0, 15, 300)
# t = np.linspace(0, 15, 300)
# plt.title(r'$T_{eq} = 5 fm/c$')

# fug = griddata((rvals[::50],tvals[::50]),fugvals[::50],(r[None,:],t[:,None]),method='linear')
# cs = plt.contour(r, t, fug, levels=[0.25, 0.5, 0.75, 0.9], colors='k', linewidths = 0.5, extend='both')
# plt.clabel(cs, inline=0, fontsize=10)


# inter = griddata((rvals,tvals),intersect,(r[None,:],t[:,None]),method='linear')
# cs = plt.contour(r, t, inter, levels=[1], colors='w', linewidths = 0.5)
# plt.clabel(cs, inline=0, fontsize=10)


# plt.scatter(rvals, tvals, c = ((results['v']**2).sum(axis=1))**(1/2), cmap=plt.cm.jet, vmin = 0, vmax = 1)
# plt.colorbar(label=r'$|v|$')

# temp = griddata((rvals,tvals),tempvals,(r[None,:],t[:,None]),method='linear')
# cs2 = plt.contour(r, t, temp, levels=[155], colors='w', linewidths = 0.5)
# plt.clabel(cs2, inline=0, fontsize=10)


""" 1D plot of initial energy density over x: """
### Saved to plots/InitialEnergyDensity.png
#plt.xlabel('x (fm)')
#plt.ylabel(r'$\epsilon (GeV/fm^3)$')
#plt.scatter(rvals, evals, s = 1.0)
    
""" Plot of energy density in central cell: """
### Saved to plots/CentralCell.png
# plt.xlabel(r'$\tau$ (fm/c)'$)
# plt.ylabel(r'$\epsilon (GeV/fm^3)$')
# rvals_original = rvals
# rvals = rvals[rvals<2*grid_step]
# tvals = tvals[rvals_original<2*grid_step]
# evals = evals[rvals_original<2*grid_step]
# plt.plot(tvals, evals)


""" Determining and plotting total entropy per space-time rapidity """
# tspace = np.linspace(0.1, 12.1, 121)
# indices = []

# for t in tspace:
#     newindex = int(np.argmax(results['x'][:,0]>t))
#     if (newindex == 0 and t > tspace[0]):
#         break
#     indices.append(newindex)

# S = np.ones(len(indices))
# sd = results['sd']
# gamma = (1-(results['v']**2).sum(axis=1))**(-1/2)
# dtspace = tspace[1]-tspace[0]
# r = np.linspace(0, 15, 100)
# t = np.linspace(0, 15, 100)
# #entropy = griddata((rvals,tvals),sd*gamma*tvals,(r[None,:],t[:,None]),method='linear')

# for i in range(1,len(indices)):
#     S[i-1] = (0.001*gamma[indices[i-1]:indices[i]]*tvals[indices[i-1]:indices[i]]*sd[indices[i-1]:indices[i]]/(indices[i]-indices[i-1])).sum()
# S[len(indices)-1] = (0.001*gamma[indices[-1]:]*tvals[indices[-1]:]*sd[indices[-1]:]/(tvals.size - indices[-1])).sum()
    
# plt.xlabel(r'$\tau$ (fm/c)')
# plt.ylabel(r'$10^{-3} dS/d\eta$')
# plt.plot(tspace[:len(indices)], S)