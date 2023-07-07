#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
 
QCD_flag = 0 # Toggle for comparison to equilibrium QCD result

""" Read data from MUSIC results """

results_file = '../../tests/pce/0001/evolution_xyeta.dat'

if (QCD_flag):
    results_file_QCD = '../../evolution_xyeta_QCDtestalt.dat'

xmax = 25
ymax = 25
t0 = 0.6
tmax = 10.6
dx = 0.1*2
dy = 0.1*2
dt = 0.0025*500
nx = int(2*xmax/dx)
ny = int(2*ymax/dy)
nt = int((tmax-t0)/dt + 1.01)

tvals = np.array([])
yvals = np.array([])

for t in np.arange(t0, tmax+dt, dt):
    if t > tmax*1.001:
        break
    tvals = np.append(tvals, np.repeat(t,nx*ny))

for y in np.arange(-ymax, ymax, dy):
    yvals = np.append(yvals, np.repeat(y, nx))
yvals = np.tile(yvals, nt)
        
xvals = np.tile(np.arange(-xmax,xmax,dx),ny*nt)

results = np.loadtxt(results_file)
tempvals = 1000 * results[:,0]
vxvals = results[:,3]
vyvals = results[:,4]
vzvals = results[:,5]
evals = results[:,1]
pvals = results[:,7]
Tpropvals = results[:,6]
svals = results[:,2]
   
rvals = np.sqrt(xvals**2 + yvals**2)
vvals = np.sqrt(vxvals**2 + vyvals**2)

if (QCD_flag):
    resultsQCD = np.loadtxt(results_file_QCD)
    tempvalsQCD = 1000 * resultsQCD[:,0]
    vxvalsQCD = resultsQCD[:,3]
    vyvalsQCD = resultsQCD[:,4]
    vzvalsQCD = resultsQCD[:,5]
    evalsQCD = resultsQCD[:,1]
    pvalsQCD = resultsQCD[:,7]
    TpropvalsQCD = resultsQCD[:,6]
       
    vvalsQCD = np.sqrt(vxvalsQCD**2 + vyvalsQCD**2)

""" For MUSIC mode 666 """

# results_file = '../../../music_mod/evolution_bjorken.dat'
# results = np.loadtxt(results_file)
# tvals = results[:,0]
# xvals = results[:,1]
# evals = results[:,3]
# tempvals = results[:,4]
# uxvals = results[:,5]
# uetavals = results[:,6]
# utvals = np.sqrt(1+uxvals**2)
# vxvals = np.abs(uxvals/utvals)

""" Filter to get one slice of y, t values """

grid_step = 0.00001
step_fraction = 0.000001
xtarget = 0.0
ytarget = 0.0
ttarget = 0.85

newyvals = yvals[abs(yvals-ytarget)<grid_step]
newxvals = xvals[abs(yvals-ytarget)<grid_step]
newtvals = tvals[abs(yvals-ytarget)<grid_step]
newtempvals = tempvals[abs(yvals-ytarget)<grid_step]
newvvals = vvals[abs(yvals-ytarget)<grid_step]
newvxvals = vxvals[abs(yvals-ytarget)<grid_step]
newvyvals = vyvals[abs(yvals-ytarget)<grid_step]
newvzvals = vzvals[abs(yvals-ytarget)<grid_step]
newevals = evals[abs(yvals-ytarget)<grid_step]
newpvals = pvals[abs(yvals-ytarget)<grid_step]
newTpropvals = Tpropvals[abs(yvals-ytarget)<grid_step]

newtempvals = newtempvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvvals = newvvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvxvals = newvxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvyvals = newvyvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newvzvals = newvzvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newevals = newevals[abs(newtvals-ttarget)<grid_step*step_fraction]
newpvals = newpvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newTpropvals = newTpropvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newxvals = newxvals[abs(newtvals-ttarget)<grid_step*step_fraction]
newtvals = newtvals[abs(newtvals-ttarget)<grid_step*step_fraction]

if (QCD_flag):
    newtempvalsQCD = tempvalsQCD[abs(yvals-ytarget)<grid_step]
    newvvalsQCD = vvalsQCD[abs(yvals-ytarget)<grid_step]
    newvxvalsQCD = vxvalsQCD[abs(yvals-ytarget)<grid_step]
    newvyvalsQCD = vyvalsQCD[abs(yvals-ytarget)<grid_step]
    newvzvalsQCD = vzvalsQCD[abs(yvals-ytarget)<grid_step]
    newevalsQCD = evalsQCD[abs(yvals-ytarget)<grid_step]
    newpvalsQCD = pvalsQCD[abs(yvals-ytarget)<grid_step]
    newTpropvalsQCD = TpropvalsQCD[abs(yvals-ytarget)<grid_step]
    
    newtvals = tvals[abs(yvals-ytarget)<grid_step]
    newtempvalsQCD = newtempvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newvvalsQCD = newvvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newvxvalsQCD = newvxvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newvyvalsQCD = newvyvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newvzvalsQCD = newvzvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newevalsQCD = newevalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newpvalsQCD = newpvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
    newTpropvalsQCD = newTpropvalsQCD[abs(newtvals-ttarget)<grid_step*step_fraction]
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


def fugacity(Tprop):
    return np.ones(Tprop.size) - np.exp((0.6-Tprop)/1.)

""" 2D scatter plot of temperature (tempvals) over space and time: """
plt.ylabel(r'$\tau$ (fm/c)')
plt.xlabel(r'$x$ (fm)')
plt.scatter(xvals[abs(yvals-ytarget)<grid_step], tvals[abs(yvals-ytarget)<grid_step], c=tempvals[abs(yvals-ytarget)<grid_step], cmap=plt.cm.jet, vmin=0, vmax=400,s=100)
plt.colorbar(label=r'$T$ (MeV)', extend='both')
#plt.scatter(xvals[abs(yvals-ytarget)<grid_step], tvals[abs(yvals-ytarget)<grid_step], c=vvals[abs(yvals-ytarget)<grid_step], cmap=plt.cm.jet, vmin=0, vmax=1,s=100)
#plt.colorbar(label=r'$v$', extend='both')
mathematica_data = np.fromfile('../../../ideal_hydro_cylindrical/mathematica_plot.dat').reshape(-1,2)
#plt.plot(mathematica_data[:,0],mathematica_data[:,1], c='g')
#plt.plot(-mathematica_data[:,0],mathematica_data[:,1], c='g')
#plt.legend()
#plt.title(r'$\tau_{eq} =$ 5 fm/c')
#plt.yscale('log')
#plt.ylim(0,0.01)
#plt.savefig('plots/VelocityComparison.png')

""" Plotting contours """
r = np.arange(-xmax,xmax,dx)
t = np.arange(t0, tmax+dt, dt)

#cs = plt.contour(r, t, tempvals[abs(yvals-ytarget)<grid_step].reshape(t.size,r.size), levels=[155., 200., 270., 350.], colors='k', linewidths = 1, extend='both')
#plt.clabel(cs, inline=0, fontsize=10, fmt = '%d')

#cs = plt.contour(r, t, fugacity(Tpropvals[abs(yvals-ytarget)<grid_step].reshape(t.size,r.size)), levels=[0.25,0.5,0.75,0.9], colors='k', linewidths = 1, extend='both')
#plt.clabel(cs, inline=0, fontsize=10)

if (QCD_flag):
     csQCD = plt.contour(r, t, tempvalsQCD[abs(yvals-ytarget)<grid_step].reshape(t.size,r.size), levels=[155., 200., 270., 350.], colors='k', linewidths = 0.5, linestyles='dashed', extend='both')
     plt.clabel(csQCD, inline=0, fontsize=0)
    
plt.plot(0,1,'k-',label=r'$\tau_{eq}$ = 5 fm/c')
plt.plot(0,1,'k--',label=r'$\tau_{eq}$ = 0 fm/c')
plt.legend(loc = 2, prop = {'size': 12})
#plt.xlim(-10,10)    

""" Entropy calculation """

# neweos = np.fromfile('../../EOS/PCE/PCE_eos_test4_e0p25spacing.dat').reshape(100000,100,4)

# import scipy.integrate
# from scipy.interpolate import PchipInterpolator
# fugvals = fugacity(Tpropvals)
# eosp1 = PchipInterpolator(neweos[10:,99,3],neweos[10:,99,1])
# eosp2 = PchipInterpolator(neweos[10:,0,3],neweos[10:,0,1])

# fugvals2 = fugvals
# fugvals2[fugvals2<1e-10] = 1e-10
# svals = (evals + pvals)/(tempvals/1000) - fugvals * np.log(fugvals2) * (eosp1(tempvals/1000) - eosp2(tempvals/1000))/(tempvals/1000)
# gammavals = (1-vvals**2)**(-1/2)
# grid_step = 0.001
# step_fraction = 1

# trange = np.arange(t0, tmax-t0+dt, dt)
# ds_deta = np.array([])
# for ttarget in trange:
#     print(ttarget)
#     newyvals = yvals[abs(tvals-ttarget)<grid_step*step_fraction]
#     newxvals = xvals[abs(tvals-ttarget)<grid_step*step_fraction]
#     newsvals = svals[abs(tvals-ttarget)<grid_step*step_fraction]
#     newgammavals = gammavals[abs(tvals-ttarget)<grid_step*step_fraction]
#     newtempvals = tempvals[abs(tvals-ttarget)<grid_step*step_fraction]

#     res = np.array([])
#     for x in np.arange(-xmax, xmax, dx):
#         newnewsvals = newsvals[abs(newxvals - x) < grid_step*step_fraction]
#         newnewgammavals = newgammavals[abs(newxvals - x) < grid_step*step_fraction]
#         res = np.append(res, scipy.integrate.trapz(newnewsvals*newnewgammavals,np.arange(-xmax, xmax, dx)))

#     ds_deta_val = ttarget*scipy.integrate.trapz(res,np.arange(-xmax, xmax, dx))

#     ds_deta = np.append(ds_deta, ds_deta_val)

# plt.plot(np.arange(t0, tmax-t0+dt, dt), ds_deta)