#!/usr/bin/env python3

""" This plots any two quantities from the gluon equation of state. """

import numpy as np
import matplotlib.pyplot as plt

eos_data = np.fromfile('../../EOS/PCE/gluon_eos_binary_e0p25spacing.dat').reshape(-1,5)
#cs2_data = np.fromfile('../../EOS/PCE/gluon_cs2_binary.dat').reshape(-1,2)

e = eos_data[:,0] # GeV/fm^3
p = eos_data[:,1] # GeV/fm^3
s = eos_data[:,2] # 1/fm^3
temp = eos_data[:,3] # GeV
cs2 = eos_data[:,4]

HBARC = 0.1973269718
var1 = e**4*HBARC
var2 = p/(e**4*HBARC)

plt.scatter(var1, var2, s=5)
plt.xlabel(r'$e$')
plt.ylabel(r'$p/e$')
plt.xlim(0,1)
#plt.ylim(0,.1)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0,1)
#plt.ylim(0.0001,100)
#plt.xlim(1e-7,100)
#plt.savefig('plots/p_over_e.png')