# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 19:20:28 2014

@author: lab
"""
from pylab import *
ion()
clf()

int_mot = array([5.5,5.2,5.,4.8,4.6,4.5,4.4,4.3,4.2,4.1,4.])
pd_on = array([1.20,1.19,1.15,1.12,1.07,1.04,1.02,0.98,0.95,0.91,0.88])
pd_off = array([.6,.606,.6,.586,.566,.56,.547,.54,.53,.52,.505])

power = array([201,200,195,189,182,177,173,169,163,158,152])
#power = power * 4.5/201
# From ipynb "ATom number, 25/02/2014
atom = [671303975.673934, 655400191.007973, 626939533.837211, 620655746.064637, 599891013.199139, 581569497.250034, 581584683.944451, 549286780.774455, 536898597.773270, 508960443.451625, 502267996.142725]
atom2 = (pd_on - pd_off) * 1.118840E+09

int_rep = arange(10,4.5,-.5)
pd_on_rep = array([1.46,1.45,1.44,1.43,1.42,1.4,1.39,1.37,1.34,1.3,1.22])
pd_off_rep = array([.71,.7,.69,.7,.69,.69,.68,.67,.66,.66,.63])

atom_rep = (pd_on_rep - pd_off_rep)* 1.118840E+09
power_rep = array([2.45,2.39,2.32,2.25,2.15,2.04,1.91,1.72,1.49,1.24,0.98])

fig = figure(1)
p1 = fig.add_subplot(221)
p2 = fig.add_subplot(222)
p3 = fig.add_subplot(223)
p4 = fig.add_subplot(224)

p1.scatter(int_mot, power)
p1.set_xlabel('AOM amplitude control (V)')
p1.set_ylabel('Power available for MOT 6-beam (mW)')

p2.scatter(power,atom, color='red')
p2.scatter(power,atom2)
p2.set_xlabel('Power available for MOT 6-beam (mW)')
p2.set_ylabel('Number of atoms (empiric)')

p3.scatter(int_rep, power_rep)
p3.set_xlabel('AOM amplitude control (V)')
p3.set_ylabel('Power available for Rep 6-beam (mW)')

p4.scatter(power_rep,atom_rep)
p4.set_xlabel('Power available for Rep 6-beam (mW)')
p4.set_ylabel('Number of atoms (empiric)')
p4.set_ybound(lower=0)
p4.set_xbound(lower=0)


fit_mot = polyfit(power,atom,1)
power_fit = power
np.put(power_fit,len(power_fit)-1,0)

for i in range(0,len(atom_fit)):
    atom_fit[i] = fit_mot[0]*power_fit[i] + fit_mot[1]

p2.plot(power_fit,atom_fit,color='red')
p2.set_ybound(lower=0)
p2.set_xbound(lower=0,upper=amax(power))
