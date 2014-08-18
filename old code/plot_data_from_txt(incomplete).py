# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:06:56 2014

@author: pedro
"""

#from pylab import *
#import scipy as sp

dname = '/home/pedro/LAB/'
fpump = 'pump8_profile.txt'
fprobe = 'probe8_profile.txt'

pump = np.loadtxt(dname+fpump, skiprows=1)
probe = np.loadtxt(dname+fprobe, skiprows=1)

pump = np.transpose(pump)
probe = np.transpose(probe)

fig = figure()
#fig.subplots_adjust(left=0.1,right=0.9)
ax1 = fig.add_subplot(211)
ax1.set_ylabel('$s(x)$')

ax2 = fig.add_subplot(212)
ax2.set_xlabel('$x/\Lambda (\mu m)$')
ax2.set_ylabel('$n(x)$')

ax1.plot(pump[0],pump[1])