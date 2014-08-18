# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 12:29:14 2014

@author: pedro
"""

from libphys import *


close()
fig = figure()
#fig.subplots_adjust(left=0.1,right=0.9)
ax1 = fig.add_subplot(211)
ax1.set_ylabel(r'$q_2$ position',fontsize=16)
ax2 = fig.add_subplot(212)
ax2.set_xlabel('image number',fontsize=18)#/\Lambda$")
ax2.set_ylabel(r'$q_2$ FWHM',fontsize=16)
x = np.linspace(1,104,103)
ax1.scatter(x,q1_peaks[:,1])
ax1.plot(x,q1_peaks[:,1])
ax2.scatter(x,q1_peaks[:,2])
ax2.plot(x,q1_peaks[:,2])