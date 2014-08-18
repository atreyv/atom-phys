# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:38:14 2013

@author: Pedro Gomes
"""

import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)

import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface
ion()                            # Turned on Matplotlib's interactive mode

from libphys import *
from scipy.signal import find_peaks_cwt
from scipy.interpolate import InterpolatedUnivariateSpline

num_files = 1
plot_x_offset = 0
half_interval = 100

plt.clf()
fig = figure(1)
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Time (s)')
#ax.set_xlabel('points')
ax.set_ylabel('PD Voltage (V)')
satabs_point_min=4200
satabs_point_max=6500
gmax=-0.0845
gfw=0.0134
goffset=.2


for i in range(0,num_files):
    f = np.loadtxt("Rb87-F=1-20130820000"+str(i)+".dat")
    f2 = FourierFilter(f[:,1], half_interval)
#    ax.plot(f[:,0] + plot_x_offset*i, f2)
#    ax.plot(f2)
    doppler_abs = np.append(f[0:satabs_point_min],f[satabs_point_max:],axis=0)
#    s = InterpolatedUnivariateSpline(doppler_abs[:,0],doppler_abs[:,1])
    gx0 = f[np.argmin(doppler_abs[:,1]),0]
#    ax.plot(doppler_abs[:,0],-doppler_abs[:,1])
#   plt.plot(f[:,0] + plot_x_offset*i, f[:,1])
    p = spleastsq(residuals_gauss,[gmax,gx0,gfw,goffset],args=(doppler_abs[:,1],doppler_abs[:,0]),xtol=10**-16)[0]
    peaks = f[:,1]-gauss(f[:,0],p)    
    f3 = FourierFilter(peaks,half_interval)    
#    ax.plot(peaks)
    ax.plot(f[:,0],peaks)
#    ax.plot(f3)
    p1max, p1x0, p1fw, p1offset = 0.02, f[np.argmax(peaks[4200:4500])+4200,0], 0.0005, 0.
    peak1 = spleastsq(residuals_gauss,[p1max,p1x0,p1fw,p1offset],args=(peaks[4200:4500],f[4200:4500,0]),xtol=10**-16)[0]
    print peak1
    ax.plot(f[4000:4700,0],gauss(f[4000:4700,0],peak1))

    p2max, p2x0, p2fw, p2offset = 0.05, f[np.argmax(peaks[4600:4900])+4600,0], 0.0005, 0.
    peak2 = spleastsq(residuals_gauss,[p2max,p2x0,p2fw,p2offset],args=(peaks[4600:4900],f[4600:4900,0]),xtol=10**-16)[0]
    print peak2
    ax.plot(f[4400:5100,0],gauss(f[4400:5100,0],peak2))

