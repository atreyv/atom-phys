# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:24:43 2014

@author: lab
"""


from libphys import *
import mayavi.mlab


#def expansion(x0):
#    return t: lambda np.sqrt

time = np.array([11,12,13,14,2,3,4,5,6,7,8,9,10])

def tof(time,Data):
    tof_x = np.array([])
    tof_y = np.array([])
    
    for i in range(0, len(Data)):
        p = fitgaussian2d(Data[i])
        tof_y = np.append(tof_y,gaussian_hwhm_to_stdev(np.abs(p[3])))
        tof_x = np.append(tof_x,gaussian_hwhm_to_stdev(np.abs(p[4])))
        print "Data point", i
    
    times = np.array([])
    tof_ys = np.array([])
    tof_xs = np.array([])
    for i in time.argsort():
        tof_ys = np.append(tof_ys,tof_y[i])
        tof_xs = np.append(tof_xs,tof_x[i])
    times = np.sort(time)
    return times, tof_x, tof_y