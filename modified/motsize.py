# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 13:24:37 2014

@author: lab
"""

from libphys import *
#import mayavi.mlab
import os


# name of directory
dname = 'D:/Data/2014/October/MOT/31_10_14/patt_cam_calibr/'

# load all files
files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
print 'Found %d files' %len(files)

#load the files
moton = np.array(plt.imread(dname+files[1]),dtype=float64)
motref = np.array(plt.imread(dname+files[0]),dtype=float64)
patton = np.array(plt.imread(dname+files[3]),dtype=float64)
pattref = np.array(plt.imread(dname+files[2]),dtype=float64)

res = np.zeros((960, 1280))
Data=np.array([res],dtype=float64)
Data+=[moton-motref]
Data=np.append(Data,[patton-pattref],axis=0)

plt.clf()
tof_x = np.array([])
tof_y = np.array([])
cm_x = np.array([])
cm_y = np.array([])
yavg_data=np.array([])
yavg=np.float

for j in range(0,len(Data)): 
#bin in x and y direction
        x=linspace(0,len(Data[j][1,:]),10000)
        y=linspace(0,len(Data[j][:,1]),10000)
        
        for i in range(0,len(Data[j][1,:])):
            yavg=Data[j][:,i]
            yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data[j][:,1])) #two different ways of summing were used just because
        
        fit_x=fitgaussian1d(None,yavg_data)
        
        xavg_data=np.zeros(len(Data[j][:,1]))
    
        for k in range(0,len(Data[j][:,1])):
            xavg_data[k]=np.sum(Data[j][k,:])/ len(Data[j][1,:]) 
        
        fit_y=fitgaussian1d(None,xavg_data) 
        
        figure(1)        
        
        plot(Data[j][fit_y[1],:],label='x')
        plot(yavg_data,label='y_avg')
        plot(x,gaussian1d(fit_x[0],fit_x[1],fit_x[2],fit_x[3])(x),label='xfit')
        plt.xlabel('x (pixel)')
        plt.ylabel('Fluorescence (a.u.)')
                
        figure(2)

        plot(Data[j][:,fit_x[1]],label='y')
        plot(xavg_data,label='x_avg')
        plot(y,gaussian1d(fit_y[0],fit_y[1],fit_y[2],fit_y[3])(y),label='yfit')
        plt.xlabel('z (pixel)')
        plt.ylabel('Fluorescence (a.u.)')

        tof_y = np.append(tof_y,np.abs(fit_y[2]))
        tof_x = np.append(tof_x,np.abs(fit_x[2]))
        cm_x = np.append(cm_x,np.abs(fit_x[1]))
        cm_y = np.append(cm_y,np.abs(fit_y[1]))
        yavg_data=np.array([])
