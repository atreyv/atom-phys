# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 19:59:46 2014

@author: pedro
"""
import os
from libphys import *

dname = '/home/pedro/LAB/DATA/2014/Aug/PF/14_08_14/pf03/'

files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
files.sort()
print 'Found %d files' %len(files)

pump_V = np.array([6,5.5,5,4.5,4,3.8,3.5,3.4,3.3,3.2,3.1,3,2.9,2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.1,2,1.9,1.8,1.7,1.6,1.5])

j=0
for i in range(20,290,10):
    if (len(files[i])>=0):
        if (files[i][0:3]=='sig'):
            for k in range(0,10):
                os.rename(dname+files[i+k],dname+files[i+k][:21]+'-'+str(pump_V[j])+files[i+k][21:])
                print 'success'
            j+=1
print j