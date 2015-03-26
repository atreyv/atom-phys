# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 19:59:46 2014

@author: pedro
"""
import os
from libphys import *

dname = '/home/pedro/LAB/DATA/2014/Dec/MOT/09_12_14/abs15/'

files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
files.sort()
print 'Found %d files' %len(files)

j=0
for i in range(0,60):
    if (len(files[i])>=30):
        if (files[i][0:8]=='abs_2014'):
            os.rename(dname+files[i],dname+files[i][0:3]+'04.2'+files[i][3:])
            print 'success'
            j+=1
print j